//! **poolgen**: quantitative and population genetics on pool sequencing (Pool-seq) data

use clap::{Parser, Subcommand, Args, ValueEnum};
use std::env;
use ndarray::prelude::*;
#[allow(warnings)]
use std::io;
mod base;
mod gp;
mod gwas;
mod popgen;
mod tables;
use base::{ChunkyReadAnalyseWrite, CrossValidation, LoadAll, Parse, SaveCsv, helpers::parse_valid_freq};
use gp::{
    ols, penalise_glmnet, penalise_lasso_like, penalise_lasso_like_with_iterative_proxy_norms,
    penalise_ridge_like, penalise_ridge_like_with_iterative_proxy_norms,
};
use gwas::*;
use popgen::*;

#[derive(Args)]
#[clap(next_help_heading = "General", display_order = 0)]
struct GeneralArgs {
    /// Filename of the input pileup or synchronised pileup file (i.e. *.pileup, *.sync, *.syncf, or *.syncx)
    #[clap(short, long, global = true)]
    fname: String,
    /// Output filename
    #[clap(short, long, global = true, default_value = "")]
    output: String,
    /// Number of threads to use for parallel processing
    #[clap(long, default_value_t = 1, global = true)]
    n_threads: usize,
}

#[derive(Args)]
#[clap(next_help_heading = "Phenotype", display_order = 2)]
struct PhenotypeArgs {
    /// Input phenotype file: csv or tsv or any delimited file
    #[clap(short, long, global = true)]
    phen_fname: String,
    /// Delimiter of the input phenotype file: comma, tab, etc...
    #[clap(long, default_value = ",", global = true)]
    phen_delim: String,
    /// Column index containing the names or IDs of the indivudals in the input phenotype file: 0, 1, 2, ...
    #[clap(long, default_value_t = 0, global = true)]
    phen_name_col: usize,
    /// Column index containing the the sizes of each pool or population: 0, 1, 2, ...
    #[clap(long, default_value_t = 1, global = true)]
    phen_pool_size_col: usize,
    /// Column indexes containing the phenotype values in the input phenotype file, e.g. 1 or 1,2,3 or 1,2,3,4 etc ...
    #[clap(
        long,
        use_value_delimiter = true,
        value_delimiter = ',',
        value_parser = clap::value_parser!(usize),
        default_value = "2",
        global = true
    )]
    phen_value_col: Vec<usize>,
}


#[derive(Args)]
#[clap(next_help_heading = "Filtering", display_order = 1)]
struct FilterArgs {
    /// Categorize lowercase reference reads in pileup as unclassified ('N')
    #[clap(long, action, global = true,  help_heading = "Filtering options")]
    keep_lowercase_reference: bool,
    /// Keep ambiguous reads during SNP filtering, i.e. keep them coded as Ns
    #[clap(long, action, global = true, help_heading = "Filtering options")]
    keep_ns: bool,
    /// Sync to csv file conversion to include all alleles or just p-1 excluding the minimum allele
    #[clap(long, action, global = true, help_heading = "Filtering options")]
    keep_p_minus_1: bool,
    /// Maximum base sequencing error rate
    #[clap(long, default_value_t = 0.01, value_parser = parse_valid_freq, global = true, help_heading = "Filtering options")]
    max_base_error_rate: f64,
    /// Minimum breadth of coverage (loci with less than this proportion of pools below min_coverage_depth will be omitted)
    #[clap(long, default_value_t = 1.0, value_parser = parse_valid_freq, global = true, help_heading = "Filtering options")]
    min_coverage_breadth: f64,
    /// Minimum depth of coverage (loci with less than min_coverage_breadth pools below this threshold will be omitted)
    #[clap(long, default_value_t = 1, global = true, help_heading = "Filtering options")]
    min_coverage_depth: u64,
    /// Minimum allele frequency (per locus, alleles which fail to pass this threshold will be omitted allowing control over multiallelic loci)
    #[clap(long, default_value_t = 0.001, value_parser = parse_valid_freq, global = true, help_heading = "Filtering options")]
    min_allele_frequency: f64,
    /// Maximum missingness rate (loci with missing data beyond this threshold will be omitted)
    #[clap(long, default_value_t = 0.0, value_parser = parse_valid_freq, global = true, help_heading = "Filtering options")]
    max_missingness_rate: f64,
}

#[derive(Args)]
#[clap(next_help_heading = "Window", display_order = 3)]
struct WindowArgs {
    /// Estimation of population genetics parameters per window, i.e. fst, pi, Watterson's theta, and Tajima's D per population per window: window size in terms of number of bases
    #[clap(long, default_value_t = 100)]
    window_size_bp: u64,
    /// Number of bases to slide the window (a good start will be half the window size)
    #[clap(long, default_value_t = 50)]
    window_slide_size_bp: u64,
    /// Estimation of population genetics parameters per window, i.e. fst, pi, Watterson's theta, and Tajima's D per population per window: minimum number of loci per window
    #[clap(long, default_value_t = 10)]
    min_loci_per_window: u64,
}

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "UPPER")]
enum GwalphaMethod {
    LS,
    ML,
}

#[derive(Subcommand)]
enum Utility {
    /// Convert a pileup file into a synchronised pileup file with a header
    #[clap(name = "pileup2sync")]
    Pileup2Sync { },
    /// Convert a vcf file into a synchronised pileup file with a header
    #[clap(name = "vcf2sync")]
    Vcf2Sync { },
    /// Convert a synchronised pileup file (sync format) into a smple comma-separated file (csv format)
    #[clap(name = "sync2csv")]
    Sync2Csv { },
    /// Perform Fisher's exact test per locus (sync format)
    #[clap(name = "fisher_exact_test")]
    FisherExactTest { },
    /// Perform Chi-squared test per locus (sync format)
    #[clap(name = "chisq_test")]
    ChisqTest { },
    /// Find the correlation between phenotypes (csv format) and allele frequencies (sync format) per locus
    #[clap(name = "pearson_corr")]
    PearsonCorr { },
    /// Find the strength of association between phenotypes (csv format) and allele frequencies (sync format) per locus using ordinary least squares (OLS) regression
    #[clap(name = "ols_iter")]
    OlsIter { 
        /// Generate plots (relevant to analysis tool)
        #[clap(long, action)]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long)]
        fname_gff: Option<String>,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action)]
        output_sig_snps_only: bool,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500)]
        gff_window_size: u64,
    },
    /// Similar to ols_iter but controls for the effects of kinship
    #[clap(name = "ols_iter_with_kinship")]
    OlsIterWithKinship {
        /// Generate plots (relevant to analysis tool)
        #[clap(long, action)]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long)]
        fname_gff: Option<String>,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500)]
        gff_window_size: u64,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action)]
        output_sig_snps_only: bool,
        /// GWAS iterative OLS using some of the kinship matrix's PCs as covariate
        #[clap(short, long, default_value_t = 0.75, value_parser = parse_valid_freq)]
        xxt_eigen_variance_explained: f64,
    },
    /// Similar to ols_iter but uses maximum likelihood estimation
    #[clap(name = "mle_iter")]
    MleIter {
        /// Generate plots (relevant to analysis tool)
        #[clap(long, action)]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long)]
        fname_gff: Option<String>,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action)]
        output_sig_snps_only: bool,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500)]
        gff_window_size: u64,
    },
    /// Similar to ols_iter_with_kinship but uses maximum likelihood estimation
    #[clap(name = "mle_iter_with_kinship")]
    MleIterWithKinship {
        /// Generate plots (relevant to analysis tool)
        #[clap(long, action)]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long)]
        fname_gff: Option<String>,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500)]
        gff_window_size: u64,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action)]
        output_sig_snps_only: bool,
        /// GWAS iterative OLS using some of the kinship matrix's PCs as covariate
        #[clap(short, long, default_value_t = 0.75, value_parser = parse_valid_freq)]
        xxt_eigen_variance_explained: f64,
    },
    /// Parametric allele effect estimation using Pool-seq data
    #[clap(name = "gwalpha")]
    Gwalpha {
        /// GWAlpha inference method to use: "LS" for least squares or "ML" for maximum likelihood estimation
        #[clap(long, value_enum, default_value = "ML")]
        gwalpha_method: GwalphaMethod,
    },
    /// Perform genomic prediction with cross-validation
    #[clap(name = "genomic_prediction_cross_validation")]
    GenomicPredictionCrossValidation {
        /// Genomic prediction cross-validation: number of replicates of k-fold cross-validation
        #[clap(long, default_value_t = 3)]
        n_reps: usize,
        /// Genomic prediction cross-validation: number of k-fold validations, i.e. number of time the data will be partitioned for training and testing each model
        #[clap(long, default_value_t = 10)]
        k_folds: usize,
    },
    /// Find the pairwise differentiation between populations using multiple/genome-wide allele frequencies
    #[clap(name = "fst")]
    Fst {
        #[clap(flatten)]
        window: WindowArgs,
    },
    /// Find the heterozygosity or nucleotide diversity (π) of each population
    #[clap(name = "heterozygosity")]
    Heterozygosity {
        #[clap(flatten)]
        window: WindowArgs,
    },
    /// Calculate Watterson's estimator
    #[clap(name = "watterson_estimator")]
    WattersonEstimator {
        #[clap(flatten)]
        window: WindowArgs,
    },
    /// Calculate Tajima's D
    #[clap(name = "tajima_d")]
    TajimaD {
        #[clap(flatten)]
        window: WindowArgs,
    }
}

#[derive(Parser)]
#[clap(
    author = "Jeff Paril",
    version = "0.1.0",
    about = "Quantitative and population genetics analyses using pool sequencing data."
)]
struct Cli {
    #[clap(flatten)]
    general_args: GeneralArgs,

    #[clap(flatten)]
    filter_args: FilterArgs,

    #[clap(flatten)]
    phenotype_args: PhenotypeArgs,

    #[clap(subcommand)]
    utility: Utility,
}

/// # poolgen: quantitative and population genetics on pool sequencing (Pool-seq) data
/// - *pileup2sync* - convert a pileup file into a synchronised pileup file with a header
/// - *vcf2sync* - convert a vcf file into a synchronised pileup file with a header
/// - *sync2csv* - convert a synchronised pileup file (sync format) into a smple comma-separated file (csv format)
/// - *fisher_exact_test* - perform Fisher's exact test per locus (sync format), i.e. using allele counts matrix of n-pools x p-alleles (Note: phenotype data is only required for the pool sizes and/or pool names)
/// - *chisq_test* - perform $\Chi^2$ test per locus (sync format), i.e. using allele counts matrix of n-pools x p-alleles (Note: likewise phenotype data is only required for the pool sizes and/or pool names)
/// - *fst* - find the pairwise differentiation between populations using multiple/genome-wide allele frequencies (sync format) with unbiased Fst (fixation index) estimates from multiallelic loci (similar to [Gautier et al, 2019](https://doi.org/10.1111/1755-0998.13557))
/// - *heterozygosity* - find the heterozygosity or nucleotide diversity ($\pi$) of each population using multiple/genome-wide allele frequencies (sync format) with unbiased Fst (fixation index) estimates from multiallelic loci (similar to [Korunes & Samuk, 2021](https://doi.org/10.1111/1755-0998.13326))
/// - *pearson_corr* - find the correlation between phenotypes (csv format) and allele frequencies (sync format) per locus
/// - *ols_iter* - find the strength of association between phenotypes (csv format) and allele frequencies (sync format) per locus using ordinary least squares (OLS) regression
/// - *ols_iter_with_kinship* - similar to *ols_iter* but controls for the effects of kinship
/// - *mle_iter* - similar to *ols_iter* but uses maximum likelihood estimation
/// - *mle_iter_with_kinship* - similar to *ols_iter_with_kinship* but uses maximum likelihood estimation
/// - *gwalpha* - parametric allele effect estimation using Pool-seq data in cases where the number of pools is small, e.g. 5 pools genotyped across thousands to hundred of thousands of loci. See [Fournier-Level et al, 2017](https://doi.org/10.1093/bioinformatics/btw805) for details.
/// - *genomic_prediction_cross_validation* - perform *k*-fold cross-validation with *r* replicates using genomic prediction models (i.e. OLS and various penalised regression models)
///  
/// Please refer to the documentation of each module for more details.
///
/// ## Examples
/// ```shell
/// cargo run -- pileup2sync -f ./tests/test.pileup -p ./tests/test.csv
/// cargo run -- fisher_exact_test -f ./tests/test.sync -p ./tests/test.csv --n-threads 32 --min-coverage-depth 10 --min-allele-frequency 0.01
/// cargo run -- chisq_test -f ./tests/test.sync -p ./tests/test.csv --n-threads 32 --min-coverage-depth 10 --min-allele-frequency 0.01
/// cargo run -- pearson_corr -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --min-coverage-depth 10 --min-allele-frequency 0.01
/// cargo run -- fst -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32
/// cargo run -- heterozygosity -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32
/// cargo run -- ols_iter -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --min-coverage-depth 10 --min-allele-frequency 0.01
/// cargo run -- mle_iter -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --min-coverage-depth 10 --min-allele-frequency 0.01
/// cargo run -- gwalpha  -f ./tests/test.sync -p ./tests/test.py --n-threads 32 --gwalpha-method ML
/// cargo run -- sync2csv -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --keep-p-minus-1
/// # cargo run -- genomic_prediction_cross_validation -f ./tests/test_MORE_POOLS.sync -p ./tests/test_MORE_POOLS.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32
/// ```
fn main() {
    env::set_var("RUST_BACKTRACE", "1");
    let cli = Cli::parse();

    let general_args = &cli.general_args;
    let phenotype_args = &cli.phenotype_args;
    let filter_args = &cli.filter_args;

    let mut output = String::new();

    let phen_format = if matches!(cli.utility, Utility::Gwalpha { .. }) {
        "gwalpha_fmt"
    } else {
        "default"
    }.to_string();

    let phen_col = phenotype_args.phen_value_col.clone();

    let file_phen = base::FilePhen {
        filename: phenotype_args.phen_fname.clone(),
        delim: phenotype_args.phen_delim.clone(),
        names_column_id: phenotype_args.phen_name_col,
        sizes_column_id: phenotype_args.phen_pool_size_col,
        trait_values_column_ids: phen_col,
        format: phen_format,
    };

    let phen = file_phen.lparse().unwrap();

    let filter_stats = base::FilterStats {
        remove_ns: !filter_args.keep_ns,
        keep_lowercase_reference: filter_args.keep_lowercase_reference,
        max_base_error_rate: filter_args.max_base_error_rate,
        min_coverage_breadth: filter_args.min_coverage_breadth,
        min_coverage_depth: filter_args.min_coverage_depth,
        min_allele_frequency: filter_args.min_allele_frequency,
        max_missingness_rate: filter_args.max_missingness_rate,
        pool_sizes: phen.pool_sizes.clone(),
    };

    match cli.utility {
        Utility::Pileup2Sync { } => {
            let file_pileup = base::FilePileup {
                filename: general_args.fname.clone(),
                pool_names: phen.pool_names,
            };
            output = format!("FILE CREATED: {}", file_pileup
                .read_analyse_write(
                    &filter_stats,
                    &general_args.output,
                    &general_args.n_threads,
                    base::pileup_to_sync,
                )
                .unwrap());
        }
        Utility::Vcf2Sync { } => {
            let file_vcf = base::FileVcf {
                filename: general_args.fname.clone(),
            };
            output = format!("FILE CREATED: {}", file_vcf
                .read_analyse_write(
                    &filter_stats,
                    &general_args.output,
                    &general_args.n_threads,
                    base::vcf_to_sync,
                )
                .unwrap());
        }
        Utility::Sync2Csv { } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "sync2csv".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = format!("FILE CREATED: {}", file_sync_phen
                .write_csv(
                    &filter_stats,
                    filter_args.keep_p_minus_1,
                    &general_args.output,
                    &general_args.n_threads,
                )
                .unwrap());
        }
        Utility::FisherExactTest { } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "fisher_exact_test".to_string()
            };
            output = format!("FILE CREATED: {}", file_sync
                .read_analyse_write(
                    &filter_stats, 
                    &general_args.output, 
                    &general_args.n_threads, 
                    tables::fisher)
                .unwrap());
        }
        Utility::ChisqTest { } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "chisq_test".to_string()
            };
            output = format!("FILE CREATED: {}", file_sync
                .read_analyse_write(
                    &filter_stats, 
                    &general_args.output, 
                    &general_args.n_threads, 
                    tables::chisq)
                .unwrap());
        }
        Utility::PearsonCorr { } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "pearson_corr".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = format!("FILE CREATED: {}", file_sync_phen
                .read_analyse_write(
                    &filter_stats,
                    &general_args.output,
                    &general_args.n_threads,
                    gwas::correlation,
                )
                .unwrap());
        }
        Utility::OlsIter { generate_plots, fname_gff, gff_window_size, output_sig_snps_only } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "ols_iter".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = format!("FILE CREATED: {}", file_sync_phen
                .read_analyse_write(
                    &filter_stats,
                    &general_args.output,
                    &general_args.n_threads,
                    gwas::ols_iterate,
                )
                .unwrap());

            let mut python_scripts: Vec<(&str, Vec<String>)> = Vec::new();

            if generate_plots {
                python_scripts.push(("plot_manhattan.py", vec![]));
                python_scripts.push(("plot_qq.py", vec![]));
            }
            if let Some(gff_filename) = fname_gff {
                let window_size_str = gff_window_size.to_string();
                python_scripts.push(("extract_snps_from_gff.py", vec![gff_filename.clone(), window_size_str]));
            }
            if output_sig_snps_only {
                python_scripts.push(("remove_insig_snps.py", vec![]));
            }
            if !python_scripts.is_empty() {
                output = base::run_python_and_append(&output.clone(), &python_scripts);
            }
        }
        Utility::OlsIterWithKinship {  generate_plots, fname_gff, gff_window_size, output_sig_snps_only, xxt_eigen_variance_explained } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "ols_iter_with_kinship".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let mut genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, filter_args.keep_p_minus_1, &general_args.n_threads)
                .unwrap();
            output = format!("FILE CREATED: {}", ols_with_covariate(
                &mut genotypes_and_phenotypes,
                xxt_eigen_variance_explained,
                &general_args.fname,
                &general_args.output,
            )
            .unwrap());
            let mut python_scripts: Vec<(&str, Vec<String>)> = Vec::new();

            if generate_plots {
                python_scripts.push(("plot_manhattan.py", vec![]));
                python_scripts.push(("plot_qq.py", vec![]));
            }
            if let Some(gff_filename) = fname_gff {
                let window_size_str = gff_window_size.to_string();
                python_scripts.push(("extract_snps_from_gff.py", vec![gff_filename.clone(), window_size_str]));
            }
            if output_sig_snps_only {
                python_scripts.push(("remove_insig_snps.py", vec![]));
            }
            if !python_scripts.is_empty() {
                output = base::run_python_and_append(&output.clone(), &python_scripts);
            }
        }
        Utility::MleIter { generate_plots, fname_gff, gff_window_size, output_sig_snps_only} => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "mle_iter".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = format!("FILE CREATED: {}", file_sync_phen
                .read_analyse_write(
                    &filter_stats,
                    &general_args.output,
                    &general_args.n_threads,
                    gwas::mle_iterate,
                )
                .unwrap());
            let mut python_scripts: Vec<(&str, Vec<String>)> = Vec::new();

            if generate_plots {
                python_scripts.push(("plot_manhattan.py", vec![]));
                python_scripts.push(("plot_qq.py", vec![]));
            }
            if let Some(gff_filename) = fname_gff {
                let window_size_str = gff_window_size.to_string();
                python_scripts.push(("extract_snps_from_gff.py", vec![gff_filename.clone(), window_size_str]));
            }
            if output_sig_snps_only {
                python_scripts.push(("remove_insig_snps.py", vec![]));
            }
            if !python_scripts.is_empty() {
                output = base::run_python_and_append(&output.clone(), &python_scripts);
            }
        }
        Utility::MleIterWithKinship {  generate_plots, fname_gff, gff_window_size, output_sig_snps_only, xxt_eigen_variance_explained } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "mle_iter_with_kinship".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let mut genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, filter_args.keep_p_minus_1, &general_args.n_threads)
                .unwrap();
            output = format!("FILE CREATED: {}", mle_with_covariate(
                &mut genotypes_and_phenotypes,
                xxt_eigen_variance_explained,
                &general_args.fname,
                &general_args.output,
            )
            .unwrap());
            let mut python_scripts: Vec<(&str, Vec<String>)> = Vec::new();

            if generate_plots {
                python_scripts.push(("plot_manhattan.py", vec![]));
                python_scripts.push(("plot_qq.py", vec![]));
            }
            if let Some(gff_filename) = fname_gff {
                let window_size_str = gff_window_size.to_string();
                python_scripts.push(("extract_snps_from_gff.py", vec![gff_filename.clone(), window_size_str]));
            }
            if output_sig_snps_only {
                python_scripts.push(("remove_insig_snps.py", vec![]));
            }
            if !python_scripts.is_empty() {
                output = base::run_python_and_append(&output.clone(), &python_scripts);
            }
        }
        Utility::Gwalpha { gwalpha_method } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "gwalpha".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            match gwalpha_method {
                GwalphaMethod::LS => {
                    output = format!("FILE CREATED: {}", file_sync_phen
                        .read_analyse_write(
                            &filter_stats,
                            &general_args.output,
                            &general_args.n_threads,
                            gwas::gwalpha_ls,
                        )
                        .unwrap())
                },
                GwalphaMethod::ML => {
                    output = format!("FILE CREATED: {}", file_sync_phen
                        .read_analyse_write(
                            &filter_stats,
                            &general_args.output,
                            &general_args.n_threads,
                            gwas::gwalpha_ml,
                        )
                        .unwrap())
                },
            }
        }
        Utility::GenomicPredictionCrossValidation { n_reps, k_folds} => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "genomic_prediction_cross_validation".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, filter_args.keep_p_minus_1, &general_args.n_threads)
                .unwrap();
            let functions: Vec<
                fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
            > = vec![
                ols,
                penalise_glmnet,
                penalise_lasso_like,
                penalise_ridge_like,
                penalise_lasso_like_with_iterative_proxy_norms,
                penalise_ridge_like_with_iterative_proxy_norms,
            ];
            let prediction_performances = genotypes_and_phenotypes
                .cross_validate(k_folds, n_reps, functions.clone())
                .unwrap();
            let (tabulated, _pred_v_expe, predictor_files) = genotypes_and_phenotypes
                .tabulate_predict_and_output(
                    &prediction_performances,
                    functions,
                    &general_args.fname,
                    &general_args.output,
                )
                .unwrap();
            output = tabulated;
            let message = "Predictors for each model are here:\n-".to_owned()
                + &predictor_files.join("\n-")[..];
            println!("{:?}", message);
        }
        Utility::Fst { window } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "fst".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, filter_args.keep_p_minus_1, &general_args.n_threads)
                .unwrap();
            let (genome_wide, per_window) = fst(
                &genotypes_and_phenotypes,
                &window.window_size_bp,
                &window.window_slide_size_bp,
                &window.min_loci_per_window,
                &general_args.fname,
                &general_args.output,
            )
            .unwrap();
            output = format!("FILE CREATED: {}", genome_wide + " and " + &per_window[..]);
        }
        Utility::Heterozygosity { window } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "heterozygosity".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &general_args.n_threads)
                .unwrap();
            output = format!("FILE CREATED: {}", pi(
                &genotypes_and_phenotypes,
                &window.window_size_bp,
                &window.window_slide_size_bp,
                &window.min_loci_per_window,
                &general_args.fname,
                &general_args.output,
            )
            .unwrap());
        }
        Utility::WattersonEstimator { window } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "watterson_estimator".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &general_args.n_threads)
                .unwrap();
            output = format!("FILE CREATED: {}", watterson_estimator(
                &genotypes_and_phenotypes,
                &file_sync_phen.pool_sizes,
                &window.window_size_bp,
                &window.window_slide_size_bp,
                &window.min_loci_per_window,
                &general_args.fname,
                &general_args.output,
            )
            .unwrap());
        }
        Utility::TajimaD { window } => {
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "tajima_d".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &general_args.n_threads)
                .unwrap(); // we need all alleles in each locus
            output = format!("FILE CREATED: {}", tajima_d(
                &genotypes_and_phenotypes,
                &file_sync_phen.pool_sizes,
                &window.window_size_bp,
                &window.window_slide_size_bp,
                &window.min_loci_per_window,
                &general_args.fname,
                &general_args.output,
            )
            .unwrap());
        }
    }

    println!("{}", output);
}

