//! **poolgen**: quantitative and population genetics on pool sequencing (Pool-seq) data

#![allow(warnings)]
use clap::{Parser, Subcommand};
use base::{GeneralArgs, PhenotypeArgs, FilterArgs, WindowArgs, prepare_phen, prepare_filterstats, parse_valid_freq};
use std::env;
use ndarray::prelude::*;
use std::io;
mod base;
mod gp;
mod gwas;
mod popgen;
mod tables;
use base::{ChunkyReadAnalyseWrite, CrossValidation, LoadAll, Parse, SaveCsv};
use gp::{
    ols, penalise_glmnet, penalise_lasso_like, penalise_lasso_like_with_iterative_proxy_norms,
    penalise_ridge_like, penalise_ridge_like_with_iterative_proxy_norms,
};
use gwas::*;
use popgen::*;

#[derive(Subcommand)]
enum Utility {
    /// Convert a pileup file into a synchronised pileup file (with a header row)
    #[command(name = "pileup2sync")]
    Pileup2Sync { 
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
    },
    /// Convert a vcf file into a synchronised pileup file (with a header row)
    #[command(name = "vcf2sync")]
    Vcf2Sync { 
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
    },
    /// Convert a synchronised pileup file (sync format) into a csv of allele frequencies
    #[command(name = "sync2csv")]
    Sync2Csv { 
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
    },
    /// Perform Fisher's exact test per locus
    #[command(name = "fisher_exact_test")]
    FisherExactTest { 
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
    },
    /// Perform Chi-squared test per locus
    #[command(name = "chisq_test")]
    ChisqTest { 
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
    },
    /// Compute Pearson's correlation between phenotypes and allele frequencies per locus
    #[command(name = "pearson_corr")]
    PearsonCorr { 
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
    },
    /// Compute genome-wide association (GWAS) per locus using ordinary least squares (OLS) regression
    #[command(name = "ols_iter")]
    OlsIter {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        /// Generate GWAS manhattan plot and QQ plot
        #[clap(long, action, help_heading = "GWAS")]
        generate_plots: bool,
        /// Path to GFF file for extracting potential causative genes from GWAS
        #[clap(long, help_heading = "GWAS")]
        fname_gff: Option<String>,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action, help_heading = "GWAS")]
        output_sig_snps_only: bool,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500, help_heading = "GWAS")]
        gff_window_size: u64,
    },
    /// Compute GWAS with OLE, but controlling for kinship
    #[command(name = "ols_iter_with_kinship")]
    OlsIterWithKinship {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        /// Generate GWAS manhattan plot and QQ plot
        #[clap(long, action, help_heading = "GWAS")]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long, help_heading = "GWAS")]
        fname_gff: Option<String>,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500, help_heading = "GWAS")]
        gff_window_size: u64,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action, help_heading = "GWAS")]
        output_sig_snps_only: bool,
        /// GWAS iterative OLS using some of the kinship matrix's PCs as covariate
        #[clap(short, long, default_value_t = 0.75, value_parser = parse_valid_freq, help_heading = "GWAS")]
        xxt_eigen_variance_explained: f64,
    },
    /// Compute genome-wide association (GWAS) per locus using maximum likelihood estimation (MLE)
    #[command(name = "mle_iter")]
    MleIter {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        /// Generate GWAS manhattan plot and QQ plot
        #[clap(long, action, help_heading = "GWAS")]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long, help_heading = "GWAS")]
        fname_gff: Option<String>,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action, help_heading = "GWAS")]
        output_sig_snps_only: bool,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500, help_heading = "GWAS")]
        gff_window_size: u64,
    },
    /// Compute GWAS with MLE, but controlling for kinship
    #[command(name = "mle_iter_with_kinship")]
    MleIterWithKinship {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        /// Generate GWAS manhattan plot and QQ plot
        #[clap(long, action, help_heading = "GWAS")]
        generate_plots: bool,
        /// GFF file for extracting potential causative genes from GWAS
        #[clap(long, help_heading = "GWAS")]
        fname_gff: Option<String>,
        /// GFF window size (look for genes within gff_window_size of a significant SNP)
        #[clap(long, default_value_t = 500, help_heading = "GWAS")]
        gff_window_size: u64,
        /// Retain only SNPs that have significant p-value (less than 0.05 bonferonni corrected)
        #[clap(long, action, help_heading = "GWAS")]
        output_sig_snps_only: bool,
        /// GWAS iterative OLS using some of the kinship matrix's PCs as covariate
        #[clap(short, long, default_value_t = 0.75, value_parser = parse_valid_freq, help_heading = "GWAS")]
        xxt_eigen_variance_explained: f64,
    },
    /// Parametric allele effect estimation using Pool-seq data
    #[command(name = "gwalpha")]
    Gwalpha {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        /// GWAlpha inference method to use: "LS" for least squares or "ML" for maximum likelihood estimation
        #[clap(long, default_value = "ML", help_heading = "GWAlpha")]
        gwalpha_method: String,
    },
    /// Perform genomic prediction with cross-validation
    #[command(name = "genomic_prediction_cross_validation")]
    GenomicPredictionCrossValidation {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        /// Number of replicates of k-fold cross-validation
        #[clap(long, default_value_t = 3, help_heading = "Genomic Prediction Cross Validation")]
        n_reps: usize,
        /// Genomic prediction cross-validation: number of k-fold validations, i.e. number of time the data will be partitioned for training and testing each model
        #[clap(long, default_value_t = 10, help_heading = "Genomic Prediction Cross Validation")]
        k_folds: usize,
    },
    /// Compute Fst/fixation index between populations
    #[command(name = "fst")]
    Fst {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        #[command(flatten)]
        window: WindowArgs,
    },
    /// Compute heterozygosity/nucleotide diversity (π)
    #[command(name = "heterozygosity")]
    Heterozygosity {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        #[command(flatten)]
        window: WindowArgs,
    },
    /// Compute Watterson's estimator of θ
    #[command(name = "watterson_estimator")]
    WattersonEstimator {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        #[command(flatten)]
        window: WindowArgs,
    },
    /// Compute Tajima's D
    #[command(name = "tajima_d")]
    TajimaD {
        #[command(flatten)]
        general_args: GeneralArgs,
        #[command(flatten)]
        phenotype_args: PhenotypeArgs,
        #[command(flatten)]
        filter_args: FilterArgs,
        #[command(flatten)]
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
    #[clap(subcommand)]
    utility: Utility,
}

/// # poolgen: quantitative and population genetics on pool sequencing (Pool-seq) data
/// - *pileup2sync* - convert a pileup file into a synchronised pileup file (with a header row)
/// - *vcf2sync* - convert a vcf file into a synchronised pileup file (with a header row)
/// - *sync2csv* - convert a synchronised pileup file into a csv of allele frequencies
/// - *fisher_exact_test* - perform Fisher's exact test per locus
/// - *chisq_test* - perform Chi-squared test per locus
/// - *pearson_corr* - compute Pearson's correlation between phenotypes and allele frequencies per locus
/// - *ols_iter* - compute genome-wide association (GWAS) using ordinary least squares (OLS) regression
/// - *ols_iter_with_kinship* - compute GWAS with OLE, but controlling for kinship
/// - *mle_iter* - compute genome-wide association (GWAS) using maximum likelihood estimation (MLE)
/// - *mle_iter_with_kinship* - compute GWAS with MLE, but controlling for kinship
/// - *gwalpha* - parametric allele effect estimation using Pool-seq data
/// - *genomic_prediction_cross_validation* - perform genomic prediction with cross-validation
/// - *fst* - compute Fst/fixation index between populations
/// - *heterozygosity* - compute heterozygosity/nucleotide diversity (π)
/// - *watterson_estimator* - compute Watterson's estimator of θ
/// - *tajima_d* - compute Tajima's D
///
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

    let mut output = String::new();

    match cli.utility {
        Utility::Pileup2Sync { general_args, phenotype_args, filter_args } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::Vcf2Sync { general_args, phenotype_args, filter_args } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::Sync2Csv { general_args, phenotype_args, filter_args } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::FisherExactTest { general_args, phenotype_args, filter_args } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::ChisqTest { general_args, phenotype_args, filter_args } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::PearsonCorr { general_args, phenotype_args, filter_args } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::OlsIter { general_args, phenotype_args, filter_args, generate_plots, fname_gff, gff_window_size, output_sig_snps_only } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::OlsIterWithKinship { general_args, phenotype_args, filter_args, generate_plots, fname_gff, gff_window_size, output_sig_snps_only, xxt_eigen_variance_explained } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::MleIter { general_args, phenotype_args, filter_args, generate_plots, fname_gff, gff_window_size, output_sig_snps_only} => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::MleIterWithKinship { general_args, phenotype_args, filter_args, generate_plots, fname_gff, gff_window_size, output_sig_snps_only, xxt_eigen_variance_explained } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::Gwalpha { general_args, phenotype_args, filter_args, gwalpha_method } => {
            let file_phen = prepare_phen(&phenotype_args, "gwalpha_fmt".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "gwalpha".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            if gwalpha_method == "LS".to_owned() {
                output = format!("FILE CREATED: {}", file_sync_phen
                    .read_analyse_write(
                        &filter_stats,
                        &general_args.output,
                        &general_args.n_threads,
                        gwas::gwalpha_ls,
                    )
                    .unwrap())
            } else {
                output = format!("FILE CREATED: {}", file_sync_phen
                    .read_analyse_write(
                        &filter_stats,
                        &general_args.output,
                        &general_args.n_threads,
                        gwas::gwalpha_ml,
                    )
                    .unwrap())
            }
        }
        Utility::GenomicPredictionCrossValidation { general_args, phenotype_args, filter_args, n_reps, k_folds} => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::Fst { general_args, phenotype_args, filter_args, window } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::Heterozygosity { general_args, phenotype_args, filter_args, window } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::WattersonEstimator { general_args, phenotype_args, filter_args, window } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
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
        Utility::TajimaD { general_args, phenotype_args, filter_args, window } => {
            let file_phen = prepare_phen(&phenotype_args, "default".to_string());
            let phen = file_phen.lparse().unwrap();
            let filter_stats = prepare_filterstats(&filter_args, &phen);
            let file_sync = base::FileSync {
                filename: general_args.fname.clone(),
                test: "tajima_d".to_string()
            };
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &general_args.n_threads)
                .unwrap();
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

