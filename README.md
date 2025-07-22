# poolgen

Quantitative and population genetics analyses using pool sequencing data (i.e. SNP data  where each sample is a pool or group of individuals, a population or a single polyploid individual).

|**Build Status**|**License**|**Documentation**|
|:---:|:---:|:---:|
| <a href="https://github.com/jeffersonfparil/poolgen/actions"><img src="https://github.com/jeffersonfparil/poolgen/actions/workflows/rust.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) | [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffersonfparil.github.io/poolgen/poolgen/) |

## Quickstart

1. Install rustup from [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install).

2. Download this repository

    ```shell
    git clone https://github.com/jeffersonfparil/poolgen.git
    ```

3. Compile and run

    ```shell
    cd poolgen/
    cargo build --release
    ```

4. View available utilities and options

    ```shell
    ./target/release/poolgen -h
    ```

## File formats

> NOTE: Header line/s and comments should be prepended by '#'.

### Pileup

Summarised or piled up base calls of aligned reads to a reference genome.

- *Column 1*:       name of chromosome, scaffold or contig
- *Column 2*:       locus position
- *Column 3*:       reference allele
- *Column 4*:       coverage, i.e. number of times the locus was sequenced
- *Column 5*:       read codes, i.e. "." ("," for reverse strand) reference allele; "A/T/C/G" ("a/t/c/g" for reverse strand) alternative alleles; "`\[+-][0-9]+[ACGTNacgtn]`" insertions and deletions; "^[" start of read including the mapping quality score; "$" end of read; and "*" deleted or missing locus.
- *Column 6*:       base qualities encoded as the `10 ^ -((ascii value of the character - 33) / 10)`
- *Columns 7 - 3n*: coverages, reads, and base qualities of *n* pools (3 columns per pool).

Pileup from alignments can be generated similar to below:
```shell
samtools mpileup \
    -b /list/of/samtools/-/indexed/bam/or/cram/files.txt \
    -l /list/of/SNPs/in/tab/-/delimited/format/or/bed/-/like.txt \
    -d 100000 \
    -q 30 \
    -Q 30 \
    -f /reference/genome.fna \
    -o /output/file.pileup
```

### Variant call format (vcf)

Canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and genotype calls are not required since allele frequencies from allele depth will be used. The input vcf file can be generated with bcftools mpileup like: `bcftools mpileup -a AD...`. The [`vcf2sync`](#vcf2sync) utility is expected to work with vcf versions 4.2 and 4.3. See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.

### Sync

An extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.

- *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
- *Column 1*:       chromosome or scaffold name
- *Column 2*:       locus position 
- *Column 3*:       reference allele, e.g. A, T, C, G 
- *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.

### Phenotypes (required)

1. A simple delimited file, e.g. "csv" or "tsv" with a column for individual IDs, pool sizes, and at least one column for the phenotypic values.

2. GWAlpha-compatible text file (i.e. "py"):

- *Line 1*: phenotype name
- *Line 2*: standard deviation of the phenotype across pools or for the entire population
- *Line 3*: minimum phenotype value
- *Line 4*: maximum phenotype value
- *Line 5*: cummulative pool sizes percentiles (e.g. 0.2,0.4,0.6,0.8,1.0)
- *Line 6*: phenotype values corresponding to each percentile (e.g. 0.16,0.20,0.23,0.27,0.42)

## Utilities

### pileup2sync

Convert pileup (`*.pileup`, `*.mpileup`) into a synchronised pileup format (`*.sync`). 

### vcf2sync

Convert variant call format (`*.vcf`) into a synchronised pileup format (`*.sync`), using allele depths to estimate allele frequencies and omitting genotype classes information including genotype likelihoods.

### sync2csv

Convert synchronised pileup format (`*.sync`) into a matrix ($n$ pools x $p$ alleles across loci) and write into a comma-delimited (csv) file.

### fisher_exact_test

Perform Fisher's exact test per locus.

### chisq_test

Perform Chi-square test per locus.

### pearson_corr

Calculate correlations between allele frequencies per locus and phenotype data.

### ols_iter

Perform ordinary linear least squares regression between allele frequencies and phenotypes per locus, independently.

### ols_iter_with_kinship

Perform ordinary linear least squares regression between allele frequencies and phenotypes using a kinship matrix ($XX' \over p$) as a covariate per locus, independently.

### mle_iter

Perform linear regression between allele frequencies and phenotypes using maximum likelihood estimation per locus, independently.

### mle_iter_with_kinship

Perform linear regression between allele frequencies and phenotypes using maximum likelihood estimation a kinship matrix ($XX' \over p$) as a covariate per locus, independently.

### gwalpha

Perform parametric genomewide association study using pool sequencing data, i.e. pool-GWAS. Refer to [Fournier-Level, et al, 2017](https://academic.oup.com/bioinformatics/article/33/8/1246/2729762) for more details.

### genomic_prediction_cross_validation

Perform genomic prediction cross-validation using various models including ordinary least squares (OLS), ridge regression (RR), least absolute shrinkage and selection operator (LASSO), and elastic-net ([glmnet](https://glmnet.stanford.edu/articles/glmnet.html)).

### fst

Estimate $F_{st}$/fixation index between pools using unbiased estimates of heterozygosity. Mean genome-wide estimates, and per sliding (overlapping/non-overlapping) window estimates are generated.

### heterozygosity

Estimates $\pi$ or heterozygosity within populations using an unbiased method ($\pi$ or $\theta_{\pi}=4N_{e}\mu$ - similar to [Korunes & Samuk 2019](https://doi.org/10.1111/1755-0998.13326). Mean genome-wide estimates, and per sliding (overlapping/non-overlapping) window estimates are generated. 

### watterson_estimator

Estimates of Watterson's estimator of $\theta$. Mean genome-wide estimates, and per sliding (overlapping/non-overlapping) window estimates are generated. 

### tajima_d

Computes [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D) per sliding (overlapping/non-overlapping) window. Mean genome-wide estimates, and per sliding (overlapping/non-overlapping) window estimates are generated. 
