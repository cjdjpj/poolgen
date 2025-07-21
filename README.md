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
    ./target/release/poolgen -h
    ```

4. Detailed documentation

    ```shell
    cargo doc --open
    ```

## File formats

## Note

**Header line/s and comments should be prepended by '#'.**

### Pileup

Summarised or piled up base calls of aligned reads to a reference genome.

- *Column 1*:       name of chromosome, scaffold or contig
- *Column 2*:       locus position
- *Column 3*:       reference allele
- *Column 4*:       coverage, i.e. number of times the locus was sequenced
- *Column 5*:       read codes, i.e. "." ("," for reverse strand) reference allele; "A/T/C/G" ("a/t/c/g" for reverse strand) alternative alleles; "`\[+-][0-9]+[ACGTNacgtn]`" insertions and deletions; "^[" start of read including the mapping quality score; "$" end of read; and "*" deleted or missing locus.
- *Column 6*:       base qualities encoded as the `10 ^ -((ascii value of the character - 33) / 10)`
- *Columns 7 - 3n*: coverages, reads, and base qualities of *n* pools (3 columns per pool).

### Variant call format (vcf)

Canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and genotype calls are not required since allele frequencies from allele depth will be used. The input vcf file can be generated with bcftools mpileup like: `bcftools mpileup -a AD...`. The [`vcf2sync`](#vcf2sync) utility is expected to work with vcf versions 4.2 and 4.3. See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.

### Sync

An extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.

- *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
- *Column 1*:       chromosome or scaffold name
- *Column 2*:       locus position 
- *Column 3*:       reference allele, e.g. A, T, C, G 
- *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.

### Phenotypes

1. A simple delimited file, e.g. "csv" and "tsv" with a column for the individual IDs, and at least one column for the phenotypic values. Header line/s and comments should be prepended by '#'.

2. GWAlpha-compatible text file (i.e. "py"):

- *Line 1*: phenotype name
- *Line 2*: standard deviation of the phenotype across pools or for the entire population
- *Line 3*: minimum phenotype value
- *Line 4*: maximum phenotype value
- *Line 5*: cummulative pool sizes percentiles (e.g. 0.2,0.4,0.6,0.8,1.0)
- *Line 6*: phenotype values corresponding to each percentile (e.g. 0.16,0.20,0.23,0.27,0.42)

## Utilities

### pileup2sync

Convert pileup from `samtools mpileup` into a [synchronised pileup format](#Sync). Pileup from alignments can be generated similar to below:

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

### vcf2sync

Convert the most widely used genotype data format, [variant call format (`*.vcf`)](https://samtools.github.io/hts-specs/VCFv4.3.pdf) into a [synchronised pileup format](#Sync), making use of allele depths to estimate allele frequencies and omitting genotype classes information including genotype likelihoods. This utility should be compatible with vcf versions 4.2 and 4.3.

### sync2csv

Convert [synchronised pileup format](#Sync) into a matrix ($n$ pools x $p$ alleles across loci) and write into a comma-delimited (csv) file.

<!-- ### impute (redacted for now 2023-11-10)

Impute allele frequencies set to missing according to another minimum depth parameter, i.e. `--min-depth-set-to-missing`. Two imputation algorithms are currently available (a third one is in the works):

1. computationally efficient mean value imputation, and
2. adaptive linkage disequilibrium-based k-nearest neighbour imputation (an extension of [LD-kNNi](https://doi.org/10.1534/g3.115.021667)). -->

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

### ridge_iter

Perform ridge regression between allele frequencies and phenotypes per locus, independently.

### genomic_prediction_cross_validation

Perform genomic prediction cross-validation using various models including ordinary least squares (OLS), ridge regression (RR), least absolute shrinkage and selection operator (LASSO), and elastic-net ([glmnet](https://glmnet.stanford.edu/articles/glmnet.html)).

### fst

Estimate pairwise genetic differentiation between pools using unbiased estimates of heterozygosity ($\pi$ or $\theta_{\pi}=4N_{e}\mu$ - similar to [Korunes & Samuk 2019](https://doi.org/10.1111/1755-0998.13326) which assumes biallelic loci). Mean genomewide estimates, and per sliding (overlapping/non-overlapping) window estimates are generated.

### heterozygosity

Estimates per sliding window heterozygosities within populations using the unbiased method discussed [above](#fst).

### watterson_estimator

Estimates of Watterson's estimator of $\theta$ which is the expected heterozygosity given the number of polymorphic loci per sliding window (overlappining/non-overlapping).

### tajima_d

Computes [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D) per sliding (overlapping/non-overlapping) window.
