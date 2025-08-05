# poolgen

Quantitative and population genetics analyses using pool sequencing data (i.e. SNP data  where each sample is a pool or group of individuals, a population or a single polyploid individual).

|**Build Status**|**License**|
|:---:|:---:|
| <a href="https://github.com/jeffersonfparil/poolgen/actions"><img src="https://github.com/jeffersonfparil/poolgen/actions/workflows/rust.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

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

## Genotype data formats

> Note: Header line/s and comments should be prepended by '#'.

### Synchronised pileup (`.sync`)
- an extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.
- tab-delimited, one row per locus.
- *Header line*: header line including the names of the samples or pools, i.e. `#chr\tpos\tref\t<id_1>\t<id_2>\t<id_3>\t<id_4>\t<id_5>`
- *Column 1*: chromosome or scaffold name
- *Column 2*: locus position 
- *Column 3*: reference allele, e.g. A, T, C, G 
- *Column/s 4 to n*: colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.
- See [`tests/test.sync`](./tests/test.sync) for an example.

### Pileup (`.mpileup`, `.pileup`)
- summarized base calls from aligned reads at each reference genome position.
- tab-delimited, one row per locus.
- *Column 1*: Chromosome, scaffold, or contig name.
- *Column 2*: Locus position (1-based).
- *Column 3*: Reference allele.
- *Column 4*: Coverage (number of reads aligned to the locus).
- *Column 5*: Read codes:
  - `"."` / `","`: Reference allele on forward/reverse strand.
  - `"A/T/C/G"` / `"a/t/c/g"`: Alternate alleles on forward/reverse strand.
  - `+/-[0-9]+[ACGTNacgtn]`: Insertions/deletions.
  - `"^"`: Start of a read (followed by mapping quality).
  - `"$"`: End of a read.
  - `"*"`: Deleted or missing base.
- *Column 6*: Base qualities (ASCII-encoded, calculated as `10^(-((ascii - 33) / 10))`).
- *Columns 7–3n*: Per-pool coverage, reads, and base qualities (3 columns per pool).
- See [`tests/test.pileup`](./tests/test.pileup) for an example.

### Variant call format (`.vcf`)

- canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and may or may not have genotypes called (e.g. generated via bctools mpileup -a AD,DP ...). If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool.
- See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.
- The allele depth information (`AD`; i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters) is used to calculate allele frequencies.
- If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. If both `GT` and `AD` fields are present, then the `AD` field takes priority.
- See [`tests/test.vcf`](./tests/test.vcf) for an example.

## Phenotypes (required)

Standard format
-  A simple delimited file, e.g. "csv" or "tsv"
-  *Column 1*: Pool/sample IDs
-  *Column 2*: Pool sizes (or 1 for individuals)
-  *Column 3+*: Phenotypic values. (must have at least one)
-  See [`tests/test.csv`](./tests/test.csv) for an example.

GWAlpha format 
- GWAlpha compatible text file (i.e. "py"):
- *Line 1*: phenotype name
- *Line 2*: standard deviation of the phenotype across pools or for the entire population
- *Line 3*: minimum phenotype value
- *Line 4*: maximum phenotype value
- *Line 5*: cummulative pool sizes percentiles (e.g. 0.2,0.4,0.6,0.8,1.0)
- *Line 6*: phenotype values corresponding to each percentile (e.g. 0.16,0.20,0.23,0.27,0.42)

## Utilities

| Utility | Description |
| --- | --- |
|pileup2sync                        |Convert a pileup file into a synchronised pileup file (with a header row)|
|vcf2sync                           |Convert a vcf file into a synchronised pileup file (with a header row)|
|sync2csv                           |Convert a synchronised pileup file (sync format) into a csv of allele frequencies|
|fisher_exact_test                  |Perform Fisher's exact test per locus|
|chisq_test                         |Perform Chi-squared test per locus|
|pearson_corr                       |Compute Pearson's correlation between phenotypes and allele frequencies per locus|
|ols_iter                           |Compute genome-wide association (GWAS) per locus using ordinary least squares (OLS) regression|
|ols_iter_with_kinship              |Compute GWAS with OLE, but controlling for kinship|
|mle_iter                           |Compute genome-wide association (GWAS) per locus using maximum likelihood estimation (MLE)|
|mle_iter_with_kinship              |Compute GWAS with MLE, but controlling for kinship|
|gwalpha                            |Parametric allele effect estimation using Pool-seq data|
|genomic_prediction_cross_validation|Perform genomic prediction with cross-validation|
|fst                                |Compute Fst/fixation index between populations|
|heterozygosity                     |Compute heterozygosity/nucleotide diversity (π)|
|watterson_estimator                |Compute Watterson's estimator of θ|
|tajima_d                           |Compute Tajima's D|
