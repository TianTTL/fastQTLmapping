# fastQTLmapping
`FastQTLmapping` is a computationally efficient, exact, and generic solver for exhaustive multiple regression analysis involving extraordinarily large numbers of dependent and explanatory variables with covariates, which is particularly helpful in QTL-like analysis. `FastQTLmapping` can afford omics data containing tens of thousands of individuals and billions of molecular loci.



## Implementation

`FastQTLmapping` accepts input files in text format and in Plink binary format. The output file is in text format and contains all test statistics for all regressions, with the ability to control the volume of the output at preset significance thresholds. 

Different thresholds can be specified according to physical distances between the markers under investigation, which facilitates the analysis of cis- and trans-mQTLs. 

Z- and rank-normalizations are optional for pre-processing certain or all input variables. 

FastQTLmapping is deployed on Linux using [MKL](https://software.intel.com/tools/onemkl) and [GSL](http://www.gnu.org/software/gsl/) library, and is run from the command line. 

If user does not enter any parameters, or if user enters invalid parameters, `fastQTLmapping` will print a help file on the screen.



## Parameter setting

| Interface                             | Description                                                  |
| ------------------------------------- | ------------------------------------------------------------ |
| `--omics1 <omics1FileName> ['bfile']` | The first omics data file path. If there is a parameter `bfile`, it means that the first omics data is genomic data in [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed); in this case only the file path without extension is given. |
| `--omics2 <omics2FileName>`           | The second omics data file path.                             |
| `--out <outputFileName>`              | Output file path.                                            |
| `[-p <globalP>`]                      | Global significance threshold. By default the significance threshold is corrected by Bonferroni correction. |
| `[--cov <covarFileName>]`             | Covariate file path. By default there is no covariate correction. |
| `[--na <NASign>]`                     | A string identifying the missing value. The default is `NA`. |
| `[--missing-rate <missingRateThd>]`   | Missing rate threshold of concatenated dependent and explanatory variables. The default value is 10%. Loci-pairs above this threshold will be filtered out. |
| `[--dl <distLv> ...]`                 | A series of increasing numbers used to divide loci-pairs into different distance levels based on physical distance. Assuming the value is set to $d_1, d_2, d_3 ... d_m$​​​​​, then the loci-pairs with physical distance belongs to $[0, d1)$​​​​ is considered as level 1， $[d_2, d_3)$​​​ is level 2, ..., $[d_m, +\infty)$​​​​ is level (m+1)。If not set, all results are categorized as level 1. |
| `[--dlp <distLvP> ...]`               | A series of numbers indicate the significance thresholds for the 1st to mth distance levels, which need to be the same length as `distLv`. The significance threshold for the (m+1)st level is `globalP`. If not set, all distance levels are set to the same significant threshold of `globalP`. |
| `[--threads <threadMaxN>]`            | Maximum number of parallelism. The default value is 1.       |
| `[--omics1norm zscore\|rank]`         | Set the normalization method for the first omics data. `zscore` is Z-value normalization. `rank` is the rank-base normalization. The default is no normalization. |
| `[--omics2norm zscore\|rank]`         | Set the normalization method for the second omics data. `zscore` is Z-value normalization. `rank` is the rank-base normalization. The default is no normalization. |

 

##   Input

Omics data is a space- or tab-delimited matrix without table headers. Each row represents a molecular loci, the first three columns record the locus name, the chromosome  and the bp, and subsequent columns record the molecular-level in each individuals. In particular, when the first omics data represents genomic, the [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed) can be used as input data.

Covariate data is a space- or tab-delimited matrix without table headers. Each column represents a sample; each row represents a covariate.

The individual order of the omics data and the covariate data must be consistent. 

Missing values need to be marked by a unique string. 



##   Output

If the program runs normally, two files are output: the file `outputFileName.log` for the log file and the file `outputFileName` for the result file.



## Operating Environment

Software environment: x64 processor-based, 64-bit operating system, Linux v3.10.0 and above.

Programming language: C++.

Library: MKL 2019.0 or above , GSL v2.6 or above, OpenMP.
