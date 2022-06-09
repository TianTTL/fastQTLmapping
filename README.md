# fastQTLmapping

`FastQTLmapping` is a computationally efficient, exact, and generic solver for exhaustive multiple regression analysis involving extraordinarily large numbers of dependent and explanatory variables with covariates, which is particularly helpful in QTL-like analysis. `FastQTLmapping` can afford omics data containing tens of thousands of individuals and billions of molecular locus.



## Operating Environment

Software environment: x64 processor-based, 64-bit operating system, Linux v3.10.0 and above.

Programming language: C++.

Library: MKL 2019.0 or above 

              GSL v2.6 or above

              OpenMP.



## Installation

1. Download the source code and unzip it.

2. Enter the source code directory.

3. Edit `Makefile` to specify the GSL and MKL paths.

4. Type `make` to install fastQTLmapping.
   
   After successful installation, an executable file `fastQTLmapping` will be generated in the current path.

## Running fastQTLmapping



FastQTLmapping exhaustively test correlations between quantitative molecular-level  of locus from two omics, respectively, by using linear regression models.

Taking eQTL analysis as an example, for each gene-SNP pair, the association between gene expression $y$ and genotype $x$ is assumed to be linear under covariable matrix $C$: 

$y=\alpha + \beta x + \gamma C + \epsilon, \epsilon \sim i.i.d. N(0, \sigma ^2)$

FastQTLmapping accepts input files in text format and in Plink binary format. The output file is in text format and contains test statistics of regression test across all locus-pairs, with the ability to control the volume of the output at preset significance thresholds.

```bash
fastQTLmapping --omics1 <omics1FileName> [bfile] --omics2 <omics2FileName>
                   --out <outputFileName> [--outPcs <outPcs>] [-p <globalP>]
                   [--ploose <PLooseMarg>] [--cov <covarFileName>] [--categ
                   <categFlag>...] [--na <NASign>] [--missing-rate
                   <missingRateThd>] [--dl <distLv>...] [--dlp <distLvP>...]
                   [--threads <threadMaxN>] [--omics1norm (zscore|rank)]
                   [--omics2norm (zscore|rank)] [--rpl <rplFileName>]
Name>]
```

**Parameter list**

| Interface                             | Description                                                                                                                                                                                                                                                |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--omics1 <omics1FileName> ['bfile']` | The first omics data file path. If there is a parameter `bfile`, it means that the first omics data is genomic data in [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed); in this case only the file path without extension is given. |
| `--omics2 <omics2FileName>`           | The second omics data file path.                                                                                                                                                                                                                           |
| `--out <outputFileName>`              | Output file path.                                                                                                                                                                                                                                          |
| `[--outPcs <outPcs>]`                 | The number of significant digits.                                                                                                                                                                                                                          |
| `[-p <globalP>`]                      | Global significance threshold. By default the significance threshold is corrected by Bonferroni correction.                                                                                                                                                |
| `[--cov <covarFileName>]`             | Covariate file path. By default there is no covariate correction.                                                                                                                                                                                          |
| `[--categ <categFlag>]`               | A series of integer indicate which rows of covariates are categorical, starting from 1.                                                                                                                                                                    |
| `[--na <NASign>]`                     | A string identifying the missing value. The default is `NA`.                                                                                                                                                                                               |
| `[--missing-rate <missingRateThd>]`   | Missing rate threshold of concatenated dependent and explanatory variables. The default value is 10%. Locus-pairs above this threshold will be filtered out.                                                                                               |
| `[--dl <distLv> ...]`                 | A series of increasing numbers used to divide locus-pairs into different distance levels based on physical distance.                                                                                                                                       |
| `[--dlp <distLvP> ...]`               | A series of numbers indicate the significance thresholds for the 1st to mth distance levels, which need to be the same length as `distLv`.                                                                                                                 |
| `[--threads <threadMaxN>]`            | Maximum number of parallelism. The default value is 1.                                                                                                                                                                                                     |
| `[--omics1norm zscore\|rank]`         | Set the normalization method for the first omics data. `zscore` is Z-value normalization. `rank` is the rank-base normalization. The default is no normalization.                                                                                          |
| `[--omics2norm zscore\|rank]`         | Set the normalization method for the second omics data. `zscore` is Z-value normalization. `rank` is the rank-base normalization. The default is no normalization.                                                                                         |
| `[--rpl  <rplFileName>]`              | Set the path of replication list file.                                                                                                                                                                                                                     |

## Input

### omics data

Omics data is a space- or tab-delimited matrix without table headers. Each row represents a molecular loci, the first three columns record the locus ID, the chromosome and the bp, and subsequent columns record the quantitative molecular-level in each individuals.  Missing values need to be marked by a unique string.

In particular, when the first omics data represents genomic, the [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed) can be used as input data.

### covariate data

It is common to include covariates in an eQTL model to account for such effects as population stratiﬁcation, gender, age, white blood count and other clinical variables.

Covariate data is a space- or tab-delimited matrix with no table headers recorded in `covarFileName` file. Each column represents a sample; each row represents a covariate. User can specify which columns of covariates are categorical variables by setting `--categ`. Numeric covariates can be represented as numbers, while categorical covariates can be represented  as numbers or strings. 

The individual order of the omics data and the covariate data must be consistent. Missing values need to be marked by a unique string. 

## Output

If the fastQTLmapping runs normally, two files are output: the `outputFileName.log` file for the log file and the `outputFileName` file for the result file.

Each line of the `outputFileName` file records the statistics of regression analysis, including the IDs of the locus(`omics1` `omics2`), estimated coefficients(`BETA`), standard errors(`SE`), t-statistics(`T`), significance levels(`P`), degrees of freedom(`NMISS`), and distance level(`distance_level`).

## List of features

**Manuscript**

If the user does not enter any parameters, or if the user enters invalid parameters, `fastQTLmapping` will print the help file on the screen.

**Normalization**

To make the method more robust to outliers in omics data, fastQTLmapping has an option that allows to normalize the quantitative molecular-level prior to any analysis.  

FastQTLmapping provides two normalization methods, Z-value normalization(`zscore`) and rank-based normalization(`rank`). The user can specify the normalization method for each omics separately by parameters `--omics1norm` and `--omics2norm`.

**Covariates correction**

FastQTLmapping reduces the covariates of the regression model by orthogonalizing $x$ and $y$ respect to $C$ using Gram-Schmidt orthogonalization.([Longley J W, Longley R D, 1997](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199707/08)4:4%3C295::AID-NLA102%3E3.0.CO;2-D)).

Categorical covariates are converted into dummy variables in the data pre-processing. Dummy Variables act as indicators of the presence or absence of a category in a categorical variable: 0 represents absence while 1 represents presence.  The conversion of categorical variables into dummy variables leads to the formation of the binary matrix where each row represents a particular category.

**Exact results**

FastQTLmapping conducts QTL mapping in two steps. In the first step, fastQTLmapping imputes all missing values by mean values. Then, under a relaxed threshold (`--outPcs`, by default 100 folds relaxed than the study-wide significance threshold), fastQTLmapping identifies all candidate meQTLs. In the second step, fastQTLmapping conducts standard multiple regressions using the original data without imputation to get exact results. In this way, fastQTLmapping guarantees exact results.

**Replacation mode**

If user sets `--rpl`, only locus-pairs in the `rplFileName` will be test.

The `rplFileName` file contains two columns of text separated by spaces, each line represents a locus-pair, and the two columns represent the locus IDs in the first and second omics, respectively.

**Distance level**

Different thresholds can be specified according to physical distances between the markers under investigation, which facilitates the analysis of cis- and trans-mQTLs.

The user can divide the locus-pairs into a series of distance levels by setting the parameter `--dl`. Assuming the value is set to $d_1, d_2, d_3 ... d_m​$​​​​, then the locus-pairs with physical distance belongs to $[0, d1)$​​​​ is considered as level 1，$[d_2, d_3)$​​​ is level 2, ..., $[d_m, +\infty)​$​​​ is level (m+1)。If not set, all results are categorized as level 1.

The user can specify the significance threshold for each distance level by setting the parameter `--dlp`. The significance threshold for the (m+1)st level is `globalP`. If not set, all distance levels are set to the same significant threshold of `globalP`.

## example

```bash
fastQTLmapping \
 --omics1 test.omics1.data \
 --omics2 test.omics2.data \
 --out test.rlt \
 --cov test.covar.data \
 --na NA --missing-rate 0.1 \
 --dl 1000000 2000000 --dlp 1 0.9 -p 0.8 \ 
 --threads 1 --omics1norm zscore --omics2norm rank
```

## Possible upcoming features

* Counting total number of regression test and calculating Bonferroni threshold for each distance level.

* FDR procedure described by Storey and Tibshirani (ST) ([Storey and Tibshirani, 2003](https://www.pnas.org/doi/abs/10.1073/pnas.1530509100)).

* Batch and online  input and processing.

* Permutation schedule for each molecular loci.

* Estimating independent test number by MLE ([Galwey N W, 2009](https://pubmed.ncbi.nlm.nih.gov/19217024/)).
