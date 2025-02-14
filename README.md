# fastQTLmapping

FastQTLmapping is a computationally efficient, exact, and generic solver for exhaustive multiple regression analysis involving extraordinarily large numbers of dependent and explanatory variables with covariates, which is particularly helpful in QTL-like analysis. FastQTLmapping can afford omics data containing tens of thousands of individuals and billions of molecular loci.

## Current Version

0.9.9

## Counting Mode

In counting mode, fastQTLmapping will count the number of loci-pairs in each distance level and calculate the significant threshold that controls FWER using bonferroni method. The distance level threshold and FWER(Family-Wise Error Rate) are user specified.

Counting mode is computationally fast and is usually performed before discovery mode.  This process can help users evaluate computational scale and estimate significance thresholds for studies.

## Discovery Mode

In discovery mode, fastQTLmapping exhaustively test correlations between quantitative molecular-level  of loci from two omics, respectively, by using linear regression models.

Taking eQTL analysis as an example, for each gene-SNP pair, the association between gene expression $y$ and genotype $x$ is assumed to be linear under covariable matrix $C$: 

$y=\alpha + \beta x + \gamma C + \epsilon, \epsilon \sim i.i.d. N(0, \sigma ^2)$

FastQTLmapping accepts input files in text format and in Plink binary format. The output file is in text format and contains test statistics of regression test across all loci-pairs, with the ability to control the volume of the output at preset significance thresholds.

## Operating Environment

Software environment: x64 processor-based, 64-bit operating system, Linux v3.10.0 and above.

Programming language: C++

Library: MKL (>= 2019.0)

              GSL (>= v2.6)

              OpenMP

## Installation

1. Download the source code and unzip it.

2. Enter the source code directory.

3. Edit `Makefile` to specify the GSL and MKL paths.

4. Type `make` to install fastQTLmapping.
   
   After successful installation, an executable file `fastQTLmapping` will be generated in the current path.

## Running fastQTLmapping

```bash
fastQTLmapping --omics1 <omics1FileName> [bfile] --omics2 <omics2FileName>
               [bfile] --out <outputFileName> [--outPcs <outPcs>] [--threads
               <threadMaxN>] [-h] count [--dl <distLv>...] [--FWER <FWER>]

fastQTLmapping --omics1 <omics1FileName> [bfile] --omics2 <omics2FileName>
               [bfile] --out <outputFileName> [--outPcs <outPcs>] [--threads
               <threadMaxN>] [-h] discovery [-p <globalP>] [--ploose
               <PLooseMarg>] [--cov <covarFileName>] [--categ
               <categFlag>...] [--na <NASign>] [--MR <msRtThd>] [--SD
               <sdThd>] [--dl <distLv>...] [--dlp <distLvP>...]
               [--omics1norm (zscore|rank)] [--omics2norm (zscore|rank)]
               [--chunk <chunkSize>]
```

**Mode-independent Parameter List**

| Interface                           | Description                                                  |
| ----------------------------------- | ------------------------------------------------------------ |
| `--omics1 <omics1FileName> [bfile]` | The first omics data file path. The `bfile` flag indicates that the input file is genomic data in [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed) format. |
| `--omics2 <omics2FileName> [bfile]` | The second omics data file path. The `bfile` flag indicates that the input file is genomic data in [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed) format. |
| `--out <outputFileName>`            | Output file path.                                            |
| `[--outPcs <outPcs>]`               | The number of significant digits. The default value is 4.    |
| `[--threads <threadMaxN>]`          | Maximum number of parallelism. The default value is 1.       |
| `[--chunk <chunkSize>]`             | The size of chunks that divides the computational load. Users can reduce the memory usage by reducing this parameter. The default value is 5000. |
| `[-h], [--help]`                    | Show manuscript.                                             |

**Counting Mode Parameter List**

| Interface            | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| `count`              | Counting mode command.                                       |
| `[--dl <distLv>...]` | A series of increasing numbers used to divide loci-pairs into different distance levels based on physical distance. The default setting is NULL. |
| `[--FWER <FWER>]`    | Family-wise error rate that expected by user. The default value is 0.05. |

**Discovery Mode Parameter List**

| Interface                     | Description                                                  |
| ----------------------------- | ------------------------------------------------------------ |
| `discovery`                   | Discovery mode command.                                      |
| `[-p <globalP>`]              | Global significance threshold. The default value is 1, which means keep all results. Note that the results of QTL mapping study usually result in huge volume, and this parameter can be used to reduce the size of the results significantly. |
| `[--ploose <PLooseMarg>]`     | Margin for calculating loosed significance thresholds. The default value is 100. |
| `[--cov <covarFileName>]`     | Covariate file path. The default setting is NULL.            |
| `[--categ <categFlag>]`       | A series of integer indicate which rows of covariates are categorical, starting from 1. |
| `[--na <NASign>]`             | A string identifying the missing value. The default is `NA`. |
| `[--MR <msRtThd>]`            | Missing rate threshold of concatenated dependent and explanatory variables. Loci-pairs above this threshold will be filtered out. The default value is 0.1. |
| `[--SD <sdThd>]`              | Standard deviation threshold of quantitative loci. Loci with a standard deviation below `sdThd` will be considered constant and filtered out. The default value is 1e-6. |
| `[--dl <distLv> ...]`         | A series of increasing numbers used to divide loci-pairs into different distance levels based on physical distance.  The default setting is NULL. |
| `[--dlp <distLvP> ...]`       | A series of numbers indicate the significance thresholds for the 1*st* to m*th* distance levels, which need to be the same length as `distLv`.  The default setting is NULL. |
| `[--omics1norm zscore\|rank]` | Set the normalization method for the first omics data. `zscore` is Z-value normalization. `rank` is rank-base normalization. If not set, no normalization will be done. The default is no normalization. |
| `[--omics2norm zscore\|rank]` | Set the normalization method for the second omics data. `zscore` is Z-value normalization. `rank` is rank-base normalization. If not set, no normalization will be done. The default is no normalization. |

## Input

### Omics Data

FastQTLmapping accepts two widely used data formats.

For quantitative omics data, input file is a space- or tab-delimited matrix without table headers. Each row represents a molecular loci. The first four columns record the loci ID, the chromosome, the start position and the end position, and subsequent columns record the quantitative molecular-level in each individuals.  Missing values need to be marked by a unique string that specified by `NASign`.

For genomics data in [Plink Binary File](http://www.cog-genomics.org/plink/1.9/formats#bed) format, use `bfile` flag to make fastQTLmapping input binary file set: `omics1FileName.bed`+`omics1FileName.bim`+`omics1FileName.fam`.

### Covariate Data

It is common to include covariates in an eQTL model to account for such effects as population stratiﬁcation, gender, age, white blood count and other clinical variables.

Covariate data is a space- or tab-delimited matrix with no table headers recorded in `covarFileName` file. Each column represents a sample; each row represents a covariate. User can specify which rows of covariates are categorical variables by setting `--categ`. Numeric covariates can be represented as numbers, while categorical covariates can be represented  as numbers or strings. 

The individual order of the omics data and the covariate data must be consistent. Missing values need to be marked by a unique string. 

## Output

If the fastQTLmapping runs normally, two files are output: the `outputFileName.log` file for the log file, the `outputFileName.cnt` file for counting mode results and the `outputFileName` file for discovery mode results.

Each line of the `outputFileName` file records the statistics of regression analysis, including the IDs of the loci(`omics1` `omics2`), number of samples without missing (`NMISS`), distance level(`distance_level`), estimated coefficients(`BETA`), standard errors(`SE`), t-statistics(`T`), significance levels(`P-value`) and Q-value(`Q-value`). 

The output file retains `outPcs` significant digits.

## Example

```bash
# counting mode
fastQTLmapping \
 --omics1 testdata/test.omics1.data \
 --omics2 testdata/test.omics2.data \
 --out testdata/test.rlt \
 count \
 --dl 1000000 2000000 --FWER 0.01

# discovery mode
fastQTLmapping \
 --omics1 testdata/test.omics1.data \
 --omics2 testdata/test.omics2.data \
 --out testdata/test.rlt \
 --threads 1 \
 discovery \
 --cov testdata/test.covar.data --categ 1 \
 --na NA --MR 0.1 \
 --dl 1000000 2000000 --dlp 1 0.9 -p 0.8 \
 --omics1norm zscore --omics2norm rank
```

## List of Features

**Manuscript**

If the user enter `-h, --help` parameters or invalid parameters, `fastQTLmapping` will print a brief help file on the screen.

**Normalization**

To make the method more robust to outliers in omics data, fastQTLmapping has an option that allows to normalize the quantitative molecular-level prior to any analysis.  

FastQTLmapping provides two normalization methods, Z-value normalization(`zscore`) and rank-based normalization(`rank`). The user can specify the normalization method for each omics separately by parameters `--omics1norm` and `--omics2norm`.

**Covariates Adjustment**

FastQTLmapping reduces the covariates of the regression model by orthogonalizing $x$ and $y$ respect to $C$ using Gram-Schmidt orthogonalization.([Longley J W, Longley R D, 1997](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199707/08)4:4%3C295::AID-NLA102%3E3.0.CO;2-D)).

Categorical covariates are converted into dummy variables in the data pre-processing. Dummy Variables act as indicators of the presence or absence of a category in a categorical variable: 0 represents absence while 1 represents presence.  The conversion of categorical variables into dummy variables leads to the formation of the binary matrix where each row represents a particular category.

When processing covariates, FastQTLmapping will perform QR decomposition on the covariate matrix. This process requires that the number of covariates (including dummy variables) does not exceed the number of samples; otherwise, an error will occur. Generally, users do not need to be concerned about such errors, but they may be triggered when there are a large number of categories for categorical covariates.

**Exact Results**

FastQTLmapping conducts QTL mapping in two steps. In the first step, fastQTLmapping imputes all missing values by mean values. Then, under a relaxed threshold (`PLooseMarg`, by default 100 folds relaxed than the study-wide significance threshold), fastQTLmapping identifies all candidate QTLs. In the second step, fastQTLmapping conducts standard multiple regressions using the original data without imputation to get exact results. In this way, fastQTLmapping guarantees exact results.

**Distance Level**

Different thresholds can be specified according to physical distances between the markers under investigation, which facilitates the analysis of cis- and trans-mQTLs.

The user can divide the loci-pairs into a series of distance levels by setting the parameter `--dl`. Assuming the value is set to $d_1, d_2, d_3 ... d_m$​​​​, then the loci-pairs with physical distance belongs to $[0, d1)$​​​​ is considered as level 1，$[d_2, d_3)$​​​ is level 2, ..., $[d_m, +\infty)$​​​ is level (m+1)。If not set, all results are categorized as level 1.

The user can specify the significance threshold for each distance level by setting the parameter `--dlp`. The significance threshold for the (m+1)*th* level is `globalP`. If not set, all distance levels are set to the same significant threshold of `globalP`.

**FDR Procedure**

FastQTLmapping calculates the q-values for all of p-values using FDR method([Storey and Tibshirani, 2003](https://doi.org/10.1073/pnas.1530509100)). A faster solution of the pi0 estimation is used by fastQTLmapping([Storey J D, Taylor J E, Siegmund D. Strong control, 2004](https://doi.org/10.1111/j.1467-9868.2004.00439.x)).

**Peak Memory Estimation**

In the preprocessing stage, fastQTLmapping will estimate memory requirements of current job and print it on the screen. Noted that when the result is massively inflated, the actual memory consumption will be larger than the estimated value.

## Possible Upcoming Features

* Replication mode
* Permutation schedule for each molecular loci.
* Estimating independent test number by MLE ([Galwey N W, 2009](https://pubmed.ncbi.nlm.nih.gov/19217024/)).

