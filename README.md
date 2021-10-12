# xQTLmapping
`xQTLmapping`能够在组学之间鉴定位点间相关性，同时表现出极高的运行效率和极小的内存开销。本软件能够支持的最大组学规模为万个样本及十亿个位点。

本软件也支持协变量校正，多线程加速，缺失值填补，数据标准化，自适应显著性阈值，可定制结果格式等功能。



## 使用方法

本软件基于Linux开发，使用命令 `fastQTLmapping`调用软件，通过命令行参数传入信息。

在不输入任何参数，或输入参数有误时，`fastQTLmapping`将会在屏幕上打印帮助文档。



## 命令行参数

| Interface                             | 程序内变量名                  | 描述                                                         |
| ------------------------------------- | ----------------------------- | ------------------------------------------------------------ |
| `--omics1 <omics1FileName> ['bfile']` | `omics1FileName`, `bfileFlag` | 第一个组学数据文件名，字符串。如果有参数`bfileFlag`，则说明第一个组学为plink binary format格式的基因组数据，此时只给出不带扩展名的文件路径。 |
| `--omics2 <omics2FileName>`           | `omics2FileName`              | 第二个组学数据文件名，字符串。                               |
| `--out <outputFileName>`              | `outputFileName`              | 输出文件文件名，字符串。                                     |
| `[-p <globalP>`]                      | `globalP`                     | 全局显著性阈值。默认为使用 Bonferroni correction 校正后的显著性阈值。 |
| `[--cov <covarFileName>]`             | `covarFileName`               | 协变量数据文件名，字符串。如不设置，则表示计算中不涉及协变量。 |
| `[--na <NASign>]`                     | `NASign`                      | 缺失值标记，字符串。默认值为 `NA`。                          |
| `[--missing-rate <missingRateThd>]`   | `missingRateThd`              | missing rate 阈值，浮点型。默认值为 10%。低于此阈值的位点将被剔除。 |
| `[--dl <distLv> ...]`                 | `distLv`                      | 用于将位点对之间根据物理距离划分为不同的水平。可以设置依次递增的多个值，如$d_1, d_2, d_3 ... d_m$​​​​​，则物理位置距离属于 $[0, d1)$​​​​ 为 level 1， $[d_2, d_3)$​​​​ 为 level 2, $[d_m, +\infty)$​​​​ 为level m+1。如不设置，则所有结果均归为 level 1。 |
| `[--dlp <distLvP> ...]`               | `distLvP`                     | 对于不同的距离水平分别设置的显著性阈值。数量等于`distLv`，分别表示第1到第m个level 的显著性阈值。第m+1个level的显著性阈值为`globalP`。如不设置，则所有距离水品的显著性阈值均为``globalP`。 |
| `[--threads <threadMaxN>]`            | `threadMaxN`                  | 最大并行数。默认值为1。                                      |
| `[--omics1norm zscore|rank]`          | `omics1Norm`                  | 设置第一个组学数据的标准化方式。`zscore`为Z值标准化。`rank`为rank-base标准化。不设置则不做标准化。 |
| `[--omics2norm zscore|rank]`          | `omics2Norm`                  | 设置第二个组学数据的标准化方式。`zscore`为Z值标准化。`rank`为rank-base标准化。不设置则不做标准化。 |

 

##   输入描述

组学数据指参与相关性检验的两个生命组学，格式为以空格或制表符分隔的文本型二维矩阵，无表头。前三列分别表示位点名称，位点所在染色体，位点位置，之后的每列表示一个样本；每行表示一个组学位点。

特别地，当第一个组学表示基因组时，其格式也可以为plink软件生成的二进制基因组文件，由三个扩展名分别为.bed .bim .fam的文件共同构成。

协变量数据指需要被校正的协变量信息，格式为以空格或制表符分隔的数值型二维矩阵，无表头。每列表示一个样本；每行表示一个协变量。

组学数据和协变量数据的样本顺序需保持一致。缺失值需用唯一字符串标识。plink二进制基因组数据的缺失值遵循plink软件的默认定义，



##   输出描述

如果程序正常运行，则会输出两个文件。设用户在命令行指定的输出文件路径为`$outputFileName`，则文件`$outputFileName.log`为日志文件，文件`$outputFileName`为结果文件。



## 运行环境



## 运行环境

硬件环境：处理器24核或以上；内存(RAM)：128 GB或以上；硬盘256GB或以上。

软件环境：基于x64处理器，64位操作系统，Linux v3.10.0及以上。

编程语言：C++。

库：MKL 2019.0或以上 ,

​		GSL v2.6或以上，

​		OpenMP。
