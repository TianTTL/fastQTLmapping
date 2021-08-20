# xQTLmapping
 fast xQTL mapping tools



##   命令行参数

| 序号 | 程序内变量名   | 描述                                   |
| ---- | -------------- | -------------------------------------- |
| 1    | omics1FileName | 基因组数据文件名，不包括扩展名，字符串 |
| 2    | omics2FileName | trait数据文件名，包括扩展名，字符串    |
| 3    | outputFileName | 输出文件文件名，字符串                 |
| 4    | NASign         | 缺失值标记，字符串                     |
| 5    | missingRateThd | missing rate 阈值，浮点型              |
| 6    | MAFThd         | MAF 阈值，浮点型                       |
| 7    | cisDist        | cis- / trans- 划分阈值                 |
| 8    | cisP           | cis- 显著性阈值，浮点型                |
| 9    | transP         | trans- 显著性阈值，浮点型              |
| 10   | thread_cout    | 并行数                                 |

 

##   输入描述

基因组数据的格式为plink软件生成的二进制基因组文件，由三个扩展名分别为.bed .bim .fam的文件共同构成。

trait数据指待鉴定QTL的任意生命组学，格式为以空格或制表符分隔的文本型二维矩阵。前三列分别表示位点名称，位点所在染色体，位点位置，之后的每列表示一个样本；每行表示一个组学位点。

基因组数据和trait数据的样本顺序需保持一致。

基因组数据的缺失值遵循plink格式的默认定义，triat文件的缺失值用唯一字符串标识，该标识需通过第四个命令行参数指明。

显著性阈值都是(0,1]之间的浮点数，可以用科学计数法表示。



##   输出描述

如果程序正常运行，则会生成三个输出文件。设用户在命令行指定的输出文件路径为`$output`，则文件`$output`为日志文件，文件`$output.cis.rlt`为 cis- 结果文件，文件`$output.trans.rlt`为 trans- 结果文件。结果包括线性回归分析中的β, se, r2, t, -lg(P)。



## 运行环境

硬件环境：处理器24核或以上；内存(RAM)：128 GB或以上；硬盘256GB或以上。

软件环境：基于x64处理器，64位操作系统，Linux v3.10.0及以上。

编程语言：C++。

库：MKL 2019.0或以上 ,GSL v2.6或以上，OpenMP。
