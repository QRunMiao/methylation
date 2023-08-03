# 目的 
学习基于bismark和DSS软件包的甲基化数据分析流程 #Bismark可高效率的分析亚硫酸氢盐序列(BS-Seq)数据。Bismark用参考基因组进行比对，同时进行胞嘧啶甲基化提取。Bismark是用Perl编写的，在命令行运行。用亚硫酸氢盐处理的reads用短Reads对准软件Bowtie 2或HISAT2进行比对。

熟悉Linux服务器工作空间和基本的生物信息数据分析协议
# 测序数据

Open browser, vist EBI-ENA search page and search for the GEO accession. On the result page, we could get the expected FTP address. Use aria2 to download the data and we could get .fastq.gz files.
```
cd /mnt/c/shengxin
mkdir -p NBT_repeat/data/seq_data/WT_mESC_rep1
mkdir -p NBT_repeat/data/seq_data/TetTKO_mESC_rep1
cd NBT_repeat/data/seq_data
aria2c -d ./WT_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/001/SRR7368841/SRR7368841.fastq.gz
aria2c -d ./TetTKO_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/005/SRR7368845/SRR7368845.fastq.gz
#aria2c c--continue[=true|false]，使用此选项可恢复下载，由网络浏览器或其他程序启动，按顺序从开始。目前只有这个选项适用于http(s)/ftp下载
-d 存储下载文件的目录
```
# 参考基因数据
Open browser, visit ensembl download page, and find the DNA .fasta file. On the result page, we could get the expected FTP address.
```
cd /mnt/c/shengxin
mkdir -p NBT_repeat/data/genome_data
cd NBT_repeat/data/genome_data

# Basic command line download
#Aria2 是一款多线程下载工具，可以同时下载多个文件 -d指定下载文件的目录，-d ./：指定下载目录为当前目录（./）；-Z：启用压缩支持。aria2c 将自动解压下载的文件。
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
```
和直接从ensembl进去选择大鼠是一样的

# 质控和修剪
```
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/

# 质控，查看原始seq数据的质量
fastqc --t 3 ./WT_mESC_rep1/*.fastq.gz ./TetTKO_mESC_rep1/*.fastq.gz
#fastqc 是一款用于分析原始测序数据质量的工具。它可以生成关于测序数据质量、GC含量、序列重复性等方面的统计报告。-t 3：指定使用 3 个线程并发处理分析任务。通配符 *.fastq.gz 表示匹配

# 质量和适配器修剪，随后fastqc操作；即依次对每个匹配到的 fastq.gz 文件进行质量修剪和过滤，并将处理后的结果保存到指定的输出目录中。同时，根据 --fastqc 选项，trim_galore 还会为每个处理后的文件生成相应的质量控制报告
trim_galore -o ./WT_mESC_rep1/trimmed_data/ --fastqc ./WT_mESC_rep1/*.fastq.gz
#trim_galore 是一个用于自动化质量修剪和去除低质量序列的工具。它可以检测和去除低质量的碱基、读长不匹配、接头污染等。

# Trim Galore used options:
-o/--output_dir <DIR>: If specified all output will be written to this directory instead of the current directory. If the directory doesn't exist it will be created for you.
--fastqc: Run FastQC in the default mode on the FastQ file once trimming is complete.
```
打开fastqc.html文件，可以看到C含量几乎没有，本样本测序测定了处理后的甲基化水平，C已经被变成U/T，只有5hmC仍被识别为C，所以是正常的
# 甲基化分析

在质量控制和适配器修剪后，可以按照基于Bismark的BS-seq数据分析协议对ACE-seq数据进行甲基化分析。

用法参考网址 

https://www.jianshu.com/p/877a5716c24a
https://www.zxzyl.com/archives/759/

bismark的分析分为三步：

基因组索引构建：需要将基因组进行相应的转换，模拟bisulfite处理后基因组，这样才能用于比对。在每个新物种比对之前，都需要进行该处理。bismark内部默认调用bowite软件进行比对。

序列比对：输入为下机序列、基因组和比对参数，bismark会输出比对结果和甲基化检测的结果，默认格式为BAM，同时还会有一个统计报告。

甲基化位点提取：这一步是可选的，会利用bismark的比对结果来获得甲基化信息。这一步会把甲基化分为不同的类别（CG，CHG和CHH），获得链特异性信息，并且提供过滤参数；也可以对甲基化信息进行调整，或者进行更多的深入分析。

# 1. 基因组索引

在进行比对之前，需要对感兴趣的基因组进行亚硫酸盐转换和索引。基于Bismark，我们可以同时使用Bowtie和Bowtie2来完成索引工作。

#下载bismark

通过官网下载 https://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark

或者brew 安装

```brew install bismark```

```
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/
#使用bismark_genome_preparation来对基因组进行处理，输入是给定的目录，目录里面是.fa或者.fasta文件。
命令：bismark_genome_preparation [options] <path_to_genome_folder> #genome_folder 基因组文件夹
# bismark_genome_preparation --path_to_bowtie这里使用bowtie2 <path_to_genome_folder>   索引路径 
结果：在基因组目录下产生Bisulfite_Genome目录

bismark_genome_preparation --bowtie2 /mnt/c/shengxin/NBT_repeat/data/genome_data/

```
```bismark_genome_preparation``` used options:

```--bowtie2```: This will create bisulfite indexes for Bowtie 2. (Default: ON).

# 2. 比对
甲基化数据分析程序的核心是将测序读数与参考基因组对齐，并且假设所有数据都已经过质量和适配器修剪。对齐报告可以在子目录bismark-align -reports/中找到。

比对的过程需要提供两个信息：
输入目录：已经建立好的索引目录，即基因组索引步骤中的索引目录
输入文件：测序reads，fastq格式

Bismark支持：

reads可以是序列格式：fastq或者fasta，也可以是压缩文件.gz

命令：USAGE: bismark [options] <genome_folder> {-1 <mates1> -2<mates2> | <singles>}
```
genome_path="/mnt/c/shengxin/NBT_repeat/data/genome_data/"
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/

# read alignment
bismark -o ./WT_mESC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./WT_mESC_rep1/SRR7368841.fastq.gz
bismark -o ./TetTKO_mESC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./TetTKO_mESC_rep1/SRR7368845.fastq.gz

结果：产生两个文件

①文件名_bismark.sam（序列比对结果的详细信息）

②文件名_bismark_mapping_report.txt（比对结果的总结）

--genome_folder: Arguments, The path to the folder containing the unmodified reference genome as well as the subfolders created by the bismark_genome_preparation script.
samtools used command and options:

cat: Command concatenate BAMs
-o: output file format
```

# Aligned reads deduplication-去重

哺乳动物的基因组是如此之大，以至于不太可能遇到几个真正独立的片段，它们排列在同一个基因组位置。这种解读更有可能是PCR扩增的结果。因此，对于大型基因组来说，去除重复读取是一条有效的途径。

重复数据删除步骤可以通过deduplicate_bismark脚本完成，重复数据删除报告可以在子目录deduplicate_reports /中找到

```
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/
mkdir -p ./WT_mESC_rep1/deduplicated_result/
mkdir -p ./TetTKO_mESC_rep1/deduplicated_result/

# aligned reads deduplication
deduplicate_bismark --bam [options] <filenames> 
--bam: The output will be written out in BAM format instead of the default SAM format.
--output_dir 输出目录 --bam  比对后的BAM输出文件；bam与sam是同一文件的不同压缩方式，sam占据的内存更小，压缩更高，二进制形式

deduplicate_bismark --bam --output_dir ./WT_mESC_rep1/deduplicated_result/ ./WT_mESC_rep1/bismark_result/SRX4241790_trimmed_bismark_bt2.bam

deduplicate_bismark --bam --output_dir ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/bismark_result/*.bam
```

# Methylation information extracting-甲基化信息提取
Bismark包中的bimark_methylation_extractor脚本可以从比对结果文件中提取甲基化信息，并充当Bismark包的端点。此外，甲基化信息可以很容易地转换成其他格式，便于脚本中的下游分析。

bismark_methylation_extractor
是一个非常重要的输出文件为全基因组胞嘧啶报告文件（*CX_report.txt），该文件是一个非常核心的记录胞嘧啶深度信息的文件，可进行甲基化水平的计算、区域甲基化信号文件的转化、差异甲基化水平的计算等。
使用bismark_methylation_extractor来提取每个C位点的甲基化的情况，可以使用--bedGraph，来产生bedgraph和coverage的文件；--count可以产生count文件，输出每个CpG位点的深度文件。

bismark_methylation_extractor用法：

    --pair-end         指定双端数据
    --comprehensive    输出CHG CHH CpG的甲基化信息
    --no-overlap       去除reads重叠区域的bias
    --bedGraph         输出bedGraph文件
    --counts           每个C上甲基化reads和非甲基化reads的数目
    --report           一个甲基化summay
    --cytosine_report  输出全基因组所有CpG
    --genome_folder    <path_to_reference_genome>
    input.bam          输入文件
    -o >output_dir>    输出路径
```
genome_path="/mnt/c/shengxin/NBT_repeat/data/genome_data/"
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/

# methylation information extracting
bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \#生信文件格式 | BedGraph（基因组浏览器绘制），甲基化提取的输出可以使用选项--bedGraph转换为一个bedGraph和coverage文件
--cytosine_report --genome_folder ${genome_path} \
-o ./WT_mESC_rep1/deduplicated_result/ ./WT_mESC_rep1/deduplicated_result/*.bam

bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \
--cytosine_report --genome_folder ${genome_path} \
-o ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/deduplicated_result/*.bam

--gzip: The methylation extractor files (CpG_OT_..., CpG_OB_... etc) will be written out in a GZIP compressed form to save disk space. This option is also passed on to the genome-wide cytosine report.
--cytosine_report指报道全基因组所有的CpG。只有当指定--cytosine_report时才需要genome_folder
--genome_folder <path>: Enter the genome folder you wish to use to extract sequences from (full path only).
```
bismark对不同状态的胞嘧啶进行了约定并记录到了BAM文件和中间文件中：
```
X 代表CHG中甲基化的C
x  代表CHG中非甲基化的C
H 代表CHH中甲基化的C
h  代表CHH中非甲基化的C
Z  代表CpG中甲基化的C
z  代表CpG中非甲基化的C
U 代表其他情况的甲基化C(CN或者CHN)
u  代表其他情况的非甲基化C (CN或者CHN)
```
基于该文件可进行单个胞嘧啶位点的甲基化水平（C/C+T）计算，比如第一行：甲基化是水平为1/(1+6)；针对指定区域的甲基化水平计算则取区域内的胞嘧啶的总C/(总C+总T)即可。

结果：
```
Bismark report for: ./TetTKO_mESC_rep1/SRR7368845.fastq.gz (version: v0.19.1)
Option '--directional' specified (default mode): alignments to complementary strands (CTOT, CTOB) were ignored (i.e. not performed)
Bismark was run with Bowtie 2 against the bisulfite genome of /mnt/c/shengxin/NBT_repeat/data/genome_data/ with the specified options: -q --score-min L,0,-0.2 --ignore-quals
Final Cytosine Methylation Report
=================================
Total number of C's analysed:	1197395

Total methylated C's in CpG context:	12638
Total methylated C's in CHG context:	64471
Total methylated C's in CHH context:	243730
Total methylated C's in Unknown context:	254

Total unmethylated C's in CpG context:	52857
Total unmethylated C's in CHG context:	162711
Total unmethylated C's in CHH context:	660988
Total unmethylated C's in Unknown context:	1098

C methylated in CpG context:	19.3%
C methylated in CHG context:	28.4%
C methylated in CHH context:	26.9%
C methylated in Unknown context (CN or CHN):	18.8%
```

# Downstream analysis-下游分析 
Based on the methylation information result got from the methylation analysis procedure, we could make some basic downstream analysis work including finding specific locus and detecting differential methylation loci (DML) or differential methylation regions (DMR). Here we use R package DSS to make the differential methylation analysis.

根据甲基化分析过程中得到的甲基化信息，我们可以进行一些基本的下游分析工作，包括寻找特定位点和检测差异甲基化位点(DML)或差异甲基化区域(DMR)。这里我们使用R包DSS进行差异甲基化分析。

## Input data preparation
DSS (Dispersion Shrinkage for Sequencing data)，为基于高通量测序数据的差异分析而设计的Bioconductor包。主要应用于BS-seq（亚硫酸氢盐测序）中计算不同组别间差异甲基化位点（DML）和差异甲基化区域（DMR）即Call DML or DMR。

DSS要求每个BS-seq样实验的数据被总结为每个CG位置的以下信息:染色体数、基因组坐标、总的读取数和显示甲基化的读取数。

所需的输入数据可以从bismark result ```.cov```文件传输，因为count文件包含以下列:chr、开始、结束、甲基化%、甲基化计数、未甲基化计数。

DSS上游的数据来自于bismark软件的结果，这里通过转换bismark_methylation_extractor的cov结果文件来获得输入文件。

可以看一下上步得到的.cov文件
```
less SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz

结果
17      867841  867841  100     1       0
17      903527  903527  100     1       0
17      903531  903531  100     1       0
17      903571  903571  100     1       0
17      1241484 1241484 100     1       0
17      1241496 1241496 100     1       0
17      1241555 1241555 100     1       0
17      1241569 1241569 100     1       0
17      1241613 1241613 100     1       0
17      1361928 1361928 0       0       1
17      1400100 1400100 0       0       1
17      1400113 1400113 0       0       1

每一行代表一个CpG site
第一列为染色体
第二列为位置
第三列为total reads
第四列为甲基化的reads
第五列为甲基化数目
第六列为未甲基化数目
```
```
cd /mnt/c/shengxin/NBT_repeat
mkdir -p R_analysis/WT_data
mkdir -p R_analysis/TetTKO_data

# store the file path information
file_WT_path="/mnt/c/shengxin/NBT_repeat/data/seq_data/WT_mESC_rep1/deduplicated_result/SRR7368841_bismark_bt2.deduplicated.bismark.cov.gz"

file_TetTKO_path="/mnt/c/shengxin/NBT_repeat/data/seq_data/TetTKO_mESC_rep1/deduplicated_result/SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz"

# copy the result file to the R analysis folder
cp ${file_WT_path} ./R_analysis/WT_data
cp ${file_TetTKO_path} ./R_analysis/TetTKO_data

# unzip
gunzip -d ./R_analysis/WT_data/SRR7368841_bismark_bt2.deduplicated.bismark.cov.gz #gunzip 是个使用广泛的解压缩程序；-d：解压缩文件
gunzip -d ./R_analysis/TetTKO_data/SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz

# transfer the .cov file to .txt file
cp ./R_analysis/WT_data/SRR7368841_bismark_bt2.deduplicated.bismark.cov ./R_analysis/WT_data/SRR7368841_methylation_result.txt
cp ./R_analysis/TetTKO_data/SRR7368845_bismark_bt2.deduplicated.bismark.cov ./R_analysis/TetTKO_data/SRR7368845_methylation_result.txt

# basic command line environment procedure 
R

library(tidyr) #tidyr 包提供了用于数据整理和重塑的函数
library(dplyr) #dplyr 包提供了用于数据处理和操作的函数
library(DSS) #提供了用于甲基化数据分析的函数
报错：没有DSS包
解决方案：
R
install.packages("BiocManager")
BiocManager::install("DSS", force = TRUE)
library(DSS)

cd /mnt/c/shengxin/NBT_repeat/R_analysis
first_file <- "./WT_data/SRR7368841_methylation_result.txt"
second_file <- "./TetTKO_data/SRR7368845_methylation_result.txt"
file_prefix <- "mm_all_chr" #prefix选项是配置安装的路径
file_save_path <- "./"

# import data
#用R的 read.table 函数来读取名为 first_file 的文件，并将其保存到一个名为 first_raw_data 的数据框中
first_raw_data <- read.table(first_file, header = T, stringsAsFactors = F)#read.table 函数用于从文本文件中读取数据，并将其存储为数据框的形式;header = T，true是数据是含表头（列名，也称变量名）,即指定第一行为列名，将其视为文件的头部。
second_raw_data <- read.table(second_file, header = T, stringsAsFactors = F) #stringsAsFactors = F意味着，“在读入数据时，遇到字符串之后，不将其转换为factors，仍然保留为字符串格式”

# data manipulation to prepare for the BSseq objection
#得到一个名为 DSS_first_input_data 的数据框，其中包含了经过处理后的数据，只保留了 chr、pos、N 和 X 四列。
DSS_first_input_data <- first_raw_data %>%
	mutate(chr = paste("chr", chr, sep = "")) %>%    
    # mutate 创建新列；将 chr 列的值转换为以 "chr" 开头的形式，例如将 "1" 转换为 "chr1"
	mutate(pos = start, N = methyled + unmethyled, X = methyled) %>%
    #使用 mutate 函数将指定列的值重新赋给新的列名，并对其中的一些列进行计算和操作。pos 列的值将与 start 列相同，N 列的值将是 methyled 列和 unmethyled 列之和，X 列的值将与 methyled 列相同。
 #paste函数是比较常用字符串处理函数，可以连接不同类型的变量及常量。基本语法如下：
paste(..., sep = " ", collapse = NULL)
其中，…表示一个或多个R可以被转化为字符型的对象；sep表示分隔符，默认为空格；
	select(chr, pos, N, X)
#select 函数选择包含指定列的数据框。这里选取了 chr、pos、N 和 X 列

DSS_second_input_data <- second_raw_data %>% #％>％来自dplyr包的管道函数，是将前一步的结果直接传参给下一步的函数
	mutate(chr = paste("chr", chr, sep = "")) %>%
	mutate(pos = start, N = methyled + unmethyled, X = methyled) %>%
	select(chr, pos, N, X)
# 报错，找不到命令%>%
解决方案
R
library(dplyr)

# create BSseq object and make the statical test
#对两个输入数据进行甲基化差异分析
bsobj <- makeBSseqData(list(DSS_first_input_data, DSS_second_input_data), c("S1", "S2"))#使用了标识符 "S1" 表示第一个样本，"S2" 表示第二个样本。将结果保存在 bsobj 中
dmlTest <- DMLtest(bsobj, group1 = c("S1"), group2 = c("S2"), smoothing = T)#DMLtest:对亚硫酸酯测序(BS-seq)数据的两组比较进行差异甲基化位点(DML)的统计检验功能。该函数接受一个BSseq对象和两个组标签，然后对每个CpG位点的差异甲基化进行统计检验。DMLtest 函数对 bsobj 执行甲基化差异检验。group1 参数指定属于第一组的样本标识符，这里为 "S1"。group2 参数指定属于第二组的样本标识符，这里为 "S2"。smoothing 参数设置为 TRUE，表示进行平滑处理

# dmls detection #DMLtest函数call DML
dmls <- callDML(dmlTest, p.threshold = 0.001) #p.threshold = 0.001 参数来指定一个阈值，即 p 值的阈值。只有具有显著差异的甲基化位点，其 p 值低于或等于这个阈值时，才会被调用为 DML
# dmrs detection
dmrs <- callDMR(dmlTest, p.threshold=0.01)

# output the results
write.table(dmlTest, paste(file_save_path, file_prefix, "_DSS_test_result.txt", sep = ""), row.names = F)
#第二个参数是要保存文件的路径和文件名。你使用 paste 函数将 file_save_path（保存路径）、file_prefix（文件前缀）和 "_DSS_test_result.txt" 进行拼接，生成完整的文件名
#row.names = F 参数表示不包含行号
write.table(dmls, paste(file_save_path, file_prefix, "_DSS_dmls_result.txt", sep = ""), row.names = F)
write.table(dmrs, paste(file_save_path, file_prefix, "_DSS_dmrs_result.txt", sep = ""), row.names = F)









