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

# quality control to see the quality of the raw seq-data
fastqc --t 3 ./WT_mESC_rep1/*.fastq.gz ./TetTKO_mESC_rep1/*.fastq.gz

# quality and adapter trimming, followd by fastqc operation
trim_galore -o ./WT_mESC_rep1/trimmed_data/ --fastqc ./WT_mESC_rep1/*.fastq.gz

# Trim Galore used options:
-o/--output_dir <DIR>: If specified all output will be written to this directory instead of the current directory. If the directory doesn't exist it will be created for you.
--fastqc: Run FastQC in the default mode on the FastQ file once trimming is complete.
```
打开fastqc.html文件，可以看到C含量几乎没有，本样本测序测定了处理后的甲基化水平，C已经被变成U/T，只有5hmC仍被识别为C，所以是正常的
# 甲基化分析

After the quality control and adapter trimming, BS-seq data analysis protocal based on Bismark could be followed to do the methylation analysis for the ACE-seq data as well.

用法参考网址 

https://www.jianshu.com/p/877a5716c24a
https://www.zxzyl.com/archives/759/

bismark的分析分为三步：

基因组索引构建：需要将基因组进行相应的转换，模拟bisulfite处理后基因组，这样才能用于比对。在每个新物种比对之前，都需要进行该处理。bismark内部默认调用bowite软件进行比对。

序列比对：输入为下机序列、基因组和比对参数，bismark会输出比对结果和甲基化检测的结果，默认格式为BAM，同时还会有一个统计报告。

甲基化位点提取：这一步是可选的，会利用bismark的比对结果来获得甲基化信息。这一步会把甲基化分为不同的类别（CG，CHG和CHH），获得链特异性信息，并且提供过滤参数；也可以对甲基化信息进行调整，或者进行更多的深入分析。

# 1. 基因组索引

Before alignments can be carried out, the genome of interest needs to be bisulfite converted in-silico and indexed. Based on Bismark, we can use both Bowtie and Bowtie2 to do the indexing work.

#下载bismark

通过官网下载 https://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark

或者brew 安装

```brew install bismark```

```
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/

命令：bismark_genome_preparation [options] <path_to_genome_folder>
# bismark_genome_preparation --path_to_bowtie这里使用bowtie2 <path_to_genome_folder>   索引路径 
结果：在基因组目录下产生Bisulfite_Genome目录

bismark_genome_preparation --bowtie2 /mnt/c/shengxin/NBT_repeat/data/genome_data/

```
```bismark_genome_preparation``` used options:

```--bowtie2```: This will create bisulfite indexes for Bowtie 2. (Default: ON).

# 2. 比对
The core of the methylation data analysis procedure is to align the sequencing reads to the reference genome, and it is assumed that all data have been quality and adapter trimmed. The alignment reports could be found in sub-directory bismark-alignment-reports/.

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

Mammalian genomes are so huge that it is rather unlikely to encounter several genuinely independent fragments which align to the very same genomic position. It is much more likely that such reads are a result of PCR amplification. For large genomes, removing duplicate reads is therefore a valid route to take.

The deduplication step could be finished through the deduplicate_bismark script, and the deduplication report could be found in sub-directory deduplication-reports/.

```
cd /mnt/c/shengxin/NBT_repeat/data/seq_data/
mkdir -p ./WT_mESC_rep1/deduplicated_result/
mkdir -p ./TetTKO_mESC_rep1/deduplicated_result/

# aligned reads deduplication
deduplicate_bismark --bam [options] <filenames> 
--bam: The output will be written out in BAM format instead of the default SAM format.
--output_dir 输出目录 --bam  比对后的BAM输出文件

deduplicate_bismark --bam --output_dir ./WT_mESC_rep1/deduplicated_result/ ./WT_mESC_rep1/bismark_result/SRX4241790_trimmed_bismark_bt2.bam

deduplicate_bismark --bam --output_dir ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/bismark_result/*.bam
```

# Methylation information extracting-甲基化信息提取
The bimark_methylation_extractor script in the Bismark package could extract the methylation information from the alignment result files and act as the endpoint of the Bismark package.  In addition, the methylation information could be easily transfered into other format facilitating the downstream analysis in this script.

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
bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \
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
gunzip -d ./R_analysis/WT_data/SRR7368841_bismark_bt2.deduplicated.bismark.cov.gz
gunzip -d ./R_analysis/TetTKO_data/SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz

# transfer the .cov file to .txt file
cp ./R_analysis/WT_data/SRR7368841_bismark_bt2.deduplicated.bismark.cov ./R_analysis/WT_data/SRR7368841_methylation_result.txt
cp ./R_analysis/TetTKO_data/SRR7368845_bismark_bt2.deduplicated.bismark.cov ./R_analysis/TetTKO_data/SRR7368845_methylation_result.txt

# basic command line environment procedure 
R
```
```
library(tidyr)
library(dplyr)
file_names <- c("mnt/c/shengxin/NBT_repeat/R_analysis/WT_data/SRR7368841_methylation_result.txt", "./TetTKO_data/SRR7368845_methylation_result.txt")
func_read_file <- function(file_name){
	dir_vec <- strsplit(file_name, split = "/")[[1]]




