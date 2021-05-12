#!/usr/bin/env bash
#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-02-22, 16:12:38
# @ Modified By: Chen Jun
# @ Last Modified: 2021-05-13, 01:23:58
#############################################

# 该脚本用于随机挑选 1000条 fastq reads, 检测nt数据库比对情况，确定是否存在异源污染

# ref：
# 随机读取reads使用ncbi在线工具比对nt数据库_穆易青的博客-CSDN博客
# https://blog.csdn.net/yangl7/article/details/111994049

# fastq=./HX-A549-BTF3_R1.fq.gz
fastq=$1
fastq_name=`basename $fastq .gz`
num=1000
if [ "$2" ]; then num=$2; fi

source /home/chenjun/.conda_bashrc_my; conda activate nt
seqtk sample -s 100 $fastq ${num}  >${num}reads_${fastq_name}
cat ${num}reads_${fastq_name} | awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' >${num}reads_${fastq_name%.fq}.fasta

echo "then, you can open link :  https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
# then, open link :  https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
