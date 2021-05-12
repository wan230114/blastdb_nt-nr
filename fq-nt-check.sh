#!/usr/bin/env bash
#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-05-11, 17:40:03
# @ Modified By: Chen Jun
# @ Last Modified: 2021-05-13, 01:34:01
#############################################

fastq=$1  # xxx.fq.gz
fastq_name=`basename $fastq .gz`
num=1000

source /home/chenjun/.conda_bashrc_my; conda activate nt

date "+%F %T"
echo step01. filter 1000 reads.
/home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/get_${num}reads.sh $fastq

# 比对
date "+%F %T"
echo step02. blast.
db=/home/chenjun/dataBase/blast_db_FASTA/nt
fasta=${num}reads_${fastq_name%.fq}.fasta
out=${num}reads_${fastq_name%.fq}.blast
time blastn -query $fasta -out $out -outfmt 6 -db $db -num_threads 10 -evalue 1e-5  -qcov_hsp_perc 50.0 -num_alignments 5
# time blastn -query $fasta -out $out -max_target_seqs 1 -outfmt 6 -db $db -num_threads 10 -evalue 1e-3

# 统计, plot
date "+%F %T"
echo step03. tongji.
/home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py $out
# /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py 9HTF2HZF013-Alignment-HitTable.txt
