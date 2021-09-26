#!/usr/bin/env bash
#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-05-11, 17:40:03
# @ Modified By: Chen Jun
# @ Last Modified: 2021-09-18, 12:18:11
#############################################

fastq=$1  # xxx.fq.gz
fastq_name=`basename $fastq .gz`
num=1000
if [ "$2" ]; then num=$2; fi
source /home/chenjun/.conda_bashrc_my; conda activate nt

date "+%F %T"
echo step01. filter $num reads.
/home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/get_1000reads.sh $fastq $num

# 比对
date "+%F %T"
echo step02. blast.
db=/home/chenjun/dataBase/blast_db_FASTA/nt
fasta=${num}reads_${fastq_name%.fq}.fasta
out=${num}reads_${fastq_name%.fq}.blast
time blastn -query $fasta -out $out -outfmt 6 -db $db -num_threads 10 -num_alignments 5
# time blastn -query $fasta -out $out -outfmt 6 -db $db -num_threads 10 -evalue 1e-5  -qcov_hsp_perc 50.0 -num_alignments 5
# time blastn -query $fasta -out $out -max_target_seqs 1 -outfmt 6 -db $db -num_threads 10 -evalue 1e-3

# 统计, plot
date "+%F %T"
echo step03. tongji.
/home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py $out
# /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py 9HTF2HZF013-Alignment-HitTable.txt
