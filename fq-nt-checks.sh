#!/usr/bin/env bash
#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-09-09, 17:55:51
# @ Modified By: Chen Jun
# @ Last Modified: 2021-10-29, 11:17:03
#############################################

Usage="Usage: sh xxx.sh  -i f1.fq.gz [f2.fq.gz ...]  [-n <int>] [-r]"

num=1000
rnaseq=false

while getopts ":i::n:r" opt; do
    case $opt in
    n) num=$OPTARG ;;
    r) rnaseq=true ;;
    i)
        infqs=("$OPTARG")
        until [[ $(eval "echo \${$OPTIND}") =~ ^-.* ]] || [ -z $(eval "echo \${$OPTIND}") ]; do
            infqs+=($(eval "echo \${$OPTIND}"))
            OPTIND=$((OPTIND + 1))
        done
        ;;
    *)
        echo $Usage
        exit 1
        ;;
    esac
done


if [ "${infqs}" ]; then
    echo infqs: ${infqs[@]}
else
    echo $Usage
    exit 1
fi

echo num: $num
echo rna: $rnaseq


#############################################
source /home/chenjun/.conda_bashrc_my; conda activate nt

date "+%F %T"
echo step01. filter $num reads.

fas=""
for fastq in ${infqs[@]}; do
    echo deal fastq: $fastq
    fastq_name=`basename $fastq .gz`
    /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/get_1000reads.sh $fastq $num
    fasta=${num}reads_${fastq_name%.fq}.fasta
    fas="$fas $fasta"
done
echo $fas


#############################################
date "+%F %T"
echo step02. blast.

if $rnaseq; then
    echo /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/blast-all/blast.py $fas -r
    echo /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/blast-all/blast.py $fas -r| sh
else
    echo /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/blast-all/blast.py $fas
    echo /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/blast-all/blast.py $fas | sh
fi

#############################################
# 统计, plot
date "+%F %T"
echo step03. tongji.

for fa in $fas; do
    out=${fa%.fa}.blast
    /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py $out
done
