#!/usr/bin/env bash
#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-09-09, 17:55:51
# @ Modified By: Chen Jun
# @ Last Modified: 2021-10-21, 09:19:28
#############################################

while getopts ":i::n:" opt; do
    case $opt in
    n) num=$OPTARG ;;
    i)
        infqs=("$OPTARG")
        until [[ $(eval "echo \${$OPTIND}") =~ ^-.* ]] || [ -z $(eval "echo \${$OPTIND}") ]; do
            infqs+=($(eval "echo \${$OPTIND}"))
            OPTIND=$((OPTIND + 1))
        done
        ;;
    *)
        echo "Usage: sh xxx.sh  -i f1.fq.gz [f2.fq.gz ...]  -n <int>"
        exit 1
        ;;
    esac
done


if [ "${infqs}" ]; then
    echo infqs: ${infqs[@]}
else
    echo "Usage: sh xxx.sh  -i f1.fq.gz [f2.fq.gz ...]  [-n <int>]"
    exit 1
fi
if [ "$num" ]; then echo -n ; else num=1000; fi

echo num: $num


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

echo /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/blast-all/blast.py $fas
echo /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/blast-all/blast.py $fas | sh

#############################################
# 统计, plot
date "+%F %T"
echo step03. tongji.

for fa in $fas; do
    out=${fa%.fa}.blast
    /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py $out
done
