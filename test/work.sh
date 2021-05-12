#!/usr/bin/env bash
#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-05-12, 21:40:18
# @ Modified By: Chen Jun
# @ Last Modified: 2021-05-12, 23:51:03
#############################################

time /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/map_taxid.py /home/chenjun/pipeline/others/dataBase/blastdb_nt-nr/test.blast
# reading database ...
# searching...
# plotting...
# real    0m12.926s
# user    0m9.968s
# sys     0m2.923s

time {
    cat test.blast | sort -k1,1 -u >test.filter.blast
    cat test.filter.blast |cut -f 2 >keys.txt
    # 下面这一步占了99%的时间，为1m20s
    awk 'NR==FNR{a[$1]=$2;next}NR!=FNR{if($1 in a)print $0}'  keys.txt  /home/chenjun/dataBase/blast_db_FASTA/db_NCBI_Taxonomy/nucl_gb.accession2taxid.filter.sorted.txt >find.txt
    cut -f 2 find.txt  >find.txt.col1
    awk 'NR==FNR{a[$1]=$0;next}NR!=FNR{if($1 in a)print a[$0]}'  /home/chenjun/dataBase/blast_db_FASTA/db_NCBI_Taxonomy/taxdump/names.dmp.db find.txt.col1 >find_res.xls
    # \rm test.filter.blast keys.txt find.txt.col1
}
# real    1m22.207s
# user    1m20.908s
# sys     0m0.939s

