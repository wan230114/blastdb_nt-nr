#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-09-17, 10:57:55
# @ Modified By: Chen Jun
# @ Last Modified: 2021-09-17, 11:49:04
#############################################

# %%
from Bio.Blast import NCBIXML
from copy import deepcopy
from datetime import datetime as dt
from Bio.Blast import NCBIWWW
fasta_string = open("./1000reads_BGC823-con-P65_R1.fasta").read()
t1 = dt.now()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)
print("耗时：", dt.now() - t1)

# %%
result_handle2 = deepcopy(result_handle)

res = result_handle2.read()
print(res, file=open("out.xml", "w"))
res
# %%
result_handle2 = deepcopy(result_handle)
# blast_record = NCBIXML.read()
i = 0
for blast_record in NCBIXML.parse(result_handle2):
    # print(dir(blast_record))
    i += 1
    print(blast_record.__dict__)
    # input()
    # if i > 4:
    #     break
    break
