#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-09-17, 10:57:55
# @ Modified By: Chen Jun
# @ Last Modified: 2021-12-31, 13:26:49
#############################################

# %%
from Bio import SeqIO
from datetime import datetime as dt
from collections import OrderedDict
import os
import pandas as pd
import argparse

#%%

parser = argparse.ArgumentParser(description='Process introduction.')
parser.add_argument('fas', type=str, nargs="+",
                    help='输入需要运行的fa')
parser.add_argument('-r', '--refseq', action='store_true',
                    default=False,
                    help='当给定此参数时，将使用refseq_rna作为blast的库(通常用于RNA检测)，否则默认使用nt库')
args = parser.parse_args()


# fas = ["./f1.fa", "./f2.fa", "./f3.fa"]
# refseq = False
fas = args.fas
refseq = args.refseq


D = OrderedDict()

db_nt = "/home/chenjun/dataBase/blast_db_FASTA/nt"
db_refseq = "/home/chenjun/dataBase/blast_db_refseq/refseq_rna"
fasta = "merge_infastas.fa"
out = "merge_infastas.fa" + ".blast"

with open(fasta, "w") as fo:
    for fa in fas:
        for i, seq_record in enumerate(SeqIO.parse(fa, "fasta"), start=1):
            newname = fa + "_" + str(i)
            D[newname] = seq_record.description
            print(">" + newname, file=fo)
            print(seq_record.seq, file=fo)
            # print(repr(seq_record.seq))
            # print(len(seq_record))

with open(f"{fasta}.info.txt", "w") as fo:
    for x, y in D.items():
        print(x, y, sep="\t", file=fo)

# %%

t0 = dt.now()
print("start blast", t0)
db = db_refseq if refseq else db_nt
CMD = f"blastn -query {fasta} -out {out} -outfmt 6 -db {db} -num_threads 10 -num_alignments 5  1>log.blast.txt 2>&1"
print(CMD)
stat = os.system(CMD)
if stat != 0:
    raise Exception("blast erro. please check file: `log.blast.txt`.")
t1 = dt.now()
print("end blast", t1)
print("Time used:", t1 - t0)


# %%
# 此部分单独用于merge文件分割

fasta = "merge_infastas.fa"
out = "merge_infastas.fa" + ".blast"


def pd_read_table_str(infile, dtype={}, **kwargs):
    # ! 二次读取的解决方案：先读第一行，将所有列指定为str，根据需求手动改类型，二次读取。
    colnames = pd.read_table(
        infile, nrows=1, keep_default_na=False, **kwargs).columns
    dtypes = {x: "str" for x in colnames}
    dtypes.update(dtype)
    # dtype.update({"c1": "int"})  # 根据需要更改
    df = pd.read_table(infile, dtype=dtypes, keep_default_na=False, **kwargs)
    # 如需将NA也识别为字符串，需指定参数，keep_default_na=False
    return df


df_tmp = pd_read_table_str(f"{fasta}.info.txt", header=None)
D2 = dict(df_tmp.values)
try:
    print("D2 == D: ", D2 == D)
except Exception:
    pass

df = pd_read_table_str(out, header=None)
df.insert(df.shape[1], "name",
          df[0].map(lambda x: "_".join(x.split("_")[:-1])))
df[0] = df[0].map(lambda x: D2[x])
for x, df_tmp in df.groupby(["name"]):
    print("-->", x + ".blast")
    df_tmp.drop("name", axis=1).to_csv(x+".blast", sep="\t", index=False, header=False)
