#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-05-10, 11:07:39
# @ Modified By: Chen Jun
# @ Last Modified: 2021-12-31, 15:03:32
#############################################

##########################################################################################
# Description：将物种对应名的数据库分割，按需读取查找，极大减少资源消耗与读取速度
#              直接输入 blast f6 结果，脚本将自动筛选第一次出现的最优比对结果，
#              将ID去判断所属物种，并进行饼图绘制
# Version： v1.2  修复比对到的物种数量较少时报错的BUG
# Version History：
#   v1.0  初版完成
#   v1.1  加入参数可以定义输入的原始reads数量，多余将视为unmap; 规范饼图最多显示条目为7
##########################################################################################


# %%
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
import pickle
import gzip
import os
import argparse


def fargv():
    parser = argparse.ArgumentParser(description='Process introduction.')
    parser.add_argument('blastfiles', type=str, nargs="+",
                        help='输入需要统计的 blast 结果')
    parser.add_argument('-n', '--numrawreads', type=int, default=None,
                        help='原始的reads数，多余的将视为unmap')
    args = parser.parse_args()
    return args


pydb = "/home/chenjun/dataBase/blast_db_FASTA/db_NCBI_Taxonomy/pydb"


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


def dict_split(d_tmp):
    """用于切分字典的函数
    Args:
        d_tmp (dict): {x:1, y:2, y2:3}
    Returns:
        name1, dx1
    """
    keys = list(d_tmp)
    # endpos = len(name)+1
    endpos = 2
    name1 = name2 = keys[0][0:endpos]
    dx1 = {}
    while not dx1 or len(dx1) == len(d_tmp):
        dx1 = {}
        p = 0
        name1 = name2 = keys[0][0:endpos]
        # print("--->", len(dx1), len(d_tmp))
        for k in keys:
            name2 = k[0:endpos]
            # print("-->", name1, name2)
            if name1 != name2:
                p = 1
                break
            dx1[k] = d_tmp[k]
        if p:
            break
        endpos += 1
        name1 = name2
    # print("split:", name1, dx1)
    return name1, dx1


def pickle_load(path):
    with gzip.GzipFile(path, "rb") as f:
        d_tmp = pickle.load(f)
    return d_tmp


def pickle_dump(var, out_db):
    # print(*(var), file=open(out_db + ".tmp", "w"), sep="\n")
    print("[dumping var]: %s" % out_db, len(var))
    if os.path.exists(out_db + ".gz"):
        raise Exception("%s文件已存在，请检查是否重复运行，或输入文件未排序导致" % out_db)
    with open(out_db, "wb") as f:
        pickle.dump(var, file=f)
    os.system("cat %s|gzip >%s.gz && rm %s" % (out_db, out_db, out_db))


def find_spe(IDs):
    # IDs = ["XM_014017284", "A00144"]
    res_dict = {}
    db_select = {}
    for ID in IDs:
        for i in range(len(ID), 0, -1):
            kw = ID[:i]
            if kw in db_info:
                break
        db_select.setdefault(kw, []).append(ID)
    # L_time = []
    # print("find index ok")
    import datetime
    t00 = datetime.datetime.now()
    for kw in db_select:
        print(f"search in {kw}", end="  ")
        t0 = datetime.datetime.now()
        try:
            db_tmp = pickle_load(f"{pydb}/%s.pydb.gz" % kw)
        except FileNotFoundError:
            print(f"  ---  Warning: {pydb}/%s.pydb.gz" % kw, "not found. this ID not mapped:", db_select[kw])
            db_tmp = {}
        t1 = datetime.datetime.now()
        print(t1 - t0)
        # L_time.append(t1 - t0)
        for ID in db_select[kw]:
            if ID in db_tmp:
                res_dict[ID] = d1[int(db_tmp[ID])]
            else:
                res_dict[ID] = "Unknown"
    t11 = datetime.datetime.now()
    print("used time:", t11 - t00)
    # x = L_time[0]
    # for xx in L_time[1:]:
    #     x += xx
    # print(x)
    return list(zip(IDs, [res_dict[x] for x in IDs]))
# find_spe(["XM_014017284", "A00144"])
# %%


def blast_res_deal(inblast, numrawreads=None):
    # inblast = "./test/test.filter.blast"
    names = [
        "Query_id",  # 查询序列ID标识
        "Subject_id",  # 比对上的目标序列ID标识
        "%identity",  # 序列比对的一致性百分比
        "alignment_length",  # 符合比对的比对区域的长度
        "mismatches",  # 比对区域的错配数
        "gap_openings",  # 比对区域的gap数目
        "q_start",  # 比对区域在查询序列(Query id)上的起始位点
        "q_end",  # 比对区域在查询序列(Query id)上的终止位点
        "s_start",  # 比对区域在目标序列(Subject id)上的起始位点
        "s_end",  # 比对区域在目标序列(Subject id)上的终止位点
        "e-value",  # 比对结果的期望值，将比对序列随机打乱重新组合，和数据库进行比对，如果功能越保守，则该值越低；该E值越高说明比对的高得分值是由GC区域，重复序列导致的。对于判断同源性是非常有意义的几个参数。
        "bit_score",  # 比对结果的bit score值
    ]
    df = pd_read_table_str(inblast, names=names)
    df = df.drop_duplicates(["Query_id"])

    df["Subject_id_clean"] = df["Subject_id"].map(
        lambda x: x.split(".")[0])
    df["Subject_id_clean"]
    df["scientific_name"] = [x[1] for x in find_spe(df["Subject_id_clean"])]
    df.to_csv(f"{inblast}.HighestScore.xls", sep="\t", index=False)

    tongji = df["scientific_name"].value_counts()
    if numrawreads and tongji.sum() < numrawreads:
        tongji.loc["unmapped reads"] = numrawreads - tongji.sum()

    # tongji.plot(kind="pie")
    df_tongji = pd.DataFrame(tongji)
    df_tongji.insert(1, "per", tongji / tongji.sum())

    df_tongji2 = df_tongji.copy()
    df_tongji2["per"] = (df_tongji2["per"]*100).map(lambda x: round(x, 3))
    # df_tongji2["per"] = df_tongji2["per"].map(lambda x: "%6.3f" % x)
    df_tongji2.rename({"scientific_name": "num", "per": "%per"},
                      axis=1, inplace=True)
    tmp = df_tongji2.copy()
    tmp.loc["Sum of reads:"] = [tmp["num"].sum(), 100]
    tmp["num"] = tmp["num"].astype("int")
    tmp.to_csv(f"{inblast}.HighestScore.Description.xls",
               sep="\t", index_label="sep")
    tmp.to_html(f"{inblast}.HighestScore.Description.html")
    if df_tongji.shape[0] > 8:
        # other_select = df_tongji["per"] < 0.01
        other_select = np.array([False]*df_tongji.shape[0])
        other_select[8:] = True
        df_tongji_other = df_tongji[other_select]
        min_num = df_tongji_other.iloc[0, 1]*100
        min_num = 0.1 if min_num < 0.1 else min_num
        othersname = "others reads (sum of <%.1f%%)" % min_num
        # df_tongji_other = df_tongji[tongji <= 2]
        # df_tongji.sum()
        if "unmapped reads" in df_tongji_other.index:
            unmap = df_tongji.loc["unmapped reads"]
            df_tongji.drop(df_tongji_other.index, inplace=True)
            df_tongji.loc[othersname] = df_tongji_other.sum()
            df_tongji.loc["unmapped reads"] = unmap
        else:
            df_tongji.drop(df_tongji_other.index, inplace=True)
            df_tongji.loc[othersname] = df_tongji_other.sum()
    # df_tongji["scientific_name"].plot(kind="pie")

    df_tongji["name"] = (
        df_tongji.index
        + "\n("
        + (df_tongji["per"]*100).map(lambda x: "%6.3f" % x)
        + "%, "
        + df_tongji["scientific_name"].map(lambda x: "%d" % x)
        + ")"
    )
    return df_tongji


def plot_pie(data, ingredients, outname):
    # data = [1,2,3]
    # ingredients = ["A","B","C"]
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    wedges, texts = ax.pie(data,
                           #   autopct=lambda pct: '%.2f%%' % pct,
                           #   textprops=dict(color="w")
                           colors=cm.get_cmap("Paired")(
                               np.arange(len(ingredients)))
                           )
    ax.legend(wedges, ingredients,
              #   title="Ingredients",
              #   loc="center left",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize='x-small'
              #   bbox_to_anchor=(0, 0, 0.5, -0.5)
              )
    # plt.setp(outname, size=8)
    ax.set_title(os.path.basename(outname), size=8, loc="left")
    # fig.suptitle(os.path.basename(outname), size=8, loc="right")
    # plt.show()
    plt.savefig(outname+".HighestScore.Description.pie.png",
                dpi=200, bbox_inches='tight')


# df_tongji = blast_res_deal("./test/test.filter.blast", 1000)
# plot_pie(data=df_tongji.scientific_name,
#          ingredients=df_tongji.name,
#          outname="./test/test.filter.blast")

# %%


# %%
print("reading database ...")


if os.path.exists(pydb):
    d1 = pickle_load(f"{pydb}/taxdum_names.pydb.gz")
    db_info = dict(pd.read_table(f"{pydb}/numbers.tsv", header=None).values)
else:
    print("no database, building database ...")
    os.makedirs(pydb, exist_ok=True)
    taxdump_names = "%s/taxdump/names.dmp" % (os.path.dirname(pydb))
    print(f"building {taxdump_names}")
    df1 = pd_read_table_str(taxdump_names, header=None)
    df1 = df1[df1[6] == "scientific name"][[0, 2]].rename(
        {0: "taxid", 2: "scientific_name"}, axis=1)
    df1["taxid"] = df1["taxid"].astype("int")
    # df[df["taxid"] == 9606]
    d1 = dict(df1.values)
    del df1
    pickle_dump(d1, f"{pydb}/taxdum_names.pydb")

    print("build gb.accession2taxid.")
    os.makedirs(pydb, exist_ok=True)
    # finame_taxid = "XR.txt"   # tail -10000000  nucl_gb.accession2taxid.filter.sorted.txt >XR.txt
    # finame_taxid = "test.txt"
    # finame_taxid = "test.sorted.txt"
    # finame_taxid = "nucl_gb.accession2taxid.filter.txt"
    finame_taxid = "nucl_gb.accession2taxid.filter.sorted.txt"
    print("building %s" % finame_taxid)
    fi = open(finame_taxid)

    split_length = {}
    # 分1000份文件，平均每份100W
    # 分500份文件，平均每份200W
    # 分100份文件，平均每份1000W
    # pd.DataFrame(d_tmp.items())[1].sum()/300
    cutoff = 10 ** 6
    # cutoff = 10
    num = 0
    # cengji = 1
    d_tmp = {}
    # def split_dict(in_tmp):
    # split_dict(in_tmp)
    old_name = ""
    while True:
        line = fi.readline()
        if not line:
            break
        k, v = line.strip().split("\t")
        name = k[0]
        if num > cutoff:
            # TODO - [ ] 逻辑日后还可进一步优化，因为当大字段遍历结束后，会将后续的新字段不断添加进来
            # TODO - [ ] 最坏的情况是每次只分割出一个字段，然后后面每次都添加一个字段
            # 此处有200w，分割之后必然少于200w
            name1, dx1 = dict_split(d_tmp)
            # cengji = len(name1)
            for x in dx1:
                d_tmp.pop(x)
            num = len(d_tmp)
            pickle_dump(dx1, os.path.join(pydb, name1 + ".pydb"))
            split_length[name1] = len(dx1)
        elif old_name and old_name != name:  # ! 该如何判断进入下一个, 大字段为节, 统一写入最后一个
            pickle_dump(d_tmp, os.path.join(pydb, old_name + ".pydb"))
            split_length[old_name] = len(d_tmp)
            d_tmp = {}
            num = 0
        # print("old_name:", old_name, "name:", name)
        # print("-->", k, v)
        d_tmp[k] = int(v)
        num += 1
        old_name = name
    pickle_dump(d_tmp, os.path.join(pydb, old_name + ".pydb"))
    split_length[old_name] = len(d_tmp)

    fi.close()
    pd.DataFrame(list(split_length.items())).to_csv(
        f"{pydb}/numbers.tsv", sep="\t", index=False, header=False)
    # 验证：
    # cat ../XR.txt |cut -f 1 >t1
    # ls -rt *tmp|xargs cat  >t2
    # diff t1 t2

# %%

if __name__ == "__main__":
    # inblast = "blast_result.blast"
    # inblast = "/home/chenjun/test/blast/test-fq/9HTF2HZF013-Alignment-HitTable.txt"
    # inblast = sys.argv[1]
    # infaNum = sys.argv[2]
    # sys.argv = ["", "1000reads_KO-C-Jun-k27me3.rep2.chip.ip_R1.fasta.blast"]
    args = fargv()
    for inblast in args.blastfiles:
        print(f"\ndeal {inblast}")
        print("searching...")
        df_tongji = blast_res_deal(inblast, args.numrawreads)
        print("plotting...")
        plot_pie(data=df_tongji.scientific_name,
                 ingredients=df_tongji.name,
                 outname=inblast)
        print("Done.")
