#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/2018 1:32 PM
# @Author  : Runsheng     
# @File    : tracklist.py

"""
Functions to handel the IO of bigg list
"""

# self import
from track import bigGenePred
from collections import OrderedDict
from utils import myexe, set_tmp, del_files
import pandas


def add_sw(bigg_file, sw_file, out="bigg_sw.bed"):
    """
    To add the splicing leader mapping score
    :param bigg_list:
    :return:
    """
    sw_dic={}
    with open(sw_file, "r") as f:
        for line in f.readlines():
            line_l=line.split(",")
            name=line_l[0]
            score=int(line_l[1])
            sw_dic[name]=score

    bigg_list= read_bigg(bigg_file)

    bigg_str=[]
    for bigg_one in bigg_list:
        score=sw_dic[bigg_one.name]
        bigg_one.score=score
        bigg_str.append(bigg_one.to_str())

    with open(out, "w") as fw:
        fw.write("\n".join(bigg_str))


def read_bigg(bigg_file):
    """

    :param bigg_file:
    :return: bigg_list
    """
    bigg_list=[]
    with open(bigg_file, "r") as f:
        for line in f.readlines():
            bigg_one=bigGenePred()
            bigg_one.from_string(line.strip())
            bigg_list.append(bigg_one)

    return bigg_list


def write_bigg(bigg_list, out="bigg_new.bed"):

    bigg_str=[]
    for bigg_one in bigg_list:
        bigg_str.append(bigg_one.to_str())

    with open(out, "w") as fw:
        fw.write("\n".join(bigg_str))


def list_to_dic(bigg_list):

    bigg_dic=OrderedDict()

    for i in bigg_list:
        bigg_dic[i.name]=i
    return bigg_dic


def bigglist_to_bedfile(bigg_list,prefix=None, dir=None):

    bigg0=bigg_list[0]
    if prefix is None:
        prefix=bigg0.name
    if dir is None:
        dir=set_tmp()

    out_exon=dir+"/{prefix}_exon.bed".format(prefix=prefix)
    out_intron=dir+"/{prefix}_intron.bed".format(prefix=prefix)

    f_exon=open(out_exon, "w")
    f_intron=open(out_intron, "w")

    for bigg in bigg_list:
        bigg.to_bedstr(gene_start=0)

        f_exon.write(bigg.exon_str)
        f_exon.write("\n")
        f_intron.write(bigg.intron_str)
        f_intron.write("\n")

    return (out_exon, out_intron)


def get_file_prefix(filepath):
    return filepath.split("/")[-1].split(".")[0]


def get_file_location(filepath):
    return "/".join(filepath.split("/")[0:-1])


def wrapper_bedtools_intersect2(bedfile1,bedfile2,outfile=None):
    """
    Using two bedfile to get the intsersection of pairs
    :param bigg_one:
    :param bigg_two:
    :return:
    """
    if outfile is None:
        prefix1=get_file_prefix(bedfile1)
        prefix2=get_file_prefix(bedfile2)
        location=get_file_location(bedfile1)

        outfile=location+"/"+"_".join([prefix1, prefix2])+".bed"

    sort_cmd1="bedtools sort -i {bed} > {bed}_s".format(bed=bedfile1)
    sort_cmd2="bedtools sort -i {bed} > {bed}_s".format(bed=bedfile2)

    _ = myexe(sort_cmd1)
    _ = myexe(sort_cmd2)

    # generate the bedfile

    cmd="bedtools intersect -wa -wb -a {bedfile1}_s -b {bedfile2}_s>{out}".format(
        bedfile1=bedfile1, bedfile2=bedfile2, out=outfile)

    _=myexe(cmd)

    ### cleanup
    bed1s=bedfile1+"_s"
    bed2s=bedfile2+"_s"
    del_files([bedfile1, bedfile2, bed1s, bed2s])

    return outfile


def count_file(thefile):
    count = 0
    for line in open(thefile).xreadlines(  ):
        count += 1
    return count


def pandas_summary(bed8file):
    """
    The bef8file is chr start end name *2 format
    :param bed8file:
    :return: the dict with (read1, read2): intersection
    """
    if count_file(bed8file)==0:
        return {}

    df=pandas.read_csv(bed8file, sep="\t", header=None)
    df["start_max"] = df[[1, 5]].max(axis=1)
    df["end_min"] = df[[2, 6]].min(axis=1)
    df["sub"] = df["end_min"] - df["start_max"]

    dfs = df[[3, 7, "sub"]]
    # debug
    #dfs.to_csv("aa.csv")
    #print dfs
    #print(len(df))
    dfs.drop_duplicates(subset=[3,7])
    #print(len(df))
    groupdfs = dfs.groupby([3, 7])
    aa = groupdfs.sum()

    intersection_dic=aa.to_dict()["sub"]

    return intersection_dic


def boundary_correct(isoform_list, read_list):
    """
    use the subread information to correct the isoforms
    """
    read_dic=list_to_dic(read_list)

    for isoform in isoform_list:
        isoform.get_exon()
        isoform.get_coverage_from_str()
        isoform.get_subread_from_str()

        for name in isoform.subread:
            read=read_dic[name]
            read.get_exon()
            for pair in read.exon:
                exon_s,exon_e=pair

        ## assume that the read have the same 3' as the isoform
        ## map the last ones or the isoforms with equal length to correct the boundary




def bigg_count_write(bigg_list, out=None):
    """
    parser the output of cluster, get the count for each isoform
    :param bigg_list:
    :return:
    """
    # store sub-read name and number
    name_dic=OrderedDict()

    for bigg in bigg_list:
        bigg.get_subread_from_str()
        if len(bigg.subread)>0:
            for name in bigg.subread:
                if "-" in name: # judge if it is a read or a isoform
                    try:
                        name_dic[name]+=1
                    except KeyError:
                        name_dic[name]=1

    for bigg in bigg_list:
        coverage=0
        for name in bigg.subread:
            if "-" in name:
                coverage+=1.0/name_dic[name]

        bigg.coverage=coverage
        bigg.write_coverage()

    # debug
    #print name_dic
    write_bigg(bigg_list,out)






