#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 8/24/2018 1:32 PM
# @Author  : Runsheng     
# @File    : tracklist.py

"""
Functions to handel the IO of bigg list
"""
import glob
import re
# self import
from trackcluster.track import bigGenePred
from collections import OrderedDict
from trackcluster.utils import myexe, set_tmp, del_files, count_file, get_file_prefix, get_file_location
from random import randint
import pandas


def add_sw(bigg_file, sw_file, out="bigg_sw.bed"):
    """
    To add the splicing leader mapping score
    sw_file is a two col file: {readname, score}
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


def bigg_addvalue(bigg_file, value_str, poskey, out="bigg_add.bed"):
    """
    add a single value to a pos, for each of the bigg in the file,write the new file out
    pos is designed to be score, GeneName, GeneName2(group), ttype(isoform_anno, read_nano, region_mark), name2 (subreads)
    :param bigg_list:
    :return:
    """
    bigg_list= read_bigg(bigg_file)
    bigg_new=[]
    for bigg_one in bigg_list:
        bigg_one.poskey=value_str
        bigg_new.append(bigg_one)
    write_bigg(bigg_new, out)


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
        fw.write("\n")


def list_to_dic(bigg_list):

    bigg_dic=OrderedDict()

    for i in bigg_list:
        bigg_dic[i.name]=i
    return bigg_dic


def remove_special_chars(my_str):
    my_new_string = re.sub('[^a-zA-Z0-9 \n\.]', '_', my_str)
    return my_new_string


def bigglist_to_bedfile(bigg_list,prefix=None, dir=None):
    """
    write all intron and all exon region to a bed file for one gene
    :param bigg_list:
    :param prefix:
    :param dir:
    :return:
    """
    surfix=str(randint(100, 999))
    bigg0=bigg_list[0]
    if prefix is None:
        prefix=bigg0.name
    if dir is None:
        dir=set_tmp()

    prefix=remove_special_chars(prefix)

    out_exon=dir+"/{prefix}_{surfix}_exon.bed".format(prefix=prefix, surfix=surfix)
    out_intron=dir+"/{prefix}_{surfix}_intron.bed".format(prefix=prefix, surfix=surfix)

    f_exon=open(out_exon, "w")
    f_intron=open(out_intron, "w")

    for bigg in bigg_list:
        bigg.to_bedstr(gene_start=0)

        f_exon.write(bigg.exon_str)
        f_exon.write("\n")
        f_intron.write(bigg.intron_str)
        f_intron.write("\n")

    f_exon.close()
    f_intron.close()

    return (out_exon, out_intron)


def wrapper_bedtools_intersect2(bedfile1,bedfile2,outfile=None):
    """
    Using two bedfile to get the intsersection of pairs
    :param bigg_one:
    :param bigg_two:
    :return:
    """
    if outfile is None:
        prefix1= get_file_prefix(bedfile1)
        prefix2= get_file_prefix(bedfile2)
        location= get_file_location(bedfile1)

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


def pandas_summary(bed8file):
    """
    The bef8file is chr start end name *2 format
    :param bed8file:
    :return: the dict with (read1, read2): intersection
    """
    if count_file(bed8file)==0:
        return {}

    df=pandas.read_csv(bed8file, sep="\t", header=None, dtype={
        0:object, 1:int, 2:int,3:object,4:object,5:int,6:int,7:object
    })

    df.drop_duplicates()
    df["start_max"] = df[[1, 5]].max(axis=1)
    df["end_min"] = df[[2, 6]].min(axis=1)
    df["sub"] = df["end_min"] - df["start_max"]

    dfs = df[[3, 7, "sub"]]
    dfs.drop_duplicates(subset=[3,7])
    groupdfs = dfs.groupby([3, 7])
    aa = groupdfs.sum()

    intersection_dic=aa.to_dict()["sub"]

    return intersection_dic


def add_subread_bigg(bigg_raw):
    """
    add two bigglist together, with its subread count
    merge the reads containing in other reads' sub reads
    :param bigg_list1:
    :param bigg_list2:
    :return:
    """

    nameset=set()
    bigg_dic=OrderedDict()

    # merge the subreads with same name
    for bigg in bigg_raw:

        if bigg.name not in nameset:
            bigg_dic[bigg.name]=bigg
            nameset.add(bigg.name)
        else:
            subread_1=bigg_dic[bigg.name].subread
            subread_2=bigg.subread
            bigg_dic[bigg.name].subread=subread_1.union(subread_2)

    return list(bigg_dic.values())


def _merge_subread_bigg(bigg_raw):
    """
    sanity check for the read number
    :param bigg_raw:
    :return:
    """

    drop=set()
    names=[x.name for x in bigg_raw]
    subreads=set()
    for bigg in bigg_raw:
        subreads=subreads.union(bigg.subread)
    for name in names:
        if name in subreads:
            print(name)

    print((sorted(names)))

    #bigg_dic=bigg
    for bigg in bigg_raw:
        pass


def get_readall_bigg(bigg_list):
    """
    sum up the read names in the bigglist, include the sub read and names
    :param bigg_list:
    :return:
    """
    name_set=set()

    name_isoform=set([x.name for x in bigg_list])
    subread=set()
    for bigg in bigg_list:
        subread=subread.union(bigg.subread)

    name_set=name_isoform.union(subread)

    return name_set


def bigg_get_namedic(bigg_list):
    """
    from a bigg list with subreads
    get a isoform:readname dic
    :param bigg_list:
    :return:
    """
    name_dic={}
    for bigg in bigg_list:
        bigg.get_subread_from_str()
        if len(bigg.subread)>0:
            for name in bigg.subread:
                try:
                    name_dic[name].append(name)
                except KeyError:
                    name_dic[name]=[]
                    name_dic[name].append(name)

    return name_dic


def __get_namedic_unique_ratio(name_dic):
    iso_unique_count={}
    iso_all_count={}

    unique_name=set()
    dup_name=set()
    full_name=set()

    for name_iso, readname in list(name_dic.items()):
        full_name=full_name.union(set(readname))


def is_a_read(name, refname_set=None):
    """
    change to a method with not in set(refname_set)

    from the string of the name, judge if it is a read or a referece name
    Note: Now use a lot of "-" as the indicator of nanopore read, but this do not work for pacbio
    So add a indicattor of long read name, usually the reference name is not too long,
    the name of raw nanopore is 36 and raw pacbio is also the same

    For the reference gene name, the ensemble name is 15, gene predictor

    :param name:
    :return:
    """
    if refname_set is None:
        if len(name.split("-")) >= 5 or len(name) >= 24:
            return True
        else:
            return False
    else:
        if name in refname_set:
            return False
        else:
            return True


def bigg_count_write_native(isoform_list, gff_bed, out=None):
    """
    parser the output of cluster, get the count for each isoform
    :param isoform_list: the isoform_list
    :return:
    """
    # getthe refname_set:
    refname_set=set([bigg.name for bigg in gff_bed])

    # store sub-read name and number
    name_dic=OrderedDict()

    for bigg in isoform_list:
        bigg.get_subread_from_str()
        if len(bigg.subread)>0:
            for name in bigg.subread:
                if is_a_read(name, refname_set=refname_set): # judge if it is a read or a isoform
                    try:
                        name_dic[name]+=1
                    except KeyError:
                        name_dic[name]=1

    for bigg in isoform_list:
        coverage=0
        for name in bigg.subread:
            if is_a_read(name, refname_set=refname_set):
                coverage+=1.0/name_dic[name]

        bigg.coverage=coverage
        bigg.write_coverage()

    # debug
    #print name_dic
    write_bigg(isoform_list, out)


def __bigg_count_write_unique(bigg_list, out=None):
    """
    use the group information from the geneName2
    :param bigg_list:
    :param out:
    :return:
    """
    pass


def cat_bed(keyword):
    """
    locate all bed file with keyword in current wkdir, read them and return a full bed file
    used in the wkdir to cat all gene output together
    :param keyword: *_simple_coveragej.bed
    :return:
    """
    bigg_full=[]
    file_l=glob.glob(pathname=keyword, recursive=True)

    for file_one in file_l:
        try:
            bigg_one=read_bigg(file_one)
            bigg_full.extend(bigg_one)
        except Exception as e: # usually for the novel gene, the region may actually contain no track
           print (file_one, e)

    print("total number of isoform in all bed files", len(bigg_full))
    return bigg_full


def junction_pre(bigg_list, bigg_ref):
    """
    make a sanity filter, rm the bigg from different chro and reverse strand
    :param bigg_list:
    :param bigg_ref:
    :return:
    """

    bigg_new=[]
    chrom=bigg_ref[0].chrom
    strand=bigg_ref[0].strand

    for bigg in bigg_list:
        if bigg.chrom ==chrom and bigg.strand==strand:
            bigg_new.append(bigg)

    return bigg_new
