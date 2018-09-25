#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:01 AM
# @Author  : Runsheng     
# @File    : cluster.py
"""
Generating similarity matrix using differ between gene models
make dendrogram for further plotting
The input is a list of
"""
# self import
from track import bigGenePred

# third part import
import scipy
import os
import operator
from tracklist import wrapper_bedtools_intersect2, bigglist_to_bedfile, pandas_summary
from utils import parmap, set_tmp


def flow_cluster(bigg_nano, bigg_gff, by="ratio_all", intronweight=0.5):

    bigg_nano.sort(key=operator.attrgetter("chromStart"))

    if by=="ratio_all":
        by1="ratio"
        by2="ratio_short"
    else:
        by1=by
        by2=by

    # hard code first filter of overalpping of 50 bp
    bigg_l1=prefilter_smallexon(bigg_nano, bigg_gff, cutoff=50) # using default cutoff 0.95
    bigg_list=bigg_l1+bigg_gff

    # can be change filters
    D1, bigg_list_by1=cal_distance(bigg_list, intronweight=intronweight, by=by1)
    _, bigg_l2=filter_D(D1, bigg_list_by1, by=by1) # using default cutoff 0.95

    D2, bigg_list_by2=cal_distance(bigg_l2, intronweight=intronweight, by=by2)
    D_remain, bigg_l3=filter_D(D2, bigg_list_by2, by=by2) # using default cutoff 0.95

    print(len(bigg_list), len(bigg_l2), len(bigg_l3))
    return D_remain, bigg_l3


def getij(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))
    return ij_list


def get_pos_dic(bigg_list):
    pos_dic={}
    for n, bigg in enumerate(bigg_list):
        pos_dic[bigg.name]=n
    return pos_dic


def select_list(bigg_list, keep):
    # re_order D and bigg_list
    bigg_list_new=[]
    for i in keep:
        bigg_list_new.append(bigg_list[i])
    return bigg_list_new


def select_D(D, keep):
    D=D[keep,:]
    D=D[:,keep]

    return D


def prefilter_smallexon(bigg_list,bigg_list_gff, cutoff=50
                        ):
    """
    remove two kind of reads:

    1. not same strand as in current annotation
    2. with less than cutoff intersection with current annotation

    :param bigg_list:
    :param cutoff: at least 50bp intersection with current annotation
    :return: retained list
    """
    if len(bigg_list_gff)==0:
        return bigg_list

    # filter 1
    strand=bigg_list_gff[0].strand
    #print strand
    bigg_list_strand=[x for x in bigg_list if x.strand==strand]

    # filter 2
    nano_exon, nano_intron=bigglist_to_bedfile(bigg_list_strand)
    gff_exon, gff_intron=bigglist_to_bedfile(bigg_list_gff)

    exonfile=wrapper_bedtools_intersect2(nano_exon, gff_exon)
    out_d=pandas_summary(exonfile)
    # debug
    #print len(out_d)
    #print out_d
    keep_name=set()
    for k, intersection in out_d.items():
        nano_name, gff_name=k
        if intersection > cutoff:
            keep_name.add(nano_name)

    bigg_list_new=[]
    for bigg in bigg_list:
        if bigg.name in keep_name:
            bigg_list_new.append(bigg)

    ### clean up
    try:
        os.remove(exonfile)
        os.remove(nano_exon)
        os.remove(nano_intron)
        os.remove(gff_exon)
        os.remove(gff_intron)
    except OSError:
        pass

    return bigg_list_new


def cal_distance(bigg_list, intronweight=0.5, by="ratio"):
    """
    :param bigg_list:
    :param intronweight: if 0, do not cal the intron to save time
    :param by: used to cal the distance between two bigg object, can be "ratio", "ratio_short", "length", "length_short"
    :return: D: distance matrix
    """
    #wkdir=set_tmp()
    #os.chdir(wkdir)

    bigg_list.sort(key=operator.attrgetter("chromStart"))

    for i in bigg_list:
        i.get_exon()
        i.to_bedstr()

    length=len(bigg_list)
    D_exon=scipy.zeros([length, length])
    D_intron=scipy.zeros([length, length])

    # get an pos combination and the name of bigg for each i
    # ij_list=getij(bigg_list)
    pos_dic=get_pos_dic(bigg_list)

    # flow begin
    file_exon, file_intron = bigglist_to_bedfile(bigg_list)

    exon_out=wrapper_bedtools_intersect2(file_exon, file_exon)
    exon_i=pandas_summary(exon_out)

    intron_out=wrapper_bedtools_intersect2(file_intron, file_intron)
    intron_i=pandas_summary(intron_out)

    for k, intersection in exon_i.items():
        name1, name2=k
        i=pos_dic[name1]
        j=pos_dic[name2]
        min_length = min(bigg_list[i].exonlen, bigg_list[j].exonlen)
        union = bigg_list[i].exonlen + bigg_list[j].exonlen - intersection
        # debug
        if union <=0:
            print "exon", name1, name2, bigg_list[i].exonlen,  bigg_list[j].exonlen, union, intersection

        if by == "ratio":
            # intron could be 0
            if min_length == 0:
                D_exon[i, j] = 0
            else:
                similar = float(intersection) / union
                D_exon[i, j] = 1 - similar

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_exon[i, j] = 0
            else:
                D_exon[i, j] = 1 - float(intersection) / min_length


    for k, intersection in intron_i.items():
        name1, name2 = k
        i = pos_dic[name1]
        j = pos_dic[name2]
        min_length = min(bigg_list[i].intronlen, bigg_list[j].intronlen)
        union = bigg_list[i].intronlen + bigg_list[j].intronlen - intersection
        # debug
        if union <=0:
            print name1, name2, bigg_list[i].intronlen,  bigg_list[j].intronlen, union, intersection
        if by == "ratio":
            # intron could be 0
            if min_length == 0:
                D_intron[i, j] = 0
            else:
                #print union
                similar = float(intersection) / union
                D_intron[i, j] = 1 - similar

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_intron[i, j] = 0
            else:
                D_intron[i, j] = 1 - float(intersection) / min_length


    D=(D_exon+intronweight*D_intron)/float(1+intronweight)

    # cleanup
    try:
        os.remove(exon_out)
        os.remove(intron_out)
    except OSError:
        pass

    # debug:
    #print("D_exon",D_exon)
    #print("D_intron", D_intron)
    #print("D",D)

    #cleanup(remove_all=True)

    return D, bigg_list


def write_D(D, bigg_list_new, outfile="./test/d.csv"):
    bigg_name=[x.name for x in bigg_list_new]

    if outfile is None:
        pass
    else:
        with open(outfile, "w") as fw:
            fw.write(",".join(bigg_name))
            fw.write("\n")
            for i in D:
                str_l=[str(x) for x in i]
                fw.write(",".join(str_l))
                fw.write("\n")



def filter_D(D, bigg_list, by="ratio", cutoff="auto"):

    """
    cutoff selection:
    learn from unc52, <0.025 in ratio_short
    return: index of the matrix that can be retained
    """
    if cutoff=="auto":
        if by=="ratio" or by=="ratio_short":
            cutoff=0.05
        if by=="length" or by=="length_short":
            cutoff=100

    # hard code a cutoff for sw score of SL
    sw_score=11

    fullset=set(range(len(D)))
    drop=set()

    for bigg in bigg_list:
        bigg.get_exon()


    # same list
    ij_list=getij(D)

    # add a filter function for D, if the distance between exon is 0, merge the small one with
    # unless need to parer the intronD and exonD separately, or else the filter should be outer function
    for i,j in ij_list:
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<=bigg_list[j].exonlen:
                        drop.add(i)
                        bigg_list[j].subread.append(bigg_list[i].name)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        drop.add(j)
                        bigg_list[i].subread.append(bigg_list[j].name)

                if by=="ratio_short":
                    if bigg_list[i].exonlen<=bigg_list[j].exonlen and bigg_list[i].score<sw_score:
                        drop.add(i)
                        bigg_list[j].subread.append(bigg_list[i].name)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen and bigg_list[j].score<sw_score:
                        drop.add(j)
                        bigg_list[i].subread.append(bigg_list[j].name)

    keep=fullset-drop
    # change the default score of gene, no need to add
    for n, bigg in enumerate(bigg_list):
        if bigg.ttype!="nanopore_read":
            keep.add(n)
    keep=sorted(list(keep))


    #### debug
    #print(len(fullset), len(drop), len(keep))
    #print(keep)

    # re_order D and bigg_list
    bigg_list_new=select_list(bigg_list, keep)
    D=select_D(D, keep)

    return D, bigg_list_new


