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
import pylab
import scipy.cluster.hierarchy as sch
import operator

from utils import parmap, set_tmp


def flow_cluster(bigg_list, by="ratio_all", intronweight=0.5, core=40):

    bigg_list.sort(key=operator.attrgetter("chromStart"))

    if by=="ratio_all":
        by1="ratio"
        by2="ratio_short"
    else:
        by1=by
        by2=by

    # hard code first filter for intron ll
    D0=cal_distance(bigg_list,intronweight=0, by="ratio", core=core)
    _, bigg_l1=prefilter_smallexon(D0, bigg_list) # using default cutoff 0.95

    # can be change filters
    D1=cal_distance(bigg_l1, intronweight=intronweight, by=by1, core=core)
    _, bigg_l2=filter_D(D1, bigg_l1, by=by1) # using default cutoff 0.95

    D2=cal_distance(bigg_l2, intronweight=intronweight, by=by2, core=core)
    D_remain, bigg_l3=filter_D(D2, bigg_l2, by=by2) # using default cutoff 0.95

    return D_remain, bigg_l3



def getij(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))
    return ij_list


def select_list(D, bigg_list, keep):
    # re_order D and bigg_list
    D=D[keep,:]
    D=D[:,keep]
    bigg_list_new=[]
    for i in keep:
        bigg_list_new.append(bigg_list[i])

    return D, bigg_list_new


def cal_distance(bigg_list, intronweight=0.5, by="ratio_short", core=40):
    """
    speed record:
    50: 177 s
    10: 0.452 s
    after update:
    50: 6s
    10: 0.2s
    100: 65S
    323:

    after bedtools:
    5:0.237
    10: 1.340
    50: 26S
    100: 97S
    323: 938s

    :param bigg_list:
    :param intronweight: if 0, do not cal the intron to save time
    :param by: used to cal the distance between two bigg object, can be "ratio", "ratio_short", "length", "length_short"
    :return: D: distance matrix
    """
    set_tmp()

    length=len(bigg_list)
    D_exon=scipy.zeros([length, length])
    D_intron=scipy.zeros([length, length])

    # use the first one as gene_start
    gene_start=bigg_list[0].chromStart
    # thr order of the bigg_list is used as key to fetch the element further

    # get an pos combination
    ij_list=getij(bigg_list)

    def run_one(n):
        i,j=n
        distance_exon=bigg_list[i].pyrange_cal_distance_exon(bigg_list[j], gene_start, by)
        if intronweight != 0:
            distance_intron=bigg_list[i].pyrange_cal_distance_intron(bigg_list[j], gene_start, by)
        else:
            distance_intron=0
        return (i,j,distance_exon, distance_intron) # need to store i,j as the order may change in the parmap

    dis_l=parmap(run_one, ij_list, core)

    for pack in dis_l:
        i,j,distance_exon, distance_intron=pack
        D_exon[i,j]=distance_exon
        D_exon[j,i]=distance_exon
        D_intron[j,i]=distance_intron
        D_intron[i,j]=distance_intron

    D=(D_exon+intronweight*D_intron)/float(1+intronweight)

    # debug:
    #print("D_exon",D_exon)
    #print("D_intron", D_intron)
    #print("D",D)

    return D


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


def prefilter_smallexon(D, bigg_list, cutoff=0.95):
    """
    only accept the D and bigg_list using exon only
    cal_distance(bigg_list, filter=False,intronweight=0,by="ratio", core=40, outfile="./test/d.csv"):

    :param D:
    :param bigg_list:
    :return: keep
    """
    drop=set()
    fullset=set(range(len(D)))

    for i in range(len(D)):
        if D[i,].mean()>=cutoff:
            drop.add(i)

    keep=sorted(list(fullset-drop))

    D, bigg_list_new=select_list(D, bigg_list, keep)
    return D, bigg_list_new


def filter_D(D, bigg_list, by="ratio", cutoff="auto"):

    """
    cutoff selection:
    learn from unc52, <0.025 in ratio_short
    return: index of the matrix that can be retained
    """
    if cutoff=="auto":
        if by=="ratio" or by=="ratio_short":
            cutoff=0.025
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
    D, bigg_list_new=select_list(D, bigg_list, keep)

    return D, bigg_list_new


