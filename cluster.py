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

from utils import parmap, set_tmp
from pybedtools import cleanup


def cal_distance(bigg_list,filter=False,intronweight=0.5,by="ratio_short", core=40):
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
    :return:
    """
    set_tmp()

    length=len(bigg_list)
    D_exon=scipy.zeros([length, length])
    D_intron=scipy.zeros([length, length])

    # use the first one as gene_start
    gene_start=bigg_list[0].chromStart
    # thr order of the bigg_list is used as key to fetch the element further

    # get an pos combination
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))

    def run_one(n):
        i,j=n
        distance_exon=bigg_list[i].bedtool_cal_distance_exon(bigg_list[j], gene_start, by)
        if intronweight != 0:
            distance_intron=bigg_list[i].bedtool_cal_distance_intron(bigg_list[j], gene_start, by)
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
    #print(D_exon)
    #print(D_intron)
    #print(D)

    # add a filter function for D, if the distance between exon is 0, merge the small one with
    # unless need to parer the intronD and exonD seperatedly, or else the filter should be outer function
    keep=range(len(D))
    if filter:
        keep=filter_D(D, bigg_list, by, cutoff="auto")
        D=D[keep,:]
        D=D[:,keep]

    bigg_list_new=[]
    for i in keep:
        bigg_list_new.append(bigg_list[i])

    with open("./test/d.csv", "w") as fw:
        for i in D:
            str_l=[str(x) for x in i]
            fw.write(",".join(str_l))
            fw.write("\n")
    cleanup(remove_all=True)

    return D, bigg_list_new


def filter_D(D, bigg_list, by="ratio", cutoff="auto"):

    """
    may need to pre-filter the bigg_list using the current gff based gene model
    cutoff selection:
    intersection < 20bp?
    sum differ in exon and intron when using the "_short" method is < 10 bp

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
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))

    for i,j in ij_list:
        if D[i,j]<cutoff:
            # debug
            #print i,j

            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<=bigg_list[j].exonlen:
                        drop.add(i)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        drop.add(j)
                if by=="ratio_short":
                    if bigg_list[i].exonlen<=bigg_list[j].exonlen and bigg_list[i].score<sw_score:
                        drop.add(i)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen and bigg_list[j].score<sw_score:
                        drop.add(j)
    keep=sorted(list(fullset-drop))
    # debug
    print(len(fullset), len(drop), len(keep))
    print(keep)
    return keep



def __plot_cluster(D):
    """
    test only
    :param D:
    :return:
    """
    fig = pylab.figure(figsize=(8, 8))
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='left')
    print Z1["leaves"]
    fig.savefig("dend.pdf")


