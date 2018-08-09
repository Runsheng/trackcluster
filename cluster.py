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
from ncls import NCLS
import scipy
import pylab
import scipy.cluster.hierarchy as sch

from utils import parmap

# NCLS is a quick method to find pyrange intersection


def cal_distance(bigg_list, core=40):
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
    length=len(bigg_list)
    D = scipy.zeros([length, length])
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
        distance=bigg_list[i].bedtool_cal_distance(bigg_list[j], gene_start)
        return distance

    dis_l=parmap(run_one, ij_list, core)

    for n, distance in zip(ij_list, dis_l):
        i,j=n
        D[i,j]=distance
        D[j,i]=distance

    with open("./test/d.csv", "w") as fw:
        for i in D:
            str_l=[str(x) for x in i]
            fw.write(",".join(str_l))
            fw.write("\n")
    return D


def filter_bigg(bigg_list):
    pass




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


