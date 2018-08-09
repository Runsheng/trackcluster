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

# NCLS is a quick method to find pyrange intersection


def cal_distance(bigg_list):
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
    for i, bigg_one in enumerate(bigg_list):
        for j, bigg_two in enumerate(bigg_list):
            D[i, j]=bigg_one.bedtool_cal_distance(bigg_two, gene_start)

    with open("d.csv", "w") as fw:
        for i in D:
            str_l=[str(x) for x in i]
            fw.write(",".join(str_l))
            fw.write("\n")
    return D


def plot_cluster(D):
    fig = pylab.figure(figsize=(8, 8))
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='left')
    print Z1["leaves"]
    fig.savefig("dend.pdf")


def cal_distance2(bigg_list):
    """
    :param bigg_list:
    :return:
    """
    for n, i in enumerate(bigg_list):
        if n==0: # use fist one to align the gene tracks
            start=i.chromStart
    pass