#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/8/2018 3:17 PM
# @Author  : Runsheng     
# @File    : plots.py

"""
Plotting function for gene track
"""
# self import
from track import bigGenePred
# third part import
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from cluster import cal_distance

def line_plot(bigg_list, color=None, out="test.pdf"):
    """

    :param bigg_list: a list contain multiple bigGenePred class
    :return:
    """

    D=cal_distance(bigg_list)

    # init a figure size
    i=bigg_list[0]

    start=i.chromStart
    i.get_exon(start)
    width=(i.chromEnd-start)/400.0 # 100 bp=1 inch width
    height=(len(bigg_list))/5.0
    fig=plt.figure(figsize=(width,height))

    print(width, height)

    #### first the dendrogram
    ax1 = fig.add_axes([0.09, 0.09, 0.2, 0.90]) # set region

    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='left')


    # second, the track for gene models

    bigg_list_new=[]

    for x in Z1["leaves"]:
        bigg_list_new.append(bigg_list[x])

    ax2=fig.add_axes([0.3,0.09,0.7,0.90])
    for n, i in enumerate(bigg_list_new):
        for exon in i.exon:
            x_start, x_end = exon
            # debug
            #print n, x_start, x_end
            l=plt.hlines(y=n,  xmin=x_start, xmax=x_end, linewidth=8)

        for intron in i.intron:
            x_start, x_end = intron
            l=plt.hlines(y=n,  xmin=x_start, xmax=x_end, linewidth=1)


    plt.axis("off")
    plt.savefig(out)
    #plt.show()
