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
    height=(len(bigg_list))/5.0+2
    fig=plt.figure(figsize=(width,height))

    print(width, height)

    #### first the dendrogram
    ax1 = fig.add_axes([0.02, 0.09, 0.2, 0.85]) # set region

    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, orientation='left')

    #ax1.set_xticks([])
    ax1.set_yticks([])

    # second, the track for gene models

    # get the order of the gene models in dendrogram
    bigg_list_new=[]
    for x in Z1["leaves"]:
        bigg_list_new.append(bigg_list[x])

    # get the real colour for each leaves (need to +1)
    color_l=[]
    dcorrd_l=[]
    for n, line in enumerate(zip(Z1["color_list"], Z1["icoord"])):
        x, y=line
        if n==0:
            pass
        else:
            if x==(Z1['color_list'])[n-1]:
                pass
            else:
                color_l.append((Z1['color_list'])[n-1])
                dcorrd_l.append(Z1["icoord"][n-1][1])
        color_l.append(x)
        dcorrd_l.append(y[1])

    dcorrd_l_new=[]
    for x in Z1["leaves"]:
        dcorrd_l_new.append(dcorrd_l[x])

    print(len(bigg_list_new), len(color_l), len(dcorrd_l))

    ax2=fig.add_axes([0.25,0.09,0.6,0.85]) # adjust the ax to fit figure
    for n,i in enumerate(bigg_list_new):
        i.get_exon(start)
        for exon in i.exon:
            x_start, x_end = exon
            # debug
            #print n, x_start, x_end
            l=plt.hlines(y=n+0.8,  xmin=x_start, xmax=x_end, linewidth=4, colors=color_l[n])

        for intron in i.intron:
            x_start, x_end = intron
            l=plt.hlines(y=n+0.8,  xmin=x_start, xmax=x_end, linewidth=1, colors=color_l[n])


    ax2.set_yticks([])
    plt.savefig(out)
    #plt.show()
