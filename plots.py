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
import numpy as np
import scipy.cluster.hierarchy as sch
from cluster import cal_distance
from copy import deepcopy


def line_plot(bigg_list, color=None, out="test.pdf"):
    """
    :param bigg_list: a list contain multiple bigGenePred class
    :return:
    """
    bigg_list_bk=deepcopy(bigg_list)

    D=cal_distance(bigg_list)

    # init a figure size
    i=bigg_list[0]

    start=i.chromStart
    i.get_exon(start)
    width=(i.chromEnd-start)/400.0 # 100 bp=1 inch width
    height=(len(bigg_list))/5.0+2  # 1 gene 0.2 inch height
    fig=plt.figure(figsize=(width,height))

    print(width, height)

    #### first the dendrogram
    ax1 = fig.add_axes([0.02, 0.09, 0.2, 0.91]) # set region

    Y = sch.linkage(D, method='single')
    Z1 = sch.dendrogram(Y, orientation='left')

    def flatten(l):
        return [item for sublist in l for item in sublist]

    x_list = flatten(Z1['icoord'])
    y_list = flatten(Z1['dcoord'])

    color_xy_list = []
    for i in Z1["color_list"]:
        color_xy_list.extend(i * 4)

    leave_cords = [x for x,y in zip(x_list, y_list) if y == 0]
    leave_color=[c for y, c in zip(y_list, color_xy_list) if y==0]

    # in the dendogram data structure,
    # leave ids are listed in ascending order according to their x-coordinate
    order = np.argsort([x for x in leave_cords])
    id_cord = [leave_cords[idx] for idx in order]  # <- main data structure
    id_color = [leave_color[idx] for idx in order]

    ymax=max(id_cord)

    #id_to_coord = dict(zip(Z1['leaves'], [leave_coords[idx] for idx in order]))

    #ax1.set_xticks([])
    #ax1.set_yticks([])

    # second, the track for gene models

    # get the order of the gene models in dendrogram
    bigg_list_new=[]
    for x in Z1["leaves"]:
        bigg_list_new.append(bigg_list_bk[x])

    ax2=fig.add_axes([0.25,0.09,0.65,0.91]) # adjust the ax to fit figure
    for n,i in enumerate(bigg_list_new):
        i.get_exon(start)
        for exon in i.exon:
            x_start, x_end = exon
            # debug
            #print n, x_start, x_end
            plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=4, colors=id_color[n])

        for intron in i.intron:
            x_start, x_end = intron
            plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=1, colors=id_color[n], linestyles="dotted")

    plt.ylim(0, ymax)
    #ax2.set_yticks([])
    plt.savefig(out)
    #plt.show()
