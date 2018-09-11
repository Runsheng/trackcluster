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
from clusteri import flow_cluster, write_D
from tracklist import write_bigg

from copy import deepcopy
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist
import operator

# static
# set colour for the dendrogram
brew_11 = ["#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd","#5e4fa2"]



def line_plot_merge(bigg_nano,
                    bigg_gff,
                    by="ratio_all",
                    intronweight=0.5,core=40,
                    out="./test/test.pdf",
                    biggout="./test/test.bed",
                    Dout="./test/d.csv",
                    color=None):
    """
    :param bigg_list: a list contain multiple bigGenePred class
    :return:
    """
    D,bigg_list_by2=flow_cluster(bigg_nano,bigg_gff, by, intronweight=intronweight, core=core )
    write_D(D, bigg_list_by2, Dout)

    # init a figure size
    bigg=bigg_gff[0]

    start=bigg.chromStart
    bigg.get_exon(0)
    width=(bigg.chromEnd-start)/100.0 # 100 bp=1 inch width
    height=(len(bigg_list_by2))/5.0+2  # 1 gene 0.2 inch height
    fig=plt.figure(figsize=(width,height))

    print(width, height)

    #### first the dendrogram
    ax1 = fig.add_axes([0.02, 0.09, 0.2, 0.91]) # set region
    Y = sch.linkage(D, method='single')


    if color is None:
        color=brew_11
    sch.set_link_color_palette(color)

    Z1 = sch.dendrogram(Y, orientation='left', color_threshold=0.5, above_threshold_color='#8c510a')

    #ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.margins(0,0)


    # second, the track for gene models
    ### to plot the gene
    def flatten(l):
        return [item for sublist in l for item in sublist]

    x_list = flatten(Z1['icoord'])
    y_list = flatten(Z1['dcoord'])

    color_xy_list = []
    for bigg_color in Z1["color_list"]:
        color_xy_list.append(bigg_color)
        color_xy_list.append(bigg_color)
        color_xy_list.append(bigg_color)
        color_xy_list.append(bigg_color)


    leave_cords = [x for x,y in zip(x_list, y_list) if y == 0]
    leave_color=[c for y, c in zip(y_list, color_xy_list) if y==0]
    # in the dendogram data structure,
    # leave ids are listed in ascending order according to their x-coordinate
    order = np.argsort([x for x in leave_cords])
    id_cord = [leave_cords[idx] for idx in order]  # <- main data structure
    id_color = [leave_color[idx] for idx in order]


    # get the order of the gene models in dendrogram
    bigg_list_new=[]
    for x in Z1["leaves"]:
        bigg_list_new.append(bigg_list_by2[x])

    ax2=fig.add_axes([0.225,0.09,0.65,0.91]) # adjust the ax to fit figure
    ax2.set_yticks([])
    ax2.set_ylim(min(id_cord)-1, max(id_cord)+1)
    ax2.margins(0, 0)
    #ax2.set_ylim(top=ymax+5)


    for n,bigg in enumerate(bigg_list_new):
        bigg.get_exon(start)
        for exon in bigg.exon:
            x_start, x_end = exon
            # debug
            #print n, x_start, x_end
            if bigg.ttype == "nanopore_read":
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=4, colors=id_color[n])
            else:
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=5, colors="black")

        for intron in bigg.intron:
            x_start, x_end = intron
            if bigg.ttype == "nanopore_read":
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=1, colors=id_color[n], linestyles="dotted")
            else:
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=1, colors=id_color[n])

    ax2.margins(0, 0)
    ax2.set_yticks([])

    plt.savefig(out)

    ### save nessary files
    for bigg in bigg_list_new:
        bigg.write_subread()
    write_bigg(bigg_list_new, out=biggout)
