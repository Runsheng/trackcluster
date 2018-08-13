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
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist
import operator

def plot_tree(bigg_list, out="test.pdf", pos=None):

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

    Z = sch.linkage(D, method='single')
    dend= sch.dendrogram(Z)

    # ----------------------------------------
    # get leave coordinates, which are at y == 0

    def flatten(l):
        return [item for sublist in l for item in sublist]

    X = flatten(dend['icoord'])
    Y = flatten(dend['dcoord'])
    leave_coords = [(x, y) for x, y in zip(X, Y) if y == 0]

    # in the dendogram data structure,
    # leave ids are listed in ascending order according to their x-coordinate
    order = np.argsort([x for x, y in leave_coords])
    id_to_coord = dict(zip(dend['leaves'], [leave_coords[idx] for idx in order]))  # <- main data structure

    # map endpoint of each link to coordinates of parent node
    children_to_parent_coords = dict()
    for i, d in zip(dend['icoord'], dend['dcoord']):
        x = (i[1] + i[2]) / 2
        y = d[1]  # or d[2]
        parent_coord = (x, y)
        left_coord = (i[0], d[0])
        right_coord = (i[-1], d[-1])
        children_to_parent_coords[(left_coord, right_coord)] = parent_coord

    # traverse tree from leaves upwards and populate mapping ID -> (x,y)
    root_node, node_list = to_tree(Z, rd=True)
    ids_left = range(len(dend['leaves']), len(node_list))

    while len(ids_left) > 0:

        for ii, node_id in enumerate(ids_left):
            node = node_list[node_id]
            if (node.left.id in id_to_coord) and (node.right.id in id_to_coord):
                left_coord = id_to_coord[node.left.id]
                right_coord = id_to_coord[node.right.id]
                id_to_coord[node_id] = children_to_parent_coords[(left_coord, right_coord)]

        ids_left = [node_id for node_id in range(len(node_list)) if not node_id in id_to_coord]

    # plot result on top of dendrogram
    ax = plt.gca()
    for node_id, (x, y) in id_to_coord.iteritems():
        if not node_list[node_id].is_leaf():
            ax.plot(y, x, 'ro')
            ax.annotate(str(node_id), (y, x), xytext=(0, -8),
                        textcoords='offset points',
                        va='top', ha='center')
        # add the lines to the x
        else:
            pass


    dend['node_id_to_coord'] = id_to_coord

    plt.savefig(out)

    return dend


def line_plot(bigg_list, out="test.pdf"):
    """
    :param bigg_list: a list contain multiple bigGenePred class
    :return:
    """
    bigg_list.sort(key=operator.attrgetter("chromStart"))

    D=cal_distance(bigg_list)

    # init a figure size
    bigg=bigg_list[0]

    start=bigg.chromStart
    bigg.get_exon(start)
    width=(bigg.chromEnd-start)/400.0 # 100 bp=1 inch width
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
    for bigg in Z1["color_list"]:
        color_xy_list.extend(bigg * 4)

    leave_cords = [x for x,y in zip(x_list, y_list) if y == 0]
    leave_color=[c for y, c in zip(y_list, color_xy_list) if y==0]

    # in the dendogram data structure,
    # leave ids are listed in ascending order according to their x-coordinate
    order = np.argsort([x for x in leave_cords])
    id_cord = [leave_cords[idx] for idx in order]  # <- main data structure
    id_color = [leave_color[idx] for idx in order]


    #id_to_coord = dict(zip(Z1['leaves'], [leave_coords[idx] for idx in order]))

    #ax1.set_xticks([])
    ax1.set_yticks([])

    # second, the track for gene models

    # get the order of the gene models in dendrogram
    bigg_list_new=[]
    for x in Z1["leaves"]:
        bigg_list_new.append(bigg_list[x])

    ax2=fig.add_axes([0.225,0.09,0.65,0.91]) # adjust the ax to fit figure
    for n,bigg in enumerate(bigg_list_new):
        bigg.get_exon(start)
        for exon in bigg.exon:
            x_start, x_end = exon
            # debug
            #print n, x_start, x_end
            if bigg.ttype == "nanopore_read":
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=4, colors=id_color[n])
            else:
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=4, colors="#756bb1")

        for intron in bigg.intron:
            x_start, x_end = intron
            if bigg.ttype == "nanopore_read":
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=1, colors=id_color[n], linestyles="dotted")
            else:
                plt.hlines(y=id_cord[n],  xmin=x_start, xmax=x_end, linewidth=1, colors=id_color[n])


    #plt.ylim(ymin=-5)
    plt.margins(0, 0)
    ax2.set_yticks([])
    plt.savefig(out)
    #plt.show()
