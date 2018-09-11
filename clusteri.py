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

from utils import parmap, set_tmp, wrapper_bedtools_intersect, wrapper_bedtools_intersection_muti


def flow_cluster(bigg_nano, bigg_gff, by="ratio_all", intronweight=0.5, core=40):

    bigg_nano.sort(key=operator.attrgetter("chromStart"))

    if by=="ratio_all":
        by1="ratio"
        by2="ratio_short"
    else:
        by1=by
        by2=by

    # hard code first filter of overalpping of 50 bp
    bigg_l1=prefilter_smallexon(bigg_nano, bigg_gff, cutoff=50, core=core) # using default cutoff 0.95
    bigg_list=bigg_l1+bigg_gff

    # can be change filters
    D1, bigg_list_by1=cal_distance(bigg_list, intronweight=intronweight, by=by1, core=core)
    _, bigg_l2=filter_D(D1, bigg_list_by1, by=by1) # using default cutoff 0.95

    D2, bigg_list_by2=cal_distance(bigg_l2, intronweight=intronweight, by=by2, core=core)
    D_remain, bigg_l3=filter_D(D2, bigg_list_by2, by=by2) # using default cutoff 0.95

    print(len(bigg_list), len(bigg_l2), len(bigg_l3))
    return D_remain, bigg_l3



def getij(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))
    return ij_list


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


def cal_distance(bigg_list, intronweight=0.5, by="ratio", core=40):
    """
    :param bigg_list:
    :param intronweight: if 0, do not cal the intron to save time
    :param by: used to cal the distance between two bigg object, can be "ratio", "ratio_short", "length", "length_short"
    :return: D: distance matrix
    """
    set_tmp()
    bigg_list.sort(key=operator.attrgetter("chromStart"))

    for i in bigg_list:
        i.to_bedfile()
        i.get_exon()

    length=len(bigg_list)
    D_exon=scipy.zeros([length, length])
    D_intron=scipy.zeros([length, length])

    # use the first one as gene_start
    gene_start= 0 #bigg_list[0].chromStart
    # thr order of the bigg_list is used as key to fetch the element further

    # get an pos combination and the name of bigg for each i
    ij_list=getij(bigg_list)
    pos_dic={}
    for n, bigg in enumerate(bigg_list):
        bigg.get_exon()
        bigg.to_bedfile()
        pos_dic[bigg.name]=n

    #par_list=[]
    #for n, bigg_n in enumerate(bigg_list):
    #    if n<len(bigg_list):
    #        bigg_rest=bigg_list[n+1:]
    #        par_list.append(bigg_n, bigg_rest)

    exon_l=wrapper_bedtools_intersection_muti(bigg_list, bigg_list, use="exon", core=40)
    intron_l=wrapper_bedtools_intersection_muti(bigg_list, bigg_list, use="intron", core=40)

    for one in exon_l:
        for two in one:
            name_bed1, name_bed2, intersection=two
            i=pos_dic[name_bed1]
            j=pos_dic[name_bed2]
            min_length = min(bigg_list[i].exonlen, bigg_list[j].exonlen)
            union=bigg_list[i].exonlen+bigg_list[j].exonlen-intersection

            if by == "ratio":
                # intron could be 0
                if min_length == 0:
                    D_exon[i,j]=0
                else:
                    similar = float(intersection) / union
                    D_exon[i,j]=1 - similar

            elif by == "ratio_short":
                # intron could be 0
                if min_length == 0:
                    D_exon[i,j]=0
                else:
                    D_exon[i,j]=1 - float(intersection) / min_length

        for one in intron_l:
            for two in one:
                name_bed1, name_bed2, intersection = two
                i = pos_dic[name_bed1]
                j = pos_dic[name_bed2]
                min_length = min(bigg_list[i].intronlen, bigg_list[j].intronlen)
                union = bigg_list[i].intronlen + bigg_list[j].intronlen - intersection

                if by == "ratio":
                    # intron could be 0
                    if min_length == 0:
                        D_intron[i, j] = 0
                    else:
                        similar = float(intersection) / union
                        D_intron[i, j] = 1 - similar

                elif by == "ratio_short":
                    # intron could be 0
                    if min_length == 0:
                        D_intron[i, j] = 0
                    else:
                        D_intron[i, j] = 1 - float(intersection) / min_length

    D=(D_exon+intronweight*D_intron)/float(1+intronweight)

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


def prefilter_smallexon(bigg_list,bigg_list_gff, cutoff=30, core=40):
    """
    only accept the D and bigg_list using exon only
    cal_distance(bigg_list, filter=False,intronweight=0,by="ratio", core=40, outfile="./test/d.csv"):

    :param D:
    :param bigg_list:
    :param cutoff: at least 50bp intersction with current annotation
    :return: keep
    """
    drop=set()
    fullset=set(range(len(bigg_list)))

    out_l=wrapper_bedtools_intersection_muti(bigg_list, bigg_list_gff, use="exon", core=core)

    for n, intersection_list in enumerate(out_l):
        numbers=[x[2] for x in intersection_list]
        if max(numbers)<=cutoff:
            drop.add(n)

    keep=sorted(list(fullset-drop))

    bigg_list_new=select_list(bigg_list, keep)
    return bigg_list_new



def filter_D(D, bigg_list, by="ratio", cutoff="auto"):

    """
    cutoff selection:
    learn from unc52, <0.025 in ratio_short
    return: index of the matrix that can be retained
    """
    if cutoff=="auto":
        if by=="ratio" or by=="ratio_short":
            cutoff=0.01
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


