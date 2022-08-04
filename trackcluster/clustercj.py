#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2/13/2021 1:44 PM
# @Author  : Runsheng     
# @File    : clustercj.py

"""
Use the self-corrected junction mode to speed up the process of clusterj
Aims:
1. The design should contain the 5' leader score, together with the filter function
2. Use only one mode, cluster the full and merge the truncated tracks, no mode selection retained
3. full junction for a gene would be a union set
"""

# self import
from utils import group_site


# std lib import
import logging
import numpy
import itertools

from collections import OrderedDict
logger = logging.getLogger('summary')


def get_junction_dic(bigg_list, ref_weight=5, read_weight=1):
    """
    get the dic showing the site: frequency
    site is as 'site'
    the bigg_list has already filtered away the reads with different chro and strand

    :param bigg_list: contain both ref and read in bigg format
    :param bigg_reads:
    :param ref_weight: the folder of the weight in
    :return: tuple(chro, pos, strand): coverage
    """

    #bigg_list.sort(key=operator.attrgetter("chromStart"))
    site_dic = {} # used to store site_cov

    for bigg in bigg_list:
        if bigg.ttype=="nanopore_read":
            track_weight=1*read_weight
        elif bigg.ttype=="isoform_anno":
            track_weight=1*ref_weight
        else:
            track_weight=1

        bigg.get_junction()
        for pos in bigg.junction:
            try:
                site_dic[pos]+=track_weight
            except KeyError:
                site_dic[pos]=track_weight

    site_cov_dic=OrderedDict(sorted(site_dic.items()))
    return site_cov_dic


def generate_binary_junction_list(junction, keys):
    """
    genearte a [0,1,1,0 ...] array to indicate if the junction is inside the keys
    length is equal to keys
    :param junction:
    :param keys:
    :return:
    """
    set_junction=set(junction)
    j_binary=[1 if i in set_junction else 0 for i in keys]
    return j_binary


def get_read_junction_dic(bigg_list, site_cov_dic):
    """
    a dataframe with readname: junction_pos
    read1: 1,0,1,0,1 (has this junction or not)

    :param bigg_dic: the dic for biggs (key is readname)
    :param site_cov_dic: the dic for junction sites, only keys are used
    :return: a pandas df
    """
    keys=list(site_cov_dic.keys())
    read_binaryj_dic=OrderedDict()

    for bigg in bigg_list:
        j_binary=generate_binary_junction_list(bigg.junction, keys)
        read_binaryj_dic[bigg.name]=numpy.array(j_binary)

    #df=pd.DataFrame.from_dict(read_binaryj_dic, orient="index",columns=keys)
    #df.index
    return read_binaryj_dic


def get_array_freq(j_del):
    """
    get the count for element number of a array
    designed for [-1,0,1] array for the junction compare (maybe others)
    :param j:
    :return: freq_dic []
    """
    nums, counts = numpy.unique(j_del, return_counts=True)
    freq_dic_raw = dict(list(zip(nums, counts)))

    # count 1 if 1 is not there
    freq_dic = {}
    for k in [-1, 0, 1]:
        try:
            freq_dic[k] = freq_dic_raw[k]
        except KeyError:
            freq_dic[k] = 0
    return freq_dic


def get_pos(junction_array, key=1):
    """
    used inside other function
    :param junction_array:
    :param key:
    :return:
    """
    pos_l=[]
    for pos, value in enumerate(junction_array):
        if value==key:
            pos_l.append(pos)
    print(pos_l)
    return pos_l


def split_site(junction_array):
    """
    [0,0,0,1,1,1,0,0,0]
    to [[0,0,0], [1,1,1], [0,0,0]]
    :param junction_array: [0,1,0] array indicate if the junction is there
    :return:
    """
    # sanity check:
    if len(junction_array)%2!=0:
        pass
        # could have some exception: 1. include the 0 or the -1 junction 2. intron retain coupled with missed exon
        # not bug, # print "Error in the number of orders !", missed_order

    # group the lists to sublist with equal values
    groups = [list(g) for k, g in itertools.groupby(junction_array)]

    return groups


def _is_junction_inside(j1, j2):
    """
    test if j1 is in j2， change into comp all junctions except single exon
    still very slow， unused
    :param j1: junction numpy array, with all aligned part
    :param j2:
    :return: 0 for contain, 3 for differ, 1 for j1 in j2, -1 for j2 in j1, "e" for error
    """
    # do not consider the single exon track/read
    # first to simplify the two unions from the large union

    if numpy.array_equal(j1,j2):
        return 0

    # select the union of the junction list
    pos_simple=[pos for pos, value in enumerate(j1+j2) if value>=1]
    j1_s=j1[pos_simple]
    j2_s=j2[pos_simple]

    #print(j1_s)
    #print(j2_s)

    if numpy.sum(j1_s)==0 or numpy.sum(j2_s)==0:
        return "is_junction_inside error: single exon"

    j12_del=j1_s-j2_s
    #freq_dic = get_array_freq(j12_del)

    # filter the easy differ ones

    #if freq_dic[-1] > 0 and freq_dic[1] > 0:  # has more and less junction, must be 3
    #    return 3

    # consider 3 types using group_l, this will be slower
    group_l=split_site(j12_del)
    #print(group_l)

    # did not consider the len=0 and len=1
    # single exon and equal junction should have been excluded
    if len(group_l)>=3:
        return 3
    elif len(group_l)==2: # should be all 0, 1 or all 0 ,-1
        # need to consider 5' or 3' missing
        # 5' is belong, 3' is differ
        a,b=group_l
        if -1<a[0]+b[0]<1: # should not happen for [-1,-1], [1,1], but no overalp
            return 3
        if a[0]+ b[0]>=1: #j1_s > j2_s
            if b[0]==0: # 5' missing
                return -1
            else: # 3' missing
                return 3

        elif a[0]+ b[0]<=-1: #j1_s < j2_s
            if b[0]==0: # 5' missing
                return 1
            else: # 3' missing
                return 3


def __compare_junction(j1, j2):
    """
    Two junctions in j1, j2, generated by get_read_junction_dic
    :param j1: numpy array, contian 0 or 1 only, the position is well aligned for all tracks,
    :param j2:
    :return: use 0:"equal", 1:"in", 2:"contain", 3:"diff"
    does not conserider single exon tracks
    """
    # do not consider the single exon track/read

    # single junction first
    pass




def get_corrected_dic(site_cov_dic, cov_cutoff=2, pos_cutoff=10):
    """
    get a dic [wrong_site: corrected_site], and use this for junction correction
    the reads contains the sites can fetched from the df
    :param site_cov_dic:
    :param cov_cutoff:
    :return: a dic used to correct the wrong junction, and the single junctions which can not be corrected
    """
    high_j=OrderedDict()
    low_j=OrderedDict()

    keys=list(site_cov_dic.keys())
    for n, k in enumerate(keys):
        if site_cov_dic[k] >= cov_cutoff:
            high_j[n]=k
        else:
            low_j[n]=k
    #print high_j, low_j

    low_j_group=group_site(list(low_j.keys()))
    #print low_j_group

    w_to_r=OrderedDict() # site can be corrected
    w_to_no=[] # site can not be corrected

    for group in low_j_group:
        start=group[0]
        end=group[-1]
        for site in group: # might be only1 or more
            try:
                # in rare cases, the previous and latter junction can both be used to correct
                # the junction (exon <20), but have not encountered this, so ignore
                if abs(low_j[site] - high_j[end+ 1]) <= pos_cutoff:
                    w_to_r[low_j[site]] = high_j[end + 1]
                elif abs(low_j[site] - high_j[start - 1]) <= pos_cutoff:
                    w_to_r[low_j[site]] = high_j[start - 1]
                else:
                    w_to_no.append(low_j[site])
            except KeyError:
                w_to_no.append(low_j[site])
    #print len(w_to_r), w_to_r
    #print "==========="
    #print len(w_to_no)
    return w_to_r, w_to_no


def get_bigg_correct(bigg, w_to_r):
    """
    bigg with junction (already run the bigg.get_junction)
    :param bigg:
    :param w_to_r:
    :return: bigg_corrected and flag (1 is corrected, 0 is untouched)
    """
    junction_new=[]
    flag=0

    if len(bigg.junction) == 0 or bigg.ttype == "isoform_anno":
        return bigg, flag

    for junction in bigg.junction:
        try:
            junction_new.append(w_to_r[junction])
            flag=1
            bigg.ttype="nanopore_read_corrected"
        except KeyError:
            junction_new.append(junction)
    #print( len(junction_new)==len(bigg.junction) )

    bigg.junction = junction_new
    bigg.write_junction_to_exon()
    bigg.exon_to_block()

    return bigg, flag


def is_junction_in_set(junction_l, set_no):
    for i in junction_l:
        if i in set_no:
            return True
    else:
        return False


def flow_junction_correct(bigg_list, cov_cutoff=2, pos_cutoff=10):
    """

    :param bigg_list:
    :param cov_cutoff:
    :param pos_cutoff: only merge sites within the distance of pos_cutoff
    :return:
    """
    site_cov_dic=get_junction_dic(bigg_list, ref_weight=5, read_weight=1) # parameters
    w_to_r, w_to_no = get_corrected_dic(site_cov_dic, cov_cutoff=cov_cutoff, pos_cutoff=pos_cutoff) # parameters

    set_no=set(w_to_no)
    bigg_rare_junction=[]
    bigg_correct=[]

    flag_count=0

    for bigg in bigg_list:
        if is_junction_in_set(bigg.junction, set_no):
            bigg_rare_junction.append(bigg)
        else: # the remains worth correction
            bigg_new, flag=get_bigg_correct(bigg, w_to_r)
            bigg_correct.append(bigg_new)
            flag_count+=flag

    logger.debug("Corrected track number is {}".format(flag_count))
    logger.debug("Number of tracks with rare junction is {}".format(len(bigg_rare_junction)))

    return bigg_correct, bigg_rare_junction


# use junction_simple_merge
