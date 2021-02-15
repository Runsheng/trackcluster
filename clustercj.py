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
from clusterj import junction_pre
from cluster import getij, get_pos_dic, select_list, select_D
from tracklist import add_subread_bigg, get_readall_bigg, list_to_dic
from utils import group_site, log_summary


# std lib import
import operator
import numpy
import pandas as pd
import itertools

from collections import OrderedDict


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
    length is equal to key
    :param junction:
    :param keys:
    :return:
    """
    set_junction=set(junction)
    j_binary=[1 if i in set_junction else 0 for i in keys]
    return j_binary


def get_read_junction_D(bigg_list, site_cov_dic):
    """
    a dataframe with readname: junction_pos
    read1: 1,0,1,0,1 (has this junction or not)

    :param bigg_dic: the dic for biggs (key is readname)
    :param site_cov_dic: the dic for junction sites, only keys are used
    :return: a pandas df
    """
    keys=site_cov_dic.keys()
    read_binaryj_dic=OrderedDict()

    for bigg in bigg_list:
        j_binary=generate_binary_junction_list(bigg.junction, keys)
        read_binaryj_dic[bigg.name]=j_binary

    df=pd.DataFrame.from_dict(read_binaryj_dic, orient="index",columns=keys)
    #df.index
    return df


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

    keys=site_cov_dic.keys()
    for n, k in enumerate(keys):
        if site_cov_dic[k] >= cov_cutoff:
            high_j[n]=k
        else:
            low_j[n]=k
    #print high_j, low_j

    low_j_group=group_site(low_j.keys())
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
        except KeyError:
            junction_new.append(junction)
    #print len(junction_new)==len(bigg.junction)

    bigg.junction = junction_new
    bigg.write_junction_to_exon()
    bigg.exon_to_block()

    return bigg, flag


def is_junction_in(junction_l, set_no):
    for i in junction_l:
        if i in set_no:
            return True
    else:
        return False


def flow_junction_correct(bigg_list):
    site_cov_dic=get_junction_dic(bigg_list, ref_weight=5, read_weight=1) # parameters
    w_to_r, w_to_no = get_corrected_dic(site_cov_dic, cov_cutoff=2, pos_cutoff=10) # parameters

    set_no=set(w_to_no)
    bigg_rare_junction=[]
    bigg_correct=[]

    flag_count=0

    for bigg in bigg_list:
        if is_junction_in(bigg.junction, set_no):
            bigg_rare_junction.append(bigg)
        else: # the remains worth correction
            bigg_new, flag=get_bigg_correct(bigg, w_to_r)
            flag_count+=flag

    logger=log_summary()
    logger.info("Corrected track number is {}".format(flag_count))
    logger.info("Number of tracks with rare junction is {}".format(len(bigg_rare_junction)))

    return bigg_correct, bigg_rare_junction



def __group_nearby_site(site_list, interval=5):
    """
    get list for all site with nearrby
    :param site_list:
    :param interval:
    :return:
    """
    # group the lists to regions
    groups = [[y[1] for y in g] for k, g in itertools.groupby(enumerate(site_list), key=lambda x: abs(x[1] - x[0])<=interval ) ]

    return groups


def get_junction_D(bigg_list, site_cov_dic):
    """
    Get the dic for sites, in pandas df form

    :param bigg_list:
    :param site_cov_dic:
    :return:
    """

    pass

