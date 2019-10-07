#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 10/2/2019 1:07 PM
# @Author  : Runsheng     
# @File    : clusterj.py
"""
Make the cluster using junction information instead of the intersection of the exon and intron
The speed is expected to be much faster
Could be useful for high accuate long reads
"""

from tracklist import list_to_dic


def get_junction_dic(bigg_list, ref_weight=5):
    """
    get the dic showing the site: frequency
    site is as 'chro_site_strand'
    :param bigg_list: contain both ref and read in bigg format
    :param bigg_reads:
    :param ref_weight: the folder of the weight in
    :return: tuple(chro, pos, strand): coverage
    """
    site_dic = {}
    for bigg in bigg_list:
        if bigg.ttype=="nanopore_read":
            track_weight=1
        elif bigg.ttype=="isoform_anno":
            track_weight=1*ref_weight
        else:
            track_weight=1
        chro=bigg.chrom
        strand=bigg.strand

        bigg.get_junction()
        for pos in bigg.junction:
            try:
                site_dic[tuple((chro, pos, strand))]+=track_weight
            except KeyError:
                site_dic[tuple((chro, pos, strand))]=track_weight

    return site_dic


def boundary_correct(bigg_isoform, site_dic, coverage_cutoff=2):
    """
    Use the high confident junctions to correcnt the low frenquent junctions

    """
    strand=bigg_isoform[0].strand
    read_dic=list_to_dic(bigg_raed)

    bigg_isoform.get_exon()
    bigg_isoform.get_junction()
    bigg_isoform.get_coverage_from_str()
    bigg_isoform.get_subread_from_str()

    sub_list=[]
    for name in bigg_isoform.subread:
        read=read_dic[name]
        read.get_junction()
        sub_list.append(read)

        ## assume that the read have the same 3' as the isoform
        ## map the last ones or the isoforms with equal length to correct the boundary

    # need to correct the start, end and junctions separately



def get_junction_dis_tab(junction_dic, coverage_cutoff=5):
    pass



def get_start_end_dic(bigg_list, type="start", ref_weight=1):
    """
    :param type: the 5'start or the 3'end of the transcript,
    :param bigg_list:
    :return: tuple(chro, pos, strand): coverage
    """
    site_dic={}
    for bigg in bigg_list:
        if bigg.ttype=="isoform_anno":
            track_weight=ref_weight
        else:
            track_weight=1

        chro=bigg.chrom
        strand=bigg.strand

        if type=="start":
            pos=bigg.chromStart if strand=="+" else bigg.chromEnd
        else:
            pos=bigg.chromStart if strand=="-" else bigg.chromEnd

        try:
            site_dic[tuple((chro, pos, strand))]+=track_weight
        except KeyError:
            site_dic[tuple((chro, pos, strand))]=track_weight

    return site_dic


def chain_site(sites_l, offset=5):
    """
    use the sorted sites_l?
    merge the sites within 5nt to a range contaning all the fields
    """

    # merge the identical one first
    sites_l = list(set(sites_l))
    sites_sorted = sorted(sites_l, key=lambda x: (x[0], x[1], x[2]))

    site_range = []

    # init the value out of the loop
    line_p = sites_sorted[0]
    chro_p, start_p, strand_p = line_p

    # init range one
    range_one = [chro_p, start_p, start_p + 1, strand_p]

    for n, line in enumerate(sites_sorted):
        if n == 0:
            pass
        else:
            line_p = sites_sorted[n - 1]
            chro_p, start_p, strand_p = line_p

            chro, start, strand = line

            if chro_p == chro and strand_p == strand and 0 < abs(start - start_p) <= offset:
                range_one[2] = start + 1
                # print range_one
                # break

            else:
                site_range.append(range_one)
                # re_init range_one
                range_one = [chro, start, start + 1, strand]

    return site_range





def simple_merge(bigg_list, intron_reteintion=True):
    """
    merge the bigg with same junction to ref list
    merge the bigg with the 5' missing junction to ref list
    leave a switch to capture the intron retention or not
    if no, merge the intron retention isoform to the nearby one
    if yes, add these back to the ref list

    and the sub-read is divided into two parts: unique and un-unique (function added to post)

    after the clustering, add the most freq used 5' and 3' to the isoform as a representitive


    :param bigg_list:
    :return:
    """
    pass

