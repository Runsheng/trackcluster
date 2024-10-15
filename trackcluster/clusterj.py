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
# std lib import

# third part import

# self import
from trackcluster.tracklist import junction_pre
from trackcluster.post import compare_ei_by_boudary, is_junction_equal
from utils import group_site
from trackcluster.cluster import select_list



def chose_read_from_list(name, bigg_list):
    for bigg in bigg_list:
        if bigg.name==name:
            return bigg
    return None


def __get_start_end_dic(bigg_list, type="start", ref_weight=1):
    """
    Todo： may need to chop the long isoform to shorter one by refine the 5' and 3'
    There are two levels of chopping， affect the junction or not
    The SL signal in nematode may work

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



def is_junction_5primer(junction):
    group_j=group_site(junction)
    if junction[0]==0 and len(group_j)==1:
        return True
    else:
        return False


def is_junction_inside(bigg1, bigg2):
    """
    merge the 5' missing but keep the 3' missing as a novel one

    :param bigg1:
    :param bigg2:
    :param ignore_5primer:
    :return:
    """
    bigg1.get_junction()

    # single exon reads, treat seperetedly
    if len(bigg1.junction)==0:
        return False

    if is_junction_equal(bigg1, bigg2, offset=0) is True:
        if bigg1.exonlen>=bigg2.exonlen:
            return False
        else:
            return True

    missed_order, extra_order=compare_ei_by_boudary(bigg1, bigg2, offset=0)

    if len(missed_order)==0:
        return False
    else:
        if is_junction_5primer(missed_order) and len(extra_order)==0:
            return True
        else:
            return False


def is_single_exon_in(bigg1, bigg2):
    """
    ignore the single exon reads with same 3'
    keep the single reads with end site before the last exon

    :param bigg1: single exon bigg
    :param bigg2: can be single exon and can be longer
    :return: boolean
    """
    bigg2.get_junction()
    # sanity check, rm the reverse and not same chro ones
    #if bigg1.chrom!=bigg2.chrom or bigg1.strand!=bigg2.strand:
    #    return False

    # single exon
    if len(bigg2.junction)==0:
        if bigg1.strand=="+":
            if bigg1.chromStart<=bigg2.chromStart or bigg1.chromEnd>=bigg2.chromEnd:
                return False
            else:
                return True
        else: # strand "-"
            if bigg1.chromStart>=bigg2.chromStart or bigg1.chromEnd<=bigg2.chromEnd:
                return False
            else:
                return True
    # reads with junction, judge if it is the 3' degradation
    # not single exon, use the last junction
    else:
        last_junction=bigg2.junction[-1]

        if bigg1.strand=="+":
            if bigg1.chromStart<=last_junction or bigg1.chromEnd>=bigg2.chromEnd:
                return False
            else:
                return True
        else: #"-"
            if bigg1.chromStart<=bigg2.chromStart or bigg1.chromEnd>=last_junction:
                return False
            else:
                return True


def junction_simple_merge(bigg_list):
    """
    merge the bigg with same junction to ref list
    merge the bigg with the 5' missing junction to ref list

    todo?:
    leave a switch to capture the intron retention or not
    if no, merge the intron retention isoform to the nearby one
    if yes, add these back to the ref list

    and the sub-read is divided into two parts: unique and un-unique (function added to post)

    after the clustering, add the most freq used 5' and 3' to the isoform as a representitive

    :param bigg_list:
    :return:
    """
    #from copy import deepcopy
    #bigg_list_bk=deepcopy(bigg_list)

    fullset=set(range(len(bigg_list)))
    drop=set()

    for i, bigg1 in enumerate(bigg_list):
        bigg1.get_junction()

        for bigg2 in bigg_list:
            if bigg1.name==bigg2.name:
                pass
            else:
                # single exon
                if len(bigg1.junction) == 0:
                    if is_single_exon_in(bigg1, bigg2):
                        drop.add(i)
                        bigg2.subread.add(bigg1.name)
                        bigg2.subread=bigg2.subread.union(bigg1.subread)
                # normal reads
                else:
                    if is_junction_inside(bigg1, bigg2):
                        drop.add(i)
                        bigg2.subread.add(bigg1.name)
                        bigg2.subread=bigg2.subread.union(bigg1.subread)

    keep=fullset-drop
    # change the default score of gene, no need to add
    for n, bigg in enumerate(bigg_list):
        if bigg.ttype=="isoform_anno":
            keep.add(n)

    keepl = sorted(list(keep))

    # write to the main
    for i in bigg_list:
        i.write_subread()

    bigg_list_new = select_list(bigg_list, keepl)

    return bigg_list_new


def flow_junction_cluster(bigg_list, bigg_ref):

    bigg_n= junction_pre(bigg_list, bigg_ref)
    bigg_merge=bigg_n+bigg_ref
    bigg_subread=junction_simple_merge(bigg_merge)

    return bigg_subread
