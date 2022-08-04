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


def __filter_site_dic(site_dic, ref=None):
    """
    filter away the sites which is in other chro or strand
    use the junction with highest freq as ref?
    :param site_dic:
    :param ref: in the format of tuple("chro", "_", "strand")
    :return: new_site_dic
    """
    site_dic_new={}
    if ref is None:
        cov_max=max(site_dic.values())
        for k, v in list(site_dic.items()):
            if v==cov_max:
                ref=k
        chro_ref, _, strand_ref=ref
    else:
        chro_ref, _, strand_ref=ref

    for k, v in list(site_dic.items()):
        chro, _, strand=k
        if chro==chro_ref and strand==strand_ref:
            site_dic_new[k]=v
    return site_dic_new


def chose_read_from_list(name, bigg_list):
    for bigg in bigg_list:
        if bigg.name==name:
            return bigg
    return None


def __get_corrected_junction(bigg, junction_dic, coverage_cutoff=2, offset=5):
    """
    too slooooow, as need to iter the full junction set, have revised in clustercj as get_bigg_correct
    Use the high confident junctions to correct the low frequent junctions
    junction dic like : (chro, pos, strand): coverage
    """
    strand=bigg.strand
    chro=bigg.chrom
    # sanity check

    bigg.get_junction()

    # get the pos that can be used
    pos_pass=[]
    for k, v in list(junction_dic.items()):
        if v>=coverage_cutoff:
            chro_d, pos_d, strand_d=k
            if chro_d==chro and strand_d==strand:
                pos_pass.append(pos_d)
    pos_pass_set=set(pos_pass)

    ###
    bigg_junction_new=[]
    corrected=[] # used for debug

    for pos in bigg.junction:
        if pos in pos_pass_set:
            bigg_junction_new.append(pos)
            corrected.append((pos, pos))
        else:
            for pos_d in pos_pass:
                if abs(pos_d-pos)<=offset:
                    bigg_junction_new.append(pos_d)
                    corrected.append((pos, pos_d))
                    break
                else: # can not be corrected, keep as it is
                    bigg_junction_new.append(pos)
    # debug
    # print("corrected", corrected)
    return bigg_junction_new


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


def __chain_site(sites_l, offset=5):
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
