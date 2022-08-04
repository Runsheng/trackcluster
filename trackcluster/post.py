#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2/15/2019 3:05 PM
# @Author  : Runsheng     
# @File    : post.py

from collections import OrderedDict

from utils import group_site


def is_junction_equal(bigg0, bigg1, offset=10):
    """

    :param bigg0:
    :param bigg1: usually the ref bigg
    :return:
    """
    bigg0.get_exon()
    bigg0.get_junction()
    bigg1.get_exon()
    bigg1.get_junction()

    #print bigg0.name, bigg0.junction
    #print bigg1.name, bigg1.junction

    match_0=[]
    match_1=[]

    for ni, i in enumerate(bigg0.junction):
        for nj, j in enumerate(bigg1.junction):
            if -offset <= i - j <= offset:
                #line=(bigg0.name, ni, i, bigg1.name, nj, j)
                match_0.append(i)
                match_1.append(j)

    if set(bigg0.junction)==set(match_0) and set(bigg1.junction)==set(match_1):
        return True
    return False


def fuzzy_intersection(junction1, junction2, offset=10):
    """
    return a dic for fuzzy intersection, keep the coord using {junction1:junction2}
    """
    match_dic=OrderedDict()

    for ni, i in enumerate(junction1):
        for nj, j in enumerate(junction2):
            if -offset <= i - j <= offset:
                match_dic[i]=j

    return match_dic


def has_new_junction(bigg0, list_ref, offset=10):
    """
    :param bigg_one:
    :param list_ref: the reference contains all ref in, list_ref can be as len=1
    :return:
    """
    junctions=[]
    bigg0.get_junction()

    for bigg_ref in list_ref:
        bigg_ref.get_exon()
        bigg_ref.get_junction()

        junctions.extend(bigg_ref.junction)

    junction_s=sorted(list(set(junctions)))

    match_dic=fuzzy_intersection(bigg0.junction,
                               junction_s, offset)

    not_in=set(bigg0.junction)-set(match_dic.keys())

    ### debug
    #print bigg0.junction
    #print sorted(not_in)
    ### debug

    if len(not_in)>0:
        #print not_in
        return True
    else:
        return False


def class_4(bigg0, list_ref, offset=10):

    #class_name=None

    is_new_junction = has_new_junction(bigg0, list_ref, offset)

    if is_new_junction is True:
        return "new_junction"

    for i in list_ref:
        is_equal=is_junction_equal(bigg0, i, offset)

        if is_equal is True:
            # compare the start and end to define the >= and <
            if abs(bigg0.chromStart-bigg0.chromEnd)>=abs(i.chromStart-i.chromEnd):
                return "all_matched>=_{}".format(i.name)
            else:
                return "all_matched_<_{}".format(i.name)

    return "new_combination"


def compare_ei_by_boudary(bigg0, bigg_ref, offset=10):
    """
    compare the exon and intron one by one by compare the junctions

    :param bigg0:
    :param bigg_ref: ref
    :return:
    """
    bigg0.get_junction()
    bigg_ref.get_junction()

    match_dic=fuzzy_intersection(bigg0.junction,
                                 bigg_ref.junction,
                                 offset)

    # correct the bigg0 junction with ref junction
    junction_new=[]
    for i in bigg0.junction:
        if i in list(match_dic.keys()):
            junction_new.append(match_dic[i])
        else:
            junction_new.append(i)

    # get the dic of the order of the boundary
    posdic_bigg0={k: v for v, k in enumerate(junction_new)}
    posdic_ref={k: v for v, k in enumerate(bigg_ref.junction)}


    # junction as end, start, end...start

    # count the missing first, in ref, not in the bigg0
    missed=set(bigg_ref.junction)-set(junction_new)
    missed_order=sorted([posdic_ref[x] for x in missed])

    #print group_site(missed_order)

    # count the new then, these can only be explained by other refs
    extra=set(junction_new)-set(bigg_ref.junction)
    extra_order=sorted([posdic_bigg0[x] for x in extra])

    #print group_site(extra_order)

    return (missed_order, extra_order)


def find_nearest_ref(bigg0, list_ref, offset=10):
    # find the nearest ref by junction description
    # nearest is: with least extra; if extra equal, least missed
    # first the regions, then the numbers of the each site

    missed_s, extra_s=compare_ei_by_boudary(bigg0, list_ref[0], offset)
    group_miss_s= group_site(missed_s)
    group_extra_s= group_site(extra_s)
    s_n=0

    for n, bigg_ref in enumerate(list_ref):
        missed, extra=compare_ei_by_boudary(bigg0, bigg_ref, offset)
        group_miss = group_site(missed)
        group_extra = group_site(extra)

        if len(group_extra) > len(group_extra_s):
            pass
        elif len(group_extra) < len(group_extra_s):
            missed_s, extra_s, group_miss_s, group_extra_s=missed, extra, group_miss, group_extra
            s_n=n
        else: # equal then compare the extra number
            if len(extra) > len(extra_s):
                pass
            elif len(extra) < len(extra_s):
                missed_s, extra_s, group_miss_s, group_extra_s = missed, extra, group_miss, group_extra
                s_n = n
            else: # equal then compare the missed
                if len(group_miss) > len(group_miss_s):
                    pass
                elif len(group_miss) < len(group_miss_s):
                    missed_s, extra_s, group_miss_s, group_extra_s = missed, extra, group_miss, group_extra
                    s_n = n
                else: # equal miss group then miss real number
                    if len(missed) > len(missed_s):
                        pass
                    elif len(missed) < len(missed_s):
                        missed_s, extra_s, group_miss_s, group_extra_s = missed, extra, group_miss, group_extra
                        s_n = n
                    else: # all equal, then use previous one
                        pass

    return list_ref[s_n]


def desc_ei_by_boundary(bigg0, bigg_ref, offset=10):
    """
    give description to the pairwise compare
    :param bigg0:
    :param bigg_ref:
    :param offset:
    :return:
    """
    missed, extra = compare_ei_by_boudary(bigg0, bigg_ref, offset)
    group_miss = group_site(missed)
    group_extra = group_site(extra)

    ref_lastexon_no= len(bigg_ref.junction)//2+1

    ### debug
    #if len(missed) % 2 != 0 or len(extra) % 2 != 0:
    #    print bigg0.name, bigg_ref.name, group_miss, len(missed), group_extra, len(extra)
    ###

    miss_desc=[]
    extra_desc=[]

    ### desc the miss first
    if len(group_miss)==0:
        miss_desc.append("No miss exon.")

    for junctions in group_miss:
        # start and end, the middle can be two
        if junctions[0]==0: # start
            miss_desc.append("5 primer miss: exon 1 to {}".format(junctions[-1]//2+1))

        if junctions[-1]==len(bigg_ref.junction)-1: # end
            miss_desc.append("3 primer miss: exon {} to {}".format(
                junctions[0] // 2 + 2, ref_lastexon_no))

        if junctions[0]!=0 and junctions[-1]!=len(bigg_ref.junction)-1: # in the middle
            if junctions[0] % 2==0: # intron miss/retention
                miss_desc.append("Intron retention: intron {} to {}".format(
                    junctions[0]//2+1, junctions[-2]//2+1))

            if junctions[0] % 2==1: # exon miss
                miss_desc.append("Exon miss: exon {} to {}".format(
                    junctions[0]//2+2, junctions[-2]//2+2))

    ### the extra part

    if len(group_extra)==0:
        extra_desc.append("No extra exon.")

    for junctions in group_extra:
        # start and end, the middle can be two
        if junctions[0]==0: # start
            extra_desc.append("5 primer extra: exon 1 to {}".format(junctions[-1]//2+1))

        if junctions[-1]==len(bigg0.junction)-1: # end
            extra_desc.append("3 primer extra: exon {} to {}".format(
                junctions[0] // 2 + 2, junctions[-2] // 2 + 2))

        if junctions[0]!=0 and junctions[-1]!=len(bigg0.junction)-1: # in the middle
            if junctions[0] % 2==0: # intron miss/retention, rare in the extra part
                extra_desc.append("Intron extra: intron {} to {}".format(
                    junctions[0]//2+1, junctions[-2]//2+1))

            if junctions[0] % 2==1: # exon miss
                extra_desc.append("Exon extra: exon {} to {}".format(
                    junctions[0]//2+2, junctions[-2]//2+2))
    # debug
    #print len(bigg_ref.junction)
    #print group_miss, group_extra
    #print miss_desc, extra_desc

    return miss_desc, extra_desc


def desc_to_text(desc):
    """
    use ; to sep different lines of desc
    :param desc:
    :return:
    """

    return ";".join(desc)


def flow_desc(bigg0, list_ref, offset=10):
    """
    write a 5 col tab to describe the exon and intron usage of the bigg and its nearest ref

    :param bigg0:
    :param list_ref:
    :param offset:
    :return:
    """
    bigg_ref=find_nearest_ref(bigg0, list_ref, offset)
    miss_desc, extra_desc= desc_ei_by_boundary(bigg0, bigg_ref, offset)
    return [bigg0.name, bigg_ref.name, bigg_ref.geneName, desc_to_text(miss_desc), desc_to_text(extra_desc)]


def flow_class4(bigg0, list_ref, offset=10):
    class4 = class_4(bigg0, list_ref, offset)
    return [bigg0.name, class4]



