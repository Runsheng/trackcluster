#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/2018 1:32 PM
# @Author  : Runsheng     
# @File    : tracklist.py

"""
Functions to handel the IO of bigg list
"""

# self import
from track import bigGenePred


def add_sw(bigg_file, sw_file, out="bigg_sw.bed"):
    """
    To add the
    :param bigg_list:
    :return:
    """
    sw_dic={}
    with open(sw_file, "r") as f:
        for line in f.readlines():
            line_l=line.split(",")
            name=line_l[0]
            score=int(line_l[1])
            sw_dic[name]=score

    bigg_list= read_bigg(bigg_file)

    bigg_str=[]
    for bigg_one in bigg_list:
        score=sw_dic[bigg_one.name]
        bigg_one.score=score
        bigg_str.append(bigg_one.to_str())

    with open(out, "w") as fw:
        fw.write("\n".join(bigg_str))


def read_bigg(bigg_file):
    """

    :param bigg_file:
    :return: bigg_list
    """
    bigg_list=[]
    with open(bigg_file, "r") as f:
        for line in f.readlines():
            bigg_one=bigGenePred()
            bigg_one.from_string(line.strip())
            bigg_list.append(bigg_one)

    return bigg_list

def write_bigg(bigg_list, out="bigg_new.bed"):

    bigg_str=[]
    for bigg_one in bigg_list:
        bigg_str.append(bigg_one.to_str())

    with open(out, "w") as fw:
        fw.write("\n".join(bigg_str))
