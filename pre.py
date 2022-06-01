#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 12/4/2018 3:43 PM
# @Author  : Runsheng     
# @File    : pre.py
"""
Functions to build the folders to start a batch run
contain tests in __main__
"""
from .utils import del_files, myexe
from .tracklist import get_file_prefix, get_file_location, count_file
import pandas


def wrapper_bedtools_intersect2_select(bedfile1,bedfile2,outfile=None):
    """
    Using two bedfile to get the intsersection of pairs
    :param bigg_one:
    :param bigg_two:
    :return:
    """
    if outfile is None:
        prefix1=get_file_prefix(bedfile1)
        prefix2=get_file_prefix(bedfile2)
        location=get_file_location(bedfile1)

        outfile=location+"/"+"_".join([prefix1, prefix2])+".bed"

    # generate the bedfile

    cmd="bedtools intersect -wa -wb -a {bedfile1}_s -b {bedfile2}_s>{out}".format(
        bedfile1=bedfile1, bedfile2=bedfile2, out=outfile)

    _=myexe(cmd)

    ### cleanup
    bed1s=bedfile1+"_s"
    bed2s=bedfile2+"_s"
    del_files([bedfile1, bedfile2, bed1s, bed2s])

    return outfile


def pandas_select(bed32file):
    """
    The bef8file is chr start end name *2 format
    :param bed8file:
    :return: the dict with (read1, read2): intersection
    """
    if count_file(bed32file) == 0:
        return {}

    df = pandas.read_csv(bed32file, sep="\t", header=None)

    df.drop_duplicates()

    readname = list(set(df[23]))

    return readname

