#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/29/2018 2:34 PM
# @Author  : Runsheng     
# @File    : batch_test.py


"""
insert a full function key value test for the plotlib
"""
import os
import unittest
from tracklist import read_bigg, write_bigg
from track import bigGenePred
from plots import line_plot_merge
from collections import OrderedDict

def process_one(key):
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    figout = "./" + key + "/" + key + ".pdf"
    biggout = "./" + key + "/" + key + "_simple.bed"
    Dout = "./" + key + "/" + key + "_simple.csv"

    if os.stat(nano_file).st_size == 0:  # no bigg nano file
        return 0
    if os.path.isfile(biggout):  # already processed
        return 0

    bigg = []
    with open(gff_file) as f:
        for line_one in f.readlines():
            bigg_one = bigGenePred()
            bigg_one.from_string(line_one)
            bigg.append(bigg_one)
    with open(nano_file) as f:
        for line_one in f.readlines():
            bigg_one = bigGenePred()
            bigg_one.from_string(line_one)
            bigg.append(bigg_one)

    line_plot_merge(bigg, out=figout,
                    biggout=biggout,
                    Dout=Dout,
                    intronweight=0.5,
                    by="ratio_all", core=40)
    return 1



def count_file(thefile):
    count = 0
    for line in open(thefile).xreadlines(  ):
        count += 1
    return count


def get_len(key):
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    figout = "./" + key + "/" + key + ".pdf"
    biggout = "./" + key + "/" + key + "_simple.bed"
    Dout = "./" + key + "/" + key + "_simple.csv"

    line_gff = count_file(gff_file)
    line_nano = count_file(nano_file)

    return (line_gff, line_nano)


if __name__ == '__main__':

    # get the full gene names
    os.chdir("/home/zhaolab1/data/nanorna/trackall2")
    gene_bigg = read_bigg("./ce10_gff.bed")
    gene_s = set()
    gene = []
    for bigg in gene_bigg:
        if bigg.geneName not in gene_s:
            gene.append(bigg.geneName)
            gene_s.add(bigg.geneName)
    print len(gene), len(gene_s)

    gene_lendic = OrderedDict()
    gene_lenlist = []
    for key in gene:
        line_gff, line_nano = get_len(key)
        gene_lendic[key] = (line_gff, line_nano)
        gene_lenlist.append((key, line_gff, line_nano))

    key_300 = [x for x in gene_lenlist if x[2] <= 100]
    print(len(key_300))
    key_300_gene = [x[0] for x in key_300]

    for key in key_300_gene:
        print key
        if key=="Y23H5A.8":
            pass
        else:
            try:
                process_one(key)
            except Exception:
                pass
