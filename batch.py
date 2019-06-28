#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 11/16/2018 2:50 PM
# @Author  : Runsheng     
# @File    : batch.py
"""
make batch runs for different genes
"""

##std lib import
import random
import os


## self import
from tracklist import read_bigg, write_bigg, add_subread_bigg, bigg_count_write, merge_subread_bigg
from plots import line_plot_merge
from utils import count_file
from cluster import flow_cluster, prefilter_smallexon, write_D

#random.seed(1234)

def process_one_subsample(key, batchsize=500, intronweight=0.5, by="ratio_all", full=False):
    # print key
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    biggout = "./" + key + "/" + key + "_simple_coverage.bed"

    if full is False:
        if os.stat(nano_file).st_size == 0:  # no bigg nano file
            return 0
        if os.path.isfile(biggout):  # already processed
            return 2

    bigg_gff = read_bigg(gff_file)
    bigg_nano_raw = read_bigg(nano_file)

    bigg_nano=prefilter_smallexon(bigg_nano_raw, bigg_gff, cutoff=50)
    n_count=50
    n=0
    #bigg_nano.sort(key=operator.attrgetter("chromStart"))

    while n < n_count and len(bigg_nano)>batchsize:
        print "n=", n
        bigg_1 = bigg_nano[:batchsize]
        bigg_2 = bigg_nano[batchsize:]
        _, bigg_list_by1 = flow_cluster(bigg_1, bigg_gff, by, intronweight=intronweight)
        bigg_nano = add_subread_bigg(bigg_list_by1 + bigg_2)
        n+=1

    #print len(bigg_nano)
    D,bigg_nano_new=flow_cluster(bigg_nano,bigg_gff, by, intronweight=intronweight)
    bigg_nano_new=add_subread_bigg(bigg_nano_new)

    #print len(bigg_nano_new)
    #write_D(D, bigg_nano_new,Dout)

    ### save nessary files
    for bigg in bigg_nano_new:
        bigg.write_subread()

    bigg_count_write(bigg_nano_new, out=biggout)
    #merge_subread_bigg(bigg_nano_new)

    return 1 # processed in this run


def process_one_subsample_try(key, batchsize=1000, intronweight=0.5, by="ratio_all", full=False):
    # print key
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    biggout = "./" + key + "/" + key + "_simple_coverage.bed"

    if full is False:
        if os.stat(nano_file).st_size == 0:  # no bigg nano file
            return 0
        if os.path.isfile(biggout):  # already processed
            return 0

    bigg_gff = read_bigg(gff_file)
    bigg_nano_raw = read_bigg(nano_file)

    try:
        bigg_nano = prefilter_smallexon(bigg_nano_raw, bigg_gff, cutoff=50)
        n_count = 100
        n = 0
        # bigg_nano.sort(key=operator.attrgetter("chromStart"))

        if bigg_nano is None:
            return 0

        try:
            while n < n_count and len(bigg_nano) > batchsize:
                # print "n=", n
                bigg_1 = bigg_nano[:batchsize]
                bigg_2 = bigg_nano[batchsize:]
                _, bigg_list_by1 = flow_cluster(bigg_1, bigg_gff, by, intronweight=intronweight)
                bigg_nano = add_subread_bigg(bigg_list_by1 + bigg_2)
                n += 1

            # print len(bigg_nano)
            D, bigg_nano_new = flow_cluster(bigg_nano, bigg_gff, by, intronweight=intronweight)
            bigg_nano_new = add_subread_bigg(bigg_nano_new)

            # print len(bigg_nano_new)
            # write_D(D, bigg_nano_new,Dout)

            ### save files
            for bigg in bigg_nano_new:
                bigg.write_subread()

            bigg_count_write(bigg_nano_new, out=biggout)
        except Exception as e:
            print e
    except Exception as e:
        print e
    # merge_subread_bigg(bigg_nano_new)
    return 1


def get_len(key):
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    figout = "./" + key + "/" + key + ".pdf"
    biggout = "./" + key + "/" + key + "_simple.bed"
    Dout = "./" + key + "/" + key + "_simple.csv"

    line_gff = count_file(gff_file)
    line_nano = count_file(nano_file)

    return (line_gff, line_nano)


if __name__=="__main__":
    os.chdir("./test")
    key = "unc52"
    process_one_subsample(key, batchsize=200, full=True)