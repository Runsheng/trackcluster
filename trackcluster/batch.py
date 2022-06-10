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
from trackcluster.tracklist import read_bigg, write_bigg, add_subread_bigg, bigg_count_write_native, merge_subread_bigg
from trackcluster.utils import count_file
from trackcluster.cluster import flow_cluster, prefilter_smallexon, write_D
from trackcluster.utils import log_summary, log_detail_file
from trackcluster.clusterj import flow_junction_cluster

#random.seed(1234)

# run the clustering for each gene


def process_one_junction_try(key, full=False, batchsize=500):
    print(key)
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    # figout = "./" + key + "/" + key + "_coverage.pdf"
    biggout = "./" + key + "/" + key + "_simple_coveragej1000.bed"
    # Dout = "./" + key + "/" + key + "_simple_coverage.csv"

    if full is False:
        if os.stat(nano_file).st_size == 0:  # no bigg nano file
            return 0
        if os.path.isfile(biggout):  # already processed
            return 2

    bigg_gff = read_bigg(gff_file)
    bigg_nano = read_bigg(nano_file)

    error_ll = []

    if bigg_nano is None:
        return 0

    ### add the bactsize part
    n_count = 100
    n = 0

    try:
        while n < n_count and len(bigg_nano) > batchsize:
            bigg_1 = bigg_nano[:batchsize]
            bigg_2 = bigg_nano[batchsize:]
            bigg_subread_by1 = flow_junction_cluster(bigg_1, bigg_gff)
            bigg_2.extend(bigg_subread_by1)
            n += 1

        bigg_subread = flow_junction_cluster(bigg_2, bigg_gff)
        write_bigg(bigg_subread, biggout)
    except Exception as e:
        error_ll.append(e)
        print(("Error", e))

    if len(bigg_subread)>=batchsize:
        return 3  # unfinished run
    return 1  # real run


def process_one_subsample_try(key, batchsize=500, intronweight=0.5, by="ratio_all", full=False):
    """
    batch run script for original trackcluster using exon/intron itersections
    :param key:
    :param batchsize:
    :param intronweight:
    :param by:
    :param full:
    :return:
    """
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    biggout = "./" + key + "/" + key + "_simple_coverage.bed"
    log_file="./" + key + "/" + key + "_run.log"

    if full is False:
        if os.stat(nano_file).st_size == 0:  # no bigg nano file
            return 0
        if os.path.isfile(biggout):  # already processed
            return 0

    bigg_gff = read_bigg(gff_file)
    bigg_nano_raw = read_bigg(nano_file)

    logger = log_detail_file(log_file)
    logger.info(key)

    try:
        bigg_nano = prefilter_smallexon(bigg_nano_raw, bigg_gff, cutoff=50)
        n_count = 100 # hard code n count for 50 times, if the reads > batchsize*time, will be
        # processed within 1000*100=100K, and the final isoform level can not be
        # more than 1000(batchsize)
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

            ### save files
            for bigg in bigg_nano_new:
                bigg.write_subread()

            bigg_count_write_native(bigg_nano_new, out=biggout)
            logger.info("Write bigg outfile to {}".format(biggout))

        except Exception as e:
            logger.debug("Error in cluster level:" + e)
    except Exception as e:
        logger.debug("Error in prefilter level:" + e)
    # merge_subread_bigg(bigg_nano_new)
    return 1


def get_len(key):
    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"

    line_gff = count_file(gff_file)
    line_nano = count_file(nano_file)

    return (line_gff, line_nano)


if __name__=="__main__":
    pass
    #os.chdir("./genes/")
    #key = "unc52"
    #print((process_one_subsample_try(key, batchsize=500, full=True)))