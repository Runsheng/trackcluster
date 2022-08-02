#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 11/16/2018 2:50 PM
# @Author  : Runsheng     
# @File    : batch.py
"""
make batch runs for different genes
"""

##std lib import
import os

## self import
from trackcluster.tracklist import read_bigg, write_bigg, add_subread_bigg, bigg_count_write_native, junction_pre
from trackcluster.utils import count_file
from trackcluster.cluster import flow_cluster, prefilter_smallexon
from trackcluster.utils import log_detail_file
from trackcluster.clusterj import junction_simple_merge
from trackcluster.clustercj import flow_junction_correct


# bar for run


#random.seed(1234)


# run the clustering for each gene
def process_one_junction_corrected_try(key, full=False, batchsize=500):

    gff_file = "./" + key + "/" + key + "_gff.bed"
    nano_file = "./" + key + "/" + key + "_nano.bed"
    bigg_out = "./" + key + "/" + key + "_simple_coveragej.bed"
    bigg_unused = "./" + key + "/" + key + "_unused.bed"
    log_file="./" + key + "/" + key + "_jrun.log"
    #print(key)

    logger = log_detail_file(log_file)
    logger.info(key+",start")

    if full is False:
        if os.stat(nano_file).st_size == 0:  # no bigg nano file
            return 0
        if os.path.isfile(bigg_out):  # already processed
            return 2

    bigg_gff = read_bigg(gff_file)
    bigg_nano_raw = read_bigg(nano_file)

    error_ll = []

    if bigg_nano_raw is None:
        return 0

    ### add the bactsize part
    n_count = 100
    n = 0
    bigg_nano = junction_pre(bigg_nano_raw, bigg_gff)

    if len(bigg_gff)==1 and bigg_gff[0].ttype=="region_mark": # exclude the land_mark bed file
        bigg_gff=[]
    bigg_merge=bigg_gff+bigg_nano # keep the ref on top of the list
    bigg_nano, bigg_rare = flow_junction_correct(bigg_merge)


    try:
        while n < n_count and len(bigg_nano) > batchsize:
            # print "n=", n
            bigg_1 = bigg_nano[:batchsize]
            bigg_2 = bigg_nano[batchsize:]
            bigg_new = junction_simple_merge(bigg_1)
            bigg_nano = add_subread_bigg(bigg_new+ bigg_2)
            n += 1

        # less than batchsize and the last collect
        bigg_new = junction_simple_merge(bigg_nano)

        logger.info([key, "isoform called:",len(bigg_new),
                               ";read with rare junction:",len(bigg_rare)])

        bigg_count_write_native(bigg_new, bigg_gff,bigg_out)
        write_bigg(bigg_rare, bigg_unused)
        logger.info(key + ",end")

    # need to add a sanity check if all bigg are there

    except Exception as e:
        error_ll.append(e)
        logger.info(("Error", e))

    return 1  # real run



def process_one_subsample_try(key, batchsize=2000, intronweight=0.5, by="ratio_all", full=False,
                              cutoff1=0.05, cutoff2=0.001,  scorecutoff=11):
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
        bigg_nano = junction_pre(bigg_nano_raw, bigg_gff)
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
                _, bigg_list_by1 = flow_cluster(bigg_1, bigg_gff, by=by, intronweight=intronweight,
                                                )
                bigg_nano = add_subread_bigg(bigg_list_by1 + bigg_2)
                n += 1

            # print len(bigg_nano)
            D, bigg_nano_new = flow_cluster(bigg_nano, bigg_gff, by, intronweight=intronweight)
            bigg_nano_new = add_subread_bigg(bigg_nano_new)

            ### save files
            for bigg in bigg_nano_new:
                bigg.write_subread()

            bigg_count_write_native(bigg_nano_new,bigg_gff, out=biggout)
            logger.info("Write bigg outfile to {}".format(biggout))

        except Exception as e:
            if e is ResourceWarning:
                pass
            else:
                logger.info("Error in cluster level:" + str(e))
    except Exception as e:
        logger.info("Error in prefilter level:" + str(e))
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