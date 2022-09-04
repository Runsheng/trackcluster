#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/8/2018 2:47 PM
# @Author  : Runsheng     
# @File    : convert.py
"""
The processing functions used to change the format of tracks
"""

# self import
from trackcluster.track import bigGenePred
from trackcluster.gff import GFF
from trackcluster.utils import get_file_prefix

# third part import
from pysam import AlignmentFile


def sam_to_bigGenePred(record, samfile):
    """
    sam is 0 based and bigg is 0 based

    :param record:
    :param samfile: the opened Alignment file from pysam
    :return:
    """

    bigg = bigGenePred()

    # rename the values
    bigg.chrom = samfile.getrname(record.reference_id)

    bigg.name = record.query_name
    bigg.strand = "-" if record.is_reverse else "+"

    # give colour to reverse and forward strand, though unused in ucsc
    bigg.reserved = [64, 224, 208] if record.is_reverse else [250, 128, 114]

    bigg.name2 = record.query_name

    # the unchanged fields
    # bigg.core= 1000
    # bigg.cdsStartStat="none"
    # bigg.cdsEndStat="none"

    # bigg.ttype="nanopore_reads"
    # bigg.geneName=""
    # bigg.geneName2="" # used to store group info, can use the bam file name as indicator
    # bigg.geneType="none"


    # using the cigar, affact the following field
    cigaryield = cigar_count(record.cigartuples)

    bigg.chromStart = record.reference_start  # -cigaryield['lclipping']
    bigg.chromEnd = record.reference_end  # +cigaryield['rclipping']

    # bigg.thickStart=0
    # bigg.thickEnd=0
    bigg.blockSizes = cigaryield["len_site"]
    bigg.chromStarts = cigaryield["start_site"]

    try:
        assert len(bigg.blockSizes) == len(bigg.chromStarts)
    except AssertionError:
        print(len(bigg.blockSizes), len(bigg.chromStarts))
    bigg.blockCount = len(bigg.blockSizes)

    bigg.exonFrames = [-1 for i in range(0, bigg.blockCount)]

    return bigg


def cigar_count(cigar_tuple):
    """
    use the cigar tuple to get
    S in L and in R
    The position of N in the reletive chro
    """
    ## 4S 0M 1I 2D 3N, N is the RNAseq gap
    ## store the start and len just as in
    cigaryield = {}

    cigaryield["len_site"] = []
    cigaryield["start_site"] = [0]

    cigaryield["lclipping"] = 0
    cigaryield["rclipping"] = 0

    region_offset = 0
    region_start = 0

    for n, i in enumerate(cigar_tuple):
        tag, number = i
        if n == 0 and tag == 4:
            cigaryield["lclipping"] = number

        if tag == 0 or tag == 2:
            region_offset += number
        if tag == 1:
            pass
        if tag == 3:  # catch a splicing event, add the region into
            cigaryield["len_site"].append(region_offset)
            cigaryield["start_site"].append(region_start + region_offset + number)

            region_offset = 0
            region_start = cigaryield["start_site"][-1]

        if n == (len(cigar_tuple) - 1):  # catch the end event
            if tag == 4:
                cigaryield["rclipping"] = number

            cigaryield["len_site"].append(region_offset)

    return cigaryield


#### another function

def gff_to_bigGenePred(gff, indicator="ID"):
    """
    Note gff is 1 start and bigGenePred is 0
    :param gff: a GFF class with only one gene inside
    :return:
    """
    # may need to return multiple bigg
    bigg_list=[]

    if gff.transcript_d is None:
        gff.transcript_format(indicator=indicator)
    #print(list(gff.transcript_d.keys()))

    for key in list(gff.transcript_d.keys()):
        bigg = bigGenePred()
        try:
            gene=gff.transcript_to_gene[key]
        except KeyError:
            print("KeyERROR", key)

        for n, record in enumerate(gff.transcript_d[key]):
            #use the first line to get basic infor
            if n==0:
                # rename the values
                bigg.chrom = record.seqid
                bigg.name = key
                bigg.strand = record.strand

                # give colour to reverse and forward strand
                bigg.reserved = [64, 224, 208] if bigg.strand=="+" else [250, 128, 114]
                bigg.name2 = key

                # the unchanged fields
                bigg.score= 100 # give a high score
                # bigg.reserved=255,128,0
                # bigg.cdsStartStat="none"
                # bigg.cdsEndStat="none"

                bigg.ttype = "isoform_anno"
                bigg.geneName = gene
                # bigg.geneName2=""
                # bigg.geneType="none"
                # bigg.thickStart=0
                # bigg.thickEnd=0
                bigg.chromStart = record.start-1 # the ucsc 0 based is [start, end), while gff is [start, end]

            bigg.chromStarts.append(record.start-bigg.chromStart-1)
            bigg.blockSizes.append(record.end-record.start+1) # [1,2] is len2

        # use the last end as end
        bigg.chromEnd = record.end # the ucsc 0 based is [start, end), while gff is [start, end]

        try:
            assert len(bigg.blockSizes) == len(bigg.chromStarts)
        except AssertionError:
            print("Error in gff convert", len(bigg.blockSizes), len(bigg.chromStarts))
        bigg.blockCount = len(bigg.blockSizes)
        bigg.exonFrames = [-1 for i in range(0, bigg.blockCount)] # ignore frame

        bigg_list.append(bigg)

    return bigg_list


if __name__ == "__main__":
    pass