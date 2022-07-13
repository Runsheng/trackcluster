#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 12/4/2018 3:43 PM
# @Author  : Runsheng     
# @File    : pre.py
"""
Functions to build the folders to start a batch run
contain tests in __main__
"""
from trackcluster.utils import del_files, myexe, get_file_prefix, get_file_location
from trackcluster.tracklist import count_file

import pandas


def wrapper_bedtools_intersect2_select(bedfile1,bedfile2,outfile=None,fraction=0.2):
    """
    Using two bedfile to get the intsersection of pairs
    :param bigg_one: the read
    :param bigg_two: the ref
    :return:
    """
    if outfile is None:
        prefix1=get_file_prefix(bedfile1)
        prefix2=get_file_prefix(bedfile2)
        location=get_file_location(bedfile1)

        outfile=location+"/"+"_".join([prefix1, prefix2])+".bed"

    # generate the bedfile, -r is reciprocal overlap, -s is strandedness/keep the same strand

    cmd="bedtools intersect -wa -wb -r  -s -f {fraction} -a {bedfile1} -b {bedfile2}>{out}".format(
        bedfile1=bedfile1, bedfile2=bedfile2, out=outfile, fraction=fraction)

    _=myexe(cmd)

    ### cleanup
    bed1s=bedfile1+"_s"
    bed2s=bedfile2+"_s"
    #del_files([bedfile1, bedfile2, bed1s, bed2s])

    return outfile


def wrapper_bedtools_merge(bigg_file, out):
    """
    run a bedtools merge to get the un-annotated region, which contains reads
    the region will use the merged region as reference
    give a name and a GeneName as novelgene_chr_start_end

    This novogene bed file will be used to create novel gene folders

    :param bigglist:
    :return: a bed filename, with col5: [chro, start, end, strand, count]
    """

    # merge with strandness, report to col5
    cmd="bedtools merge -s  -c 4 -o count -i {biggfile} > {out}".format(
        biggfile=bigg_file, out=out
    )
    _=myexe(cmd)

    return out


def get_gendic_bedinter(interfile):
    """
    read the bedtools intersection file for the annotation bed and read bed
    # consdier the read is bed1 and gff is bed2
    :param interfile:
    :return:[readname:(gene1, gene2...)]
    """
    read_gene = {}
    f = open(interfile,"r")
    for line in f.readlines():
        line_l = line.split("\t")
        name = line_l[3]
        gene = line_l[-3]
        try:
            read_gene[name].add(gene)
        except KeyError:
            read_gene[name]=set()
            read_gene[name].add(gene)

    f.close()
    return read_gene


def group_bigg_by_gene(bigglist):
    """
    reverse the process of add gene
    :param bigglist:
    :return: {genename:bigglist} [bigglist_novelgene]
    """

    gene_bigg = {}
    bigg_novelgene=[]

    for bigg in bigglist:

        # also collect the track which did not have gene annotation
        if bigg.geneName=="none":
            bigg_novelgene.append(bigg)
        else:
            for gene in bigg.geneName.split("||"):
                try:
                    gene_bigg[gene].append(bigg)
                except KeyError:
                    gene_bigg[gene] = []
                    gene_bigg[gene].append(bigg)
    return gene_bigg, bigg_novelgene


def tracklist_add_gene(bigg_nano, read_gene):
    """
    :param bigg_nano: read track to bed added
    :param read_gene:
    :return: bigg_nano with genName has mutiple gene1||gene2||gene3
    """
    bigg_new=[]
    for bigg in bigg_nano:
        try:
            gene_set=read_gene[bigg.name]
            genename_str="||".join( list(gene_set) )
            bigg.geneName=genename_str
        except KeyError:
            pass
        bigg_new.append(bigg)
    return bigg_new






def flow_write_gene_bigg(bigg_nano, bigg_ref, fraction=0.2):

    out=wrapper_bedtools_intersect2_select(bigg_nano, bigg_ref, outfile=None, fraction=fraction)
    read_gene=get_gendic_bedinter(out)

    return read_gene



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

