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
from trackcluster.track import bigGenePred

import pandas


def wrapper_bedtools_intersect2_select(bedfile1,bedfile2,outfile=None,fraction_bed1=0.05, fraction_bed2=0.1):
    """
    Using two bedfile to get the intsersection of pairs
    use simular fraction as substract instead of 0.2, the bigg can contain large intron in first intersect

    # -f 0.01= 1% from bigg, -F 0.05, means 5% from gene annotation

    :param bigg_one: the read
    :param bigg_two: the ref
    :param: fraction, may need to inlcude some read track with very long
    :return:
    """
    if outfile is None:
        prefix1=get_file_prefix(bedfile1)
        prefix2=get_file_prefix(bedfile2)
        location=get_file_location(bedfile1)

        outfile=location+"/"+"_".join([prefix1, prefix2])+".bed"

    # generate the bedfile, -r is reciprocal overlap, -s is strandedness/keep the same strand

    # sometime the name check for an intact bed and an non-intact one can result in a blank file
    # add -nonamecheck
    cmd="bedtools intersect -nonamecheck -wa -wb  -s -f {f1} -F {f2} -a {bedfile1} -b {bedfile2}>{out}".format(
        bedfile1=bedfile1, bedfile2=bedfile2, out=outfile, f1=fraction_bed1, f2=fraction_bed2)

    _=myexe(cmd)

    ### cleanup
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
    cmd="bedtools merge -nonamecheck -s  -c 4 -o count -i {biggfile} > {out}".format(
        biggfile=bigg_file, out=out
    )
    _=myexe(cmd)

    return out


def wrapper_bedtools_subtract(bigg_regionmark, bigg_gff, bigg_regionmark_f, f1=0.01, f2=0.01):
    """

    compare the region mark with original gff to remove
    :param bigg_regionmark:
    :param bigg_gff:
    :return:
    """
    # strandness, -A remove entire

    # the substract will change the original novel file
    # the output could be 0 if re-run the order
    sorted1=bigg_regionmark+"_s"
    sorted2=bigg_gff+"_s"


    def get_cmd_sort(infile, outfile):
        cmd_sort="bedtools sort -i {infile}> {outfile}".format(infile=infile, outfile=outfile)
        return cmd_sort

    cmd1=get_cmd_sort(bigg_regionmark, sorted1)
    cmd2=get_cmd_sort(bigg_gff,sorted2)

    myexe(cmd1)
    myexe(cmd2)

    # cutoff, if the novel region contain the following, then exclude
    # -f 0.01= 1% from region mark, -F 0.05, means 5% from gene annotation
    cmd="bedtools subtract -nonamecheck -s -A -f {f1} -F {f2} -a {bigg_regionmark} -b {bigg_gff}> {bigg_regionmark_f}".format(
        bigg_regionmark=sorted1, bigg_gff=sorted2, bigg_regionmark_f=bigg_regionmark_f, f1=f1, f2=f2
    )
    #print(cmd)
    _=myexe(cmd)

    del_files([sorted1, sorted2])
    return bigg_regionmark_f


def mergedbed2bigg(merge_bedfile, count_cutoff=5):
    """
    read the output from bedtools merge
    create a pseudo reference bigg for further use
    :param merge_bedfile:
    :return: a region_mark bigg used to get novel genes to new region
    """
    bigglist=[]
    with open(merge_bedfile, "r") as f:
        for line in f.readlines():
            line_l=line.strip().split("\t")
            chro, start, end, strand, count=line_l

            if int(count)>=count_cutoff:
                # generate a chro_start_end string for name
                str_l=["NOVEL",chro,start,end, "f"] if strand=="+" else ["NOVEL",chro, start, end, "r"]
                name_str="_".join(str_l)

                bigg=bigGenePred()
                bigg.chrom=chro
                bigg.chromStart=start
                bigg.chromEnd=end
                bigg.strand=strand

                # create a single exon gene
                bigg.chromStarts=[0]
                bigg.blockSizes=[int(end)-int(start)]
                bigg.blockCount=1
                bigg.exonFrames=[-1]
                # end for single exon

                bigg.name=name_str
                bigg.geneName=name_str # gene will be added to geneName

                bigg.ttype="region_mark"

                bigglist.append(bigg)
    return bigglist



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




def flow_write_gene_bigg(bigg_nano, bigg_ref, f1=0.05, f2=0.1):

    out=wrapper_bedtools_intersect2_select(bigg_nano, bigg_ref, outfile=None, fraction_bed1=f1, fraction_bed2=f2)
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

