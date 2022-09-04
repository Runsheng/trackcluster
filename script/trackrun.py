#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 3/10/2020 2:45 PM
# @Author  : Runsheng     
# @File    : trackrun.py

"""
A wrapper for the general processing of the trackcluster jobs
"""

import argparse
import sys
import os
import inspect
import logging
from collections import OrderedDict

# non-std lib
from tqdm import tqdm

from trackcluster.tracklist import read_bigg, write_bigg
from trackcluster.utils import is_bin_in, is_package_installed, get_file_prefix
from trackcluster import __version__
from trackcluster.flow import flow_clusterj_all_gene_novel, flow_cluster_all_gene_novel, flow_count

logger = logging.getLogger('summary')
logger.setLevel(logging.INFO)

class CMD(object):

    def __init__(self):
        parser=argparse.ArgumentParser(
            description="Trackcluster cmd lines",
            usage=""" trackrun.py <command> [<args>]

-------
The command contains:
pre: prepare the run folder by separating the tracks from the same locus
clusterj:run the trackcluster in junction mode to get the isoforms and the counting
cluster:run the trackcluster eon/intron intersection to get the isoforms and the counting
desc: compare the novel isoforms with the existing annotations to give an description for the min edit distance between a novel isoform with its nearest reference annotation
test: run test for installation 
------
Examples for running command:
trackrun.py clusterj -s reads.bed -r refs.bed -t 40 # run in junction mode
trackrun.py cluster -s reads.bed -r refs.bed -t 40 # run in exon/intron intersection mode， slower
trackrun.py desc --isoform isoform.bed --reference ref.bed > desc.bed 

# test if all dependencies are installed
trackrun.py test --install

version {version}
            """.format(version =__version__),
            formatter_class = argparse.RawDescriptionHelpFormatter)
        parser.add_argument("command", help="Subcommand to run")
        args=parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)
        getattr(self, args.command)()



    def cluster(self):
        parser=argparse.ArgumentParser(
            description="Original cluster, using intersection of intron and exon as matrix"
        )
        parser.add_argument("-d", "--folder", default=os.getcwd(),
                            help="the folder contains all the seperated tracks in different locus/genes, default is the current dir")

        parser.add_argument("-s", "--sample", help="the bigg format of the read track, with the key of GeneName")
        parser.add_argument("-r", "--reference", help="the bigg format of the reference annotation track")
        parser.add_argument("-t", "--thread", default=32, type=int,
                            help="the max thread used to run some of the process")
        parser.add_argument("-f1", "--intersect1", default=0.01, type=float,
                            help="the min intersection in read track, used to assign gene name to a read using bedtools")
        parser.add_argument("-f2", "--intersect2", default=0.05, type=float,
                            help="the min intersection in isoform track, used to assign gene name to a read using bedtools")
        parser.add_argument("-c", "--count", default=5, type=int,
                            help="the min cutoff for a novel isoform be retained in counting")
        parser.add_argument("-p", "--prefix", default=None, type=int,
                            help="the min cutoff for a novel isoform be retained in counting")
        parser.add_argument("-b", "--batchsize", default=2000, type=int,
                            help="the max reads can be processed in one batch")

        ###cluster specific para
        parser.add_argument("--tmp", default=None,
                            help="the tmp dir, default is the current dir/tmp")
        parser.add_argument("-c1", "--cutoff1", default=0.05, type=float,
                            help="the dissimilarity cutoff to merge reads in round1 (equal) merge")
        parser.add_argument("-c2", "--cutoff2", default=0.05, type=float,
                            help="the dissimilarity cutoff to merge reads in round2 (within) merge")
        parser.add_argument("-w", "--intronweight", default=0.5, type=float,
                            help="the weight for the intron in read cluster")

        arg_use = sys.argv[2:]
        if len(arg_use) >= 4:
            args = parser.parse_args(arg_use)
        else:
            parser.print_help()
            sys.exit(1)
        if args.prefix is None:
            args.prefix=get_file_prefix(args.sample,sep=".")

        #args = parser.parse_args(args=None if sys.argv[2:] else ['--help'])

        flow_cluster_all_gene_novel(wkdir=args.folder,
                                 prefix=args.prefix,
                                 nano_bed=args.sample,
                                 gff_bed=args.reference,
                                 core=args.thread,
                                 f1=args.intersect1,
                                 f2=args.intersect2,
                                 count_cutoff=args.count,
                                 batchsize=args.batchsize,
                                 intronweight=args.intronweight,
                                 cutoff1=args.cutoff1,
                                 cutoff2=args.cutoff2)


    def clusterj(self):
        parser = argparse.ArgumentParser(
            description="cluster using junction information"
        )
        parser.add_argument("-d", "--folder", default=os.getcwd(),
                            help="the folder contains all files, default is the current dir")
        parser.add_argument("-s", "--sample", help="the bigg format of the read track, with group information in GeneName2")
        parser.add_argument("-r", "--reference", help="the bigg format of the reference annotation track")
        parser.add_argument("-t", "--thread", default=32, type=int,
                            help="the max thread used to run some of the process")
        parser.add_argument("-f1", "--intersect1", default=0.01, type=float,
                            help="the min intersection in read track, used to assign gene name to a read using bedtools")
        parser.add_argument("-f2", "--intersect2", default=0.05, type=float,
                            help="the min intersection in isoform track, used to assign gene name to a read using bedtools")
        parser.add_argument("-c", "--count", default=5, type=int,
                            help="the min cutoff for a novel isoform be retained in counting")
        parser.add_argument("-p", "--prefix", default=None, type=int,
                            help="the min cutoff for a novel isoform be retained in counting")
        parser.add_argument("-b", "--batchsize", default=2000, type=int,
                            help="the max reads can be processed in one batch")
        #args = parser.parse_args(args=None if sys.argv[2:] else ['--help'])
        arg_use=sys.argv[2:]
        if len(arg_use)>=4:
            args = parser.parse_args(arg_use)
        else:
            parser.print_help()
            sys.exit(1)

        if args.prefix is None:
            args.prefix=get_file_prefix(args.sample,sep=".")

        flow_clusterj_all_gene_novel(wkdir=args.folder,
                                     prefix=args.prefix,
                                     nano_bed=args.sample,
                                     gff_bed=args.reference,
                                     core=args.thread,
                                     f1=args.intersect1,
                                     f2=args.intersect2,
                                     count_cutoff=args.count,
                                     batchsize=args.batchsize)

    def count(self):
        parser = argparse.ArgumentParser(
            description="Counting the cluster result isoform file to get the expression csv"
        )
        parser.add_argument("-d", "--folder", default=os.getcwd(),
                            help="the folder contains all files, default is the current dir")
        parser.add_argument("-s", "--sample",
                            help="the bigg format of the read track, with group information in GeneName2")
        parser.add_argument("-r", "--reference", help="the bigg format of the reference annotation track")
        parser.add_argument("-i", "--isoform",
                            help="the isoform bed file from clustering")
        parser.add_argument("-p", "--prefix", default=None, type=int,
                            help="the min cutoff for a novel isoform be retained in counting")

        arg_use=sys.argv[2:]
        if len(arg_use)>=4:
            args = parser.parse_args(arg_use)
        else:
            parser.print_help()
            sys.exit(1)

        if args.prefix is None:
            args.prefix=get_file_prefix(args.sample,sep=".")

        flow_count(wkdir=args.folder,
                   prefix=args.prefix,
                   nano_bed=args.sample,
                   isoform_bed=args.isoform,
                   gff_bed=args.reference)

    def desc(self):
        parser = argparse.ArgumentParser(
            description="Counting the cluster result isoform file to get the expression csv"
        )
        parser.add_argument("-d", "--folder", default=os.getcwd(),
                            help="the folder contains all files, default is the current dir")
        parser.add_argument("-s", "--sample",
                            help="the bigg format of the read track, with group information in GeneName2")
        parser.add_argument("-r", "--reference", help="the bigg format of the reference annotation track")
        parser.add_argument("-i", "--isoform",
                            help="the isoform bed file from clustering")
        parser.add_argument("-p", "--prefix", default=None, type=int,
                            help="the min cutoff for a novel isoform be retained in counting")

        arg_use=sys.argv[2:]
        if len(arg_use)>=4:
            args = parser.parse_args(arg_use)
        else:
            parser.print_help()
            sys.exit(1)

        if args.prefix is None:
            args.prefix=get_file_prefix(args.sample,sep=".")

    def test(self):
        """
        run simple tests using the data files contained in the test
        :return:
        """
        currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        parentdir = os.path.dirname(os.path.dirname(currentdir))
        #sys.path.insert(0,parentdir)
        ### need to use the file from the test path to run the test
        parser = argparse.ArgumentParser(
            description="Test functions"
        )
        parser.add_argument("--install", action='store_true',
                            help="test the install of all needed packages")
        parser.add_argument("--pre", action="store_true",
                            help="test the pre and makedir functions")
        args = parser.parse_args(sys.argv[2:])

        if args.install: ## test installed packages, samtools, bedtools
            if is_bin_in("samtools") and is_bin_in("bedtools") and is_bin_in("minimap2"):
                logger.info("Pass")
            else:
                logger.info("Check if samtools, bedtools and minimap2 are in $PATH")

            for package_name in ["pysam", "Bio", "numpy", "pandas", "tqdm"]:
                if is_package_installed(package_name):
                    logger.info("Package {} installed".format(package_name) )
                else:
                    logger.info("Import Error: Package {} not installed".format(package_name) )


        if args.pre: ### test the prepare function using test
           pass

if __name__ =="__main__":
    CMD()
