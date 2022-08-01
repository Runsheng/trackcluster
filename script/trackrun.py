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
from trackcluster.flow import flow_clusterj_all_gene_novel

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
cluster:run the trackcluster main function to get the isoforms and the counting
desc: compare the novel isoforms with the existing annotations to give an description for the min edit distance between a novel isoform with its nearest reference annotation

test: run the test code for some gene models in the test folder
------
A example for running annotation command:
trackrun.py pre --reference ref.bed --tracks in.bed --out trackall
trackrun.py cluster --folder trackall
trackrun.py clusterj --folder trackall # run in junction mod
trackrun.py desc --isoform isoform.bed --reference ref.bed > desc.bed 

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




    def __pre(self):
        """
        prepare the data folder used for the
        :return:
        """
        parser=argparse.ArgumentParser(
            description="prepare the run dir for each genes"
        )
        parser.add_argument("-d", "--wkdir", default=os.getcwd(),
                            help="the working dir, default is the current dir")
        parser.add_argument("-t", "--tmp", default=None,
                            help="the working dir, default is the current dir")

        parser.add_argument("-s", "--sample", help="the bigg format of the read track, with the key of GeneName")
        parser.add_argument("-r", "--reference", help="the bigg format of the reference annotation track")
        args = parser.parse_args(sys.argv[2:])

        #### main flow
        logger.info("Start file reading")
        anno_bigg=read_bigg(args.reference)
        gene_anno=self.group_bigg_by_gene(anno_bigg)
        nano_bigg=read_bigg(args.sample)
        gene_nano = self.group_bigg_by_gene(nano_bigg)
        logger.info("End of file reading")

        ### change dir
        os_origin=os.getcwd()
        os.chdir(args.wkdir)

        logger.info("Start of dir making")
        for gene, nano_bigg in gene_nano.items():
            anno_bigg = gene_anno[gene]
            try:
                os.mkdir(gene)
            except OSError:
                pass

            anno_out = "./{gene}/{gene}_gff.bed".format(gene=gene)
            nano_out = "./{gene}/{gene}_nano.bed".format(gene=gene)

            write_bigg(anno_bigg, anno_out)
            write_bigg(nano_bigg, nano_out)

        os.chdir(os_origin)
        logging.info("END of dir making")

    def __cluster(self):
        parser=argparse.ArgumentParser(
            description="Original cluster, using intersection of intron and exon as matrix"
        )
        parser.add_argument("-d", "--folder", default=os.getcwd(),
                            help="the folder contains all the seperated tracks in different locus/genes, default is the current dir")
        parser.add_argument("-t", "--tmp", default=None,
                            help="the tmp dir, default is the current dir/tmp")
        parser.add_argument("-s", "--sample", help="the bigg format of the read track, with the key of GeneName")
        parser.add_argument("-r", "--reference", help="the bigg format of the reference annotation track")

        args = parser.parse_args(sys.argv[2:])
        #print(args.fasta)


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
        parser.add_argument("-b", "--batchsize", default=1000, type=int,
                            help="the max reads can be processed in one batch")
        #args = parser.parse_args(args=None if sys.argv[2:] else ['--help'])
        args = parser.parse_args(sys.argv[2:] if sys.argv[2:] is not None else ["--help"])

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

    def desc(self):
        pass


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
