#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 3/10/2020 2:45 PM
# @Author  : Runsheng     
# @File    : trackcluster.py

"""
A wrapper for the general processing of the trackcluster jobs
"""

import argparse
import sys
import os
import inspect
from collections import OrderedDict
from .tracklist import read_bigg, write_bigg
import logging
from .utils import log_detail_file, log_summary, is_bin_in, is_package_installed

from .clusterj import flow_junction_cluster
logger = logging.getLogger('summary')
logger.setLevel(logging.INFO)

class CMD(object):

    def __init__(self):
        parser=argparse.ArgumentParser(
            description="Mitovar cmd lines",
            usage=""" trackcluster.py <command> [<args>]

-------
The command contains:
pre: prepare the run folder by separating the tracks from the same locus
cluster:run the trackcluster main function to get the isoforms and the counting
desc: compare the novel isoforms with the existing annotations to give an description for the min edit distance between a novel isoform with its nearest reference annotation

test: run the test code for some gene models in the test folder
------
A example for running annotation command:
trackcluster.py pre --reference ref.bed --tracks in.bed --out trackall
trackcluster.py cluster --folder trackall
trackcluster.py clusterj --folder trackall # run in junction mod
trackcluster.py desc --isoform isoform.bed --reference ref.bed > desc.bed 

trackcluster.py test --install

version 0.1.01
            """,
            formatter_class = argparse.RawDescriptionHelpFormatter)
        parser.add_argument("command", help="Subcommand to run")
        args=parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)
        getattr(self, args.command)()


    @staticmethod
    def group_bigg_by_gene(bigglist):
        gene_bigg = OrderedDict()

        for bigg in bigglist:
            try:
                gene_bigg[bigg.geneName].append(bigg)
            except KeyError:
                gene_bigg[bigg.geneName] = []
                gene_bigg[bigg.geneName].append(bigg)
        return gene_bigg


    def pre(self):
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

    def cluster(self):
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
            if is_bin_in("samtools") and is_bin_in("bedtools"):
                logger.info("Pass")
            else:
                logger.info("Check samtools and bedtools installion")

            for package_name in ["pysam", "Bio", "numpy", "pandas"]:
                if is_package_installed(package_name):
                    logger.info("Package {} installed".format(package_name) )
                else:
                    logger.info("Import Error: Package {} not installed".format(package_name) )



        if args.pre: ### test the prepare function using test
           pass


if __name__=="__main__":
    CMD()

