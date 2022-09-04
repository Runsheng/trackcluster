#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 6:56 PM
# @Author  : Runsheng     
# @File    : utils.py

# std library import
import itertools
import multiprocessing
import subprocess
import signal
import os
import sys
import logging
from collections import OrderedDict

# third part import
from Bio import SeqIO,Seq
import gzip


def count_file(thefile):
    """
    count bed file line number have to have at least 4 column (chro, start, end)
    :param thefile:
    :return:
    """
    count = 0
    f=open(thefile, "r")
    for line in f.readlines():
        try:
            if len(line.split("\t"))>=3:
                count += 1
        except Exception:
            pass
    f.close()
    return count


def get_name_from_bedfile(bedfile1):
    return bedfile1.split("_")[0].split("/")[-1]


def fasta2dic(fastafile):
    """
    Give a fastq file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large readfile may take a lot of RAM
    """
    if ".gz" in fastafile:
        handle=gzip.open(fastafile, "r")
    else:
        handle=open(fastafile, "r")
    record_dict=SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict


def chr_select(record_dict, chro, start,end):
    """
    Note the start and end is 0 based and [)
    give the name of refdic, and the chr, start and end to be used
    return the name and sequence (both as str)
    for example, chrcut(record_dict, "I", 100,109) returns
     ("I:100_109","AAAAAAAAAA")
    """
    name=chro+ ":"+str(start)+"_"+str(end)
    seq=str(record_dict[chro][start:end].seq)
    return name,seq


def myexe(cmd, timeout=10):
    """
    a simple wrap of the shell
    mainly used to run the bwa mem mapping and samtool orders
    """
    def setupAlarm():
        signal.signal(signal.SIGALRM, alarmHandler)
        signal.alarm(timeout)

    def alarmHandler(signum, frame):
        sys.exit(1)

    proc=subprocess.Popen(cmd, shell=True, preexec_fn=setupAlarm,
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=os.getcwd())
    #print("Running:  ",  cmd)
    out, err=proc.communicate()
    #print("stderr out:  ",err, "\n","return code:  ",proc.returncode)
    return out


def is_bin_in(cmd_name):
    """
    used to test if bedtools and samtools are in the bin
    :param cmd_name:
    :return:
    """
    out=myexe(cmd_name+" --version")
    # could return name or just version number with 2.1.1
    if cmd_name in str(out) or "." in str(out): # py3, out is byte
        return True
    else:
        return False


def is_package_installed(package_name):
    try:
        import package_name
    except ImportError:
        return False
    return True


def set_tmp(wkdir=None):
    """
    set the the tmp dir to the mem
    """
    if wkdir is None:
        aa=os.popen("echo $XDG_RUNTIME_DIR")
        dirpath=aa.read().strip()
        aa.close()
        return dirpath
    else:
        return wkdir


def reverse_complement(seq):
    """
        Given: A DNA string s of length at most 1000 bp.
        Return: The reverse complement sc of s.
        due to the complement_map,
        the symbol such as \n and something else is illegal
        the input need to be pure sequence
    """
    return Seq.reverse_complement(seq)


def del_files(filename_l):
    for filename in filename_l:
        try:
            os.remove(filename)
        except OSError:
            pass


def parmap(f, X, nprocs=multiprocessing.cpu_count()):
    """
    a function to use muti map inside a function
    modified from stackoverflow, 3288595
    :param f:
    :param X:
    :param nprocs: core, if not given, use all core
    :return:
    """
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
            for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]


def fun(f, q_in, q_out):
    """
    for parmap
    :param f:
    :param q_in:
    :param q_out:
    :return:
    """
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))


#### add two logger to reduece the use of print
def log_summary():
    """
    # define only the summary logger first
    """
    logger = logging.getLogger('summary')
    logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def log_detail_file(filename):
    """
    # define only the summary logger first
    """
    logger = logging.getLogger('details')
    logger.setLevel(logging.DEBUG)

    ch = logging.FileHandler(filename=filename)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def group_site(missed_order):
    # sanity check:
    if len(missed_order)%2!=0:
        pass
        # could have some exception: 1. include the 0 or the -1 junction 2. intron retain coupled with missed exon
        # not bug, # print "Error in the number of orders !", missed_order

    # group the lists to regions
    groups = [[y[1] for y in g] for k, g in itertools.groupby(enumerate(missed_order), key=lambda x: x[0] - x[1])]

    return groups


def __group_nearby_site(site_list, interval=5):
    """
    get list for all site with nearby
    :param site_list:
    :param interval:
    :return:
    """
    # group the lists to regions
    groups = [[y[1] for y in g] for k, g in itertools.groupby(enumerate(site_list), key=lambda x: abs(x[1] - x[0])<=interval ) ]

    return groups


def get_file_prefix(filepath, sep="_"):
    return filepath.split("/")[-1].split(sep)[0]


def get_file_location(filepath):
    return "/".join(filepath.split("/")[0:-1])

def list2file(ll, genename_file):
    """
    IO function
    write a list to a file
    :param ll:
    :param genename_file:
    :return:
    """
    with open(genename_file, "w") as fw:
       for gene in ll:
           fw.write(gene)
           fw.write("\n")
    return genename_file

def file2list(genename_file):
    """
    IO function
    reverse of file2list
    :param genename_file:
    :return:
    """
    gene_l=[]
    with open(genename_file, "r") as f:
        for line in f.readlines():
            gene_l.append(line.strip())
    return gene_l

def print_dic(dic, n=100):
    """
    test code, print the first 100 line of dic
    """
    for k in list(dic.keys())[:n]:
        print(k, dic[k])