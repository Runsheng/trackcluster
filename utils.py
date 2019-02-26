#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 6:56 PM
# @Author  : Runsheng     
# @File    : utils.py

# std library import
import multiprocessing
import subprocess
import signal
import os
import sys
from collections import OrderedDict

# third part import
from Bio import SeqIO,Seq
import gzip


def count_file(thefile):
    count = 0
    for line in open(thefile).xreadlines(  ):
        count += 1
    return count


def get_name_from_bedfile(bedfile1):
    return bedfile1.split("_")[0].split("/")[-1]


def fasta2dic(fastafile):
    """
    Give a fastq file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large readfile may take a lot of RAM
    """
    if ".gz" in fastafile:
        handle=gzip.open(fastafile, "rU")
    else:
        handle=open(fastafile, "rU")
    record_dict=SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict


def chr_select(record_dict, chro, start,end):
    """
    Note the start and end is 0 based
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
    out, err=proc.communicate()
    #print(err, proc.returncode)
    return out


def set_tmp(wkdir=None):
    """
    set the the tmp dir to the mem
    """
    if wkdir is None:
        aa=os.popen("echo $XDG_RUNTIME_DIR")
        dirpath=aa.read().strip()
        #pybedtools.set_tempdir(dirpath)
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