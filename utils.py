#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 6:56 PM
# @Author  : Runsheng     
# @File    : utils.py

import time
from functools import wraps
import multiprocessing
import subprocess
import signal
import os
import sys
from collections import OrderedDict
import functools

def count_file(thefile):
    count = 0
    for line in open(thefile).xreadlines(  ):
        count += 1
    return count


def get_name_from_bedfile(bedfile1):
    return bedfile1.split("_")[0].split("/")[-1]


def __wrapper_bedtools_jaccard(bedfile1, bedfile2):

    cmd="bedtools jaccard -a {bedfile1} -b {bedfile2}".format(bedfile1=bedfile1, bedfile2=bedfile2)
    out=myexe(cmd)

    foo=out.split("\n")[1].split("\t")

    jaccard={"intersection":float(foo[0]),
             "union-intersection":float(foo[1]),
             "jaccard":float(foo[2]),
             "n_intersections":float(foo[3])}
    return jaccard


def __wrapper_bedtools_intersection_muti(bigg_list1, bigg_list2, use="exon", core=40):

    pair_list = []

    for i in bigg_list1:
        pair_list.append((i, bigg_list2))

    def run_one(pair):
        bigg_one, bigg_list=pair
        return __wrapper_bedtools_intersect(bigg_one, bigg_list, use)

    out_l = parmap(run_one, pair_list, core)
    # flatten the result
    #out=[]
    #for i in out_l:
    #    for j in i:
    #        out.append(j)
    return out_l


def __wrapper_bedtools_intersect(bigg_one, bigg_list, use="exon"):

    # generate the bedfile
    bigg_one.to_bedfile()
    for i in bigg_list:
        i.to_bedfile()

    if use== "exon":
        bedfile1=bigg_one.exon_file
        bedfile_list=[x.exon_file for x in bigg_list]
    elif use== "intron":
        # intron could be 0
        if count_file(bigg_one.intron_file)==0:
            return []
        else:
            bedfile1=bigg_one.intron_file
            bedfile_list=[x.intron_file for x in bigg_list]

    bedfile_str=" ".join(bedfile_list)

    cmd="bedtools intersect -wa -wb -a {bedfile1} -b {bedfile2}".format(bedfile1=bedfile1, bedfile2=bedfile_str)

    out=myexe(cmd)
    #print cmd
    out_l=out.strip().split("\n")

    i_d=[]
    name_bed1=get_name_from_bedfile(bedfile1)
    target_intersection=OrderedDict()

    for i in bedfile_list:
        name_bed2=get_name_from_bedfile(i)
        target_intersection[name_bed2]=0

    if out_l[0]!="":
        for line in out_l:
            line_l=line.split("\t")
            #print line_l
            if len(line_l)==7:
                _,q_start,q_end,t_name,_,t_start,t_end=line_l
                intersection=min(int(q_end), int(t_end))-max(int(q_start), int(t_start))
                name_bed2=get_name_from_bedfile(t_name)
                target_intersection[name_bed2]+=intersection
            # bedtools -filenames do not show the name of bed2 if bed2 is only one file
            elif len(line_l)==6: # just one in bed2
                _,q_start,q_end,_,t_start,t_end=line_l
                intersection=min(int(q_end), int(t_end))-max(int(q_start), int(t_start))
                name_bed2=get_name_from_bedfile(bedfile_list[0])
                target_intersection[name_bed2]+=intersection

    else: # all intersection is 0
        pass

    for name_bed2,intersection in target_intersection.items():
        i_d.append((name_bed1, name_bed2, intersection))
    return i_d


def myexe(cmd, timeout=0):
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