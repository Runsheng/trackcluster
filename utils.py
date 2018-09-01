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
import pybedtools



def wrapper_bedtools_jaccard(bedfile1, bedfile2):

    cmd="bedtools jaccard -a {bedfile1} -b {bedfile2}".format(bedfile1=bedfile1, bedfile2=bedfile2)
    out=myexe(cmd)

    foo=out.split("\n")[1].split("\t")

    jaccard={"intersection":float(foo[0]),
             "union-intersection":float(foo[1]),
             "jaccard":float(foo[2]),
             "n_intersections":float(foo[3])}
    return jaccard

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


def set_tmp():
    """
    set the the tmp dir to the mem
    """
    aa=os.popen("echo $XDG_RUNTIME_DIR")
    dirpath=aa.read().strip()
    pybedtools.set_tempdir(dirpath)

    return dirpath


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