#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/30/2018 6:05 PM
# @Author  : Runsheng     
# @File    : pyrange_test.py


from pyranges.pyranges import PyRanges

a_start=[10, 30, 50, 60, 70, 80]
a_end=[20, 40, 60,70,80,90]

#a_start=[14647320, 14647919, 14648551, 14648885, 14649494, 14649843, 14650060, 14650325, 14650593, 14651002, 14652682, 14653660, 14654218, 14654570, 14656089, 14656327, 14657212, 14657637, 14658474]
#a_end=[14647857, 14648142, 14648827, 14649435, 14649798, 14650013, 14650279, 14650538, 14650898, 14651210, 14652809, 14653815, 14654458, 14654777, 14656270, 14657097, 14657491, 14657919, 14658715]

#b_start=[5,32,60]
#b_end=[11,41,70]

b_start=[1, 5]
b_end=[4, 7]

a=PyRanges(seqnames="1",  starts=a_start, ends=a_end, strands="+")
b=PyRanges(seqnames="1", starts=b_start, ends=b_end, strands="+")


try:
    c=a.jaccard(b, False)
except ValueError:
    c={"jacard": 0}
#print c
