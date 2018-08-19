#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/18/2018 1:48 PM
# @Author  : Runsheng     
# @File    : ssw.py


from LR_toolkit.ssw.ssw_wrap import Aligner


def ssw_wrapper(seq1, seq2, match=2, mismatch=2,
                gap_open=3, gap_extend=1):
    """
    parameter are write inside the function
    seq1 is ref and seq2 is query
    # todo : leave a api to change matrix
    """

    seq1_str = str(seq1)
    seq2_str = str(seq2)
    # reduce the gap open score from 3 to 1 for nanopore reads
    aligner = Aligner(seq1,
                      match, mismatch, gap_open, gap_extend,
                      report_cigar=True, )
    aln = aligner.align(seq2)  # min_score=20, min_len=10)

    return aln


# copy the library from primer_RFLP(runsheng, 2016 version)
import re


def consensus(query, ref, cigar_string, ref_shift=0):
    """
    to get a concensus sequence as the longest cons of both sequence,
    and fill all disagreement with N
    for example,
    query:AAATA-TAGAA
    ref:  AAACACTA-AA
    cons: AAANANCANAA

    use query to get the con
    then get two chain file, 0 based
    """

    pattern = re.compile('([0-9]*)([MIDNSHP=X])')

    query_out = []
    ref_out = []
    cons = []

    query_pos = 0
    ref_pos = ref_shift

    for length, code in pattern.findall(cigar_string):
        length = int(length)

        if code == "M":
            for i in range(length):
                q_s = query[query_pos].upper()
                r_s = ref[ref_pos].upper()

                query_out.append(q_s)
                ref_out.append(r_s)

                cons.append(q_s) if q_s == r_s else cons.append("N")

                ref_pos += 1
                query_pos += 1

        elif code in "D":
            for i in range(length):
                r_s = ref[ref_pos]

                ref_out.append(r_s)
                query_out.append("-")
                cons.append("N")

                ref_pos += 1

        elif code in "IHS":
            for i in range(length):
                q_s = query[query_pos]

                query_out.append(q_s)
                ref_out.append("-")
                cons.append("N")
                query_pos += 1

    return "".join(cons), "".join(query_out), "".join(ref_out)


def get_chain(seq_out):
    """
    The postion with "-" as a gap have to be used to change the cord
    The postion is so called "chain" file for one sequence
    """
    return [i for i, nucl in enumerate(seq_out) if nucl == "-"]