#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/8/2018 2:44 PM
# @Author  : Runsheng
# @File    : gff.py

"""
Used to parser the gff gene terms and convert to bigGenePred object
GFF is a collection of one large gff file
"""

# standard library import
from collections import namedtuple, OrderedDict
from operator import attrgetter

# set nametuple for gff class
GFF_FIELD = ['seqid', 'source', 'type', 'start', 'end',
                 'score', 'strand', 'frame', 'attributes']

GFF_record = namedtuple("gff", GFF_FIELD)


def parse_gff_line(line):
    """
    :return: gff_records list, contains 9 col nametuple
    """

    line_l = line.strip().split("\t")

    # test if the line is 9 cols
    if len(line_l) != len(GFF_FIELD):
        print("{0} != {1}".format(len(line_l), len(GFF_FIELD)))
        pass

    # in case some gff have quotes in the cols, replace them
    else:
        data = {
            'seqid': None if line_l[0] == '.' else line_l[0].replace('"', ''),
            'source': None if line_l[1] == '.' else line_l[1].replace('"', ''),
            'type': None if line_l[2] == '.' else line_l[2].replace('"', ''),
            'start': None if line_l[3] == '.' else int(line_l[3]),
            'end': None if line_l[4] == '.' else int(line_l[4]),
            'score': None if line_l[5] == '.' else float(line_l[5]),
            'strand': None if line_l[6] == '.' else line_l[6].replace('"', ''),
            'frame': None if line_l[7] == '.' else line_l[7].replace('"', ''),
            'attributes': parse_attributes(line_l[8])
        }

    return GFF_record(**data)


def parse_attributes(attr_string):
    """
    used to parse the GFF attributes field (col 9) into a dict
    example: ID=ctg7180000006767:hsp:104246:3.10.0.0;Parent=ctg7180000006767:hit:31434:3.10.0.0; \
    Target=CBG07519 47 108;Length=292;Gap=F2 M62 R1
    return: [('ID', 'ctg7180000006767:hsp:104246:3.10.0.0'), ('Parent', 'ctg7180000006767:hit:31434:3.10.0.0'), \
     ('Target', 'CBG07519 47 108'), ('Length', '292'), ('Gap', 'F2 M62 R1')]))
    """
    attr_dic = OrderedDict()

    if attr_string == ".":
        return attr_dic
    else:
        attr_dic = OrderedDict()
        for attribute in attr_string.strip().split(";"):
            if len(attribute) > 0:
                elems = attribute.strip().split('=')
                key = elems[0]
                value = ' '.join(elems[1:])
                if value[0] == '"':
                    value = value[1:]
                if value[-1] == '"':
                    value = value[0:-1]
                attr_dic[key] = value
    return attr_dic


class GFF(object):
    """
    A general class to parse the gff file
    For speed issue, just read all lines into mem
    """

    def __init__(self, gfffile, indicator="sequence_name"):

        # init in file reading
        self.filename=gfffile
        self.gff_string_list = []
        self.gff_to_list() # here use list for the flexibility instead of generator

        self.gff_records=[] # hold multiple GFF_record in a list

        # other
        self.gene_d=None
        self.transcript_to_gene=None
        self.transcript_d=None

        self.indicator=indicator

    def gff_to_list(self):
        """
        :return: gfflines list
        """
        with open(self.filename, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line[0] != "#" and len(line) > 18: # ignore the # and small lines
                    self.gff_string_list.append(line)

    def gene_format(self):
        """
        format the gff line to blocks,
        {gene name, [mRNAline, cdsline...]}

        the gff need to be sorted and all gene/mRNA/exon/cds/UTR is together
        """
        indicator=self.indicator

        gene_l=self.gff_string_list

        gene_d = OrderedDict()
        for line in gene_l:
            line_l = line.split("\t")

            gff_type = line_l[2]
            attributes = line_l[-1].strip()

            if gff_type == "gene":
                for attribute in attributes.split(";"):
                    if indicator in attribute: # the key "sequence_name" only for wormbase gff
                        kk = attribute.split("=")[1]
                        gene_d[kk] = []
                        gene_d[kk].append(line)

            else:
                for attribute in attributes.split(";"):
                    if "ID" in attribute:
                        parent = attribute.split("=")[1]
                        if kk in parent:
                            gene_d[kk].append(line)
                            break
                    if "Parent" in attribute:
                        parent = attribute.split("=")[1]
                        if kk in parent:
                            gene_d[kk].append(line)
                            break
        self.gene_d=gene_d

    def transcript_format(self, keys=None):
        transcript_to_gene={} # {mrna_name: genename}
        transcript_d=OrderedDict() # {mrna_name: gff_line}

        if self.gene_d is None:
            self.gene_format()

        if keys is None: # if no chose ones, use all mrna
            keys=self.gene_d.keys()


        for k in keys:
            gff_list=self.gene_d[k]

            # first test the type of the gene by test if CDS is in the gff_line
            type_set=set([x.split("\t")[2] for x in gff_list])
            is_protein_coding=True if "CDS" in type_set else False
            exon_junction=set()

            for line in gff_list:
                record=parse_gff_line(line)

                if record.type == "exon":
                    try:
                        transcript_name = record.attributes["Parent"].split(":")[1]  # only for wormbase
                    except IndexError:
                        transcript_name = record.attributes["Parent"]  # general
                    transcript_to_gene[transcript_name] = k

                    try:
                        transcript_d[transcript_name].append(record)
                    except KeyError:
                        transcript_d[transcript_name] = []
                        transcript_d[transcript_name].append(record)

                    exon_junction.add(record.start)
                    exon_junction.add(record.end)

            for line in gff_list:
                # exon can be used to record genes from non-coding genes
                if is_protein_coding: # need to add the not in exon cds to the key
                    record=parse_gff_line(line)
                    if record.type=="CDS" and (record.start not in exon_junction) and (record.end not in exon_junction): # cds can have multiple parents
                        parent_l= record.attributes["Parent"].split(",")
                        for parent in parent_l:
                            try:
                                transcript_name=parent.split(":")[1] # only for wormbase
                            except IndexError:
                                transcript_name=parent # only for wormbase
                            #transcript_to_gene[transcript_name]=k
                            try:
                                transcript_d[transcript_name].append(record)
                            except KeyError:
                                transcript_d[transcript_name]=[]
                                transcript_d[transcript_name].append(record)

        for k, v in transcript_d.items():
            v.sort(key=attrgetter('start'))

        self.transcript_to_gene=transcript_to_gene
        self.transcript_d=transcript_d

    def gff_write(self, out, keys=None):

        def write_line():
            fw.write(line)

        fw=open(out,"w")

        if keys is None:
           keys=self.gene_d.keys()

        for k in keys:
            print k
            v=self.gene_d[k]
            for line in v:
                if "gene" in line.split("\t")[2]:
                    write_line()
            for line in v:
                if "mRNA" in line.split("\t")[2]:
                    write_line()
            for line in v:
                if "exon" in line.split("\t")[2]:
                    write_line()
            for line in v:
                if "CDS" in line.split("\t")[2]:
                    write_line()
            for line in v:
                if "UTR" in line.split("\t")[2]:
                    write_line()
        fw.close()




class WormBaseGFF(GFF):
    """
    specific parser for WormBase GFF file, include:
    un-sorted CDS and exon
    mutiple parents:
    mRNA (transcript):Parent=Gene:WBxxxxx
    exon:Parent=Transcript:Y74C9A.3
    CDS: same
    intron: same
    UTR: same

    """
    pass





