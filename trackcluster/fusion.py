#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/4/2025 4:45
# @Author  : Runsheng
# @File    : fusion.py

"""
add a single file to handle the fusion/readthrough reads, using different methods and filters

The first aim is to locate the operons in C. elegans genome with fusion reads
filters: >=2 or more reads support
"""

import os
import sys
import argparse
from collections import Counter, defaultdict
from functools import partial

# Import from trackcluster package
from trackcluster.flow import flow_add_gene
from trackcluster.tracklist import read_bigg, write_bigg
from trackcluster.utils import get_file_prefix


def flow_fusion(wkdir, prefix, bigg_gff_file, bigg_nano_file, f1=0.1, f2=0.1, min_support=2):
    """
    The function used to capture the fusion reads which can indicate operons (readthrough)
    :param wkdir: working directory
    :param prefix: output file prefix
    :param bigg_gff_file: reference annotation file in bigg format
    :param bigg_nano_file: nanopore reads file in bigg format
    :param f1: minimum intersection fraction for read track
    :param f2: minimum intersection fraction for isoform track
    :param min_support: minimum number of reads supporting a fusion event
    :return: dictionary of fusion events and their supporting reads
    """
    os.chdir(wkdir)
    print(f"Processing fusion reads in {wkdir}")
    print(f"Parameters: f1={f1}, f2={f2}, min_support={min_support}")

    # Add gene annotations to reads
    bigg_new = flow_add_gene(wkdir, prefix, bigg_gff_file, bigg_nano_file, f1=f1, f2=f2)
    print(f"Processed {len(bigg_new)} reads with gene annotations")

    # Identify fusion reads by checking gene names
    fusion_d = {}
    fusion_gene_pairs = defaultdict(list)
    
    for bigg in bigg_new:
        genename_l = bigg.geneName.split("||")
        if len(genename_l) > 1:
            # This is a fusion read spanning multiple genes
            fusion_d[bigg.name] = genename_l
            # Create sorted gene pair for counting
            gene_pair = "||".join(sorted(genename_l))
            fusion_gene_pairs[gene_pair].append(bigg.name)

    print(f"Found {len(fusion_d.keys())} potential fusion reads")

    # Write individual fusion reads
    fw_name_l = [prefix, str(f1), str(f2), "fusion.txt"]
    fw_name = "_".join(fw_name_l)

    with open(fw_name, "w") as fw:
        fw.write("read_name\tgenes\tgene_count\n")
        for k, v in fusion_d.items():
            v_str = ";".join(v)
            fw.write(f"{k}\t{v_str}\t{len(v)}\n")

    # Count fusion events and filter by minimum support
    gene_dic = parser_fusion(fw_name)
    filtered_fusions = {k: v for k, v in gene_dic.items() if v >= min_support}
    
    # Write fusion event counts
    fw_name_l = [prefix, str(f1), str(f2), "count.txt"]
    fw_name = "_".join(fw_name_l)
    with open(fw_name, "w") as fw:
        fw.write("gene_combination\tread_count\n")
        for k, v in gene_dic.items():
            fw.write(f"{k}\t{v}\n")

    # Write filtered high-confidence fusion events
    fw_name_l = [prefix, str(f1), str(f2), f"filtered_min{min_support}.txt"]
    fw_name = "_".join(fw_name_l)
    with open(fw_name, "w") as fw:
        fw.write("gene_combination\tread_count\tsupporting_reads\n")
        for gene_pair, read_count in filtered_fusions.items():
            supporting_reads = ";".join(fusion_gene_pairs[gene_pair])
            fw.write(f"{gene_pair}\t{read_count}\t{supporting_reads}\n")

    print(f"Found {len(filtered_fusions)} high-confidence fusion events (>={min_support} reads)")
    
    return fusion_d, filtered_fusions


def parser_fusion(file_fusion, sep="\t"):
    """
    Parse fusion file and count gene combinations
    :param file_fusion: fusion file path
    :param sep: separator character
    :return: dictionary of gene combinations and their counts
    """
    gene_line = []
    with open(file_fusion, "r") as f:
        # Skip header
        next(f)
        for line in f:
            line_l = line.strip().split(sep)
            if len(line_l) >= 2:
                gene_line.append(line_l[1])
    
    # Count gene combination occurrences
    gene_dic = dict(Counter(gene_line))
    return gene_dic


def get_fusion_from_bigg(bigg_list, prefix, min_support=2):
    """
    Get fusion reads from a list of bigg objects
    :param bigg_list: list of bigGenePred objects
    :param prefix: output file prefix
    :param min_support: minimum number of reads supporting a fusion event
    :return: dictionary of fusion reads and filtered high-confidence events
    """
    fusion_d = {}
    fusion_gene_pairs = defaultdict(list)
    
    for bigg in bigg_list:
        genename_l = bigg.geneName.split("||")
        if len(genename_l) > 1:
            fusion_d[bigg.name] = genename_l
            # Create sorted gene pair for counting
            gene_pair = "||".join(sorted(genename_l))
            fusion_gene_pairs[gene_pair].append(bigg.name)

    # Count and filter fusion events
    gene_combinations = [";".join(genes) for genes in fusion_d.values()]
    gene_counts = Counter(gene_combinations)
    filtered_fusions = {k: v for k, v in gene_counts.items() if v >= min_support}

    # Write fusion results
    fusion_file = f"{prefix}_fusion.txt"
    with open(fusion_file, "w") as fw:
        fw.write("read_name\tgenes\tgene_count\n")
        for k, v in fusion_d.items():
            v_str = ";".join(v)
            fw.write(f"{k}\t{v_str}\t{len(v)}\n")

    # Write high-confidence fusion events
    filtered_file = f"{prefix}_filtered_fusion.txt"
    with open(filtered_file, "w") as fw:
        fw.write("gene_combination\tread_count\tsupporting_reads\n")
        for gene_comb, count in filtered_fusions.items():
            # Find supporting reads for this combination
            supporting_reads = []
            for read_name, genes in fusion_d.items():
                if ";".join(genes) == gene_comb:
                    supporting_reads.append(read_name)
            fw.write(f"{gene_comb}\t{count}\t{';'.join(supporting_reads)}\n")

    return fusion_d, filtered_fusions


def analyze_fusion_patterns(fusion_dict, output_prefix):
    """
    Analyze patterns in fusion events
    :param fusion_dict: dictionary of fusion events
    :param output_prefix: prefix for output files
    """
    # Analyze gene pair frequencies
    gene_pairs = []
    gene_frequencies = defaultdict(int)
    
    for read_name, genes in fusion_dict.items():
        # Count individual gene frequencies in fusions
        for gene in genes:
            gene_frequencies[gene] += 1
        
        # Create all possible gene pairs
        if len(genes) >= 2:
            for i in range(len(genes)):
                for j in range(i+1, len(genes)):
                    pair = tuple(sorted([genes[i], genes[j]]))
                    gene_pairs.append(pair)
    
    pair_counts = Counter(gene_pairs)
    
    # Write gene frequency analysis
    freq_file = f"{output_prefix}_gene_frequencies.txt"
    with open(freq_file, "w") as fw:
        fw.write("gene\tfusion_frequency\n")
        for gene, freq in sorted(gene_frequencies.items(), key=lambda x: x[1], reverse=True):
            fw.write(f"{gene}\t{freq}\n")
    
    # Write gene pair analysis
    pair_file = f"{output_prefix}_gene_pairs.txt"
    with open(pair_file, "w") as fw:
        fw.write("gene1\tgene2\tpair_frequency\n")
        for (gene1, gene2), freq in sorted(pair_counts.items(), key=lambda x: x[1], reverse=True):
            fw.write(f"{gene1}\t{gene2}\t{freq}\n")
    
    print(f"Analysis complete. Found {len(gene_frequencies)} genes in fusion events")
    print(f"Found {len(pair_counts)} unique gene pairs")
    
    return gene_frequencies, pair_counts


def flow_fusion_analysis(wkdir, prefix, bigg_gff_file, bigg_nano_file, 
                        f1=0.1, f2=0.1, min_support=2, analyze_patterns=True):
    """
    Complete fusion analysis workflow
    :param wkdir: working directory
    :param prefix: output file prefix
    :param bigg_gff_file: reference annotation file
    :param bigg_nano_file: nanopore reads file
    :param f1: minimum intersection fraction for read track
    :param f2: minimum intersection fraction for isoform track
    :param min_support: minimum number of reads supporting a fusion event
    :param analyze_patterns: whether to perform pattern analysis
    :return: tuple of (fusion_dict, filtered_fusions, gene_frequencies, pair_counts)
    """
    print("=== Fusion Analysis Workflow ===")
    print(f"Working directory: {wkdir}")
    print(f"Output prefix: {prefix}")
    print(f"Reference annotation: {bigg_gff_file}")
    print(f"Nanopore reads: {bigg_nano_file}")
    print(f"Parameters: f1={f1}, f2={f2}, min_support={min_support}")
    print()
    
    # Run fusion detection
    fusion_dict, filtered_fusions = flow_fusion(
        wkdir=wkdir,
        prefix=prefix,
        bigg_gff_file=bigg_gff_file,
        bigg_nano_file=bigg_nano_file,
        f1=f1,
        f2=f2,
        min_support=min_support
    )
    
    gene_frequencies = {}
    pair_counts = {}
    
    # Perform pattern analysis if requested
    if analyze_patterns and fusion_dict:
        print("\n=== Performing Pattern Analysis ===")
        gene_frequencies, pair_counts = analyze_fusion_patterns(fusion_dict, prefix)
    
    # Summary
    print("\n=== Summary ===")
    print(f"Total fusion reads: {len(fusion_dict)}")
    print(f"High-confidence fusion events: {len(filtered_fusions)}")
    
    if filtered_fusions:
        print("\nTop fusion events:")
        sorted_fusions = sorted(filtered_fusions.items(), key=lambda x: x[1], reverse=True)
        for i, (genes, count) in enumerate(sorted_fusions[:10]):
            print(f"  {i+1}. {genes}: {count} reads")
    
    return fusion_dict, filtered_fusions, gene_frequencies, pair_counts


def main():
    """
    Main function for fusion read detection script
    """
    parser = argparse.ArgumentParser(
        description="Detect fusion/readthrough events in long read RNA sequencing data"
    )
    parser.add_argument("-d", "--workdir", default=os.getcwd(),
                       help="Working directory (default: current directory)")
    parser.add_argument("-p", "--prefix", required=True,
                       help="Output file prefix")
    parser.add_argument("-g", "--gff", required=True,
                       help="Reference annotation file in bigg format")
    parser.add_argument("-n", "--nano", required=True,
                       help="Nanopore reads file in bigg format")
    parser.add_argument("-f1", "--fraction1", type=float, default=0.1,
                       help="Minimum intersection fraction for read track (default: 0.1)")
    parser.add_argument("-f2", "--fraction2", type=float, default=0.1,
                       help="Minimum intersection fraction for isoform track (default: 0.1)")
    parser.add_argument("-m", "--min-support", type=int, default=2,
                       help="Minimum number of reads supporting a fusion event (default: 2)")
    parser.add_argument("--analyze", action="store_true",
                       help="Perform additional pattern analysis")
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.gff):
        print(f"Error: GFF file {args.gff} not found")
        sys.exit(1)
    if not os.path.exists(args.nano):
        print(f"Error: Nanopore file {args.nano} not found")
        sys.exit(1)
    
    # Run fusion analysis
    fusion_dict, filtered_fusions, gene_frequencies, pair_counts = flow_fusion_analysis(
        wkdir=args.workdir,
        prefix=args.prefix,
        bigg_gff_file=args.gff,
        bigg_nano_file=args.nano,
        f1=args.fraction1,
        f2=args.fraction2,
        min_support=args.min_support,
        analyze_patterns=args.analyze
    )
    
    print(f"\nOutput files written with prefix: {args.prefix}")


if __name__ == "__main__":
    main()


