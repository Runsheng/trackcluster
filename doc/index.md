# TrackCluster
![PyPI](https://img.shields.io/pypi/v/trackcluster?color=green)

Trackcluster is an isoform calling and quantification pipeline for long RNA/cDNA reads.


Walkthrough for the using of trackcluster in isofrom calling in worms/mouse/human/virus datasets.

## preprocessing 
1. Long read QC

Some of the basic characteristics need to be known before the further analysis. For example, the read length and mapped
read length; the read estimated quality, mapping quality and the read mapping rate. 

We recommend to use the following tools for the long read QC:

Giraffe: https://github.com/lrslab/Giraffe_View


2. Mapping Nanopore long reads to genome

The mapping process is the first step for the isoform calling. 
We recommend to use minimap2 for the mapping of the long reads to the genome.
    - Mapping the reads for eukaryotic genomes
    ```bash
    minimap2 -ax splice -uf -k14 --secondary=no --MD -t 8 ref.fa read.fq.gz | samtools view -bS - | samtools sort -o read.bam -
    ```
