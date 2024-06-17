# TrackCluster
![PyPI](https://img.shields.io/pypi/v/trackcluster?color=green)

Trackcluster is an isoform calling and quantification pipeline for long RNA/cDNA reads.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Scripts](#scripts)
- [Walkthrough](#walkthrough)

####
The ongoing development can be found in the ["dev"](https://github.com/Runsheng/trackcluster/tree/dev) branch.
#### TODO: 
1. fix "fusion" classification. 2. speed for clusterj. 3. splicing leader/5' indicator finding using pyssw. 
4.add function to get the CDS start (maybe ATG) and CDS end.   

## <a name="overview"></a>Overview
A pipeline for reference-based isoform identification and quantification using long reads. This pipeline was designed to use **only** long and nosisy reads to make a valid transcriptome. An indicator for the intact 5' could be very helpful to the pipeline, i.e, the splicing leader in the mRNA of nematodes. 

The major input/output for this pipeline is "bigg"--["bigGenePred"](https://github.com/Runsheng/trackcluster/blob/master/test/bigGenePred.as) format. 

## <a name="requirements"></a>Requirements

1. developed on python 3.9, tested on python 3.6 and above (or 2.7.10+), should work with most of the py3 versions
2. samtools V2.0+ , bedtools V2.24+  and minimap2 V2.24+ in your $PATH
```bash
# install the external bins with conda
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c bioconda minimap2
```

## Installation
```bash
# use pip from pypi
pip install trackcluster
# or pip from source code for the latest version
git clone https://github.com/Runsheng/trackcluster.git
pip install ./trackcluster
```

## Recommendations
1. UCSC Kent source tree (for generating binary track), used only in bigg2b.py

## Scripts
All scripts can be run directly from shell after pip installation.
- **trackrun.py**: the main script for trackcluser run
- **bam2bigg.py**: convert the mapped read from the bam file, to bigg track format
- **gff2bigg.py**: convert the isoform annotation in gff3 to bigg format 
- bigg2b.py: convert the bigg track into binary format for better loading in IGV/UCSC
- biggmutant.py: change the value of one column in a bigglist

## <a name="walkthrough"></a>Walkthrough
```bash
# test if all dependencies are installed
trackrun.py test --install

# prepare the reference annotation bed file from gff file
# tested on Ensembl, WormBase and Arapost gff
gff2bigg.py -i ensemblxxxx.gff3 -o ref.bed 
# WormBase full gff contains too many information, need to extract the lines from WormBase only
cat c_elegans.PRJNA13758.WS266.annotations.gff3 |grep WormBase > ws266.gff
gff2bigg.py -i ws266.gff -o ref.bed
# the ref.bed can be sorted to speed up the analysis
bedtools sort -i ref.bed > refs.bed # refs.bed contains the sorted, know transcripts from gff annotation

# generate the read track from minimap2 bam file
bam2bigg.py -b group1.bam -o group1.bed
bam2bigg.py -b group2.bam -o group2.bed

# merge the bed file and sort
cat group1.bed group2.bed > read.bed
bedtools sort -i read.bed > reads.bed

# Examples for running commands:
trackrun.py clusterj -s reads.bed -r refs.bed -t 40 # run in junction mode, will generate the isoform.bed
trackrun.py count -s reads.bed -r refs.bed -i isoform.bed # generate the csv file for isoform expression
# alternative for cluster
trackrun.py cluster -s reads.bed -r refs.bed -t 40 # run in exon/intron intersection mode， slower, will generate the isoform.bed

# the post analysis could include the classification of novel isoforms
trackrun.py desc --isoform isoform.bed --reference ref.bed # generate the description for each novel isoform
# this part can be run directly on reads, to count the frequency of splicing events in reads, like intron_retention
trackrun.py addgene -r ref.bed -s reads.bed # will generate reads_gene.bed
trackrun.py desc --isoform reads_gene.bed --reference ref.bed # will generated reads_desc.txt and reads_class12.txt 

```


## Citation
Please kindly cite our paper for using trackcluster in your work.

Li, R., Ren, X., Ding, Q., Bi, Y., Xie, D. and Zhao, Z., 2020. Direct full-length RNA sequencing reveals unexpected transcriptome complexity during *Caenorhabditis elegans* development. **Genome research**, 30(2), pp.287-298.