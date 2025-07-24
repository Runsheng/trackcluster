# The protocol to run the fusion detection in read/isoform levels.

## Fusion read types:
1. In cancer reasearch, the fusion read are generally generated from the fusion of genomes. Two independent part of 
genome has been linked together, and the exons are cated together to form a long fusion read. The gene fusion will cause 
a lot of consequences in cancer research.
2. The other type of fusion is the readthrough of one long read, spanning adjudictive genes. This fusion would indicate 
a possible operon in _C. elegans_ or _E. coli_ genome. And we could infer the fusion events by directly using reads or 
using the high expressed isoforms. 


## the usage in C. elegans genome
In C. elegans, fusion reads can be used to directly for the detection of operon genes.But we need to apply different
filters before we can get reliable results. 

### filters need to be applied
1. The genes assigned is not close to each other, which is a long range fusion, which could be false positive. Even the 
case is true positive, this type of fusion is not likely to indicate an operon.
2. Gene direction filter, the operon holds genes in one direction.
3. Fusion length cutoff, need to ensure at least one exon is full contained, or maybe for single exon gene, >50%
of the later gene's read length is covered in the read.


#####test data
```bash
-rw-rw-r--  1 li   li   1021885460 Apr  8 16:28 reads_sorted.bed
-rw-rw-r--  1 li   li     11610085 Apr  8 16:34 refs.bed
```

```bash
trackrun.py addgene -r ref.bed -s reads.bed -f1 0.1 -f2 0.1
 
# will generate reads_gene.bed
trackrun.py desc --isoform reads_gene.bed --reference ref.bed # will generated reads_desc.txt and reads_class12.txt 


```