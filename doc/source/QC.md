# long read QC for RNA reads
The read length and the relative read length/real length (full length ratio) are essential for the downstream analysis. 
For instance, we will not suggest isoform calling with full length ratio lower than 10%. However, the read counting can still
be done.

### Basic statistics

Some of the basic characteristics need to be known before the further analysis. For example, the read length and mapped
read length; the read estimated quality, mapping quality and the read mapping rate. 

We recommend to use the following tools for the long read QC:

Giraffe: https://github.com/lrslab/Giraffe_View

### Full length read ratio estimation

1. Common case, with no 5' indicator. 

The sequencing of most of the long reads are started from the 3' end (PolyA site), so 5' indicator like the splicing leader
or artificial sequence added after de-capping could indicate if one read is likely to be full length. But for most of the
sequencing reads, we do not have this resource. As a result, we will use >95% to estimate the full length ratio.

2. Special case, with 5' indicator.
-   Splicing leaders: Some species using both cis and trans splicing, like C. elegans. The 80% trans-spliced transcripts have the splicing leader sequence at the 5' end of 
mRNA reads. In this case, we can use the splicing leader sequence to estimate the full length ratio. The remained splicing 
leader sequences are ~22nt long short sequences hanging at the 5' end of the reads.
-  Artificial sequence: Some of the long reads are added with artificial sequence after de-capping. 
The artificial sequence could be used to indicate full length read. For example, the Cappble-seq reads are generated by 
decapping the 5' G cap and adding the 5' artificial sequence. Some old fashion ways like 5'RACE will also give the users
some 5' sequences to indicate the start of a transcript.



