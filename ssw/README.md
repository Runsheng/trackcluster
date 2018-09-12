# This is a fork of SSW library. And only the python wrapper is retained. 
### The build-in python wrapper for SSW library is forked from pyDNA package [https://github.com/a-slide/pyDNA]
### Please refer to [https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library] for details.


##Python interface

###How to use the python wrapper ssw_wrap.py

A Python wrapper partially implements the c library functionality (only DNA sequences can be aligned for now). c libraries are completely integrated in a simple module that do not require any C programming knowledge to be used.
Briefly, An aligner object can be initialized with alignment parameters and a reference subject sequence, then the object method *align* can be called with a query sequence and filtering conditions (min score and min length) as many time as desired.
Depending of the score and length requested an python object PyAlignRes will be eventually returned.

To use the python wrapper, please:

1. Compile the src folder by either using the makefile or by compiling a dynamic shared library with gcc ```gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h```
2. libssw.so and ssw_wrap.py can them be put in the same folder of your own program files.
3. Depending of the LINUX OS version installed it may be required to modify the LD_LIBRARY_PATH environment variable to use the dynamic library libssw.so by one of the 2 following possibilities :

    * Export the path or the directory containing the library ```LD_LIBRARY_PATH=path_of_the_library```
    * For a definitive inclusion edit /etc/ld.so.conf and add the path of the lib directory. Then, update the cache by using ```/sbin/ldconfig``` as root
4. In a python script or in a interactive interpreter the main class can be imported with : ```from ssw_wrap import Aligner```
5. Instantiate the Aligner class with initial parameters, including the reference subject sequence. ```Aligner(myref, match=2, mismatch=2, gap_open=3, gap_extension=1, report_secondary=False, report_cigar=False)```
6. Call the object align method with a query sequence, the minimal score and length for the alignment to be reported ```res = ssw.align(myquery, min_score=10, min_len=20)```
7. Parse the returned PyAlignRes object for alignment result description 

###Run pyssw standalone 
```
Usage: pyssw.py -s subject.fasta -q fastq (or fasta) [Facultative options]

Options:
 --version             show program's version number and exit
 -h, --help            show this help message and exit
 -s SUBJECT, --subject=SUBJECT Path of the fasta file containing the subject genome sequence [REQUIRED]
 -q QUERY, --query=QUERY   Path of the fastq or fasta file containing the short read to be aligned [REQUIRED]
 -t QTYPE, --qtype=QTYPE   Type of the query file = fastq or fasta. [default:fastq]
 -m MATCH, --match=MATCH   positive integer for weight match in genome sequence alignment. [default: 2]
 -x MISMATCH, --mismatch=MISMATCH  positive integer. The negative value will be used as weight mismatch in genome sequence alignment.[default: 2]
 -o GAP_OPEN, --gap_open=GAP_OPEN  positive integer. The negative value will be used as weight for the gap opening. [default: 3]
 -e GAP_EXTEND, --gap_extend=GAP_EXTEND    positive integer. The negative value will be used as weight for the gap opening. [default: 1]
 -f MIN_SCORE, --min_score=MIN_SCORE   integer. Consider alignments having a score <= MIN_SCORE as not aligned. [default: 0]
 -l MIN_LEN, --min_len=MIN_LEN integer. Consider alignments having a length <= as not aligned. [default: 0]
 -r, --reverse         Flag. Align query in forward and reverse orientation and choose the best alignment. [Set by default]
 -u, --unaligned       Flag. Write unaligned reads in sam output [Unset by default]
```

###pyssw output

For now, the program only output a SAM like file. The file encoding is sightly more complete than the version created by the C/C++ software and is fully compatible with downstream NGS utilities such as samtools. The MAPQ score field is set to 0, but the ssw score is reported in the AS optional tag. 

```
@HD	VN:1.0	SO:unsorted
@SQ	SN:Virus_genome	LN:2216
@PG	ID:Striped-Smith-Waterman	PN:pyssw	VN:0.1
@CO	Score_values = match 2, mismatch 2, gap_open 3, gap_extend 1
@CO	Filter Options = min_score 250, min_len 0
M00785:2:000000000-A60JC:1:1101:14942:1516	16	Virus_genome	1773	0	151M	*	0	0    TTGTTT...TATTAC >A1AAF...@1@@2@	    AS:i:302
M00785:2:000000000-A60JC:1:1101:13644:1543	0	Virus_genome	1985	0	21M1I129M	*	0	0    TTTAAA...CTCGCT	AAAA11.../<</</     AS:i:289
M00785:2:000000000-A60JC:1:1101:13883:1547	16	Virus_genome	716	0	36M2D115M	*	0	0	CTTACC...GGAATT	>AA1>C...////?F    AS:i:294
...
```
###Speed and memory usage of pyssw

Intel Core i5 3320PM CPU @ 2.60GHz x 2

* E. coli K12 U0096.3 (4,641,652 nucleotides) vs 10,000 illumina reads (50pb) from [CLoVR public repository](http://data.clovr.org/d/17/e-coli-illumina-inputs)

    ```pyssw.py -s Ecoli_K12_U0096.3.fa -q Ecoli_illumina_10k.fastq -r -f 80```

    Time : ~55m 40s

* Short virus reference (2,216 nt) vs 100,000 illumina Miseq reads (150pb) (see demo dataset)

    ```pyssw.py -s Virus_genome.fa.gz -q 100k_illumina1.fastq.gz -f 250```

    Time : ~60s

