# nanopore-readmapper
Readmapper for RNA secondary structure determination data containing full length reads of RNA. 

# Usage

## read_mapper_2.py
``` shell
python3 read_mapper_2.py ref.fasta reads.fastq output.sam -m 1 -x 0 -o (-2) -e (-1) 
```
This script takes a reference genome in fasta format and a fastq file with reads and creates a SAM file with the alignments of the reads to the reference genome. The SAM file is outputted to the file specified by the user. The SAM file can then be sorted and indexed using samtools.


## SAM_creator.py
Contains the class SAM (complete alignment in SAM format) and SAMline (single read alignment) which is used in read_mapper_2.py to create a SAM file from a given fastq file with reads and a given reference genome in fasta format. The SAM file is created by aligning the reads to the reference genome using the Biopython Align package.

