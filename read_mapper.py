from Bio import Align
from Bio import SeqIO
import argparse

argparser = argparse.ArgumentParser(description='Read mapper')
argparser.add_argument('fasta', help='fasta file')
argparser.add_argument('fastq', help='fastq file')
argparser.add_argument('sam', help='sam file')
args = argparser.parse_args()

fasta = args.fasta
fastq = args.fastq
sam_file = args.sam

reference = SeqIO.read(fasta, "fasta")
reads = list(SeqIO.parse(fastq, "fastq"))

def makeCIGAR(alignment):
    cigar = ""
    for i in range(len(alignment[0])):
        if alignment[0][i] == "-":
            cigar += "I"
        elif alignment[1][i] == "-":
            cigar += "D"
        else:
            cigar += "M"
    return cigar

def runLengthEncoding(cigar):
    cigar_rle = ""
    i = 0
    while i < len(cigar):
        count = 1
        j = i
        while j < len(cigar)-1 and cigar[j] == cigar[j+1]:
            count += 1
            j += 1
        cigar_rle += str(count) + cigar[i]
        i += count
    return cigar_rle

def runLengthEncodedCIGAR(alignment):
    cigar = makeCIGAR(alignment)
    cigar_rle = runLengthEncoding(cigar)
    return cigar_rle

def makeHeader(alignment):
    header = ""
    header += "@HD" + "\t" + "VN:1.6" + "\t" + "SO:unsorted" + "\n"
    header += "@SQ" + "\t" + "SN:" + alignment.target.name + "\t" + "LN:" + str(len(alignment.target.seq)) + "\n"
    return header

def makeSAM(alignment):
    READ_SEQ=alignment.query.seq
    READ_NAME=alignment.query.name
    READ_ID=alignment.query.id
    SEQ_LEN=len(alignment.query.seq)
    POS_count = 0
    CIGAR=runLengthEncodedCIGAR(alignment)
    for c in CIGAR:
        if c.isdigit():
            POS_count += 1
        elif c == "D":
            POS = CIGAR[0:POS_count]
            POS = str(int(POS)+1)
            CIGAR = CIGAR[POS_count+1:]
            break
        else:
            POS = "1"
            break
    POS_int = int(POS)

    SAM = str(READ_NAME + "\t" + "0" + "\t" + alignment.target.name + "\t" + POS + "\t" + "0" + "\t" + CIGAR + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + READ_SEQ + "\t" + "*" + "\n")
    return SAM

def alignToReference(fasta, read):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -1.0
    aligner.query_end_gap_score = 0.0
    allignments = aligner.align(fasta, read)
    alignment = allignments[0]
    return alignment

alignments = []
for read in reads:
    alignments.append(alignToReference(reference, read))

# print(alignments[10][0])
# print(alignments[10][1])
# print(makeCIGAR(alignments[10]))
# print(runLengthEncoding(makeCIGAR(alignments[10])))
# print(makeHeader(alignments[0]))

# print(runLengthEncodedCIGAR(alignments[0]))
# print(alignments[0])

# counter = 0
# for alignment in alignments[0:101]:
#     counter += 1
#     print(str(counter)+": "+runLengthEncodedCIGAR(alignment))

with open(sam_file, 'w') as f:
    f.write(makeHeader(alignments[0]))
    progress = 0
    for alignment in alignments:
        if progress % 100 == 0:
            print("Progress: "+ str(progress/len(alignments)))
        f.write(makeSAM(alignment))
        progress += 1
