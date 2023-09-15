from Bio import Align
from Bio import SeqIO
import sys
import argparse
import SAM_creater

argparser = argparse.ArgumentParser(description='Read mapper for RNA nanopore full-length reads')
argparser.add_argument('fasta', help='fasta file with reference sequence')
argparser.add_argument('fastq', help='fastq file with reads')
argparser.add_argument('sam', help='sam file output')
argparser.add_argument('-m','--match', required=False, default=1.0, help='match score (0 or positive)', type=float)
argparser.add_argument('-x','--mismatch', required=False, default=0.0, help='mismatch score (0 or negative)', type=float)
argparser.add_argument('-o','--gapopen', required=False, default=2.0, help='gap open score (Should be positive, as it will be treated as negative)', type=float)
argparser.add_argument('-e','--gapext', required=False, default=1.0, help='gap extension score (Should be positive, as it will be treated as negative)', type=float)
args = argparser.parse_args()

fasta = args.fasta
fastq = args.fastq
sam_file = args.sam

reference = SeqIO.read(fasta, "fasta")
reads = list(SeqIO.parse(fastq, "fastq"))

def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.3+
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print("{}[{}{}] {}/{}".format(prefix, "#"*x, "."*(size-x), j, count), 
                end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)

print("Aligning reads to reference...")
alignments = SAM_creater.alignToReference(reference, reads, match=args.match, mismatch=args.mismatch, gapopen=-args.gapopen, gapext=-args.gapext,verbose=False)
print("Finished aligning!")

SAM = SAM_creater.SAM()
SAM.header = SAM_creater.makeHeader(reference.name, str(len(reference.seq)))

for alignment in progressbar(range(len(alignments)), "Writing SAM file: ", 40):
    S = SAM_creater.SAMline()
    S.fromAlignment(alignments[alignment])
    SAM.append(S)

SAM.WriteSAM(sam_file)