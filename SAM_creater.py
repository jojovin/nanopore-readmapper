from Bio import Align



def makeCIGAR(alignment) -> str:
    cigar = ""
    for i in range(len(alignment[0])):
        if alignment[0][i] == "-":
            cigar += "I"
        elif alignment[1][i] == "-":
            cigar += "D"
        else:
            cigar += "M"
    return cigar

def runLengthEncodedCIGAR(alignment) -> str:
    def runLengthEncoding(cigar) -> str:
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
    cigar = makeCIGAR(alignment)
    cigar_rle = runLengthEncoding(cigar)
    return cigar_rle

def makeHeader(reference: str, length="0") -> str:
    header = ""
    header += "@HD" + "\t" + "VN:1.6" + "\t" + "SO:unsorted" + "\n"
    header += "@SQ" + "\t" + "SN:" + reference + "\t" + "LN:" + length + "\n"
    return header

def alignToReference(fasta, reads, match=1.0, mismatch=0.0, gapopen=-2.0, gapext=-1.0,verbose=False):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gapopen
    aligner.extend_gap_score = gapext
    aligner.query_end_gap_score = 0.0
    alignments = []
    if verbose:
        print("Aligning reads to reference...")
    count = 0
    for read in reads:
        if verbose:
            print("Aligning read " + str(count) + " of " + str(len(reads)-1) + "...")
        alignment = aligner.align(fasta, read)
        alignments.append(alignment[0])
        count += 1
    if verbose:
        print("Finished aligning!")
    return alignments

class SAMline:
    def __init__(self, QNAME="*", FLAG=0, RNAME="*", POS=0, MAPQ=0, CIGAR="*", RNEXT="*", PNEXT=0, TLEN=0, SEQ="*", QUAL="*") -> None:
        self.READ_NAME=QNAME
        self.FLAG=FLAG
        self.REF_NAME=RNAME
        self.POS=POS
        self.MAPQ=MAPQ
        self.CIGAR=CIGAR
        self.RNEXT=RNEXT
        self.PNEXT=PNEXT
        self.TLEN=TLEN
        self.READ_SEQ=SEQ
        self.QUAL=QUAL
        self.line=str(self.READ_NAME + "\t" + str(self.FLAG) + "\t" + self.REF_NAME + "\t" + str(POS) + "\t" + str(self.MAPQ) + "\t" + self.CIGAR + "\t" + self.RNEXT + "\t" + str(self.PNEXT) + "\t" + str(self.TLEN) + "\t" + self.READ_SEQ + "\t" + self.QUAL + "\n")
        return

    def fromAlignment(self, alignment, verbose=False) -> None:
        self.REF_NAME=alignment.target.name
        self.READ_SEQ=alignment.query.seq
        self.READ_NAME=alignment.query.name
        self.READ_ID=alignment.query.id
        self.SEQ_LEN=len(alignment.query.seq)
        POS_count = 0
        self.CIGAR=runLengthEncodedCIGAR(alignment)
        for c in self.CIGAR:
            if c.isdigit():
                POS_count += 1
            elif c == "D":
                POS = self.CIGAR[0:POS_count]
                POS = str(int(POS)+1)
                self.CIGAR = self.CIGAR[POS_count+1:]
                break
            else:
                POS = "1"
                break
        POS_int = int(POS)
        self.line=str(self.READ_NAME + "\t" + str(self.FLAG) + "\t" + self.REF_NAME + "\t" + str(POS) + "\t" + str(self.MAPQ) + "\t" + self.CIGAR + "\t" + self.RNEXT + "\t" + str(self.PNEXT) + "\t" + str(self.TLEN) + "\t" + self.READ_SEQ + "\t" + self.QUAL + "\n")
        if verbose:
            print("Made SAM line for read " + self.READ_NAME + "from alignment.")
        return
    
    def fromSAM(self, line) -> None:
        self.line = line
        line = line.split("\t")
        self.READ_NAME=line[0]
        self.FLAG=int(line[1])
        self.REF_NAME=line[2]
        self.POS=int(line[3])
        self.MAPQ=int(line[4])
        self.CIGAR=line[5]
        self.RNEXT=line[6]
        self.PNEXT=int(line[7])
        self.TLEN=int(line[8])
        self.READ_SEQ=line[9]
        self.QUAL=line[10]
        return
    
    def __str__(self) -> str:
        return self.line
    
    def __repr__(self) -> str:
        return self.line
    
    def __len__(self) -> int:
        return len(self.READ_SEQ)
    
class SAM:
    def __init__(self, SAMlines=[]) -> None:
        self.header = ""
        if len(SAMlines) != 0:
            #Make header from first SAM line
            self.header += makeHeader(SAMlines[0].REF_NAME)
        else:
            #No header available
            self.header = ""
        self.SAMlines = SAMlines
        return
    
    def makeHeader(self, reference: str, length="0") -> str:
        self.header = ""
        self.header += "@HD" + "\t" + "VN:1.6" + "\t" + "SO:unsorted" + "\n"
        self.header += "@SQ" + "\t" + "SN:" + reference + "\t" + "LN:" + length + "\n"
        return self.header
    
    def __str__(self) -> str:
        return self.header + "".join(str(self.SAMlines))

    def append(self, line: SAMline) -> None:
        if self.header == "":
            self.header += makeHeader(line.REF_NAME)
        self.SAMlines.append(line)
        return

    def __len__(self) -> int:
        return len(self.SAMlines)

    def WriteSAM(self, sam_file: str) -> None:
        if self.header == "":
            print("No header available. SAM file will NOT be valid! (Continuing anyway...)")
        with open(sam_file, 'w') as f:
            f.write(self.header)
            for line in self.SAMlines:
                f.write(line.line)
        return