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

def runLengthEncodedCIGAR(alignment) -> str:
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
    def __init__(self, alignment, verbose=False) -> None:
        self.RNAME=alignment.target.name
        self.READ_SEQ=alignment.query.seq
        self.READ_NAME=alignment.query.name
        self.READ_ID=alignment.query.id
        self.SEQ_LEN=len(alignment.query.seq)
        self.FLAG="0"
        self.MAPQ="0"
        self.RNEXT="*"
        self.PNEXT="0"
        self.TLEN="0"
        self.QUAL="*"
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
        self.line = str(self.READ_NAME + "\t" + self.FLAG + "\t" + self.RNAME + "\t" + POS + "\t" + self.MAPQ + "\t" + self.CIGAR + "\t" + self.RNEXT + "\t" + self.PNEXT + "\t" + self.TLEN + "\t" + self.READ_SEQ + "\t" + self.QUAL + "\n")
        if verbose:
            print("Made SAM line for read " + self.READ_NAME)
        return
    
    def __str__(self) -> str:
        return self.line
    
    def __repr__(self) -> str:
        return self.line
    
    def __len__(self) -> int:
        return len(self.alignment)
    
class SAM:
    def __init__(self, SAMlines=[]) -> None:
        self.header = ""
        if len(SAMlines) != 0:
            #Make header from first SAM line
            self.header += makeHeader(SAMlines[0].RNAME)
        else:
            #No header available
            self.header = ""
        self.SAMlines = SAMlines
        return
    
    def __str__(self) -> str:
        return self.header + "".join(self.SAMlines)

    def append(self, SAMline: SAMline) -> None:
        if self.header == "":
            self.header += makeHeader(SAMline.RNAME)
        self.SAMlines.append(SAMline)
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