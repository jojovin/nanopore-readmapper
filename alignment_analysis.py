import SAM_creater
import re
import argparse

def main():
    argparser = argparse.ArgumentParser(description='Analysis of insertions in SAM file')
    argparser.add_argument('sam', help='sam file')
    argparser.add_argument('out', help='output file')
    args = argparser.parse_args()

    SAM = readSAM(args.sam)

    insmaps = []
    for line in SAM.SAMlines:
        insmaps.append(mapInsertionsFromSAMline(line))

    out_list=insertionDict(insmaps)
    with open(args.out, "w") as f:
        f.write("position,insertions\n")
        for i in range(len(out_list)):
            f.write(str(i)+","+str(out_list[i])+"\n")

def readSAM(sam_file):
    SAM = SAM_creater.SAM()
    with open(sam_file, "r") as f:
        for line in f:
            if line[0] == "@":
                continue
            else:
                S = SAM_creater.SAMline()
                S.fromSAM(line)
                SAM.append(S)
    return SAM

def mapInsertionsFromSAMline(SAMline: SAM_creater.SAMline):
    CIGAR = SAMline.CIGAR
    pos = SAMline.POS
    sCIGAR = re.split("(M|D|I)", CIGAR)
    map = []
    for i in range(len(sCIGAR)-1):
        if sCIGAR[i+1] == "I":
            if int(i)>2:
                map.append(pos)
        elif sCIGAR[i+1] == "D" or sCIGAR[i+1] == "M":
            pos+=int(sCIGAR[i])
        else:
            continue
    return map

def insertionDict(maps):
    insertions = {}
    for map in maps:
        for i in map:
            if i not in insertions:
                insertions[i] = 1
            else:
                insertions[i] += 1
    positions=sorted(insertions)
    insertions_list = []
    for position in positions:
        insertions_list.append(insertions[position])         
    return insertions_list