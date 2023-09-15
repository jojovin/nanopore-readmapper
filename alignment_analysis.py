import SAM_creater
import re
import argparse

def main():
    argparser = argparse.ArgumentParser(description='Analysis of insertions in SAM file')
    argparser.add_argument('sam', help='sam file')
    argparser.add_argument('out', help='output file')
    args = argparser.parse_args()

    SAM, reflen = readSAM(args.sam)

    insmaps = []
    for line in SAM.SAMlines:
        insmaps.append(mapInsertionsFromSAMline(line))

    out_list=insertionDict(insmaps,reflen)
    with open(args.out, "w") as f:
        f.write("position,insertions\n")
        for i in range(len(out_list)):
            f.write(str(i)+","+str(out_list[i])+"\n")

def readSAM(sam_file):
    '''
    Reads a SAM file and returns a SAM object and the reference length
    '''
    SAM = SAM_creater.SAM()
    reflen = 0
    with open(sam_file, "r") as f:
        for line in f:
            if line[0:3] == "@SQ":
                reference = line.split("\t")
                for i in reference:
                    if i[0:2] == "LN":
                        reflen = int(i[3:])
                        break
                continue
            if line[0] == "@":
                continue
            else:
                S = SAM_creater.SAMline()
                S.fromSAM(line)
                SAM.append(S)
    return SAM, reflen

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

def insertionDict(maps,reflen=0):
    insertions = {}
    for map in maps:
        for i in map:
            if i not in insertions:
                insertions[i] = 1
            else:
                insertions[i] += 1
    for i in range(reflen):
        if i not in insertions:
            insertions[i] = 0
    positions=sorted(insertions)
    insertions_list = []
    for position in positions:
        insertions_list.append(insertions[position])         
    return insertions_list

if __name__ == "__main__":
    main()