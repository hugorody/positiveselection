#!/usr/bin/python3

import sys

msafile = sys.argv[1]
#output *.nogaps.fas

#parse msa FASTA file
seqs1 = {}
with open(msafile,"r") as set_fasta:
    for line in set_fasta:
        line = line.rstrip()
        if line != '':
            if line[0] == '>':
                words=line.split()
                name=words[0][1:]
                seqs1[name]=''
            else:
                seqs1[name] = seqs1[name] + line


gaps = [] #list of the positions having gaps "-"
for i in seqs1.items():
    id = i[0]
    seq = list(i[1])
    countgaps = 0
    for j in seq:
        countgaps += 1
        if j == "-": #gaps are computed for positions having "-"
            if countgaps not in gaps:
                gaps.append(countgaps)

fill = []
outfile1 = open(msafile+".nogaps.fas","w")
newseqs = {}
for i in seqs1.items():
    id = i[0]
    seq = list(i[1])
    countgaps = 0
    newseq = []
    for j in seq:
        countgaps += 1
        if countgaps not in gaps:
            newseq.append(j)
            if countgaps not in fill:
                fill.append(countgaps)
    outfile1.write(">"+id+"\n"+"".join(newseq)+"\n")
    newseqs[id] = "".join(newseq)


outcoordinates = open("alignment_coordinates.tsv","w")
outcoordinates.write("Position\tMaintained sites\tRemoved sites\n")

gaps_size = len(gaps)
fill_size = len(fill)

#define sizes
if gaps_size >= fill_size:
    myrange = gaps_size
else:
    myrange = fill_size


#create table
for i in range(0,myrange):
    position = i + 1

    if i <= fill_size - 1:
        letterfill = fill[i]
    else:
        letterfill = "NA"

    if i <= gaps_size - 1:
        lettergaps = gaps[i]
    else:
        lettergaps = "NA"

    outcoordinates.write (str(position)+"\t"+str(letterfill)+"\t"+str(lettergaps)+"\n")

#print stop codons if they are present in the alignment
stopcodons = ["TAA","TAG","TGA"]
for i in newseqs.items():
    seq = list(i[1])
    codonn = 0
    for j in range(0,len(seq),3):
        codonn += 1
        codon = "".join(seq[j:j+3])
        if codon in stopcodons:
            print (i[0],"".join(seq[j:j+3]),codonn*3)
