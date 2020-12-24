#!/usr/bin/env python3

import sys
import re

File = sys.argv[1]
FastaDict = {}
SeqName = ''
Seq = ''

with open(File,'r') as fasta_file:
        for line in fasta_file:
                line = line.rstrip()
                #print(line)
                if(line.startswith('>')):
                        SeqName = line.split(' ')[0]
                        SeqName = SeqName.replace('>','')
                        Seq = ''
                else:
                        Seq += line
                        FastaDict[SeqName] = Seq

#print(FastaDict)
#print('')

seq = [FastaDict[i] for i in FastaDict]
seqNames = list(FastaDict.keys())
#print(seqNames)
#print(seq)

def Translate(seq):
        change1 = seq.replace('A', 'u')
        change2 = change1.replace('C', 'g')
        change3 = change2.replace('T', 'a')
        change4 = change3.replace('G', 'c')
        final = change4.upper()
        return final

seq_t = [Translate(i) for i in seq]
#print(seq_t)

#Definir os frames
def Frames(seq):
	codon1 = [seq[i:i+3] for i in range(0, len(seq),3)]
	Frame1 = ' '.join(codon1)
	codon2 = [seq[i:i+3] for i in range(1, len(seq), 3)]
	Frame2 = ' '.join(codon2)
	codon3 = [seq[i:i+3] for i in range(2, len(seq), 3)]
	Frame3 = ' '.join(codon3)
	seq_inv = seq[::-1]
	codon4 = [seq_inv[i:i+3] for i in range(0, len(seq_inv), 3)]
	Frame4 = ' '.join(codon4)
	codon5 = [seq_inv[i:i+3] for i in range(1, len(seq_inv), 3)]
	Frame5 = ' '.join(codon5)
	codon6 = [seq_inv[i:i+3] for i in range(2, len(seq_inv), 3)]
	Frame6 = ' '.join(codon6)
	Frames = [Frame1, Frame2, Frame3, Frame4, Frame5, Frame6]
	Codons = [codon1, codon2, codon3, codon4, codon5, codon6]
	return Frames 

value = 0
Seq_Frame = []
for i in seq_t:
        frame = Frames(seq_t[value])
        value += 1
        Seq_Frame.append(frame)
#print(Seq_Frame)


#-------------------------------------------------------------
lenght = len(Seq_Frame[0])
frameName = []
frameName1 = []
frameName2 = []

for i in range(1, lenght+1):
	names = f'Frame_{i}'
	frameName1.append(names)
	frameName2.append(names)
frameName.append(frameName1)
frameName.append(frameName2)

dictFrame = {}
#------------------------------------

teste = Seq_Frame[1][1]
#print(teste)

def maxORF(teste):
	ORFmax = max(re.findall(r'AUG(?:(?!UAA|UAG|UGA)...)*(?:UAA|UAG|UGA)',teste), key = len)
	return ORFmax

#print(maxORF(teste))

ORFs = []
for i in Seq_Frame:
	for a in i:
		if 'AUG' in a:
			ORF = maxORF(a)
			ORFs.append(ORF)
		else:
			result = 'NONE'
			ORFs.append(results)
print(ORFs)
