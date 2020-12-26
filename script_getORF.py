#!/usr/bin/env python3

import sys
import re

#QUESTÃO 2.1 -----------------------------------------------------
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

#QUESTÃO 2.2 ----------------------------------------------------

seq = [FastaDict[i] for i in FastaDict]
seqNames = list(FastaDict.keys())
#print(seqNames)
#print(seq)

def Transcription(seq):  
        change1 = seq.replace('A', 'u')
        change2 = change1.replace('C', 'g')
        change3 = change2.replace('T', 'a')
        change4 = change3.replace('G', 'c')
        final = change4.upper()
        return final

seq_t = [Transcription(i) for i in seq]
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

seqFrame2 = [val for sublist in Seq_Frame for val in sublist]
#print(seqFrame2)

seq_interesse = re.compile(r'AUG(?:(?!UAA|UAG|UGA)...)*(?:UAA|UAG|UGA)')

ORFs = []
value = 0
for i in seqFrame2:
	seq_int = seq_interesse.findall(seqFrame2[value])
	value += 1
	ORFs.append(seq_int)
#print(ORFs)

for i in ORFs:
	if not i:
		i.append('NONE')


cut1 = [i for i in range(0, len(ORFs), 6)]
cut2 = [i for i in range(6, len(ORFs)+1, 6)]
cuts = list(zip(cut1, cut2))
#print(cuts)
#print(cut2)

genesORFs = [ORFs[s:e] for s, e in cuts]
#print(genesORFs)

dictFrameORF = dict(zip(seqNames, genesORFs))
#print(dictFrameORF)


ORFgenes = []
for i in dictFrameORF:
	orfs = dictFrameORF[i]
	newlist = [val for sublist in orfs for val in sublist]
	ORFgenes.append(newlist)
#print(ORFgenes)	

dictORFs = dict(zip(seqNames, ORFgenes))
#print(dictORFs)


def maxLenght(ORFlist):
	lenght = []
	for i in range(0, len(ORFlist)):
		seq = ORFlist[i]
		size = len(seq)
		lenght.append(size)
		index = lenght.index(max(lenght))
	return ORFlist[index]
 
teste3 = (ORFgenes[0])
#print(len(maxLenght(teste3)))
#print(len(ORFgenes[1]))

maxORFs = []
for i in dictORFs:
	orfs = dictORFs[i]
	maxseq = maxLenght(orfs)
	maxORFs.append(maxseq)
#print(len(maxORFs[0]))

#QUESTÃO 2.3 ------------------------------------------------------------------

geneticCode = {"UUU":"F", "UUC": "F", "UUA":"L", "UUG":"L",
             "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
             "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
             "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
       	     "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
             "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
             "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
             "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
             "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
             "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
             "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
             "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
             "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
             "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
             "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
             "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


maxORFs_ = []
for i in maxORFs:
	seq = i.replace(' ','')
	maxORFs_.append(seq)

remover = 'NONE'
while remover in maxORFs_:
	maxORFs_.remove(remover)
print(maxORFs_)
#print(maxORFs_[102])

# ------------------------Verificando se é miltiplo de 3
a = [len(i) for i in maxORFs_]
#print(a)

for i in a:
	if not(i%3):
		print('S')
	else:
		print(a.index(i))

def Translate(ORF, genetic_code):
	protein = ''
	for i in range(0, len(ORF), 3):
		codon = ORF[i:i+3]
		protein += genetic_code[codon]
	return protein

Proteins = []
for i in maxORFs_:
	protein = Translate(i, geneticCode)
	Proteins.append(protein)
print(Proteins)
	

#QUESTÃO 2.4 ------------------------------------------------------------------



