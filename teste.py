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

def Translate(seq):  #MUDAR O NOME
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
#print(frameName)
dictFrame = {}
#------------------------------------

teste = Seq_Frame[1][1]
#print(teste)

def maxORF(teste):
	ORFmax = re.findall(r'AUG(?:(?!UAA|UAG|UGA)...)*(?:UAA|UAG|UGA)',teste)
	return ORFmax

seqFrame2 = [val for sublist in Seq_Frame for val in sublist]
#print(seqFrame2)
#print(maxORF(teste))

seq_interesse = re.compile(r'AUG(?:(?!UAA|UAG|UGA)...)*(?:UAA|UAG|UGA)')

ORFs = []
value = 0
for i in seqFrame2:
	seq_int = seq_interesse.findall(seqFrame2[value])
	value += 1
	ORFs.append(seq_int)
#print(ORFs)

ORFsGene1 = ORFs[0:6]  #TENTAR TRANSFORMAR NUM LOOP
ORFsGene1_ = [valores for sublista in ORFsGene1 for valores in sublista]
ORFsGene2 = ORFs[6:12]
ORFsGene2_ = [valores for sublist in ORFsGene2 for valores in sublist]

value = 0
lenght1 = []
for i in ORFsGene1_:
	seq = (i.replace(' ','')) 
	lenght = (len(seq))
	lenght1.append(lenght)

value = 0
lenght2 = []
for i in ORFsGene2_:
	seq = (i.replace(' ', ''))
	lenght = (len(seq))
	lenght2.append(lenght)

lenghtMaxORF = []
lenghtMaxORF.append(max(lenght1))
lenghtMaxORF.append(max(lenght2))
#print(type(lenghtMaxORF[0]))

index1 = lenght1.index(lenghtMaxORF[0])
index2 = lenght2.index(lenghtMaxORF[1])

maxORFs = []
maxORFs.append(ORFsGene1_[index1])
maxORFs.append(ORFsGene2_[index2])
#print(maxORFs)

#QUESTÃO 2.3

Peptideos = {"UUU":"F", "UUC": "F", "UUA":"L", "UUG":"L",
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

def Translate2(ORF, genetic_code):
	protein = ''
	for i in range(0, len(ORF), 3):
		codon = ORF[i:i+3]
		protein += genetic_code[codon]
	return protein

Sera = maxORFs_[0]
print(Translate2(Sera, Peptideos))
	





