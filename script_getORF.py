#!/usr/bin/env python3

import sys
import re

#QUESTÃO 2.1 -----------------------------------------------------
#Abrindo o arquivo multifasta a partir da linha de comando usando sys e salvando em um dicionário
File = sys.argv[1]
FastaDict = {}
SeqName = ''
Seq = ''

with open(File,'r') as fasta_file:
        for line in fasta_file:
                line = line.rstrip()
                #print(line)
                if(line.startswith('>')):
                        SeqName = line.split('\n')[0]
                        SeqName = SeqName.replace('>','')
                        Seq = ''
                else:
                        Seq += line
                        FastaDict[SeqName] = Seq
#print(FastaDict)

#QUESTÃO 2.2 ----------------------------------------------------

seq = [FastaDict[i] for i in FastaDict] #Lista com as sequências de DNA
seqNames = list(FastaDict.keys()) #Lista com os nomes das sequências de DNA
#print(seq)
#print(seqNames)

#Função para trascrever a sequêmcia de DNA em mRNA, que posteriormente será traduzida
def Transcription(seq):  
        change1 = seq.replace('A', 'u')
        change2 = change1.replace('C', 'g')
        change3 = change2.replace('T', 'a')
        change4 = change3.replace('G', 'c')
        final = change4.upper()
        return final

seq_t = [Transcription(i) for i in seq] #Lista com as sequências transcritas
#print(seq_t)

#Função para definir os seis diferentes frames
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
	return Frames

#Loop para obter os frames de cada sequência -> Gera uma lista onde cada valor é composto por 6 listas (uma para cada frame)
value = 0
Seq_Frame = []
for i in seq_t:
        frame = Frames(seq_t[value])
        value += 1
        Seq_Frame.append(frame)
#print(Seq_Frame)

#Transforma a lista de lista em apenas uma lista, onde cada 6 valores corresponde a um fasta
seqFrame2 = [val for sublist in Seq_Frame for val in sublist]
#print(seqFrame2)

#Objeto para armazenar o código de busca dos ORFs
seq_interesse = re.compile(r'AUG(?:(?!UAA|UAG|UGA)...)*(?:UAA|UAG|UGA)')

#Loop para procurar dentro de cada frames todos os ORFs presentes
ORFs = []
value = 0
for i in seqFrame2:
	seq_int = seq_interesse.findall(seqFrame2[value])
	value += 1
	ORFs.append(seq_int)
#print(ORFs)

#Caso algum frame nãocodifique nenhum ORF o valor 'NONE' é inserido
for i in ORFs:
	if not i:
		i.append('NONE')

#Sequências para colocar selecionar os frames de cada ORF e armazena-los em uma lista 
cut1 = [i for i in range(0, len(ORFs), 6)] #Valores iniciais do corte
cut2 = [i for i in range(6, len(ORFs)+1, 6)] #Valores finais do corte
cuts = list(zip(cut1, cut2)) 
#print(cuts)

genesORFs = [ORFs[s:e] for s, e in cuts] 
#Os valores da primeira lista representam cada ID e dentro de cada valor há uma lista com os ORFs de cada frame 
#print(genesORFs) 

#Dicionário onde o nome da ID são é a key e os values são os ORFs armazenados em listas   
dictFrameORF = dict(zip(seqNames, genesORFs))
#print(dictFrameORF)

#Loop para criar uma lista onde cada valor é composto pelos ORFs de cada ID
ORFgenes = []
for i in dictFrameORF:
	orfs = dictFrameORF[i]
	newlist = [val for sublist in orfs for val in sublist]
	ORFgenes.append(newlist)
#print(ORFgenes)	

#Dicionário onde a key é a ID e o value os ORFs respectivos
dictORFs = dict(zip(seqNames, ORFgenes))
#print(dictORFs)

#Função para encontrar o ORF mais comprido dentre os seis frames
def maxORF(ORFlist):
	lenght = []
	for i in range(0, len(ORFlist)):
		seq = ORFlist[i]
		size = len(seq)
		lenght.append(size)
		index = lenght.index(max(lenght))
	return ORFlist[index]
 
#Lista com o maior ORF de cada ID
maxORFs = []
for i in dictORFs:
	orfs = dictORFs[i]
	maxseq = maxORF(orfs)
	maxORFs.append(maxseq)
#print(maxORFs)

#Dicionário com os maiores ORFs e seus ID
dictmaxORFs = dict(zip(seqNames, maxORFs))
for k, v in dictmaxORFs.items():
	print('O ORF máximo para a sequência: ' + str(k) + ' é:' + '\n' + str(v))
print(' ')

#QUESTÃO 2.3 ------------------------------------------------------------------

#Código de tradução:
Legenda = 'F = Fenilalanina; L = Leucina; S = Serina; Y = Tirosina; C = Cisteína; W = Triptofano; P = Prolina; H = Histidina; Q = Glutamina; R = Arginina; I = Isoleucina; M = Metionina; T = Treonina; N = Aspargina; K = Lisina; V = Valina; A = Alanina; D = Ác. Aspártico; E = Ác. Glutâmico; G = Glicina'

geneticCode = {"UUU":"F", "UUC": "F", "UUA":"L", "UUG":"L",
             "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
             "UAU":"Y", "UAC":"Y", "UAA":"stop", "UAG":"stop",
             "UGU":"C", "UGC":"C", "UGA":"stop", "UGG":"W",
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

#Rretirando o ' ' entre os códons
maxORFs_ = []
for i in maxORFs:
	seq = i.replace(' ','')
	maxORFs_.append(seq)

#Função para traduzir o código em um polipeptídeo 
def Translate(ORF, genetic_code):
	protein = ''
	for i in range(0, len(ORF), 3):
		codon = ORF[i:i+3]
		protein += genetic_code[codon]
	return protein

#Lista com a proteína referente a cada ORF máximo
Proteins = []
for i in maxORFs_:
	protein = Translate(i, geneticCode)
	Proteins.append(protein)
#print(Proteins)

#Dicionário onde as proteínas são os value e as ID as keys
dictProteins = dict(zip(seqNames, Proteins))
for k, v in dictProteins.items():
	print('O polipeptídeo codificado pelo max ORF de: ' + str(k) + ' é: ' + '\n' + str(v))
print(Legenda)

#QUESTÃO 2.4 ------------

#Criando o arquivo ORF.nfa:
fileORFs = 'ORF.fna'
File = open(fileORFs, 'w')

for k, v in dictmaxORFs.items():
	File.write('>' + str(k) + '\n' + str(v) + '\n\n')

File.close()

#QUESTÃO 2.5 ----------------

#Criando o arquivo ORF.faa:
dictProteins = dict(zip(seqNames, Proteins))
#print(dictProteins)

fileProtein = 'ORF.faa'
FileP = open(fileProtein, 'w')

for k, v in dictProteins.items():
	FileP.write('>' + str(k) + '\n' + str(v) + '\n\n')

FileP.close()
