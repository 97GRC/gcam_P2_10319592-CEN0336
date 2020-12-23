#!/usr/bin/env python3

import numpy as np

seq = "ATGCGATATGCGATCGATCGATCGCGCGCGCAGCTATAAACGTTGACGTTACTAGAC"
#Traduzir os genes
def Translate(seq):
	change1 = seq.replace('A', 'u')
	change2 = change1.replace('C', 'g')
	change3 = change2.replace('T', 'a')
	change4 = change3.replace('G', 'c')
	final = change4.upper()
	return final

seq_t = Translate(seq)

#Definir os frames
def Frames(seq):
	frame1 = [seq[i:i+3] for i in range(0, len(seq),3)]
	frame2 = [seq[i:i+3] for i in range(1, len(seq), 3)]
	frame3 = [seq[i:i+3] for i in range(2, len(seq), 3)]
	frame4 = [seq[i:i+3] for i in range(3, len(seq), 3)]
	frame5 = [seq[i:i+3] for i in range(4, len(seq), 3)]
	frame6 = [seq[i:i+3] for i in range(5, len(seq), 3)]
	Frames = [frame1, frame2, frame3, frame4, frame5, frame6]
	return Frames 

frames = Frames(seq_t)
#print(frames)

#Encontrar os ORFs

def StartPosition(frames):
	startPosition = []
	for i in frames:
	        if 'AUG' in i:
        	        startPosition.append(i.index('AUG'))
        	else:
                	startPosition.append('x')
	return startPosition

print(StartPosition(frames))


def StopPosition(frames):
	stopPosition = []
	for i in frames:
		if 'UAA' in i:
			stopPosition.append(i.index('UAA'))
		elif 'UGA' in i:
			stopPosition.append(i.index('UGA'))
		elif 'UAG' in i:
			stopPosition.append(i.index('UAG'))
		else:
			stopPosition.append('x')	
	return stopPosition

print(StopPosition(frames))


