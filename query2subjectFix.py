#!/usr/bin/env python3
import csv
from Bio import SeqIO
import os
import re

def createFile(outputFile, inputFile):

	try:
			
		if os.path.exists(outputFile):
			os.remove(outputFile)
		## Python will convert \n to os.linesep

		if not os.path.exists(outputFile):
		    with open(outputFile, 'w'): pass

		f = open(outputFile,'w')

		x = 0;
		for record in SeqIO.parse(inputFile, "fasta"):
			
			# print(record)
			desc = record.description
			if "[" not in record.description:
				continue;

			genome = desc[desc.find('[')+1: desc.find(']') ]
			genome = genome.replace(' ', '_');
			line = record.id + ',' + genome.split(".")[0];

			f.write(line + '\n') #Give your csv text here.
		pass

	except Exception as e:

		print( str(e) + ' Script ended')
		exit();

	return True;


