#!/usr/bin/env python

import sys
import os

inFile = open(sys.argv[1],'r')
FileName = os.path.basename(sys.argv[1])
FileName1 = FileName.replace('.mpileup','')

if os.path.getsize(sys.argv[1])>0:
	
	for line in inFile:

		############################
		#Initiate types for output
		types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

		############################
		# Read file if not empty
		if len(line) > 0:

			############################
			#Check healthy data 
			data = line.strip().split('\t')
			doctor = data[0]
			if "chr" in doctor and (len(data)>=5):  #len(data)>5 because some .mpileup files are corrupted and only contain 4 fields
	
				############################
				#Obtain all info from line
				chromosome= data[0]
				bp = data[1]
				bases = data[4].upper()

				############################
				#Parse information
				i = 0
				jump = 0
				while i < len(bases):
					base = bases[i]
					if base == '^' or base == '$':
						i += 1
					elif base == '-':
						i += 1
						types['-'] += 1
						if bases[i+1].isdigit():    # Addition for deletions longer than 9 bp
							jump = int(bases[i] + bases[i+1])
							i += 1
						else:
							jump = int(bases[i])
						i = i + jump
						i += 1   # Addition to fix MPileup bug - an additional character is inserted on every addition/deletion 
					elif base == '*':
						i += 1
					elif base == '+':
						i += 1
						if bases[i+1].isdigit():    # Addition for insertions longer than 9 bp
							addNum = int(bases[i] + bases[i+1])
							i += 1
						else:
							addNum = int(bases[i])	
						addSeq = ''
						for a in range(addNum):
							i += 1
							addSeq += bases[i]
						i += 1  # Addition to fix MPileup bug - an additional character is inserted on every addition/deletion 
						types['+'].append(addSeq)
					else:
						if types.has_key(base):
							types[base] += 1
						else:
							types['X'].append(base)

					i += 1

				adds = '.'
				if len(types['+']) > 0:
						adds = ','.join(types['+'])

				amb = '.'
				if len(types['X']) > 0:
						amb = ','.join(types['X'])
		
			############################
			#Write error on file if corrupted data
			else:
				chromosome = "CD"
				bp = "CD"
				adds = "."
				amb = "."
		
		else:
			chromosome = "ND"
			bp = "ND"
			adds = "."
			amb = "."
		############################
		# End of file analysis

		############################
		# Print output data - even if file is empty

		out = [FileName1,chromosome,bp,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb]
		print '\t'.join([str(x) for x in out])
else:
	############################
	# Write line of empty output
	types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}	
	chromosome = "ND"
	bp = "ND"
	adds = "."
	amb = "."
	out = [FileName1,chromosome,bp,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb]
	print '\t'.join([str(x) for x in out])
