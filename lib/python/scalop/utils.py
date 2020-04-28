from .region_definitions import Accept
from .anarci import run_anarci

import os, sys, time, re
import pickle
import numpy as np

a = Accept()

resns = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H','I','L','K','M','F','P','S','T','W','Y','V']
cutoffs={'L1':0.82,'L2':1,'L3':0.91,'H1':0.8,'H2':0.63} # TODO: clustering cut-off

def getnumberedCDRloop(numberedAb,cdr,scheme,definition):
	a.numbering_scheme = scheme
	a.definition = definition
	a.set_regions(['CDR'+cdr])
	cdrloop = [pos for pos in numberedAb if a.accept(pos[0],cdr[-2]) and pos[1] != '-']
	if cdrloop == []: return [],[]
	# find 5 residues of anchors
	for pos in numberedAb: 
		if pos[1] == '-': numberedAb.remove(pos)
	startposi = sorted(numberedAb).index(sorted(cdrloop)[0])
	endposi = sorted(numberedAb).index(sorted(cdrloop)[len(cdrloop)-1])
	anchor = [x[0][0] for x in sorted(numberedAb)[startposi-5:startposi]+sorted(numberedAb)[endposi:endposi+5]]
	return sorted(cdrloop), anchor

def getnumberedCDRloop_builddb(_numberedAb,cdr,scheme,definition):
	a.numbering_scheme = scheme
	a.definition = definition
	a.set_regions([cdr])
	
	cdrloop = [pos for pos in _numberedAb if a.accept(pos[0],cdr[-2]) and pos[1] != '-']
	if cdrloop == []: return [],[]	
	numberedAb = [pos for pos in _numberedAb if pos[1] != '-']
	# find 5 residues of anchors
	startposi = sorted(numberedAb).index(sorted(cdrloop)[0])
	endposi = sorted(numberedAb).index(sorted(cdrloop)[len(cdrloop)-1])
	anchor = [sorted(numberedAb)[startposi-5],sorted(numberedAb)[endposi+5]]
	return sorted(cdrloop), anchor


def anarci_fun(list_of_inputs,numbering_scheme): # for parallel pool
	output = {chain:{Abentry[0]:{} for Abentry in list_of_inputs[chain]} for chain in list_of_inputs}
	for chain in list_of_inputs:
		anarcioutput = run_anarci(list_of_inputs[chain], scheme=numbering_scheme, assign_germline = False,ncpu=8)
		# matching Abid to the numbered Ab
		for eachAbi in range(len(anarcioutput[0])):
			output[chain][anarcioutput[0][eachAbi][0]] = dict(anarcioutput[1][eachAbi][0][0])
	return output

def get_PSSM(storepos):
	nseq={pos:sum(storepos[pos]) for pos in storepos}
	for pos in sorted(storepos):
		for resi in range(20):
			if storepos[pos][resi] == 0:
				storepos[pos][resi] = np.float(0.001)
			else:
				storepos[pos][resi] /=float(nseq[pos])
	pssm = {pos: np.log2(np.array(storepos[pos])/float(0.05)) for pos in sorted(storepos)}
	return pssm
	
# From BioPython
def SimpleFastaParser(handle):
	"""Generator function to iterate over Fasta records (as string tuples).
	
	For each record a tuple of two strings is returned, the FASTA title
	line (without the leading '>' character), and the sequence (with any
	whitespace removed). The title line is not divided up into an
	identifier (the first word) and comment or description.
	
	>>> with open("Fasta/dups.fasta") as handle:
		...     for values in SimpleFastaParser(handle):
		...         print(values)
		...
	('alpha', 'ACGTA')
	('beta', 'CGTC')
	('gamma', 'CCGCC')
	('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
	('delta', 'CGCGC')
	
	"""
	# Skip any text before the first record (e.g. blank lines, comments)
	while True:
		line = handle.readline()
		if line == "":
			return  # Premature end of file, or just empty?
		if line[0] == ">":
			break
	while True:
		if line[0] != ">":
			raise ValueError(
				"Records in Fasta files should start with '>' character")
		title = line[1:].rstrip()
		lines = []
		line = handle.readline()
		while True:
			if not line:
				break
			if line[0] == ">":
				break
			lines.append(line.rstrip())
			line = handle.readline()
			
		# Remove trailing whitespace, and any internal spaces
		# (and any embedded \r which are possible in mangled files
		# when not opened in universal read lines mode)
		yield title, "".join(lines).replace(" ", "").replace("\r", "")
		if not line:
			return  # StopIteration
	assert False, "Should not reach this line"
	
def write_txt(assignresults,outputdir, scheme, definition, graftedstructs=None, opf=""):
	lineout = ['Numbering scheme: {0} \nCDR definition: {1}'.format(scheme,definition)]
	lineout.append('\t'.join(['Input','CDR','Sequence','Canonical','Median']))
	if graftedstructs: lineout[-1]+='\tGrafted_Loop'
	if opf == "": opf = os.path.join(outputdir,'{0}.txt'.format(time.time()))
	for i in range(len(assignresults)):
		for cdr in sorted(assignresults[i]['outputs']):
			lineout.append('\t'.join([assignresults[i]['seqname']] + assignresults[i]['outputs'][cdr]))
			if graftedstructs and cdr in graftedstructs: lineout[-1]+='\t'+graftedstructs[cdr]
	if not os.path.exists(outputdir): os.mkdir(outputdir) 
	with open(opf,'w') as f:
		f.writelines('\n'.join(lineout))
	return opf

def print_result(assignresults, scheme, definition):
	lineout = ['Numbering scheme: {0} \nCDR definition: {1}'.format(scheme,definition)]
	lineout.append('\t'.join(['Input','CDR','Sequence','Canonical','Median']))
	
	for i in range(len(assignresults)):
		for cdr in sorted(assignresults[i]['outputs']):
			lineout.append('\t'.join([assignresults[i]['seqname']] + assignresults[i]['outputs'][cdr]))
	print(('\n'.join(lineout)))
