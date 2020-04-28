#!/usr/bin/env python
"""
graft_loop.py
Description:    Follows on from FREAD but is a simpler, loop grafting method
Jun 23, 2016
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys, os, json, re, datetime, pickle, math
import numpy as np
from array import array
helpdoc = """
LoopGraft
Requirements:
    FREAD -- (c) Sebastian Kelm, Charlotte M. Deane

"""
epilogue = """
(c) 2016,   Authors:    Jinwoo Leem
            Supervisor: Charlotte M. Deane
            Group:      Oxford Protein Informatics Group, University of Oxford
"""

from scalop.prosci.util.pdb import Pdb, residueCode
from scalop.prosci.util.residue import ResidueList, Residue
from scalop.prosci.loops.fread_meld import meld, add_oxygens
from scalop.prosci.loops import ccd
from scalop.prosci.esst import load_esst


### Define functions ###

def rename_chain(residuelist, target_chain):
	for r in residuelist:
		for a in r:
			a.chain = target_chain

def cut_side_chains(residuelist):
	"""
		Cut out the side chains and just keep the backbone and CB
	"""
	tmp = ResidueList([])
	for r in residuelist:
		l = list(r.iter_backbone())
		l.append(r.CB)
		tmp.append(Residue(l))
	return tmp
		
# mini-FREAD
def miniFREAD(inputseq,esst,cmlist,freaddbcdr,blacklist=[]):
	tables = esst.tables
	ascii2index = esst.ascii2index
	seqmap = tuple([ascii2index[ord(s)] for s in inputseq])

	prevdihed=["", [None]*len(seqmap)]
	def score_sequence(loopseq, dihed):
		if prevdihed[0] != dihed:
			prevdihed[0] = dihed
			for i,x in enumerate(dihed):
				prevdihed[1][i] = tables[int(x)][seqmap[i]]
		score=0
		for i,x in enumerate(loopseq):
			score += prevdihed[1][i][ascii2index[ord(x)]]  # Speed-optimised version
		return score
	essstore = []
	for cm,cmseq,cmdih,a1,a2,a3,a4 in freaddbcdr:
		if cm not in cmlist or cm in blacklist or len(cmseq)-4!=len(inputseq): continue
		essstore.append((score_sequence(cmseq[2:-2],cmdih[2:-2]),cm,l2obj([np.asarray(a1),np.asarray(a2),np.asarray(a3),np.asarray(a4)])))
	essstore.sort(reverse=True)
	return [cminfo for cminfo in essstore if cminfo[0] == essstore[0][0]] # cluster members with the maximum ess

# Decoding
class emptyatom(object):
	def __init__(self):
		self.xyz = np.zeros(3)
class emptyrl(object):
	def __init__(self):
		self.C = self.CA = self.N = self.O = emptyatom()
def l2obj(anchorres):
	anchorsobj = []
	for a in anchorres:
		aobj = emptyrl()
		aobj.C.xyz,aobj.CA.xyz,aobj.N.xyz,aobj.O.xyz = a         
		anchorsobj.append(aobj)
	return anchorsobj



### Parsing step ###
# Parse arguments
ess = load_esst()
def export_structure(args, assignmentresults, output=''):
	anchor_length = 2
	scalop_path  = os.path.split(__file__)[0] # from scalop database
	# Parse files
	is_model_generated = 0
	p_istructure = Pdb(args.structuref)
	if (args.hc!='' and args.hc not in p_istructure.get_chain_codes()) or (args.lc!='' and args.lc not in p_istructure.get_chain_codes()): 
		print("Input chain IDs do not match")
		return None
	tlist = []
	if args.dbv == 'latest':
		for n in os.listdir(os.path.join(scalop_path,'database')):
			if re.match(r'{0}_{1}_v(\d+)-(\d+).pickle'.format(args.scheme,args.definition),n) == None: continue
			nt = datetime.datetime.strptime(''.join(re.match(r'{0}_{1}_v(\d+)-(\d+).pickle'.format(args.scheme,args.definition),n).groups()),'%Y%m')
			tlist.append((nt,n))
		fname = tlist[tlist.index(max(tlist))][1]
	else:
		fname = '{0}_{1}_v{2}.pickle'.format(args.scheme,args.definition,args.dbv)	

	with open(os.path.join(scalop_path,'database',fname)) as f:
		_,_,cm = pickle.load(f) 
	with open(os.path.join(scalop_path,'database','{0}_{1}.freaddb'.format(args.scheme,args.definition))) as f:
		freaddb = json.load(f) 
	
	loopstograft = {}
	graftedstructs = {}
	for scalopres in assignmentresults:
		if 'outputs' not in scalopres: continue
		c = None
		if 'H' in list(scalopres['outputs'].keys())[0]:
			c = args.hc
		elif 'L' in list(scalopres['outputs'].keys())[0]:
			c = args.lc
		if c=='' or c == None: 
			print('Missing input chain ID')
			continue
		iscompletestruc = scalopres['input'][1]==p_istructure.get_chain(c).get_seq()
		localpdb = p_istructure.get_chain(c)
		if not iscompletestruc:
			cdrlens = [len(n[1]) for n in sorted(scalopres['outputs'].values())]
			if len(localpdb.get_seq()) + sum(cdrlens) != len(scalopres['input'][1]): 
				print('The input structure ({}) and CDR sequences ({}) do not match with input sequence ({})'.format(len(localpdb.get_seq()) , sum(cdrlens) , len(scalopres['input'][1])))
				continue
			cdri = 0
		for cdr,loop,cl,_ in sorted(scalopres['outputs'].values()):
                        if cl == 'None': 
                                if not iscompletestruc: cdri += 1
                                continue
                        if iscompletestruc:	
                                lsi,lei = scalopres['input'][1].index(loop), scalopres['input'][1].index(loop) + len(loop)
                        else:
                                lsi=lei = scalopres['input'][1].index(loop)-sum(cdrlens[:cdri])
                                cdri += 1
					
			
                        ri = ResidueList(localpdb).deep_copy()
                        i_N = ResidueList(localpdb).deep_copy()[lsi-2:lsi]
                        i_C = ResidueList(localpdb).deep_copy()[lei:lei+2]
                        if hasattr(args,'blacklist'):
                                maxessanchor = miniFREAD(loop,ess,cm[cdr][cl]['Cluster members'],freaddb[cdr],args.blacklist)
                        else:
                                maxessanchor = miniFREAD(loop,ess,cm[cdr][cl]['Cluster members'],freaddb[cdr])
                        rmsdcm = []
                        for essscore,loopname,anchors in maxessanchor:
                                if essscore < 25:continue
                                rmsdcm.append((ccd.superimpose(anchors, i_N+i_C),loopname, anchors))
                        if rmsdcm == []:
                                print('No loop with ESS > 25 fits')
                                continue
                        loopname,anchors = sorted(rmsdcm)[0][1:]
                        graftedstructs[cdr]=str(loopname)
                        p_loop = Pdb(str(os.path.join(args.loopdb,'{0}_{1}'.format(args.scheme,args.definition),cdr,loopname+'.pdb')))
                        rloop = ResidueList(p_loop).deep_copy()[3:-3]
                        ccd.superimpose_withloop(anchors,i_N+i_C,rloop)
                        i_Nn = i_N.deep_copy()
                        i_Cn = i_C.deep_copy()
                        meld(i_Nn, rloop[:2])
                        meld(i_Cn, rloop[-2:], invertWeights = True)
                        rloop = rloop[2:-2]
                        rloop[-1].O = None # This oxygen is wrong anyway, so delete it
                        add_oxygens(i_Nn+rloop[:1], force=True) # Add oxygen to the N-anchor
                        add_oxygens(rloop[-1:]+i_Cn+[ri[lei+3]], force=True) # Add oxygen to C-anchor, and add the i+1th residue to the C-terminus where i is anchor length
                        rename_chain(rloop, c)

                        for loopi in range(len(rloop)):
                                rloop[loopi].set_type(residueCode(loop[loopi]))
                        if c not in loopstograft: loopstograft.update({c:[]})
                        loopstograft[c].append((lsi,cut_side_chains(i_Nn.deep_copy()),cut_side_chains(rloop.deep_copy()),cut_side_chains(i_Cn.deep_copy()),lei,ri.deep_copy()))
		if output != '':
			model = ResidueList([])	
			for c in loopstograft:
				localpdb = p_istructure.get_chain(c)
				starti = 0
				numcdr = len(loopstograft[c])
				cloopi = 0
				for lsi,i_Nn,rloop,i_Cn,lei,ri in sorted(loopstograft[c]):
					cloopi += 1
					### ASSEMBLE MODEL ###
				 	
					# Copy everything preceding the loop
					for i in range(starti,lsi-2):
						model.append(ri[i].copy())
				 
					# Add the newly formed anchor and loop
					model.extend(i_Nn)
					model.extend(rloop)
					model.extend(i_Cn)
				 
					# Add C-terminal end of structure
					end_of_segment = len(ri) if cloopi == numcdr else sorted(loopstograft[c])[cloopi][0]-2
					
					for i in range(lei+2, end_of_segment):
						model.append(ri[i].copy())
				 
					# Replace main chain oxygen before the loop, to ensure the peptide bond there is planar
					model[lsi-1].O = None
					add_oxygens(model, start=lsi-1, end=lsi+1, force=True)
					if cloopi != numcdr: starti = sorted(loopstograft[c])[cloopi][0]-2
			model.renumber_atoms()
			with open(output, 'w') as modelfile:
				print(model, file=modelfile)
			is_model_generated = 1
	return graftedstructs, is_model_generated
