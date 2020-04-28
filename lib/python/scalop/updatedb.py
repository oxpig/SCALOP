import os, sys, re, operator,json, argparse, subprocess, time, datetime
import pickle as pickle
import numpy as np
import pandas as pd
import multiprocessing
import multiprocessing.pool
from string import ascii_uppercase
import subprocess

from .anarci import run_anarci
from .utils import *


a = Accept()
cdrs = ["H1","H2","L1","L2","L3"]
scalop_path  = os.path.split(__file__)[0]
thism = datetime.date.today().replace(day=1)
lastm = thism - datetime.timedelta(days=1)
res_dict = {
		 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
		 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
		 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
		 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'TYS': 'Y', 
		 'MET': 'M', 'MSE': 'M'
}
def loop_parse(s, lanchor,lcdr, ucdr,uanchor, cdr, Abchain,bfactor):

	residuelist = []

	for r in s.get_residues():
		if hasattr(r, 'chain_type') and r.chain_type == cdr[-2] and r.id[1] >= lanchor and r.id[1] <= uanchor and r.get_full_id()[3].lower() == Abchain.lower() and r.get_resname() in res_dict:
			for atom in r:
				if r.id[1]>=lcdr and r.id[1]< ucdr and atom.name in ["CA","C","N","O"] and atom.bfactor >= bfactor: # only those in the CDR loop backbone
					return []
			residuelist.append(r)
	return residuelist

class NoDaemonProcess(multiprocessing.Process):
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)
class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def runbuilddb(args):
	global allVanarcioutput, Abs
	dbdir = args.dbdir
	print('dbdir', dbdir)
	struc_cutoff = args.struc_cutoff
	
	
	chaintypes = ['H','L']
	allVanarciinput = {x:[] for x in chaintypes}
	allVanarcioutput = {}
	bannedlist = ['1a14','3u4e','4k3d','4k3e','4ocs','4od1','4od3','4org','4ud3','5drn','5dt1','5e99','5ihu','5ijv','5ilt','5uxq','7fab'] # TODO: optional banned list
	
	db.set_numbering_scheme(args.scheme)
	db.set_region_definition(args.definition)
	a.numbering_scheme = args.scheme
	a.definition = args.definition
	numnewchain = 0
	print('Screening through SAbDab to find antibody structures with resolution <= %1.1f angstrom' %struc_cutoff)
	# import Antibodies database (all)	
	if os.path.exists(os.path.join(dbdir,'Abs.obj')):
		with open(os.path.join(dbdir,'Abs.obj'),'r') as f:
			Abs = pickle.load(f)
	else:
		Abs = {}
	for entry in sorted(db):
		if entry in bannedlist: continue
		p = db.fetch(entry)

		# ignore if resolution > struc_cutoff
		try:
			if p.get_resolution().isalpha() or float(p.get_resolution()) > struc_cutoff: continue
		except ValueError:
			continue
		for tempAb in p.fabs:
			if tempAb.VH != 'NA': 
				Abid = '_'.join((tempAb.id[0],tempAb.VH))
				allVanarciinput['H'].append((Abid,tempAb.get_sequence().split('/')[0]))
				Abs[(Abid,'H')]=tempAb
				numnewchain += 1
			if tempAb.VL != 'NA': 
				Abid = '_'.join((tempAb.id[0],tempAb.VL))
				allVanarciinput['L'].append((Abid,tempAb.get_sequence().split('/')[1]))
				Abs[(Abid,'L')]=tempAb
				numnewchain += 1
	if numnewchain == 0: return 1
	if os.path.exists(dbdir) == False:
		os.mkdir(dbdir)
	# Dump to save time
	with open(os.path.join(dbdir,'Abs.obj'),'wb') as f:
		pickle.dump(Abs,f)
	print('Found {0} high-quality antibodies'.format(len(Abs)))

	# Run anarci on all Abs to number the Abs
	print('Running anarci')
	allVanarcioutput = anarci_fun(allVanarciinput,args.scheme)

	# Dump to save time
	with open(os.path.join(dbdir,'allVanarcioutput.obj'),'wb') as f:
		pickle.dump(allVanarcioutput,f)
	
	print('Building CDR structures database')

	for cdr in cdrs:
		cdrdbdir = os.path.join(dbdir,cdr)
		if os.path.exists(cdrdbdir) == False:
			os.mkdir(cdrdbdir)
	maxproc = 4 
	chunksize = int(np.ceil(float(len(Abs))/float(maxproc)))
	joblist = [(ab,Abs[ab],allVanarcioutput,args) for ab in sorted(Abs)]
	print("Creating %d (non-daemon) workers and jobs in main process." %(min(maxproc,len(joblist))))
	pool = MyPool(min(maxproc,len(joblist)))
	results = pool.map(write_strucs, joblist)
	pool.close()
	pool.join()
	return 0
	

def write_strucs(xxx_todo_changeme):
	(Abid,Abinfo,allVanarcioutput,args) = xxx_todo_changeme
	problems = []
	
	try:
		structure = Abinfo.get_structure(scheme=args.scheme,definition=args.definition)
		for cdr in cdrs:

			if Abid[1] != cdr[-2]: continue
			outfdir = os.path.join(args.dbdir,'%s/%s.pdb'%(cdr,Abid[0]))
			if os.path.exists(outfdir): continue
			
			try:
				# ignore if Ab has missing residues in the CDR loops
				a = Accept()
				a.numbering_scheme = args.scheme
				a.definition = args.definition
				a.set_regions(['CDR'+cdr])
				
				try:
					missingr = Abinfo.get_missing_residues(scheme=args.scheme)
				except IndexError:
					missingr = []
				if sum([[a.accept(num,chain) for num in missingr[chain]].count(1) for chain in missingr]) != 0:	continue
				cdrseq,anchor = getnumberedCDRloop_builddb(list(allVanarcioutput[cdr[-2]][Abid[0]].items()),'CDR'+cdr,args.scheme,args.definition)
				if cdrseq == []: continue
				
				startpos = sorted(cdrseq)[0][0][0]
				endpos = sorted(cdrseq)[len(cdrseq)-1][0][0]
				lanchor = anchor[0][0][0]
				uanchor = anchor[1][0][0]
				residuelist = loop_parse(structure, lanchor, startpos, endpos, uanchor, cdr, Abid[0][-1],args.bfactor) # Higher quality loops only.
				if residuelist == [] or len(residuelist) != len(cdrseq)+10: 
					continue
		
				# Write to pdb
				pdboutf = open(outfdir,'wb')
				n = 1
				for eachresidue in residuelist:
					aapdb, n = eachresidue._get_output_string(select_all(), n)
					pdboutf.write(aapdb)
				pdboutf.write('TER\nEND\n')
				pdboutf.close()
				if os.path.getsize(outfdir) == 0: os.system('rm -rf {0}'.format(outfdir))
			except Exception as e: 
				print((e), Abid[0], cdr)
				pass
	except Exception as e: 
		print((e), Abid[0])
		problems.append(Abid[0])
		
	return 0

def clustering(args):
	# Compile the DTW algorithm to calculate the structural differences between loops
	cmd = []
	if 'ReadDir_DTW.o' not in os.listdir('builddb_code'):
		if args.arma_inc!='':
			cmd.append('g++ --std=c++11 builddb_code/ReadDir_DTW.cpp -o builddb_code/ReadDir_DTW.o -I{0} -lboost_system -lboost_filesystem -DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_DONT_USE_WRAPPER -llapack -lblas'.format(args.arma_inc))
		else:
			cmd.append('g++ --std=c++11 builddb_code/ReadDir_DTW.cpp -o builddb_code/ReadDir_DTW.o -lboost_system -lboost_filesystem -DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_DONT_USE_WRAPPER -llapack -lblas')

	if os.path.exists('%s/results/'%(args.dbdir)) == False:os.mkdir('%s/results/'%(args.dbdir))

	# Parallelise the process of calculating DTW-RMSD
	
	for cdr in cdrs:
		if cdr != cdrs[-1]:
			cmd.append('(./builddb_code/ReadDir_DTW.o {1}/{0} {1}/results/{0}filelist {1}/results/{0}distmatrix ; python builddb_code/DTW_vis_fixedcut-offs.py --directory={1}/{0} --filelist={1}/results/{0}filelist --distfile={1}/results/{0}distmatrix )&'.format(cdr,args.dbdir))
		else: # Last one
			cmd.append('(./builddb_code/ReadDir_DTW.o {1}/{0} {1}/results/{0}filelist {1}/results/{0}distmatrix ; python builddb_code/DTW_vis_fixedcut-offs.py --directory={1}/{0} --filelist={1}/results/{0}filelist --distfile={1}/results/{0}distmatrix );'.format(cdr,args.dbdir))
	for cdr in cdrs:
		cmd.append('mv ' +os.path.join(args.dbdir,'results/Cluster_data%s_%sA.txt'%(cdr,cutoffs[cdr])) + ' ' + os.path.join(args.dbdir,'results/Cluster_data%s_%sA_%s_%s.txt'%(cdr,cutoffs[cdr],lastm.strftime("%Y"),lastm.strftime("%m"))))
	# bash file
	bashfn = 'compile_and_run_{}.sh'.format(time.time())
	with open(bashfn,'wb') as bashf:
		bashf.write('\n'.join(cmd))
	sprocess = subprocess.Popen(['bash',bashfn], stdout=subprocess.PIPE)
	sprocess.wait()
	os.remove(bashfn)

def extractcluster(args):

	print('Loading length-independent CDR canonical form clusters')

	folderdir = args.dbdir

	Ablicluster = {cdr:{} for cdr in cdrs}
	Abtempstore = {}
	clustercenters = {cdr:{} for cdr in cdrs}
	clusterabstore = {cdr:{} for cdr in cdrs}
	clusterloopstore = {cdr:{} for cdr in cdrs}
	
	clusterout = {cdr:{} for cdr in cdrs}
	pssm = {cdr:{} for cdr in cdrs}
	with open(os.path.join(folderdir,'allVanarcioutput.obj'),'rb') as f:
		allVanarcioutput = pickle.load(f) # format allVanarcioutput[chain][Ab-chainids]
	
	for cdr in cdrs:
		
		chain = cdr[0]
		ClusterInfo = pd.read_csv(os.path.join(folderdir,'results/Cluster_data%s_%sA_%s_%s.txt'%(cdr,cutoffs[cdr],lastm.strftime("%Y"),lastm.strftime("%m"))),sep='\t')
		allclusters = set(ClusterInfo['Cluster'])
		clusterabstore[cdr].update({'UC':[]})
		for clno in allclusters:
			tempclustername = cdr + '-' + str(clno) if clno != -1 else 'UC'
			allloops = [Ab.split('.')[0] for Ab in list(ClusterInfo[ClusterInfo['Cluster']==clno]['Loop_name'])]
			lengthstore = []
			if tempclustername == 'UC':
				clusterabstore[cdr]['UC']+=allloops
			else:
				clustercenter = [Ab.split('.')[0] for Ab in list(ClusterInfo[(ClusterInfo['Cluster']==clno) & (ClusterInfo['Cluster_center'] == 'Yes')]['Loop_name'])]
				clustercenters[cdr][tempclustername]=clustercenter[0]
				clusterabstore[cdr].update({tempclustername:allloops})
				for Abid in allloops:
					if Abid not in allVanarcioutput[chain]:
						print(Abid, ' is not found in the updated entry')
						continue
					if tempclustername not in clusterloopstore[cdr]: 
						clusterloopstore[cdr].update({tempclustername:{}})
						
					if Abid not in clusterloopstore[cdr][tempclustername]:		
						[loop,_] = getnumberedCDRloop(list(allVanarcioutput[chain][Abid].items()),cdr,args.scheme,args.definition)
						if sum([identity({chain:dict(loop)},{chain:dict(clusterloopstore[cdr][tempclustername][Abother])},scheme=args.scheme,definition=args.definition) == 1.00 for Abother in clusterloopstore[cdr][tempclustername]]) != 0: continue # skip if not unique
						clusterloopstore[cdr][tempclustername].update({Abid:loop}) # only store unique loops
						if len(loop) not in lengthstore: lengthstore.append(len(loop))
				if len(clusterloopstore[cdr][tempclustername]) < 6:
					clusterloopstore[cdr].pop(tempclustername,None)
					clustercenters[cdr].pop(tempclustername,None)
					clusterabstore[cdr].pop(tempclustername,None)
					tempclustername = 'UC'
					clusterabstore[cdr][tempclustername]+=allloops
			if tempclustername not in clusterout[cdr]: clusterout[cdr].update({tempclustername:{'Cluster members':[],'Lengths':[]} })
			clusterout[cdr][tempclustername]['Cluster members']+=allloops
			clusterout[cdr][tempclustername]['Lengths']=lengthstore


		# Rename clusters
		clusters = {}
		for cl in clusterout[cdr]:
			if cl == 'UC': continue
			lengths = str(sorted(clusterout[cdr][cl]['Lengths']))[1:-1].replace(' ','')
			if lengths not in clusters: clusters.update({lengths:{}})
			clusters[lengths].update({cl:len(clusterloopstore[cdr][cl])})
		for lengths in clusters:
			asciiindex = 0
			for cl,size in sorted(list(clusters[lengths].items()), key=operator.itemgetter(1),reverse=True):
				newcl = '-'.join([cdr,lengths,ascii_uppercase[asciiindex]])
				clusterabstore[cdr][newcl]=clusterabstore[cdr].pop(cl)
				clusterloopstore[cdr][newcl]=clusterloopstore[cdr].pop(cl)
				clusterout[cdr][newcl]=clusterout[cdr].pop(cl)
				clustercenters[cdr][newcl]=clustercenters[cdr].pop(cl)
				asciiindex+=1

		# Build PSSM
		for clustername in list(clusterabstore[cdr].keys()):
			if clustername == 'UC': continue 
			storepos = {} # local storage for the cluster
			lengthstore = clusterout[cdr][clustername]['Lengths']
			
			for Abid in clusterloopstore[cdr][clustername]:
				for pos,residue in clusterloopstore[cdr][clustername][Abid]:
					residuei = resns.index(residue)
					if pos not in storepos: storepos[pos]=np.zeros(20)
					storepos[pos][residuei] += float(1)	
		
			pssm[cdr].update({clustername:{'PSSM':get_PSSM(storepos),'Lengths':lengthstore}})
		if cdr == 'L2':
			clsizes = {n:len(clusterloopstore[cdr][n]) for n in clusterloopstore[cdr]}
			largestcl = [n for n,v in list(clsizes.items()) if v == max(clsizes.values())][0]
			for cl in list(clusterout[cdr].keys()):
				if cl != largestcl and cl != 'UC': 
					clusterout[cdr]['UC']['Cluster members']+=clusterout[cdr][cl]['Cluster members']
					clusterout[cdr].pop(cl,None)
					pssm[cdr].pop(cl,None)
					clustercenters[cdr].pop(cl,None)

	with open(os.path.normpath(os.path.join(folderdir+'/ClusterMember.json')),'wb') as f: # a little weird hack
		json.dump(clusterout,f)


	newdbname = os.path.join(scalop_path,'database','{0}_{1}_v{2}-{3}.pickle'.format(args.scheme,args.definition,lastm.strftime("%Y"),lastm.strftime("%m")))
	with open(newdbname,'wb') as f:
		pickle.dump((pssm,clustercenters,clusterout),f)
	
	# Cluster analysis
	print('CDR \t Cluster \t Size')
	for cdr in clusterout:
		size = 0
		for cl in clusterout[cdr]:
			size += len(clusterout[cdr][cl]['Cluster members'])
			print('\t'.join([cdr,cl,str(len(clusterout[cdr][cl]['Cluster members']))]))
		print('\t'.join([cdr,'Total',str(size)]))
	return 0

def update(args):
	args.dbdir = os.path.join(args.dbdir,'{0}_{1}'.format(args.scheme,args.definition)) 
	try:
		sys.path.append(args.sabdabdir)
		sys.path.append(args.sabdabpydir)
	except subprocess.CalledProcessError: 
		if input('Have you installed SAbDAb?'):
			args.sabdabdir = input('Please enter the location of SAbDab bin script (e.g. ~/bin/):')
			args.sabdabpydir = input('Please enter the location of ABDB python module (e.g. ~/bin/Python/ABDB):')
			sys.path.append(args.sabdabdir)
			sys.path.append(args.sabdabpydir)
		else:
			print('Please install SAbDab before performing the update!')
			sys.exit(1)
	if args.sabdabu == 'yes': os.system(os.path.join(args.sabdabdir,'SAbDab') + ' -u')
	global identity, db, select_all  
	from ABDB.AB_Utils.calculations import identity
	from ABDB.AB_Utils import Accept
	from ABDB import database as db
	from ABDB.AbPDB.Select import select_all
	newchain = runbuilddb(args)
	if newchain == 1: 
		print('No new structure to be processed')
	else:
		clustering(args)
		extractcluster(args)
	
