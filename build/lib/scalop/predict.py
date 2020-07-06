# Assign canonical forms based on sequence-based length-independent clustering results

import numpy as np
import json, pickle, os, time, datetime, re, sys, multiprocessing
from .anarci import run_anarci
from .utils import getnumberedCDRloop, resns, SimpleFastaParser, write_txt, print_result
from .LoopGraft import export_structure

def ip(scheme,definition,dbv):

    global outm
    global clustercenters
    scalop_path  = os.path.split(__file__)[0] # from scalop database
    
    tlist = []
    if dbv == 'latest':
        for n in os.listdir(os.path.join(scalop_path,'database')):
            if re.match(r'{0}_{1}_v(\d+)-(\d+).pickle'.format(scheme,definition),n) == None: continue
            nt = datetime.datetime.strptime(''.join(re.match(r'{0}_{1}_v(\d+)-(\d+).pickle'.format(scheme,definition),n).groups()),'%Y%m')
            tlist.append((nt,n))
        if tlist == []: 
            sys.stderr.write('Database is missing. Aborting!\n')
            sys.exit(1)
            
        fname = tlist[tlist.index(max(tlist))][1]
    else:
        fname = '{0}_{1}_v{2}.pickle'.format(scheme,definition,dbv)
    if not os.path.exists(os.path.join(scalop_path,'database',fname)):
        sys.stderr.write('Database {0} not found. Please review the database directory.'%dbv)
        sys.exit(1)
    with open(os.path.join(scalop_path,'database',fname), 'rb') as f:
        outm,clustercenters,_ = pickle.load(f, encoding='latin1') 
def _score_licanonicalliPSSM(sequence, pssm):
    score = float(0)
    for pos,res in sequence:
        if pos not in pssm or pssm[pos][resns.index(res.upper())] == np.nan: continue
        
        score += pssm[pos][resns.index(res)]
    return float(score)/float(len(sequence))
    
assignmentthresholds={'H1':0.5,'H2':-0.5,'L1':-0.5,'L2':-1,'L3':-1} # legacy from 2017 PSSM
cdrs = {'H':['H1','H2'],'L':['L1','L2','L3']}

def _assign(sequence, cdr):
    assignmentthreshold =assignmentthresholds[cdr]
    loopseq = ''.join([x[1] for x in sequence]).upper()
    cdr = cdr.upper()
    if any([n not in resns for n in loopseq]): return [cdr,'None','None','None']
    seqlen = len(sequence)
    _clusters = [cluster for cluster in sorted(outm[cdr]) if seqlen in outm[cdr][cluster]['Lengths'] ] # clusters that have the sequence length
    # assign all length-8 L2 sequences (North et al. definition; or loops of the majority length(s) ) to 1 cluster
    if cdr == 'L2':
        charac_seqlen = list(outm['L2'].values())[0]['Lengths']
        if seqlen in charac_seqlen:
            ID = list(outm['L2'].keys())[0]
            center = clustercenters[cdr][ID]
            return [cdr,loopseq, ID, center]
        else:
            return [cdr,loopseq, 'None','None']
    scores = [ _score_licanonicalliPSSM(sequence, outm[cdr][cluster]['PSSM']) for cluster in sorted(outm[cdr]) if seqlen in outm[cdr][cluster]['Lengths'] ]

    if any(np.greater(np.array(scores),assignmentthreshold)):
        pos = scores.index(max(scores))
        ID = _clusters[pos]
        center = clustercenters[cdr][ID]
        return [cdr,loopseq, ID, center]
    else: # all were 0
        return [cdr,loopseq, 'None', 'None']

def assignweb(args):

    seqs = {}
    permitbuildstruct = False
    if os.path.isfile(args.seq):
        with open(args.seq) as f:
            _seqs = list(SimpleFastaParser(f))
        seqs = []
        for seqid,seq in _seqs:
            if '/' in seq:
                seqh,seql = seq.split('/')
                seqh = re.sub('[\W]+','',seqh)
                seql = re.sub('[\W]+','',seql)
                seqs.append((seqid+'_1',seqh))
                seqs.append((seqid+'_2',seql))
            else:
                seqs.append((seqid,re.sub('[\W\/]+','',seq)))
    else: # Single sequence input
        permitbuildstruct = True
        if '/' in args.seq: # must be Heavy-chain/Light-chain 
            seqh,seql = args.seq.split('/')
            seqs = [('input_1',seqh),('input_2',seql)]
        else:
            seqs = [('input',args.seq)]
    
    assignresults = [{} for i in range(len(seqs))]

    # Number the heavy / light chain using ANARCI
    ncpu = multiprocessing.cpu_count()
    numberedseqs = run_anarci( seqs ,scheme=args.scheme,ncpu=ncpu,assign_germline=False) 
    
    # Import the right version of pssm database
    ip(args.scheme,args.definition,args.dbv)    
    
    for i in range(len(numberedseqs[0])):
        seqname = numberedseqs[0][i][0]
        assert(seqname == seqs[i][0])
        assignresults[i].update({'seqname':seqname})
        assignresults[i].update({'input':seqs[i]})
        assignresults[i].update({'outputs':{}})
        # Find if ANARCI can number this chain (whether or not it is an antibody chain)
        if numberedseqs[2][i] == None: continue
        # Find whether heavy or light chain
        chain = 'H' if numberedseqs[2][i][0]['chain_type'] == 'H' else 'L'

        for cdr in cdrs[chain]:
            # Frame the CDR region
            loop,_ = getnumberedCDRloop(numberedseqs[1][i][0][0],cdr,args.scheme,args.definition) 
    
            # Assign canonical form
            assignresults[i]['outputs'].update({cdr:_assign(loop,cdr)})
    opf = write_txt(assignresults, args.outputdir, args.scheme, args.definition)
    outpdbf = False

    if permitbuildstruct and args.structuref!='':
    	_outpdbf = opf.replace('.txt','.pdb')
    	graftedstructs,is_model_generated = export_structure(args, assignresults, _outpdbf)
    	if is_model_generated == 1: 
    		opf = write_txt(assignresults, args.outputdir, args.scheme, args.definition, graftedstructs)
    		outpdbf = _outpdbf
    return assignresults,opf,outpdbf
    
def assigncmd(args):

    seqs = {}
    permitbuildstruct = False
    if os.path.isfile(args.seq):
        with open(args.seq) as f:
            _seqs = list(SimpleFastaParser(f))
        seqs = []
        for seqid,seq in _seqs:
            if '/' in seq:
                seqh,seql = seq.split('/')
                seqh = re.sub('[\W]+','',seqh)
                seql = re.sub('[\W]+','',seql)
                seqs.append((seqid+'_1',seqh))
                seqs.append((seqid+'_2',seql))
            else:
                seqs.append((seqid,re.sub('[\W\/]+','',seq)))
    else: # Single sequence input
        permitbuildstruct = True
        if '/' in args.seq: # must be Heavy-chain/Light-chain 
            seqh,seql = args.seq.split('/')
            seqs = [('input_1',seqh),('input_2',seql)]
        else:
            seqs = [('input',args.seq)]
    
    assignresults = [{} for i in range(len(seqs))]

    # Number the heavy / light chain using ANARCI
    ncpu = multiprocessing.cpu_count()
    numberedseqs = run_anarci( seqs ,scheme=args.scheme,ncpu=ncpu,assign_germline=False)  
    
    # Import the right version of pssm database
    ip(args.scheme,args.definition,args.dbv)    
    
    for i in range(len(numberedseqs[0])):
        seqname = numberedseqs[0][i][0]
        assert(seqname == seqs[i][0])
        assignresults[i].update({'seqname':seqname})
        assignresults[i].update({'input':seqs[i]})
        assignresults[i].update({'outputs':{}})
        # Find if ANARCI can number this chain (whether or not it is an antibody chain)
        if numberedseqs[2][i] == None: continue
        # Find whether heavy or light chain
        chain = 'H' if numberedseqs[2][i][0]['chain_type'] == 'H' else 'L'
        
        for cdr in cdrs[chain]:
            # Frame the CDR region
            loop,_ = getnumberedCDRloop(numberedseqs[1][i][0][0],cdr,args.scheme,args.definition) 
    
            # Assign canonical form
            assignresults[i]['outputs'].update({cdr:_assign(loop,cdr)})
    
    # Output format
    if args.outputdir == False:
        print_result(assignresults, args.scheme, args.definition)
    elif args.outputdir == '':
        pass # No output
    elif args.outputformat in ['txt','csv']:
        opf = write_txt(assignresults, args.outputdir, args.scheme, args.definition)
        if permitbuildstruct and args.structuref!='':
            outpdbf = opf.replace('.txt','.pdb')
            graftedstructs,is_model_generated = export_structure(args, assignresults, outpdbf)
            if is_model_generated: opf = write_txt(assignresults, args.outputdir, args.scheme, args.definition, graftedstructs)
        print('Write results to %s' % opf)
    elif args.outputformat =='json':
        opf = os.path.join(args.outputdir,'{0}.{1}'.format(time.time(),args.outputformat))
        if not os.path.exists(args.outputdir): os.mkdir(args.outputdir) 
        with open(opf,'wb') as f:
            json.dump(assignresults,f)
        print('Write results to %s' % opf)
    return assignresults
    
def assign(seq,scheme='imgt',definition='north',dbv='latest',structuref='',loopdb='',hc='',lc='',blacklist=[]):
    '''
    args: <sequence(s)> <numbering scheme> <cdr definition> <db version> <structure file> <loop database directory> <heavy chain ID> <light chain ID> <blacklisted PDB_CHAIN IDs>
    '''

    seqs = []
    permitbuildstruct = False
    if type(seq) == str and os.path.isfile(seq):
        with open(seq) as f:
            _seqs = list(SimpleFastaParser(f))
        
        for seqid,seq in _seqs:
            if '/' in seq:
                seqh,seql = seq.split('/')
                seqh = re.sub('[\W]+','',seqh)
                seql = re.sub('[\W]+','',seql)
                seqs.append(('{}_1'.format(seqi),seqh))
                seqs.append(('{}_2'.format(seqi),seql))
            else:
                seqs.append((seqid,re.sub('[\W\/]+','',seq)))
    elif type(seq) == dict: # in dict
        for seqn in seq:
            if '/' in seq[seqn]: # must be Heavy-chain/Light-chain 
                seqh,seql = seq[seqn].split('/')
                seqs.append(('{}_1'.format(seqn),seqh))
                seqs.append(('{}_2'.format(seqn),seql))
            else:
                seqs.append((seqn,seq[seqn]))

    elif type(seq) == list: # in unnamed list
        for seqi in range(len(seq)):
            if '/' in seq[seqi]: # must be Heavy-chain/Light-chain 
                seqh,seql = seq[seqi].split('/')
                seqs.append(('{}_1'.format(seqi),seqh))
                seqs.append(('{}_2'.format(seqi),seql))
            else:
                seqs.append((str(seqi),seq[seqi]))

    else: # Single sequence input
        permitbuildstruct = True
        if '/' in seq: # must be Heavy-chain/Light-chain 
            seqh,seql = seq.split('/')
            seqs = [('input_1',seqh),('input_2',seql)]
        else:
            seqs = [('input',seq)]
    assignresults = [{} for i in range(len(seqs))]
    # Number the heavy / light chain using ANARCI

    ncpu = multiprocessing.cpu_count()

    numberedseqs = run_anarci( seqs,scheme=scheme,ncpu=ncpu,assign_germline=False)
    
    # Import the right version of pssm database
    ip(scheme,definition,dbv)       
            
    for i in range(len(numberedseqs[0])):
        # Find if ANARCI can number this chain (whether or not it is an antibody chain)
        if numberedseqs[2][i] == None: continue
        # Find whether heavy or light chain
        chain = 'H' if numberedseqs[2][i][0]['chain_type'] == 'H' else 'L'
        seqname = numberedseqs[0][i][0]
        assert(seqname == seqs[i][0])
        assignresults[i].update({'seqname':seqname})
        assignresults[i].update({'input':seqs[i]})
        assignresults[i].update({'outputs':{}})
        for cdr in cdrs[chain]:
            # Frame the CDR region
            loop,_ = getnumberedCDRloop(numberedseqs[1][i][0][0],cdr,scheme,definition) 
    
            # Assign canonical form
            assignresults[i]['outputs'].update({cdr:_assign(loop,cdr)})
    if permitbuildstruct and structuref!='':
        from collections import namedtuple
        args = namedtuple('args',['scheme','definition','dbv','structuref','loopdb','hc','lc','blacklist'])(scheme=scheme,definition=definition,dbv=dbv,structuref=os.path.abspath(structuref),loopdb=loopdb,hc=hc,lc=lc,blacklist=blacklist)
        graftedstructs, is_model_generated = export_structure(args, assignresults, '')
        
        for i in range(len(assignresults)):
            if 'outputs' not in assignresults[i]: continue
            for cdr in assignresults[i]['outputs']:
                if cdr in graftedstructs: 
                    assignresults[i]['outputs'][cdr]+=[graftedstructs[cdr]]
    return assignresults
