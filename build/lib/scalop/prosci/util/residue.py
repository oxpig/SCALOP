#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Residue class represents a single amino acid in a PDB file.
#
# Author: Sebastian Kelm
# Created: 07/05/2009
#
# Revisions:
# 07/05/2009    Sebastian Kelm    Moved from mescol.py into its own script, and generalised to all-atom data.
#

from array import array

import numpy

from scalop.prosci.common import join, OverridableField, SettableField
from scalop.prosci.util.pdb import Pdb, Atom, residueLetter, residueCode


def is_residue_boundary(atomlist, index):
  if index < 0:
    index = len(atomlist)+index
  assert index <= len(atomlist)
  if index == 0 or index == len(atomlist):
    return True
  return (atomlist[index].ires    != atomlist[index-1].ires) \
      or (atomlist[index].chain   != atomlist[index-1].chain) \
      or (atomlist[index].inscode != atomlist[index-1].inscode)


def iterate_residue_boundaries(atomlist):
  istart = 0
  iend   = 1
  while istart < len(atomlist):
    # search forward, for the end of the residue
    while iend < len(atomlist) and not is_residue_boundary(atomlist, iend):
      iend += 1
    yield (istart, iend)
    istart = iend
    iend += 1


class ResidueList(list):
  "A list of Residue objects."
  
  chain = SettableField("chain")
  code = OverridableField("chain")
  
  def __init__(self, pdb, code=None):
    "Takes a list or tuple or Pdb object as argument. Otherwise passes the argument to Pdb() first."
    if not pdb:
      pass
    elif isinstance(pdb, list) or isinstance(pdb, tuple):
      for r in pdb:
        self.append(r)
    else:
      if not isinstance(pdb, Pdb):
        pdb = Pdb(pdb)
      for r in pdb.xresidues():
        self.append(Residue(r))
      if pdb.code:
        self.code = pdb.code
    if code is not None:
      self.code = code
  
  
  def __repr__(self):
    return "ResidueList(%s)"%(list.__repr__(self))
  
  def __str__(self):
    out=""
    for r in self:
      if r:
        out += str(r)
    return out
  
  def __getslice__(self, start=None, end=None):
    return ResidueList(list.__getslice__(self, start, end))
  
  def __add__(self, other):
    return ResidueList(list.__add__(self, other))
  def __radd__(self, other):
    return ResidueList(list.__radd__(self, other))
  def __iadd__(self, other):
    self.extend(other)
    return self
  
  
  def get_seq(self):
    "Returns a FASTA string representation of the sequence."
    seq = array('c')
    for r in self:
      if r is not None:
        seq.append(r.get_seq())
    return seq.tostring()
  
  
  def get_gapped_seq(self, gapchar='-', collapse_gaps=False):
    "Returns a FASTA string representation of the sequence, with gaps in the residue numbering represented by special characters (default: '-')."
    seq = array('c')
    prevres = None
    for r in self:
      if r is not None:
        if prevres and r.ires > prevres.ires + 1:
          if collapse_gaps:
            seq.append(gapchar)
          else:
            for i in range(r.ires - prevres.ires - 1):
              seq.append(gapchar)
        seq.append(r.get_seq())
      else:
        seq.append(gapchar)
      prevres = r
    return seq.tostring()
  
  
  def get_coords(self, atomfilter=lambda atom: True):
    "Get a numpy array of coordinates"
    coords=[]
    for r in self:
      if r is not None:
        for a in r:
          if atomfilter(a):
            coords.append(a.xyz)
    return numpy.array(coords)
  
  
  def deep_copy(self):
    rl = ResidueList([])
    for r in self:
      if r is not None:
        rl.append(r.copy())
      else:
        rl.append(None)
    return rl
  
  
  def get_gap_list(self, ignore_residue_numbers=False):
    """Returns a list of gaps, as determined using residue numbers. None elements count as missing residues.
    
    Format: [[index, gaplength], ...] , where index is the array index of the gap and gaplength is the number of missing residues."""
    results=[]
    
    start=-1
    for i, r in enumerate(self):
      if r is None:
        if start < 0:
          start = i
      elif start >= 0:
        results.append([start, i])
        start = -1
    if start >= 0:
      results.append([start, len(self)])
    
    if not ignore_residue_numbers:
      len_before = len(results)
      for i in range(1, len(self)):
        if self[i-1] is not None and self[i] is not None:
          gaplen = self[i].CA.ires - self[i-1].CA.ires - 1
          if gaplen > 0:
            results.append([i, gaplen])
      if len_before and len_before != len(results):
        results.sort()
    
    return results
  
  
  def to_pdb(self, atomfilter=lambda atom: True):
    "Returns a Pdb object containing the Atom objects within this ResidueList"
    p = Pdb(self.code, [])
    for res in self:
      for atm in res:
        if atomfilter(atm):
          p.data.append(atm)
    return p
  
  
  def renumber_atoms(self):
    "Renumbers all atoms, starting from 1 upwards."
    i = 1
    for r in self:
      for a in r:
        a.iatom = i
        i += 1
  
  
  def iter_backbone(self):
    "Iterate over all backbone atoms that aren't None"
    for r in self:
      for a in r.iter_backbone():
        yield a
  
  
  def split_chains(self):
    "Returns a list of ResidueList objects, split at the positions where chain codes differ between subsequent residues."
    chains = [ResidueList([])]
    for r in self:
      if chains[-1] and chains[-1][-1].chain != r.chain:
        chains.append(ResidueList([]))
      chains[-1].append(r)
    for c in chains:
      c.code = self.code
    return chains
  
  
  def remove_residues(self, func):
    i=0
    while i<len(self):
      if func(self[i]):
        del self[i]
      else:
        i += 1
  
  
  def find_residue(self, ires, inscode="", chain=""):
    if isinstance(ires, Residue):
      inscode = ires.inscode
      chain = ires.chain
      ires = ires.ires
    
    for i, res in enumerate(self):
      if (res.ires == ires) and (inscode == res.inscode):
        if (not chain) or (chain == res.chain):
          return i
    return -1



class Residue(object):
  "Class representing a single amino acid. Has fields for all main chain and CB atoms. Remaining atoms are in the list self.rest"

  res = SettableField("res")
  ires = SettableField("ires")
  inscode = SettableField("inscode")
  chain = SettableField("chain")
  hetatm = SettableField("hetatm")
  
  
  def __init__(self, atoms=None):
    "Expects a list of Atom objects as input, or a Pdb object."
    
    self.N = None
    self.CA = None
    self.C = None
    self.O = None
    self.CB = None
    self.rest = []
    
    if atoms is not None:
      for a in atoms:
        if a is None:
          continue
        self.add(a)
      
      if self.res in ("GLY", "ALA", "PRO"):
        del self.rest[:]
        if self.res == "GLY":
          self.CB = None
  
  def __bool__(self):
    return bool(self.N or self.CA or self.C or self.O or self.CB or self.rest)
  
  def __iter__(self):
    if self.N is not None:
      yield self.N
    if self.CA is not None:
      yield self.CA
    if self.C is not None:
      yield self.C
    if self.O is not None:
      yield self.O
#    if self.CA.res != "GLY":
    if self.CB is not None:
      yield self.CB
    for a in self.rest:
      yield a
  
  
  def add(self, a):
    if not isinstance(a, Atom):
      a = Atom(a)
    if a.atom in ("N", "CA", "C", "O", "CB"):
      curratm = getattr(self, a.atom)
      # Only save atom with lowest altloc
      if curratm is None or curratm.altloc < a.altloc:
        setattr(self, a.atom, a)
    else:
      # TODO: do something about multiple altlocs
      self.rest.append(a)
  
  
  def iter_backbone(self):
    "Iterate over all backbone atoms that aren't None"
    if self.N is not None:
      yield self.N
    if self.CA is not None:
      yield self.CA
    if self.C is not None:
      yield self.C
    if self.O is not None:
      yield self.O
  
  
  def __str__(self):
    s=array('c')
    for a in self:
      s.extend(str(a))
      s.extend("\n")
    return s.tostring()
  
  
  def __repr__(self):
    return "Residue(%s)"%(repr([self.N, self.CA, self.C, self.O, self.CB]+self.rest))
  
  
  def get_seq(self):
    "Get the amino acid letter representing this Residue"
    return residueLetter(self.res)
  
  def get_atom(self, atype):
    "Get atom of a specified type, if it exists. Else raise ValueError."
    atom = None
    if atype in ("N", "CA", "C", "O", "CB"):
      try:
        atom = getattr(self, atype)
      except AttributeError:
        pass
    else:
      for a in self.rest:
        if a.atom == atype:
          atom = a
          break
    
    if atom is not None:
      return atom
    
    raise ValueError("Residue '%s/%d%s/%s' does not contain an atom of type '%s'" % (self.chain, self.ires, self.inscode, self.res, atype))
  
  
  def set_type(self, newtype):
    """Set a new residues type.
    
    Deletes the side chain unless newtype == self.type. Keeps CB unless newtype == 'GLY'."""
    if newtype == self.res:
      rebuild_CB(self)
      return
    
    assert len(newtype) == 3
    assert newtype.isupper()
    
    if newtype == "GLY":
      self.CB = None
    del self.rest[:]
    
    #print "Relabeling residue (%s.%d%s) from %s to %s" % (self.chain, self.ires, self.inscode, self.res, newtype)
    
    self.res = newtype
    assert self.CA.res == newtype
    
    rebuild_CB(self)
  
  
  def copy(self):
    "Make a copy of this object."
    return Residue([a.copy() for a in self])
  
  def equals_in_name(other):
    return self.chain == other.chain and self.ires == other.ires and self.res == other.res and self.inscode == other.inscode
  
  def is_left_neighbour_of(other, minlen_sq=1.0, maxlen_sq=4.0):
    return is_residue_consecutive(self, other, minlen_sq, maxlen_sq)
    
  def is_right_neighbour_of(other, minlen_sq=1.0, maxlen_sq=4.0):
    return is_residue_consecutive(other, self, minlen_sq, maxlen_sq)




def is_residue_consecutive(res_a, res_b, minlen_sq=1.0, maxlen_sq=4.0):
    a = res_a.C
    b = res_b.N
    if a is None or b is None:
      return False
    d_sq = (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2
    if (d_sq < minlen_sq or d_sq > maxlen_sq):
      return False
    return True

def is_residue_numbering_consecutive(res_a, res_b):
  if res_a.chain != res_b.chain:
    return False
  if res_a.ires == res_b.ires - 1:
    return res_b.inscode in ("", "A")
  if res_a.ires == res_b.ires:
    return ord(res_a.inscode) == ord(res_b.inscode) - 1
  return False


cb_dictionary = {'CYS': [-0.3494607343363032, 0.79155529324285911, 0.49818968247044454, 1.5296357057135819], 'ASP': [-0.34996888099893109, 0.78680052199417105, 0.50545821158245419, 1.5318055409680786], 'SER': [-0.34908005233451428, 0.79092525569532701, 0.50008816225520469, 1.5298749212857365], 'GLN': [-0.35044189807274562, 0.78789868002634589, 0.50391579356824856, 1.5307608842551144], 'LYS': [-0.35014881263058767, 0.78862863785937809, 0.50338903814841085, 1.5311461087301732], 'ILE': [-0.36317297924619868, 0.78202773684664761, 0.50393645313205515, 1.5451861881417324], 'ASN': [-0.34954774900564267, 0.78213004101991257, 0.5126836623595018, 1.5321433083090938], 'THR': [-0.36088247834694348, 0.787447497345662, 0.49708105577718109, 1.5401087763522436], 'PHE': [-0.35040956754722258, 0.78912631290685253, 0.50162366982130213, 1.5339592272837324], 'ALA': [-0.34777316024597776, 0.78875689088188949, 0.50566390742614242, 1.5253623347586374], 'MET': [-0.35070773060772492, 0.78911407016265778, 0.50143527870114502, 1.530863273877535], 'HIS': [-0.3478745720232258, 0.78875770121288913, 0.50381793689150922, 1.53204409902269], 'LEU': [-0.34927313668815713, 0.7911706392471316, 0.49990908813589296, 1.5313231292876908], 'ARG': [-0.35073600249773451, 0.78838749120654794, 0.50304768909818021, 1.5309172528392483], 'TRP': [-0.3496482864883515, 0.78907538478814376, 0.50231116231957951, 1.5332689020009584], 'PRO': [-0.22910526303159487, 0.84414806305461676, 0.48310967532320609, 1.5324482100664505], 'VAL': [-0.36481633408183051, 0.78221455833922093, 0.50264044284361142, 1.5453102405120656], 'GLU': [-0.35008028783676887, 0.78714165176127382, 0.50553324134098787, 1.5316362178823235], 'TYR': [-0.34979587900927023, 0.78983362800133428, 0.5011455140152965, 1.5342036619521204]}


cb_dictionary['UNK'] = cb_dictionary['XXX'] = cb_dictionary['ALA']


def rebuild_CB(res, replace=False, ignore_glycine=True):
  "Rebuild a given Residue's missing C-beta atom co-ordinates"
  
  if not (res.CA and res.C and res.N):
    # Can't build CB
    #print "Missing main chain atoms - CA,C,N needed to rebuild CB"
    return False
  
  if res.CB and not replace:
    # There is already a CB atom present and we were told not to replace it
    #print "CB already present"
    return False
  
  if res.res == "GLY" and ignore_glycine:
    # Glycines don't have CB atoms
    #print "Ignoring glycine"
    return False
  
  CB = res.CA.copy()
  CB.atom = "CB"
  
  n_direction = res.N.xyz - res.CA.xyz
  n_direction /= numpy.linalg.norm(n_direction)
  c_direction = res.C.xyz - res.CA.xyz
  c_direction /= numpy.linalg.norm(c_direction)
  norm_direction = numpy.cross(n_direction, c_direction)
  norm_direction /= numpy.linalg.norm(norm_direction)
  norm2_direction = numpy.cross(n_direction, norm_direction)
  norm2_direction /= numpy.linalg.norm(norm_direction)
  
  n_factor, norm_factor, norm2_factor, cb_dist = cb_dictionary[residueCode(residueLetter(res.res))]
  
  cb_vector = (n_direction * n_factor) + (norm_direction * norm_factor) + (norm2_direction * norm2_factor)
  cb_vector /= numpy.linalg.norm(cb_vector)
  cb_vector *= cb_dist
  
  CB.xyz += cb_vector
  
  res.CB = CB
  return True
