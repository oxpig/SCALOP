#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Handle protein structure files in PDB format
#
# Author: Sebastian Kelm
# Created: 07/05/2009
#
# Revisions:
# 07/05/2009    Sebastian Kelm    Moved from mescol.py into its own script, and generalised to all-atom data.
#


import sys
import subprocess
import gzip
from array import array
from io import StringIO

import numpy

from scalop.prosci.common import *


# Dictionary for the conversion from 3-letter to 1-letter amino acid code
#
residues = {}
residues['ALA'] = 'A'
residues['CYS'] = 'C'
residues['ASP'] = 'D'
residues['GLU'] = 'E'
residues['PHE'] = 'F'
residues['GLY'] = 'G'
residues['HIS'] = 'H'
residues['ILE'] = 'I'
residues['LYS'] = 'K'
residues['LEU'] = 'L'
residues['MET'] = 'M'
residues['ASN'] = 'N'
residues['PRO'] = 'P'
residues['GLN'] = 'Q'
residues['ARG'] = 'R'
residues['SER'] = 'S'
residues['THR'] = 'T'
residues['VAL'] = 'V'
residues['TRP'] = 'W'
residues['TYR'] = 'Y'

residues['SEC'] = 'U' # Selenocysteine
residues['ASX'] = 'B' # ASP or ASN
residues['XLE'] = 'J' # LEU or ILE
residues['GLX'] = 'Z' # GLU or GLN
residues['XXX'] = 'X' # Unknown amino acid


residues_inverted = {}
for code in residues:
  residues_inverted[residues[code]] = code


# NOTE: Duplicate mappings
residues['MSE'] = 'M' # Selenomethionine - pretend it's just methionine
residues['XAA'] = 'X' # Unknown amino acid
residues['UNK'] = 'X' # Unknown amino acid


def residueLetter(code):
  try:
    return residues[code]
  except KeyError:
    return 'X'

def residueCode(letter):
  try:
    return residues_inverted[letter]
  except KeyError:
    return 'XXX'

def residue_is_amino_acid(code):
  return code in residues


# Global function to make sure a PDB code is ok
def validate_pdb_code(pdbcode):
    import re
    p = re.compile('[<>&\|!]')
    assert isinstance(pdbcode, str) and len(pdbcode)==4, "Invalid PDB code: "+pdbcode
    assert not p.search(pdbcode), "PDB code contains illegal characters: "+pdbcode


def get_pdb_file(pdb_code):
    import os
    from glob import glob
    try:
      import getpdb
    except ImportError:
      sys.stderr.write("ERROR importing module 'getpdb'. PDB code arguments to the Pdb() constructor will not work.\n")
      return None
    
    validate_pdb_code(pdb_code)
    filepattern = getpdb.get_pdb_path(pdb_code)
    if not glob(filepattern):
      return None
    output = subprocess.Popen('zcat '+filepattern+' | egrep "^ATOM|^HETATM"', shell=True, stdout=subprocess.PIPE).stdout
    if output == None:
      sys.stderr.write("ERROR reading coordinate file for pdb "+pdb_code+"\n")
    return output


def diffPdb(bigModel, smallModel):
  "diffPdb(bigModel, smallModel) : return the atoms present in bigModel but not smallModel"
  
  excludeMe=set([])
  for start, end in smallModel.xresidue_boundaries():
    atm = smallModel[start]
    excludeMe.add((atm.ires, atm.inscode, atm.res, atm.chain))
  
  difference = Pdb("difference", [])
  for start, end in bigModel.xresidue_boundaries():
    atm = bigModel[start]
    if (atm.ires, atm.inscode, atm.res, atm.chain) not in excludeMe:
      difference.append_atoms(bigModel[start:end])
  return difference


def intersectPdb(models):
  "intersectPdb(models) : Gives back the atoms, for each model, corresponding to the intersection of all models."
  
  if not models:
    return list()
  
  intersection = set([])
  
  for start, end in models[0].xresidue_boundaries():
    atm = models[0][start]
    intersection.add((atm.ires, atm.inscode, atm.res, atm.chain))
  
  for m in models[1:]:
    res_set=set([])
    for start, end in m.xresidue_boundaries():
      atm = m[start]
      res_set.add((atm.ires, atm.inscode, atm.res, atm.chain))
    intersection &= res_set
    
  newmodels=[]
  for m in models:
    newmodel = Pdb(m.code + "_subset", [])
    for start, end in m.xresidue_boundaries():
      atm = m[start]
      if (atm.ires, atm.inscode, atm.res, atm.chain) in intersection:
        newmodel.append_atoms(m[start:end])
    newmodels.append(newmodel)
  
  return newmodels



atom_order_list = ["N", "H", "CA", "sidechain atoms go here", "C", "O"]
def atom_type_key(a):
    try:
      i=atom_order_list.index(a.atom)
      return i
    except ValueError:
      return 3


def intersectAtomTypes(rS, rP):
    "intersectAtomTypes(rS, rP) : for both given residues, return the set of atoms present in the other"
    
    atoms_S = set([x.atom for x in rS])
    atoms_P = set([x.atom for x in rP])
    
    onlyInS = atoms_S - atoms_P
    onlyInP = atoms_P - atoms_S
    
    rS=rS[:]
    rP=rP[:]
    
    if not isinstance(rS, Pdb):
      rS = Pdb("residueS", rS)
    if not isinstance(rP, Pdb):
      rP = Pdb("residueP", rP)
    
    for atmtype in onlyInS:
      rS.delete_atoms_by_type(atmtype)
    for atmtype in onlyInP:
      rP.delete_atoms_by_type(atmtype)
    
    rS.sort(key=atom_type_key)
    rP.sort(key=atom_type_key)
    
    return rS, rP



class NoAtomCoordinatesError(ParsingError):
  pass


class Atom(object):
  "The Atom class represents a single ATOM line in a PDB file."
  
  # field names
  fields = ("iatom", "atom", "altloc", "res", "chain", "ires", "inscode", "x", "y", "z", "occup", "b", "element", "charge")
  
  class Coord(object):
    "Descriptor to access individual x,y,z coordinates"
    def __init__(self, i):
        self.i = i
    def __get__(self, obj, cls=None): 
        return obj.xyz[self.i]
    def __set__(self, obj, val): 
        obj.xyz[self.i] = val
    def __delete__(self, obj): 
        pass
  
  x = Coord(0)
  y = Coord(1)
  z = Coord(2)
  
  
  def __init__(self, line=None):
    if line != None:
      self.init(str(line))
    else:
      self.hetatm    = False
      self.iatom     = -1
      self.atom      = "XXX"
      self.altloc    = ""
      self.res       = "XXX"
      self.chain     = ""
      self.ires      = -1
      self.inscode   = ""
      self.xyz       = numpy.array([0.0, 0.0, 0.0])
      self.occup     = 1.0
      self.b         = 1.0
      self.element   = ""
      self.charge    = ""
  
  
  def init(self, line):
    if line.startswith("HETATM"):
      self.hetatm = True
    else:
      self.hetatm = False
    
    try:
      self.iatom = int(line[ 6:11].strip())
    except ValueError:
      self.iatom = -1
    
    self.atom      =       line[12:16].strip().upper()
    self.altloc    =       line[16:17].strip().upper()
    self.res       =       line[17:20].upper()
    self.chain     =       line[21:22].strip()
    self.ires      =   int(line[22:26].strip())
    self.inscode   =       line[26:27].strip()
    self.xyz       = numpy.array([0.0, 0.0, 0.0])
    try:
      self.xyz[0]    = float(line[30:38].strip())
      self.xyz[1]    = float(line[38:46].strip())
      self.xyz[2]    = float(line[46:54].strip())
    except ValueError:
      self.xyz     = None
    
    self.occup     = 1.0
    self.b         = 1.0
    self.element   = ""
    self.charge    = ""
    try:
      self.occup     = float(line[54:60].strip())
      self.b         = float(line[60:66].strip())
      self.element   =       line[76:78].strip()
      self.charge    =       line[78:80].strip()
    except ValueError:
      pass
    
    if self.xyz is None:
      raise NoAtomCoordinatesError()
  
  
  def copy(self, other=None):
    "If other is None, returns a copy of this Atom. Otherwise, changes this Atom by taking on the values of other."
    if other is None:
      a = Atom()
      a.hetatm    = self.hetatm
      a.iatom     = self.iatom
      a.atom      = self.atom
      a.altloc    = self.altloc
      a.res       = self.res
      a.chain     = self.chain
      a.ires      = self.ires
      a.inscode   = self.inscode
      a.xyz       = self.xyz.copy()
      a.occup     = self.occup
      a.b         = self.b
      a.element   = self.element
      a.charge    = self.charge
      return a
    else:
      self.hetatm    = other.hetatm
      self.iatom     = other.iatom
      self.atom      = other.atom
      self.altloc    = other.altloc
      self.res       = other.res
      self.chain     = other.chain
      self.ires      = other.ires
      self.inscode   = other.inscode
      self.xyz       = other.xyz.copy()
      self.occup     = other.occup
      self.b         = other.b
      self.element   = other.element
      self.charge    = other.charge
  
  
  def is_amino(self):
    return residue_is_amino_acid(self.res)
  
  def is_main_chain(self):
    return self.atom in ("CA", "C", "N", "O")
  
  
  # Make Atom object seem like a list
  #
  def __iter__(self):
      yield self.iatom
      yield self.atom
      yield self.altloc
      yield self.res
      yield self.chain
      yield self.ires
      yield self.inscode
      yield self.x
      yield self.y
      yield self.z
      yield self.occup
      yield self.b
      yield self.element
      yield self.charge
  
  def __getitem__(self, i):
      if isinstance(i, int):
          return self.__dict__[Atom.fields[i]]
      else:
          return self.__dict__[str(i)]
  
  def __len__(self):
      return len(Atom.fields)
  
  def __str__(self):
      if self.hetatm:
        datatype = "HETATM"
      else:
        datatype = "ATOM  "
      if self.atom[0:1].isdigit():
        atom = "%-4s" % (self.atom)
      else:
        atom = "%1s%-3s" % (self.atom[-4:-3], self.atom[-3:])
      return "%s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%-2s" % (datatype, self.iatom, atom, self.altloc, self.res, self.chain, self.ires, self.inscode, self.x, self.y, self.z, self.occup, self.b, self.element, self.charge)
  def __repr__(self):
      return 'Atom("%s")' % (str(self))
  
  def __cmp__(self, other):
      v = cmp(type(self), type(other))
      if v != 0:
        return v
      v = cmp((self.chain, self.iatom, self.altloc), (other.chain, other.iatom, other.altloc))
      if v != 0:
        return v
      return cmp((self.atom, self.res, self.ires, self.inscode, self.x, self.y, self.z, self.occup, self.b, self.element, self.charge), (other.atom, other.res, other.ires, other.inscode, other.x, other.y, other.z, other.occup, other.b, other.element, other.charge))
  
  def __hash__(self):
      return hash(self.chain) ^ hash(self.res) ^ hash(self.iatom) ^ hash(self.altloc)



class Pdb(object):
    "The Pdb class represents a protein structure as parsed from a *.pdb or *.atm file."
    
    def __init__(self, pdb, data=None, keepBefore1=True, deleteInsertions=False, nofilter=False, cryst1=""):
        self.parent = None
        self.code   = ""
        self.data   = []
        self.ligands= [] # hetero atoms
        self.cryst1 = cryst1
        
        # indeces for faster data access
        self.chainlimits = None
        self.hashindex = None
        
        if not pdb and not data:
          # We're a dummy object. Don't bother trying to parse anything.
          return
        
        if pdb is None:
          pdb = ""
        
        if isinstance(pdb, Pdb):
          self.code   = pdb.code
          self.parent = pdb
          if data is not None:
            self.data   = [d for d in data]
          else:
            self.data   = [d for d in pdb.data]
        elif isinstance(pdb, str):
            if data is None:
              self.code = basename_noext(pdb)
              if len(self.code) >= 7 and self.code.startswith("pdb"):
                self.code = self.code[3:]
              
              try:
                #p = subprocess.Popen(["zcat", pdb], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                #f = StringIO(p.communicate()[0])
                #if p.returncode != 0:
                #  raise IOError()
                f = gzip.open(pdb)
                try:
                  self.parse_file(f)
                finally:
                  f.close()
                self.code = basename_noext(self.code)
              except IOError:
                f = open(pdb)
                try:
                  self.parse_file(f)
                finally:
                  f.close()
            else:
              self.code = pdb
              self.parse_file(data)
        else:
          if data is not None:
            raise ArgumentError("Pdb() takes exactly 2 arguments (name, data), or a single argument of class Pdb, or a valid filename!")
          data = pdb
          pdb = ""
          self.parse_file(data)
        
        self.unindex()
        
        if data and not self.data:
          raise ParsingError("Given PDB data/file was not successfully parsed: "+repr(pdb))
        
        if self.parent is None and not nofilter:
          self.filter(negres=(not keepBefore1), inscode=deleteInsertions)
          """
          self.filter_missing_CA()
          if not keepBefore1: 
            self.filter_residues_before_1()
          self.filter_duplicate_unlabelled_atoms()
          self.filter_insertions()
          self.filter_alternative_locations()
          """
        
        if not self.cryst1 and self.parent:
          root = self.get_root_object()
          self.cryst1 = root.cryst1
    
    
    # Make Pdb object seem like a list - content is the atoms list
    #
    def __iter__(self):
        return self.data.__iter__()
    
    def __getitem__(self, i):
        if isinstance(i, tuple):
          if not len(i)==2:
            raise TypeError("Pdb.__getitem__(i) : i must be of type int or tuple(len=2)")
          return self[i[0]:i[1]]
        return self.data[i]
    def __getslice__(self, i, j):
        return Pdb(self, self.data[i:j])
    def __len__(self):
        return len(self.data)
    
    def __str__(self):
        if self.code:
          return "TITLE     %-70s\n%s%s\n" % (self.code, self.cryst1, join("\n", self.data))
        else:
          return "%s%s\n" % (self.cryst1, join("\n", self.data))
    
    def __repr__(self):
        #return "{Pdb:%s:%s}" % (self.code, self.data.__str__())
        return 'Pdb("%s", %s, cryst1=%s)' % (self.code, repr(self.data), repr(self.cryst1))
    
    def __cmp__(self, other):
        v = cmp(type(self), type(other))
        if v != 0:
          return v
        v = cmp(self.code, other.code)
        if v != 0:
          return v
        return cmp(self.data, other.data)
    
    #def __eq__(self, other):
        #if type(other) == type(None):
          #return False
        #return 0 == cmp(self.code, other.code)
    
    #def __ne__(self, other):
        #if type(other) == type(None):
          #return True
        #return 0 != cmp(self.code, other.code)
    
    def __hash__(self):
      return hash(self.code) ^ hash(self.data)
    
    
    # Instance methods
    #
    
    
    def get_ligand_string(self):
        return join("\n", self.ligands) + "\n"
    
    def parse_file(self, data):
        """parse_file(data) : sets this Pdb object to contain the given data.
        
        data can be any iterable object, containing objects whose string representation corresponds to a line in a PDB file."""
        if not data:
            return None
        
        atoms = []
        ligands = []
        chaincodes = set([])
        must_renumber = False
        for item in data:
          item = str(item)
          for line in item.splitlines():
            if line.startswith('CRYST1'):
              self.cryst1 = line.rstrip('\n\r')+"\n"
            elif line.startswith('ATOM  ') or line.startswith('HETATM'):
              try:
                atm = Atom(line)
              except NoAtomCoordinatesError:
                sys.stderr.write("WARNING: Pdb '%s' contains atoms without coordinates\n"%self.code)
                continue
              
              if atm.res == "MSE":
                #sys.stderr.write("WARNING: Pdb '%s' contains MSE. Rewriting as MET.\n"%self.code)
                atm.res = "MET"
                if atm.atom == "SE":
                  atm.atom = "S"
                  atm.element = "S"
              
              if atm.hetatm and atm.is_amino():
                if not chaincodes or atm.chain in chaincodes:
                  #sys.stderr.write("WARNING: Pdb '%s' contains amino acids labelled as HETATM. Relablling as ATOM.\n"%self.code)
                  atm.hetatm = False
              
              if not atm.hetatm:
                atoms.append(atm)
                chaincodes.update(atm.chain)
                must_renumber |= atm.iatom < 0
              elif atm.res != "HOH": # Ignore waters
                ligands.append(atm)
                must_renumber |= atm.iatom < 0
              
        if not atoms:
            sys.stderr.write("WARNING: Failed to parse ATM file for pdb '%s' (file: %s). No valid ATOM coodinates?\n" % (self.code, file))
        
        self.data = atoms
        self.ligands = ligands
        if must_renumber:
          self.renumber_atoms()
        self.unindex()
    
    
    def append_atoms(self, data):
        """append(data) : appends the data to this Pdb object.
        """
        if data:
          assert isinstance(data[0], Atom), "Can only append() Atom objects to Pdb"
          self.data.extend(data)
        self.unindex()
    
    
    def sort(self, key=None):
        """sort() : sorts the atoms in this Pdb, in place, by (chain id, atom id, altloc).
        
        As a last-resort comparison, all ATOM data is compared alpha-numerically in the order they appear in the PDB file format (left to right).
        
        If the PDB is formatted correctly, calling this function on a Pdb object should not change the object at all. NOTE: This messes up duplicate chains though!"""
        if not self.data:
          return
        self.data.sort(key=key)
        self.unindex()
    
    
    def unindex(self):
        self.hashindex=None
        self.chainlimits=None


    def build_hashindex(self):
        self.hashindex = {}
        for i,atom in enumerate(self):
            self.hashindex[atom] = i
     
    
    def build_chainindex(self):
        self.chainlimits = {}
        
        if not self:
          return
        
        prev = self[0]
        
        # Start first chain
        self.chainlimits[prev.chain] = [0]
        
        for i,atom in enumerate(self):
            if (prev.chain != atom.chain) or (prev > atom):
                # Terminate previous chain
                if prev.chain not in self.chainlimits:
                    self.chainlimits[prev.chain] = [i]
                else:
                    self.chainlimits[prev.chain].append(i)
                
                # Start new chain
                if atom.chain not in self.chainlimits:
                    self.chainlimits[atom.chain] = [i]
                else:
                    self.chainlimits[atom.chain].append(i)
                prev = atom
        
        # Terminate last chain
        assert prev.chain in self.chainlimits, "Error calculating chain limits for PDB %s: Last chain limit belongs to unknown chain!?" % (self.code)
        self.chainlimits[prev.chain].append(len(self))
        
        for k in self.chainlimits:
            self.chainlimits[k] = tuple(self.chainlimits[k])
    
    
    def get_chain_limits(self, c):
        c = c.strip()
        if not self.chainlimits:
          self.build_chainindex()
        if self.chainlimits and c in self.chainlimits:
            return self.chainlimits[c]
        else:
            return None
    
    
    def xchains(self):
        "Generates one chain Pdb object at a time, ignoring duplicate chains."
        for c in self.get_chain_codes():
          yield self.get_chain(c, False)
    
    def chains(self):
        return list(self.xchains())
    
    
    def get_chain(self, chains, include_duplicates=True):
        data = []
        for c in chains:
            limits = self.get_chain_limits(c)
            if limits:
                if include_duplicates:
                    for i in range(0, len(limits), 2):
                        data.extend(self.data[limits[i]:limits[i+1]])
                else:
                    data.extend(self.data[limits[0]:limits[1]])
        return Pdb(self, data)
    
    
    def chaincount(self, include_duplicates=True):
        if not self.chainlimits:
          self.build_chainindex()
        if not include_duplicates:
          return len(self.chainlimits)
        n = 0
        for k,v in list(self.chainlimits.items()):
          n+=len(v)/2
        return n


    def get_first_chain(self):
        if not len(self):
            return self
        c = self[0].chain
        return self.get_chain(c, False)
    
    
    def get_chain_codes(self):
        """get_chain_codes() : Returns a tuple of chain identifiers, corresponding to all peptide chains in this Pdb. If the Pdb contains duplicate chains, they will still appear only once in the returned tuple."""
        if not len(self):
            return ()
        
        if not self.chainlimits:
          self.build_chainindex()
        
        c = list(self.chainlimits.keys())
        c.sort()
        return tuple(c)
    
    
    def get_structure_lign(self):
        """get_structure_lign() : Returns a string containing the description line, in ALI format, (beginning with "structure:") corresponding to this Pdb."""
        startnum   = self[0].ires
        startchain = self[0].chain
        endnum     = self[-1].ires
        endchain   = self[-1].chain
        
        return "structure:%4s:%4d :%1s:%4d :%1s::::" % (self.code, startnum, startchain, endnum, endchain)
    
    
    def rescount(self):
        """rescount() : Returns the number of residues in this protein (including all chains)."""
        n=0
        for rb in self.xresidue_boundaries():
          n+=1
        return n
    
    
    def get_CA(self, allowAlternatives=False):
        """get_CA() : Returns a Pdb object containing all C-alpha atoms."""
        return self.get_atoms_by_type("CA", allowAlternatives)
    
    def get_backbone(self, allowAlternatives=False):
        """get_backbone() : Returns a Pdb object containing all backbone (N, CA, C, O) atoms."""
        return self.get_atoms_by_type(("N", "CA", "C", "O"), allowAlternatives)
    
    def get_backbone_and_CB(self, allowAlternatives=False):
        """get_backbone_and_CB() : Returns a Pdb object containing all backbone (N, CA, C, O) and CB atoms."""
        return self.get_atoms_by_type(("N", "CA", "C", "O", "CB"), allowAlternatives)
    
    def get_CB(self, allowAlternatives=False):
        """get_CB() : Returns a Pdb object containing all CB atoms and CA atoms for Glycine."""
        newatoms = []
        for start,end in self.xresidue_boundaries():
          selected = None
          for i in range(start,end):
            a = self[i]
            if a.atom in ("CA", "CB"):
              selected = a
              if (allowAlternatives or not a.altloc or a.altloc == 'A') and (a.atom == "CB" or a.res == "GLY"):
                break
          assert selected is not None
          newatoms.append(selected)
        return Pdb(self, newatoms)
    
    
    def get_atoms_by_type(self, types, allowAlternatives=False):
        """
        get_atoms_by_type(types, allowAlternatives=False) : Returns a Pdb object containing all atoms of the specified type(s).
        
        types can be a string or a list of strings. allowAlternatives determines whether more than one version of a single atom is allowed.
        """
        if isinstance(types, str):
          types = (types,)
        
        newatoms = []
        for a in self:
          if a.atom in types:
            if allowAlternatives or not a.altloc or a.altloc == 'A':
              newatoms.append(a)
        return Pdb(self, newatoms)
    
    
    def get_atoms(self, slice=None, atom_types=None):
        """get_atoms() : Returns a list of Atom objects, matching the given criteria."""

        if slice:
          start, end = slice
        else:
          start=0
          end=len(self)
        
        if not atom_types:
          return self.data[start:end]

        result=[]
        for i in range(start, end):
          a = self.data[i]
          if a.atom in atom_types:
            result.append(a)
        return result
    
    
    def index(self, atom):
        if not len(self):
          raise KeyError("Pdb.index() : Atom not found. Pdb is empty.")
        elif not self.hashindex:
          self.build_hashindex()
        
        try:
          return self.hashindex[atom]
        except KeyError as e:
          raise KeyError("Pdb.index() : Given atom not found in index: %s" % (str(atom)))
    
    
    def get_nth_residue_boundaries(self, n):
        """get_nth_residue() : Returns the boundaries of the Nth residue (where 0 is the first residue)."""
        if n == 0:
          return self.get_residue_boundaries(0)
        
        for i, pair in enumerate(self.xresidue_boundaries()):
            if i == n:
              return pair
        raise ArgumentError("Nth residue does not exist. N (%d) >= rescount (%d)." % (n, self.rescount()))
    
    
    def get_residues(self, indeces):
        """get_residues(indeces) : Returns a new Pdb object containing the specified residues (where 0 is the first residue)."""
        
        atoms = Pdb(self, [])
        n=0
        for i, res in enumerate(self.xresidue_boundaries()):
          if n >= len(indeces):
            break
          if indeces[n] > i:
            continue
          if indeces[n] < i:
            raise ArgumentError("Residue indeces not unique and sorted:" + str(indeces))
          atoms.append_atoms(self.data[res[0]:res[1]])
          n+=1
        
        if n < len(indeces):
          raise ArgumentError("Reached end of Pdb before extracting residue %d. Actual residue count: %d" % (indeces[n], self.rescount()))
        
        return atoms
    
    
    def is_residue_boundary(self, index):
        if index < 0:
          index = len(self)+index
        assert index <= len(self)
        if index == 0 or index == len(self):
          return True
        return (self[index].ires    != self[index-1].ires) \
            or (self[index].chain   != self[index-1].chain) \
            or (self[index].inscode != self[index-1].inscode)
    
    
    def get_residue_slice_boundaries(self, offset, length, beginHere=None):
        """get_residue_slice_boundaries(offset, length) : Returns the start and end indices of the specified residue slice with the specified length and starting offset residues from beginHere (which is an atom index)."""
        
        assert length>=0
        
        if beginHere == None:
          beginHere = 0
        
        if offset<0 and -offset>=length:
          return self.rget_residue_slice_boundaries(-(offset+length), length, beginHere)
        
        # Make sure we always include the entire residue, if beginHere is within a residue
        if not self.is_residue_boundary(beginHere):
          beginHere = self.get_residue_boundaries(beginHere)
          beginHere = beginHere[0]
        
        if length==0:
          return (beginHere, beginHere)
        
        start=len(self)
        end=len(self)
        l=0
        
        if offset < 0:
          for i, res in enumerate(self.xrresidue_boundaries(beginHere)):
            start=res[0]
            l+=1
            if i+1 == -offset:
              break
        
        for i, res in enumerate(self.xresidue_boundaries(beginHere)):
          if i>=offset:
            if i==offset:
              start=res[0]
            end=res[1]
            l+=1
            if l>=length:
              break

        return (start, end)
    
    
    def get_residue_slice(self, offset, length, beginHere=None):
        """get_residue_slice(offset, length [, beginHere]) : Returns a new Pdb object containing the specified residue slice with the specified length and starting offset residues from beginHere."""
        s, e = self.get_residue_slice_boundaries(offset, length, beginHere)
        return self[s:e]
    
    
    def delete_residue_slice(self, offset, length, beginHere=None):
        s, e = self.get_residue_slice_boundaries(offset, length, beginHere)
        del self.data[s:e]
        self.unindex()
    
    
    def __delitem__(self, i):
        del self.data[i]
        self.unindex()


    def __delslice__(self, i, j):
        del self.data[i:j]
        self.unindex()
    
    
    def rget_residue_slice_boundaries(self, offset, length, beginHere=None):
        """rget_residue_slice(offset, length) : Returns a new Pdb object containing the specified residue slice with the specified length and ending offset residues from beginHere."""
        
        assert length>=0
        
        if beginHere == None:
          beginHere = len(self)
        
        if offset<0 and -offset>=length:
          return self.get_residue_slice_boundaries(-(offset+length), length, beginHere)
        
        # Make sure we always include the entire residue, if beginHere is within a residue
        if not self.is_residue_boundary(beginHere):
          beginHere = self.get_residue_boundaries(beginHere)
          beginHere = beginHere[0]
        
        if length==0:
          return (beginHere, beginHere)
        
        start=0
        end=0
        l=0
        
        if offset < 0:
          for i, res in enumerate(self.xresidue_boundaries(beginHere)):
            end=res[1]
            l+=1
            if i+1 == -offset:
              break
        
        for i, res in enumerate(self.xrresidue_boundaries(beginHere)):
          if i>=offset:
            if i==offset:
              end=res[1]
            start=res[0]
            l+=1
            if l>=length:
              break
        
        return (start, end)
    
    
    def rget_residue_slice(self, offset, length, beginHere=None):
        """rget_residue_slice(offset, length [, beginHere]) : Returns a new Pdb object containing the specified residue slice with the specified length and ending offset residues from beginHere."""
        s, e = self.rget_residue_slice_boundaries(offset, length, beginHere)
        return self[s:e]
    
    
    def rdelete_residue_slice(self, offset, length, beginHere=None):
        s, e = self.rget_residue_slice_boundaries(offset, length, beginHere)
        del self.data[s:e]
        self.unindex()
    
    
    def get_residue_boundaries(self, atom):
        """get_residue_boundaries(atom) : atom can be an integer index or an Atom object. Returns a pair of indeces (usable in slice syntax), delimiting the enclosing residue."""
        if isinstance(atom, Atom):
            iatom = self.index(atom)
        elif isinstance(atom, int):
            iatom = atom
            atom = self[iatom]
        else:
            raise ValueError("get_residue : argument must be of type Atom or int. Given : %s : %s" % (type(atom), repr(atom)))
        
        istart = iatom
        iend = iatom + 1  # slice syntax. self[i:i+1] returns a Pdb object containing a single atom.
        
        # search backwards, for the beginning of the residue
        while istart > 0 and not self.is_residue_boundary(istart):
          istart -= 1
        
        # search forwards, for the end of the residue
        while iend < len(self) and not self.is_residue_boundary(iend):
          iend += 1
        
        return (istart, iend)
    
    
    def residue_boundaries(self):
        return list(self.xresidue_boundaries())
    
    def xresidue_boundaries(self, beginHere=None):
        if None==beginHere:
          beginHere=0
        assert beginHere>=0 and beginHere<=len(self), "index out of range: %d (len=%d)"%(beginHere,len(self))
        
        istart = beginHere
        iend   = beginHere+1  # slice syntax. self[i:i+1] returns a Pdb object containing a single atom.
        
        while istart < len(self):
          # search forward, for the end of the residue
          while iend < len(self) and not self.is_residue_boundary(iend):
            iend += 1
          yield (istart, iend)
          istart = iend
          iend += 1
    
    
    def xrresidue_boundaries(self, beginHere=None):
        if None==beginHere:
          beginHere=len(self)
        assert beginHere>=0 and beginHere<=len(self), "index out of range: %d (len=%d)"%(beginHere,len(self))
        
        istart = beginHere-1
        iend   = beginHere  # slice syntax. self[i:i+1] returns a Pdb object containing a single atom.
        
        while iend > 0:
          # search backward, for the start of the residue
          while istart > 0 and not self.is_residue_boundary(istart):
            istart -= 1
          yield (istart, iend)
          iend = istart
          istart -= 1
    
    
    def residues(self):
        return list(self.xresidues())
    
    
    def xresidues(self):
        for start, end in self.xresidue_boundaries():
          yield self[start:end]
    
    
    def xrresidues(self):
        for start, end in self.xrresidue_boundaries():
          yield self[start:end]
    
    
    def get_gap_list(self):
      """get_gap_list(self) : returns a list of gaps, as determined using residue numbers.
      
      Format: [[gapstart, gaplength, index], ...] , where gapstart is the first missing residue's residue number and index is the list index where missing residues should be inserted."""
      results=[]
      for i in range(1, len(self)):
        if self.is_residue_boundary(i):
          gaplen = self[i].ires - self[i-1].ires - 1
          if gaplen > 0:
            results.append([self[i-1].ires+1, gaplen, i])
      return results
    

    def split(self):
      gaplist = self.get_gap_list()
      if not gaplist:
        return [self[:]]
      result=[]
      prev=0
      for g in gaplist:
        i=g[2]
        result.append(self[prev:i])
        prev=i
      result.append(self[prev:])
      return result

    
    def renumber(self):
        "renumber() : Renumbers all atoms/residues such that their numbers correspond to their position in the Pdb's Atom list"
        if not len(self):
          return
        resnum=1
        atmnum=1
        prev=self.data[0]
        for start, end in self.xresidue_boundaries():
          if prev.chain != self.data[start].chain or prev > self.data[start]:
            resnum = 1
          for i in range(start, end):
            self.data[i].ires = resnum
            self.data[i].iatom = atmnum
            self.data[i].inscode = ""
            atmnum+=1
          resnum+=1
          prev = self.data[end-1]
    
    def tostring_renumbered(self):
        "Returns a string representing a renumber()-ed version of this Pdb object"
        if not len(self):
          return ""
        
        output = array('c', "TITLE     %-70s\n"%(self.code))
        tempatm = self[0].copy()
        
        resnum=1
        atmnum=1
        prev=self.data[0]
        for start, end in self.xresidue_boundaries():
          if prev.chain != self.data[start].chain or prev > self.data[start]:
            resnum = 1
          for i in range(start, end):
            tempatm.copy(self.data[i])
            tempatm.ires = resnum
            tempatm.iatom = atmnum
            tempatm.inscode = ""
            output.extend(str(tempatm))
            output.append("\n")
            atmnum+=1
          resnum+=1
          prev = self.data[end-1]
        
        return output.tostring()
    
    def write_renumbered(self, out):
        "Writes a renumber()-ed version of this Pdb object to the given file object"
        if not len(self):
          return
        
        tempatm = self[0].copy()
        
        resnum=1
        atmnum=1
        prev=self.data[0]
        for start, end in self.xresidue_boundaries():
          if prev.chain != self.data[start].chain or prev > self.data[start]:
            resnum = 1
          for i in range(start, end):
            tempatm.copy(self.data[i])
            tempatm.ires = resnum
            tempatm.iatom = atmnum
            tempatm.inscode = ""
            out.write(str(tempatm))
            out.write("\n")
            atmnum+=1
          resnum+=1
          prev = self.data[end-1]
    
    
    def renumber_atoms(self):
        "renumber_atoms() : Renumbers all atoms, starting from 1 upwards, in the order they appear in the Atom list."
        
        offset = 1
        
        if len(self) > 0:
          for i,a in enumerate(self):
            a.iatom = i+1
          offset = self[-1].iatom + 1
        
        for i,a in enumerate(self.ligands):
          a.iatom = i+offset
    
    
    def filter(self, negres=False, noca=True, inscode=False, altloc=True, duperes=True):
        if not self.data:
          return
        
        newatoms = []
        atmhashset = set([])
        
        for start,end in self.xresidue_boundaries():
          if negres:
            if self[start].ires < 1:
              continue
          
          if noca:
            ok=False
            for i in range(start,end):
              if self[i].atom == "CA":
                ok=True
                break
            if not ok:
              continue
          
          resatms = []
          for ix in range(start, end):
            atm = self[ix]
            
            atmhash = "%s %s %s"%(atm.ires, atm.atom, atm.chain)
            
            if not inscode:
              atmhash += " %s"%(atm.inscode)
            
            if not altloc:
              atmhash += " %s"%(atm.altloc)
            else:
              atm.altloc = ""
            
            if atmhash in atmhashset:
              continue
            
            atmhashset.add(atmhash)
            resatms.append(atm)
          
          if noca:
            ok=False
            for atm in resatms:
              if atm.atom == "CA":
                ok=True
                break
            if not ok:
              continue
          
          if altloc:
            # Only keep the first residue type
            # This avoids keeping bits of several alternative residue types,
            # in cases where altloc is used on a per-residue level.
            #
            for atm in resatms:
              if atm.res == resatms[0].res:
                newatoms.append(atm)
          else:
            newatoms.extend(resatms)
          
        
        if len(newatoms) != len(self.data):
          self.data = newatoms
          self.unindex()
          return True
        return False
    
    
    def deep_copy(self):
        p = Pdb("", [])
        p.parent = self.parent
        p.code = self.code
        p.data = [atm.copy() for atm in self.data]
        p.ligands = [atm.copy() for atm in self.ligands]
        return p
    
    
    def get_residue(self, atom):
        """get_residue(atom) : atom can be an integer index or an Atom object. Returns a Pdb object containing all atoms belonging to the enclosing residue."""
        istart, iend = self.get_residue_boundaries(atom)
        return self[istart:iend]
    
    
    def get_seq(self, gapped=False, firstres=None, lastres=None):
        """get_seq() : Returns a FASTA string representation of the sequence.
        
        If gapped==True, '-' is inserted for missing residues and '/' in between chains.
        One letter is assigned to each residue returned by Pdb.xresidues()."""
        assert firstres==None or firstres>0
        assert lastres==None or (firstres==None and lastres>0) or lastres>=firstres
        
        if len(self) == 0:
          return ""
        
        seq = array('c')
        
        prev = self[0]
        
        if gapped and firstres and prev.ires > firstres:
          seq.extend('-'*(prev.ires-firstres))
        
        seq.append(residueLetter(prev.res))
        
        iter = self.xresidue_boundaries()
        next(iter) # skip first residue
        for start, end in iter:
            curr = self[start]
            
            if gapped:
              if   curr.chain != prev.chain:  # start of new chain
                  seq.append('/')
              elif not (curr.ires  == prev.ires+1 or (curr.ires == prev.ires and curr.inscode != prev.inscode)): # gap in same chain
                  seq.extend('-'*(curr.ires-prev.ires-1))
            
            seq.append(residueLetter(curr.res))
            prev = curr
        
        if gapped and lastres and prev.ires < lastres:
          seq.extend('-'*(lastres-prev.ires))
        
        return seq.tostring()
    
    
    def get_root_object(self):
        obj = self
        while obj.parent:
          obj = obj.parent
          if id(obj) == id(self):
            raise RuntimeError("Pdb.get_root_object() : Circular parent/child relationship. Aborting.") # avoid infinite loops
        return obj
    
    
    def delete_atoms_by_type(self, atmtype):
        newatoms = []
        for atm in self:
          if atm.atom != atmtype:
            newatoms.append(atm)
        if len(newatoms) < len(self):
          self.data = newatoms
          self.unindex()
          return True
        return False


    def has_atom_type(self, atype):
        for atm in self:
          if atm.atom == atype:
            return True
        return False
   
    def has_residue_type(self, rtype):
        for atm in self:
          if atm.res == rtype:
            return True
        return False
    
    
    # geometry methods
    
    def rotate(self, rotation_matrix):
        for atm in self:
            atm.xyz = numpy.dot(atm.xyz, rotation_matrix)
    
    def translate(self, translation_vector):
        for atm in self:
            atm.xyz += translation_vector
    
    def centroid(self):
        centre=numpy.array([0.0, 0.0, 0.0])
        for atm in self:
          centre += atm.xyz
        centre /= len(self)
        return centre
