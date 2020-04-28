import scalop.prosci.loops.ccd as ccd
import numpy
from scalop.prosci.util.residue import is_residue_consecutive as is_consecutive
def meld(fragScaffold, fragPrediction, invertWeights=False):
    """Averages the coordinates of the two Pdb arguments. Move the first object's co-ordinates onto the averaged position.
    
    By default, the first object is assumed to be part of the N-terminal fragment, the second is part of the C-terminal fragment. This can be reversed by setting invertWeights=True."""
    
    resS = fragScaffold
    resP = fragPrediction
    L = len(resS)
    
    assert len(resS) == len(resP)
    
    def averageCoord(P, S):
      # P = coordinate of prediction
      # S = coordinate of scaffold
      # D = distance (in residues) from (loop+anchor)-fragment end
      # L = anchor length
      return 1.0/(L+1) * (D*P + (L+1-D)*S)
    
    for i, (rS, rP) in enumerate(zip(resS, resP)):
      if invertWeights:
        D = len(resS)-i
      else:
        D = i+1
      
      newN = averageCoord(rP.N.xyz, rS.N.xyz)
      newCA = averageCoord(rP.CA.xyz, rS.CA.xyz)
      newC = averageCoord(rP.C.xyz, rS.C.xyz)
      
      T_from, T_to, rotmat = ccd.get_rotmat([rS.N.xyz, rS.CA.xyz, rS.C.xyz], [newN, newCA, newC])
      
      rS.O = None # Remove main chain oxygen - can regenerate this later
      
      for a in rS:
        a.xyz = numpy.dot(a.xyz - T_from, rotmat) + T_to

def add_oxygens(residues, start=0, end=None, dO=1.23, force=False):
  if end is None:
    end = len(residues) - 1
  for i in range(start, end):
    # add missing mainchain oxygens, making them planar
    # with CA, C and the next residue's N
    #
    r = residues[i]
    if r is None:
      continue
    
    if r.O is None:
      q = residues[i+1]
      if q is None or q.N is None or r.C is None:
        continue
      
      if not force and not is_consecutive(r, q):
        continue
      
      O = r.C.copy()
      O.atom = "O"
      if r.C.element:
        O.element = "O"
      O.iatom = r.C.iatom+1

      vO_1 = r.C.xyz - r.CA.xyz
      vO_2 = r.C.xyz - q.N.xyz
      vO_3 = (vO_1 / numpy.linalg.norm(vO_1)) + (vO_2 / numpy.linalg.norm(vO_2))
      O.xyz = r.C.xyz + (vO_3 / numpy.linalg.norm(vO_3)) * dO

      r.O = O
      