#! /usr/bin/env python3

import re
import sys
import numpy as np
from copy import deepcopy
import xml.etree.ElementTree as et



## reads in a qmcpack particle file
## places a single atom-centered gaussian around each atom

def indent(elem, level=0):
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if elem.text:
      t = re.search(r"(.+)\n",elem.text,flags=re.DOTALL)
      if t is not None:
        elem.text = t.group(1) + i
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i

# main gaussian class
class gaussian:

  def _D_to_B(self, A, D):
    B = np.dot(A,D)
    return B
  
  def _B_to_D(self, A, B):
    try:
      D = np.linalg.solve(A,B)
    except:
      return np.nan
    return D
    # raise Exception("_B_to_D: conversion not implemented.")

  def _K_to_C(self, A, D, K):
    C = reduce(np.dot, [D,A,D]) + K
    return C

  def _C_to_K(self, A, B, C):
    try:
      D = np.linalg.solve(A,B)
    except:
      return np.nan
    K = C - np.dot(D,B)
    return K
    #raise Exception("_C_to_K: conversion not implemented.")

  # arguments:
  #   name
  #   A should be a 3x3 np.array
  #   B, D should be a 3 element np.array
  #   C, K should just be a floating point number
  def __init__(self, name, A, B, C, D = None, K = None):
    self.A = np.ndarray(buffer=deepcopy(A),shape=(3,3),dtype=np.float64)
    self.D = None if D is None else np.array(deepcopy(D),dtype=np.float64)
    self.K = None if K is None else np.float64(deepcopy(K))
    self.B = None if B is None else np.array(deepcopy(B),dtype=np.float64)
    self.C = None if C is None else np.float64(deepcopy(C))
    if B is not None and C is not None:
      self.D = self._B_to_D(self.A, self.B)
      self.K = self._C_to_K(self.A, self.B, self.C)
    elif D is not None and K is not None: # prefer to use D, K
      self.B = self._D_to_B(self.A, self.D)
      self.C = self._K_to_C(self.A, self.D, self.K)
    # assign variables
    self.name = deepcopy(name)
  
  def __str__(self):
    A = self.A
    B = self.B
    C = self.C
    D = self.D
    K = self.K
    name = self.name
    if B is not None and C is not None:
      s = ("gaussian: name: {}\n"
          "  A: {} {} {}\n"
          "     {} {} {}\n"
          "     {} {} {}\n"
          "  B: {} {} {}\n"
          "  C: {}\n")
      val_array = [ela for row in A for ela in row] + [elb for elb in B] + [C]
      s = s.format(name, *val_array)
      return s
    elif D is not None and D is not np.nan and K is not None and K is not np.nan:
      s = ("gaussian: name: {}\n"
          "  A: {} {} {}\n"
          "     {} {} {}\n"
          "     {} {} {}\n"
          "  D: {} {} {}\n"
          "  K: {}\n")
      val_array = [ela for row in A for ela in row] + [elb for elb in D] + [K]
      s = s.format(name, *val_array) 
      return s
    else:
      return "incomplete gaussian definition: A: {}, B: {}, C: {}, D: {}, K: {}".format(
          A is not None, B is not None, C is not None, D is not None, K is not None)

  # build and return xml tag structure 
  def xml_tag(self, var_list = ["A", "B", "C"], opt=True):
    # serialize A, B
    function_tag = et.Element("function",{"id":self.name})
    if "A" in var_list:
      A_serial = []
      for i in range(3):
        for j in range(3):
          A_serial.append(self.A[i,j]) if j >= i else None
      A_tag = et.Element("var",{ "name": "A", "opt": "(opt_a)" if opt else "False" } )
      A_tag.text = " ".join( map(str, A_serial ) )
      function_tag.append(A_tag)
    if "B" in var_list:
      B_tag = et.Element("var",{ "name": "B", "opt": "(opt_b)" if opt else "False" } )
      B_tag.text = " ".join( map(str, list(self.B)) )
      function_tag.append(B_tag)
    if "C" in var_list:
      C_tag = et.Element("var",{ "name": "C", "opt": "(opt_c)" if opt else "False" } )
      C_tag.text = str(self.C)
      function_tag.append(C_tag)
    if "D" in var_list:
      D_tag = et.Element("var",{ "name": "D", "opt": "(opt_d)" if opt else "False" } )
      D_tag.text = " ".join( map(str, list(self.D)) )
      function_tag.append(D_tag)
    if "K" in var_list:
      K_tag = et.Element("var",{ "name": "K", "opt": "(opt_k)" if opt else "False" } )
      K_tag.text = str(self.K)
      function_tag.append(K_tag)
    return function_tag

  def __deepcopy__(self,memo):
    gaussian_copy = gaussian(self.name, self.A, self.B, self.C)
    return gaussian_copy

  # g*c or g1*g2
  def __mul__(self, other):
    # other is another gaussian
    try:
      A = self.A + other.A
      B = self.B + other.B
      C = self.C + other.C
      D = self._B_to_D(self.A, self.B)
      K = self._C_to_K(self.A, self.B, self.C)
      return gaussian("none",A,B,C)
    # other is a number
    except:
      A = self.A
      B = self.B
      C = self.C + np.log(other)
      return gaussian("none",A,B,C)
  
  # c*g
  def __rmul__(self, other):
    A = self.A
    B = self.B
    C = self.C + np.log(other)
    return gaussian("none",A,B,C)

  # g *= c or g1 *= g2
  def __imul__(self, other):
    # other is another gaussian
    try:
      self.A += other.A
      self.B += other.B
      self.C += other.C
      self.D = self._B_to_D(self.A, self.B)
      self.K = self._C_to_K(self.A, self.B, self.C)
    # other is a number
    except:
      self.C += np.log(other)
    return self

  # g / c or g1 / g2
  def __div__(self,other):
    try:
      A = self.A - other.A
      B = self.B - other.B
      C = self.C - other.C
      self.D = self._B_to_D(self.A, self.B)
      self.K = self._C_to_K(self.A, self.B, self.C)
      return gaussian("none",A,B,C)
    except:
      A = self.A
      B = self.B
      C = self.C - np.log(other)
      self.D = self._B_to_D(self.A, self.B)
      self.K = self._C_to_K(self.A, self.B, self.C)
      return gaussian("none",A,B,C)

  # c / g
  def __rdiv__(self,other):
    A = self.A
    B = self.B
    C = self.C - np.log(other)
    self.D = self._B_to_D(self.A, self.B)
    self.K = self._C_to_K(self.A, self.B, self.C)
    return gaussian("none",A,B,C)


  # g /= c or g1 /= g2
  def __idiv__(self,other):
    # other is another gaussian
    try:
      self.A -= other.A
      self.B -= other.B
      self.C -= other.C
      self.D = self._B_to_D(self.A, self.B)
      self.K = self._C_to_K(self.A, self.B, self.C)
    # other must be a number
    except:
      self.C -= np.log(other)
    return self


  def logval(self,r):
    logval = reduce(np.dot,[r,self.A,r]) - 2*np.dot(self.B,r) + self.C
    return logval

  def gradlog(self,r):
    gradlog = np.array(2*np.dot(self.A,r) - 2*self.B)
    return gradlog

  def val(self,r):
    return np.exp(self.logval(r))

  def gradval(self,r):
    val = self.val(r)
    gradlog = self.gradlog(r)
    gradval = gradlog * val
    return gradval

  def hessval(self,r):
    A = self.A
    x,y,z = r
    val = self.val(r)
    gradval = self.gradval(r)
    gradlog = self.gradlog(r)
    gh = (2*A + np.outer(gradlog,gradlog))*val
    return gh

  # cf functions correspond to the counting function 1/1+g
  def cf_val(self,r):
    cf_val = 1/(1+self.val(r))
    return cf_val

  def cf_gradval(self,r):
    val = self.val(r)
    cf_gradval = -(2*np.dot(self.A,r) - 2*self.B)*val/(1 + val)**2
    return cf_gradval

  def cf_hessval(self,r):
    val = self.val(r)
    gradval = self.gradval(r)
    cf_val = self.cf_val(r)
    cf_gradval = self.cf_gradval(r)
    gh = self.hessval(r)
    cf_hessval = np.outer(gradval, gradval)*cf_val**2 - gh*cf_val
    return cf_hessval

def atomic_coords(ptclfile):
  ptcltree = et.parse(ptclfile)
  ptclroot = ptcltree.getroot()
  posnode = ptclroot.find(".//particleset[@name='ion0']/attrib[@name='position']")
  posblock =  posnode.text
  pts = []
  re_num = "(-?\d\.\d+e[+-]\d\d)"
  re_pos = "{re_num}\s+{re_num}\s+{re_num}".format(re_num=re_num)
  for pos in re.finditer(re_pos, posblock):
    pts.append( np.array(pos.groups(), dtype=np.float64) )
  return pts

def gaussian_generator(A, B, C, prefix="g"):
  n = len(A)
  g = []
  if n == 0:
    return g
  for i, (A,B,C) in enumerate(zip(A,B,C)):
    name = "{}{:0>{width}}".format(prefix, i, width=int(np.log10(n)))
    g.append(gaussian(name,A=A,B=B,C=C))
  return g

def build_cartesian_gaussians(pts):
  alpha = 1
  A = [-alpha*np.eye(3,3) for pt in pts]
  B = [np.dot(Apt, pt) for Apt, pt in zip(A, pts)]
  C = [reduce(np.dot, [pt, Apt, pt]) for Apt, pt in zip(A, pts)]
  g = gaussian_generator(A, B, C)
  return g

# writes an ncjf to the specified outfile
def write_ncjf(outpath, g_list, F, ncjf_name, gref = None):
  ng = len(g_list)
  cjf_tag = et.Element("ncjf")
  jast_tag = et.Element("jastrow",{'type':"Counting",'name':ncjf_name, 'region':'normalized_gaussian'})
  cjf_tag.append(jast_tag)
  # build debug tag
  #debug_tag = et.Element("debug", {'name':'debug', 'seqlen':'5', 'period':'10000'})
  #jast_tag.append(debug_tag)
  # build F matrix tag
  F_tri = F[np.triu_indices(ng)]
  F_tag = et.Element("var", {'name':"F", "opt":"(opt_F)"})
  # pretty print F-matrix as an upper-triangle
  F_str = ""
  for i in range(ng):
    F_str += "\n    " + " "*15*i + "{:>15.5e}"*(ng-i) + " "
  F_tag.text = F_str.format(*F_tri)
  F_tag.text += "\n    "
  jast_tag.append(F_tag)
  # build region tag
  if gref is None:
    gref = g_list[0]
  region_tag = et.Element("region", {"reference_id":gref.name, "opt":"(opt_C)"})
  jast_tag.append(region_tag)
  # gaussian function tags
  for gaussian in g_list:
    gaussian_tag = gaussian.xml_tag(var_list=["A","B","C"], opt=(gref.name != gaussian.name))
    region_tag.append(gaussian_tag)
  # indent and write to file
  indent(cjf_tag)
  print "  writing cjf tags to file: {}".format(outpath)
  et.ElementTree(cjf_tag).write(outpath)

if __name__ == "__main__":
  if len(sys.argv) < 3:
    print "usage: ncjf_gen.py qmc.ptcl.xml cjf.xml"
    sys.exit(0)
  ptclfile = sys.argv[1]
  cjffile = sys.argv[2]
  acoords = atomic_coords(ptclfile)
  g_list = build_cartesian_gaussians(acoords)
  Fdim = len(g_list)
  write_ncjf(cjffile, g_list, np.zeros((Fdim, Fdim)), "atom-centered")


