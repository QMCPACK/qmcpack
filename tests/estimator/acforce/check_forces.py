#! /usr/bin/env python3

import sys
import numpy as np
import math


reference_key = { "Index" : 0,
  "LocalEnergy" : 1,
  "LocalPotential" : 3,
  "Kinetic" : 4,
  "NonLocalECP" : 8,
  "ACForce_hf_0_0" : 9,       
  "ACForce_pulay_0_0" : 10,    
  "ACForce_Ewfgrad_0_0" : 11,
  "ACForce_wfgrad_0_0" : 12, 
  "ACForce_hf_0_1" : 13,       
  "ACForce_pulay_0_1" : 14,    
  "ACForce_Ewfgrad_0_1" : 15,
  "ACForce_wfgrad_0_1" : 16, 
  "ACForce_hf_0_2" : 17,       
  "ACForce_pulay_0_2" : 18,    
  "ACForce_Ewfgrad_0_2" : 19,
  "ACForce_wfgrad_0_2" : 20, 
  "ACForce_hf_1_0" : 21,       
  "ACForce_pulay_1_0" : 22,    
  "ACForce_Ewfgrad_1_0" : 23,
  "ACForce_wfgrad_1_0" : 24, 
  "ACForce_hf_1_1" : 25,       
  "ACForce_pulay_1_1" : 26,    
  "ACForce_Ewfgrad_1_1" : 27,
  "ACForce_wfgrad_1_1" : 28, 
  "ACForce_hf_1_2" : 29,       
  "ACForce_pulay_1_2" : 30,    
  "ACForce_Ewfgrad_1_2" : 31,
  "ACForce_wfgrad_1_2" : 32 } 

reference_vals = { "LocalEnergy" : -1.5586292800e+01,
  "LocalPotential" : -2.8321193825e+01,
  "Kinetic" : 1.2734901025e+01,
  "NonLocalECP" : 1.8497959759e+00,
  "ACForce_0_0"  : 4.9434111234e-01,
  "ACForce_Ewfgrad_0_0" : -4.2095283359e+00,
  "ACForce_wfgrad_0_0" : 3.0004600095e-01,
  "ACForce_0_1"  : 3.5819786542e-01,
  "ACForce_Ewfgrad_0_1" : 2.3852191199e+00,
  "ACForce_wfgrad_0_1" : -1.8114010291e-01,
  "ACForce_0_2"  : 5.2292922210e+00,
  "ACForce_Ewfgrad_0_2" : 2.4541551533e+01,
  "ACForce_wfgrad_0_2" : -1.5614748618e+00,
  "ACForce_1_0"  : 1.1896691601e+00,
  "ACForce_Ewfgrad_1_0" : -8.1704406106e+00,
  "ACForce_wfgrad_1_0" : 3.3676305176e-01,
  "ACForce_1_1"  : 4.7617264236e+00,
  "ACForce_Ewfgrad_1_1" : -1.7346902278e+01,
  "ACForce_wfgrad_1_1" : 8.3752812817e-01,
  "ACForce_1_2"  : -4.1298777683e-02,
  "ACForce_Ewfgrad_1_2" : -4.6648162310e+01,
  "ACForce_wfgrad_1_2" : 2.6949237554e+00,
}  

#zero indexed.  
def grab_data_line(fname,stepnum):
  f=open(fname,'r')
  f.readline() #get rid of the header.
  for line in f:
    sl=line.split()
    if int(sl[0])==stepnum:
      myarray=list(map(float,sl))
      return np.array(myarray)
  return -1
    
if __name__ == "__main__":
  #CN molecule, so Natoms = 2
  natom=2
  #3D system, so.
  ndim=3

  relative_tol=1e-5

  result=grab_data_line("vmc.s000.scalar.dat",4)
  all_pass=True
  if not math.isclose(reference_vals["LocalEnergy"],result[reference_key["LocalEnergy"]],rel_tol=relative_tol):
    print("Error.  LocalEnergy Ref = ",reference_vals["LocalEnergy"], " Val = ",result[reference_key["LocalEnergy"]])
    all_pass=False
  if not math.isclose(reference_vals["Kinetic"],result[reference_key["Kinetic"]],rel_tol=relative_tol):
    print("Error.  Kinetic Ref = ",reference_vals["Kinetic"], " Val = ",result[reference_key["Kinetic"]])
    all_pass=False
  if not math.isclose(reference_vals["NonLocalECP"],result[reference_key["NonLocalECP"]],rel_tol=relative_tol):
    print("Error.  NonLocalECP Ref = ",reference_vals["NonLocalECP"], " Val = ",result[reference_key["NonLocalECP"]])
    all_pass=False

  for iat in range(0,natom):
    for idim in range(0,ndim):
      totforce = result[reference_key["ACForce_hf_%d_%d"%(iat,idim)]] + result[reference_key["ACForce_pulay_%d_%d"%(iat,idim)]]
      if not math.isclose(reference_vals["ACForce_%d_%d"%(iat,idim)], totforce, rel_tol=relative_tol):
        all_pass=False
        ref=reference_vals["ACForce_%d_%d"%(iat,idim)]
        val=totforce
        print("Error. ACForce_%d_%d Ref = "%(iat,idim),ref," Val = ",val)
      if not math.isclose(reference_vals["ACForce_Ewfgrad_%d_%d"%(iat,idim)],result[reference_key["ACForce_Ewfgrad_%d_%d"%(iat,idim)]], rel_tol=relative_tol):
        all_pass=False
        ref=reference_vals["ACForce_Ewfgrad_%d_%d"%(iat,idim)]
        val=result[reference_key["ACForce_Ewfgrad_%d_%d"%(iat,idim)]]
        print("Error. ACForce_Ewfgrad_%d_%d Ref = "%(iat,idim),ref," Val = ",val)
      if not math.isclose(reference_vals["ACForce_wfgrad_%d_%d"%(iat,idim)],result[reference_key["ACForce_wfgrad_%d_%d"%(iat,idim)]], rel_tol=relative_tol):
        all_pass=False
        ref=reference_vals["ACForce_wfgrad_%d_%d"%(iat,idim)]
        val=result[reference_key["ACForce_wfgrad_%d_%d"%(iat,idim)]]
        print("Error. ACForce_wfgrad_%d_%d Ref = "%(iat,idim),ref," Val = ",val)

        
  if(all_pass):
    exit(0)
  else:
    exit(1)
