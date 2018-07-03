######################################################################################
## This file is distributed under the University of Illinois/NCSA Open Source License.
## See LICENSE file in top directory for details.
##
## Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
##
## File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
##                    Thomas Applencourt, applencourt@anl.gov,  Argonne National Laboratory
##
## File created by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
#######################################################################################



def savetoqmcpack(cell,mf,title="Default",kpts=[]):
  import h5py, re, sys
  from collections import defaultdict
  from pyscf.pbc import gto, scf, df, dft
  from numpy import empty
  

  PBC=False
  Gamma=False
  UnRestricted=False
  Complex=False
  Python3=False
  Python2=False

  if sys.version_info >= (3, 0):
     sys.stdout.write("Using Python 3.x\n") 
     Python3=True
  else:
     sys.stdout.write("Using Python 2.x\n") 
     Python2=True
 

  val=str(mf)
  ComputeMode= re.split('[. ]',val)

  SizeMode=len(ComputeMode)
  for i in range(SizeMode):
     if ComputeMode[i] in ("UHF","KUHF","UKS"):
           UnRestricted=True
     if ComputeMode[i]=="pbc":
           PBC=True

  if PBC and len(kpts) == 0:
        #sys.exit("You need to specify explicit the list of K-point (including gamma)")
        Gamma=True



  IonName=dict([('H',1),  ('He',2),  ('Li',3),('Be',4),  ('B', 5),  ('C', 6),  ('N', 7),('O', 8),  ('F', 9),   ('Ne',10),   ('Na',11),('Mg',12),   ('Al',13),   ('Si',14),   ('P', 15),   ('S', 16),('Cl',17),   ('Ar',18),   ('K', 19),   ('Ca',20),   ('Sc',21),   ('Ti',22),   ('V', 23),   ('Cr',24),   ('Mn',25),   ('Fe',26),   ('Co',27),   ('Ni',28),   ('Cu',29),   ('Zn',30),   ('Ga',31),   ('Ge',32),   ('As',33),   ('Se',34),   ('Br',35),   ('Kr',36),   ('Rb',37),   ('Sr',38),   ('Y', 39),  ('Zr',40),   ('Nb',41),   ('Mo',42),   ('Tc',43),   ('Ru',44),   ('Rh',45),   ('Pd',46),   ('Ag',47),   ('Cd',48),   ('In',49),   ('Sn',50),   ('Sb',51),   ('Te',52),   ('I', 53),   ('Xe',54),   ('Cs',55),   ('Ba',56),   ('La',57),   ('Ce',58), ('Pr',59),   ('Nd',60),   ('Pm',61),   ('Sm',62),   ('Eu',63),   ('Gd',64),   ('Tb',65),   ('Dy',66),   ('Ho',67),  ('Er',68),   ('Tm',69),   ('Yb',70),   ('Lu',71),   ('Hf',72),   ('Ta',73),   ('W', 74),   ('Re',75),   ('Os',76),   ('Ir',77),   ('Pt',78),   ('Au',79),   ('Hg',80), ('Tl',81),   ('Pb',82),  ('Bi',83),   ('Po',84),   ('At',85),   ('Rn',86),   ('Fr',87),   ('Ra',88),   ('Ac',89),   ('Th',90),   ('Pa',91),   ('U', 92),   ('Np',93)]) 



  H5_qmcpack=h5py.File(title+'.h5','w')
  groupApp=H5_qmcpack.create_group("application")
  if Python3: 
     strList=['PySCF']
     asciiList = [n.encode("ascii", "ignore") for n in strList]
     groupApp.create_dataset('code', (1,),'S5', asciiList)
  else:
     CodeData  = groupApp.create_dataset("code",(1,),dtype="S5")
     CodeData[0:] = "PySCF"

  CodeVer  = groupApp.create_dataset("version",(3,),dtype="i4")
  CodeVer[0:] = 1
  CodeVer[1:] = 4
  CodeVer[2:] = 2

  GroupPBC=H5_qmcpack.create_group("PBC")
  GroupPBC.create_dataset("PBC",(1,),dtype="b1",data=PBC)
  natom=cell.natm

  dt = h5py.special_dtype(vlen=bytes)
  #Group Atoms
  groupAtom=H5_qmcpack.create_group("atoms")

  #Dataset Number Of Atoms
  groupAtom.create_dataset("number_of_atoms",(1,),dtype="i4",data=natom)

  #Dataset Number Of Species 
  #Species contains (Atom_Name, Atom_Number,Atom_Charge,Atom_Core)
  l_atoms = [ (cell.atom_symbol(x),IonName[cell.atom_symbol(x)],cell.atom_charge(x),cell.atom_nelec_core(x)) for x in  range(natom)  ] 


  d = defaultdict(list)
  for i,t in enumerate(l_atoms):
        d[t].append(i)


  idxSpeciestoAtoms = dict()
  uniq_atoms= dict()
  for i, (k,v) in enumerate(d.items()):
        idxSpeciestoAtoms[i] = v
        uniq_atoms[i] = k

  idxAtomstoSpecies = dict()
  for k, l_v in idxSpeciestoAtoms.items():
        for v in l_v:
            idxAtomstoSpecies[v] = k
 
  NbSpecies=len(idxSpeciestoAtoms.keys())

  groupAtom.create_dataset("number_of_species",(1,),dtype="i4",data=NbSpecies)

  #Dataset positions 
  MyPos=groupAtom.create_dataset("positions",(natom,3),dtype="f8")
  for x in range(natom): 
    MyPos[x:]=cell.atom_coord(x)

  #Group Atoms
  for x in range(NbSpecies):
    atmname=str(uniq_atoms[x][0])
    groupSpecies=groupAtom.create_group("species_"+str(x))
    groupSpecies.create_dataset("atomic_number",(1,),dtype="i4",data=uniq_atoms[x][1])
    mylen="S"+str(len(atmname))
    if Python3:
       strList=[atmname]
       asciiList = [n.encode("ascii", "ignore") for n in strList]
       groupSpecies.create_dataset('name', (1,),mylen, asciiList)
    else:
       AtmName=groupSpecies.create_dataset("name",(1,),dtype=mylen)
       AtmName[0:]=atmname

    groupSpecies.create_dataset("charge",(1,),dtype="f8",data=uniq_atoms[x][2])
    groupSpecies.create_dataset("core",(1,),dtype="f8",data=uniq_atoms[x][3])
  SpeciesID=groupAtom.create_dataset("species_ids",(natom,),dtype="i4")

  for x in range(natom):
      SpeciesID[x:]  = idxAtomstoSpecies[x]



  #Parameter Group
  GroupParameter=H5_qmcpack.create_group("parameters")
  GroupParameter.create_dataset("ECP",(1,),dtype="b1",data=bool(cell.has_ecp()))
  bohrUnit=True
  Spin=cell.spin 

  GroupParameter.create_dataset("Unit",(1,),dtype="b1",data=bohrUnit) 
  GroupParameter.create_dataset("NbAlpha",(1,),dtype="i4",data=cell.nelec[0]) 
  GroupParameter.create_dataset("NbBeta",(1,),dtype="i4",data=cell.nelec[1]) 
  GroupParameter.create_dataset("NbTotElec",(1,),dtype="i4",data=cell.nelec[0]+cell.nelec[1])
  GroupParameter.create_dataset("spin",(1,),dtype="i4",data=Spin) 
   

  #basisset Group
  GroupBasisSet=H5_qmcpack.create_group("basisset")
  #Dataset Number Of Atoms
  GroupBasisSet.create_dataset("NbElements",(1,),dtype="i4",data=NbSpecies)
  if Python3:
     strList=['LCAOBSet']
     asciiList = [n.encode("ascii", "ignore") for n in strList]
     GroupBasisSet.create_dataset('name', (1,),'S8', asciiList)
  else:
     LCAOName=GroupBasisSet.create_dataset("name",(1,),dtype="S8")
     LCAOName[0:]="LCAOBSet"

  #atomicBasisSets Group
  for x in range(NbSpecies):

    MyIdx=idxAtomstoSpecies[x]
    atomicBasisSetGroup=GroupBasisSet.create_group("atomicBasisSet"+str(x))
    mylen="S"+str(len(uniq_atoms[x][0]))

    if Python3:
       strList=[uniq_atoms[x][0]]
       asciiList = [n.encode("ascii", "ignore") for n in strList]
       atomicBasisSetGroup.create_dataset('elementType', (1,),mylen, asciiList)
       if cell.cart==True:
          strList=['cartesian']
          asciiList = [n.encode("ascii", "ignore") for n in strList]
          atomicBasisSetGroup.create_dataset('angular', (1,),'S9', asciiList)

          strList=['Gamess']
          asciiList = [n.encode("ascii", "ignore") for n in strList]
          atomicBasisSetGroup.create_dataset('expandYlm', (1,),'S6', asciiList)

       else:
          strList=['spherical']
          asciiList = [n.encode("ascii", "ignore") for n in strList]
          atomicBasisSetGroup.create_dataset('angular', (1,),'S9', asciiList)

          strList=['pyscf']
          asciiList = [n.encode("ascii", "ignore") for n in strList]
          atomicBasisSetGroup.create_dataset('expandYlm', (1,),'S5', asciiList)
    else:
       elemtype=atomicBasisSetGroup.create_dataset("elementType",(1,),dtype=mylen)
       elemtype[0:]=uniq_atoms[x][0]
       if cell.cart==True:
          Angular=atomicBasisSetGroup.create_dataset("angular",(1,),dtype="S9")
          ExpandYLM=atomicBasisSetGroup.create_dataset("expandYlm",(1,),dtype="S6")
          Angular[0:]="cartesian"
          ExpandYLM[0:]="Gamess"
       else:
          Angular=atomicBasisSetGroup.create_dataset("angular",(1,),dtype="S9")
          Angular[0:]="spherical"
          ExpandYLM=atomicBasisSetGroup.create_dataset("expandYlm",(1,),dtype="S5")
          ExpandYLM[0:]="pyscf"



    atomicBasisSetGroup.create_dataset("grid_npts",(1,),dtype="i4",data=1001)
    atomicBasisSetGroup.create_dataset("grid_rf",(1,),dtype="i4",data=100)
    atomicBasisSetGroup.create_dataset("grid_ri",(1,),dtype="f8",data=1e-06)

    mylen="S"+str(len(cell.basis))
    if Python3:
      strList=['log']
      asciiList = [n.encode("ascii", "ignore") for n in strList]
      atomicBasisSetGroup.create_dataset('grid_type', (1,),'S3', asciiList)
      if (len(cell.basis)<=2):
         strList=['gaussian']
         asciiList = [n.encode("ascii", "ignore") for n in strList]
         atomicBasisSetGroup.create_dataset('name', (1,),'S8', asciiList)
      else:
         strList=[cell.basis]
         asciiList = [n.encode("ascii", "ignore") for n in strList]
         atomicBasisSetGroup.create_dataset('name', (1,),mylen, asciiList)
      strList=['no']
      asciiList = [n.encode("ascii", "ignore") for n in strList]
      atomicBasisSetGroup.create_dataset('normalized', (1,),'S2', asciiList)
    else:
      gridType=atomicBasisSetGroup.create_dataset("grid_type",(1,),dtype="S3")
      gridType[0:]="log"
      if (len(cell.basis)<=2):
        nameBase=atomicBasisSetGroup.create_dataset("name",(1,),dtype="S8")
        nameBase[0:]="gaussian"
      else:
        nameBase=atomicBasisSetGroup.create_dataset("name",(1,),dtype=mylen)
        nameBase[0:]=cell.basis

      Normalized=atomicBasisSetGroup.create_dataset("normalized",(1,),dtype="S2")
      Normalized[0:]="no"


       


 

    nshell = cell.atom_shell_ids(MyIdx)
    n=0
    for i in nshell:
        l = cell.bas_angular(i)   
        contracted_coeffs = cell.bas_ctr_coeff(i)
        contracted_exp =cell.bas_exp(i)
        for line in zip(*contracted_coeffs):
          BasisGroup=atomicBasisSetGroup.create_group("basisGroup"+str(n))


          mylen="S"+str(len((uniq_atoms[x][0]+str(n)+str(l))))
          if Python3:
            strList=['Gaussian'] 
            asciiList = [n.encode("ascii", "ignore") for n in strList]
            BasisGroup.create_dataset('type',(1,),'S8',asciiList)


            strList=[uniq_atoms[x][0]+str(n)+str(l)] 
            asciiList = [n.encode("ascii", "ignore") for n in strList]
            BasisGroup.create_dataset('rid', (1,),mylen, asciiList)
          else:
            basisType=BasisGroup.create_dataset("type",(1,),dtype="S8")
            basisType[0:]="Gaussian"
            RID=BasisGroup.create_dataset("rid",(1,),dtype=mylen)
            RID[0:]=(uniq_atoms[x][0]+str(n)+str(l))


          BasisGroup.create_dataset("Shell_coord",(3,),dtype="f8",data=cell.bas_coord(i))
          BasisGroup.create_dataset("NbRadFunc",(1,),dtype="i4",data=cell.bas_nprim(i))
          Val_l=BasisGroup.create_dataset("l",(1,),dtype="i4",data=l)
          Val_n=BasisGroup.create_dataset("n",(1,),dtype="i4",data=n)
          RadGroup=BasisGroup.create_group("radfunctions")
          #print "<basisGroup",n," rid=",uniq_atoms[x][0]+str(n)+str(l)," n=",n,"  l=",l ,"NbRadFunc=",cell.bas_nprim(i),"type=Gaussian>"
          IdRad=0

          for e,c in zip(contracted_exp,line):
              DataRadGrp=RadGroup.create_group("DataRad"+str(IdRad))
              DataRadGrp.create_dataset("exponent",(1,),dtype="f8",data=e)
              DataRadGrp.create_dataset("contraction",(1,),dtype="f8",data=c)
              #print  "<radfunc exponent=",e," contraction=",c, "DataRad=",n,"IdRad=",IdRad,"/>"
              IdRad+=1
          n+=1

    atomicBasisSetGroup.create_dataset("NbBasisGroups",(1,),dtype="i4",data=n)

  def is_complex(l):
      try:
              return is_complex(l[0])
      except:
              return bool(l.imag)

    



  if cell.cart==True:
    # Generated from read_order.py in Numerics/codegen
    d_gms_order = {
        0:[""],
        1:["x","y","z"],
        2:["xx","yy","zz","xy","xz","yz"],
        3:["xxx","yyy","zzz","xxy","xxz","yyx","yyz","zzx","zzy","xyz"],
        4:["xxxx","yyyy","zzzz","xxxy","xxxz","yyyx","yyyz","zzzx","zzzy","xxyy","xxzz","yyzz","xxyz","yyxz","zzxy"],
        5:["xxxxx","yyyyy","zzzzz","xxxxy","xxxxz","yyyyx","yyyyz","zzzzx","zzzzy","xxxyy","xxxzz","yyyxx","yyyzz","zzzxx","zzzyy","xxxyz","yyyxz","zzzxy","xxyyz","xxzzy","yyzzx"],
        6:["xxxxxx","yyyyyy","zzzzzz","xxxxxy","xxxxxz","yyyyyx","yyyyyz","zzzzzx","zzzzzy","xxxxyy","xxxxzz","yyyyxx","yyyyzz","zzzzxx","zzzzyy","xxxxyz","yyyyxz","zzzzxy","xxxyyy","xxxzzz","yyyzzz","xxxyyz","xxxzzy","yyyxxz","yyyzzx","zzzxxy","zzzyyx","xxyyzz"],
    }
    
    d_l = {'s':0, 'p':1, 'd':2, 'f':3, 'g':4, 'h':5, 'i':6}
  
    def n_orbital(n):
      if n==0:
          return 1
      elif n==1:
          return 3
      else:
          return 2*n_orbital(n-1)-n_orbital(n-2)+1
  
  
    def compare_gamess_style(item1, item2):
      # Warning:
      # 	- d_gms_order is a global variable
      n1,n2 = map(len,(item1,item2))
      assert (n1 == n2)
      try:
          l = d_gms_order[n1]
      except KeyError:
          return 0
      else:
          a = l.index(item1)
          b = l.index(item2)
          return ((a>b) - (a<b)) #cmp( a, b )
 
    def compare_python3(item1, item2):
         return compare_gamess_style(item1[0],item2[0])
 
    ao_label = cell.ao_labels(False)
  
  
    # Create a list of shell
    l_l = []
    for label, name, t, l in ao_label:
          # Change yyx -> xyy "
          q =  "".join(sorted(l, key=l.count, reverse=True))
          l_l.append(q)
  
    # Pyscf ordering of shell
    l_order = list(range(len(l_l)))
  
    # Shell ordering indexed
    n = 1
    l_order_new = []
    for  i,(label, name, t, l) in enumerate(ao_label):
          r = d_l[t[-1]]
          # print r,n_orbital(r)
          if n != 1:
                  n-=1
          else:
                  from functools import cmp_to_key
                  n = n_orbital(r)
                  unordered_l = l_l[i:i+n]
                  unordered = l_order[i:i+n]
                  #print i,n,unordered 
                  ordered = [x for _,x in sorted(zip(unordered_l,unordered),key=cmp_to_key(compare_python3))]
                  l_order_new.extend(ordered)
  
    
    def order_mo_coef(ll):
          # Order a list of transposed mo_coeff (Ao,Mo) -> (Mo,Ao) ordered
          # Warning:
          #	- l_order_new is used as global variable
          #	- gamess order
  
          ll_new= []
          for l in zip(*ll):
              ll_new.append([l[i] for i in l_order_new])
          return ll_new
  
  mo_coeff = mf.mo_coeff
  if len(kpts)==0:
     Complex=False
  else:
     Complex=True
  GroupParameter.create_dataset("IsComplex",(1,),dtype="b1",data=Complex)

 
  GroupParameter.create_dataset("SpinUnResticted",(1,),dtype="b1",data=UnRestricted)
  GroupNbkpts=H5_qmcpack.create_group("Nb_KPTS")
  if not PBC:
    Nbkpts=1
    GroupNbkpts.create_dataset("Nbkpts",(1,),dtype="i4",data=Nbkpts)
    
    GroupDet=H5_qmcpack.create_group("KPTS_0")
    if UnRestricted==False:
      NbMO=len(mo_coeff)
      NbAO=len(mo_coeff[0])
      if cell.cart==True:
        eigenset=GroupDet.create_dataset("eigenset_0",(NbMO,NbAO),dtype="f8",data=order_mo_coef(mo_coeff))
      else:
        eigenset=GroupDet.create_dataset("eigenset_0",(NbMO,NbAO),dtype="f8",data=list(zip(*mo_coeff)))
    else:
      NbMO=len(mo_coeff[0])
      NbAO=len(mo_coeff[0][0])
      eigenset_up=GroupDet.create_dataset("eigenset_0",(NbMO,NbAO),dtype="f8",data=order_mo_coef(mo_coeff[0]))
      eigenset_dn=GroupDet.create_dataset("eigenset_1",(NbMO,NbAO),dtype="f8",data=order_mo_coef(mo_coeff[1]))
  else:
    #Cell Parameters
    GroupCell=H5_qmcpack.create_group("Cell")
    GroupCell.create_dataset("LatticeVectors",(3,3),dtype="f8",data=cell.lattice_vectors())


    if Gamma:  
       if not UnRestricted:
         NbMO=len(mo_coeff)
         NbAO=len(mo_coeff[0])
       else:
         NbMO=len(mo_coeff[0])
         NbAO=len(mo_coeff[0][0])
    else:
       if not UnRestricted:
         NbMO=len(mo_coeff[0])
         NbAO=len(mo_coeff[0][0])
       else:
         NbMO=len(mo_coeff[0][0])
         NbAO=len(mo_coeff[0][0][0])



    def get_mo(mo_coeff, cart):
        return order_mo_coef(mo_coeff) if cart else zip(*mo_coeff)
    if Gamma:
      Nbkpts=1
    else:
      Nbkpts=len(kpts)
    GroupNbkpts.create_dataset("Nbkpts",(1,),dtype="i4",data=Nbkpts)

    for i in range(Nbkpts):
      GroupDet=H5_qmcpack.create_group("KPTS_"+str(i))
      if Gamma:
         GroupDet.create_dataset("Coord",(1,3),dtype="f8",data=[0.0,0.0,0.0])
         if not UnRestricted:
           mo_coeff_ = get_mo(mo_coeff, cell.cart) 
            
           eigenset=GroupDet.create_dataset("eigenset_0",(NbMO,NbAO),dtype="f8",data=mo_coeff_) 
           eigenvalue=GroupDet.create_dataset("eigenval_0",(1,NbMO),dtype="f8",data=mf.mo_energy)
 
         else:
 
           mo_coeff_up = get_mo(mo_coeff[0], cell.cart) 
           mo_coeff_down = get_mo(mo_coeff[1], cell.cart)
 
           GroupDet.create_dataset("eigenset_0",(NbMO,NbAO),dtype="f8",data=mo_coeff_up)
           GroupDet.create_dataset("eigenset_1",(NbMO,NbAO),dtype="f8",data=mo_coeff_down)
 
           GroupDet.create_dataset("eigenval_0",(1,NbMO),dtype="f8",data=mf.mo_energy[0])
           GroupDet.create_dataset("eigenval_1",(1,NbMO),dtype="f8",data=mf.mo_energy[1])
      else:
         GroupDet.create_dataset("Coord",(1,3),dtype="f8",data=kpts[i])
         if not UnRestricted:
           mo_coeff_real = get_mo(mo_coeff[i].real, cell.cart) 
           mo_coeff_imag = get_mo(mo_coeff[i].imag, cell.cart) 
    #       moc_pack = empty((NbMO,NbAO,2),dtype=float)
    #       moc_pack[:,:,0] = mo_coeff_real
    #       moc_pack[:,:,1] = mo_coeff_imag
            

    #       GroupDet.create_dataset("eigenset_0",(NbMO,NbAO,2),dtype="f8",data=moc_pack) 

           GroupDet.create_dataset("eigenset_0_real",(NbMO,NbAO),dtype="f8",data=mo_coeff_real) 
           GroupDet.create_dataset("eigenset_0_imag",(NbMO,NbAO),dtype="f8",data=mo_coeff_imag) 
           GroupDet.create_dataset("eigenval_0",(1,NbMO),dtype="f8",data=mf.mo_energy[i])
 
         else:
 
           mo_coeff_up_real = get_mo(mo_coeff[0][i].real, cell.cart) 
           mo_coeff_up_imag = get_mo(mo_coeff[0][i].imag, cell.cart) 
           mo_coeff_down_real = get_mo(mo_coeff[1][i].real, cell.cart)
           mo_coeff_down_imag = get_mo(mo_coeff[1][i].imag, cell.cart)
 
           GroupDet.create_dataset("eigenset_0_real",(NbMO,NbAO),dtype="f8",data=mo_coeff_up_real)
           GroupDet.create_dataset("eigenset_0_imag",(NbMO,NbAO),dtype="f8",data=mo_coeff_up_imag)
           GroupDet.create_dataset("eigenset_1_real",(NbMO,NbAO),dtype="f8",data=mo_coeff_down_real)
           GroupDet.create_dataset("eigenset_1_imag",(NbMO,NbAO),dtype="f8",data=mo_coeff_down_imag)
 
           GroupDet.create_dataset("eigenval_0",(1,NbMO),dtype="f8",data=mf.mo_energy[0][i])
           GroupDet.create_dataset("eigenval_1",(1,NbMO),dtype="f8",data=mf.mo_energy[1][i])


  GroupParameter.create_dataset("numMO",(1,),dtype="i4",data=NbMO)
  GroupParameter.create_dataset("numAO",(1,),dtype="i4",data=NbAO)
  
  H5_qmcpack.close()

  print ('Wavefunction successfuly saved to QMCPACK HDF5 Format')
  print ('Use: "convert4qmc -pyscf  {}.h5" to generate QMCPACK input files'.format(title))
  # Close the file before exiting


