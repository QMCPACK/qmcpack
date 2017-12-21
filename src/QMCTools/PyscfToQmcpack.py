
def savetoqmcpack(cell,mf,Title):
  import h5py, re
  from collections import defaultdict


  PBC=False
  RHF=False
  UHF=False
  ROHF=False

  val=str(mf)
  ComputeMode= re.split('[. ]',val)

  SizeMode=len(ComputeMode)
  print 'MySize=',SizeMode
  for i in range(SizeMode):
     if ComputeMode[i]=="RHF":
           RHF=True
     if ComputeMode[i]=="ROHF":
           ROHF=True
     if ComputeMode[i]=="UHF":
           UHF=True
     if ComputeMode[i]=="pbc":
           PBC=True
     
  print 'PBC=',PBC,' RHF=',RHF,' UHF=',UHF,' ROHF=',ROHF


  IonName=dict([('H',1),  ('He',2),  ('Li',3),('Be',4),  ('B', 5),  ('C', 6),  ('N', 7),('O', 8),  ('F', 9),   ('Ne',10),   ('Na',11),('Mg',12),   ('Al',13),   ('Si',14),   ('P', 15),   ('S', 16),('Cl',17),   ('Ar',18),   ('K', 19),   ('Ca',20),   ('Sc',21),   ('Ti',22),   ('V', 23),   ('Cr',24),   ('Mn',25),   ('Fe',26),   ('Co',27),   ('Ni',28),   ('Cu',29),   ('Zn',30),   ('Ga',31),   ('Ge',32),   ('As',33),   ('Se',34),   ('Br',35),   ('Kr',36),   ('Rb',37),   ('Sr',38),   ('Y', 39),  ('Zr',40),   ('Nb',41),   ('Mo',42),   ('Tc',43),   ('Ru',44),   ('Rh',45),   ('Pd',46),   ('Ag',47),   ('Cd',48),   ('In',49),   ('Sn',50),   ('Sb',51),   ('Te',52),   ('I', 53),   ('Xe',54),   ('Cs',55),   ('Ba',56),   ('La',57),   ('Ce',58), ('Pr',59),   ('Nd',60),   ('Pm',61),   ('Sm',62),   ('Eu',63),   ('Gd',64),   ('Tb',65),   ('Dy',66),   ('Ho',67),  ('Er',68),   ('Tm',69),   ('Yb',70),   ('Lu',71),   ('Hf',72),   ('Ta',73),   ('W', 74),   ('Re',75),   ('Os',76),   ('Ir',77),   ('Pt',78),   ('Au',79),   ('Hg',80), ('Tl',81),   ('Pb',82),  ('Bi',83),   ('Po',84),   ('At',85),   ('Rn',86),   ('Fr',87),   ('Ra',88),   ('Ac',89),   ('Th',90),   ('Pa',91),   ('U', 92),   ('Np',93)]) 


  H5_qmcpack=h5py.File(Title+'.h5','w')
  groupApp=H5_qmcpack.create_group("application")
  CodeData  = groupApp.create_dataset("code",(1,),dtype="S5")
  CodeData[0:] = "PySCF"
  CodeVer  = groupApp.create_dataset("version",(3,),dtype="i4")
  CodeVer[0:] = 1
  CodeVer[1:] = 4
  CodeVer[2:] = 2

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



  #print 'l_atoms', l_atoms
  #print 'uniq_atoms', uniq_atoms.items()
  #print 'idxS2A', idxSpeciestoAtoms.items()
  #print 'idxA2S', idxAtomstoSpecies.items()
 
  NbSpecies=len(idxSpeciestoAtoms.keys())


  groupAtom.create_dataset("number_of_species",(1,),dtype="i4",data=NbSpecies)


  #Dataset positions 
  MyPos=groupAtom.create_dataset("positions",(natom,3),dtype="f8")
  for x in range(natom): 
    MyPos[x:]=cell.atom_coord(x)
    print 'Positions:',cell.atom_coord(x)

  #Group Atoms
  for x in range(NbSpecies):
    atmname=str(uniq_atoms[x][0])
    groupSpecies=groupAtom.create_group("species_"+str(x))
    groupSpecies.create_dataset("atomic_number",(1,),dtype="i4",data=uniq_atoms[x][1])
    mylen="S"+str(len(atmname))
    AtmName=groupSpecies.create_dataset("name",(1,),dtype=mylen)
    AtmName[0:]=atmname
    groupSpecies.create_dataset("charge",(1,),dtype="f8",data=uniq_atoms[x][2])
    groupSpecies.create_dataset("core",(1,),dtype="f8",data=uniq_atoms[x][3])
  SpeciesID=groupAtom.create_dataset("species_ids",(natom,),dtype="i4")
  for x in range(natom):
  	SpeciesID[x:]  = idxAtomstoSpecies[x]


   



  #print 'Total number of contracted GTO', cell.nao_cart()
  #print 'Total number of NR', cell.nao_nr()
  #print 'offset of every Shell', cell.ao_loc_2c()
  #print 'aoslice NR by atom',cell.aoslice_nr_by_atom()
  #print 'Nb of Shells per atom: 1',cell.atom_nshells(0)
  #print 'Nb of Shells per atom: 2',cell.atom_nshells(1)
 
  #for j in range(natom):
  #  for i in range(cell.atom_nshells(j)):
       #print 'shell ID per Atom; Atom=',j, 'ShellId=',i, cell.atom_shell_ids(j)[i]
       #print ' Contraction Coeeficient for Atom ',j,' ShellID=',i,cell.bas_exp(cell.atom_shell_ids(j)[i])
       #print ' Atom ',j,' ShellID=',i,'length of Cartesian',cell.bas_nprim(cell.atom_shell_ids(j))


  #Parameter Group
  GroupParameter=H5_qmcpack.create_group("parameters")
   
  if cell.ecp=={}:
     ECP=False
  else:
     ECP=True
 
  GroupParameter.create_dataset("ECP",(1,),dtype="b1",data=ECP)


  bohrUnit=True
  Spin=cell.spin 

  GroupParameter.create_dataset("Unit",(1,),dtype="b1",data=bohrUnit) 
  GroupParameter.create_dataset("NbAlpha",(1,),dtype="i4",data=cell.nelec[0]) 
  GroupParameter.create_dataset("NbBeta",(1,),dtype="i4",data=cell.nelec[1]) 
  GroupParameter.create_dataset("NbTotElec",(1,),dtype="i4",data=cell.nelec[0]+cell.nelec[1])
  GroupParameter.create_dataset("spin",(1,),dtype="i4",data=Spin) 

  

#  print 'Spin=',Spin
#  print 'Total Nb of Electrons=',cell.nelec[0]+cell.nelec[1]
#  print 'Total Nb of Electrons Alpha=',cell.nelec[0]
#  print 'Total Nb of Electrons Beta=',cell.nelec[1]





  #basisset Group
  GroupBasisSet=H5_qmcpack.create_group("basisset")
  #Dataset Number Of Atoms
  GroupBasisSet.create_dataset("NbElements",(1,),dtype="i4",data=NbSpecies)

  LCAOName=GroupBasisSet.create_dataset("name",(1,),dtype="S8")
  LCAOName[0:]="LCAOBset"

  #atomicBasisSets Group
  for x in range(NbSpecies):
    atomicBasisSetGroup=GroupBasisSet.create_group("atomicBasisSet"+str(x))
    #Dataset NbBasisGroups
    atomicBasisSetGroup.create_dataset("NbBasisGroups",(1,),dtype="i4",data=cell.atom_nshells(idxAtomstoSpecies[x]))
    Angular=atomicBasisSetGroup.create_dataset("angular",(1,),dtype="S9")
    Angular[0:]="cartesian"
    for i in range(cell.atom_nshells(idxAtomstoSpecies[x])):
      BasisGroup=atomicBasisSetGroup.create_group("basisGroup"+str(i))
      BasisGroup.create_dataset("NbRadFunc",(1,),dtype="i4",data=cell.bas_nprim(cell.atom_shell_ids(idxAtomstoSpecies[x])[i]))
      Val_l=BasisGroup.create_dataset("l",(1,),dtype="i4",data=cell.bas_angular(cell.atom_shell_ids(idxAtomstoSpecies[x])[i]))
      Val_n=BasisGroup.create_dataset("n",(1,),dtype="i4",data=i)
      RadGroup=BasisGroup.create_group("radfunctions")
      for j in range(cell.bas_nprim(cell.atom_shell_ids(idxAtomstoSpecies[x])[i])):
         DataRadGrp=RadGroup.create_group("DataRad"+str(j))
         DataRadGrp.create_dataset("contraction",(1,),dtype="f8",data=cell.bas_ctr_coeff(cell.atom_shell_ids(idxAtomstoSpecies[x])[i])[j])
         DataRadGrp.create_dataset("exponent",(1,),dtype="f8",data=cell.bas_exp(cell.atom_shell_ids(idxAtomstoSpecies[x])[i])[j])
      mylen="S"+str(len((uniq_atoms[x][0]+str(i)+str(cell.bas_angular(cell.atom_shell_ids(idxAtomstoSpecies[x])[i])))))
      RID=BasisGroup.create_dataset("rid",(1,),dtype=mylen)
      RID[0:]=(uniq_atoms[x][0]+str(i)+str(cell.bas_angular(cell.atom_shell_ids(idxAtomstoSpecies[x])[i])))
      basisType=BasisGroup.create_dataset("type",(1,),dtype="S8")
      basisType[0:]="Gaussian"
    mylen="S"+str(len(uniq_atoms[x][0]))
    elemtype=atomicBasisSetGroup.create_dataset("elementType",(1,),dtype=mylen)
    elemtype[0:]=uniq_atoms[x][0]
    atomicBasisSetGroup.create_dataset("grid_npts",(1,),dtype="i4",data=1001)
    atomicBasisSetGroup.create_dataset("grid_rf",(1,),dtype="i4",data=100)
    atomicBasisSetGroup.create_dataset("grid_ri",(1,),dtype="f8",data=1e-06)
    #atomicBasisSetGroup.create_dataset("grid_type",(1,),dtype=dt,data="log")
    gridType=atomicBasisSetGroup.create_dataset("grid_type",(1,),dtype="S3")
    gridType[0:]="log"
     
    mylen="S"+str(len(cell.basis))
    nameBase=atomicBasisSetGroup.create_dataset("name",(1,),dtype=mylen)
    nameBase[0:]=cell.basis
    Normalized=atomicBasisSetGroup.create_dataset("normalized",(1,),dtype="S2")
    Normalized[0:]="no"
 
    



  if PBC==False:
    GroupDet=H5_qmcpack.create_group("determinant")

    if (RHF==True) or (ROHF==True): 
      GroupParameter.create_dataset("SpinResticted",(1,),dtype="b1",data=True)
      mo_coeff = mf.mo_coeff
      NbMO=len(mo_coeff)
      NbAO=len(mo_coeff[0])
      MO=GroupParameter.create_dataset("numMO",(1,),dtype="i4",data=NbMO)
      AO=GroupParameter.create_dataset("numAO",(1,),dtype="i4",data=NbAO)
      eigenset=GroupDet.create_dataset("eigenset_0",(NbMO,NbAO),dtype="f8",data=mo_coeff)
    else:
      GroupParameter.create_dataset("SpinResticted",(1,),dtype="b1",data=False)
      mo_coeff = mf.mo_coeff
      NbMO_up=len(mo_coeff[0])
      NbAO_up=len(mo_coeff[0][0])
      eigenset_up=GroupDet.create_dataset("eigenset_0",(NbMO_up,NbAO_up),dtype="f8",data=mo_coeff[0])

      NbMO_dn=len(mo_coeff[1])
      NbAO_dn=len(mo_coeff[1][0])
      eigenset_dn=GroupDet.create_dataset("eigenset_1",(NbMO_dn,NbAO_dn),dtype="f8",data=mo_coeff[1])
      
      
      MO_UP=GroupParameter.create_dataset("numMO_up",(1,),dtype="i4",data=NbMO_up)
      AO_UP=GroupParameter.create_dataset("numAO_up",(1,),dtype="i4",data=NbAO_up)

      MO_DN=GroupParameter.create_dataset("numMO_dn",(1,),dtype="i4",data=NbMO_dn)
      AO_DN=GroupParameter.create_dataset("numAO_dn",(1,),dtype="i4",data=NbAO_dn)
      
  else:
    #Electrons Group
    GroupElectron=H5_qmcpack.create_group("electrons")
   

  # Close the file before exiting
  H5_qmcpack.close()
