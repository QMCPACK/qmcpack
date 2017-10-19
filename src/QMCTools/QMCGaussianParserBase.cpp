//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "QMCTools/QMCGaussianParserBase.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Utilities/OhmmsInfo.h"
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include <iterator>
#include <algorithm>
#include <numeric>
#include "io/hdf_archive.h"
#include <set>
#include <map>
#include "QMCTools/GTO2GridBuilder.h"
#include "QMCApp/InitMolecularSystem.h"

//std::vector<std::string> QMCGaussianParserBase::IonName;
const int OhmmsAsciiParser::bufferSize;
std::map<int,std::string> QMCGaussianParserBase::IonName;
std::vector<std::string> QMCGaussianParserBase::gShellType;
std::vector<int> QMCGaussianParserBase::gShellID;

QMCGaussianParserBase::QMCGaussianParserBase():
  Title("sample"),basisType("Gaussian"),basisName("generic"),DoCusp(false),
  Normalized("no"),gridPtr(0),multideterminant(false),ci_threshold(0.01)
  ,usingCSF(false),readNO(0),readGuess(0),zeroCI(false),target_state(0)
  ,orderByExcitation(false), addJastrow(true), addJastrow3Body(false),QP(false)
{
}

QMCGaussianParserBase::QMCGaussianParserBase(int argc, char** argv):
  BohrUnit(true),SpinRestricted(false),NumberOfAtoms(0),NumberOfEls(0),DoCusp(false),
  SpinMultiplicity(0),NumberOfAlpha(0),NumberOfBeta(0),SizeOfBasisSet(0),
  Title("sample"),basisType("Gaussian"),basisName("generic"),numMO(0),numMO2print(-1),
  Normalized("no"),gridPtr(0),multideterminant(false),ci_threshold(0.01),target_state(0),
  angular_type("spherical"),usingCSF(false),readNO(0),readGuess(0),zeroCI(false)
  ,orderByExcitation(false), addJastrow(true), addJastrow3Body(false), QP(false)
{
  //IonSystem.setName("i");
  IonChargeIndex=IonSystem.getSpeciesSet().addAttribute("charge");
  ValenceChargeIndex=IonSystem.getSpeciesSet().addAttribute("valence");
  AtomicNumberIndex=IonSystem.getSpeciesSet().addAttribute("atomicnumber");
  std::cout << "Index of ion charge " << IonChargeIndex << std::endl;
  std::cout << "Index of valence charge " << ValenceChargeIndex << std::endl;

  createGridNode(argc,argv);
}

void QMCGaussianParserBase::init()
{
  IonName[1] = "H";
  IonName[2] = "He";
  IonName[3] = "Li";
  IonName[4] = "Be";
  IonName[5] = "B";
  IonName[6] = "C";
  IonName[7] = "N";
  IonName[8] = "O";
  IonName[9] = "F";
  IonName[10] = "Ne";
  IonName[11] = "Na";
  IonName[12] = "Mg";
  IonName[13] = "Al";
  IonName[14] = "Si";
  IonName[15] = "P";
  IonName[16] = "S";
  IonName[17] = "Cl";
  IonName[18] = "Ar";
  IonName[19] = "K";
  IonName[20] = "Ca";
  IonName[21] = "Sc";
  IonName[22] = "Ti";
  IonName[23] = "V";
  IonName[24] = "Cr";
  IonName[25] = "Mn";
  IonName[26] = "Fe";
  IonName[27] = "Co";
  IonName[28] = "Ni";
  IonName[29] = "Cu";
  IonName[30] = "Zn";
  IonName[31] = "Ga";
  IonName[32] = "Ge";
  IonName[33] = "As";
  IonName[34] = "Se";
  IonName[35] = "Br";
  IonName[36] = "Kr";
  IonName[37] = "Rb";
  IonName[38] = "Sr";
  IonName[39] = "Y";
  IonName[40] = "Zr";
  IonName[41] = "Nb";
  IonName[42] = "Mo";
  IonName[43] = "Tc";
  IonName[44] = "Ru";
  IonName[45] = "Rh";
  IonName[46] = "Pd";
  IonName[47] = "Ag";
  IonName[48] = "Cd";
  IonName[49] = "In";
  IonName[50] = "Sn";
  IonName[51] = "Sb";
  IonName[52] = "Te";
  IonName[53] = "I";
  IonName[54] = "Xe";
  IonName[55] = "Cs";
  IonName[56] = "Ba";
  IonName[57] = "La";
  IonName[58] = "Ce";
  IonName[59] = "Pr";
  IonName[60] = "Nd";
  IonName[61] = "Pm";
  IonName[62] = "Sm";
  IonName[63] = "Eu";
  IonName[64] = "Gd";
  IonName[65] = "Tb";
  IonName[66] = "Dy";
  IonName[67] = "Ho";
  IonName[68] = "Er";
  IonName[69] = "Tm";
  IonName[70] = "Yb";
  IonName[71] = "Lu";
  IonName[72] = "Hf";
  IonName[73] = "Ta";
  IonName[74] = "W";
  IonName[75] = "Re";
  IonName[76] = "Os";
  IonName[77] = "Ir";
  IonName[78] = "Pt";
  IonName[79] = "Au";
  IonName[80] = "Hg";
  IonName[81] = "Tl";
  IonName[82] = "Pb";
  IonName[83] = "Bi";
  IonName[84] = "Po";
  IonName[85] = "At";
  IonName[86] = "Rn";
  IonName[87] = "Fr";
  IonName[88] = "Ra";
  IonName[89] = "Ac";
  IonName[90] = "Th";
  IonName[91] = "Pa";
  IonName[92] = "U";
  IonName[93] = "Np";
  gShellType.resize(10);
  gShellType[1]="s";
  gShellType[2]="sp";
  gShellType[3]="p";
  gShellType[4]="d";
  gShellType[5]="f";
  gShellType[6]="g";
  gShellType[7]="h";
  gShellType[8]="h1";
  gShellType[9]="h2";
  gShellID.resize(10);
  gShellID[1]=0;
  gShellID[2]=0;
  gShellID[3]=1; //gShellID[4]=2; gShellID[5]=3; gShellID[6]=4; gShellID[7]=5;
  for(int i=4,l=2; i<gShellID.size(); ++i,++l)
    gShellID[i]=l;
}

void QMCGaussianParserBase::setOccupationNumbers()
{
  int ds=SpinMultiplicity-1;
  NumberOfBeta= (NumberOfEls-ds)/2;
  NumberOfAlpha= NumberOfEls-NumberOfBeta;
  if(!SpinRestricted)
    //UHF
  {
    std::multimap<value_type,int> e;
    //for(int i=0; i<SizeOfBasisSet; i++) e.insert(std::pair<value_type,int>(EigVal_alpha[i],0));
    //for(int i=0; i<SizeOfBasisSet; i++) e.insert(std::pair<value_type,int>(EigVal_beta[i],1));
    for(int i=0; i<numMO; i++)
      e.insert(std::pair<value_type,int>(EigVal_alpha[i],0));
    for(int i=0; i<numMO; i++)
      e.insert(std::pair<value_type,int>(EigVal_beta[i],1));
    int n=0;
    std::multimap<value_type,int>::iterator it(e.begin());
    LOGMSG("Unrestricted HF. Sorted eigen values")
    while(n<NumberOfEls && it != e.end())
    {
      LOGMSG(n << " " << (*it).first << " " << (*it).second)
      //if((*it).second == 0) {NumberOfAlpha++;}
      //else {NumberOfBeta++;}
      ++it;
      ++n;
    }
  }
  //}
  LOGMSG("Number of alpha electrons " << NumberOfAlpha)
  LOGMSG("Number of beta electrons " << NumberOfBeta)
  //Occ_alpha.resize(SizeOfBasisSet,0);
  //Occ_beta.resize(SizeOfBasisSet,0);
  Occ_alpha.resize(numMO,0);
  Occ_beta.resize(numMO,0);
  for(int i=0; i<NumberOfAlpha; i++)
    Occ_alpha[i]=1;
  for(int i=0; i<NumberOfBeta; i++)
    Occ_beta[i]=1;
}

xmlNodePtr QMCGaussianParserBase::createElectronSet()
{
//  const double ang_to_bohr=1.0/0.529177e0;
//   if(!BohrUnit) IonSystem.R *= ang_to_bohr;
//   SpeciesSet& ionSpecies(IonSystem.getSpeciesSet());
//   for(int i=0; i<NumberOfAtoms; i++) {
//     ionSpecies.addSpecies(GroupName[i]);
//   }
//   for(int i=0; i<NumberOfAtoms; i++) {
//     ionSpecies(IonChargeIndex,IonSystem.GroupID[i])=Qv[i];
//   }
  ParticleSet els;
  els.setName("e");
  std::vector<int> nel(2);
  nel[0]=NumberOfAlpha;
  nel[1]=NumberOfBeta;
  els.create(nel);
  int iu=els.getSpeciesSet().addSpecies("u");
  int id=els.getSpeciesSet().addSpecies("d");
  int ic=els.getSpeciesSet().addAttribute("charge");
  els.getSpeciesSet()(ic,iu)=-1;
  els.getSpeciesSet()(ic,id)=-1;
  //Create InitMolecularSystem to assign random electron positions
  InitMolecularSystem m(0,"test");
  if(IonSystem.getTotalNum()>1)
  {
    std::cout <<"Total number atoms "<<IonSystem.getTotalNum()<< std::endl;
    m.initAtom(&IonSystem,&els);
    //m.initMolecule(&IonSystem,&els);
  }
  else
  {
    m.initAtom(&IonSystem,&els);
  }
  XMLSaveParticle o(els);
  return o.createNode(false);
}

xmlNodePtr QMCGaussianParserBase::createESPSet(int iesp)
{
  const double ang_to_bohr=1.0/0.529177e0;
  if(!BohrUnit)
    ESPSystem[iesp].R *= ang_to_bohr;
  SpeciesSet& ionSpecies(ESPSystem[iesp].getSpeciesSet());
  for(int i=0; i<ESPSystem[iesp].getTotalNum(); i++)
  {
    ionSpecies(ESPValenceChargeIndex[iesp],i)=1;
  }
  XMLSaveParticle o(ESPSystem[iesp]);
  return o.createNode(Periodicity);
}

xmlNodePtr QMCGaussianParserBase::createIonSet()
{
  const double ang_to_bohr=1.0/0.529177e0;
  if(!BohrUnit)
    IonSystem.R *= ang_to_bohr;
  double CoreTable[] =
  {
    0, /* index zero*/
    1,2,                           /*H He */
    2,2,2,2,2,2,2,10,              /*Li-Ne*/
    10,10,10,10,10,10,10,18,       /*Na-Ar*/
    18,18,18,18,18,18,18,18,18,18, /*N-Zn*/
    28,28,28,28,28,36,             /*Ga-Kr*/
    36,36,36,36,36,36,36,36,36,36, /*Rb-Cd*/
    46,46,46,46,46,54              /*In-Xe*/
  };
  SpeciesSet& ionSpecies(IonSystem.getSpeciesSet());
  for(int i=0; i<ionSpecies.getTotalNum(); i++)
  {
    int z = static_cast<int>(ionSpecies(AtomicNumberIndex,i));
    double valence = ionSpecies(IonChargeIndex,i);
    if(valence>CoreTable[z]&& FixValence!=true)
      valence-=CoreTable[z];
    ionSpecies(ValenceChargeIndex,i)=valence;
  }
  XMLSaveParticle o(IonSystem);
  return o.createNode(Periodicity);
}



xmlNodePtr QMCGaussianParserBase::createBasisSet()
{
  xmlNodePtr bset = xmlNewNode(NULL,(const xmlChar*)"basisset");
  xmlNewProp(bset,(const xmlChar*)"name",(const xmlChar*)"LCAOBSet");
  /*
  xmlNodePtr cur = xmlAddChild(bset,xmlNewNode(NULL,(const xmlChar*)"distancetable"));
  xmlNewProp(cur,(const xmlChar*)"source",(const xmlChar*)"i");
  xmlNewProp(cur,(const xmlChar*)"target",(const xmlChar*)"e");
  */
  xmlNodePtr cur=NULL;
  std::map<int,int> species;
  int gtot = 0;
  for(int iat=0; iat<NumberOfAtoms; iat++)
  {
    int itype = IonSystem.GroupID[iat];
    int ng = 0;
    std::map<int,int>::iterator it=species.find(itype);
    if(it == species.end())
    {
      for(int ig=gBound[iat]; ig<gBound[iat+1]; ig++)
      {
        ng += gNumber[ig];
      }
      species[itype] = ng;
      if(cur)
      {
        cur = xmlAddSibling(cur,createCenter(iat,gtot));
      }
      else
      {
        cur = xmlAddChild(bset,createCenter(iat,gtot));
      }
    }
    else
    {
      ng = (*it).second;
    }
    gtot += ng;
  }
  return bset;
}



xmlNodePtr QMCGaussianParserBase::createBasisSetWithHDF5()
{
  int counter=0;

  xmlNodePtr bset = xmlNewNode(NULL,(const xmlChar*)"basisset");
  hdf_archive hout(0); 
  hout.create(h5file.c_str(),H5F_ACC_TRUNC);
  hout.push("basisset",true);
  std::string BasisSetName("LCAOBSet");
  hout.write(BasisSetName,"name");
  hout.push("atomicBasisSet",true);
  
  std::map<int,int> species;
  int gtot = 0;
  for(int iat=0; iat<NumberOfAtoms; iat++)
  {
    int itype = IonSystem.GroupID[iat];
    int ng = 0;
    std::map<int,int>::iterator it=species.find(itype);
    if(it == species.end())
    {
      for(int ig=gBound[iat]; ig<gBound[iat+1]; ig++)
      {
        ng += gNumber[ig];
      }
      species[itype] = ng;
      createCenterH5(iat,gtot,counter);
      counter++;
    }
    else
    {
      ng = (*it).second;
    }
    gtot += ng;
  }
  hout.write(counter,"NumElemets");
  hout.close();
  return bset;
}


xmlNodePtr
QMCGaussianParserBase::createDeterminantSetWithHDF5()
{
  setOccupationNumbers();
  xmlNodePtr slaterdet = xmlNewNode(NULL,(const xmlChar*)"slaterdeterminant");
  std::ostringstream up_size, down_size, b_size;
  up_size <<NumberOfAlpha;
  down_size << NumberOfBeta;
  b_size<<numMO;
  //create a determinant Up
  xmlNodePtr udet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(udet,(const xmlChar*)"id",(const xmlChar*)"updet");
  xmlNewProp(udet,(const xmlChar*)"size",(const xmlChar*)up_size.str().c_str());
  if (DoCusp==true)
     xmlNewProp(udet,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/updet.cuspInfo.xml");

  //add occupation
  xmlNodePtr occ_data = xmlNewNode(NULL,(const xmlChar*)"occupation");
  xmlNewProp(occ_data,(const xmlChar*)"mode",(const xmlChar*)"ground");
  xmlAddChild(udet,occ_data);
  //add coefficients
  xmlNodePtr coeff_data = xmlNewNode(NULL,(const xmlChar*)"coefficient");
  xmlNewProp(coeff_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
  xmlNewProp(coeff_data,(const xmlChar*)"spindataset",(const xmlChar*)"0");
  xmlAddChild(udet,coeff_data);
  //add udet to slaterdet
  xmlNodePtr cur = xmlAddChild(slaterdet,udet);
  std::vector<int> dim(2, SizeOfBasisSet);
  dim[numMO];
  
  hdf_archive hout(0);
  hout.open(h5file.c_str(),H5F_ACC_RDWR);
  hout.push("basisset");
  hout.push("determinant_0",true);

  Matrix<double> Ctemp(SizeOfBasisSet,SizeOfBasisSet);

  int n=0;  
  for (int i=0; i<SizeOfBasisSet; i++)
     for (int j=0; j<SizeOfBasisSet; j++){
         Ctemp[i][j]=EigVec[n];
         n++;
     }

  hout.write(Ctemp,"eigenset_0");


  xmlNodePtr ddet;
  if(SpinRestricted)
  {
    ddet = xmlCopyNode(udet,1);
    xmlSetProp(ddet,(const xmlChar*)"id",(const xmlChar*)"downdet");
    xmlSetProp(ddet,(const xmlChar*)"size",(const xmlChar*)down_size.str().c_str());
    if (DoCusp==true)
       xmlNewProp(udet,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/updet.cuspInfo.xml");
  }
  else
  {
    ddet = xmlCopyNode(udet,2);
    xmlSetProp(ddet,(const xmlChar*)"id",(const xmlChar*)"downdet");
    xmlSetProp(ddet,(const xmlChar*)"size",(const xmlChar*)down_size.str().c_str());
    if (DoCusp==true)
       xmlNewProp(udet,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/updet.cuspInfo.xml");
    xmlNodePtr o= xmlAddChild(ddet,xmlCopyNode(occ_data,1));
    xmlNodePtr c= xmlCopyNode(coeff_data,1);
    xmlSetProp(c,(const xmlChar*)"spindataset",(const xmlChar*)"1");
    o = xmlAddSibling(o,c);

  n=numMO*SizeOfBasisSet;  
  for (int i=0; i<SizeOfBasisSet; i++)
     for (int j=0; j<SizeOfBasisSet; j++){
         Ctemp[i][j]=EigVec[n];
         n++;
     }

   hout.write(Ctemp,"eigenset_1");

  }
  cur = xmlAddSibling(cur,ddet);
    hout.close();
  //return slaterdeterminant node
  return slaterdet;
}




xmlNodePtr
QMCGaussianParserBase::createDeterminantSet()
{
  setOccupationNumbers();
  xmlNodePtr slaterdet = xmlNewNode(NULL,(const xmlChar*)"slaterdeterminant");
  //check spin-dependent properties
  //int nup = NumberOfEls/2;
  //int ndown = NumberOfEls-nup;
  std::ostringstream up_size, down_size, b_size, occ;
  //up_size <<NumberOfAlpha; down_size << NumberOfBeta; b_size<<SizeOfBasisSet;
  up_size <<NumberOfAlpha;
  down_size << NumberOfBeta;
  b_size<<numMO;
  //create a determinant Up
  xmlNodePtr adet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(adet,(const xmlChar*)"id",(const xmlChar*)"updet");
  xmlNewProp(adet,(const xmlChar*)"size",(const xmlChar*)up_size.str().c_str());
  if (DoCusp==true)
     xmlNewProp(adet,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/updet.cuspInfo.xml");

  //occ<<"\n";
  //vector<int>::iterator it(Occ_alpha.begin());
  //int i=0;
  //while(i<SizeOfBasisSet) {
  //  int n = (i+10<SizeOfBasisSet)? 10 : SizeOfBasisSet-i;
  //  copy(it, it+n, std::ostream_iterator<int>(occ," "));
  //  occ << "\n"; it += 10; i+=10;
  //}
  //xmlNodePtr occ_data
  //  = xmlNewTextChild(adet,NULL,(const xmlChar*)"occupation",(const xmlChar*)occ.str().c_str());
  //xmlNewProp(occ_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
  xmlNodePtr occ_data = xmlNewNode(NULL,(const xmlChar*)"occupation");
  xmlNewProp(occ_data,(const xmlChar*)"mode",(const xmlChar*)"ground");
  xmlAddChild(adet,occ_data);
  //int btot=SizeOfBasisSet*SizeOfBasisSet;
  int btot=numMO*SizeOfBasisSet;
  int n=btot/4, b=0;
  int dn=btot-n*4;
  std::ostringstream eig;
  eig.setf(std::ios::scientific, std::ios::floatfield);
  eig.setf(std::ios::right,std::ios::adjustfield);
  eig.precision(14);
  eig << "\n";
  for(int k=0; k<n; k++)
  {
    eig << std::setw(22) << EigVec[b] << std::setw(22) << EigVec[b+1] << std::setw(22) << EigVec[b+2] << std::setw(22) <<  EigVec[b+3] << "\n";
    b += 4;
  }
  for(int k=0; k<dn; k++)
  {
    eig << std::setw(22) << EigVec[b++];
  }
  if(dn)
    eig << std::endl;
  xmlNodePtr det_data
  = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eig.str().c_str());
  xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
  xmlNewProp(det_data,(const xmlChar*)"id",(const xmlChar*)"updetC");
  xmlNodePtr cur = xmlAddChild(slaterdet,adet);
  adet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(adet,(const xmlChar*)"id",(const xmlChar*)"downdet");
  xmlNewProp(adet,(const xmlChar*)"size",(const xmlChar*)down_size.str().c_str());
  if (DoCusp==true)
     xmlNewProp(adet,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/downdet.cuspInfo.xml");

  {
    //std::ostringstream occ_beta;
    //occ_beta<<"\n";
    //it=Occ_beta.begin();
    //int i=0;
    //while(i<SizeOfBasisSet) {
    //  int n = (i+10<SizeOfBasisSet)? 10 : SizeOfBasisSet-i;
    //  copy(it, it+n, std::ostream_iterator<int>(occ_beta," "));
    //  occ_beta << "\n"; it += 10; i+=10;
    //}
    //occ_data=xmlNewTextChild(adet,NULL,(const xmlChar*)"occupation",(const xmlChar*)occ_beta.str().c_str());
    //xmlNewProp(occ_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
    occ_data = xmlNewNode(NULL,(const xmlChar*)"occupation");
    xmlNewProp(occ_data,(const xmlChar*)"mode",(const xmlChar*)"ground");
    xmlAddChild(adet,occ_data);
    std::ostringstream eigD;
    eigD.setf(std::ios::scientific, std::ios::floatfield);
    eigD.setf(std::ios::right,std::ios::adjustfield);
    eigD.precision(14);
    eigD << "\n";
    //b=SizeOfBasisSet*SizeOfBasisSet;
    b=numMO*SizeOfBasisSet;
    for(int k=0; k<n; k++)
    {
      eigD << std::setw(22) << EigVec[b] << std::setw(22) << EigVec[b+1] << std::setw(22) << EigVec[b+2] << std::setw(22) <<  EigVec[b+3] << "\n";
      b += 4;
    }
    for(int k=0; k<dn; k++)
    {
      eigD << std::setw(22) << EigVec[b++];
    }
    if(dn)
      eigD << std::endl;
    if(SpinRestricted)
      det_data
      = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eig.str().c_str());
    else
      det_data
      = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eigD.str().c_str());
    xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
    xmlNewProp(det_data,(const xmlChar*)"id",(const xmlChar*)"downdetC");
  }
  cur = xmlAddSibling(cur,adet);
  return slaterdet;
}

void
QMCGaussianParserBase::createSPOSets(xmlNodePtr spoUP, xmlNodePtr spoDN)
{
  setOccupationNumbers();
  std::ostringstream up_size, down_size, b_size, occ, nstates_alpha,nstates_beta;
  //up_size <<NumberOfAlpha; down_size << NumberOfBeta; b_size<<SizeOfBasisSet;
  up_size <<NumberOfAlpha;
  down_size << NumberOfBeta;
  b_size<<numMO;
  nstates_alpha <<ci_nstates+ci_nca;;
  nstates_beta <<ci_nstates+ci_ncb;

  xmlNewProp(spoUP,(const xmlChar*)"name",(const xmlChar*)"spo-up");
  xmlNewProp(spoDN,(const xmlChar*)"name",(const xmlChar*)"spo-dn");
  xmlNewProp(spoUP,(const xmlChar*)"size",(const xmlChar*)nstates_alpha.str().c_str());
  xmlNewProp(spoDN,(const xmlChar*)"size",(const xmlChar*)nstates_beta.str().c_str());
  if (DoCusp==true){
     xmlNewProp(spoUP,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/spo-up.cuspInfo.xml");
     xmlNewProp(spoDN,(const xmlChar*)"cuspInfo",(const xmlChar*)"../CuspCorrection/spo-dn.cuspInfo.xml");
  }
  xmlNodePtr occ_data = xmlNewNode(NULL,(const xmlChar*)"occupation");
  xmlNewProp(occ_data,(const xmlChar*)"mode",(const xmlChar*)"ground");
  xmlAddChild(spoUP,occ_data);
  //int btot=SizeOfBasisSet*SizeOfBasisSet;
  int btot=numMO*SizeOfBasisSet;
  int n=btot/4, b=0;
  int dn=btot-n*4;
  std::ostringstream eig;
  eig.setf(std::ios::scientific, std::ios::floatfield);
  eig.setf(std::ios::right,std::ios::adjustfield);
  eig.precision(14);
  eig << "\n";
  for(int k=0; k<n; k++)
  {
    eig << std::setw(22) << EigVec[b] << std::setw(22) << EigVec[b+1] << std::setw(22) << EigVec[b+2] << std::setw(22) <<  EigVec[b+3] << "\n";
    b += 4;
  }
  for(int k=0; k<dn; k++)
  {
    eig << std::setw(22) << EigVec[b++];
  }
  if(dn)
    eig << std::endl;
  xmlNodePtr det_data
  = xmlNewTextChild(spoUP,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eig.str().c_str());
  xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
  xmlNewProp(det_data,(const xmlChar*)"id",(const xmlChar*)"updetC");
  {
    occ_data = xmlNewNode(NULL,(const xmlChar*)"occupation");
    xmlNewProp(occ_data,(const xmlChar*)"mode",(const xmlChar*)"ground");
    xmlAddChild(spoDN,occ_data);
    std::ostringstream eigD;
    eigD.setf(std::ios::scientific, std::ios::floatfield);
    eigD.setf(std::ios::right,std::ios::adjustfield);
    eigD.precision(14);
    eigD << "\n";
    //b=SizeOfBasisSet*SizeOfBasisSet;
    b=numMO*SizeOfBasisSet;
    for(int k=0; k<n; k++)
    {
      eigD << std::setw(22) << EigVec[b] << std::setw(22) << EigVec[b+1] << std::setw(22) << EigVec[b+2] << std::setw(22) <<  EigVec[b+3] << "\n";
      b += 4;
    }
    for(int k=0; k<dn; k++)
    {
      eigD << std::setw(22) << EigVec[b++];
    }
    if(dn)
      eigD << std::endl;
    if(SpinRestricted)
      det_data
      = xmlNewTextChild(spoDN,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eig.str().c_str());
    else
      det_data
      = xmlNewTextChild(spoDN,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eigD.str().c_str());
    xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
    xmlNewProp(det_data,(const xmlChar*)"id",(const xmlChar*)"downdetC");
  }
}

xmlNodePtr
QMCGaussianParserBase::createMultiDeterminantSetQP()
{
    xmlNodePtr multislaterdet = xmlNewNode(NULL,(const xmlChar*)"multideterminant");
    xmlNewProp(multislaterdet,(const xmlChar*)"optimize",(const xmlChar*)"yes");
    xmlNewProp(multislaterdet,(const xmlChar*)"spo_up",(const xmlChar*)"spo-up");
    xmlNewProp(multislaterdet,(const xmlChar*)"spo_dn",(const xmlChar*)"spo-dn");
    xmlNodePtr detlist = xmlNewNode(NULL,(const xmlChar*)"detlist");
    std::ostringstream nstates,cisize,cinca,cincb,cinea,cineb,ci_thr;
    cisize <<ci_size;
    nstates <<ci_nstates;
    cinca <<ci_nca;
    cincb <<ci_ncb;
    cinea <<ci_nea;
    cineb <<ci_neb;
    ci_thr <<ci_threshold;
    xmlNewProp(detlist,(const xmlChar*)"size",(const xmlChar*)cisize.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"type",(const xmlChar*)"DETS");
    xmlNewProp(detlist,(const xmlChar*)"nca",(const xmlChar*)cinca.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"ncb",(const xmlChar*)cincb.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"nea",(const xmlChar*)cinea.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"neb",(const xmlChar*)cineb.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"nstates",(const xmlChar*)nstates.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"cutoff",(const xmlChar*)ci_thr.str().c_str());
    if(CIcoeff.size() == 0)
    {
      std::cerr <<" CI configuration list is empty. \n";
      exit(101);
    }
    if(CIcoeff.size() != CIalpha.size() || CIcoeff.size() != CIbeta.size())
    {
      std::cerr <<" Problem with CI configuration lists. \n";
      exit(102);
    }
    int iv=0;
    for(int i=0; i<CIcoeff.size(); i++)
    {
        xmlNodePtr ci = xmlNewNode(NULL,(const xmlChar*)"ci");
        std::ostringstream coeff;
        std::ostringstream qc_coeff;
        qc_coeff<<CIcoeff[i];
        coeff<<CIcoeff[i];
       /*
        if(zeroCI && i==0)
        {
          coeff<<1.0;
        }
        else
          if(zeroCI && i>0)
          {
            coeff<<0.0;
          }
          else
          {
            coeff<<CIcoeff[i];
          }
       */
        std::ostringstream tag;
        tag<<"CIcoeff_" <<iv++;
        xmlNewProp(ci,(const xmlChar*)"id",(const xmlChar*) tag.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"coeff",(const xmlChar*) coeff.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"qc_coeff",(const xmlChar*) qc_coeff.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"alpha",(const xmlChar*) CIalpha[i].substr(0,ci_nstates).c_str());
        xmlNewProp(ci,(const xmlChar*)"beta",(const xmlChar*) CIbeta[i].substr(0,ci_nstates).c_str());
        xmlAddChild(detlist,ci);
      }
    xmlAddChild(multislaterdet,detlist);
  return multislaterdet;
  }

xmlNodePtr
QMCGaussianParserBase::createMultiDeterminantSetVSVB()
{
  xmlNodePtr multislaterdet = xmlNewNode(NULL,(const xmlChar*)"multideterminant");
  xmlNewProp(multislaterdet,(const xmlChar*)"optimize",(const xmlChar*)"yes");
  xmlNewProp(multislaterdet,(const xmlChar*)"spo_up",(const xmlChar*)"spo-up");
  xmlNewProp(multislaterdet,(const xmlChar*)"spo_dn",(const xmlChar*)"spo-dn");
  xmlNodePtr detlist = xmlNewNode(NULL,(const xmlChar*)"detlist");
  std::ostringstream nstates,cisize,cinca,cincb,cinea,cineb,ci_thr;
  cisize <<ci_size;
  nstates <<ci_nstates;
  cinca <<ci_nca;
  cincb <<ci_ncb;
  cinea <<ci_nea;
  cineb <<ci_neb;
  ci_thr <<ci_threshold;
  xmlNewProp(detlist,(const xmlChar*)"size",(const xmlChar*)cisize.str().c_str());
  xmlNewProp(detlist,(const xmlChar*)"type",(const xmlChar*)"CSF");
  xmlNewProp(detlist,(const xmlChar*)"nca",(const xmlChar*)cinca.str().c_str());
  xmlNewProp(detlist,(const xmlChar*)"ncb",(const xmlChar*)cincb.str().c_str());
  xmlNewProp(detlist,(const xmlChar*)"nea",(const xmlChar*)cinea.str().c_str());
  xmlNewProp(detlist,(const xmlChar*)"neb",(const xmlChar*)cineb.str().c_str());
  xmlNewProp(detlist,(const xmlChar*)"nstates",(const xmlChar*)nstates.str().c_str());
  xmlNewProp(detlist,(const xmlChar*)"cutoff",(const xmlChar*)ci_thr.str().c_str());
  CIexcitLVL.clear();
  for(int i=0; i<CSFocc.size(); i++)
  {
    CIexcitLVL.push_back(0);
    //cout<<CSFocc[i] <<" " <<CIexcitLVL.back() << std::endl;
  }
  // order dets according to ci coeff
  std::vector<std::pair<double,int> > order;
  for(int i=0; i<coeff2csf.size(); i++)
  {
    std::pair<double,int> cic(std::abs(coeff2csf[i].second),i);
    order.push_back(cic);
  }
  std::vector<std::pair<double,int> >::reverse_iterator it(order.rbegin());
  std::vector<std::pair<double,int> >::reverse_iterator last(order.rend());
  int iv=0;
  while(it != last)
  {
    int nq = (*it).second;
    xmlNodePtr csf = xmlNewNode(NULL,(const xmlChar*)"csf");
    std::ostringstream qc_coeff;
    qc_coeff<<coeff2csf[nq].second;
    std::ostringstream coeff;
    std::ostringstream exct;
    exct<<CIexcitLVL[nq];
    coeff<<coeff2csf[nq].second;
    std::ostringstream tag;
    tag<<"CSFcoeff_" <<iv;
    xmlNewProp(csf,(const xmlChar*)"id",(const xmlChar*) tag.str().c_str());
    xmlNewProp(csf,(const xmlChar*)"exctLvl",(const xmlChar*) exct.str().c_str());
    xmlNewProp(csf,(const xmlChar*)"coeff",(const xmlChar*) coeff.str().c_str());
    xmlNewProp(csf,(const xmlChar*)"qchem_coeff",(const xmlChar*) qc_coeff.str().c_str());
    xmlNewProp(csf,(const xmlChar*)"occ",(const xmlChar*) CSFocc[nq].substr(0,ci_nstates).c_str());
    for(int i=0; i<CSFexpansion[nq].size(); i++)
    {
      xmlNodePtr ci = xmlNewNode(NULL,(const xmlChar*)"det");
      std::ostringstream coeff0;
      coeff0<<CSFexpansion[nq][i];
      std::ostringstream tag0;
      tag0<<"csf_" <<iv <<"-" <<i;
      xmlNewProp(ci,(const xmlChar*)"id",(const xmlChar*) tag0.str().c_str());
      xmlNewProp(ci,(const xmlChar*)"coeff",(const xmlChar*) coeff0.str().c_str());
      xmlNewProp(ci,(const xmlChar*)"alpha",(const xmlChar*) CSFalpha[nq][i].substr(0,ci_nstates).c_str());
      xmlNewProp(ci,(const xmlChar*)"beta",(const xmlChar*) CSFbeta[nq][i].substr(0,ci_nstates).c_str());
      xmlAddChild(csf,ci);
    }
    xmlAddChild(detlist,csf);
    it++;
    iv++;
  }
  xmlAddChild(multislaterdet,detlist);
  return multislaterdet;
}
xmlNodePtr
QMCGaussianParserBase::createMultiDeterminantSet()
{
  xmlNodePtr multislaterdet = xmlNewNode(NULL,(const xmlChar*)"multideterminant");
  xmlNewProp(multislaterdet,(const xmlChar*)"optimize",(const xmlChar*)"yes");
  xmlNewProp(multislaterdet,(const xmlChar*)"spo_up",(const xmlChar*)"spo-up");
  xmlNewProp(multislaterdet,(const xmlChar*)"spo_dn",(const xmlChar*)"spo-dn");
  if(usingCSF)
  {
    xmlNodePtr detlist = xmlNewNode(NULL,(const xmlChar*)"detlist");
    std::ostringstream nstates,cisize,cinca,cincb,cinea,cineb,ci_thr;
    cisize <<ci_size;
    nstates <<ci_nstates;
    cinca <<ci_nca;
    cincb <<ci_ncb;
    cinea <<ci_nea;
    cineb <<ci_neb;
    ci_thr <<ci_threshold;
    xmlNewProp(detlist,(const xmlChar*)"size",(const xmlChar*)cisize.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"type",(const xmlChar*)"CSF");
    xmlNewProp(detlist,(const xmlChar*)"nca",(const xmlChar*)cinca.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"ncb",(const xmlChar*)cincb.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"nea",(const xmlChar*)cinea.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"neb",(const xmlChar*)cineb.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"nstates",(const xmlChar*)nstates.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"cutoff",(const xmlChar*)ci_thr.str().c_str());
    //std::vector<std::pair<int,double> >::iterator it(coeff2csf.begin());
    //std::vector<std::string>::iterator occit(CSFocc.begin());
    //std::vector<std::pair<int,double> >::iterator last(coeff2csf.end());
    CIexcitLVL.clear();
    for(int i=0; i<CSFocc.size(); i++)
    {
      CIexcitLVL.push_back(numberOfExcitationsCSF(CSFocc[i]));
      //cout<<CSFocc[i] <<" " <<CIexcitLVL.back() << std::endl;
    }
    // order dets according to ci coeff
    std::vector<std::pair<double,int> > order;
    if(orderByExcitation)
    {
      std::cout <<"Ordering csfs by excitation level. \n";
      int maxE =  *max_element(CIexcitLVL.begin(),CIexcitLVL.end());
      std::vector<int> pos(maxE);
      int ip1,ip2,cnt=0,cnt2;
      // add by excitations, and do partial sorts of the list
      // messy but I dont want std::pair< std::pair<> > types right now
      for(int i=maxE; i>=0; i--)
      {
        ip1 = ip2 = cnt;
        cnt2=0;
        for(int k=0; k<CIexcitLVL.size(); k++)
        {
          if(CIexcitLVL[k] == i)
          {
            std::pair<double,int> cic(std::abs(coeff2csf[k].second),k);
            order.push_back(cic);
            cnt2++;
            cnt++;
          }
        }
        if(cnt2 > 0)
          sort(order.begin()+ip1,order.end());
      }
    }
    else
    {
      for(int i=0; i<coeff2csf.size(); i++)
      {
        std::pair<double,int> cic(std::abs(coeff2csf[i].second),i);
        order.push_back(cic);
      }
      sort(order.begin(),order.end());
    }
    std::vector<std::pair<double,int> >::reverse_iterator it(order.rbegin());
    std::vector<std::pair<double,int> >::reverse_iterator last(order.rend());
    int iv=0;
    while(it != last)
    {
      int nq = (*it).second;
      xmlNodePtr csf = xmlNewNode(NULL,(const xmlChar*)"csf");
      std::ostringstream qc_coeff;
      qc_coeff<<coeff2csf[nq].second;
      std::ostringstream coeff;
      std::ostringstream exct;
      exct<<CIexcitLVL[nq];
      if(zeroCI && iv==0)
      {
        coeff<<1.0;
      }
      else
        if(zeroCI && iv>0)
        {
          coeff<<0.0;
        }
        else
        {
          coeff<<coeff2csf[nq].second;
        }
      std::ostringstream tag;
      tag<<"CSFcoeff_" <<iv;
      xmlNewProp(csf,(const xmlChar*)"id",(const xmlChar*) tag.str().c_str());
      xmlNewProp(csf,(const xmlChar*)"exctLvl",(const xmlChar*) exct.str().c_str());
      xmlNewProp(csf,(const xmlChar*)"coeff",(const xmlChar*) coeff.str().c_str());
      xmlNewProp(csf,(const xmlChar*)"qchem_coeff",(const xmlChar*) qc_coeff.str().c_str());
      xmlNewProp(csf,(const xmlChar*)"occ",(const xmlChar*) CSFocc[nq].substr(0,ci_nstates).c_str());
      for(int i=0; i<CSFexpansion[nq].size(); i++)
      {
        xmlNodePtr ci = xmlNewNode(NULL,(const xmlChar*)"det");
        std::ostringstream coeff0;
        coeff0<<CSFexpansion[nq][i];
        std::ostringstream tag0;
        tag0<<"csf_" <<iv <<"-" <<i;
        xmlNewProp(ci,(const xmlChar*)"id",(const xmlChar*) tag0.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"coeff",(const xmlChar*) coeff0.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"alpha",(const xmlChar*) CSFalpha[nq][i].substr(0,ci_nstates).c_str());
        xmlNewProp(ci,(const xmlChar*)"beta",(const xmlChar*) CSFbeta[nq][i].substr(0,ci_nstates).c_str());
        xmlAddChild(csf,ci);
      }
      xmlAddChild(detlist,csf);
      it++;
      iv++;
    }
    xmlAddChild(multislaterdet,detlist);
  }
  else
    // usingCSF
  {
    xmlNodePtr detlist = xmlNewNode(NULL,(const xmlChar*)"detlist");
    ci_size=0;
    for(int i=0; i<CIcoeff.size(); i++)
      if(std::abs(CIcoeff[i]) > ci_threshold)
        ci_size++;
    std::ostringstream nstates,cisize,cinca,cincb,cinea,cineb,ci_thr;
    cisize <<ci_size;
    nstates <<ci_nstates;
    cinca <<ci_nca;
    cincb <<ci_ncb;
    cinea <<ci_nea;
    cineb <<ci_neb;
    ci_thr <<ci_threshold;
    xmlNewProp(detlist,(const xmlChar*)"size",(const xmlChar*)cisize.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"type",(const xmlChar*)"DETS");
    xmlNewProp(detlist,(const xmlChar*)"nca",(const xmlChar*)cinca.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"ncb",(const xmlChar*)cincb.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"nea",(const xmlChar*)cinea.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"neb",(const xmlChar*)cineb.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"nstates",(const xmlChar*)nstates.str().c_str());
    xmlNewProp(detlist,(const xmlChar*)"cutoff",(const xmlChar*)ci_thr.str().c_str());
    if(CIcoeff.size() == 0)
    {
      std::cerr <<" CI configuration list is empty. \n";
      exit(101);
    }
    if(CIcoeff.size() != CIalpha.size() || CIcoeff.size() != CIbeta.size())
    {
      std::cerr <<" Problem with CI configuration lists. \n";
      exit(102);
    }
    int iv=0;
    /*
    while(iv < CIcoeff.size() && std::abs(CIcoeff[iv]) < ci_threshold) iv++;

    {
      xmlNodePtr ci = xmlNewNode(NULL,(const xmlChar*)"ci");
      std::ostringstream coeff; coeff<<CIcoeff[iv];
      xmlNewProp(ci,(const xmlChar*)"coeff",(const xmlChar*) coeff.str().c_str());
      xmlNewProp(ci,(const xmlChar*)"alpha",(const xmlChar*) CIalpha[iv].c_str());
      xmlNewProp(ci,(const xmlChar*)"beta",(const xmlChar*) CIbeta[iv].c_str());
      xmlAddChild(detlist,ci);
      iv++;
    }
    */
    for(int i=0; i<CIcoeff.size(); i++)
    {
      if(std::abs(CIcoeff[i]) > ci_threshold)
      {
        xmlNodePtr ci = xmlNewNode(NULL,(const xmlChar*)"ci");
        std::ostringstream coeff;
        std::ostringstream qc_coeff;
        qc_coeff<<CIcoeff[i];
        if(zeroCI && i==0)
        {
          coeff<<1.0;
        }
        else
          if(zeroCI && i>0)
          {
            coeff<<0.0;
          }
          else
          {
            coeff<<CIcoeff[i];
          }
        std::ostringstream tag;
        tag<<"CIcoeff_" <<iv++;
        xmlNewProp(ci,(const xmlChar*)"id",(const xmlChar*) tag.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"coeff",(const xmlChar*) coeff.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"qc_coeff",(const xmlChar*) qc_coeff.str().c_str());
        xmlNewProp(ci,(const xmlChar*)"alpha",(const xmlChar*) CIalpha[i].substr(0,ci_nstates).c_str());
        xmlNewProp(ci,(const xmlChar*)"beta",(const xmlChar*) CIbeta[i].substr(0,ci_nstates).c_str());
        xmlAddChild(detlist,ci);
      }
    }
    xmlAddChild(multislaterdet,detlist);
  } //usingCSF
  return multislaterdet;
}

void QMCGaussianParserBase::createCenterH5(int iat, int off_,int numelem)
{

  CurrentCenter=GroupName[iat];
  int numbasisgroups(0);
  std::stringstream tempcenter;
  std::string CenterID;
  tempcenter<<CurrentCenter<<"";
  CenterID=tempcenter.str();
  std::stringstream tempElem;
  std::string ElemID0="Element",ElemID;
  tempElem<<ElemID0<<numelem;
  ElemID=tempElem.str();

  hdf_archive hout(0);
  hout.open(h5file.c_str(),H5F_ACC_RDWR);
  hout.push("basisset");
  hout.push("atomicBasisSet");
  hout.push(ElemID.c_str(),true);
  hout.write(basisName,"name");
  hout.write(angular_type,"angular");
  hout.write(basisType,"type");
  hout.write(CenterID,"elementType");
  hout.write(Normalized,"normalized");

  double ValgridFirst(1.e-6), ValgridLast(1.e2);
  int Valnpts(1001);
  std::string gridType("log");
  double gridFirst(ValgridFirst);
  double gridLast(ValgridLast);
  int gridSize(Valnpts);
  
  hout.write(gridType,"grid_type");
  hout.write(gridFirst,"grid_ri");
  hout.write(gridLast,"grid_rf");
  hout.write(gridSize,"grid_npts");

  for(int ig=gBound[iat], n=0; ig< gBound[iat+1]; ig++,n++)
  {
    createShellH5(n, ig, off_,numelem);
    off_ += gNumber[ig];
    numbasisgroups=n+1;
  }
  //std::cout<<"numbasisgroups="<<numbasisgroups<<std::endl; 
  
  hout.write(numbasisgroups,"NbBasisGroups");
  hout.close();
}


xmlNodePtr QMCGaussianParserBase::createCenter(int iat, int off_)
{
  //CurrentCenter = IonName[GroupID[iat]];
  //CurrentCenter = IonSystem.Species.speciesName[iat];

  CurrentCenter=GroupName[iat];

  xmlNodePtr abasis = xmlNewNode(NULL,(const xmlChar*)"atomicBasisSet");
  xmlNewProp(abasis,(const xmlChar*)"name",(const xmlChar*)basisName.c_str());
  //xmlNewProp(abasis,(const xmlChar*)"angular",(const xmlChar*)"spherical");
  xmlNewProp(abasis,(const xmlChar*)"angular",(const xmlChar*)angular_type.c_str());
  xmlNewProp(abasis,(const xmlChar*)"type",(const xmlChar*)basisType.c_str());
  xmlNewProp(abasis,(const xmlChar*)"elementType",(const xmlChar*)CurrentCenter.c_str());
  xmlNewProp(abasis,(const xmlChar*)"normalized",(const xmlChar*)Normalized.c_str());
  xmlAddChild(abasis,xmlCopyNode(gridPtr,1));
  for(int ig=gBound[iat], n=0; ig< gBound[iat+1]; ig++,n++)
  {
    createShell(n, ig, off_,abasis);
    off_ += gNumber[ig];
  }
  return abasis;
}


void
QMCGaussianParserBase::createShellH5(int n, int ig, int off_,int numelem)
{
  int gid(gShell[ig]);
  int ng(gNumber[ig]);



  char l_name[4],n_name[4],a_name[32];
  sprintf(a_name,"%s%d%d",CurrentCenter.c_str(),n,gShellID[gid]);
  sprintf(l_name,"%d",gShellID[gid]);
  sprintf(n_name,"%d",n);

  std::string aa_name(a_name);
  std::string an_name(n_name);
  std::string al_name(l_name);
  std::string at_name("Gaussian");
  //std::string basisGroupID="basisGroup"+aa_name;
  std::string basisGroupID="basisGroup"+an_name;

  //std::cout << "basisGroupID="<<basisGroupID<<std::endl;
  std::stringstream tempElem;
  std::string ElemID0="Element",ElemID;
  tempElem<<ElemID0<<numelem;
  ElemID=tempElem.str();

  hdf_archive hout(0);
  hout.open(h5file.c_str(),H5F_ACC_RDWR);
  hout.push("basisset");
  hout.push("atomicBasisSet");
  hout.push(ElemID.c_str());
  hout.push(basisGroupID.c_str(),true);
  hout.write(aa_name,"rid");
  hout.write(n,"n");
  hout.write(gShellID[gid],"l");
  hout.write(at_name,"type");
  hout.write(ng,"NbRadFunc");

  hout.push("radfunctions",true);

  if(gid == 2)
  {
    std::string basisGroupID2="basisGroup2"+aa_name;
    std::string valac("1");
    hout.push(basisGroupID2.c_str(),true);
    hout.write(aa_name,"rid");
    hout.write(an_name,"n");
    hout.write(valac,"l");
    hout.write(at_name,"type");
    hout.write(ng,"NbRadFunc");
    hout.push("radfunctions2",true);
    hout.pop();
    hout.pop();
  }

  for(int ig=0, i=off_; ig<ng; ig++, i++)
  {
    std::stringstream tempdata; 
    std::string dataradID0="DataRad",dataradID;
    tempdata<<dataradID0<<ig;    
    dataradID=tempdata.str();
    hout.push(dataradID.c_str(),true);
    hout.write(gExp[i],"exponent");
    hout.write(gC0[i],"contraction");
    hout.pop();
    if(gid == 2)
    {
       std::string basisGroupID2="basisGroup2"+aa_name;
       hout.push(basisGroupID2.c_str());
       hout.push("radfunctions2");
       std::stringstream tempdata2; 
       std::string datarad2ID0="DataRad2",datarad2ID;
       tempdata2<<datarad2ID0<<ig;    
       datarad2ID=tempdata2.str();
       hout.push(datarad2ID.c_str(),true);
       hout.write(gExp[i],"exponent");
       hout.write(gC1[i],"contraction");
       hout.pop();
       hout.pop();
    }
  }
  hout.close();
 
}


void
QMCGaussianParserBase::createShell(int n, int ig, int off_, xmlNodePtr abasis)
{
  int gid(gShell[ig]);
  int ng(gNumber[ig]);
  xmlNodePtr ag = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
  xmlNodePtr ag1 = 0;
  char l_name[4],n_name[4],a_name[32];
  sprintf(a_name,"%s%d%d",CurrentCenter.c_str(),n,gShellID[gid]);
  sprintf(l_name,"%d",gShellID[gid]);
  sprintf(n_name,"%d",n);
  xmlNewProp(ag,(const xmlChar*)"rid", (const xmlChar*)a_name);
  xmlNewProp(ag,(const xmlChar*)"n", (const xmlChar*)n_name);
  xmlNewProp(ag,(const xmlChar*)"l", (const xmlChar*)l_name);
  xmlNewProp(ag,(const xmlChar*)"type", (const xmlChar*)"Gaussian");
  if(gid == 2)
  {
    sprintf(a_name,"%s%d1",CurrentCenter.c_str(),n);
    ag1 = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
    xmlNewProp(ag1,(const xmlChar*)"rid", (const xmlChar*)a_name);
    xmlNewProp(ag1,(const xmlChar*)"n", (const xmlChar*)n_name);
    xmlNewProp(ag1,(const xmlChar*)"l", (const xmlChar*)"1");
    xmlNewProp(ag1,(const xmlChar*)"type", (const xmlChar*)"Gaussian");
  }
  for(int ig=0, i=off_; ig<ng; ig++, i++)
  {
    std::ostringstream a,b,c;
    a.setf(std::ios::scientific, std::ios::floatfield);
    b.setf(std::ios::scientific, std::ios::floatfield);
    a.precision(12);
    b.precision(12);
    a<<gExp[i];
    b<<gC0[i];
    xmlNodePtr anode = xmlNewNode(NULL,(const xmlChar*)"radfunc");
    xmlNewProp(anode,(const xmlChar*)"exponent", (const xmlChar*)a.str().c_str());
    xmlNewProp(anode,(const xmlChar*)"contraction", (const xmlChar*)b.str().c_str());
    xmlAddChild(ag,anode);
    if(gid ==2)
    {
      c.setf(std::ios::scientific, std::ios::floatfield);
      c.precision(12);
      c <<gC1[i];
      anode = xmlNewNode(NULL,(const xmlChar*)"radfunc");
      xmlNewProp(anode,(const xmlChar*)"exponent", (const xmlChar*)a.str().c_str());
      xmlNewProp(anode,(const xmlChar*)"contraction", (const xmlChar*)c.str().c_str());
      xmlAddChild(ag1,anode);
    }
  }
  xmlAddChild(abasis,ag);
  if(gid == 2)
    xmlAddChild(abasis,ag1);
}


xmlNodePtr QMCGaussianParserBase::createJ3()
{
  xmlNodePtr j3 = xmlNewNode(NULL,(const xmlChar*)"jastrow");
  xmlNewProp(j3,(const xmlChar*)"name", (const xmlChar*)"J3");
  xmlNewProp(j3,(const xmlChar*)"type", (const xmlChar*)"eeI");
  xmlNewProp(j3,(const xmlChar*)"function", (const xmlChar*)"polynomial");
  xmlNewProp(j3,(const xmlChar*)"source", (const xmlChar*)"ion0");
  xmlNewProp(j3,(const xmlChar*)"print", (const xmlChar*)"yes");
  SpeciesSet& ionSpecies(IonSystem.getSpeciesSet());
  for(int i=0; i<ionSpecies.getTotalNum(); i++)
  {
    xmlNodePtr uuc = xmlNewNode(NULL,(const xmlChar*)"correlation");
    xmlNewProp(uuc,(const xmlChar*)"ispecies", (const xmlChar*)ionSpecies.speciesName[i].c_str());
    xmlNewProp(uuc,(const xmlChar*)"especies", (const xmlChar*)"u");
    xmlNewProp(uuc,(const xmlChar*)"isize", (const xmlChar*)"3");
    xmlNewProp(uuc,(const xmlChar*)"esize", (const xmlChar*)"3");
    xmlNewProp(uuc,(const xmlChar*)"rcut", (const xmlChar*)"5");

    xmlNodePtr a= xmlNewTextChild(uuc,NULL,(const xmlChar*)"coefficients",(const xmlChar*)"\n        ");
    std::ostringstream o1;
    o1<< "uu"<<ionSpecies.speciesName[i];
    xmlNewProp(a,(const xmlChar*)"id", (const xmlChar*)o1.str().c_str());
    xmlNewProp(a,(const xmlChar*)"type", (const xmlChar*)"Array");
    xmlNewProp(a,(const xmlChar*)"optimize", (const xmlChar*)"yes");
    xmlAddChild(j3,uuc);

    xmlNodePtr udc = xmlNewNode(NULL,(const xmlChar*)"correlation");
    xmlNewProp(udc,(const xmlChar*)"ispecies", (const xmlChar*)ionSpecies.speciesName[i].c_str());
    xmlNewProp(udc,(const xmlChar*)"especies1", (const xmlChar*)"u");
    xmlNewProp(udc,(const xmlChar*)"especies2", (const xmlChar*)"d");
    xmlNewProp(udc,(const xmlChar*)"isize", (const xmlChar*)"3");
    xmlNewProp(udc,(const xmlChar*)"esize", (const xmlChar*)"3");
    xmlNewProp(udc,(const xmlChar*)"rcut", (const xmlChar*)"5");

    xmlNodePtr b= xmlNewTextChild(udc,NULL,(const xmlChar*)"coefficients",(const xmlChar*)"\n        ");
    std::ostringstream o2;
    o2<< "ud"<<ionSpecies.speciesName[i];
    xmlNewProp(b,(const xmlChar*)"id", (const xmlChar*)o2.str().c_str());
    xmlNewProp(b,(const xmlChar*)"type", (const xmlChar*)"Array");
    xmlNewProp(b,(const xmlChar*)"optimize", (const xmlChar*)"yes");
    xmlAddChild(j3,udc);

  }
  return j3;
}

xmlNodePtr QMCGaussianParserBase::createJ2()
{
  xmlNodePtr j2 = xmlNewNode(NULL,(const xmlChar*)"jastrow");
  xmlNewProp(j2,(const xmlChar*)"name", (const xmlChar*)"J2");
  xmlNewProp(j2,(const xmlChar*)"type", (const xmlChar*)"Two-Body");
  xmlNewProp(j2,(const xmlChar*)"function", (const xmlChar*)"Bspline");
  xmlNewProp(j2,(const xmlChar*)"print", (const xmlChar*)"yes");
  {
    xmlNodePtr uu = xmlNewNode(NULL,(const xmlChar*)"correlation");
    xmlNewProp(uu,(const xmlChar*)"rcut", (const xmlChar*)"10");
    xmlNewProp(uu,(const xmlChar*)"size", (const xmlChar*)"10");
    xmlNewProp(uu,(const xmlChar*)"speciesA", (const xmlChar*)"u");
    xmlNewProp(uu,(const xmlChar*)"speciesB", (const xmlChar*)"u");
    xmlNodePtr a= xmlNewTextChild(uu,NULL,(const xmlChar*)"coefficients",(const xmlChar*)"0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0");
    xmlNewProp(a,(const xmlChar*)"id", (const xmlChar*)"uu");
    xmlNewProp(a,(const xmlChar*)"type", (const xmlChar*)"Array");
    xmlAddChild(j2,uu);
  }
  {
    xmlNodePtr uu = xmlNewNode(NULL,(const xmlChar*)"correlation");
    xmlNewProp(uu,(const xmlChar*)"rcut", (const xmlChar*)"10");
    xmlNewProp(uu,(const xmlChar*)"size", (const xmlChar*)"10");
    xmlNewProp(uu,(const xmlChar*)"speciesA", (const xmlChar*)"u");
    xmlNewProp(uu,(const xmlChar*)"speciesB", (const xmlChar*)"d");
    //xmlNodePtr a = xmlNewNode(NULL,(const xmlChar*)"coefficients");
    xmlNodePtr a= xmlNewTextChild(uu,NULL,(const xmlChar*)"coefficients",(const xmlChar*)"0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0");
    xmlNewProp(a,(const xmlChar*)"id", (const xmlChar*)"ud");
    xmlNewProp(a,(const xmlChar*)"type", (const xmlChar*)"Array");
    // xmlAddChild(uu,a);
    xmlAddChild(j2,uu);
  }
  return j2;
}

xmlNodePtr QMCGaussianParserBase::createJ1()
{
  xmlNodePtr j1 = xmlNewNode(NULL,(const xmlChar*)"jastrow");
  xmlNewProp(j1,(const xmlChar*)"name", (const xmlChar*)"J1");
  xmlNewProp(j1,(const xmlChar*)"type", (const xmlChar*)"One-Body");
  xmlNewProp(j1,(const xmlChar*)"function", (const xmlChar*)"Bspline");
  xmlNewProp(j1,(const xmlChar*)"source", (const xmlChar*)"ion0");
  xmlNewProp(j1,(const xmlChar*)"print", (const xmlChar*)"yes");
  SpeciesSet& ionSpecies(IonSystem.getSpeciesSet());
  for(int i=0; i<ionSpecies.getTotalNum(); i++)
  {
    xmlNodePtr c = xmlNewNode(NULL,(const xmlChar*)"correlation");
    xmlNewProp(c,(const xmlChar*)"rcut", (const xmlChar*)"10");
    xmlNewProp(c,(const xmlChar*)"size", (const xmlChar*)"10");
    xmlNewProp(c,(const xmlChar*)"cusp", (const xmlChar*)"0");
    xmlNewProp(c,(const xmlChar*)"elementType", (const xmlChar*)ionSpecies.speciesName[i].c_str());
    xmlNodePtr a= xmlNewTextChild(c,NULL,(const xmlChar*)"coefficients",(const xmlChar*)"0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0");
    std::ostringstream o;
    o<< 'e'<<ionSpecies.speciesName[i];
    xmlNewProp(a,(const xmlChar*)"id", (const xmlChar*)o.str().c_str());
    xmlNewProp(a,(const xmlChar*)"type", (const xmlChar*)"Array");
    xmlAddChild(j1,c);
  }
  return j1;
}

void QMCGaussianParserBase::map2GridFunctors(xmlNodePtr cur)
{
  using namespace qmcplusplus;
  xmlNodePtr anchor = cur;
  //xmlNodePtr grid_ptr = 0;
  std::vector<xmlNodePtr> phi_ptr;
  std::vector<QuantumNumberType> nlms;
  int Lmax = 0;
  int current = 0;
  std::string acenter("none");
  const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"elementType");
  if(aptr)
    acenter = (const char*)aptr;
  xmlNodePtr grid_ptr=0;
  cur = anchor->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "grid")
      grid_ptr = cur;
    else
      if(cname == "basisGroup")
      {
        int n=1,l=0,m=0;
        const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"n");
        if(aptr)
          n = atoi((const char*)aptr);
        aptr = xmlGetProp(cur,(const xmlChar*)"l");
        if(aptr)
          l = atoi((const char*)aptr);
        Lmax = std::max(l,Lmax);
        phi_ptr.push_back(cur);
        nlms.push_back(QuantumNumberType());
        nlms[current][0]=n;
        nlms[current][1]=l;
        nlms[current][2]=m;
        ++current;
      }
    cur = cur->next;
  }
  if(grid_ptr == 0)
  {
    LOGMSG("Grid is not defined: using default")
    //xmlAddChild(anchor,gridPtr);
    grid_ptr = xmlCopyNode(gridPtr,1);
    xmlAddChild(anchor,grid_ptr);
  }
  RGFBuilderBase::CenteredOrbitalType aos(Lmax);
  bool normalized(Normalized=="yes");
  GTO2GridBuilder rbuilder(normalized);
  rbuilder.setOrbitalSet(&aos,acenter);
  rbuilder.addGrid(grid_ptr);
  for(int i=0; i<nlms.size(); i++)
  {
    rbuilder.addRadialOrbital(phi_ptr[i],nlms[i]);
  }
  rbuilder.print(acenter,1);
}

void QMCGaussianParserBase::createGridNode(int argc, char** argv)
{
  gridPtr = xmlNewNode(NULL,(const xmlChar*)"grid");
  std::string gridType("log");
  std::string gridFirst("1.e-6");
  std::string gridLast("1.e2");
  std::string gridSize("1001");
  int iargc=0;
  while(iargc<argc)
  {
    std::string a(argv[iargc]);
    if(a == "-gridtype")
    {
      gridType=argv[++iargc];
    }
    else
      if(a == "-frst")
      {
        gridFirst=argv[++iargc];
      }
      else
        if(a == "-last")
        {
          gridLast=argv[++iargc];
        }
        else
          if(a == "-size")
          {
            gridSize=argv[++iargc];
          }
          else
            if(a == "-nojastrow")
            {
              addJastrow=false;
            }
            else
              if(a == "-add3BodyJ")
              {
                addJastrow3Body=true;
              }
              else
                if(a == "-numMO")
                {
                  numMO2print = atoi(argv[++iargc]); 
                }
    ++iargc;
  }
  xmlNewProp(gridPtr,(const xmlChar*)"type",(const xmlChar*)gridType.c_str());
  xmlNewProp(gridPtr,(const xmlChar*)"ri",(const xmlChar*)gridFirst.c_str());
  xmlNewProp(gridPtr,(const xmlChar*)"rf",(const xmlChar*)gridLast.c_str());
  xmlNewProp(gridPtr,(const xmlChar*)"npts",(const xmlChar*)gridSize.c_str());
}

void QMCGaussianParserBase::dump(const std::string& psi_tag,
                                 const std::string& ion_tag)
{
  std::cout << " QMCGaussianParserBase::dump " << std::endl;
  {
    xmlDocPtr doc_p = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr qm_root_p = xmlNewNode(NULL, BAD_CAST "qmcsystem");
    xmlAddChild(qm_root_p,createElectronSet());
    xmlAddChild(qm_root_p,createIonSet());
    xmlDocSetRootElement(doc_p, qm_root_p);
    std::string fname = Title+"."+basisName+".ptcl.xml";
    xmlSaveFormatFile(fname.c_str(),doc_p,1);
    xmlFreeDoc(doc_p);
  }
  xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  {
    //wavefunction
    xmlNodePtr wfPtr = xmlNewNode(NULL,(const xmlChar*)"wavefunction");
    xmlNewProp(wfPtr,(const xmlChar*)"name",(const xmlChar*)psi_tag.c_str());
    xmlNewProp(wfPtr,(const xmlChar*)"target",(const xmlChar*)"e");
    {
      xmlNodePtr detPtr = xmlNewNode(NULL, (const xmlChar*) "determinantset");
      xmlNewProp(detPtr,(const xmlChar*)"type",(const xmlChar*)"MolecularOrbital");
      xmlNewProp(detPtr,(const xmlChar*)"name",(const xmlChar*)"LCAOBSet");
      xmlNewProp(detPtr,(const xmlChar*)"source",(const xmlChar*)ion_tag.c_str());
      xmlNewProp(detPtr,(const xmlChar*)"transform",(const xmlChar*)"yes");

      if (DoCusp==true)
         xmlNewProp(detPtr,(const xmlChar*)"cuspCorrection",(const xmlChar*)"yes");
      if(UseHDF5)
      {
         xmlNewProp(detPtr,(const xmlChar*)"href",(const xmlChar*)h5file.c_str());
      }

      {
        if(UseHDF5)
        {
          xmlNodePtr bsetPtr = createBasisSetWithHDF5();
        }
        else
        {
          xmlNodePtr bsetPtr = createBasisSet();
          xmlAddChild(detPtr,bsetPtr);
        }
        if(multideterminant)
        {
          xmlNodePtr spoupPtr = xmlNewNode(NULL,(const xmlChar*)"sposet");
          xmlNodePtr spodnPtr = xmlNewNode(NULL,(const xmlChar*)"sposet");
          xmlNewProp(spoupPtr,(const xmlChar*)"basisset",(const xmlChar*)"LCAOBSet");
          xmlNewProp(spodnPtr,(const xmlChar*)"basisset",(const xmlChar*)"LCAOBSet");
          createSPOSets(spoupPtr,spodnPtr);
          xmlAddChild(detPtr,spoupPtr);
          xmlAddChild(detPtr,spodnPtr);
          xmlNodePtr multislaterdetPtr=NULL;


          if (VSVB!=true && QP!=true) 
             multislaterdetPtr = createMultiDeterminantSet();
          else
          {
             if (VSVB)        
                multislaterdetPtr = createMultiDeterminantSetVSVB();
             else
                multislaterdetPtr = createMultiDeterminantSetQP();

          }
  
          xmlAddChild(detPtr,multislaterdetPtr);
        }
        else
        {
          xmlNodePtr slaterdetPtr=NULL;
          if(UseHDF5)
          {
            slaterdetPtr = createDeterminantSetWithHDF5();
          }
          else
          {
            slaterdetPtr = createDeterminantSet();
          }
          xmlAddChild(detPtr,slaterdetPtr);
        }
      }
      xmlAddChild(wfPtr,detPtr);
      if(addJastrow)
      {
        std::cout << "Adding Two-Body and One-Body jastrows with rcut=\"10\" and size=\"10\"" << std::endl;
        xmlAddChild(wfPtr,createJ2());
        xmlAddChild(wfPtr,createJ1());
      }
      if(addJastrow3Body)
      {
        std::cout << "Adding Three-Body rcut=\"4\". " << std::endl;
        xmlAddChild(wfPtr,createJ3());
      }
    }
    xmlAddChild(qm_root,wfPtr);
  }
  xmlDocSetRootElement(doc, qm_root);
  xmlXPathContextPtr m_context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr result
  = xmlXPathEvalExpression((const xmlChar*)"//atomicBasisSet",m_context);
  if(!xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    for(int ic=0; ic<result->nodesetval->nodeNr; ic++)
    {
      xmlNodePtr cur = result->nodesetval->nodeTab[ic];
      map2GridFunctors(cur);
    }
  }
  xmlXPathFreeObject(result);
  std::string fname = Title+"."+basisName+".xml";
  xmlSaveFormatFile(fname.c_str(),doc,1);
  xmlFreeDoc(doc);
  if (numMO*SizeOfBasisSet>=500 && !UseHDF5)
     std::cout<<"Consider using HDF5 via -hdf5 for higher performance and smaller wavefunction files"<<std::endl;
}
void QMCGaussianParserBase::Fmodump(const std::string& psi_tag,
                                 const std::string& ion_tag,
                                 std::string Mytag)
{
  std::cout << " QMCGaussianParserBase::Fmodump " << std::endl;
   if (DoCusp==true){
      // Creating a ptcl file with no esp for Cusp computation // 
      {
   
       xmlDocPtr doc_p = xmlNewDoc((const xmlChar*)"1.0");
       xmlNodePtr qm_root_p = xmlNewNode(NULL, BAD_CAST "qmcsystem");
       xmlAddChild(qm_root_p,createElectronSet());

     for(int i=0 ; i<NumberOfAtoms; i++)
       xmlAddChild(qm_root_p,createIonSet());
     for(int i=0 ; i<NumberOfAtoms; i++)
       xmlDocSetRootElement(doc_p, qm_root_p);
       std::string fname = Mytag+"-cusp.ptcl.xml";
       xmlSaveFormatFile(fname.c_str(),doc_p,1);
       xmlFreeDoc(doc_p);
      }
   }
   {
    xmlDocPtr doc_p = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr qm_root_p = xmlNewNode(NULL, BAD_CAST "qmcsystem");
    xmlAddChild(qm_root_p,createElectronSet());
    xmlAddChild(qm_root_p,createIonSet());
    for (int iesp=0;iesp<TotNumMonomer;iesp++)
    {
       if(FMO1)
          if (iesp!=FMOIndexI)
             xmlAddChild(qm_root_p,createESPSet(iesp));

       if(FMO2)
          if (iesp!=FMOIndexI && iesp!=FMOIndexJ)
             xmlAddChild(qm_root_p,createESPSet(iesp));

       if(FMO3)
          if (iesp!=FMOIndexI && iesp!=FMOIndexJ && iesp!=FMOIndexK)
             xmlAddChild(qm_root_p,createESPSet(iesp));


    }
    xmlDocSetRootElement(doc_p, qm_root_p);
    std::string fname = Mytag+".ptcl.xml";
    xmlSaveFormatFile(fname.c_str(),doc_p,1);
    xmlFreeDoc(doc_p);
  }
  xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  {
    //wavefunction
    xmlNodePtr wfPtr = xmlNewNode(NULL,(const xmlChar*)"wavefunction");
    xmlNewProp(wfPtr,(const xmlChar*)"id",(const xmlChar*)psi_tag.c_str());
    xmlNewProp(wfPtr,(const xmlChar*)"target",(const xmlChar*)"e");
    {
      xmlNodePtr detPtr = xmlNewNode(NULL, (const xmlChar*) "determinantset");
      xmlNewProp(detPtr,(const xmlChar*)"name",(const xmlChar*)"LCAOBSet");
      xmlNewProp(detPtr,(const xmlChar*)"type",(const xmlChar*)"MolecularOrbital");
      xmlNewProp(detPtr,(const xmlChar*)"transform",(const xmlChar*)"yes");
      xmlNewProp(detPtr,(const xmlChar*)"source",(const xmlChar*)ion_tag.c_str());
      if(DoCusp==true)
          xmlNewProp(detPtr,(const xmlChar*)"cuspCorrection",(const xmlChar*)"yes");
      {
        xmlNodePtr bsetPtr = createBasisSet();
        xmlAddChild(detPtr,bsetPtr);

        if(multideterminant)
        {
          xmlNodePtr spoupPtr = xmlNewNode(NULL,(const xmlChar*)"sposet");
          xmlNodePtr spodnPtr = xmlNewNode(NULL,(const xmlChar*)"sposet");
          xmlNewProp(spoupPtr,(const xmlChar*)"basisset",(const xmlChar*)"LCAOBSet");
          xmlNewProp(spodnPtr,(const xmlChar*)"basisset",(const xmlChar*)"LCAOBSet");
          createSPOSets(spoupPtr,spodnPtr);
          xmlAddChild(detPtr,spoupPtr);
          xmlAddChild(detPtr,spodnPtr);
          xmlNodePtr multislaterdetPtr=NULL;
          multislaterdetPtr = createMultiDeterminantSet();
          xmlAddChild(detPtr,multislaterdetPtr);
        }
        else
        {
          xmlNodePtr slaterdetPtr=NULL;
          if(UseHDF5)
          {
            slaterdetPtr = createDeterminantSetWithHDF5();
          }
          else
          {
            slaterdetPtr = createDeterminantSet();
          }
          xmlAddChild(detPtr,slaterdetPtr);
        }
      }
      xmlAddChild(wfPtr,detPtr);
      if(addJastrow)
      {
        std::cout << "Adding Two-Body and One-Body jastrows with rcut=\"10\" and size=\"10\"" << std::endl;
        xmlAddChild(wfPtr,createJ2());
        xmlAddChild(wfPtr,createJ1());
      }
      if(addJastrow3Body)
      {
        std::cout << "Adding Three-Body rcut=\"10\". " << std::endl;
        xmlAddChild(wfPtr,createJ3());
      }
    }
    xmlAddChild(qm_root,wfPtr);
  }

  xmlDocSetRootElement(doc, qm_root);
  xmlXPathContextPtr m_context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr result
  = xmlXPathEvalExpression((const xmlChar*)"//atomicBasisSet",m_context);
  if(!xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    for(int ic=0; ic<result->nodesetval->nodeNr; ic++)
    {
      xmlNodePtr cur = result->nodesetval->nodeTab[ic];
      map2GridFunctors(cur);

    }
  }
  xmlXPathFreeObject(result);
  std::string fname = Mytag+".wfs.xml";
  xmlSaveFormatFile(fname.c_str(),doc,1);
  xmlFreeDoc(doc);

  //Hamiltonian

  xmlDocPtr doc_h = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qm_root_h = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  std::string fname_h = Mytag+".ham.xml";

  xmlNodePtr hamPtr = xmlNewNode(NULL,(const xmlChar*)"hamiltonian");
  xmlNewProp(hamPtr,(const xmlChar*)"name",(const xmlChar*)"h0");
  xmlNewProp(hamPtr,(const xmlChar*)"type",(const xmlChar*)"generic");
  xmlNewProp(hamPtr,(const xmlChar*)"target",(const xmlChar*)"e");
  {
    xmlNodePtr pairpot1 = xmlNewNode(NULL,(const xmlChar*)"pairpot");
    xmlNewProp(pairpot1,(const xmlChar*)"name", (const xmlChar*)"ElecElec");
    xmlNewProp(pairpot1,(const xmlChar*)"type", (const xmlChar*)"coulomb");
    xmlNewProp(pairpot1,(const xmlChar*)"source", (const xmlChar*)"e");
    xmlNewProp(pairpot1,(const xmlChar*)"target", (const xmlChar*)"e");
    xmlAddChild(hamPtr,pairpot1);

    xmlNodePtr pairpot2 = xmlNewNode(NULL,(const xmlChar*)"pairpot");
    xmlNewProp(pairpot2,(const xmlChar*)"name", (const xmlChar*)"IonElec");
    xmlNewProp(pairpot2,(const xmlChar*)"type", (const xmlChar*)"coulomb");
    xmlNewProp(pairpot2,(const xmlChar*)"source", (const xmlChar*)ion_tag.c_str());
    xmlNewProp(pairpot2,(const xmlChar*)"target", (const xmlChar*)"e");
    xmlAddChild(hamPtr,pairpot2);

    xmlNodePtr pairpot3 = xmlNewNode(NULL,(const xmlChar*)"pairpot");
    xmlNewProp(pairpot3,(const xmlChar*)"name", (const xmlChar*)"IonIon");
    xmlNewProp(pairpot3,(const xmlChar*)"type", (const xmlChar*)"coulomb");
    xmlNewProp(pairpot3,(const xmlChar*)"source", (const xmlChar*)ion_tag.c_str());
    xmlNewProp(pairpot3,(const xmlChar*)"target", (const xmlChar*)ion_tag.c_str());
    xmlAddChild(hamPtr,pairpot3);

    for (int iesp=0;iesp<TotNumMonomer;iesp++)
    {
       std::string name;
       if(FMO1)
          if (iesp!=FMOIndexI)
          {
             name="nm"+ESPSystem[iesp].getName();
             xmlNodePtr pairpot4 = xmlNewNode(NULL,(const xmlChar*)"pairpot");
             xmlNewProp(pairpot4,(const xmlChar*)"name", (const xmlChar*)name.c_str());
             xmlNewProp(pairpot4,(const xmlChar*)"type", (const xmlChar*)"coulomb");
             xmlNewProp(pairpot4,(const xmlChar*)"source", (const xmlChar*)ESPSystem[iesp].getName().c_str());
             xmlNewProp(pairpot4,(const xmlChar*)"target", (const xmlChar*)"e");
             xmlAddChild(hamPtr,pairpot4);
          }
          
       if(FMO2)
          if (iesp!=FMOIndexI && iesp!=FMOIndexJ)
          {
             name="nm"+ESPSystem[iesp].getName();
             xmlNodePtr pairpot4 = xmlNewNode(NULL,(const xmlChar*)"pairpot");
             xmlNewProp(pairpot4,(const xmlChar*)"name", (const xmlChar*)name.c_str());
             xmlNewProp(pairpot4,(const xmlChar*)"type", (const xmlChar*)"coulomb");
             xmlNewProp(pairpot4,(const xmlChar*)"source", (const xmlChar*)ESPSystem[iesp].getName().c_str());
             xmlNewProp(pairpot4,(const xmlChar*)"target", (const xmlChar*)"e");
             xmlAddChild(hamPtr,pairpot4);
          }

       if(FMO3)
          if (iesp!=FMOIndexI && iesp!=FMOIndexJ && iesp!=FMOIndexK)
          {
             name="nm"+ESPSystem[iesp].getName();
             xmlNodePtr pairpot4 = xmlNewNode(NULL,(const xmlChar*)"pairpot");
             xmlNewProp(pairpot4,(const xmlChar*)"name", (const xmlChar*)name.c_str());
             xmlNewProp(pairpot4,(const xmlChar*)"type", (const xmlChar*)"coulomb");
             xmlNewProp(pairpot4,(const xmlChar*)"source", (const xmlChar*)ESPSystem[iesp].getName().c_str());
             xmlNewProp(pairpot4,(const xmlChar*)"target", (const xmlChar*)"e");
             xmlAddChild(hamPtr,pairpot4);
          }
    }

  }
  xmlAddChild(qm_root_h,hamPtr);

  xmlDocSetRootElement(doc_h, qm_root_h);
  xmlSaveFormatFile(fname_h.c_str(),doc_h,1);
  xmlFreeDoc(doc_h);

}
int QMCGaussianParserBase::numberOfExcitationsCSF( std::string& occ)
{
  int res=0;
  for(int i=ci_neb; i<ci_nstates; i++)
  {
    if(i < ci_nea && occ[i]=='2')
    {
      res++;  //excitation into singly occupied alpha states in the reference
    }
    else
    {
      if(occ[i] == '1' )
        res++;
      else
        if(occ[i] == '2' )
          res+=2;
    }
  }
  return res;
}

