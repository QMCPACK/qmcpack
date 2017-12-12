//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "QMCTools/CasinoParser.h"
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>


CasinoParser::CasinoParser()
{
  basisName = "casino-G2";
  Normalized = "yes";
}

CasinoParser::CasinoParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "casino-G2";
  Normalized = "yes";
}

void CasinoParser::parse(const std::string& fname)
{
  if(multideterminant)
  {
    std::cerr <<"Multideterminant parser for Casino is not implemented. \n";
    exit(301);
  }
  std::ifstream fin(fname.c_str());
  //Grep the first word of the first line to assign a file name
  fin.getline(dbuffer,sizeof(dbuffer));
  std::istringstream a(dbuffer);
  a>>Title;
  LOGMSG("Looking for Periodicity ")
  search(fin,"Periodicity");
  int periodicity;
  getValue(fin,periodicity);
  if(periodicity > 0)
    Periodicity=true;
  std::string spin_unrestricted;
  LOGMSG("Looking for Spin ")
  search(fin,"Spin");
  getValue(fin,spin_unrestricted);
  if(spin_unrestricted.find("false")<spin_unrestricted.size())
  {
    SpinRestricted=true;
  }
  else
  {
    SpinRestricted=false;
  }
  LOGMSG("Looking for electrons ")
  search(fin,"electrons");
  getValue(fin,NumberOfEls);
  LOGMSG("Looking for GEOMETRY ")
  search(fin, "GEOMETRY");
  getGeometry(fin);
  SpeciesSet& ionSpecies(IonSystem.getSpeciesSet());
  //GroupID is the atomic number. Have to match the name and the real GroupID
  for(int i=0; i<NumberOfAtoms; i++)
  {
    if(IonSystem.GroupID[i]>200)
      IonSystem.GroupID[i]-=200;
    int atomic_number = IonSystem.GroupID[i];
    GroupName[i]=IonName[IonSystem.GroupID[i]];
    int gid=IonSystem.GroupID[i]=IonSystem.getSpeciesSet().addSpecies(GroupName[i]);
    ionSpecies(AtomicNumberIndex,gid)=atomic_number;
    ionSpecies(ValenceChargeIndex,gid)=Qv[i];
  }
  LOGMSG("Looking for BASIS ")
  search(fin, "BASIS");
  getGaussianCenters(fin);
  //search(fin, "MULTIDETERMINANT");
  EigVal_alpha.resize(SizeOfBasisSet);
  EigVal_beta.resize(SizeOfBasisSet);
  std::vector<value_type> etemp;
// morales: setting numMO to SizeOfBasisSet, for now
  numMO = SizeOfBasisSet;
  //////////////
  LOGMSG("Looking for EIGENVECTOR ")
  search(fin, "EIGENVECTOR");
  int nstates=SizeOfBasisSet;
  if(SpinRestricted)
  {
    EigVec.resize(SizeOfBasisSet*SizeOfBasisSet);
    etemp.resize(SizeOfBasisSet);
  }
  else
  {
    EigVec.resize(2*SizeOfBasisSet*SizeOfBasisSet);
    etemp.resize(2*SizeOfBasisSet);
  }
  getValues(fin,EigVec.begin(), EigVec.end());
  std::string aline;
  getline(fin,aline,'\n');
  getline(fin,aline,'\n');
  if(aline.find("EIGENVALUES")<aline.size())
  {
    LOGMSG("Looking for EIGENVALUES ")
    skiplines(fin,2);
    getValues(fin,etemp.begin(), etemp.end());
    copy(etemp.begin(), etemp.begin()+SizeOfBasisSet, EigVal_alpha.begin());
    if(SpinRestricted)
    {
      EigVal_beta = EigVal_alpha;
    }
    else
    {
      copy(etemp.begin()+SizeOfBasisSet, etemp.end(), EigVal_beta.begin());
    }
  }
  else
  {
    WARNMSG("Missing EIGENVALUES. Modify determinant/@orbitals")
  }
  //Need to make correction to the eigen vectors
  int tot=0;
  for(int i=0; i<nstates; i++)
  {
    for(int j=0; j<SizeOfBasisSet; j++, tot++)
      EigVec[tot]*=BasisCorrection[j];
  }
}

void CasinoParser::getGeometry(std::istream& is)
{
  //Number of atoms
  getNumberOfAtoms(is);
  search(is,"Atomic positions");
  getValues(is,IonSystem.R.begin(),IonSystem.R.end());
  search(is,"Atomic numbers");
  getValues(is,IonSystem.GroupID.begin(),IonSystem.GroupID.end());
  search(is,"Valence charges");
  getValues(is,Qv.begin(),Qv.end());
  if(Periodicity)
  {
    std::vector<double> lat(9);
    search(is,"Primitive lattice");
    getValues(is,lat.begin(),lat.end());
    int ij=0;
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++,ij++)
        IonSystem.Lattice.R(i,j)=lat[ij];
    IonSystem.Lattice.reset();
    std::cout << "Lattice vectors " << std::endl;
    IonSystem.Lattice.print(std::cout);
  }
}

void CasinoParser::getNumberOfAtoms(std::istream& is)
{
  search(is,"Number");
  getValue(is,NumberOfAtoms);
  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  Qv.resize(NumberOfAtoms);
  gBound.resize(NumberOfAtoms+1);
}
void CasinoParser::getAtomicPositions(std::istream& is)
{
}
void CasinoParser::getAtomicNumbers(std::istream& is)
{
}
void CasinoParser::getValenceCharges(std::istream& is)
{
}

void CasinoParser::getGaussianCenters(std::istream& is)
{
  int n=0;
  std::streampos pivot= is.tellg();
  search(is, "Number of shells");
  getValue(is,n);
  gShell.resize(n);
  gNumber.resize(n);
  LOGMSG("Number of shells " << n)
  is.seekg(pivot);//rewind it
  search(is, "Number of basis functions");
  getValue(is,SizeOfBasisSet);
  LOGMSG("Number of basis functions " << SizeOfBasisSet);
  is.seekg(pivot);//rewind it
  search(is, "Number of Gaussian primitives");
  getValue(is,n);
  gExp.resize(n);
  gC0.resize(n);
  gC1.resize(n);
  LOGMSG("Number of Gaussian primitives " << n)
  is.seekg(pivot);//rewind it
  search(is, "Code for shell types");
  getValues(is,gShell.begin(), gShell.end());
  LOGMSG("Checking shell types")
  copy(gShell.begin(), gShell.end(),std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
  //becomes true, if there is sp shell
  bool SPshell(false);
  for(int ig=0; ig<gShell.size(); ig++)
  {
    if(gShell[ig] == 2)
      SPshell=true;
  }
  LOGMSG("Checking the number of gaussians per shell")
  is.seekg(pivot);//rewind it
  search(is, "Number of primitive Gaussians");
  getValues(is,gNumber.begin(), gNumber.end());
  LOGMSG("Checking the bound of shells")
  is.seekg(pivot);//rewind it
  search(is, "Sequence number");
  getValues(is,gBound.begin(), gBound.end());
  for(int i=0; i<gBound.size(); i++)
    gBound[i]-= 1;
  is.seekg(pivot);//rewind it
  LOGMSG("Checking the gaussian exponents")
  search(is, "Exponents");
  getValues(is,gExp.begin(), gExp.end());
  LOGMSG("Checking the gaussian contractions")
  search(is, "contraction");
  getValues(is,gC0.begin(), gC0.end());
  if(SPshell)
    //only check if there is 2nd
  {
    LOGMSG("Checking the gaussian contractions (sp)")
    search(is, "2nd");
    getValues(is,gC1.begin(), gC1.end());
  }
  else
  {
    LOGMSG("Input file does not have sp shell")
  }
  makeCorrections();
}

double CasinoParser::contractionCorrection(int shell_id, double alpha)
{
  const double pi = 4.0*atan(1.0);
  const double sqrtpi=sqrt(pi);
  double fac;
  switch(shell_id)
  {
  case(1): //s, 1/(2 sqrt(pi))
    fac = 2.0*sqrtpi;
    break;
  case(2): // sp
    fac = 2.0*sqrtpi;
    break;
  case(3): // p
    fac = sqrt(4.0/3.0)*sqrtpi;
    break;
  case(4): // d
    fac = sqrt(16.0/15.0)*sqrtpi;
    break;
  case(5): // f
    fac = 1.0e0;
    //fac *= pow(2.0e0*alpha,2)/sqrt(pi); break;
  default: // others, return 1 for now
    fac = 1.0e0;
    break;
  }
  return fac;
}

/** make corrections to gC0 and gC1 and tabulate the corrections to the eigen vectors
 */
void CasinoParser::makeCorrections()
{
  int n=gC0.size();
  for(int i=0; i<n; i++)
  {
    gC0[i] *= contractionCorrection(gShell[i],gExp[i]);
    if(gShell[i] == 2)
      gC1[i] *= contractionCorrection(3,gExp[i]);
  }
  //s, sp and p do not need any correction
  BasisCorrection.resize(SizeOfBasisSet,1.0);
  std::vector<int> offset(10,0);
  offset[1]=1; //s
  offset[2]=4; //sp
  offset[3]=3; //p
  offset[4]=5; //d
  offset[5]=7; //f
  const double m40=sqrt(3.0e0);
  const double m50=sqrt(15.0e0/8.0e0);
  const double m51=1.5e0*sqrt(5.0e0);
  const double m52=15.0e0/sqrt(2.0e0);
  const double m53=15.0e0*sqrt(3.0e0);
  int basisCount=0;
  for(int i=0; i<gShell.size(); i++)
  {
    if(gShell[i] == 4)
      //d-orbital corrections
    {
      BasisCorrection[basisCount++]=m40;//m=0
      BasisCorrection[basisCount++]=0.5;//m=1
      BasisCorrection[basisCount++]=0.5;//m=-1
      BasisCorrection[basisCount++]=1.0;//m=2
      BasisCorrection[basisCount++]=0.5;//m=-2
    }
    else
      if(gShell[i] == 5)
        //f-orbital corrections
      {
        BasisCorrection[basisCount++]=m50;//m=0
        BasisCorrection[basisCount++]=m51;//m=1
        BasisCorrection[basisCount++]=m51;//m=-1
        BasisCorrection[basisCount++]=m52;//m=2
        BasisCorrection[basisCount++]=m52;//m=-2
        BasisCorrection[basisCount++]=m53;//m=3
        BasisCorrection[basisCount++]=m53;//m=-3
      }
      else
        //increase the count
      {
        basisCount += offset[gShell[i]];
      }
  }
}
