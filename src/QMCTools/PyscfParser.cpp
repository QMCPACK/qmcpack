//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//
// File created by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////





#include "QMCTools/PyscfParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>



PyscfParser::PyscfParser()
{
  basisName = "Gaussian";
  Normalized = "no";
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  FixValence=true;
  NFZC=0;
}

PyscfParser::PyscfParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "Gaussian";
  Normalized = "no";
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  SpinRestricted=true;
  readtype=0;
  FixValence=true;
  NFZC=0;
}

void PyscfParser::parse(const std::string& fname)
{


  hdf_archive hin(0);

  if(!hin.open(fname.c_str(),H5F_ACC_RDONLY))
  {
       std::cerr<<"Could not open H5 file"<<std::endl;
       abort();
  }


  hin.push("parameters");

//  hin.read(usingECP,"ECP");
  hin.read(ECP,"ECP");
  //std::cout <<"usingECP: " <<(usingECP?("yes"):("no")) << std::endl;
  std::cout <<"usingECP: " <<(ECP?("yes"):("no")) << std::endl;
  std::cout.flush();

  multideterminant=false;
  std::cout <<"Multideterminant: " <<(multideterminant?("yes"):("no")) << std::endl;
  hin.read(SpinRestricted,"SpinResticted");

  if(SpinRestricted)
  { 
     hin.read(numAO,"numAO");
     hin.read(numMO,"numMO");

  }
  else
  {
     int numAO_up,numAO_dn,numMO_up,numMO_dn;

     hin.read(numAO_up,"numAO_up");
     hin.read(numMO_up,"numMO_up");
     hin.read(numAO_dn,"numAO_dn");
     hin.read(numMO_dn,"numMO_dn");
     if (numAO_up==numAO_dn)
           numAO=numAO_up;
     else
     {
           std::cout<<"numAO_up="<<numAO_up<<"  numAO_dn="<<numAO_dn<<std::endl;
           std::cerr<<"The number of AO for the up Orbitals are different than the number of AOs for down orbitals. This is probably an error in your Pyscf input. Please contact QMCPACK developers"<<std::endl;     
           abort();
     }
     if (numMO_up==numMO_dn)
           numMO=numMO_up;
     else
     {
           std::cout<<"numMO_up="<<numMO_up<<"  numMO_dn="<<numMO_dn<<std::endl;
           std::cerr<<"The number of MO for the up Orbitals are different than the number of MOs for down orbitals. This is probably an error in your Pyscf input. Please contact QMCPACK developers"<<std::endl;     
           abort();
     }
  }



  std::cout <<"NUMBER OF AOs: " <<numAO << std::endl;
  SizeOfBasisSet = numAO; 
  std::cout <<"Size of Basis Set: " <<SizeOfBasisSet << std::endl;
  std::cout <<"NUMBER OF MOs: " <<numMO << std::endl;

  bool bohr=true;
  hin.read(bohr,"Unit");
  std::cout<<"Unit in Bohr="<<bohr<<std::endl;
  BohrUnit=bohr; 

  hin.read(NumberOfAlpha,"NbAlpha");
  hin.read(NumberOfBeta,"NbBeta");
  hin.read(NumberOfEls,"NbTotElec");
  hin.read(SpinMultiplicity,"spin");

  std::cout <<"Number of alpha electrons: " <<NumberOfAlpha << std::endl;
  std::cout <<"Number of beta electrons: " <<NumberOfBeta << std::endl;
  std::cout <<"Number of electrons: " <<NumberOfEls << std::endl;
  std::cout <<"SPIN MULTIPLICITY: " <<SpinMultiplicity << std::endl;

  hin.pop();
  hin.push("atoms");
  hin.read(NumberOfAtoms,"number_of_atoms");
  std::cout <<"NUMBER OF ATOMS: " <<NumberOfAtoms << std::endl;
  hin.close();


  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  getGeometry(fname);

}

void PyscfParser::getGeometry(const std::string& fname)
{

  hdf_archive hin(0);

  if(!hin.open(fname.c_str(),H5F_ACC_RDONLY))
  {
       std::cerr<<"Could not open H5 file"<<std::endl;
       abort();
  }
  hin.push("atoms");

  //atomic numbers
  std::vector<int> atomic_number;
  std::vector<double> q,pos;
  tags.clear();
  //read atomic info
  //MAP: (Atom Number, Species Index)
  std::vector<int> idx(NumberOfAtoms);
  std::map<int,int> AtomIndexmap;
  hin.read(idx,"species_ids");
  for(int i=0;  i<NumberOfAtoms; i++)
  {
     AtomIndexmap.insert(std::pair<int,int>(i,idx[i]) );  
  }
  for(int i=0;  i<NumberOfAtoms; i++)
  {
    std::string speciesName("species_");
    speciesName=speciesName + std::to_string(AtomIndexmap[i]);
    hin.push(speciesName);
    int zint,mycore;
    double z;
    std::string Name;
    hin.read(zint,"atomic_number"); 
    atomic_number.push_back(zint);
    hin.read(z,"charge"); 
    q.push_back(z); 
    hin.read(Name,"name");
    tags.push_back(Name);
    hin.pop();
  }

  Matrix<double> IonPos(NumberOfAtoms,3);
  hin.read(IonPos,"positions");

  SpeciesSet& species(IonSystem.getSpeciesSet());
  for(int i=0; i<NumberOfAtoms; i++)
  {
    IonSystem.R[i][0]=IonPos[i][0];
    IonSystem.R[i][1]=IonPos[i][1];
    IonSystem.R[i][2]=IonPos[i][2];
    GroupName[i]=IonName[atomic_number[i]];
    int speciesID = species.addSpecies(GroupName[i]);
    IonSystem.GroupID[i]=speciesID;
    species(AtomicNumberIndex,speciesID)=atomic_number[i];
    species(IonChargeIndex,speciesID)=q[i];
   
  }
}

