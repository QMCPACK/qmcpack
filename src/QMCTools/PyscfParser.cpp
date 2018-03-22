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
#include <sstream>


char *binpad (unsigned long int  n, size_t sz);

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


  hin.push("PBC");
  hin.read(PBC,"PBC");
  std::cout <<"Periodic Boundary Comditions: " <<(PBC?("yes"):("no")) << std::endl;
  std::cout.flush();

  hin.pop();

  hin.push("parameters");

  hin.read(ECP,"ECP");
  std::cout <<"usingECP: " <<(ECP?("yes"):("no")) << std::endl;
  std::cout.flush();

//  multideterminant=false;
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
  if (PBC){
     getCell(fname);
     getKpts(fname);
     if (debug){
         getGaussianCenters(fname);
         getMO(fname);
 
     }
   
  }
  getGeometry(fname);

  if(multideterminant)
  {
     outputFile;
     hdf_archive hin(0);

     if(!hin.open(outputFile.c_str(),H5F_ACC_RDONLY))
     {
          std::cerr<<"Could not open H5 file"<<std::endl;
          abort();
     }
     
     if(!hin.push("MultiDet"))
     {
          std::cerr<<"Could not open H5 file"<<std::endl;
          abort();
     }
     else
     {
           hin.read(ci_size,"NbDet");
           CIcoeff.clear();
           CIalpha.clear();
           CIbeta.clear();
           CIcoeff.resize(ci_size);
           CIalpha.resize(ci_size);
           CIbeta.resize(ci_size);


           hin.read(ci_nstates,"nstate");
           int N_int;
           const int bit_kind = 64; 
           hin.read(N_int,"Nbits");


           Matrix<unsigned long int> tempAlpha(ci_size,N_int);
           hin.read(tempAlpha,"CI_Alpha");

           Matrix<unsigned long int> tempBeta(ci_size,N_int);
           hin.read(tempBeta,"CI_Beta");
            
           std::string MyCItempAlpha,MyCItempBeta;
           MyCItempAlpha.resize(ci_nstates+1);
           MyCItempBeta.resize(ci_nstates+1);

           for (int ni=0; ni<ci_size;ni++)
           {
                int j=0;
                int jj=0;
                for (int k=0; k<N_int; k++)
                {
      	 	    for(int i=0; i<bit_kind;i++) 
                    {
                        if ( j <ci_nstates ) 
                        {
                     	   MyCItempAlpha[j] = binpad(tempAlpha[ni][k],bit_kind)[i]; 
             		   j++;
        		}
                        if ( jj <ci_nstates ) 
                        {
                     	   MyCItempBeta[jj] = binpad(tempBeta[ni][k],bit_kind)[i]; 
             		   jj++;
        		}
                     }
                }
                CIalpha[ni]=MyCItempAlpha;
                CIbeta[ni]=MyCItempBeta;
           }

           hin.read(CIcoeff,"Coeff");

           int ds=SpinMultiplicity-1;
           int neb= (NumberOfEls-ds)/2;
           int nea= NumberOfEls-NumberOfBeta;
           ci_nea= NumberOfAlpha;
           ci_neb= NumberOfBeta;           
           ci_nca = nea-ci_nea;
           ci_ncb = neb-ci_neb;
           std::cout <<" Done reading CIs!!"<< std::endl;
           hin.close(); 
     }     
  }  
}

void PyscfParser::getCell(const std::string& fname)
{
  X.resize(3);
  Y.resize(3);
  Z.resize(3);

  hdf_archive hin(0);

  if(!hin.open(fname.c_str(),H5F_ACC_RDONLY))
  {
       std::cerr<<"Could not open H5 file"<<std::endl;
       abort();
  }
  hin.push("Cell");
  Matrix<double> LatticeVec(3,3);                                                                                                       
  hin.read(LatticeVec,"LatticeVectors"); 
  X[0]=LatticeVec[0][0];
  X[1]=LatticeVec[0][1];
  X[2]=LatticeVec[0][2];
  Y[0]=LatticeVec[1][0];
  Y[1]=LatticeVec[1][1];
  Y[2]=LatticeVec[1][2];
  Z[0]=LatticeVec[2][0];
  Z[1]=LatticeVec[2][1];
  Z[2]=LatticeVec[2][2];
  hin.close();
  std::cout<<"Lattice parameters in Bohr:"<<std::endl;
  std::cout<<X[0]<<"  "<<X[1]<<"  "<<X[2]<<std::endl;
  std::cout<<Y[0]<<"  "<<Y[1]<<"  "<<Y[2]<<std::endl;
  std::cout<<Z[0]<<"  "<<Z[1]<<"  "<<Z[2]<<std::endl;
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
  hin.close();
}

void PyscfParser::getKpts(const std::string& fname)
{
  Matrix <double> MyVec(1,3);
  hdf_archive hin(0);

  if(!hin.open(fname.c_str(),H5F_ACC_RDONLY))
  {
       std::cerr<<"Could not open H5 file"<<std::endl;
       abort();
  }

  hin.push("Nb_KPTS");
  hin.read(NbKpts,"Nbkpts");
  hin.pop();
  Kpoints_Coord.resize(NbKpts);

  for (int i=0; i<NbKpts;i++)
  {
     
     Kpoints_Coord[i].resize(3);
     std::stringstream ss;
     ss << "KPTS_" << i;
     if(!hin.push(ss.str()))
     {
          std::cerr<<"Could not find coordinates for Kpoint Nb"<<i<<std::endl;
          abort();
     }
     hin.read(MyVec,"Coord");
     hin.pop(); 
     std::cout<<"MyCoord="<<MyVec[0][0]<<"  "<<MyVec[0][1]<<"  "<<MyVec[0][2]<<std::endl;
     Kpoints_Coord[i][0]=MyVec[0][0];
     Kpoints_Coord[i][1]=MyVec[0][1];
     Kpoints_Coord[i][2]=MyVec[0][2];
      
  }
  hin.close();



}



char *binpad (unsigned long int n, size_t sz)
{
    static char s[64 + 1] = {0};
    char *p = s + 64 ;

    for (int i = sz ; i > 0; i--)
        *--p = ( (n>> (i-1) ) & 1) ? '1' : '0';
    
    return p;
}



void PyscfParser::getMO(const std::string & fname)
{
  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  EigVec.resize(2*SizeOfBasisSet*numMO);
  
  std::string setname;
  Matrix<double> CartMat(SizeOfBasisSet,SizeOfBasisSet);

  hdf_archive hin(0);

  if(!hin.open(fname.c_str(),H5F_ACC_RDONLY))
  {
       std::cerr<<"Could not open H5 file"<<std::endl;
       abort();
  }
    char name[72];
    sprintf(name,"%s","/KPTS_0/eigenset_0");
    setname=name;
    if(!hin.read(CartMat,setname))
    {
       setname="SPOSetBase::putFromH5 Missing "+setname+" from HDF5 File.";
       APP_ABORT(setname.c_str());
    }
    hin.close();
  //hin.push("KPTS_0");
  //hin.read(CartMat,"eingenset_0"); 
  int cnt=0;
  for(int i=0; i<numMO; i++)
    for(int k=0; k<SizeOfBasisSet; k++)
      EigVec[cnt++] = CartMat[i][k];

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
  std::cout<<eig.str().c_str()<<std::endl;         
  std::cout <<"Finished reading MO." << std::endl;
  hin.close();
}

void PyscfParser::getGaussianCenters(const std::string fname)
{
}
