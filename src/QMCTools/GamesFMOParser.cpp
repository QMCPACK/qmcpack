//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////





#include "QMCTools/GamesFMOParser.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>



GamesFMOParser::GamesFMOParser()
{
  basisName = "Gaussian-G2";
  Normalized = "no";
   Mono=false;
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  psi_tag="psi0";
  ion_tag="ion0";
  SpinRestricted=true;
  TotNumAtom=0;
  NumMonomer=0;
  NumDimer=0;
  NumTrimer=0;
  FMOMethod=0;
  FMO=true;
  MonomerID=0;
  FixValence=true;
}

GamesFMOParser::GamesFMOParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "Gaussian-G2";
  Normalized = "no";
  Mono=false;
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  psi_tag="psi0";
  ion_tag="ion0";
  SpinRestricted=true;
  TotNumAtom=0;
  NumMonomer=0;
  NumDimer=0;
  NumTrimer=0;
  FMOMethod=0;
  FMO=true;
  MonomerID=0;
  FixValence=true;
}


void GamesFMOParser::parse(const std::string& fname)
{

  std::streampos pivot;
  std::ifstream fin(fname.c_str());
  pivot_begin= fin.tellg();
  std::string aline;
  int count;

  IDmer *IDDimer;
  IDmer *IDTrimer;



  fin.seekg(pivot_begin);


  search(fin,"NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS",aline);
  parsewords(aline.c_str(),currentWords);
  SizeOfBasisSet = atoi(currentWords[6].c_str());

  search(fin,"TOTAL NUMBER OF ATOMS",aline);
  parsewords(aline.c_str(),currentWords);
  TotNumAtom = atoi(currentWords[4].c_str());

  search(fin,"NUMBER OF MONOMERS",aline);
  parsewords(aline.c_str(),currentWords);
  NumMonomer = atoi(currentWords[3].c_str());

  search(fin,"NUMBER OF DIMERS",aline);
  parsewords(aline.c_str(),currentWords);
  NumDimer = atoi(currentWords[3].c_str());

  search(fin,"NUMBER OF TRIMERS",aline);
  parsewords(aline.c_str(),currentWords);
  NumTrimer = atoi(currentWords[3].c_str());

  search(fin,"FMO METHOD",aline);
  parsewords(aline.c_str(),currentWords);
  FMOMethod = atoi(currentWords[2].c_str());


  if(FMOMethod==3 && NumTrimer==0){
     std::cerr <<" Problem in the Input file! FMO level is 3 but the number of Trimers is 0"<< std::endl;
     abort();
  }
  std::cout <<" Tot Atoms="<<TotNumAtom<<" NumMonomer="<<NumMonomer<<"  NumDimer="<<NumDimer<<"  NumTrimer="<<NumTrimer<<"  FMOMethod="<<FMOMethod<< std::endl;

  MOtype = "Canonical";
  readtype=0;
  TotNumMonomer=NumMonomer;




  IDMonomer = new std::string[NumMonomer];
  IDDimer = new IDmer[NumDimer];
  ESPSystem = new ParticleSet [NumMonomer];
  ESPIonChargeIndex = new int [NumMonomer];
  ESPValenceChargeIndex = new int [NumMonomer];
  ESPAtomicNumberIndex = new int [NumMonomer];
  ESPGroupName.resize(NumMonomer);
  ESP.resize(NumMonomer);


  for (int iesp=0;iesp<NumMonomer;iesp++){
     ESPIonChargeIndex[iesp]=ESPSystem[iesp].getSpeciesSet().addAttribute("charge");
     ESPValenceChargeIndex[iesp]=ESPSystem[iesp].getSpeciesSet().addAttribute("valence");
     ESPAtomicNumberIndex[iesp]=ESPSystem[iesp].getSpeciesSet().addAttribute("atomicnumber");
  }


  pivot= fin.tellg();
  for(int i=0;i<NumMonomer;i++)
  {
     fin.seekg(pivot);
     std::ostringstream convert;
     convert <<i+1;
     IDMonomer[i]="CURRENT N-MER COORDINATES, I=      "+convert.str()+" J=      0 K=      0";
     std::cout <<IDMonomer[i]<< std::endl;
     getAllMonomer(fin,i);
  }


  fin.seekg(pivot);
  count=0;
  for(int i=2; i<=NumMonomer; i++){
     for(int j=1; j<=i-1;j++){
        std::ostringstream convertI,convertJ;
        convertI <<i;
        convertJ <<j;
        IDDimer[count].MyId="CURRENT N-MER COORDINATES, I=      "+convertI.str()+" J=      "+convertJ.str()+" K=      0";
        IDDimer[count].IndexI=i-1;
        IDDimer[count].IndexJ=j-1;
        IDDimer[count].IndexK=0;
        std::cout <<IDDimer[count].MyId<< std::endl;
        count++;

     }
  }

  if(FMOMethod==3){

    fin.seekg(pivot);
    count=0;
    IDTrimer = new IDmer[NumTrimer];
    for(int i=3; i<=NumMonomer; i++){
       for(int j=2; j<=i-1;j++){
          for(int k=1; k<=j-1;k++){
             std::ostringstream convertI,convertJ,convertK;
             convertI <<i;
             convertJ <<j;
             convertK <<k;
             IDTrimer[count].MyId="CURRENT N-MER COORDINATES, I=      "+convertI.str()+" J=      "+convertJ.str()+" K=      "+convertK.str();
             std::cout <<IDTrimer[count].MyId<< std::endl;

             IDTrimer[count].IndexI=i-1;
             IDTrimer[count].IndexJ=j-1;
             IDTrimer[count].IndexK=k-1;
             count++;


           }
       }
    }
  }



  Tot_Monomer = new Monomer[NumMonomer];

  for(int Mer=0;Mer<FMOMethod;Mer++)
  {
     int This=0;
     if (Mer==0){
        This=NumMonomer;
        std::cout <<"Running Monomers"<< std::endl;
        FMO1=true;
        FMO2=false;
        FMO3=false;
     }

     if (Mer==1){
        This=NumDimer;
        std::cout <<"Running Dimers"<< std::endl;
        FMO1=false;
        FMO2=true;
        FMO3=false;
     }

     if (Mer==2){
        This=NumTrimer;
        std::cout <<"Running Trimers"<< std::endl;
        FMO1=false;
        FMO2=false;
        FMO3=true;
     }

     for (int i=0; i<This; i++)
     {
         std::streampos pivot;
         std::ifstream fin(fname.c_str());
         pivot_begin= fin.tellg();
         std::cout <<"Working On "<<Mer<<"-"<<i<< std::endl;
         std::string temp_tag,Frag_Tag;


         if (Mer==0){
           std::ostringstream convert;
           convert <<i+1;
           Frag_Tag="Monomer"+convert.str();
           Mono=true;
           temp_tag=IDMonomer[i];

           FMOIndexI=i;
           FMOIndexJ=FMOIndexK=0;
           MonomerID++;
         }

         if (Mer==1){
           std::ostringstream convert;
           convert <<i+1;
           Frag_Tag="Dimer"+convert.str();
           temp_tag=IDDimer[i].MyId;
           FMOIndexI=IDDimer[i].IndexI;
           FMOIndexJ=IDDimer[i].IndexJ;
           FMOIndexK=IDDimer[i].IndexK;
         }

         if(FMOMethod==3){
            if (Mer==2){
              std::ostringstream convert;
              convert <<i+1;
              Frag_Tag="Trimer"+convert.str();
              temp_tag=IDTrimer[i].MyId;
              FMOIndexI=IDTrimer[i].IndexI;
              FMOIndexJ=IDTrimer[i].IndexJ;
              FMOIndexK=IDTrimer[i].IndexK;
            }
         }


         if(lookFor(fin,temp_tag,aline))
         {
            std::cout <<"temp_tag="<<temp_tag<< std::endl;
            search(fin,"NUMBER OF ATOMS IN FRAGMENT",aline);
            parsewords(aline.c_str(),currentWords);
            NumberOfAtoms= atoi(currentWords[5].c_str());
            std::cout <<"Total Number of Atoms="<<NumberOfAtoms<< std::endl;

            search(fin,"NUMBER OF CARTESIAN ATOMIC ORBITALS",aline);
            parsewords(aline.c_str(),currentWords);
            SizeOfBasisSet= atoi(currentWords[5].c_str());
            std::cout <<"Number of Cartesian atomic orbitals="<<SizeOfBasisSet<< std::endl;


            search(fin,"SPIN MULTIPLICITY",aline);
            parsewords(aline.c_str(),currentWords);
            SpinMultiplicity= atoi(currentWords[2].c_str());
            std::cout <<"SpinMultiplicity="<<SpinMultiplicity<< std::endl;
            search(fin,"NUMBER OF ELECTRONS ",aline);
            parsewords(aline.c_str(),currentWords);
            NumberOfEls= atoi(currentWords[3].c_str());
            std::cout <<"Number of electrons="<<NumberOfEls<< std::endl;
            search(fin,"NUMBER OF OCCUPIED ORBITALS (ALPHA)",aline);
            parsewords(aline.c_str(),currentWords);
            NumberOfAlpha = atoi(currentWords[5].c_str());
            std::cout <<"Number of alpha electrons: " <<NumberOfAlpha << std::endl;
            search(fin,"NUMBER OF OCCUPIED ORBITALS (BETA)",aline);
            parsewords(aline.c_str(),currentWords);
            NumberOfBeta = atoi(currentWords[5].c_str());
            std::cout <<"Number of beta electrons: " <<NumberOfBeta << std::endl;



            IonSystem.setName(ion_tag);
            IonSystem.create(NumberOfAtoms);

            GroupName.resize(NumberOfAtoms);
            getGeometry(fin);

            search(fin,"TOTAL NUMBER OF MOS IN VARIATION SPACE",aline);
            parsewords(aline.c_str(),currentWords);
            numMO= atoi(currentWords[7].c_str());
            std::cout <<"Number of Molecular Orbitals="<<numMO<< std::endl;


            fin.seekg(pivot_begin);
            getGaussianCenters(fin);
            fin.seekg(pivot_begin);
            search(fin,temp_tag,aline);
            //Looking for Eigenvectors
            lookFor(fin,"ATOM CHARGE        X                 Y                 Z (ANGST)   ESP CHARGE");

            getMO(fin,temp_tag);
            Fmodump(psi_tag, ion_tag,Frag_Tag);


            IonSystem.clear();
            GroupName.clear();
            fin.close();
         }
         else
         {
            std::cerr <<"Could not find all atoms in FMO-Games output"<< std::endl;
            std::cerr <<"	  Problem in "<<temp_tag<< std::endl;
            abort();
         }
      }
   }
   fin.close();
   for(int i=0;i<NumMonomer;i++){
      ESPSystem[i].clear();
   }
   ESPGroupName.clear();
   if (FMOMethod==3)
      delete [] IDTrimer;

   delete [] IDMonomer;
   delete [] IDDimer;
   delete [] Tot_Monomer;
   delete [] ESPIonChargeIndex;
   delete [] ESPValenceChargeIndex;
   delete [] ESPAtomicNumberIndex;
   delete [] ESPSystem;
   OHMMS::Controller->finalize();

}
void GamesFMOParser::getAllMonomer(std::istream& is, int Index)
{
  const double ang_to_bohr=1.0/0.529177e0;
  std::string aline;
  //atomic numbers
  std::vector<int> atomic_number,core;
  std::vector<double> q,pos;
  int natms=0;
  int MyAtoms=0;
  //tags.clear();
  std::string temp;


  if(search(is,IDMonomer[Index],aline))
  {
    getwords(currentWords,is);
    MyAtoms= atoi(currentWords[5].c_str());

    std::ostringstream temp;
    std::string ESP_ion_tag;
    temp<<Index;
    ESP_ion_tag="ESP_Charge"+temp.str();
    ESPSystem[Index].setName(ESP_ion_tag);
    ESPSystem[Index].create(MyAtoms);
    ESPGroupName[Index].resize(MyAtoms);


    search(is,"ATOM CHARGE        X                 Y                 Z (ANGST)   ESP CHARGE",aline);
    parsewords(aline.c_str(),currentWords);
    if(currentWords[0] == "ATOM" &&
        currentWords[1] == "CHARGE" &&
        currentWords[2] == "X" )
    {
        getwords(currentWords,is);
        while(currentWords[0]!="------------------------------------------------------------")
        {
          natms++;
          double z=atof(currentWords[1].c_str());
          int zint = (int)z;
          atomic_number.push_back(zint);
          q.push_back(z);
          tags.push_back(currentWords[0]);
          pos.push_back(atof(currentWords[2].c_str()));
          pos.push_back(atof(currentWords[3].c_str()));
          pos.push_back(atof(currentWords[4].c_str()));
          ESP[Index].push_back(atof(currentWords[5].c_str()));
          getwords(currentWords,is);
        }
      }
      else
      {
         std::cerr <<"Could not find Atomic coordinates"<< std::endl;
         abort();
      }

    if(natms != MyAtoms)
    {
      std::cerr <<"Could not find atomic coordinates for all atoms. \n";
      abort();
    }

    SpeciesSet& ESPspecies(ESPSystem[Index].getSpeciesSet());
    for(int i=0, ii=0; i<MyAtoms; i++)
     {
       std::string temp;
       std::ostringstream convert1,convert2;
       convert1<<Index;
       convert2<<i;
       ESPSystem[Index].R[i][0]=pos[ii++];
       ESPSystem[Index].R[i][1]=pos[ii++];
       ESPSystem[Index].R[i][2]=pos[ii++];
       temp=IonName[atomic_number[i]]+convert1.str()+convert2.str();
       ESPGroupName[Index].push_back(temp);
       int speciesID = ESPspecies.addSpecies(ESPGroupName[Index].back());
       ESPSystem[Index].GroupID[i]=speciesID;
       ESPspecies(ESPAtomicNumberIndex[Index],speciesID)=1;
       ESPspecies(ESPIonChargeIndex[Index],speciesID)=ESP[Index][i];
      }

    ESPSystem[Index].R *= ang_to_bohr;

  }
  else
  {
     std::cerr <<" Error Finding Monomer coordinates!!! Impossible to create ESP charges!!"<< std::endl;
     abort();

  }
}

void GamesFMOParser::getGeometry(std::istream& is)
{
  //atomic numbers
  const double ang_to_bohr=1.0/0.529177e0;
  std::vector<int> atomic_number,core;
  std::vector<double> q,pos;
  int natms=0;
  tags.clear();
  getwords(currentWords,is);
  if(currentWords[0] == "ATOM" &&
      currentWords[1] == "CHARGE" &&
      currentWords[2] == "X" )
  {
      getwords(currentWords,is);
      while(currentWords[0]!="------------------------------------------------------------")
      {
        natms++;
        double z=atof(currentWords[1].c_str());
        int zint = (int)z;  // is this dangerous???
        atomic_number.push_back(zint);
        q.push_back(z);
        tags.push_back(currentWords[0]);
        pos.push_back(atof(currentWords[2].c_str()));
        pos.push_back(atof(currentWords[3].c_str()));
        pos.push_back(atof(currentWords[4].c_str()));
        //cout<<"Atom number "<<natms<<"  Words:"<<currentWords[0]<<"   "<< currentWords[2]<<"   "<< currentWords[3]<<"   "<< currentWords[4]<<"   "<< std::endl;
        getwords(currentWords,is);

      }
    }
    else
    {
       std::cerr <<"Could not find Atomic coordinates"<< std::endl;
       abort();
    }

    std::cout <<"natms="<<natms<<"    NumberOfAtoms="<<NumberOfAtoms<< std::endl;
  if(natms != NumberOfAtoms)
  {
    std::cerr <<"Could not find atomic coordinates for all atoms. \n";
    abort();
  }
  SpeciesSet& species(IonSystem.getSpeciesSet());
  for(int i=0, ii=0; i<NumberOfAtoms; i++)
  {
    IonSystem.R[i][0]=pos[ii++];
    IonSystem.R[i][1]=pos[ii++];
    IonSystem.R[i][2]=pos[ii++];
    GroupName[i]=IonName[atomic_number[i]];
    int speciesID = species.addSpecies(GroupName[i]);
    IonSystem.GroupID[i]=speciesID;
    species(AtomicNumberIndex,speciesID)=atomic_number[i];
    species(IonChargeIndex,speciesID)=q[i];
  }

    IonSystem.R *= ang_to_bohr;
}




void GamesFMOParser::getGaussianCenters(std::istream& is)
{
  gBound.resize(NumberOfAtoms+1);
  int ng,nx;
  int OldShell=1;
  std::string aline;
  std::map<std::string,int> basisDataMap;
  int nUniqAt=0;
  for(int i=0; i<NumberOfAtoms; i++)
  {
    std::map<std::string,int>::iterator it(basisDataMap.find(tags[i]));
    if(it == basisDataMap.end())
    {
      basisDataMap[tags[i]]=nUniqAt++;
    }
  }

  std::vector<std::vector<double> > expo(nUniqAt),coef(nUniqAt),coef2(nUniqAt);
  std::vector<int> nshll(nUniqAt,0); //use this to
  std::vector<std::vector<int> > ncoeffpershell(nUniqAt);
  std::vector<std::vector<std::string> > shID(nUniqAt);
  std::map<std::string,int> gsMap;
  gsMap[std::string("S")]=1;
  gsMap[std::string("SP")]=2;
  gsMap[std::string("P")]=3;
  gsMap[std::string("D")]=4;
  gsMap[std::string("F")]=5;
  gsMap[std::string("G")]=6;
  gsMap[std::string("H")]=7;
  gsMap[std::string("I")]=8;
  is.seekg(pivot_begin);
  bool found=false;
  while(!found)
  {
    if(is.eof())
    {
      std::cerr <<"Problem with basis set data.\n";
      abort();
    }
    getwords(currentWords,is);
    if(currentWords.size() < 6)
      continue;
    if(currentWords[0] == "SHELL" &&
        currentWords[1] == "TYPE" &&
        currentWords[2] == "PRIMITIVE" &&
        currentWords[3] == "EXPONENT" &&
        currentWords[4] == "CONTRACTION" &&
        currentWords[5] == "COEFFICIENT(S)")
      found=true;
  }


  int currPos=-1;


  while(true)
  {
    getwords(currentWords,is);

    if(currentWords.empty()) continue;

    if(currentWords[0] == "TOTAL" && currentWords[1] == "NUMBER" &&
        currentWords[2] == "OF" && currentWords[3] == "BASIS")
    {
      ng=atoi(currentWords.back().c_str());
      break;
    }
    if(currentWords.size() == 1) //found Species
    {
      std::map<std::string,int>::iterator it(basisDataMap.find(currentWords[0]));
      if(it == basisDataMap.end())
      {
        std::cerr <<"Error in parser.\n";
        abort();
      }
      currPos=it->second;
      bool newgroup=(nshll[currPos]==0);
      if(newgroup)
      {
        ncoeffpershell[currPos].clear();
        ncoeffpershell[currPos].push_back(0);
        shID[currPos].clear();
        shID[currPos].push_back("NONE");
      }

      getwords(currentWords,is); //empty line after species

      while(true)
      {
        std::streampos pivot= is.tellg();
        getwords(currentWords,is);
        if(currentWords.empty()) //empty line after the shell
        {
          if(newgroup)
          {
            nshll[currPos]++;
            ncoeffpershell[currPos].push_back(0);
            shID[currPos].push_back("NONE");
          }
          continue;
        }
        if(currentWords.size()==1 || currentWords[0]=="TOTAL")
        {//use the size and TOTAL to skip to the new group
          is.seekg(pivot);
          break;
        }
        else
        {
          if(newgroup)
          {
            expo[currPos].push_back(atof(currentWords[3].c_str()));
            coef[currPos].push_back(atof(currentWords[4].c_str()));
            ncoeffpershell[currPos][nshll[currPos]]++;
            shID[currPos][nshll[currPos]] = currentWords[1];

            if(gsMap[currentWords[1]] == 2)
            {
              std::cerr <<"Can't handle SP basis states yet. Fix later.\n";
              abort();
            }
            if(gsMap[currentWords[1]] >= 7)
            {
              std::cerr <<"Can't handle H basis states or higher yet. Fix later.\n";
              abort();
            }
            if (debug){
              std::cout << currPos << ":" <<expo[currPos].back() << " " << coef[currPos].back() << " "
                << ncoeffpershell[currPos][nshll[currPos]]
                << " " << shID[currPos][nshll[currPos]] << std::endl;
            }
          }
        }
      }
    }
  }

  gShell.clear();
  gNumber.clear();
  gExp.clear();
  gC0.clear();
  gC1.clear();
  int gtot=0;
  for(int i=0; i<NumberOfAtoms; i++)
  {
    std::map<std::string,int>::iterator it(basisDataMap.find(tags[i]));
    if(it == basisDataMap.end())
    {
      std::cerr <<"Error in parser.\n";
      abort();
    }
    gBound[i] = gtot;
    int indx = it->second;
    gtot+=nshll[indx];
    for(int k=0; k<nshll[indx]; k++)
      gShell.push_back(gsMap[shID[indx][k]]);
    for(int k=0; k<nshll[indx]; k++)
      gNumber.push_back(ncoeffpershell[indx][k]);
    for(int k=0; k<expo[indx].size(); k++)
      gExp.push_back(expo[indx][k]);
    for(int k=0; k<coef[indx].size(); k++)
      gC0.push_back(coef[indx][k]);
  }
  gBound[NumberOfAtoms] = gtot;
//  ng=gShell.size();

}


void GamesFMOParser::getMO(std::istream& is,std::string temp_tag)
{
  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  EigVec.resize(2*SizeOfBasisSet*numMO);
  std::string aline;
  for(int i=0;i<NumberOfAtoms;i++)
     getwords(currentWords,is);
  getwords(currentWords,is);  // ----------------------
  getwords(currentWords,is);  //  TOTAL NUMBER OF MOS IN VARIATION SPACE
  getwords(currentWords,is);  //  Empty Line



  std::vector<double> dummy(50);
  Matrix<double> CartMat(numMO,SizeOfBasisSet);
  std::streampos pivot;
  pivot= is.tellg();
  std::vector<std::string> CartLabels(SizeOfBasisSet);
// this is not the best way, you should use the basis type (e.g. S,P,D,etc) to do this
  getwords(currentWords,is);
  getwords(currentWords,is);
  getwords(currentWords,is);
  if(readtype==2)
    getwords(currentWords,is);
  for(int k=0; k<SizeOfBasisSet; k++)
  {
    getwords(currentWords,is);
    if(currentWords.size() == 8)
    {
      CartLabels[k] = currentWords[2];
      CartLabels[k].erase(0,1); // remove
    }
    else
    {
      CartLabels[k] = currentWords[3];
    }
//cout<<"label: " <<k <<"  " <<CartLabels[k] << std::endl; std::cout.flush();
  }


  is.seekg(pivot);
  getMO_single_set(is, CartMat, EigVal_alpha,temp_tag);
  int cnt=0;
  for(int i=0; i<numMO; i++)
    for(int k=0; k<SizeOfBasisSet; k++)
      EigVec[cnt++] = CartMat[i][k];

// beta states for now
  if(!SpinRestricted)
  {
    search(is,temp_tag);
    for(int i=0;i<NumberOfAtoms;i++)
       getwords(currentWords,is);

    getwords(currentWords,is);  // ----------------------
    getwords(currentWords,is);  // empty line
    getMO_single_set(is, CartMat, EigVal_beta,temp_tag);
  }

  for(int i=0; i<numMO; i++)
    for(int k=0; k<SizeOfBasisSet; k++)
      EigVec[cnt++] = CartMat[i][k];
  std::cout <<"Finished reading MO." << std::endl;
}


void GamesFMOParser::getMO_single_set(std::istream& is, Matrix<double> &CartMat, std::vector<value_type>& EigVal,std::string temp_tag)
{
  int nq = numMO/5;
  int rem = numMO%5;
  int cnt=0;
  for(int i=0; i<nq; i++)
  {
    getwords(currentWords,is);
    if(readtype==2)
      getwords(currentWords,is);
    getwords(currentWords,is);
    EigVal[cnt] = atof(currentWords[0].c_str()) ;
    EigVal[cnt+1] = atof(currentWords[1].c_str()) ;
    EigVal[cnt+2] = atof(currentWords[2].c_str()) ;
    EigVal[cnt+3] = atof(currentWords[3].c_str()) ;
    EigVal[cnt+4] = atof(currentWords[4].c_str()) ;
    getwords(currentWords,is);
    for(int k=0; k<SizeOfBasisSet; k++)
    {
      getwords(currentWords,is);
//cout<<"i,k,size: " <<i <<" " <<k <<" " <<currentWords.size() <<" " <<currentWords[4] << std::endl;
      if(currentWords.size() == 8)
        // G basis or higher TAG gets mixed with atom id
      {
        CartMat[cnt][k] = atof(currentWords[3].c_str()) ;
        CartMat[cnt+1][k] = atof(currentWords[4].c_str()) ;
        CartMat[cnt+2][k] = atof(currentWords[5].c_str()) ;
        CartMat[cnt+3][k] = atof(currentWords[6].c_str()) ;
        CartMat[cnt+4][k] = atof(currentWords[7].c_str()) ;
      }
      else
      {
        CartMat[cnt][k] = atof(currentWords[4].c_str()) ;
        CartMat[cnt+1][k] = atof(currentWords[5].c_str()) ;
        CartMat[cnt+2][k] = atof(currentWords[6].c_str()) ;
        CartMat[cnt+3][k] = atof(currentWords[7].c_str()) ;
        CartMat[cnt+4][k] = atof(currentWords[8].c_str()) ;
      }
    }
    getwords(currentWords,is);
    cnt+=5;
//cout<<"cnt: " <<cnt << std::endl; std::cout.flush();
  }
//cout<<"done with main block, reading rem: " <<rem << std::endl; std::cout.flush();
  if(rem > 0)
  {
    getwords(currentWords,is);
    if(readtype==2)
      getwords(currentWords,is);
    getwords(currentWords,is);
    for(int i=0; i<rem; i++)
    {
      EigVal[cnt+i] = atof(currentWords[i].c_str()) ;
    }
    getwords(currentWords,is);
    for(int k=0; k<SizeOfBasisSet; k++)
    {
      getwords(currentWords,is);
      if(currentWords.size() == 3+rem)
        // G basis or higher TAG gets mixed with atom id
      {
        for(int i=0; i<rem; i++)
        {
          CartMat[cnt+i][k] = atof(currentWords[3+i].c_str()) ;
        }
      }
      else
      {
        for(int i=0; i<rem; i++)
        {
          CartMat[cnt+i][k] = atof(currentWords[4+i].c_str()) ;
        }
      }
    }
    getwords(currentWords,is);
  }
//cout<<"done with rem block, writing eigV: " << std::endl; std::cout.flush();
}

/////////////////////CLASS Atom Definition/////////

Monomer::Monomer():
    X(0),Y(0),Z(0),ESPq(0),q(0),AtomIndex(0)
{}

IDmer::IDmer():
IndexI(0),IndexJ(0),IndexK(0)
{}

void Monomer::print_Geometry(){
   std::cout <<" My Tag="<< tags<<"   X="<<X<<"   Y="<< Y<<" Z="<< Z<<"  Charge="<<q<<"  ESP Charge="<<ESPq<<"  Index="<<AtomIndex<< std::endl;
}

