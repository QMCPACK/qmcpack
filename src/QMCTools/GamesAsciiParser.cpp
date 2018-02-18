//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////





#include "QMCTools/GamesAsciiParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>


void Cartesian2Spherical(int n, double* Cart, double* Sphe);

GamesAsciiParser::GamesAsciiParser()
{
  basisName = "Gaussian";
  Normalized = "no";
  usingECP=false;
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  NFZC=0;
  ECP=false;
}

GamesAsciiParser::GamesAsciiParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "Gaussian";
  Normalized = "no";
  usingECP=false;
  ECP=false;
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  SpinRestricted=true;
  readtype=0;
  NFZC=0;
}

void GamesAsciiParser::parse(const std::string& fname)
{
  std::ifstream fin(fname.c_str());
  if (fin.fail())
  {
    std::cerr << "Error when opening file: " << fname << std::endl;
    abort();
  }
  pivot_begin= fin.tellg();
  std::string aline;
  // if basis functions are removed, this will be modified below
  search(fin,"NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS",aline);
  parsewords(aline.c_str(),currentWords);
  SizeOfBasisSet = atoi(currentWords[6].c_str());
  search(fin,"SPIN MULTIPLICITY",aline);
  parsewords(aline.c_str(),currentWords);
  SpinMultiplicity = atoi(currentWords[2].c_str());
  std::cout <<"SPIN MULTIPLICITY: " <<SpinMultiplicity << std::endl;
  search(fin,"TOTAL NUMBER OF ATOMS",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfAtoms = atoi(currentWords[4].c_str());
  std::cout <<"NUMBER OF ATOMS: " <<NumberOfAtoms << std::endl;
  if(lookFor(fin,"SCFTYP=UHF",aline))
  {
    SpinRestricted=false;
    std::cout <<"Spin Unrestricted MOs" << std::endl;
  }
  else
  {
    std::cout <<"Spin Restricted MOs" << std::endl;
  }
  if(lookFor(fin,"TOTAL NUMBER OF MOS IN VARIATION SPACE=",aline))
  {
    parsewords(aline.c_str(),currentWords);
    numMO = atoi(currentWords[7].c_str());
    std::cout <<"NUMBER OF MOs: " <<numMO << std::endl;
  }
  else
  {
    fin.close();
    fin.open(fname.c_str());
    pivot_begin= fin.tellg();
    if(lookFor(fin,"SET, THE NUMBER OF SPHERICAL HARMONICS KEPT IN THE VARIATION SPACE IS",aline))
    {
      parsewords(aline.c_str(),currentWords);
      numMO = atoi(currentWords[12].c_str());
      std::cout <<"NUMBER OF MOs: " <<numMO << std::endl;
    }
    else
    {
      fin.close();
      fin.open(fname.c_str());
      pivot_begin= fin.tellg();
      std::cout <<"Didn't find reduction of variational space, assuming cartesian number of MO's. \n";
      numMO = SizeOfBasisSet;
      //abort();
    }
  }
  if(numMO2print <= 0) numMO2print = numMO;
  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  getGeometry(fin);
  fin.seekg(pivot_begin);
  if(usingECP)
  {
    std::cout <<"Using ECP." << std::endl;
    ECP=true;
    search(fin,"NUMBER OF ELECTRONS KEPT IN THE CALCULATION IS",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfEls = atoi(currentWords[8].c_str());
    std::cout <<"Number of electrons: " <<NumberOfEls << std::endl;
    std::cout.flush();
    search(fin,"NUMBER OF OCCUPIED ORBITALS (ALPHA) KEPT IS",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfAlpha = atoi(currentWords[7].c_str());
    std::cout <<"Number of alpha electrons: " <<NumberOfAlpha << std::endl;
    search(fin,"NUMBER OF OCCUPIED ORBITALS (BETA ) KEPT IS",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfBeta = atoi(currentWords[8].c_str());
    std::cout <<"Number of beta electrons: " <<NumberOfBeta << std::endl;
  }
  else
  {
    search(fin,"NUMBER OF ELECTRONS ",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfEls = atoi(currentWords[3].c_str());
    std::cout <<"Number of electrons: " <<NumberOfEls << std::endl;
    search(fin,"NUMBER OF OCCUPIED ORBITALS (ALPHA)",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfAlpha = atoi(currentWords[5].c_str());
    std::cout <<"Number of alpha electrons: " <<NumberOfAlpha << std::endl;
    search(fin,"NUMBER OF OCCUPIED ORBITALS (BETA )",aline);
    parsewords(aline.c_str(),currentWords);
    NumberOfBeta = atoi(currentWords[6].c_str());
    std::cout <<"Number of beta electrons: " <<NumberOfBeta << std::endl;
  }
  getGaussianCenters(fin);
  fin.seekg(pivot_begin);
  if(readNO > 0)
    // look for natural orbitals
  {
// output from ALDET and GUGA CI
    std::cout <<"Reading " <<readNO <<" orbitals from file.\n";
    numMO=readNO;
    if(lookFor(fin,"NATURAL ORBITALS IN ATOMIC ORBITAL BASIS"))
    {
      MOtype = "NaturalOrbitals";
      readtype=1;
      std::cout <<"Reading Natural Orbitals from ALDET/GUGA/FSOCI run output. \n";
    }
    else
    {
      fin.close();
      fin.open(fname.c_str());
// output from MCSCF run
      if(lookFor(fin,"MCSCF NATURAL ORBITALS"))
      {
        MOtype = "NaturalOrbitals";
        readtype=2;
        std::cout <<"Reading Natural Orbitals from MCSCF run output. \n";
      }
      else
      {
        std::cerr <<"Could not find Natural Orbitals. \n";
        abort();
      }
    }
  }
  else
    if (readGuess > 0)
    {
      std::cout <<"Reading " <<readGuess <<" orbitals from file.\n";
      numMO=readGuess;
      if(lookFor(fin,"     INITIAL GUESS ORBITALS"))
      {
        MOtype = "InitialGuess";
        readtype=0;
        std::cout <<"Reading INITIAL GUESS ORBITALS output. \n";
      }
      else
      {
        std::cerr <<"Could not find INITIAL GUESS ORBITALS. \n";
        abort();
      }
    }
    else
      // look for eigenvectors
    {
      if(lookFor(fin,"   EIGENVECTORS"))
      {
        MOtype = "Canonical";
        readtype=0;
        numMO=numMO2print;
        std::cout <<"Reading RHF Canonical Orbitals from Gamess output. \n";
      }
      else
      {
        fin.close();
        fin.open(fname.c_str());
// output
        if(lookFor(fin,"MCSCF OPTIMIZED ORBITALS"))
        {
          MOtype = "Canonical";
          readtype=0;
          numMO=numMO2print;
          std::cout <<"Reading Optimized Orbitals from MCSCF run output. \n";
        }
        else
        {
          std::cerr <<"Could not find eigenstates. \n";
          abort();
        }
      }
    }
//  fin.close(); fin.open(fname.c_str());
  getMO(fin);
  fin.close();
// using a possibly different output file for ci coefficients
  if(multideterminant)
  {
    fin.open(outputFile.c_str());
    if (fin.fail())
    {
      std::cerr << "Error when opening file: " << outputFile << std::endl;
      abort();
    }
    pivot_begin= fin.tellg();
//cout<<"looking for dets " << std::endl;
//cout.flush();
    if(lookFor(fin,"GUGA DISTINCT ROW TABLE"))
    {
      std::cout <<"Found GUGA ROW TABLE, reading CSF." << std::endl;
//cout.flush();
      if(!lookFor(fin,"SYMMETRIES FOR THE",aline))
      {
        std::cerr <<"Could not find number of frozen core orbitals in output file.\n";
        abort();
      }
      else
      {
        NFZC = atoi(aline.substr(20,3).c_str());
        NAC = atoi(aline.substr(30,3).c_str());
        NEXT = atoi(aline.substr(42,3).c_str());
        NTOT=NEXT+NAC;
        std::cout <<"# core, #active, #external: " <<NFZC <<" " <<NAC <<" " <<NEXT << std::endl;
      }
//cout.flush();
      fin.seekg(pivot_begin);
      getCSF(fin);
    }
    else
    {
      std::cout <<"Could not find GUGA ROW TABLE, looking for Slater Dets." << std::endl;
//cout.flush();
      fin.close();
      fin.open(outputFile.c_str());
      pivot_begin= fin.tellg();
      if(lookFor(fin,"DIRECT DETERMINANT ORMAS-CI"))
      {
        std::cout <<"Found ORMAS-CI" << std::endl;
//cout.flush();
        fin.close();
        fin.open(outputFile.c_str());
        pivot_begin= fin.tellg();
        getORMAS(fin);
      }
      else
      {
        std::cout <<"Assuming ALDET-CI" << std::endl;
//cout.flush();
        fin.close();
        fin.open(outputFile.c_str());
        pivot_begin= fin.tellg();
        getCI(fin);
      }
    }
    fin.close();
  }
}

void GamesAsciiParser::getGeometry(std::istream& is)
{
  //atomic numbers
  std::vector<int> atomic_number,core;
  std::vector<double> q,pos;
  int natms=0;
  tags.clear();
  is.seekg(pivot_begin);
  //read atomic info
  bool notfound=true;
  do
  {
    if(is.eof())
    {
      std::cerr <<"Could not find atomic coordinates. \n";
      abort();
    }
    getwords(currentWords,is);
    if(currentWords.size() < 4 )
      continue;
    if(currentWords[0] == "ATOM" &&
        currentWords[1] == "ATOMIC" &&
        currentWords[2] == "COORDINATES" &&
        currentWords[3] == "(BOHR)"  )
    {
      getwords(currentWords,is);  // second header line
      notfound=false;
      getwords(currentWords,is);
      while(currentWords.size() != 0)
      {
        if(currentWords[0] == "INTERNUCLEAR")
          break;
        natms++;
        double z=atof(currentWords[1].c_str());
        int zint = (int)z;  // is this dangerous???
        atomic_number.push_back(zint);
        q.push_back(z);  // if using ECPs, change below
        tags.push_back(currentWords[0]);
        pos.push_back(atof(currentWords[2].c_str()));
        pos.push_back(atof(currentWords[3].c_str()));
        pos.push_back(atof(currentWords[4].c_str()));
        getwords(currentWords,is);
      }
    }
  }
  while(notfound);
  // effective charges are read from ECP section
  if(natms != NumberOfAtoms)
  {
    std::cerr <<"Could not find atomic coordinates for all atoms. \n";
    abort();
  }
// this is risky but works for now
  is.seekg(pivot_begin);
  notfound=true;
  while(notfound)
  {
    if(is.eof())
    {
      std::cerr <<"Problem looking for ECPs, this should not happen. Contact developers for help. \n";
      abort();
    }
    getwords(currentWords,is);
// this should appear below the ECP section in the output file
// so use this to avoid going all the way to the bottom
    if(currentWords.size() < 2 )
      continue;
    if( currentWords[0] == "ECP" &&
        currentWords[1] == "POTENTIALS")
      // eureka!!!
    {
      usingECP=true;
      ECP=true;
      core.resize(NumberOfAtoms);
      getwords(currentWords,is); // -------------
// this only works if all atoms have an ECP, fix later
//      for(int i=0; i<NumberOfAtoms; i++) {
// fixing this problem
      bool done=false;
      while(!done)
      {
        if(is.eof())
        {
          std::cerr <<"Found ECPs, but problem looking ZCORE data.\n";
          abort();
        }
        getwords(currentWords,is);
        if(currentWords.size() == 0)
          continue;
        if(currentWords.size() >= 4)
        {
          if(currentWords[0] == "THE" &&
              currentWords[1] == "ECP"  &&
              currentWords[2] == "RUN"  &&
              currentWords[3] == "REMOVES" )
          {
            done=true;
          }
        }
        if(currentWords[0] == "PARAMETERS" &&
            currentWords[1] == "FOR" )
        {
          //done=true;
          std::vector<std::string>::iterator it,it0;
          it = find(currentWords.begin(),currentWords.end(),"ZCORE");
          it0 = find(currentWords.begin(),currentWords.end(),"ATOM");
          if(it0 == currentWords.end())
          {
            std::cerr <<"Problem with ECP data. Didn't found ATOM tag\n";
            std::cerr << is.rdbuf() << std::endl;
            abort();
          }
          it0++;
          int nq0 = atoi(it0->c_str())-1;
          if(it != currentWords.end())
          {
            it++;
            core[nq0] = atoi(it->c_str());
            q[nq0] -= core[nq0];
            std::cout <<"Found ECP for atom " <<nq0 <<" with zcore " <<core[nq0] << std::endl;
          }
          else
          {
            it = find(currentWords.begin(),currentWords.end(),"ATOM");
            if(it == currentWords.end())
            {
              std::cerr <<"Problem with ECP data. Didn't found ATOM tag\n";
              std::cerr <<"Atom: " <<nq0 << std::endl;
              abort();
            }
            std::vector<std::string>::iterator it2=it;
            it2++;
            int nq = atoi(it2->c_str());
            if(nq != nq0+1)
            {
              std::cerr <<"Problem with ECP data. ID's don't agree\n";
              std::cerr <<"Atom: " <<nq0 << std::endl;
              abort();
            }
            it = find(it2,currentWords.end(),"ATOM");
            if(it == currentWords.end())
            {
              std::cerr <<"Problem with ECP data (2).\n";
              std::cerr <<"Atom: " <<nq0 << std::endl;
              abort();
            }
            nq = atoi((it+1)->c_str());
            core[nq0] = core[nq-1];
            q[nq0] -= core[nq0];
            std::cout <<"Found ECP for atom " <<nq0 <<" with zcore " <<core[nq0] << std::endl;
          }
        }
      }
      notfound=false;
    }
    else
    {
      if(currentWords.size() < 3 )
        continue;
      if( currentWords[0] == "1" &&
          currentWords[1] == "ELECTRON" &&
          currentWords[2] == "INTEGRALS" )
        break;
    }
  }
  std::cout <<"usingECP: " <<(usingECP?("yes"):("no")) << std::endl;
  std::cout.flush();
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
}

void GamesAsciiParser::getGaussianCenters(std::istream& is)
{
  gBound.resize(NumberOfAtoms+1);
  int ng,nx;
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
  gsMap[std::string("L")]=2;
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

  getwords(currentWords,is);  // empty line

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

            if(gsMap.find(currentWords[1]) == gsMap.end())
            {
              std::cerr <<"Unhandled primitive type " << currentWords[1] << std::endl;
              abort();
            }
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
            if(debug){
               std::cout << currPos << ":" <<expo[currPos].back() << " " << coef[currPos].back() << " "
                 << ncoeffpershell[currPos][nshll[currPos]]
                 << " " << shID[currPos][nshll[currPos]] << std::endl;
            }
          }
        }
      }
    }
  }


  /*
  getwords(currentWords,is);  // tag of first atom
  for(int i=0; i<nUniqAt-1; i++)
  {
    int currPos;
    if(currentWords.size() == 0)
    {
      std::cerr <<"Error in parser.\n";
      abort();
    }
    std::map<std::string,int>::iterator it(basisDataMap.find(currentWords[0]));
    if(it == basisDataMap.end())
    {
      std::cerr <<"Error in parser.\n";
      abort();
    }
    currPos=it->second;
    getwords(currentWords,is); // empty line
    nshll[currPos]=0;
    ncoeffpershell[currPos].clear();
    ncoeffpershell[currPos].push_back(0);
    shID[currPos].clear();
    shID[currPos].push_back("NONE");
    while(true)
    {
      getwords(currentWords,is);
      if(currentWords.size() == 0)
      {
        nshll[currPos]++;
        ncoeffpershell[currPos].push_back(0);
        shID[currPos].push_back("NONE");
        continue;
      }
      if(basisDataMap.find(currentWords[0]) != basisDataMap.end())
        break;
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
    }
  }
  {
    // one last time
    int currPos;
    if(currentWords.size() == 0)
    {
      std::cerr <<"Error in parser.\n";
      abort();
    }
    std::map<std::string,int>::iterator it(basisDataMap.find(currentWords[0]));
    if(it == basisDataMap.end())
    {
      std::cerr <<"Error in parser.\n";
      abort();
    }
    currPos=it->second;
    getwords(currentWords,is); // empty line
    nshll[currPos]=0;
    ncoeffpershell[currPos].clear();
    ncoeffpershell[currPos].push_back(0);
    shID[currPos].clear();
    shID[currPos].push_back("NONE");
    while(true)
    {
      getwords(currentWords,is);
      if(currentWords.size() == 0)
      {
        nshll[currPos]++;
        ncoeffpershell[currPos].push_back(0.0);
        shID[currPos].push_back("NONE");
        continue;
      }
      if(currentWords[0] == "TOTAL" && currentWords[1] == "NUMBER" &&
          currentWords[2] == "OF" && currentWords[3] == "BASIS")
      {
        ng=atoi(currentWords[7].c_str());
        break;
      }
      expo[currPos].push_back(atof(currentWords[3].c_str()));
      coef[currPos].push_back(atof(currentWords[4].c_str()));
      ncoeffpershell[currPos][nshll[currPos]]++;
      shID[currPos][nshll[currPos]] = currentWords[1];
      if(gsMap[currentWords[1]] == 2)
      {
        std::cerr <<"Can't handle SP basis states yet. Fix later.\n";
        abort();
      }
    }
  }
*/
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
}

void GamesAsciiParser::getMO(std::istream& is)
{
  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  EigVec.resize(2*SizeOfBasisSet*numMO);
  std::string aline;
  //if(MOtype == "Canonical")
  //  search(is,"   EIGENVECTORS");
  //else if(MOtype == "NaturalOrbitals") {
  //  if(readtype==1)
  //    search(is,"NATURAL ORBITALS IN ATOMIC ORBITAL BASIS"); // ci
  //  else if(readtype == 2)
  //    search(is,"MCSCF NATURAL ORBITALS");  // mcscf
  //}
  getwords(currentWords,is);  // ----------------------
  getwords(currentWords,is);  // empty line
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
  getMO_single_set(is, CartMat, EigVal_alpha);
  int cnt=0;
  for(int i=0; i<numMO; i++)
    for(int k=0; k<SizeOfBasisSet; k++)
      EigVec[cnt++] = CartMat[i][k];

// beta states for now
  if(!SpinRestricted)
  {
    search(is,"   EIGENVECTORS");
    getwords(currentWords,is);  // ----------------------
    getwords(currentWords,is);  // empty line
    getMO_single_set(is, CartMat, EigVal_beta);
  }

  for(int i=0; i<numMO; i++)
    for(int k=0; k<SizeOfBasisSet; k++)
      EigVec[cnt++] = CartMat[i][k];
  std::cout <<"Finished reading MO." << std::endl;
}

void GamesAsciiParser::getMO_single_set(std::istream& is, Matrix<double> &CartMat, std::vector<value_type>& EigVal)
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

void Cartesian2Spherical(int n, double* Cart, double* Sphe)
{
  switch(n)
  {
  case 1:
  {
    *Sphe = *Cart;
    break;
  }
  case 3:
  {
    // m = -1
    *(Sphe) = *(Cart+1);
    // m = 0
    *(Sphe+1) = *(Cart+2);
    // m = 1
    *(Sphe+2) = *(Cart);
    break;
  }
  case 5:
  {
    // m = -2
    *(Sphe)   = *(Cart+3);
    // m = -1
    *(Sphe+1) = *(Cart+5);
    // m = 0
    *(Sphe+2) = *(Cart+2) - 0.5*(*(Cart)+*(Cart+1));
    // m = 1
    *(Sphe+3) = *(Cart+4);
    // m = 2
    *(Sphe+4) = std::sqrt(0.75)*(*(Cart)-*(Cart+1));
    break;
  }
  case 7:
  {
    // m = -3
    *(Sphe)   = -1.0*std::sqrt(5.0/8.0)*(*(Cart+1)) + std::sqrt(9.0/8.0)*(*(Cart+3)) ;
    // m = -2
    *(Sphe+1) = *(Cart+9);
    // m = -1
    *(Sphe+2) = std::sqrt(6.0/5.0)*(*(Cart+8)) - std::sqrt(3.0/8.0)*(*(Cart+1)) - std::sqrt(6.0/5.0)*(*(Cart+3))/4.0;
    // m = 0
    *(Sphe+3) = *(Cart+2) - 3.0/std::sqrt(10.0)*(*(Cart+4)+*(Cart+6));
    // m = 1
    *(Sphe+4) = std::sqrt(6.0/5.0)*(*(Cart+7)) - std::sqrt(3.0/8.0)*(*(Cart)) - std::sqrt(6.0/5.0)*(*(Cart+5))/4.0;
    // m = 2
    *(Sphe+5) = std::sqrt(3.0/4.0)*(*(Cart+4)-*(Cart+6));
    // m = 3
    *(Sphe+6) = -1.0*std::sqrt(5.0/8.0)*(*(Cart)) + std::sqrt(9.0/8.0)*(*(Cart+5)) ;
    break;
  }
  case 9:
  {
    // m = -4
    *(Sphe)   = std::sqrt(5.0/4.0)*(*(Cart+3)-*(Cart+5));
    // m = -3
    *(Sphe+1) = -1.0*std::sqrt(5.0/8.0)*(*(Cart+6))+std::sqrt(9.0/8.0)*(*(Cart+12));
    // m = -2
    *(Sphe+2) = std::sqrt(9.0/7.0)*(*(Cart+14))-std::sqrt(5.0/28.0)*(*(Cart+3)+*(Cart+5));
    // m = -1
    *(Sphe+3) = std::sqrt(10.0/7.0)*(*(Cart+8))-0.75*std::sqrt(10.0/7.0)*(*(Cart+6))-0.75*std::sqrt(2.0/7.0)*(*(Cart+12));
    // m = 0
    *(Sphe+4) = *(Cart+2)+std::sqrt(9.0/32.0)*(*(Cart)+*(Cart+1))-3.0*std::sqrt(6.0/35.0)*(*(Cart+10)+*(Cart+11)-0.25*(*(Cart+9)));
    // m = 1
    *(Sphe+5) = std::sqrt(10.0/7.0)*(*(Cart+7))-0.75*std::sqrt(10.0/7.0)*(*(Cart+4))-0.75*std::sqrt(2.0/7.0)*(*(Cart+13));
    // m = 2
    *(Sphe+6) = 1.5*std::sqrt(3.0/7.0)*(*(Cart+10)-*(Cart+11))-std::sqrt(5.0/16.0)*(*(Cart)-*(Cart+1));
    // m = 3
    *(Sphe+7) = std::sqrt(5.0/8.0)*(*(Cart+4))-std::sqrt(9.0/8.0)*(*(Cart+13));
    // m = 4
    *(Sphe+8) = std::sqrt(35.0)/8.0*(*(Cart)+*(Cart+1)) - std::sqrt(3.0)*0.75*(*(Cart+9));
    break;
  }
  /* GAMESS doesn't allow H or higher
      case 11:
      {
        // m = -5
        *(Sphe)   = 0.75*std::sqrt(7.0/8.0)*(*(Cart)+*(Cart))+std::sqrt()*(*(Cart)+*(Cart))-std::sqrt()*(*(Cart)+*(Cart));
        // m = -4
        *(Sphe+1)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = -3
        *(Sphe+2)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = -2
        *(Sphe+3)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = -1
        *(Sphe+4)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = 0
        *(Sphe+5)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = 1
        *(Sphe+6)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = 2
        *(Sphe+7)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = 3
        *(Sphe+8)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = 4
        *(Sphe+9)   = std::sqrt()*(*(Cart)+*(Cart));
        // m = 5
        *(Sphe+10)   = std::sqrt()*(*(Cart)+*(Cart));
      }
  */
  default:
  {
    std::cerr <<"Error in Cartesian2Spherical. Invalid n: " <<n << std::endl;
    abort();
  }
  }
}

// read CSF from output file, NPRT=2 is required in the $CIDRT/DRT
// section
void GamesAsciiParser::getCSF(std::istream& is)
{
  //look for CI coefficients, only working for slater dets right now
  bool notfound=true;
  ci_size=0;
  CSFocc.clear();
  CSFalpha.clear();
  CSFbeta.clear();
  coeff2csf.clear();
  usingCSF=true;

  // set a count to check if we arrive our target state or not
  int state_num = -1;

  std::cout << "Target State Number is " << target_state << std::endl;

  do
  {
    if(is.eof())
    {
      std::cerr <<"Could not find CSF expansion. \n";
      abort();
    }
    getwords(currentWords,is);
    if(currentWords.size() < 3 )
      continue;
    if(currentWords[0] == "CSF" &&
        currentWords[1] == "COEF" &&
        currentWords[2] == "OCCUPANCY" )
    {

      // add the state number by one
      state_num++;

      // if we have not reached target state, continue
      if (state_num != target_state)
        continue;

      getwords(currentWords,is);  // --------
      notfound=false;
      getwords(currentWords,is);
      while(currentWords.size() != 0)
      {
        if(currentWords[0] == "......" || currentWords[1] == "END")
          break;
        double cof = atof(currentWords[1].c_str());
        if(std::abs(cof) > ci_threshold)
        {
          ci_size++;
          int nq = atoi(currentWords[0].c_str());
          std::pair<int,double> cic(nq,cof);
          coeff2csf.push_back(cic);
          if(NTOT < 50)
          {
            CSFocc.push_back(currentWords[2]);
          }
          else
            if(NTOT == 50)
            {
              CSFocc.push_back(currentWords[2]);
              getwords(currentWords,is);
            }
            else
            {
              std::string tmp=currentWords[2];
              getwords(currentWords,is);
              tmp+=currentWords[0];
              CSFocc.push_back(tmp);
            }
          getwords(currentWords,is);
        }
        else
        {
          if(NTOT < 50)
            getwords(currentWords,is);
          else
          {
            getwords(currentWords,is);
            getwords(currentWords,is);
          }
        }
      }
    }
  }
  while(notfound);
  std::cout <<"Done reading csf coefficients." << std::endl;
  std::cout <<"Found: " <<coeff2csf.size() <<" CSFs.\n";
  std::cout.flush();
// look for highest occupied MO to avoid using unnecesary ones
  ci_nstates = 0;
  for(int i=0; i<CSFocc.size(); i++)
  {
    int max=CSFocc[i].size();
    for(int k=CSFocc[i].size()-1; k>=0; k--)
    {
      if(CSFocc[i][k] == '1' || CSFocc[i][k] == '2')
      {
        max = k+1;
        break;
      }
    }
    if(ci_nstates < max)
      ci_nstates=max;
  }
  CSFalpha.resize(ci_size);
  CSFbeta.resize(ci_size);
  CSFexpansion.resize(ci_size);
// now rewind and look for CSF definitions
  is.seekg(pivot_begin);
  if(!lookFor(is,"DETERMINANT CONTRIBUTION TO CSF'S"))
  {
    std::cerr <<"Could not find CSF determinant contributions. Please use NPRT=2 in $CIDRT/DRT input section of gamess. \n";
    abort();
  }
  getwords(currentWords,is); // ----------------
  getwords(currentWords,is); // ----------------
  if(currentWords[0] != "CASE" || currentWords[1] != "VECTOR")
  {
    std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (1). \n";
    abort();
  }
  int ds=SpinMultiplicity-1;
  int neb= (NumberOfEls-ds)/2;
  int nea= NumberOfEls-NumberOfBeta;
  ci_nca = ci_ncb = NFZC;
  std::vector<int> csfOccup;
  bool done=false,first=true;
  int cnt=1,current=0;
  int naea=0,naeb=0;
  std::string aline;
  while(current < ci_size)
  {
    if(is.eof())
    {
      std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (2). \n";
      abort();
    }
    getwords(currentWords,is);
    getwords(currentWords,is);
    getwords(currentWords,is);
// checking
    //if(currentWords[0] != "FOR" || currentWords[1] != "MS") {
    //  std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (3). \n";
    //  abort();
    //}
    getwords(currentWords,is,aline);
    if(aline.substr(1,3) != "CSF")
    {
      std::cerr <<"aline:" <<aline << std::endl;
      std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (4). \n";
      abort();
    }
    if(coeff2csf[current].first == cnt)
      // read dets
    {
// first time the std::string is longer
      csfOccup.clear();
      {
        std::string alp(ci_nstates,'0'),beta(ci_nstates,'0');
        if(first)
        {
          first=false;
          int num=(aline.size()-33)/3;
          for(int i=0; i<num; i++)
          {
            //int nq = atoi(currentWords[i].c_str());
            int nq = atoi(aline.substr(33+i*3,3).c_str());
            csfOccup.push_back(nq);
            if(nq > 0)
            {
              if(nq-1 >= ci_nstates+ci_nca)
              {
                std::cerr <<"Problems with det std::string #,nq,i: " <<cnt <<"  " <<nq <<"  " <<i << std::endl;
                std::cout <<"line: " <<aline << std::endl;
                std::cout <<"alpha: " <<alp << std::endl;
                std::cout <<"beta: " <<beta << std::endl;
                for(int i=6; i<currentWords.size(); i++)
                  std::cerr <<currentWords[i] <<" ";
                std::cerr << std::endl;
                abort();
              }
              alp.at(nq-1-ci_nca) = '1';
              naea++;
            }
            else
            {
              if(-nq-1 >= ci_nstates+ci_ncb)
              {
                std::cerr <<"Problems with det std::string #,nq,i: " <<cnt <<"  " <<nq <<"  " <<i << std::endl;
                std::cout <<"line: " <<aline << std::endl;
                std::cout <<"alpha: " <<alp << std::endl;
                std::cout <<"beta: " <<beta << std::endl;
                for(int i=6; i<currentWords.size(); i++)
                  std::cerr <<currentWords[i] <<" ";
                std::cerr << std::endl;
                abort();
              }
              beta.at(-nq-1-ci_ncb) = '1';
              naeb++;
            }
          }
        }
        else
        {
          int na=0,nb=0;
          int num=(aline.size()-33)/3;
          for(int i=0; i<num; i++)
          {
            //int nq = atoi(currentWords[i].c_str());
            int nq = atoi(aline.substr(33+i*3,3).c_str());
            csfOccup.push_back(nq);
            if(nq > 0)
            {
              if(nq-1 >= ci_nstates+ci_nca)
              {
                std::cerr <<"Problems with det std::string #,nq,i: " <<cnt <<"  " <<nq <<"  " <<i << std::endl;
                std::cout <<"line: " <<aline << std::endl;
                std::cout <<"alpha: " <<alp << std::endl;
                std::cout <<"beta: " <<beta << std::endl;
                for(int i=6; i<currentWords.size(); i++)
                  std::cerr <<currentWords[i] <<" ";
                std::cerr << std::endl;
                abort();
              }
              alp.at(nq-1-ci_nca) = '1';
              na++;
            }
            else
            {
              if(-nq-1 >= ci_nstates+ci_ncb)
              {
                std::cerr <<"Problems with det std::string #,nq,i: " <<cnt <<"  " <<nq <<"  " <<i << std::endl;
                std::cout <<"line: " <<aline << std::endl;
                std::cout <<"alpha: " <<alp << std::endl;
                std::cout <<"beta: " <<beta << std::endl;
                for(int i=6; i<currentWords.size(); i++)
                  std::cerr <<currentWords[i] <<" ";
                std::cerr << std::endl;
                abort();
              }
              beta.at(-nq-1-ci_ncb) = '1';
              nb++;
            }
          }
          if(na != naea || nb != naeb)
          {
            std::cerr <<"Problems with det std::string #: " <<cnt << std::endl;
            std::cout <<"line: " <<aline << std::endl;
            std::cout <<"alpha: " <<alp << std::endl;
            std::cout <<"beta: " <<beta << std::endl;
            for(int i=6; i<currentWords.size(); i++)
              std::cerr <<currentWords[i] <<" ";
            std::cerr << std::endl;
            abort();
          }
        }
        double sg = getCSFSign(csfOccup);
        CSFalpha[current].push_back(alp);
        CSFbeta[current].push_back(beta);
        //CSFexpansion[current].push_back(atof(currentWords[4].c_str())*sg);
        CSFexpansion[current].push_back(atof(aline.substr(20,9).c_str())*sg);
        getwords(currentWords,is,aline);
      }
      while(currentWords.size() != 0)
      {
        if(is.eof())
        {
          std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (5). \n";
          abort();
        }
        if(currentWords[0] == "CASE" && currentWords[1] == "VECTOR")
        {
          cnt++;
          if(cnt < 10000000 && atoi(currentWords[2].c_str()) != cnt)
          {
            std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (6). \n";
            abort();
          }
          break;
        }
        csfOccup.clear();
        std::string alp(ci_nstates,'0'),beta(ci_nstates,'0');
        int num=(aline.size()-33)/3;
        for(int i=0; i<num; i++)
        {
          //int nq = atoi(currentWords[i].c_str());
          int nq = atoi(aline.substr(33+i*3,3).c_str());
          csfOccup.push_back(nq);
          if(nq > 0)
          {
            if(nq-1 >= ci_nstates+ci_nca)
            {
              std::cerr <<"Problems with det std::string #,nq,i: " <<cnt <<"  " <<nq <<"  " <<i << std::endl;
              std::cout <<"line: " <<aline << std::endl;
              std::cout <<"alpha: " <<alp << std::endl;
              std::cout <<"beta: " <<beta << std::endl;
              for(int i=6; i<currentWords.size(); i++)
                std::cerr <<currentWords[i] <<" ";
              std::cerr << std::endl;
              abort();
            }
            alp.at(nq-1-ci_nca) = '1';
          }
          else
          {
            if(-nq-1 >= ci_nstates+ci_ncb)
            {
              std::cerr <<"Problems with det std::string #,nq,i: " <<cnt <<"  " <<nq <<"  " <<i << std::endl;
              std::cout <<"line: " <<aline << std::endl;
              std::cout <<"alpha: " <<alp << std::endl;
              std::cout <<"beta: " <<beta << std::endl;
              for(int i=6; i<currentWords.size(); i++)
                std::cerr <<currentWords[i] <<" ";
              std::cerr << std::endl;
              abort();
            }
            beta.at(-nq-1-ci_ncb) = '1';
          }
        }
        double sg = getCSFSign(csfOccup);
        CSFalpha[current].push_back(alp);
        CSFbeta[current].push_back(beta);
        //CSFexpansion[current].push_back(atof(currentWords[2].c_str())*sg);
        CSFexpansion[current].push_back(atof(aline.substr(20,9).c_str())*sg);
        getwords(currentWords,is,aline);
      }
      current++;
    }
    else
      // not interested in this CSF, so read until next
    {
      while(currentWords.size() != 0)
      {
        if(is.eof())
        {
          std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (5). \n";
          abort();
        }
        if(currentWords[0] == "CASE" && currentWords[1] == "VECTOR")
        {
          cnt++;
          if(cnt < 10000000 && atoi(currentWords[2].c_str()) != cnt)
          {
            std::cerr <<"Problems reading DETERMINANT CONTRIBUTION TO CSF'S (6). \n";
            abort();
          }
          break;
        }
        getwords(currentWords,is,aline);
      }
    }
  }
  std::cout <<"Done reading csf expansion." << std::endl;
  std::cout.flush();
  ci_nea=0;
  for(int i=0; i<CSFalpha[0][0].size(); i++)
    if(CSFalpha[0][0].at(i)=='1')
      ci_nea++;
  ci_neb=0;
  for(int i=0; i<CSFbeta[0][0].size(); i++)
    if(CSFbeta[0][0].at(i)=='1')
      ci_neb++;
//  int ds=SpinMultiplicity-1;
//  int neb= (NumberOfEls-ds)/2;
//  int nea= NumberOfEls-NumberOfBeta;
//  ci_nca = nea-ci_nea;
//  ci_ncb = neb-ci_neb;
  /*
    std::cout <<"Summary. #ci: " <<ci_size << std::endl;
    for(int i=0; i<ci_size; i++) {
      std::cout <<"c: " <<coeff2csf[i].second << std::endl;
      for(int k=0; k<CSFexpansion[i].size(); k++)
        std::cout <<"    " <<k <<"  " <<CSFexpansion[i][k]
            <<"  " <<CSFalpha[i][k]
            <<"  " <<CSFbeta[i][k] << std::endl;
    }
  */
}

void GamesAsciiParser::getCI(std::istream& is)
{
  is.seekg(pivot_begin);
  //look for CI coefficients
  bool notfound=true;
  ci_size=0;
  CIcoeff.clear();
  CIalpha.clear();
  CIbeta.clear();
  do
  {
    if(is.eof())
    {
      std::cerr <<"Could not find CI expansion. \n";
      abort();
    }
    getwords(currentWords,is,0,std::string("|"));
    if(currentWords.size() < 3 )
      continue;
    if(currentWords[0].find("ALP") == 0 &&
        currentWords[1].find("BET") == 0 &&
        currentWords[2] == "COEFFICIENT" )
    {
      getwords(currentWords,is);  // --------
      notfound=false;
      getwords(currentWords,is);
      while(currentWords.size() != 0)
      {
        if(currentWords[0] == "....." || currentWords[1] == "DONE")
          break;
        ci_size++;
        CIcoeff.push_back(atof(currentWords[4].c_str()));
        CIalpha.push_back(currentWords[0]);
        CIbeta.push_back(currentWords[2]);
        getwords(currentWords,is);
      }
    }
  }
  while(notfound);
  ci_nea=ci_neb=0;
  for(int i=0; i<CIalpha[0].size(); i++)
    if(CIalpha[0].at(i) == '1')
      ci_nea++;
  for(int i=0; i<CIbeta[0].size(); i++)
    if(CIbeta[0].at(i) == '1')
      ci_neb++;
  if(CIalpha[0].size() != CIbeta[0].size())
  {
    std::cerr <<"QMCPack can't handle different number of active orbitals in alpha and beta channels right now. Contact developers for help (Miguel).\n";
    abort();
  }
  int ds=SpinMultiplicity-1;
  int neb= (NumberOfEls-ds)/2;
  int nea= NumberOfEls-NumberOfBeta;
  ci_nca = nea-ci_nea;
  ci_ncb = neb-ci_neb;
  ci_nstates = CIalpha[0].size();
}

void GamesAsciiParser::getORMAS(std::istream& is)
{
  is.seekg(pivot_begin);
  //look for CI coefficients
  bool notfound=true;
  ci_size=0;
  CIcoeff.clear();
  CIalpha.clear();
  CIbeta.clear();
  std::string aline;
  if(!lookFor(is,"NUMBER OF CORE ORBITALS",aline))
  {
    std::cerr <<"Couldn't find # of CORE ORBITALS in ORMAS.\n";
    abort();
  }
  parsewords(aline.c_str(),currentWords);
  ci_nca = ci_ncb = atoi(currentWords[4].c_str());
  if(!lookFor(is,"NUMBER OF ACTIVE ORBITALS",aline))
  {
    std::cerr <<"Couldn't find # of ACTIVE ORBITALS in ORMAS.\n";
    abort();
  }
  parsewords(aline.c_str(),currentWords);
  int nactive(atoi(currentWords[4].c_str()));
  if(!lookFor(is,"NUMBER OF ALPHA ELECTRONS",aline))
  {
    std::cerr <<"Couldn't find # of ALPHA ELECTRONS in ORMAS.\n";
    abort();
  }
  parsewords(aline.c_str(),currentWords);
  //ci_nea = atoi(currentWords[4].c_str());
  ci_nea = atoi(currentWords[6].c_str());
  if(!lookFor(is,"NUMBER OF BETA ELECTRONS",aline))
  {
    std::cerr <<"Couldn't find # of BETA ELECTRONS in ORMAS.\n";
    abort();
  }
  parsewords(aline.c_str(),currentWords);
  //ci_neb = atoi(currentWords[4].c_str());
  ci_neb = atoi(currentWords[6].c_str());
  std::cout <<"ORMAS: nea,neb,ncore,nact: "
       <<ci_nea <<" "
       <<ci_neb <<" "
       <<ci_nca <<" "
       <<nactive <<"\n";
  int ds=SpinMultiplicity-1;
  int neb= (NumberOfEls-ds)/2;
  int nea= NumberOfEls-NumberOfBeta;
  if( ci_nca != nea-ci_nea)
  {
    std::cerr <<"Inconsistent number of core electrons: " <<ci_nca <<" " <<nea-ci_nea << std::endl;
    abort();
  }
  if( ci_ncb != neb-ci_neb)
  {
    std::cerr <<"Inconsistent number of core electrons: " <<ci_ncb <<" " <<neb-ci_neb << std::endl;
    abort();
  }
  std::string dummy_alpha(nactive,'0');
  std::string dummy_beta(nactive,'0');
  int nskip = ci_nea+ci_neb+2;
  do
  {
    if(is.eof())
    {
      std::cerr <<"Could not find ORMAS CI expansion. \n";
      abort();
    }
    getwords(currentWords,is);
    if(currentWords.size() < 5 )
      continue;
    if(currentWords[0] == "ALPHA" &&
        currentWords[2] == "BETA" &&
        currentWords[4] == "COEFFICIENT" )
    {
      getwords(currentWords,is);  // 1      2
      getwords(currentWords,is);  // --------
      notfound=false;
      getwords(currentWords,is);
      while(currentWords.size() != 0)
      {
        if(currentWords[0] == "....." || currentWords[1] == "DONE")
          break;
        double cof = atof(currentWords[nskip].c_str());
        if(std::abs(cof) > ci_threshold)
        {
          ci_size++;
          CIcoeff.push_back(cof);
          CIalpha.push_back(dummy_alpha);
          CIbeta.push_back(dummy_beta);
          for(int i=0; i<ci_nea; i++)
            (CIalpha.back())[atoi(currentWords[i].c_str())-1]='1';
          for(int i=0; i<ci_neb; i++)
            (CIbeta.back())[atoi(currentWords[ci_nea+1+i].c_str())-1]='1';
        }
        getwords(currentWords,is);
      }
    }
  }
  while(notfound);
  ci_nstates = 0;
  for(int i=0; i<ci_size; i++)
  {
    int max=0; //=nactive;
    for(int k=nactive-1; k>=0; k--)
    {
      if(CIalpha[i][k] == '1' || CIbeta[i][k] == '1')
      {
        max = k+1;
        break;
      }
    }
    //cout<<i <<" " <<max << std::endl;
    if(ci_nstates < max)
      ci_nstates=max;
  }
}

double GamesAsciiParser::getCSFSign(std::vector<int> & occ)
{
// reference ordering is irrelevant as long as it is consistent
// within all determinants.
  double res=1.0;
  int n=occ.size();
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
    {
      if(occ[j] > occ[i])
      {
        res*=-1.0;
        double tmp=occ[i];
        occ[i]=occ[j];
        occ[j]=tmp;
      }
    }
  return res;
}


