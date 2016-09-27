//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "QMCTools/GaussianFCHKParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>


GaussianFCHKParser::GaussianFCHKParser()
{
  basisName = "Gaussian-G2";
  Normalized = "no";
}

GaussianFCHKParser::GaussianFCHKParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "Gaussian-G2";
  Normalized = "no";
}

void GaussianFCHKParser::parse(const std::string& fname)
{
  std::ifstream fin(fname.c_str());
  getwords(currentWords,fin); //1
  Title = currentWords[0];
  getwords(currentWords,fin);//2  SP RHF Gen
  //if(currentWords[1]=="ROHF" || currentWords[1]=="UHF") {
  if(currentWords[1]=="UHF")
  {
// mmorales: this should be determined by the existence of "Beta MO", since
// there are many ways to get unrestricted runs without UHF (e.g. UMP2,UCCD,etc)
    SpinRestricted=false;
    //std::cout << " Spin Unrestricted Calculation (UHF). " << std::endl;
  }
  else
  {
    SpinRestricted=true;
    //std::cout << " Spin Restricted Calculation (RHF). " << std::endl;
  }
  getwords(currentWords,fin);//3  Number of atoms
  NumberOfAtoms = atoi(currentWords.back().c_str());
  // TDB: THIS FIX SHOULD BE COMPATIBLE WITH MY OLD FCHK FILES
  //getwords(currentWords,fin); //4 Charge
  bool notfound = true;
  while(notfound)
  {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++)
    {
      if("Charge" == currentWords[i])
      {
        notfound = false;
      }
    }
  }
  getwords(currentWords,fin); //5 Multiplicity
  SpinMultiplicity=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //6 Number of electrons
  NumberOfEls=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //7 Number of alpha electrons
  int nup=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //8 Number of beta electrons
  int ndown=atoi(currentWords.back().c_str());
  //NumberOfAlpha=nup-ndown;
  NumberOfAlpha=nup;
  NumberOfBeta=ndown;
  getwords(currentWords,fin); //9 Number of basis functions
  SizeOfBasisSet=atoi(currentWords.back().c_str());
  getwords(currentWords,fin); //10 Number of independant functions
  int NumOfIndOrb=atoi(currentWords.back().c_str());
  std::cout <<"Number of independent orbitals: " <<NumOfIndOrb << std::endl;
  numMO = NumOfIndOrb;
  // TDB: THIS ADDITION SHOULD BE COMPATIBLE WITH MY OLD FCHK FILES
  std::streampos pivottdb = fin.tellg();
  int ng;
  notfound = true;
  while(notfound)
  {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++)
    {
      if("contracted" == currentWords[i])
      {
        ng=atoi(currentWords.back().c_str());
        notfound = false;
      }
    }
  }
  // TDB: THIS FIX SHOULD BE COMPATIBLE WITH MY OLD FCHK FILES
  // getwords(currentWords,fin); //Number of contracted shells
  // getwords(currentWords,fin); //Number of contracted shells
  // getwords(currentWords,fin); //Number of contracted shells
  // int nx=atoi(currentWords.back().c_str()); //Number of exponents
  int nx;
  notfound = true;
  while(notfound)
  {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++)
    {
      if("primitive" == currentWords[i])
      {
        nx=atoi(currentWords.back().c_str()); //Number of exponents
        notfound = false;
      }
    }
  }
  //allocate everything here
  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  gBound.resize(NumberOfAtoms+1);
  gShell.resize(ng);
  gNumber.resize(ng);
  gExp.resize(nx);
  gC0.resize(nx);
  gC1.resize(nx);
  // TDB: THIS ADDITION SHOULD BE COMPATIBLE WITH MY OLD FCHK FILES
  fin.seekg(pivottdb);//rewind it
  getGeometry(fin);
  std::cout << "Number of gaussians " << ng << std::endl;
  std::cout << "Number of primitives " << nx << std::endl;
  std::cout << "Number of atoms " << NumberOfAtoms << std::endl;
  search(fin, "Shell types");
  getGaussianCenters(fin);
  std::cout << " Shell types reading: OK" << std::endl;
// mmorales: SizeOfBasisSet > numMO always, so leave it like this
  EigVal_alpha.resize(SizeOfBasisSet);
  EigVal_beta.resize(SizeOfBasisSet);
// mmorales HACK HACK HACK, look for a way to rewind w/o closing/opening a file
  SpinRestricted = !(lookFor(fin, "Beta MO"));
  fin.close();
  fin.open(fname.c_str());
  search(fin, "Alpha Orbital"); //search "Alpha Orbital Energies"
// only read NumOfIndOrb
  getValues(fin,EigVal_alpha.begin(), EigVal_alpha.begin()+numMO);
  std::cout << " Orbital energies reading: OK" << std::endl;
  if(SpinRestricted)
  {
    EigVec.resize(2*numMO*SizeOfBasisSet);
    EigVal_beta=EigVal_alpha;
    std::vector<value_type> etemp;
    search(fin, "Alpha MO");
    getValues(fin,EigVec.begin(), EigVec.begin()+SizeOfBasisSet*numMO);
    copy(EigVec.begin(),EigVec.begin()+SizeOfBasisSet*numMO,EigVec.begin()+SizeOfBasisSet*numMO);
    std::cout << " Orbital coefficients reading: OK" << std::endl;
  }
  else
  {
    EigVec.resize(2*numMO*SizeOfBasisSet);
    std::vector<value_type> etemp;
    search(fin, "Beta Orbital");
    getValues(fin,EigVal_beta.begin(), EigVal_beta.begin()+numMO);
    std::cout << " Read Beta Orbital energies: OK" << std::endl;
    search(fin, "Alpha MO");
    getValues(fin,EigVec.begin(), EigVec.begin()+SizeOfBasisSet*numMO);
    search(fin, "Beta MO");
    getValues(fin,EigVec.begin()+numMO*SizeOfBasisSet, EigVec.begin()+2*numMO*SizeOfBasisSet);
    std::cout << " Alpha and Beta Orbital coefficients reading: OK" << std::endl;
  }
  if(multideterminant)
  {
    std::ifstream ofile(outputFile.c_str());
    if(ofile.fail())
    {
      std::cerr <<"Failed to open output file from gaussian. \n";
      exit(401);
    }
    int statesType = 0;
    std::streampos beginpos = ofile.tellg();
    bool found = lookFor(ofile, "SLATER DETERMINANT BASIS");
    if(!found)
    {
      ofile.close();
      ofile.open(outputFile.c_str());
      ofile.seekg(beginpos);
      found = lookFor(ofile, "Slater Determinants");
      if(found)
      {
        statesType=1;
        ofile.close();
        ofile.open(outputFile.c_str());
        ofile.seekg(beginpos);
        search(ofile, "Total number of active electrons");
        getwords(currentWords,ofile);
        ci_nstates = atoi(currentWords[5].c_str());
        search(ofile, "Slater Determinants");
      }
      else
      {
        std::cerr <<"Gaussian parser currenly works for slater determinant basis only. Use SlaterDet in CAS() or improve parser.\n";
        abort();
      }
    }
    beginpos = ofile.tellg();
// mmorales: gaussian by default prints only 50 states
// if you want more you need to use, iop(5/33=1), but this prints
// a lot of things. So far, I don't know how to change this
// without modifying the source code, l510.F or utilam.F
    found = lookFor(ofile, "Do an extra-iteration for final printing");
    std::map<int,int> coeff2confg;
//   290 FORMAT(1X,7('(',I5,')',F10.7,1X)/(1X,7('(',I5,')',F10.7,1X)))
    if(found)
    {
// this is tricky, because output depends on diagonalizer used (Davidson or Lanczos) for now hope this works
      std::streampos tmppos = ofile.tellg();
// long output with iop(5/33=1)
      int coeffType = 0;
      if(lookFor(ofile,"EIGENVECTOR USED TO COMPUTE DENSITY MATRICES SET"))
      {
        coeffType=1;
      }
      else
      {
// short list (50 states) with either Dav or Lanc
        ofile.close();
        ofile.open(outputFile.c_str());
        ofile.seekg(tmppos);//rewind it
        if(!lookFor(ofile,"EIGENVALUE "))
        {
          ofile.close();
          ofile.open(outputFile.c_str());
          ofile.seekg(tmppos);//rewind it
          if(!lookFor(ofile,"Eigenvalue"))
          {
            std::cerr <<"Failed to find CI voefficients.\n";
            abort();
          }
        }
      }
      std::streampos strpos = ofile.tellg();
      ofile.close();
      ofile.open(outputFile.c_str());
      ofile.seekg(beginpos);//rewind it
      std::string aline;
      int numCI;
      if(lookFor(ofile,"NO OF BASIS FUNCTIONS",aline))
      {
        parsewords(aline.c_str(),currentWords);
        numCI=atoi(currentWords[5].c_str());
      }
      else
      {
        ofile.close();
        ofile.open(outputFile.c_str());
        if(lookFor(ofile,"Number of configurations",aline))
        {
          parsewords(aline.c_str(),currentWords);
          numCI=atoi(currentWords[3].c_str());
        }
        else
        {
          std::cerr <<"Problems finding total number of configurations. \n";
          abort();
        }
      }
      std::cout <<"Total number of CI configurations in file (w/o truncation): " <<numCI << std::endl;
      ofile.close();
      ofile.open(outputFile.c_str());
      ofile.seekg(strpos);//rewind it
      int cnt=0;
      CIcoeff.clear();
      CIalpha.clear();
      CIbeta.clear();
      if(coeffType == 0)
      {
// short list
        for(int i=0; i<7; i++)
        {
          int pos=2;
          std::string aline;
          getline(ofile,aline,'\n');
//cout<<"aline: " <<aline << std::endl;
          for(int i=0; i<7; i++)
          {
            //int q = atoi( (aline.substr(pos,pos+4)).c_str() );
            //coeff2confg[q] = cnt++;
            //CIcoeff.push_back( atof( (aline.substr(pos+6,pos+15)).c_str() ) );
            //pos+=18;
            int q = atoi( (aline.substr(pos,pos+7)).c_str() );
            coeff2confg[q] = cnt++;
            CIcoeff.push_back( atof( (aline.substr(pos+9,pos+18)).c_str() ) );
            pos+=21;
//cout<<"confg, coeff: " <<q <<"  " <<CIcoeff.back() << std::endl;
          }
        }
        {
          int pos=2;
          std::string aline;
          getline(ofile,aline,'\n');
          //int q = atoi( (aline.substr(pos,pos+4)).c_str() );
          //coeff2confg[q] = cnt++;
          //CIcoeff.push_back( atof( (aline.substr(pos+6,pos+15)).c_str() ) );
          int q = atoi( (aline.substr(pos,pos+7)).c_str() );
          coeff2confg[q] = cnt++;
          CIcoeff.push_back( atof( (aline.substr(pos+9,pos+18)).c_str() ) );
//cout<<"confg, coeff: " <<q <<"  " <<CIcoeff.back() << std::endl;
        }
      }
      else
      {
// long list
        int nrows = numCI/5;
        int nextra = numCI - nrows*5;
        int indx[5];
        ofile.close();
        ofile.open(outputFile.c_str());
        ofile.seekg(strpos);//rewind it
        for(int i=0; i<nrows; i++)
        {
          getwords(currentWords,ofile);
          if(currentWords.size() != 5)
          {
            std::cerr <<"Error reading CI configurations-line: " <<i <<"\n";
            abort();
          }
          for(int k=0; k<5; k++)
            indx[k]=atoi(currentWords[k].c_str());
          getwords(currentWords,ofile);
          if(currentWords.size() != 6)
          {
            std::cerr <<"Error reading CI configurations-line: " <<i <<"\n";
            abort();
          }
          for(int k=0; k<5; k++)
          {
            // HACK HACK HACK - how is this done formally???
            for(int j=0; j<currentWords[k+1].size(); j++)
              if(currentWords[k+1][j] == 'D')
                currentWords[k+1][j]='E';
            double ci = atof(currentWords[k+1].c_str());
            if(std::abs(ci) > ci_threshold )
            {
              coeff2confg[indx[k]] = cnt++;
              CIcoeff.push_back(ci);
//cout<<"ind,cnt,c: " <<indx[k] <<"  " <<cnt-1 <<"  " <<CIcoeff.back() << std::endl;
            }
          }
        }
        if(nextra>0)
        {
          getwords(currentWords,ofile);
          if(currentWords.size() != nextra)
          {
            std::cerr <<"Error reading CI configurations last line \n";
            abort();
          }
          for(int k=0; k<nextra; k++)
            indx[k]=atoi(currentWords[k].c_str());
          getwords(currentWords,ofile);
          if(currentWords.size() != nextra+1)
          {
            std::cerr <<"Error reading CI configurations last line \n";
            abort();
          }
          for(int k=0; k<nextra; k++)
          {
            double ci = atof(currentWords[k+1].c_str());
            if(std::abs(ci) > ci_threshold )
            {
              coeff2confg[indx[k]] = cnt++;
              CIcoeff.push_back(ci);
            }
          }
        }
      }
      std::cout <<"Found " <<CIcoeff.size() <<" coeffficients after truncation. \n";
      if(statesType==0)
      {
        ofile.close();
        ofile.open(outputFile.c_str());
        ofile.seekg(beginpos);//rewind it
// this might not work, look for better entry point later
        search(ofile,"Truncation Level=");
//cout<<"found Truncation Level=" << std::endl;
        getwords(currentWords,ofile);
        while(!ofile.eof() && (currentWords[0] != "no." || currentWords[1] != "active" || currentWords[2] != "orbitals") )
        {
//        std::cout <<"1. " <<currentWords[0] << std::endl;
          getwords(currentWords,ofile);
        }
        ci_nstates=atoi(currentWords[4].c_str());
        getwords(currentWords,ofile);
// can choose specific irreps if I want...
        while(currentWords[0] != "Configuration" || currentWords[2] != "Symmetry" )
        {
//        std::cout <<"2. " <<currentWords[0] << std::endl;
          getwords(currentWords,ofile);
        }
        CIbeta.resize(CIcoeff.size());
        CIalpha.resize(CIcoeff.size());
        bool done=false;
// can choose specific irreps if I want...
        while(currentWords[0] == "Configuration" && currentWords[2] == "Symmetry" )
        {
          int pos = atoi(currentWords[1].c_str());
          std::map<int,int>::iterator it = coeff2confg.find( pos );
//cout<<"3. configuration: " <<currentWords[1].c_str() << std::endl;
          if(it != coeff2confg.end())
          {
            std::string alp(currentWords[4]);
            std::string beta(currentWords[4]);
            if(alp.size() != ci_nstates)
            {
              std::cerr <<"Problem with ci string. \n";
              abort();
            }
            for(int i=0; i<alp.size(); i++)
            {
              if(alp[i] == 'a')
                alp[i]='1';
              if(alp[i] == 'b')
                alp[i]='0';
              if(beta[i] == 'a')
                beta[i]='0';
              if(beta[i] == 'b')
                beta[i]='1';
            }
            if(done)
            {
              // check number of alpha/beta electrons
              int n1=0;
              for(int i=0; i<alp.size(); i++)
                if(alp[i] == '1')
                  n1++;
              if(n1 != ci_nea)
              {
                std::cerr <<"Problems with alpha ci string: "
                    << std::endl <<alp << std::endl <<currentWords[3] << std::endl;
                abort();
              }
              n1=0;
              for(int i=0; i<beta.size(); i++)
                if(beta[i] == '1')
                  n1++;
              if(n1 != ci_neb)
              {
                std::cerr <<"Problems with beta ci string: "
                    << std::endl <<beta << std::endl <<currentWords[3] << std::endl;
                abort();
              }
            }
            else
            {
              // count number of alpha/beta electrons
              ci_nea=0;
              for(int i=0; i<alp.size(); i++)
                if(alp[i] == '1')
                  ci_nea++;
              ci_neb=0;
              for(int i=0; i<beta.size(); i++)
                if(beta[i] == '1')
                  ci_neb++;
              ci_nca = nup-ci_nea;
              ci_ncb = ndown-ci_neb;
            }
            CIalpha[it->second] = alp;
            CIbeta[it->second] = beta;
//cout<<"alpha: " <<alp <<"  -   "  <<CIalpha[it->second] << std::endl;
          }
          getwords(currentWords,ofile);
        }
//cout.flush();
      }
      else
      {
        // coefficient list obtained with iop(4/21=110)
        ofile.close();
        ofile.open(outputFile.c_str());
        ofile.seekg(beginpos);//rewind it
        bool first_alpha = true;
        bool first_beta = true;
        std::string aline;
        getline(ofile,aline,'\n');
        CIbeta.resize(CIcoeff.size());
        CIalpha.resize(CIcoeff.size());
        for(int nst=0; nst<numCI; nst++)
        {
          getwords(currentWords,ofile);  // state number
          int pos = atoi(currentWords[0].c_str());
          std::map<int,int>::iterator it = coeff2confg.find( pos );
          if(it != coeff2confg.end())
          {
            getwords(currentWords,ofile);  // state number
            if(first_alpha)
            {
              first_alpha=false;
              ci_nea = currentWords.size();
              ci_nca = nup-ci_nea;
              std::cout <<"nca, nea, nstates: " <<ci_nca <<"  " <<ci_nea <<"  " <<ci_nstates << std::endl;
            }
            else
            {
              if(currentWords.size() != ci_nea)
              {
                std::cerr <<"Problems with alpha string: " <<pos << std::endl;
                abort();
              }
            }
            std::string alp(ci_nstates,'0');
            for(int i=0; i<currentWords.size(); i++)
            {
//              std::cout <<"i, alpOcc: " <<i <<"  " <<atoi(currentWords[i].c_str())-1 << std::endl; std::cout.flush();
              alp[ atoi(currentWords[i].c_str())-1 ] = '1';
            }
            getwords(currentWords,ofile);  // state number
            if(first_beta)
            {
              first_beta=false;
              ci_neb = currentWords.size();
              ci_ncb = ndown-ci_neb;
              std::cout <<"ncb, neb, nstates: " <<ci_ncb <<"  " <<ci_neb <<"  " <<ci_nstates << std::endl;
            }
            else
            {
              if(currentWords.size() != ci_neb)
              {
                std::cerr <<"Problems with beta string: " <<pos << std::endl;
                abort();
              }
            }
            std::string beta(ci_nstates,'0');
            for(int i=0; i<currentWords.size(); i++)
            {
//              std::cout <<"i, alpOcc: " <<i <<"  " <<atoi(currentWords[i].c_str())-1 << std::endl; std::cout.flush();
              beta[ atoi(currentWords[i].c_str())-1 ] = '1';
            }
            CIalpha[it->second] = alp;
            CIbeta[it->second] = beta;
//cout<<"alpha: " <<alp << std::endl <<"beta: " <<beta << std::endl;
            std::string aline1;
            getline(ofile,aline1,'\n');
          }
          else
          {
            getwords(currentWords,ofile);
            getwords(currentWords,ofile);
            std::string aline1;
            getline(ofile,aline1,'\n');
          }
        }
      }
      ofile.close();
    }
    else
    {
      std::cerr <<"Could not find CI coefficients in gaussian output file. \n";
      abort();
    }
    std::cout <<" size of CIalpha,CIbeta: " <<CIalpha.size() <<"  " <<CIbeta.size() << std::endl;
  }
}

void GaussianFCHKParser::getGeometry(std::istream& is)
{
  //atomic numbers
  std::vector<int> atomic_number(NumberOfAtoms);
  std::vector<double> q(NumberOfAtoms);
  //read atomic numbers
  search(is, "Atomic numbers");//search for Atomic numbers
  getValues(is,atomic_number.begin(),atomic_number.end());
  std::streampos pivot= is.tellg();
  //read effective nuclear charges
  search(is, "Nuclear");//search for Nuclear
  getValues(is,q.begin(),q.end());
  is.seekg(pivot);//rewind it
  search(is,"coordinates");
  std::vector<double> pos(NumberOfAtoms*OHMMS_DIM);
  getValues(is,pos.begin(),pos.end());
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

void GaussianFCHKParser::getGaussianCenters(std::istream& is)
{
  //map between Gaussian to Casino Shell notation
  std::map<int,int> gsMap;
  gsMap[0] =1; //s
  gsMap[-1]=2; //sp
  gsMap[1] =3; //p
  gsMap[-2]=4; //d
  gsMap[-3]=5; //f
  gsMap[-4]=6; //g
  gsMap[-5]=7; //l=5 h
  gsMap[-6]=8; //l=6 h1??
  gsMap[-7]=9; //l=7 h2??
  gsMap[-8]=10; //l=8 h3??
  std::vector<int> n(gShell.size()), dn(NumberOfAtoms,0);
  bool SPshell(false);
  getValues(is,n.begin(), n.end());
  for(int i=0; i<n.size(); i++)
  {
    gShell[i]=gsMap[n[i]];
    if(n[i] == -1)
      SPshell=true;
  }
  search(is, "Number");
  getValues(is,gNumber.begin(), gNumber.end());
  search(is, "Shell");
  getValues(is,n.begin(), n.end());
  for(int i=0; i<n.size(); i++)
    dn[n[i]-1]+=1;
  gBound[0]=0;
  for(int i=0; i<NumberOfAtoms; i++)
  {
    gBound[i+1]=gBound[i]+dn[i];
  }
  search(is, "Primitive");
  getValues(is,gExp.begin(), gExp.end());
  search(is, "Contraction");
  getValues(is,gC0.begin(), gC0.end());
  if(SPshell)
  {
    search(is, "P(S=P)");
    getValues(is,gC1.begin(), gC1.end());
  }
}

