#include "QMCTools/VSVBParser.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

void Cartesian2Spherical(int n, double* Cart, double* Sphe);

VSVBParser::VSVBParser()
{
  basisName = "Gaussian-G2";
  Normalized = "no";
  usingECP=false;
  BohrUnit=false;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  NFZC=0;
  NumSpinCoupledOrbitals=0;
  NumOpenShellOrbitals=0;
  NumDoubleOccupiedOrbitals=0;
}

VSVBParser::VSVBParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "Gaussian-G2";
  Normalized = "no";
  usingECP=false;
  BohrUnit=false;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  NFZC=0;
  NumSpinCoupledOrbitals=0;
  NumOpenShellOrbitals=0;
  NumDoubleOccupiedOrbitals=0;
}

void VSVBParser::parse(const std::string& fname)
{
  
  std::ifstream fin(fname.c_str());
  pivot_begin= fin.tellg();
  std::string aline;
  bool MDVSVB=true;
  int NbVbStructures;
  
  search(fin,"The number of atoms or centers is",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfAtoms = atoi(currentWords[7].c_str());
  cout<<"NUMBER OF ATOMS: " <<NumberOfAtoms <<endl;

  search(fin,"basis functions",aline);
  parsewords(aline.c_str(),currentWords);
  SizeOfBasisSet = atoi(currentWords[2].c_str());
  cout<<" SizeOfBasisSet = "<<SizeOfBasisSet<<endl; 


  if(lookFor(fin,"Orbital information"))
  {
      search(fin,"spin-coupled orbitals labeled",aline);
      parsewords(aline.c_str(),currentWords);
      NumSpinCoupledOrbitals= atoi(currentWords[2].c_str());
      cout<<"Number of Spin Coupled Orbitals = "<<NumSpinCoupledOrbitals<<endl;

      search(fin,"open-shell orbitals",aline);
      parsewords(aline.c_str(),currentWords);
      NumOpenShellOrbitals= atoi(currentWords[2].c_str());
      cout<<"Number of Open Shell Orbitals = "<<NumOpenShellOrbitals<<endl;

      search(fin,"double-occupied orbital",aline);
      parsewords(aline.c_str(),currentWords);
      NumDoubleOccupiedOrbitals= atoi(currentWords[2].c_str());
      cout<<"Number of doubly Occupied Orbitals = "<<NumDoubleOccupiedOrbitals<<endl;
  }
  else{
     cout<<"Could not find Orbital information!!!! Check VSVB output!"<<endl;
     abort();
  }

  numMO=NumSpinCoupledOrbitals+NumOpenShellOrbitals+NumDoubleOccupiedOrbitals;
  cout<<"NUMBER OF MOs: " <<numMO <<endl;
  SpinMultiplicity = NumOpenShellOrbitals+1; 
  cout<<"SPIN MULTIPLICITY: " <<SpinMultiplicity <<endl;

  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);

  fin.close();
  fin.open(fname.c_str());
  pivot_begin= fin.tellg();
  getGeometry(fin);
  fin.seekg(pivot_begin);

  search(fin,"The total number of electrons is",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfEls = atoi(currentWords[6].c_str());
  cout<<"Number of electrons: " <<NumberOfEls <<endl;
  /*search(fin,"NUMBER OF OCCUPIED ORBITALS (ALPHA)",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfAlpha = atoi(currentWords[5].c_str());
  cout<<"Number of alpha electrons: " <<NumberOfAlpha <<endl;
  search(fin,"NUMBER OF OCCUPIED ORBITALS (BETA )",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfBeta = atoi(currentWords[6].c_str());
  cout<<"Number of beta electrons: " <<NumberOfBeta <<endl;*/

  NumberOfAlpha = NumberOfEls/2; 
  NumberOfBeta =  NumberOfEls/2; 
  cout<<"Number of alpha electrons: " <<NumberOfAlpha <<endl;
  cout<<"Number of beta electrons: " <<NumberOfBeta <<endl;

  getGaussianCenters(fin);
  getMO(fin);

  fin.seekg(pivot_begin);

  
  if(lookFor(fin,"This is a single-determinant VSVB wave function"))
  {
      MDVSVB=false;
      multideterminant=false;
      cout<<"This is a single-determinant VSVB calculation!"<<endl<<"CONVERSION TO QMCPACK INPUT COMPLETE!!"<<endl;
  }
  else
  {
      MDVSVB=true;
      multideterminant=true;
      cout<<"This is a multi-determinant VSVB calculation!"<<endl<<"Adding the VB structures to the determinant"<<endl;

  }
  fin.close();

  if(MDVSVB)
  {
    fin.open(fname.c_str());
    pivot_begin= fin.tellg();
    search(fin,"Determinant information");
    if(search(fin,"VB structures",aline))
    {
      getwords(currentWords,fin);
      parsewords(aline.c_str(),currentWords);
      cout<<"Found "<<currentWords[2]<<" VB structures"<<endl;
      NbVbStructures=atoi(currentWords[2].c_str());
      NFZC = 0; 
      NAC = numMO; 
      NEXT =0; 
      NTOT=NEXT+NAC;
      cout<<"# core, #active, #external: " <<NFZC <<" " <<NAC <<" " <<NEXT <<endl;
      //fin.seekg(pivot_begin);
      getMDVSVB(fin,NbVbStructures);
    }
    else
    {
       cout<<" Could not find VB structures!! Check VSVB output!"<<endl;
       abort();
    }

  }
    fin.close();

}


void VSVBParser::getGeometry(std::istream& is)
{
  std::string aline;
  bool all_atoms=true;
  //atomic numbers
  vector<int> atomic_number,core;
  vector<double> q,pos;
  int natms=0;
  int test=0;
  tags.clear();
  //is.seekg(pivot_begin);
  //read atomic info
  if(search(is,"The geometry (in Angstroms) is",aline))
  {
      getwords(currentWords,is);
      do
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
        getwords(currentWords,is);
        if(currentWords[0] == "end" &&
           currentWords[1] == "of" &&
           currentWords[2] == "the" &&
           currentWords[3] == "geometry"  )
        {
           all_atoms=false;
        }
      }
      while(all_atoms);
   }
   else
   {
    cerr<<"Could not find atomic coordinates. \n";
    abort();
   }

  if(natms != NumberOfAtoms)
  {
    cerr<<"Could not find atomic coordinates for all atoms. \n";
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
}

void VSVBParser::getGaussianCenters(std::istream& is)
{
  string Shell_temp;
  gBound.resize(NumberOfAtoms+1);
  int ng,nx;
  string aline;
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
  vector<vector<double> > expo(nUniqAt),coef(nUniqAt),coef2(nUniqAt);
  vector<int> nshll(nUniqAt,0); //use this to 
  vector<vector<int> > ncoeffpershell(nUniqAt);
  vector<vector<std::string> > shID(nUniqAt);
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
      cerr<<"Problem with basis set data.\n";
      abort();
    }
    getwords(currentWords,is);
    if(currentWords.size() < 5)
      continue;
    if(currentWords[0] == "Beginning" &&
        currentWords[1] == "of" &&
        currentWords[2] == "the" &&
        currentWords[3] == "basis" &&
        currentWords[4] == "set") 
        found=true; 
  }


  is.seekg(pivot_begin);


  int currPos=-1;
  lookFor(is,"Beginning of the basis set"); 
 int NbCoeffperShell=0;
 getwords(currentWords,is);
 bool allbase=true;
 while(allbase==true)
 {
    if(currentWords.empty()) getwords(currentWords,is);;
    if(currentWords[0] == "The" &&
        currentWords[1] == "basis" &&
        currentWords[2] == "set" &&
        currentWords[3] == "for" &&
        currentWords[4] == "atom" &&
        currentWords[6] == "is") //FOUND SPECIES
    {
        std::map<std::string,int>::iterator it(basisDataMap.find(currentWords[5]));
        if(it == basisDataMap.end())
        {
           cerr<<"Error in parser.\n";
           abort();
        }
        currPos=it->second;
        nshll[currPos]=0;
        ncoeffpershell[currPos].clear();
        ncoeffpershell[currPos].push_back(0);
        shID[currPos].clear();
        shID[currPos].push_back("NONE");

        bool group=true; 
        do{ 
           getwords(currentWords,is);
           if(currentWords[0]=="S"||currentWords[0]=="P"||currentWords[0]=="D"||currentWords[0]=="F"||currentWords[0]=="G"||currentWords[0]=="H"||currentWords[0]=="I") //empty line after the shell
           {
              shID[currPos].push_back(currentWords[0]);
              Shell_temp=currentWords[0];
              NbCoeffperShell=atoi(currentWords[1].c_str());
              ncoeffpershell[currPos][nshll[currPos]]=NbCoeffperShell;
              for(int nbcoef=0;nbcoef<NbCoeffperShell;nbcoef++){
                 getwords(currentWords,is);
                 expo[currPos].push_back(atof(currentWords[0].c_str()));
                 coef[currPos].push_back(atof(currentWords[1].c_str()));
                 shID[currPos][nshll[currPos]] = Shell_temp; 
                 cout << currPos << ":" <<expo[currPos].back() << " " << coef[currPos].back() << " " 
                 << ncoeffpershell[currPos][nshll[currPos]] 
                 << " " << shID[currPos][nshll[currPos]]<<endl; 
              }              
            nshll[currPos]++;
            ncoeffpershell[currPos].push_back(0);
            shID[currPos].push_back("NONE");

           }
           else{
              if(currentWords[0] == "end" && currentWords[1] == "of")
              {
                 ng=SizeOfBasisSet;
                 group=false;
                 allbase=false;              
                 break;
              }
              else
              {
                 break;
              }
           }
        }while (group==true);
      }
      else
      {
          cerr<<"error in parser"<<endl;
          abort();
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
      cerr<<"Error in parser.\n";
      abort();
    }
    gBound[i] = gtot;
    int indx = it->second;
    gtot+=nshll[indx];
    for(int k=0; k<nshll[indx]; k++){
      gShell.push_back(gsMap[shID[indx][k]]);

    }
    for(int k=0; k<nshll[indx]; k++)
      gNumber.push_back(ncoeffpershell[indx][k]);
    for(int k=0; k<expo[indx].size(); k++)
      gExp.push_back(expo[indx][k]);
    for(int k=0; k<coef[indx].size(); k++)
      gC0.push_back(coef[indx][k]);
  }
  gBound[NumberOfAtoms] = gtot;
}

void VSVBParser::getMO(std::istream& is)
{
  int nq = SizeOfBasisSet/4;
  int rem = SizeOfBasisSet%4;
  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  EigVec.resize(2*SizeOfBasisSet*numMO);
  std::string aline;
  vector <double> SpinCoupledOrbitals;
  vector <double> OpenShellOrbitals;
  vector <double> DoubleOccupiedOrbitals;

  //cout<<"Size of BasisSet="<<SizeOfBasisSet<<" Size of nq="<<nq<<"  Size of rem="<<rem<<" nq+rm="<<4*nq+rem<<endl; 
  search(is,"Orbital information ...");
  getwords(currentWords,is); //Number of spin-coupled orbitals
  if(NumSpinCoupledOrbitals !=0){
     getwords(currentWords,is); // Orbital 0
     for (int i=0; i<NumSpinCoupledOrbitals; i++){
        for (int j=0;j<nq; j++){
           getwords(currentWords,is);
           for (int k=0;k<4;k++)
              SpinCoupledOrbitals.push_back(atof(currentWords[k].c_str()));
           
        }
        if(rem>0){
           getwords(currentWords,is);
           for (int k=0;k<rem;k++)
              SpinCoupledOrbitals.push_back(atof(currentWords[k].c_str()));
        }
        getwords(currentWords,is);//end of Orbital
        getwords(currentWords,is);
     }
  }
  cout<<"Size of elements for spin-coupled orbitals ="<<SpinCoupledOrbitals.size()<<endl;
  
  getwords(currentWords,is); //Empty line
  getwords(currentWords,is); //Number of open shell orbitals 

  if(NumOpenShellOrbitals !=0){
     getwords(currentWords,is); // Orbital 0
     for (int i=0; i<NumOpenShellOrbitals; i++){
        for (int j=0;j<nq; j++){
           getwords(currentWords,is);
           for (int k=0;k<4;k++)
              OpenShellOrbitals.push_back(atof(currentWords[k].c_str()));

        }
        if(rem>0){
           getwords(currentWords,is);
           for (int k=0;k<rem;k++)
              OpenShellOrbitals.push_back(atof(currentWords[k].c_str()));
        } 
        getwords(currentWords,is);
        getwords(currentWords,is);
     }
  }
  cout<<"Size of elements for Open Shell orbitals="<<OpenShellOrbitals.size()<<endl;

  getwords(currentWords,is); //Empty line
  getwords(currentWords,is); //Number of doubly occupied orbitals 

  if(NumDoubleOccupiedOrbitals !=0){
     getwords(currentWords,is); // Orbital 0
     for (int i=0; i<NumDoubleOccupiedOrbitals; i++){
        for (int j=0;j<nq; j++){
           getwords(currentWords,is);
           for (int k=0;k<4;k++)
              DoubleOccupiedOrbitals.push_back(atof(currentWords[k].c_str()));

        }
        if(rem>0){
           getwords(currentWords,is);
           for (int k=0;k<rem;k++)
              DoubleOccupiedOrbitals.push_back(atof(currentWords[k].c_str()));
        } 
        getwords(currentWords,is);
        getwords(currentWords,is);
     }
  }

  cout<<"Size of elements for Doubly Occupied orbitals="<<DoubleOccupiedOrbitals.size()<<endl;

    int cnt=0;
    for(int spin=0;spin<2;spin++){
       if(SpinCoupledOrbitals.size()>0){
          for (int i=0;i<SpinCoupledOrbitals.size();i++)
             EigVec[cnt++] = SpinCoupledOrbitals[i]; 
       }
   
       if(SpinCoupledOrbitals.size()>0){
          for (int i=0;i<OpenShellOrbitals.size();i++)
             EigVec[cnt++] = OpenShellOrbitals[i]; 
       }
   
       if(DoubleOccupiedOrbitals.size()>0){
          for (int i=0;i<DoubleOccupiedOrbitals.size();i++)
             EigVec[cnt++] = DoubleOccupiedOrbitals[i]; 
       }
   
    }
  cout<<"Finished reading MO." <<endl;
}

void VSVBParser::getMDVSVB(std::istream& is,int NbVbStructures)
{
  int VB_indx;
  double VB_Coeff; 

  ci_size=0;
  CSFocc.clear();
  CSFalpha.clear();
  CSFbeta.clear();
  coeff2csf.clear();
  CSFexpansion.clear(); 
  CSFexpansion.resize(NbVbStructures);
  CSFalpha.resize(NbVbStructures);
  CSFbeta.resize(NbVbStructures);
  
  for (int numVB=0; numVB<NbVbStructures;numVB++)
  {

    getwords(currentWords,is);
    if(currentWords.empty()) getwords(currentWords,is);
    if(currentWords[0] == "end" &&
        currentWords[1] == "of" &&
        currentWords[2] == "determinant" &&
        currentWords[3] == "information" )
    {
        cout<<"Done reading VB determinants"<<endl;
        break;
    }

    VB_Coeff=atof(currentWords[currentWords.size()-1].c_str());
    pair<int,double> cic(numVB,VB_Coeff);
    coeff2csf.push_back(cic);
    CSFocc.push_back(getOccup(1));
    VB_indx=atoi(currentWords[2].c_str());
    bool Val=true;
    while(Val==true){
       getwords(currentWords,is);
       if(currentWords[0] == "end" &&
           currentWords[1] == "of" &&
           currentWords[2] == "VB" &&
           currentWords[3] == "structure" )
       {
           cout<<"Done reading CI for VB "<<VB_indx<<endl;
           Val=false;
           break;
       }
       CSFexpansion[numVB].push_back(atof(currentWords[2].c_str()));
       //Alpha spin
       getwords(currentWords,is);
       getwords(currentWords,is);
       CSFalpha[numVB].push_back(getOccup(0));
       //Beta spin
       getwords(currentWords,is);
       getwords(currentWords,is);
       CSFbeta[numVB].push_back(getOccup(0));

       
    } 
    if(is.eof())
    {
       cerr<<"Could not find VB structures!. \n";
       abort();
    }
    
  }
  cout<<" Done reading VB structures!!"<<endl;
  ci_size=NbVbStructures;
  ci_nstates = numMO;
  int ds=SpinMultiplicity-1;
  ci_neb= (NumberOfEls-ds)/2;
  ci_nea= NumberOfEls-NumberOfBeta;
  ci_nca = 0;
  ci_ncb = 0;
}


/*
void VSVBParser::getMDVSVB(std::istream& is,int NbVbStructures)
{
  //is.seekg(pivot_begin);
  //look for CI coefficients
  bool notfound=true;
  string aline;
  int VB_indx;
  int Count=0;
  double VB_Coeff; 
  ci_size=0;
  CIcoeff.clear();
  CIalpha.clear();
  CIbeta.clear();

  

  for (int numVB=0; numVB<NbVbStructures;numVB++)
  {
    getwords(currentWords,is);
    if(currentWords.empty()) getwords(currentWords,is);
    if(currentWords[0] == "end" &&
        currentWords[1] == "of" &&
        currentWords[2] == "determinant" &&
        currentWords[3] == "information" )
    {
        cout<<"Done reading VB determinants"<<endl;
        break;
    }
    VB_Coeff=atof(currentWords[currentWords.size()-1].c_str());
    VB_indx=atoi(currentWords[2].c_str());
    bool Val=true;
    while(Val==true){
       getwords(currentWords,is);
       if(currentWords[0] == "end" &&
           currentWords[1] == "of" &&
           currentWords[2] == "VB" &&
           currentWords[3] == "structure" )
       {
           cout<<"Done reading CI for VB "<<VB_indx<<endl;
           Val=false;
           break;
       }
       ci_size++;
       CIcoeff.push_back(atof(currentWords[2].c_str()));
       //Alpha spin
       getwords(currentWords,is);
       getwords(currentWords,is);
       CIalpha.push_back(getOccup());
       //Beta spin
       getwords(currentWords,is);
       getwords(currentWords,is);
       CIbeta.push_back(getOccup());
       
    } 
    if(is.eof())
    {
       cerr<<"Could not find VB structures!. \n";
       abort();
    }
    
  }
  cout<<" Done reading VB structures!!"<<endl;
  ci_nea=ci_neb=0;
  
  for(int i=0; i<CIalpha[0].size(); i++)
    if(CIalpha[0].at(i) == '1')
      ci_nea++;
  for(int i=0; i<CIbeta[0].size(); i++)
    if(CIbeta[0].at(i) == '1')
      ci_neb++;
  if(CIalpha[0].size() != CIbeta[0].size())
  {
    cerr<<"QMCPack can't handle different number of active orbitals in alpha and beta channels right now. Contact developers for help (Miguel).\n";
    abort();
  }
  int ds=SpinMultiplicity-1;
  int neb= (NumberOfEls-ds)/2;
  int nea= NumberOfEls-NumberOfBeta;
  ci_nca = nea-ci_nea;
  ci_ncb = neb-ci_neb;
  ci_nstates = CIalpha[0].size();
}
*/
std::string VSVBParser::getOccup(int val)
{

    string Occup; 

    if (val==0){
       char Temp[numMO+1];
       for (int i=0;i<numMO;i++)
          Temp[i]='0';
      
       for (int i=0; i<currentWords.size();i++)
          Temp[atoi(currentWords[i].c_str())]='1';

       for (int i=0;i<numMO;i++)
          Occup+=Temp[i];
    }       
    else
    {
       for (int i=0;i<NumSpinCoupledOrbitals;i++)
          Occup+='1';
       for (int i=0;i<NumOpenShellOrbitals;i++)
          Occup+='1';
       for (int i=0;i<NumDoubleOccupiedOrbitals;i++)
          Occup+='2';
    }
    return Occup;
}
