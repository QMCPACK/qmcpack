#include "QMCTools/QPParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>

using namespace std;


QPParser::QPParser()
{
  basisName = "Gaussian-G2";
  Normalized = "no";
  usingECP=false;
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  readtype=0;
  NFZC=0;
  numAO=0;
  FixValence=true;
  
}

QPParser::QPParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "Gaussian-G2";
  Normalized = "no";
  usingECP=false;
  BohrUnit=true;
  MOtype="Canonical";
  angular_type="cartesian";
  SpinRestricted=true;
  readtype=0;
  NFZC=0;
  numAO=0;
  FixValence=true;
}

void QPParser::parse(const std::string& fname)
{
  std::ifstream fin(fname.c_str());
  pivot_begin= fin.tellg();
  std::string aline;

  search(fin,"do_pseudo",aline);
  parsewords(aline.c_str(),currentWords);
  if(currentWords[1]=="True")
     usingECP=true; 
  else
     usingECP=false;
  cout<<"usingECP: " <<(usingECP?("yes"):("no")) <<endl;
  cout.flush();

  search(fin,"multi_det",aline);
  parsewords(aline.c_str(),currentWords);
  if(currentWords[1]=="True")
     multideterminant=true;
  else
     multideterminant=false;

  cout<<"Multideterminant: " <<(multideterminant?("yes"):("no")) <<endl;



  search(fin,"ao_num",aline);
  parsewords(aline.c_str(),currentWords);
  numAO = atoi(currentWords[1].c_str());
  cout<<"NUMBER OF AOs: " <<numAO <<endl;
  SizeOfBasisSet = numAO; 
  cout<<"Size of Basis Set: " <<SizeOfBasisSet <<endl;

  search(fin,"mo_num",aline);
  parsewords(aline.c_str(),currentWords);
  numMO = atoi(currentWords[1].c_str());
  cout<<"NUMBER OF MOs: " <<numMO <<endl;





  search(fin,"elec_alpha_num",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfAlpha = atoi(currentWords[1].c_str());
  cout<<"Number of alpha electrons: " <<NumberOfAlpha <<endl;
  search(fin,"elec_beta_num",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfBeta = atoi(currentWords[1].c_str());
  cout<<"Number of beta electrons: " <<NumberOfBeta <<endl;


  search(fin,"elec_tot_num",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfEls = atoi(currentWords[1].c_str());
  cout<<"Number of electrons: " <<NumberOfEls <<endl;







  search(fin,"spin_multiplicity",aline);
  parsewords(aline.c_str(),currentWords);
  SpinMultiplicity = atoi(currentWords[1].c_str());
  cout<<"SPIN MULTIPLICITY: " <<SpinMultiplicity <<endl;

  SpinRestricted=true;


  search(fin,"nucl_num",aline);
  parsewords(aline.c_str(),currentWords);
  NumberOfAtoms = atoi(currentWords[1].c_str());
  cout<<"NUMBER OF ATOMS: " <<NumberOfAtoms <<endl;


  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  getGeometry(fin);
  fin.seekg(pivot_begin);

  getGaussianCenters(fin);
  fin.seekg(pivot_begin);
  MOtype = "Canonical";
  readtype=0;
  getMO(fin);
  fin.close();


  if(multideterminant)
  {
    fin.open(fname.c_str());
    QP=true;
    pivot_begin= fin.tellg();
    search(fin,"BEGIN_DET",aline);
    getwords(currentWords,fin);//Empty Line
    getwords(currentWords,fin);//mo_num
    getwords(currentWords,fin);//Number Of determinants
    cout<<"Found "<<currentWords[1]<<" Determinants"<<endl;
    ci_size=atoi(currentWords[1].c_str());
    NFZC = 0; 
    NAC = numMO; 
    NEXT =0; 
    NTOT=NEXT+NAC;
    cout<<"# core, #active, #external: " <<NFZC <<" " <<NAC <<" " <<NEXT <<endl;
    //fin.seekg(pivot_begin);
    getQPCI(fin);
  }

}

void QPParser::getGeometry(std::istream& is)
{
  //atomic numbers
  vector<int> atomic_number,core;
  vector<double> q,pos;
  int natms=0;
  tags.clear();
  is.seekg(pivot_begin);
  //read atomic info
  bool notfound=true;

  do
  {
    if(is.eof())
    {
      cerr<<"Could not find atomic coordinates. \n";
      abort();
    }
    getwords(currentWords,is);
    if(currentWords.size() < 4 )
      continue;
    if (currentWords.size() == 4)
    if(currentWords[0] == "Atomic" &&
        currentWords[1] == "coord" &&
        currentWords[2] == "in" &&
        currentWords[3] == "Bohr"  )
    {
      notfound=false;
      getwords(currentWords,is);
      while(currentWords.size() != 0)
      {
        if(currentWords[0] == "BEGIN_BASIS_SET")
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
  if(natms != NumberOfAtoms)
  {
    cerr<<"Could not find atomic coordinates for all atoms. \n";
    abort();
  }
//ECP PART!!!
  if(usingECP==true){
    is.seekg(pivot_begin);
    notfound=true;
 
    while(notfound)
    {
       if(is.eof())
       {
          cerr<<"Problem looking for ECPs, this should not happen. Contact developers for help. \n";
          abort();
       }
       getwords(currentWords,is);
// this should appear below the ECP section in the output file
// so use this to avoid going all the way to the bottom
       if(currentWords.size() == 0 )
           continue;
       if( currentWords[0] == "BEGIN_PSEUDO")
       {
           core.resize(NumberOfAtoms);
       //    getwords(currentWords,is); // -------------
           bool done=false;
           while(!done)
           {
              if(is.eof())
              {
                 cerr<<"2 Found ECPs, but problem looking ZCORE data.\n";
                 abort();
              }
              getwords(currentWords,is);
              if(currentWords.size() == 0)
                 continue;
              if(currentWords.size() == 1)
              {
                 if(currentWords[0] == "END_PSEUDO") 
                    done=true;
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
                    cerr<<"Problem with ECP data. Didn't found ATOM tag\n";
                    cerr<< is.rdbuf() <<endl;
                    abort();
                 }
                 it0++;
                 int nq0 = atoi(it0->c_str())-1;
                 if(it != currentWords.end())
                 {
                    it++;
                    cout<<"Avant nq0="<<nq0<<" core[nq0]="<<core[nq0]<<"   q[nq0]="<<q[nq0]<<endl;
                    core[nq0] = atoi(it->c_str());
                    q[nq0] -= core[nq0];
                       
                    cout<<"Apres nq0="<<nq0<<" core[nq0]="<<core[nq0]<<"   q[nq0]="<<q[nq0]<<endl;
                    cout<<"1 Found ECP for atom " <<nq0 <<" with zcore " <<core[nq0] <<endl;
                 }
                 else
                 {
                    it = find(currentWords.begin(),currentWords.end(),"ATOM");
                    if(it == currentWords.end())
                    {
                        cerr<<"Problem with ECP data. Didn't found ATOM tag\n";
                        cerr<<"Atom: " <<nq0 <<endl;
                        abort();
                    }
                    std::vector<std::string>::iterator it2=it;
                    it2++;
                    int nq = atoi(it2->c_str());
                    if(nq != nq0+1)
                    {
                        cerr<<"Problem with ECP data. ID's don't agree\n";
                        cerr<<"Atom: " <<nq0 <<endl;
                        abort();
                    }
                    it = find(it2,currentWords.end(),"ATOM");
                    if(it == currentWords.end())
                    {
                       cerr<<"Problem with ECP data (2).\n";
                       cerr<<"Atom: " << nq0 <<endl;
                       abort();
                    }
                    nq = atoi((it+1)->c_str());
                    core[nq0] = core[nq-1];
                    q[nq0] -= core[nq0];
                    cout<<"2 Found ECP for atom " <<nq0 <<" with zcore " <<core[nq0] <<endl;
                 }
              }
           }
           notfound=false;
       }
       else
       {
          if(currentWords.size() < 2 )
              continue;
          if( currentWords[0] == "END_PSEUDO")
               break;
       }
    }
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

void QPParser::getGaussianCenters(std::istream& is)
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
    if(currentWords.size() >2)
      continue;
    if(currentWords[0] == "BEGIN_BASIS_SET") 
        found=true; 
  }


  is.seekg(pivot_begin);


  int currPos=-1;
 lookFor(is,"BEGIN_BASIS_SET"); 
 int NbCoeffperShell=0;
 getwords(currentWords,is);
 bool allbase=true;
 while(allbase==true)
 {
    if(currentWords.empty()) getwords(currentWords,is);
    if(currentWords.size() ==1)
    {

        std::map<std::string,int>::iterator it(basisDataMap.find(currentWords[0]));
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
                 expo[currPos].push_back(atof(currentWords[1].c_str()));
                 coef[currPos].push_back(atof(currentWords[2].c_str()));
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
              if(currentWords[0] == "END_BASIS_SET")
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

void QPParser::getMO(std::istream& is)
{
  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  EigVec.resize(2*SizeOfBasisSet*numMO);
  std::string aline;
  search(is,"BEGIN_MO");
  getwords(currentWords,is);  // empty line
  
  vector<double> dummy(50);
  Matrix<double> CartMat(numMO,SizeOfBasisSet);
  streampos pivot;
  pivot= is.tellg();
  std::vector<std::string> CartLabels(SizeOfBasisSet);
  getwords(currentWords,is);
  for(int k=0; k<SizeOfBasisSet; k++)
  {
    getwords(currentWords,is);
    CartLabels[k] = currentWords[0];
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
    
    search(is,"BEGIN_MO");
    getwords(currentWords,is);  // empty line
    getMO_single_set(is, CartMat, EigVal_beta);
  }
    
  for(int i=0; i<numMO; i++)
    for(int k=0; k<SizeOfBasisSet; k++)
      EigVec[cnt++] = CartMat[i][k];
  cout<<"Finished reading MO." <<endl;
}

void QPParser::getMO_single_set(std::istream& is, Matrix<double> &CartMat, std::vector<value_type>& EigVal)
{
  int nq = numMO/4;
  int rem = numMO%4;
  int cnt=0;
  for(int i=0; i<nq; i++)
  {
    getwords(currentWords,is);
    for(int k=0; k<SizeOfBasisSet; k++)
    {
      getwords(currentWords,is);
      if(currentWords.size() == 6)
      {
        CartMat[cnt][k] = atof(currentWords[2].c_str()) ;
        CartMat[cnt+1][k] = atof(currentWords[3].c_str()) ;
        CartMat[cnt+2][k] = atof(currentWords[4].c_str()) ;
        CartMat[cnt+3][k] = atof(currentWords[5].c_str()) ;
      }
      else
      {
          cerr<<"Problem reading orbitals!!" <<endl;
          abort();
      }
    }
    getwords(currentWords,is);
    cnt+=4;
  }
  if(rem > 0)
  {
    getwords(currentWords,is);

    for(int k=0; k<SizeOfBasisSet; k++)
    {
      getwords(currentWords,is);
      for(int i=0; i<rem; i++)
      {
        CartMat[cnt+i][k] = atof(currentWords[2+i].c_str()) ;
      }
    }
    getwords(currentWords,is);
  }
}



void QPParser::getQPCI(std::istream& is)
{
  int Ci_indx;
  double Ci_Coeff; 

  CIcoeff.clear();
  CIalpha.clear();
  CIbeta.clear();
  CIcoeff.resize(ci_size);
  CIalpha.resize(ci_size);
  CIbeta.resize(ci_size);


  for (int i=0; i<ci_size;i++)
  {

    getwords(currentWords,is); //Empty Line
    if(currentWords[0] == "END_DET")
    {
        cout<<"Done reading determinants"<<endl;
        break;
    }

    getwords(currentWords,is); //Coeff
    CIcoeff[i]=atof(currentWords[0].c_str());

    //Alpha spin
    getwords(currentWords,is);
    CIalpha[i]=currentWords[0];
    //Beta spin
    getwords(currentWords,is);
    CIbeta[i]=currentWords[0];
  }
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
  cout<<" Done reading CIs!!"<<endl;
  ci_nstates = CIalpha[0].size();
}
