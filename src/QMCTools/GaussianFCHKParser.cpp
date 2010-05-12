#include "QMCTools/GaussianFCHKParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

GaussianFCHKParser::GaussianFCHKParser() {
  basisName = "Gaussian-G2";
  Normalized = "no";
}

GaussianFCHKParser::GaussianFCHKParser(int argc, char** argv): 
  QMCGaussianParserBase(argc,argv) {
  basisName = "Gaussian-G2";
  Normalized = "no";
}

void GaussianFCHKParser::parse(const std::string& fname) {

  std::ifstream fin(fname.c_str());

  getwords(currentWords,fin); //1
  Title = currentWords[0];

  getwords(currentWords,fin);//2  SP RHF Gen
  //if(currentWords[1]=="ROHF" || currentWords[1]=="UHF") {
  if(currentWords[1]=="UHF") {
// mmorales: this should be determined by the existence of "Beta MO", since
// there are many ways to get unrestricted runs without UHF (e.g. UMP2,UCCD,etc)
    SpinRestricted=false;
    std::cout << " Spin Unrestricted Calculation (UHF). " << endl;    
  } else {
    SpinRestricted=true;
    std::cout << " Spin Restricted Calculation (RHF). " << endl;    
  }

  getwords(currentWords,fin);//3  Number of atoms
  NumberOfAtoms = atoi(currentWords.back().c_str());

  // TDB: THIS FIX SHOULD BE COMPATIBLE WITH MY OLD FCHK FILES
  //getwords(currentWords,fin); //4 Charge
  bool notfound = true;
  while(notfound) {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++){
      if("Charge" == currentWords[i]){
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

  // TDB: THIS ADDITION SHOULD BE COMPATIBLE WITH MY OLD FCHK FILES 
  streampos pivottdb = fin.tellg();

  int ng; 
  notfound = true;
  while(notfound) {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++){
      if("contracted" == currentWords[i]){
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
  while(notfound) {
    std::string aline;
    getwords(currentWords,fin);
    for(int i=0; i<currentWords.size(); i++){
      if("primitive" == currentWords[i]){
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
  std::cout << " Shell types reading: OK" << endl;

// mmorales:
  EigVal_alpha.resize(SizeOfBasisSet);
  EigVal_beta.resize(SizeOfBasisSet);
  search(fin, "Alpha Orbital"); //search "Alpha Orbital Energies"
// only read NumOfIndOrb
  getValues(fin,EigVal_alpha.begin(), EigVal_alpha.begin()+NumOfIndOrb);
  std::cout << " Orbital energies reading: OK" << endl;
  if(SpinRestricted) {
    EigVec.resize(2*SizeOfBasisSet*SizeOfBasisSet);
    EigVal_beta=EigVal_alpha;
    vector<value_type> etemp;
    search(fin, "Alpha MO");

    getValues(fin,EigVec.begin(), EigVec.begin()+SizeOfBasisSet*NumOfIndOrb); 
    std::copy(EigVec.begin(),EigVec.begin()+SizeOfBasisSet*NumOfIndOrb,EigVec.begin()+SizeOfBasisSet*SizeOfBasisSet);
    std::cout << " Orbital coefficients reading: OK" << endl;
  }
  else {
    EigVec.resize(2*SizeOfBasisSet*SizeOfBasisSet);
    vector<value_type> etemp;
    search(fin, "Beta Orbital"); 
    getValues(fin,EigVal_beta.begin(), EigVal_beta.begin()+NumOfIndOrb);
    std::cout << " Read Beta Orbital energies: OK" << endl;

    search(fin, "Alpha MO");
    getValues(fin,EigVec.begin(), EigVec.begin()+SizeOfBasisSet*NumOfIndOrb); 

    search(fin, "Beta MO");
    getValues(fin,EigVec.begin()+SizeOfBasisSet*SizeOfBasisSet, EigVec.begin()+SizeOfBasisSet*SizeOfBasisSet+SizeOfBasisSet*NumOfIndOrb); 

    std::cout << " Alpha and Beta Orbital coefficients reading: OK" << endl;
  }

  if(multideterminant) {

    std::ifstream ofile(outputFile.c_str());
    if(ofile.fail()) {
      cerr<<"Failed to open output file from gaussian. \n";
      exit(401);
    }

    bool found = lookFor(ofile, "SLATER DETERMINANT BASIS");
    if(!found) {
      cerr<<"Gaussian parser currenly works for slater determinant basis only.\n";
      abort();
    }

// mmorales: gaussian by default prints only 50 states
// if you want more you need to use, iop(5/33=1), but thus prints
// a lot of things. This stuff is printed in l510.F, line ~93300
// so a simple modification can fix this, but is not general for everyone
    //found = lookFor(ofile, "EIGENVALUES AND  EIGENVECTORS OF CI MATRIX");
    found = lookFor(ofile, "Do an extra-iteration for final printing");
/*    if(!found) {
      ofile.close();
      ofile.open(outputFile.c_str());
      if(ofile.fail()) {
        cerr<<"Failed to open output file from gaussian a second time. \n";
        exit(401);
      }
      found = lookFor(ofile, "EIGENVALUES AND  EIGENVECTORS OF CI MATRIX");
    }
*/
    map<int,int> coeff2confg;
//   290 FORMAT(1X,7('(',I5,')',F10.7,1X)/(1X,7('(',I5,')',F10.7,1X)))
    if(found) {

      streampos tmppos = ofile.tellg();
      if(!lookFor(ofile,"EIGENVALUE")) {
        ofile.seekg(tmppos);//rewind it
        if(!lookFor(ofile,"Eigenvalue")) {
          cerr<<"Failed to find VI voefficients.\n";
          abort();
        }
      }
      int cnt=0;
      CIcoeff.clear();
      CIalpha.clear(); 
      CIbeta.clear(); 
      for(int i=0; i<7; i++) {
        int pos=2;
        std::string aline; 
        getline(ofile,aline,'\n');
        for(int i=0; i<7; i++) {
          int q = atoi( (aline.substr(pos,pos+4)).c_str() );  
          coeff2confg[q] = cnt++;  
          CIcoeff.push_back( atof( (aline.substr(pos+6,pos+15)).c_str() ) );
          pos+=18; 
//cout<<"confg, coeff: " <<q <<"  " <<CIcoeff.back() <<endl;
        }
      }
      {
        int pos=2;
        std::string aline;
        getline(ofile,aline,'\n');
        int q = atoi( (aline.substr(pos,pos+4)).c_str() );  
        coeff2confg[q] = cnt++;  
        CIcoeff.push_back( atof( (aline.substr(pos+6,pos+15)).c_str() ) );
//cout<<"confg, coeff: " <<q <<"  " <<CIcoeff.back() <<endl;
      } 
      ofile.close();  // can I rewind???
      ofile.open(outputFile.c_str());
      if(ofile.fail()) {
        cerr<<"Failed to open output file from gaussian. \n";
        abort();
      }

// this might not work, look for better entry point later
      search(ofile,"Truncation Level=");
//cout<<"found Truncation Level=" <<endl; 
      getwords(currentWords,ofile);
      while(!ofile.eof() && (currentWords[0] != "no." || currentWords[1] != "active" || currentWords[2] != "orbitals") )
      {
//        cout<<"1. " <<currentWords[0] <<endl;
        getwords(currentWords,ofile);
      }
      ci_nstates=atoi(currentWords[4].c_str());
     
      getwords(currentWords,ofile);
// can choose specific irreps if I want...
      while(currentWords[0] != "Configuration" || currentWords[2] != "Symmetry" ) 
      {
//        cout<<"2. " <<currentWords[0] <<endl;
        getwords(currentWords,ofile);
      }

      CIbeta.resize(CIcoeff.size());
      CIalpha.resize(CIcoeff.size());
     
      bool done=false;
// can choose specific irreps if I want...
      while(currentWords[0] == "Configuration" && currentWords[2] == "Symmetry" ) 
      {
         int pos = atoi(currentWords[1].c_str());
         map<int,int>::iterator it = coeff2confg.find( pos );

//cout<<"3. configuration: " <<currentWords[1].c_str() <<endl;
        
         if(it != coeff2confg.end()) {

           std::string alp(currentWords[4]);
           std::string beta(currentWords[4]);

           if(alp.size() != ci_nstates) {
             cerr<<"Problem with ci string. \n";
             abort();
           }
           
           for(int i=0; i<alp.size(); i++) {
             if(alp[i] == 'a') alp[i]='1';
             if(alp[i] == 'b') alp[i]='0';
             if(beta[i] == 'a') beta[i]='0';
             if(beta[i] == 'b') beta[i]='1';
           } 
           if(done) {
           // check number of alpha/beta electrons
             int n1=0;
             for(int i=0; i<alp.size(); i++)  
               if(alp[i] == '1') n1++;
             if(n1 != ci_nea) {
               cerr<<"Problems with alpha ci string: " 
                   <<endl <<alp <<endl <<currentWords[3] <<endl;
               abort();
             } 
             n1=0;
             for(int i=0; i<beta.size(); i++) 
               if(beta[i] == '1') n1++;
             if(n1 != ci_neb) {
               cerr<<"Problems with beta ci string: "
                   <<endl <<beta <<endl <<currentWords[3] <<endl;
               abort();
             } 
           } else {
           // count number of alpha/beta electrons
             ci_nea=0;
             for(int i=0; i<alp.size(); i++) 
               if(alp[i] == '1') ci_nea++;
             ci_neb=0;
             for(int i=0; i<beta.size(); i++)
               if(beta[i] == '1') ci_neb++;
             ci_nca = nup-ci_nea;
             ci_ncb = ndown-ci_neb;
           } 
           CIalpha[it->second] = alp; 
           CIbeta[it->second] = beta; 
//cout<<"alpha: " <<alp <<"  -   "  <<CIalpha[it->second] <<endl;
         } 
         getwords(currentWords,ofile);
      }
//cout.flush();
      ofile.close();
    } else {
// if iop(5/33=1), look for "FINAL EIGENVALUES AND EIGENVECTORS"
      cerr<<"Could not find CI coefficients in gaussian output file. \n";
      abort();
    }

  }

}

void GaussianFCHKParser::getGeometry(std::istream& is) {
  //atomic numbers
  vector<int> atomic_number(NumberOfAtoms);
  vector<double> q(NumberOfAtoms);

  //read atomic numbers
  search(is, "Atomic numbers");//search for Atomic numbers
  getValues(is,atomic_number.begin(),atomic_number.end());

  streampos pivot= is.tellg();
  //read effective nuclear charges
  search(is, "Nuclear");//search for Nuclear
  getValues(is,q.begin(),q.end());

  is.seekg(pivot);//rewind it
  search(is,"coordinates");
  vector<double> pos(NumberOfAtoms*OHMMS_DIM);
  getValues(is,pos.begin(),pos.end());

  SpeciesSet& species(IonSystem.getSpeciesSet());
  for(int i=0, ii=0; i<NumberOfAtoms; i++) {
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

void GaussianFCHKParser::getGaussianCenters(std::istream& is) {

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

  vector<int> n(gShell.size()), dn(NumberOfAtoms,0);

  bool SPshell(false);
  getValues(is,n.begin(), n.end());
  for(int i=0; i<n.size(); i++){
    gShell[i]=gsMap[n[i]];
    if(n[i] == -1)SPshell=true;
  }

  search(is, "Number");
  getValues(is,gNumber.begin(), gNumber.end());

  search(is, "Shell");
  getValues(is,n.begin(), n.end());
  for(int i=0; i<n.size(); i++) dn[n[i]-1]+=1;

  gBound[0]=0;
  for(int i=0; i<NumberOfAtoms; i++) {
    gBound[i+1]=gBound[i]+dn[i];
  }

  search(is, "Primitive");
  getValues(is,gExp.begin(), gExp.end());

  search(is, "Contraction");
  getValues(is,gC0.begin(), gC0.end());

  if(SPshell){
    search(is, "P(S=P)");
    getValues(is,gC1.begin(), gC1.end());
  }
}

