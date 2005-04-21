#include "QMCTools/QMCGaussianParserBase.h"
#include <iterator>
#include <algorithm>
#include <numeric>
#include <set>
#include <map>
using namespace std;
#include "QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.h"

//std::vector<std::string> QMCGaussianParserBase::IonName;
std::map<int,std::string> QMCGaussianParserBase::IonName;
std::vector<std::string> QMCGaussianParserBase::gShellType;
std::vector<int> QMCGaussianParserBase::gShellID;

QMCGaussianParserBase::QMCGaussianParserBase(): 
  Title("sample"),basisType("Gaussian"),basisName("generic"),
  Normalized("no"),gridPtr(0)
{
}

QMCGaussianParserBase::QMCGaussianParserBase(int argc, char** argv):
  Title("sample"),basisType("Gaussian"),basisName("generic"),  
  Normalized("no"),gridPtr(0)
{

  createGridNode(argc,argv);
}

void QMCGaussianParserBase::init() {
  //IonName.resize(24);
  IonName[1] = "H"; IonName[2] = "He"; IonName[3] = "Li";
  IonName[4] = "Be"; IonName[5] = "B"; IonName[6] = "C";
  IonName[7] = "N"; IonName[8] = "O"; IonName[9] = "F";
  IonName[10] = "Ne"; IonName[11] = "Na"; IonName[12] = "Mg";
  IonName[13] = "Al"; IonName[14] = "Si"; IonName[15] = "P";
  IonName[16] = "S"; IonName[17] = "Cl"; IonName[18] = "Ar";
  gShellType.resize(7);
  gShellType[1]="s"; gShellType[2]="sp"; gShellType[3]="p"; gShellType[4]="d"; gShellType[5]="f"; gShellType[6]="g";
  gShellID.resize(7);
  gShellID[1]=0; gShellID[2]=0; gShellID[3]=1; gShellID[4]=2; gShellID[5]=3; gShellID[6]=4;
}


xmlNodePtr QMCGaussianParserBase::createBasisSet() {

  //xmlNodePtr wf_root = xmlNewNode(NULL, BAD_CAST "determinantset");
  //xmlNewProp(wf_root,(const xmlChar*)"type",(const xmlChar*)"MolecularOrbital");

  xmlNodePtr bset = xmlNewNode(NULL,(const xmlChar*)"basisset");
  xmlNewProp(bset,(const xmlChar*)"ref",(const xmlChar*)"i");

  xmlNodePtr cur = xmlAddChild(bset,xmlNewNode(NULL,(const xmlChar*)"distancetable"));
  xmlNewProp(cur,(const xmlChar*)"source",(const xmlChar*)"i");
  xmlNewProp(cur,(const xmlChar*)"target",(const xmlChar*)"e");

  std::map<int,int> species;
  int gtot = 0;
  for(int iat=0; iat<NumberOfAtoms; iat++) {
    int itype = GroupID[iat];
    int ng = 0;
    std::map<int,int>::iterator it=species.find(itype);
    if(it == species.end()) {
      for(int ig=gBound[iat]; ig<gBound[iat+1]; ig++) {
        ng += gNumber[ig];
      }
      species[itype] = ng;
      cur = xmlAddSibling(cur,createCenter(iat,gtot));
      //xmlNodePtr abasis = createCenter(iat,gtot);
      //if(cur)
      //  cur = xmlAddSibling(cur,abasis);
      //else
      //  cur = xmlAddChild(bset,abasis);
    } else {
      ng = (*it).second;
    }
    gtot += ng;
  }

  return bset;
}

xmlNodePtr 
QMCGaussianParserBase::createDeterminantSet() {

  xmlNodePtr slaterdet = xmlNewNode(NULL,(const xmlChar*)"slaterdeterminant");

  //check spin-dependent properties
  int nup = NumberOfEls/2;
  int ndown = NumberOfEls-nup;
  std::ostringstream up_size, down_size, b_size, occ;
  up_size <<nup; down_size << ndown; b_size<<SizeOfBasisSet;

  //create a determinant Up
  xmlNodePtr adet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(adet,(const xmlChar*)"id",(const xmlChar*)"updet");
  xmlNewProp(adet,(const xmlChar*)"orbitals",(const xmlChar*)up_size.str().c_str());

  vector<int> iocc(SizeOfBasisSet,0);
  for(int j=0;j<nup; j++) iocc[j]=1;

  occ<<"\n";
  vector<int>::iterator it(iocc.begin()); 
  int i=0;
  while(i<SizeOfBasisSet) {
    int n = (i+10<SizeOfBasisSet)? 10 : SizeOfBasisSet-i;
    std::copy(it, it+n, ostream_iterator<int>(occ," "));
    occ << "\n"; it += 10; i+=10;
  }

  xmlNodePtr occ_data 
    = xmlNewTextChild(adet,NULL,(const xmlChar*)"occupation",(const xmlChar*)occ.str().c_str());
  xmlNewProp(occ_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());

  std::ostringstream eig;
  eig.setf(std::ios::scientific, std::ios::floatfield);
  eig.setf(std::ios::right,std::ios::adjustfield);
  eig.precision(14);
  int btot=SizeOfBasisSet*SizeOfBasisSet;
  int n=btot/4, b=0;
  int dn=btot-n*4;
  eig << "\n";
  for(int k=0; k<n; k++) {
    eig << setw(22) << EigVec[b] << setw(22) << EigVec[b+1] << setw(22) << EigVec[b+2] << setw(22) <<  EigVec[b+3] << "\n";
    b += 4;
  }
  for(int k=0; k<dn; k++) {
    eig << setw(22) << EigVec[b];
  }
  if(dn) eig << endl;
  xmlNodePtr det_data 
    = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eig.str().c_str());

  //xmlNodePtr det_data 
  //  = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)"\n");//eig.str().c_str());

  //int btot=SizeOfBasisSet*SizeOfBasisSet;
  //int n=btot/4, b=0;
  //int dn=btot-n*4;
  //char v[100];
  //for(int k=0; k<n; k++) {
  //  sprintf(v,"%20.12e %20.12e %20.12e %20.12e\n",EigVec[b],EigVec[b+1],EigVec[b+2],EigVec[b+3]);
  //  b += 4;
  //  xmlTextConcat(adet,(const xmlChar*)v,100);
  //}
  //for(int k=0; k<dn; k++) {
  //  sprintf(v,"%20.12e\n",EigVec[b++]);
  //  xmlTextConcat(adet,(const xmlChar*)v,25);
  //}
  //xmlNodePtr det_data 
  //  = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)EigVecU.c_str());
  xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());

  xmlNodePtr cur = xmlAddChild(slaterdet,adet);
  adet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(adet,(const xmlChar*)"id",(const xmlChar*)"downdet");
  xmlNewProp(adet,(const xmlChar*)"orbitals",(const xmlChar*)down_size.str().c_str());

  if(SpinRestricted)
    xmlNewProp(adet,(const xmlChar*)"ref",(const xmlChar*)"updet");
  else {
    //std::ostringstream eigD;
    //eigD.setf(std::ios::scientific, std::ios::floatfield);
    //eigD.precision(12);
    //int btot2=2*btot;
    //for(int k=btot; k<btot2;) {
    //  eig << setw(20) << EigVec[k++];
    //  if(k/4 == 0) eig << endl;
    //}
    //det_data = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)eigD.str().c_str());
    //xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());
  }

  cur = xmlAddSibling(cur,adet);
  return slaterdet;
}

xmlNodePtr QMCGaussianParserBase::createCenter(int iat, int off_) {

  CurrentCenter = IonName[GroupID[iat]];
  xmlNodePtr abasis = xmlNewNode(NULL,(const xmlChar*)"atomicBasisSet");
  xmlNewProp(abasis,(const xmlChar*)"name",(const xmlChar*)basisName.c_str());
  xmlNewProp(abasis,(const xmlChar*)"angular",(const xmlChar*)"spherical");
  xmlNewProp(abasis,(const xmlChar*)"type",(const xmlChar*)basisType.c_str());
  xmlNewProp(abasis,(const xmlChar*)"elementType",(const xmlChar*)CurrentCenter.c_str());
  xmlNewProp(abasis,(const xmlChar*)"normalized",(const xmlChar*)Normalized.c_str());
  //xmlNodePtr grid_ptr = xmlCopyNode(gridPtr,1);
  xmlAddChild(abasis,xmlCopyNode(gridPtr,1));
  //std::string hdf = CurrentCenter+"-"+basisName+".h5";
  //xmlNewProp(abasis,(const xmlChar*)"src", (const xmlChar*)(hdf.c_str()));
  for(int ig=gBound[iat], n=0; ig< gBound[iat+1]; ig++,n++) {
    createShell(n, ig, off_,abasis);
    off_ += gNumber[ig];
  }

  return abasis;
}

void
QMCGaussianParserBase::createShell(int n, int ig, int off_, xmlNodePtr abasis) {

  int gid(gShell[ig]);
  int ng(gNumber[ig]);

  xmlNodePtr ag = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
  xmlNodePtr ag1 = 0;

  char l_name[4],n_name[4],a_name[32];

  sprintf(a_name,"%s%d%d",CurrentCenter.c_str(),n,gShellID[gid]);
  cout << "*********** (species,input-gid,l) = " << CurrentCenter << " " << gid << " " << gShellID[gid] << " ^^" << a_name << "^^" << endl;
  sprintf(l_name,"%d",gShellID[gid]);
  sprintf(n_name,"%d",n);
  xmlNewProp(ag,(const xmlChar*)"rid", (const xmlChar*)a_name);
  xmlNewProp(ag,(const xmlChar*)"n", (const xmlChar*)n_name);
  xmlNewProp(ag,(const xmlChar*)"l", (const xmlChar*)l_name);
  xmlNewProp(ag,(const xmlChar*)"type", (const xmlChar*)"Gaussian");
  if(gid == 2) {
    sprintf(a_name,"%s%d1",CurrentCenter.c_str(),n);
    ag1 = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
    xmlNewProp(ag1,(const xmlChar*)"rid", (const xmlChar*)a_name);
    xmlNewProp(ag1,(const xmlChar*)"n", (const xmlChar*)n_name);
    xmlNewProp(ag1,(const xmlChar*)"l", (const xmlChar*)"1");
    xmlNewProp(ag1,(const xmlChar*)"type", (const xmlChar*)"Gaussian");
  }

  for(int ig=0, i=off_; ig<ng; ig++, i++) {
    std::ostringstream a,b,c;
    a.setf(std::ios::scientific, std::ios::floatfield);
    b.setf(std::ios::scientific, std::ios::floatfield);
    a.precision(12);
    b.precision(12);

    a<<gExp[i]; b<<gC0[i];

    xmlNodePtr anode = xmlNewNode(NULL,(const xmlChar*)"radfunc");
    xmlNewProp(anode,(const xmlChar*)"exponent", (const xmlChar*)a.str().c_str());
    xmlNewProp(anode,(const xmlChar*)"contraction", (const xmlChar*)b.str().c_str());
    xmlAddChild(ag,anode);
    if(gid ==2) {
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
  if(gid == 2) xmlAddChild(abasis,ag1);
}

void QMCGaussianParserBase::map2GridFunctors(xmlNodePtr cur) {

  using namespace ohmmsqmc;

  xmlNodePtr anchor = cur;
  //xmlNodePtr grid_ptr = 0;

  vector<xmlNodePtr> phi_ptr;
  vector<QuantumNumberType> nlms;

  int Lmax = 0;
  int current = 0;
  std::string acenter("none");

  const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"elementType");
  if(aptr) acenter = (const char*)aptr;

  xmlNodePtr grid_ptr=0;

  cur = anchor->children;
  while(cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname == "grid") 
      grid_ptr = cur;
    else if(cname == "basisGroup") {
      int n=1,l=0,m=0;
      const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"n");
      if(aptr) n = atoi((const char*)aptr);
      aptr = xmlGetProp(cur,(const xmlChar*)"l");
      if(aptr) l = atoi((const char*)aptr);

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

  if(grid_ptr == 0) {
    std::cout << "Grid is not defined: using default" << std::endl;
    //xmlAddChild(anchor,gridPtr);
    grid_ptr = xmlCopyNode(gridPtr,1);
    xmlAddChild(anchor,grid_ptr);
  }

  RGFBuilderBase::CenteredOrbitalType aos(Lmax);
  bool normalized(Normalized=="yes");
  RGFBuilderBase* rbuilder = new GTO2GridBuilder(normalized);
  rbuilder->setOrbitalSet(&aos,acenter);
  rbuilder->addGrid(grid_ptr);
  for(int i=0; i<nlms.size(); i++) {
    rbuilder->addRadialOrbital(phi_ptr[i],nlms[i]);
  }
  rbuilder->print(acenter,1);
}

void QMCGaussianParserBase::createGridNode(int argc, char** argv) {

  gridPtr = xmlNewNode(NULL,(const xmlChar*)"grid");
  string gridType("log");
  string gridFirst("1.e-6");
  string gridLast("1.e2");
  string gridSize("1001");
  int iargc=0;
  while(iargc<argc) {
    string a(argv[iargc]);
    if(a == "-gridtype") {
      gridType=argv[++iargc];
    } else if(a == "-frst") {
      gridFirst=argv[++iargc];
    } else if(a == "-last") {
      gridLast=argv[++iargc];
    } else if(a == "-size") {
      gridSize=argv[++iargc];
    }
    ++iargc;
  }
  xmlNewProp(gridPtr,(const xmlChar*)"type",(const xmlChar*)gridType.c_str());
  xmlNewProp(gridPtr,(const xmlChar*)"ri",(const xmlChar*)gridFirst.c_str());
  xmlNewProp(gridPtr,(const xmlChar*)"rf",(const xmlChar*)gridLast.c_str());
  xmlNewProp(gridPtr,(const xmlChar*)"npts",(const xmlChar*)gridSize.c_str());
}
