#include "QMCTools/QMCGaussianParserBase.h"
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
using namespace std;
#include "QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.h"

std::vector<std::string> QMCGaussianParserBase::IonName;
std::vector<std::string> QMCGaussianParserBase::gShellType;
std::vector<int> QMCGaussianParserBase::gShellID;

void QMCGaussianParserBase::init() {
  IonName.resize(24);
  IonName[1] = "H"; IonName[2] = "He"; IonName[3] = "Li";
  IonName[4] = "Be"; IonName[5] = "B"; IonName[6] = "C";
  IonName[7] = "N"; IonName[8] = "O"; IonName[9] = "F";
  IonName[10] = "Ne"; IonName[11] = "Na"; IonName[12] = "Mg";
  IonName[13] = "Al"; IonName[14] = "Si"; IonName[15] = "P";
  IonName[16] = "S"; IonName[17] = "Cl"; IonName[18] = "Ar";
  gShellType.resize(5);
  gShellType[1]="s"; gShellType[2]="sp"; gShellType[3]="p"; gShellType[4]="d";
  gShellID.resize(5);
  gShellID[1]=0; gShellID[2]=0; gShellID[3]=1; gShellID[4]=2;
}


xmlNodePtr QMCGaussianParserBase::createBasisSet() {

  //xmlNodePtr wf_root = xmlNewNode(NULL, BAD_CAST "determinantset");
  //xmlNewProp(wf_root,(const xmlChar*)"type",(const xmlChar*)"MolecularOrbital");

  xmlNodePtr bset = xmlNewNode(NULL,(const xmlChar*)"basisset");
  xmlNewProp(bset,(const xmlChar*)"ref",(const xmlChar*)"i");
  xmlNodePtr cur=0;

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
      xmlNodePtr abasis = createCenter(iat,gtot);
      if(cur)
        cur = xmlAddSibling(cur,abasis);
      else
        cur = xmlAddChild(bset,abasis);
    } else {
      ng = (*it).second;
    }
    gtot += ng;
  }

  return bset;
  //xmlAddChild(wf_root,bset);
  //return wf_root;
}

xmlNodePtr QMCGaussianParserBase::createDeterminantSet() {
  xmlNodePtr slaterdet = xmlNewNode(NULL,(const xmlChar*)"slaterdeterminant");

  //check spin-dependent properties
  int nup = NumberOfEls/2;
  int ndown = NumberOfEls-nup;
  std::ostringstream up_size, down_size, b_size;
  up_size <<nup; down_size << ndown; b_size<<SizeOfBasisSet;

  //create a determinant Up
  xmlNodePtr adet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(adet,(const xmlChar*)"id",(const xmlChar*)"updet");
  xmlNewProp(adet,(const xmlChar*)"orbitals",(const xmlChar*)up_size.str().c_str());

  //std::ostreagstream up_occ, down_occ;
  //for(int i=0; i<nup; i++) up_occ << "1 ";
  //for(int i=0; i<down; i++) up_occ << "1 ";

  xmlNodePtr det_data 
    = xmlNewTextChild(adet,NULL,(const xmlChar*)"coefficient",(const xmlChar*)EigVecU.c_str());
  xmlNewProp(det_data,(const xmlChar*)"size",(const xmlChar*)b_size.str().c_str());

  xmlNodePtr cur = xmlAddChild(slaterdet,adet);
  adet = xmlNewNode(NULL,(const xmlChar*)"determinant");
  xmlNewProp(adet,(const xmlChar*)"id",(const xmlChar*)"downdet");
  xmlNewProp(adet,(const xmlChar*)"orbitals",(const xmlChar*)down_size.str().c_str());
  xmlNewProp(adet,(const xmlChar*)"ref",(const xmlChar*)"updet");
  cur = xmlAddSibling(cur,adet);

  return slaterdet;
}
xmlNodePtr QMCGaussianParserBase::createCenter(int iat, int _off) {

  std::string aname(IonName[GroupID[iat]]);
  xmlNodePtr abasis = xmlNewNode(NULL,(const xmlChar*)"basis");
  xmlNewProp(abasis,(const xmlChar*)"type",(const xmlChar*)"GTO");
  xmlNewProp(abasis,(const xmlChar*)"species",(const xmlChar*)aname.c_str());

  xmlNodePtr agrid = xmlNewNode(NULL,(const xmlChar*)"grid");
  xmlNewProp(agrid,(const xmlChar*)"ri",    (const xmlChar*)"1.e-6");
  xmlNewProp(agrid,(const xmlChar*)"rf",    (const xmlChar*)"10.0");
  xmlNewProp(agrid,(const xmlChar*)"npts",    (const xmlChar*)"1001");
  xmlNodePtr cur = xmlAddChild(abasis,agrid);

  for(int ig=gBound[iat], n=0; ig<gBound[iat+1]; ig++,n++) {

    int gid = gShell[ig];
    xmlNodePtr ag=createShell(_off,gNumber[ig],gC0);

    //Add attributes for (n,l)
    std::ostringstream id_assigned,n_name,l_name;
    n_name << n; l_name << gShellID[gid];
    id_assigned << aname << n << "_"<<gid;
    xmlNewProp(ag,(const xmlChar*)"id",    (const xmlChar*)(id_assigned.str().c_str()));
    xmlNewProp(ag,(const xmlChar*)"symbol",(const xmlChar*)(gShellType[gid].c_str()));
    xmlNewProp(ag,(const xmlChar*)"n",     (const xmlChar*)(n_name.str().c_str()));
    xmlNewProp(ag,(const xmlChar*)"l",     (const xmlChar*)(l_name.str().c_str()));
    cur = xmlAddSibling(cur,ag);

    if(gid == 2) {//special case for 2 (sp) need to add P
      ag=createShell(_off,gNumber[ig],gC1);
      id_assigned << "1";
      xmlChar pname[]={"1"};
      xmlNewProp(ag,(const xmlChar*)"id",    (const xmlChar*)(id_assigned.str().c_str()));
      xmlNewProp(ag,(const xmlChar*)"symbol",(const xmlChar*)(gShellType[gid].c_str()));
      xmlNewProp(ag,(const xmlChar*)"n",     (const xmlChar*)(n_name.str().c_str()));
      xmlNewProp(ag,(const xmlChar*)"l",     pname);
      cur = xmlAddSibling(cur,ag);
    }
    _off += gNumber[ig];
  }

  //std::cout << "Create the orbitals " << std::endl;
  //map2GridFunctors(abasis);
  return abasis;
}

xmlNodePtr 
QMCGaussianParserBase::createShell(int _off, int ng, const std::vector<double>& coeff) {

  xmlNodePtr ag = xmlNewNode(NULL,(const xmlChar*)"phi");
  for(int ig=0; ig<ng; ig++, _off++) {
    std::ostringstream a,b;
    a.setf(std::ios::scientific, std::ios::floatfield);
    a.precision(12);
    b.setf(std::ios::scientific, std::ios::floatfield);
    b.precision(12);
    xmlNodePtr aR = xmlNewNode(NULL,(const xmlChar*)"Rnl");
    a << gExp[_off]; b<< coeff[_off];
    xmlNewProp(aR,(const xmlChar*)"alpha",(const xmlChar*)a.str().c_str());
    xmlNewProp(aR,(const xmlChar*)"c",(const xmlChar*)b.str().c_str());
    xmlAddChild(ag,aR);
  }

  return ag;
}

void QMCGaussianParserBase::map2GridFunctors(xmlNodePtr cur) {

  using namespace ohmmsqmc;

  xmlNodePtr anchor = cur;
  xmlNodePtr grid_ptr = 0;

  vector<xmlNodePtr> phi_ptr;
  vector<QuantumNumberType> nlms;

  int Lmax = 0;
  int current = 0;
  std::string acenter("none");

  const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"species");
  if(aptr) acenter = (const char*)aptr;

  cout << "Building basis set for " << acenter << endl;

  cur = anchor->children;
  while(cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname == "grid") 
      grid_ptr = cur;
    else if(cname == "phi") {
      phi_ptr.push_back(cur);
      nlms.push_back(QuantumNumberType());
      int n=1,l=0,m=0;
      const xmlChar* lptr = xmlGetProp(cur,(const xmlChar*)"l");
      if(lptr) l = atoi((const char*)lptr);
      const xmlChar* nptr = xmlGetProp(cur,(const xmlChar*)"n");
      if(nptr) n = atoi((const char*)nptr);
      const xmlChar* mptr = xmlGetProp(cur,(const xmlChar*)"m");
      if(mptr) m = atoi((const char*)mptr);
      Lmax = std::max(l,Lmax);
      nlms[current][0]=n;
      nlms[current][1]=l;
      nlms[current][2]=m;
      current++;
    }
    cur = cur->next;
  }

  if(grid_ptr == 0) {
    std::cout << "Grid is not defined" << std::endl;
    return;
  }

  RGFBuilderBase::CenteredOrbitalType aos(Lmax);

  RGFBuilderBase* rbuilder = new GTO2GridBuilder(true);
  rbuilder->setOrbitalSet(&aos,acenter);
  rbuilder->addGrid(grid_ptr);
  for(int i=0; i<nlms.size(); i++) {
    rbuilder->addRadialOrbital(phi_ptr[i],nlms[i]);
  }
  rbuilder->print(acenter,1);

}

