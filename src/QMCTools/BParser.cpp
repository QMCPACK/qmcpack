//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "QMCTools/BParser.h"
#include "QMCTools/BMakeFunc.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <set>
#include <map>

xmlNodePtr AGPLambda::createNode()
{
  xmlNodePtr aptr = xmlNewNode(NULL,(const xmlChar*)"lambda");
  std::ostringstream i,j,x;
  x.setf(std::ios::scientific, std::ios::floatfield);
  x.precision(8);
  i<<I;
  j<<J;
  x<<X;
  xmlNewProp(aptr,(const xmlChar*)"i",(const xmlChar*)i.str().c_str());
  xmlNewProp(aptr,(const xmlChar*)"j",(const xmlChar*)j.str().c_str());
  xmlNewProp(aptr,(const xmlChar*)"c",(const xmlChar*)x.str().c_str());
  return aptr;
}

BParser::BParser():
  DetShells(0),J3Shells(0),J2Index(0),
  DetSize(0), J3Size(0), DetNonZero(0), J3NonZero(0)
{
  basisName = "AGP";
  Normalized = "no";
  BMakeFuncBase::init();
}

BParser::BParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv),
  DetShells(0),J3Shells(0),J2Index(0),
  DetSize(0), J3Size(0), DetNonZero(0), J3NonZero(0)
{
  basisName = "AGP";
  Normalized = "no";
  BMakeFuncBase::init();
}

void BParser::parse(const std::string& fname)
{
  std::ifstream fin(fname.c_str());
  //# Nelup  #Nel  # Ion
  getwords(currentWords,fin);
  getwords(currentWords,fin);
  NumberOfAlpha =atoi(currentWords[0].c_str());
  NumberOfEls=atoi(currentWords[1].c_str());
  NumberOfAtoms = atoi(currentWords[2].c_str());
  NumberOfBeta=NumberOfEls-NumberOfAlpha;
  SpinRestricted = (NumberOfAlpha == NumberOfBeta);
  //# Shell Det.   # Shell Jas.
  getwords(currentWords,fin);
  getwords(currentWords,fin);
  DetShells=atoi(currentWords[0].c_str());
  J3Shells=atoi(currentWords[1].c_str());
  //# Jas 2body  # Det   #  3 body atomic par
  getwords(currentWords,fin);
  getwords(currentWords,fin);
  J2Index=atoi(currentWords[0].c_str());
  //# Det mat. =/0  # Jas mat. =/0
  getwords(currentWords,fin);
  getwords(currentWords,fin);
  DetNonZero=atoi(currentWords[0].c_str());
  J3NonZero=atoi(currentWords[1].c_str());
  // # Eq. Det atomic par.  # Eq. 3 body atomic. par.
  getwords(currentWords,fin);
  getwords(currentWords,fin);
  //# unconstrained iesfree,iessw,ieskinr,I/O flag
  getwords(currentWords,fin);
  getwords(currentWords,fin);
  //# Ion coordinates
  getGeometry(fin);
  if(J2Index != 0)
  {
    search(fin,"Jastrow two body");
    std::cout << "Found Jastrow Two Body" << std::endl;
  }
  search(fin,"atomic wf");
  getBasisSetForDet(fin);
  if(J3Shells != 0)
  {
    search(fin,"atomic Jastrow wf");
    std::cout << "Getting Atomic three-body Jastrow Wfs  " << std::endl;
    getBasisSetForJ3(fin);
  }
  search(fin,"Occupation atomic orbitals");
  getOccupationForDet(fin);
  if(J3Shells != 0)
  {
    search(fin,"Occupation atomic orbitals  Jastrow");
    getOccupationForJ3(fin);
  }
  search(fin,"detmat");
  getLambdaForDet(fin);
  if(J3Shells != 0)
  {
    search(fin,"jasmat");
    getLambdaForJ3(fin);
  }
}

void BParser::getGeometry(std::istream& is)
{
  getwords(currentWords,is);
  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  Qv.resize(NumberOfAtoms);
  for(int i=0; i< NumberOfAtoms; i++)
  {
    getwords(currentWords,is);
    double q =atof(currentWords[0].c_str());
    int atomic_number = static_cast<int>(q);
    int gid = IonSystem.GroupID[i]
              =IonSystem.getSpeciesSet().addSpecies(IonName[atomic_number]);
    IonSystem.getSpeciesSet()(IonChargeIndex,gid)=q;
    IonSystem.getSpeciesSet()(AtomicNumberIndex,gid)=q;
    GroupName[i]=IonName[atomic_number];
    int dir=0;
    for(int j=1; j<currentWords.size(); j++,dir++)
      IonSystem.R[i][dir]=atof(currentWords[j].c_str());
    if(dir < 3)
    {
      getwords(currentWords,is);
      int j=0;
      while(dir<3 && j<currentWords.size())
      {
        IonSystem.R[i][dir]=atof(currentWords[j].c_str());
        dir++;
        j++;
      }
    }
  }
}

void BParser::getBasisSetForDet(std::istream& is)
{
  int nitem =  DetShells;
  std::vector<int> rnl(IonSystem.getTotalNum(),0);
  detBasisPerAtom.resize(IonSystem.getTotalNum(),0);
  while(nitem)
  {
    getwords(currentWords,is);
    int angL=(atoi(currentWords[0].c_str())-1)/2;
    int nterms=atoi(currentWords[1].c_str());
    int iflag=atoi(currentWords[2].c_str());
    std::vector<std::string> items(nterms+1);
    //white-space at the end is killing me
    std::vector<std::string> temp;
    int inw=0;
    while(inw<nterms+1)
    {
      getwords(temp,is);
      items[inw++]=temp[0];
      for(int i=1; i<temp.size(); i++)
      {
        if(temp[i].size()>1)
          items[inw++]=temp[i];
      }
    }
    int centerID=atoi(items[0].c_str())-1;
    std::vector<BMakeFuncBase*>* b=0;
    std::map<int,std::vector<BMakeFuncBase*>*>::iterator it=detBasisSet.find(centerID);
    if(it == detBasisSet.end())
    {
      b = new std::vector<BMakeFuncBase*>;
      detBasisSet[centerID]=b;
    }
    else
    {
      b = (*it).second;
    }
    char rname[8];
    sprintf(rname,"R%d%d",centerID,rnl[centerID]);
    rnl[centerID]+=1;
    detBasisPerAtom[centerID]+=2*angL+1;
    //call the factory
    BMakeFuncBase* afunc = createBMakeFunc(iflag);
    afunc->BasisID=rname;
    afunc->L=angL;
    afunc->put(items);
    b->push_back(afunc);
    nitem--;
  }
}


void BParser::getBasisSetForJ3(std::istream& is)
{
  int nitem =  J3Shells;
  std::vector<int> rnl(IonSystem.getTotalNum(),0);
  j3BasisPerAtom.resize(IonSystem.getTotalNum(),0);
  while(nitem)
  {
    getwords(currentWords,is);
    int angL=(atoi(currentWords[0].c_str())-1)/2;
    int nterms=atoi(currentWords[1].c_str());
    int iflag=atoi(currentWords[2].c_str());
    std::vector<std::string> items(nterms+1);
    getwords(items,is);
    int centerID=atoi(items[0].c_str())-1;
    std::vector<BMakeFuncBase*>* b=0;
    std::map<int,std::vector<BMakeFuncBase*>*>::iterator it=j3BasisSet.find(centerID);
    if(it == j3BasisSet.end())
    {
      b = new std::vector<BMakeFuncBase*>;
      j3BasisSet[centerID]=b;
    }
    else
    {
      b = (*it).second;
    }
    char rname[8];
    sprintf(rname,"R%d%d",centerID,rnl[centerID]);
    rnl[centerID]+=1;
    j3BasisPerAtom[centerID]+=2*angL+1;
    //call the factory
    BMakeFuncBase* afunc = createBMakeFunc(iflag);
    afunc->BasisID=rname;
    afunc->L=angL;
    afunc->put(items);
    b->push_back(afunc);
    nitem--;
  }
}

void BParser::getOccupationForDet(std::istream& is)
{
  int tot= std::accumulate(detBasisPerAtom.begin(),detBasisPerAtom.end(),0);
  detOcc.resize(tot);
  DetSize=0;
  int item=0;
  while(item<tot)
  {
    getwords(currentWords,is);
    DetSize += detOcc[item]=atoi(currentWords[0].c_str());
    item++;
  }
  std::cout << "Occupation for determinants: size = " << DetSize << std::endl;
  copy(detOcc.begin(),detOcc.end(), std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
}

void BParser::getOccupationForJ3(std::istream& is)
{
  int tot= std::accumulate(j3BasisPerAtom.begin(),j3BasisPerAtom.end(),0);
  j3Occ.resize(tot);
  J3Size=0;
  int item=0;
  while(item<tot)
  {
    getwords(currentWords,is);
    J3Size += j3Occ[item]=atoi(currentWords[0].c_str());
    item++;
  }
  std::cout << "Occupation for J3: size= " << J3Size << std::endl;
  copy(j3Occ.begin(),j3Occ.end(), std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
}

void BParser::getLambdaForDet(std::istream& is)
{
  std::cout << "Number of non-zero determinant lambdas " << DetNonZero << std::endl;
  //detPairedLambda.reserve(DetNonZero);
  //detUnPairedLambda.reserve(DetNonZero);
  int d=0;
  while(d < DetNonZero)
  {
    getwords(currentWords,is);
    int i=atoi(currentWords[0].c_str());
    int j=atoi(currentWords[1].c_str());
    double x=atof(currentWords[2].c_str());
    if(j>DetSize)
    {
      detUnPairedLambda.push_back(AGPLambda(i,j-DetSize,x));
    }
    else
    {
      detPairedLambda.push_back(AGPLambda(i,j,x));
    }
    d++;
  }
  std::cout << "Non-zero elements of paried Lambda for the Determinant" << std::endl;
  for(int i=0; i<detPairedLambda.size(); i++)
  {
    std::cout << detPairedLambda[i].I << " " <<  detPairedLambda[i].J << " " <<  detPairedLambda[i].X << std::endl;
  }
  std::cout << "Non-zero elements of un-paried Lambda for the Determinant" << std::endl;
  for(int i=0; i<detUnPairedLambda.size(); i++)
  {
    std::cout << detUnPairedLambda[i].I << " " <<  detUnPairedLambda[i].J << " " <<  detUnPairedLambda[i].X << std::endl;
  }
}

void BParser::getLambdaForJ3(std::istream& is)
{
  j3Lambda.reserve(J3NonZero);
  int j3=0;
  while(j3<J3NonZero)
  {
    getwords(currentWords,is);
    j3Lambda.push_back(AGPLambda(currentWords));
    j3++;
  }
  std::cout << "Non-zero elements for J3 Lambda " << J3NonZero << std::endl;
  for(int i=0; i<j3Lambda.size(); i++)
  {
    std::cout << j3Lambda[i].I << " " <<  j3Lambda[i].J << " "
         <<  j3Lambda[i].X << std::endl;
  }
}

xmlNodePtr
BParser::createBasisSet(std::map<int,std::vector<BMakeFuncBase*>*>& bset,
                        std::vector<int>& basisPerAtom, std::vector<int>& occ, bool jastrow)
{
  xmlNodePtr bPtr = xmlNewNode(NULL, (const xmlChar*) "basisset");
  int boffset=0;
  std::vector<bool> newCenter(IonSystem.getSpeciesSet().getTotalNum(),true);
  std::map<int,std::vector<BMakeFuncBase*>*>::iterator it(bset.begin()), it_end(bset.end());
  while(it != it_end)
  {
    int id=(*it).first;
    int centerID=IonSystem.GroupID[id];
    if(newCenter[centerID])
    {
      std::ostringstream s;
      s << GroupName[id];
      xmlNodePtr cPtr = xmlNewNode(NULL, (const xmlChar*) "atomicBasisSet");
      xmlNewProp(cPtr,(const xmlChar*)"elementType",(const xmlChar*)s.str().c_str());
      if(jastrow)
      {
        xmlNewProp(cPtr,(const xmlChar*)"type",(const xmlChar*)"Gaussian");
        xmlNewProp(cPtr,(const xmlChar*)"normalized",(const xmlChar*)"yes");
        xmlNewProp(cPtr,(const xmlChar*)"angular",(const xmlChar*)"spherical");
        xmlNewProp(cPtr,(const xmlChar*)"expandYlm",(const xmlChar*)"no");
      }
      else
      {
        xmlNewProp(cPtr,(const xmlChar*)"type",(const xmlChar*)"STO");
        xmlNewProp(cPtr,(const xmlChar*)"normalized",(const xmlChar*)"no");
        xmlNewProp(cPtr,(const xmlChar*)"angular",(const xmlChar*)"spherical");
        xmlNewProp(cPtr,(const xmlChar*)"expandYlm",(const xmlChar*)"no");
      }
      newCenter[centerID]=false;
      std::vector<BMakeFuncBase*>& rgroup(*((*it).second));
      int b=boffset;
      for(int k=0; k<rgroup.size(); k++)
      {
        int angL=rgroup[k]->L;
        bool duplicated=false;
        switch(angL)
        {
        case(1):
          if(occ[b++])
          {
            rgroup[k]->M=1;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          if(occ[b++])
          {
            rgroup[k]->M=-1;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          if(occ[b++])
          {
            rgroup[k]->M=0;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          break;
        case(2):
          if(occ[b++])
          {
            rgroup[k]->M=0;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          if(occ[b++])
          {
            rgroup[k]->M=2;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          if(occ[b++])
          {
            rgroup[k]->M=-2;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          if(occ[b++])
          {
            rgroup[k]->M=-1;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          if(occ[b++])
          {
            rgroup[k]->M=1;
            xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
            duplicated=true;
          }
          break;
        default:
          for(int m=-angL; m<=angL; m++)
          {
            if(occ[b++])
            {
              rgroup[k]->M=m;
              xmlAddChild(cPtr,rgroup[k]->createBasisGroup(duplicated));
              duplicated=true;
            }
          }
        }//switch(angL)
      }//angL
      xmlAddChild(bPtr,cPtr);
    }
    boffset+=basisPerAtom[id];
    ++it;
  }
  return bPtr;
}

/** xmlNode for determinant set
 *
 *\xmlonly
   <determinantset type="AGP" transform="yes" source="ion0">
     <coefficients offset="1" size="28">
       <lamdba i="1" j="1" c="1.0">

     </coefficients>
   </determinantset>
 \endxmlonly
 */
xmlNodePtr BParser::createDeterminantSet()
{
  xmlNodePtr detPtr = xmlNewNode(NULL, (const xmlChar*) "determinantset");
  xmlNewProp(detPtr,(const xmlChar*)"type",(const xmlChar*)"AGP");
  xmlNewProp(detPtr,(const xmlChar*)"transform",(const xmlChar*)"yes");
  xmlNewProp(detPtr,(const xmlChar*)"source",(const xmlChar*)IonSystem.getName().c_str());
  std::cout << "Checking the basis set for the determinants " << std::endl;
  xmlAddChild(detPtr,createBasisSet(detBasisSet,detBasisPerAtom, detOcc,false));
  //add basis here
  if(detPairedLambda.size())
  {
    std::ostringstream s;
    s << DetSize;
    xmlNodePtr cptr = xmlNewNode(NULL, BAD_CAST "coefficients");
    xmlNewProp(cptr,(const xmlChar*)"offset",(const xmlChar*)"1");
    xmlNewProp(cptr,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
    for(int i=0; i<detPairedLambda.size(); i++)
    {
      xmlAddChild(cptr,detPairedLambda[i].createNode());
    }
    xmlAddChild(detPtr,cptr);
  }
  if(detUnPairedLambda.size())
  {
    std::ostringstream s;
    s << NumberOfAlpha-NumberOfBeta;
    xmlNodePtr cptr = xmlNewNode(NULL, BAD_CAST "unpaired");
    xmlNewProp(cptr,(const xmlChar*)"offset",(const xmlChar*)"1");
    xmlNewProp(cptr,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
    for(int i=0; i<detUnPairedLambda.size(); i++)
    {
      xmlAddChild(cptr,detUnPairedLambda[i].createNode());
    }
    xmlAddChild(detPtr,cptr);
  }
  return detPtr;
}

/** xmlNode for J3 set
 *
 *\xmlonly
   <jastrow name="J3G" type="Three-Body-Geminal" function="gto" source="ion0">
     <coefficients offset="1" size="28">
       <lamdba i="1" j="1" c="1.0">
     </coefficients>
   </jastrow>
 \endxmlonly
 */
xmlNodePtr BParser::createJ3()
{
  xmlNodePtr j3Ptr = xmlNewNode(NULL, (const xmlChar*) "jastrow");
  xmlNewProp(j3Ptr,(const xmlChar*)"name",(const xmlChar*)"J3G");
  xmlNewProp(j3Ptr,(const xmlChar*)"type",(const xmlChar*)"Three-Body-Geminal");
  xmlNewProp(j3Ptr,(const xmlChar*)"function",(const xmlChar*)"gto");
  xmlNewProp(j3Ptr,(const xmlChar*)"transform",(const xmlChar*)"yes");
  xmlNewProp(j3Ptr,(const xmlChar*)"source",(const xmlChar*)IonSystem.getName().c_str());
  std::cout << "Checking the basis set for the Jastrow " << std::endl;
  xmlAddChild(j3Ptr, createBasisSet(j3BasisSet,j3BasisPerAtom, j3Occ, true));
  std::ostringstream s;
  s << J3Size;
  xmlNodePtr cptr = xmlNewNode(NULL, BAD_CAST "coefficients");
  xmlNewProp(cptr,(const xmlChar*)"offset",(const xmlChar*)"1");
  xmlNewProp(cptr,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
  for(int i=0; i<j3Lambda.size(); i++)
  {
    xmlAddChild(cptr,j3Lambda[i].createNode());
  }
  xmlAddChild(j3Ptr,cptr);
  return j3Ptr;
}

void BParser::dump(const std::string& psi_tag,
                   const std::string& ion_tag)
{
  std::cout << " BParser::dump " << std::endl;
  xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  {
    //particleset
    xmlAddChild(qm_root,createElectronSet(ion_tag));
    xmlAddChild(qm_root,createIonSet());
    //wavefunction
    xmlNodePtr wfPtr = xmlNewNode(NULL,(const xmlChar*)"wavefunction");
    xmlNewProp(wfPtr,(const xmlChar*)"id",(const xmlChar*)psi_tag.c_str());
    xmlNewProp(wfPtr,(const xmlChar*)"target",(const xmlChar*)"e");
    if(DetSize)
    {
      xmlAddChild(wfPtr,createDeterminantSet());
    }
    if(J3Size)
    {
      xmlAddChild(wfPtr,createJ3());
    }
    xmlAddChild(qm_root,wfPtr);
  }
  xmlDocSetRootElement(doc, qm_root);
  std::string fname = basisName+".xml";
  xmlSaveFormatFile(fname.c_str(),doc,1);
  xmlFreeDoc(doc);
}
