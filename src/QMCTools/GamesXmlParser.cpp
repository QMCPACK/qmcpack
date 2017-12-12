//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "QMCTools/GamesXmlParser.h"
#include "QMCHamiltonians/ConservedEnergy.h"

//namespace qmcplusplus {

GamesXmlParser::GamesXmlParser()
{
  basisName = "GamesXml";
  Normalized = "no";
}

GamesXmlParser::GamesXmlParser(int argc, char** argv):
  QMCGaussianParserBase(argc,argv)
{
  basisName = "GamesXml";
  Normalized = "no";
}


void GamesXmlParser::parse(const std::string& fname)
{
  if(multideterminant)
  {
    std::cerr <<"Multideterminant parser for Gamesxml is not implemented. \n";
    exit(201);
  }
  xmlDocPtr m_doc = xmlParseFile(fname.c_str());
  if (m_doc == NULL)
  {
    ERRORMSG("File " << fname << " is invalid")
    xmlFreeDoc(m_doc);
    return;
  }
  xmlNodePtr cur = xmlDocGetRootElement(m_doc);
  if(!xmlStrEqual(cur->name,(const xmlChar*)"GAMESS"))
  {
    ERRORMSG("File " << fname << " does not have GAMESS as its root. Invalid")
    xmlFreeDoc(m_doc);
    return;
  }
  //xmlNodePtr for atoms
  std::vector<xmlNodePtr> aPtrList;
  //xmlNodePtr for eigvectors
  std::vector<xmlNodePtr> ePtrList;
  //xmlNodePtr for gaussian basis
  std::vector<xmlNodePtr> bPtrList;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "IN")
    {
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)cur1->name);
        if(cname1 == "RUN_TITLE")
        {
          std::string atitle;
          putContent(atitle,cur1);
          std::string::size_type wh=atitle.find("...");
          if(wh<atitle.size())
            atitle.erase(wh,atitle.size()-wh);
          Title = atitle;
        }
        else
          if(cname1 == "CONTRL")
          {
            getControlParameters(cur1);
          }
        cur1=cur1->next;
      }//everything within IN
    }
    else
      if(cname == "OUT")
      {
        xmlNodePtr cur1=cur->children;
        while(cur1 != NULL)
        {
          std::string cname1((const char*)cur1->name);
          if(cname1 == "SYSTEM_STATE")
          {
            //Unit needs to be generalized!!
            std::string unitL((const char*)xmlGetProp(cur1,(const xmlChar*)"UNITS"));
            if(unitL == "ANGS")
              BohrUnit=false;
            xmlNodePtr cur2 = cur1->children;
            while(cur2 != NULL)
            {
              std::string cname2((const char*)cur2->name);
              if(cname2 == "ATOM")
              {
                aPtrList.push_back(cur2);
              }
              else
                if(cname2 == "VEC")
                {
                  ePtrList.push_back(cur2);
                }
              cur2=cur2->next;
            }
          }
          else
            if(cname1 == "PDATA")
            {
              xmlNodePtr cur2 = cur1->children;
              while(cur2 != NULL)
              {
                std::string cname2((const char*)cur2->name);
                if(cname2 == "PATOMIC_BASIS_SET")
                {
                  bPtrList.push_back(cur2);
                }
                cur2=cur2->next;
              }
            }
          cur1=cur1->next;
        }//everything within OUT
      }
    cur=cur->next;
  }
  //xmlXPathContextPtr m_context = xmlXPathNewContext(m_doc);
  getGeometry(aPtrList);
  getGaussianCenters(bPtrList);
  getEigVectors(ePtrList);
  //xmlXPathFreeContext(m_context);
  xmlFreeDoc(m_doc);
}

void GamesXmlParser::getControlParameters(xmlNodePtr cur)
{
  std::string a;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "SCFTYP")
    {
      putContent(a,cur);
      if(a == "RHF" || a == "ROHF")
        SpinRestricted=true;
      else
        if(a == "URHF" || a == "UHF")
          SpinRestricted=false;
    }
    else
      if(cname == "MULT")
      {
        putContent(SpinMultiplicity,cur);
      }
    cur=cur->next;
  }
}

void GamesXmlParser::getGeometry(std::vector<xmlNodePtr>& aPtrList)
{
  NumberOfAtoms = aPtrList.size();
  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  Qv.resize(NumberOfAtoms);
  double nel=0;
  for(int i=0; i<NumberOfAtoms; i++)
  {
    xmlNodePtr cur=aPtrList[i]->children;
    std::string atName;
    double q;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == "ATOM_NAME")
      {
        //string aname;
        putContent(atName,cur);
      }
      else
        if(cname == "ATOMIC_NUMBER")
        {
          putContent(q,cur);
          Qv[i]=q;
          nel+=q;
        }
        else
          if(cname == "ATOM_POSITION")
          {
            xmlNodePtr tcur=cur->children;
            while(tcur!=NULL)
            {
              std::string tname((const char*)tcur->name);
              double x,y,z;
              if(tname == "XCOORD")
                putContent(x,tcur);
              else
                if(tname == "YCOORD")
                  putContent(y,tcur);
                else
                  if(tname == "ZCOORD")
                    putContent(z,tcur);
              IonSystem.R[i]=SingleParticlePos_t(x,y,z);
              tcur=tcur->next;
            }
          }
      cur=cur->next;
    }//loop-cur
    int atomic_number = static_cast<int>(q);
    int gid = IonSystem.GroupID[i]
              =IonSystem.getSpeciesSet().addSpecies(IonName[atomic_number]);
    IonSystem.getSpeciesSet()(IonChargeIndex,gid)=q;
    IonSystem.getSpeciesSet()(AtomicNumberIndex,gid)=q;
    GroupName[i]=IonName[atomic_number];
  }//i
  NumberOfEls=static_cast<int>(nel);
  NumberOfBeta=(NumberOfEls-NumberOfAlpha)/2;
  NumberOfAlpha=NumberOfEls-NumberOfBeta;
  std::cout << "Number of atoms " << NumberOfAtoms << std::endl;
  std::cout << "Number of electrons " << NumberOfEls << std::endl;
  std::cout << "Number of electrons (ALPHA) " << NumberOfAlpha << std::endl;
  std::cout << "Number of electrons (BETA) " << NumberOfBeta << std::endl;
  std::cout << "Group ID " << std::endl;
  copy(IonSystem.GroupID.begin(), IonSystem.GroupID.end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
}

void GamesXmlParser::getGaussianCenters(std::vector<xmlNodePtr>& bPtrList)
{
  //if(bPtrList.size() != aPtrList.size())
  gBound.push_back(0);
  int offset=0;
  double zeta,c;
  SizeOfBasisSet=0;
  for(int i=0; i<bPtrList.size(); i++)
  {
    std::string p;
    int ng_tot=0,ng;
    xmlNodePtr cur=bPtrList[i]->children;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == "PSHELL")
      {
        ng_tot++;
        xmlNodePtr cur1=cur->children;
        int gshellType=1;
        while(cur1!= NULL)
        {
          std::string tname((const char*)cur1->name);
          if(tname == "PTYPE")
          {
            putContent(p,cur1);
            if(p == "S")
            {
              gshellType=1;
              SizeOfBasisSet+=1;
            }
            else
              if(p == "P")
              {
                gshellType=3;
                SizeOfBasisSet+=3;
              }
              else
                if(p == "D")
                {
                  gshellType=4;
                  SizeOfBasisSet+=5;
                }
            gShell.push_back(gshellType);
          }
          else
            if(tname == "PNGAUSS")
            {
              putContent(ng,cur1);
              gNumber.push_back(ng);
              //  ng_tot+=ng;
            }
            else
              if(tname == "PGAUSSIAN")
              {
                xmlNodePtr cur2=cur1->children;
                while(cur2 != NULL)
                {
                  std::string cname2((const char*)cur2->name);
                  if(cname2 == "PZETA")
                  {
                    putContent(zeta,cur2);
                    gExp.push_back(zeta);
                  }
                  else
                    if(cname2 == "PCONE")
                    {
                      putContent(c,cur2);
                      gC0.push_back(c);
                    }
                  cur2=cur2->next;
                }
                std::cout << "zeta,c " << zeta << " " << c << std::endl;
              }
          cur1=cur1->next;
        }
      }
      cur=cur->next;
    }
    offset+=ng_tot;
    gBound.push_back(offset);
  }
  std::cout << "Bound of gauassians " << std::endl;
  copy(gBound.begin(), gBound.end(),std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
  std::cout << "Number of shell type " << std::endl;
  copy(gShell.begin(), gShell.end(),std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
  std::cout << "Number of gaussians per shell " << std::endl;
  copy(gNumber.begin(), gNumber.end(),std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
  gC1.resize(gC0.size(),0.0);
}

void GamesXmlParser::getEigVectors(std::vector<xmlNodePtr>& ePtrList)
{
  std::vector<xmlNodePtr> a;
  //vector<int> numorb(ePtrList.size());
  for(int i=0; i<ePtrList.size(); i++)
  {
    xmlNodePtr cur=ePtrList[i]->children;
    int n=0;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == "ORB")
      {
        a.push_back(cur);
        n++;
      }
      cur=cur->next;
    }
  }
  //adhoc
  //if(ePtrList.size()>1) SpinRestricted=false;
  std::cout << "Size of eig vectors " << a.size() << " x " << SizeOfBasisSet << std::endl;
  EigVal_alpha.resize(SizeOfBasisSet);
  EigVal_beta.resize(SizeOfBasisSet);
// mmorales: not sure if this is correct, this leaves it unmodified
  numMO = SizeOfBasisSet;
  EigVec.resize(a.size()*SizeOfBasisSet);
  int ii=0;
  double x;
  for(int i=0; i<a.size(); i++)
  {
    xmlNodePtr cur=a[i]->children;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname== "EIGENVALUE")
      {
        if(i<SizeOfBasisSet)
        {
          putContent(x,cur);
          EigVal_alpha[i]=x;
        }
        else
        {
          putContent(x,cur);
          EigVal_beta[i-SizeOfBasisSet]=x;
        }
      }
      else
        if(cname == "BASIS_COEFF")
        {
          putContent(x,cur);
          EigVec[ii]=x;
          ii++;
        }
      cur=cur->next;
    }
  }
  if(SpinRestricted)
    EigVal_beta=EigVal_alpha;
}
//}
