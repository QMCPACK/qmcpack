//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Quadrature.h"
#include "Numerics/Transform2GridFunctor.h"
#include "QMCHamiltonians/FSAtomPseudoPot.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/SimpleParser.h"
#include "Message/CommOperators.h"
#include <cmath>
#include <qmc_common.h>


namespace qmcplusplus
{

ECPComponentBuilder::ECPComponentBuilder(const std::string& aname, Communicate* c):
  MPIObjectBase(c),
  RcutMax(-1), NumNonLocal(0), Lmax(0), AtomicNumber(0), Zeff(0), Species(aname), Nrule(4),
  grid_global(0),pp_loc(0), pp_nonloc(0)
{
  angMon["s"]=0;
  angMon["p"]=1;
  angMon["d"]=2;
  angMon["f"]=3;
  angMon["g"]=4;
  angMon["0"]=0;
  angMon["1"]=1;
  angMon["2"]=2;
  angMon["3"]=3;
  angMon["4"]=4;
}

bool ECPComponentBuilder::parse(const std::string& fname, xmlNodePtr cur)
{
  const xmlChar* rptr=xmlGetProp(cur,(const xmlChar*)"cutoff");
  if(rptr != NULL)
    RcutMax = atof((const char*)rptr);

  return read_pp_file(fname);
}

int ReadFileBuffer::get_file_length(std::ifstream *f) const
{
    f->seekg (0, std::ios::end);
    int len = f->tellg();
    f->seekg (0, std::ios::beg);
    return len;
}

void ReadFileBuffer::reset()
{
  if (is_open)
  {
    if(myComm == NULL || myComm->rank() == 0)
    {
      delete fin;
      fin = NULL;
    }
    delete[] cbuffer;
    cbuffer = NULL;
    is_open = false;
    length = 0;
  }
}

bool ReadFileBuffer::open_file(const std::string &fname)
{

  reset();

  if (myComm == NULL || myComm->rank() == 0)
  {
    fin = new std::ifstream(fname.c_str());
    if (fin->is_open())
      is_open=true;
  }
  if (myComm) myComm->bcast(is_open);
  return is_open;
}

bool ReadFileBuffer::read_contents()
{
  if (!is_open)
    return false;

  if(myComm == NULL || myComm->rank() == 0)
  {
    length = get_file_length(fin);
  }
  if (myComm) myComm->bcast(length);

  cbuffer = new char[length+1];
  cbuffer[length] = '\0';

  if (myComm == NULL || myComm->rank() == 0)
    fin->read (cbuffer,length);

  if (myComm != NULL)
    myComm->bcast(cbuffer,length);

  return true;
}


bool ECPComponentBuilder::read_pp_file(const std::string &fname)
{
  ReadFileBuffer buf(myComm);
  bool okay = buf.open_file(fname);
  if(!okay)
  {
    APP_ABORT("ECPComponentBuilder::read_pp_file  Missing PP file " + fname +"\n");
  }

  okay = buf.read_contents();
  if(!okay)
  {
    APP_ABORT("ECPComponentBuilder::read_pp_file Unable to read PP file " + fname +"\n");
  }

  xmlDocPtr m_doc = xmlReadMemory(buf.contents(),buf.length,NULL,NULL,0);

  if (m_doc == NULL)
  {
    xmlFreeDoc(m_doc);
    APP_ABORT("ECPComponentBuilder::read_pp_file xml file "+fname+" is invalid");
  }
  // Check the document is of the right kind
  xmlNodePtr cur = xmlDocGetRootElement(m_doc);
  if (cur == NULL)
  {
    xmlFreeDoc(m_doc);
    APP_ABORT("Empty document");
  }
  bool success=put(cur);
  xmlFreeDoc(m_doc);
  return success;
}

bool ECPComponentBuilder::put(xmlNodePtr cur)
{
  int nk=0;
  //vector<RealType> kpts;
  std::vector<xmlNodePtr> semiPtr;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "header")
    {
      Zeff = atoi((const char*)xmlGetProp(cur,(const xmlChar*)"zval"));
      AtomicNumber = atoi((const char*)xmlGetProp(cur,(const xmlChar*)"atomic-number"));
    }
    else if(cname == "grid")
    {
      //capture the global grid
      grid_global = createGrid(cur);
    }
    else if(cname == "semilocal")
    {
      semiPtr.push_back(cur);//save the pointer
    }
    else if(cname == "local")
    {
      buildLocal(cur);
    }
    // else if(cname == "sphericalGrid")
    // {
    //  nk=atoi((const char*)xmlGetProp(cur,(const xmlChar*)"size"));
    //  kpts.resize(nk*4);
    //  putContent(kpts,cur);
    // }
    cur=cur->next;
  }
  if(semiPtr.size())
  {
    if(pp_nonloc==0)
      pp_nonloc=new NonLocalECPComponent;
    if(pp_loc)
    {
      for(int i=0; i<semiPtr.size(); i++)
        addSemiLocal(semiPtr[i]);
    }
    else
      buildSemiLocalAndLocal(semiPtr);
  }
  if(pp_nonloc)
  {
    SetQuadratureRule(Nrule);
    app_log() << "    Non-local pseudopotential parameters" << std::endl;
    pp_nonloc->print(app_log());
    app_log() << "    Maximum cutoff radius " << pp_nonloc->Rmax << std::endl;
  }
  return true;
}

void ECPComponentBuilder::printECPTable()
{
  if(!qmc_common.io_node || qmc_common.mpi_groups>1) return;

  char fname[12];
  sprintf(fname,"%s.pp.dat",Species.c_str());
  std::ofstream fout(fname);
  fout.setf(std::ios::scientific, std::ios::floatfield);
  fout.precision(12);
  int nl=pp_nonloc?pp_nonloc->nlpp_m.size():0;
  RealType d=1.7e-2;
  RealType rt=0.13*d;
  if(nl)
  {
    fout << "#  Lmax = " << Lmax+1 << " nonlocal L channels" << nl << std::endl;
    fout << "#  Units = bohr hartree " << std::endl;
    fout << "#  r  -r*V/zeff   Vnl ... " << std::endl;
    while(rt<5)
    {
      fout << rt << std::setw(25) << pp_loc->splint(rt);
      for(int l=0; l<nl; l++)
        fout << std::setw(25) << pp_nonloc->nlpp_m[l]->splint(rt);
      fout << std::endl;
      rt+=d;
    }
  }
  else
  {
    fout << "#  Units = bohr hartree " << std::endl;
    fout << "#  r  -r*V/zeff " << std::endl;
    while(rt<5)
    {
      fout << rt << std::setw(25) << pp_loc->splint(rt) << std::endl;;
      rt+=d;
    }
  }
}

void ECPComponentBuilder::SetQuadratureRule(int rule)
{
  Quadrature3D<RealType> myRule(rule);
  pp_nonloc->sgridxyz_m=myRule.xyz_m;
  pp_nonloc->sgridweight_m=myRule.weight_m;
  // Allocate storage for wave function ratios
  pp_nonloc->resize_warrays(myRule.nk,NumNonLocal,Lmax);
}

} // namespace qmcPlusPlus
