//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Jastrow/JastrowBasisBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/AtomicBasisBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/BsplineAOBuilder.h"
#include "QMCWaveFunctions/Jastrow/CBSOBuilder.h"
#include "QMCWaveFunctions/SparseLocalizedBasisSet.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

/** constructor
 * \param els reference to the electrons
 * \param ions reference to the ions
 */
JastrowBasisBuilder::JastrowBasisBuilder(ParticleSet& els, ParticleSet& ions,
    const std::string& functype, bool usespline):
  targetPtcl(els), sourcePtcl(ions), UseSpline(usespline),FuncType(functype), myBasisSet(0)
{
}

JastrowBasisBuilder::~JastrowBasisBuilder()
{
}

template<typename RFBUILDER>
void JastrowBasisBuilder::createLocalizedBasisSet(xmlNodePtr cur)
{
  typedef typename RFBUILDER::CenteredOrbitalType COT;
  typedef SparseLocalizedBasisSet<COT> ThisBasisSetType;
  ThisBasisSetType* curBasis= new ThisBasisSetType(sourcePtcl,targetPtcl);
  //create the size vector with zeros
  SizeOfBasisPerCenter.resize(sourcePtcl.getSpeciesSet().getTotalNum(),0);
  //create the basis set
  cur = cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "atomicBasisSet")
    {
      const xmlChar* eptr=xmlGetProp(cur,(const xmlChar*)"elementType");
      if(eptr == NULL)
      {
        app_error() << "   Missing elementType attribute of atomicBasisSet.\n Abort at MOBasisBuilder::put " << std::endl;
        OHMMS::Controller->abort();
      }
      std::string elementType((const char*)eptr);
      std::map<std::string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
      if(it == aoBuilders.end())
      {
        AtomicBasisBuilder<RFBUILDER>* any = new AtomicBasisBuilder<RFBUILDER>(elementType);
        any->put(cur);
        COT* aoBasis= any->createAOSet(cur);
        if(aoBasis)
          //add the new atomic basis to the basis set
        {
          int activeCenter =sourcePtcl.getSpeciesSet().findSpecies(elementType);
          curBasis->add(activeCenter, aoBasis);
          aoBuilders[elementType]=any;
          printRadialFunctors(elementType,aoBasis);
          SizeOfBasisPerCenter[activeCenter]=aoBasis->getBasisSetSize();
        }
      }
    }
    cur = cur->next;
  }
  //resize the basis set
  curBasis->setBasisSetSize(-1);
  //clean up
  if(myBasisSet) delete myBasisSet;
  myBasisSet=curBasis;
}

bool JastrowBasisBuilder::put(xmlNodePtr cur)
{
  if(myBasisSet)
    return true;
  if(UseSpline)
  {
    createLocalizedBasisSet<CBSOBuilder>(cur);
  }
  else
  {
    if(FuncType == "Bspline")
    {
      createLocalizedBasisSet<BsplineAOBuilder>(cur);
    }
    if(FuncType == "gto" || FuncType == "GTO")
      createLocalizedBasisSet<GTOBuilder>(cur);
  }
  return true;
}

template<typename COT>
void JastrowBasisBuilder::printRadialFunctors(const std::string& elementType, COT* aoBasis)
{
#if !defined(HAVE_MPI)
  std::string fname(elementType);
  fname.append(".j3.dat");
  std::ofstream fout(fname.c_str());
  int nr=aoBasis->Rnl.size();
  fout << "# number of radial functors = " << nr << std::endl;
  double r=0.0;
  while(r<20)
  {
    fout << r ;
    for(int i=0; i<nr; i++)
      fout << " " << aoBasis->Rnl[i]->evaluate(r,1.0/r);
    fout << std::endl;
    r += 0.013;
  }
#endif
}
}
