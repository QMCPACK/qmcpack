//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
    const string& functype, bool usespline):
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
    string cname((const char*)(cur->name));
    if(cname == "atomicBasisSet")
    {
      const xmlChar* eptr=xmlGetProp(cur,(const xmlChar*)"elementType");
      if(eptr == NULL)
      {
        app_error() << "   Missing elementType attribute of atomicBasisSet.\n Abort at MOBasisBuilder::put " << endl;
        OHMMS::Controller->abort();
      }
      string elementType((const char*)eptr);
      map<string,BasisSetBuilder*>::iterator it = aoBuilders.find(elementType);
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
void JastrowBasisBuilder::printRadialFunctors(const string& elementType, COT* aoBasis)
{
#if !defined(HAVE_MPI)
  string fname(elementType);
  fname.append(".j3.dat");
  ofstream fout(fname.c_str());
  int nr=aoBasis->Rnl.size();
  fout << "# number of radial functors = " << nr << endl;
  double r=0.0;
  while(r<20)
  {
    fout << r ;
    for(int i=0; i<nr; i++)
      fout << " " << aoBasis->Rnl[i]->evaluate(r,1.0/r);
    fout << endl;
    r += 0.013;
  }
#endif
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1623 $   $Date: 2007-01-18 18:01:15 -0600 (Thu, 18 Jan 2007) $
 * $Id: JastrowBasisBuilder.h 1623 2007-01-19 00:01:15Z jnkim $
 ***************************************************************************/
