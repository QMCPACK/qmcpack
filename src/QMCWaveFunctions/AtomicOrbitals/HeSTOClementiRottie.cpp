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
    
    
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

/** constructor
 *@param els the electron particleset
 *@param wfs trial wavefuntion to which determinant terms are added
 *@param ions the ion particleset
 *
 *@note A DiracDeterminant is a single determinant, a SlaterDeterminant is
 *the product of DiracDeterminants while a MultiDeterminant is a linear
 *combination of SlaterDeterminants
 */
HePresetHFBuilder::HePresetHFBuilder(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions):
  OrbitalBuilderBase(els,wfs)
{
  //the electron-ion distancetable
  DistanceTableData* d_ei = DistanceTable::add(ions,els);
  typedef DiracDeterminant<HePresetHF> Det_t;
  //pointer to the Helium single particle orbitals
  HePresetHF* heorb = new HePresetHF;
  for(int i=0; i<heorb->N; i++)
  {
    LOGMSG(" Slater Component (n,zeta,c)= 1 " << heorb->Z[i] << " " << heorb->C[i])
  }
  //set the distance table
  heorb->setTable(d_ei);
  //determinant for up-electron
  Det_t *DetU = new Det_t(*heorb,0);
  //determinant for particle 0 and size 1
  DetU->set(0,1);
  ///determinant for down-electron
  Det_t* DetD = new Det_t(*heorb,1);
  //determinant for particle 1 and size 1
  DetD->set(1,1);
  //add a SlaterDeterminant
  SlaterDeterminant<HePresetHF> *asymmpsi=new SlaterDeterminant<HePresetHF>;
  //add the DiracDeterminants to the SlaterDeterminant
  asymmpsi->add(DetU);
  asymmpsi->add(DetD);
  asymmpsi->setBasisSet(new HePresetHF::BasisSet_t);
  //add the SlaterDeterminant to the trial wavefuntion
  targetPsi.addOrbital(asymmpsi);
}

bool
HePresetHFBuilder::put(xmlNodePtr cur)
{
  //DistanceTableData* d_ei = NULL;
  //cur = cur->xmlChildrenNode;
  //while(cur != NULL) {
  //  std::string cname((const char*)(cur->name));
  //  if(cname == dtable_tag) {
  //    std::string dtable((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
  //    dtable.append((const char*)(xmlGetProp(cur,(const xmlChar *)"target")));
  //    d_ei= DistanceTable::getTable(dtable.c_str());
  //    LOGMSG("DistanceTable is selected = " << dtable)
  //  }
  //  cur = cur->next;
  //}
  //typedef DiracDeterminant<HePresetHF> Det_t;
  //HePresetHF* heorb = new HePresetHF;
  //heorb->setTable(d_ei);
  //Det_t *DetU = new Det_t(*heorb,0);
  //DetU->set(0,1);
  //Det_t* DetD = new Det_t(*heorb,1);
  //DetD->set(1,1);
  //SlaterDeterminant<HePresetHF> *asymmpsi=new SlaterDeterminant<HePresetHF>;
  //asymmpsi->add(DetU);
  //asymmpsi->add(DetD);
  //targetPsi.add(asymmpsi);
  return true;
}
}
