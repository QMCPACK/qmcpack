//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.h"

qmcplusplus::HFAtomicSTOSet::HFAtomicSTOSet(): Ylm(0)
{
  RealType C[] = {0.76838,0.22346,0.04082,-0.00994,0.00230};
  RealType screen[] = {1.41714, 2.37682, 4.39628, 6.52699, 7.94252};
  STONorm<RealType> anorm(2);
  for(int i=0; i<5; i++)
  {
    RnlPool.push_back(new RadialOrbital_t(0,screen[i],anorm(0,screen[i])));
  }
  Orbital.push_back(new SPO_t(Ylm.index(0,0), Ylm, RnlPool,C));
}

// void qmcplusplus::PresetHFBuilder(xmlNodePtr cur,TrialWaveFunction& wfs_ref) {

//   LOGMSG("Calling PresetHFBuilder to initialize a trial wave fucntion")
//   DistanceTableData* d_ei = DistanceTable::getTable(1);
//   typedef DiracDeterminant<HFAtomicSTOSet> Det_t;

//   STONorm<RealType> anorm(4);
//   Det_t *DetU = NULL, *DetD = NULL;
//   ///Test case using default constructor for He
//   HFAtomicSTOSet* he = new HFAtomicSTOSet;
//   he->setTable(d_ei);
//   DetU = new Det_t(*he,0);
//   DetU->set(0,1);
//   DetD = new Det_t(*he,1);
//   DetD->set(1,1);

//   SlaterDeterminant<HFAtomicSTOSet>
//     *asymmpsi = new SlaterDeterminant<HFAtomicSTOSet>;
//   if(DetU) asymmpsi->add(DetU);
//   if(DetD) asymmpsi->add(DetD);
//   wfs_ref.add(asymmpsi,d_ei);
// }
