//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_SHO_BASIS_BUILDER_H
#define QMCPLUSPLUS_SHO_BASIS_BUILDER_H

#include <QMCWaveFunctions/HarmonicOscillator/SHOSet.h>
#include <QMCWaveFunctions/BasisSetBase.h>
#include <QMCWaveFunctions/SPOSetInfo.h>

namespace qmcplusplus
{

  struct SHOSetBuilder : public BasisSetBuilder
  {

    //enum{DIM=OHMMS_DIM}

    ParticleSet& Ps;

    RealType length;
    RealType mass;
    RealType energy;
    PosType  center;

    int nstates;
    int nmax;
    TinyVector<int,DIM> ind_dims;

    SPOSetInfoSimple<SHOState> basis_states;

    //construction/destruction
    SHOSetBuilder(ParticleSet& P);

    ~SHOSetBuilder();

    //reset parameters
    void reset();

    //BasisSetBuilder interface
    SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);

    SPOSetBase* createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input);
    
    //unneeded BasisSetBuilder interface functions
    bool put(xmlNodePtr cur)
    { 
      return true; 
    }

    //local functions
    void update_basis_states(int smax);
    void report(const std::string& pad="");
  };

}


#endif
