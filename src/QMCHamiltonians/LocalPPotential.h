//////////////////////////////////////////////////////////////////
// (c) Copyright 2004 by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_LOCALPPOTENTIAL_H
#define QMCPLUSPLUS_LOCALPPOTENTIAL_H
#include <fstream>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
//#include "Numerics/HDFNumericAttrib.h"

namespace qmcplusplus {


  /** @ingroup hamiltonian
   * \brief Evaluate the local potentials (either pseudo or full core) around each ion.
   */

  struct LocalPPotential: public QMCHamiltonianBase {

    typedef OneDimGridBase<ValueType> GridType;
    typedef OneDimGridFunctor<ValueType> LocalPotentialType;

  /**
   *\brief Contains a set of radial grid potentials around a 
   center.  
   *
   Evaluates the grid potential given the electron-
   ion distances.
  */
  struct RadialPotentialSet {
    ///the radial potentials
    vector<LocalPotentialType*> lpp_m;
    ///the radial grids
    vector<GridType*> grid_m;

    ///destructor
    ~RadialPotentialSet();

    ///add a new grid and radial grid potential
    void add(GridType* agrid, LocalPotentialType* pp) {
      grid_m.push_back(agrid);
      lpp_m.push_back(pp);
    }

    /*!
     *\param d_table the distance table containing electron-
     ion distances  
     *\param iat the index of the nucleus
     *\return \f$ \sum_{I=1}^{N_C}\sum_{i=1}^{N_P} V(r_{iI}) \f$
     *\brief Evaluates the set of radial grid potentials around
     a single ion.
    */
    inline ValueType evaluate(DistanceTableData* d_table, int iat) {
      ValueType esum = 0.0;     
      for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++){
        //evaluate takes care of the grid 
        //vector<GridType*>::iterator gid(grid_m.begin()), gid_end(grid_m.end());
        //while(gid != gid_end) {
        //  (*gid)->locate(d_table->r(nn));++gid;
        //}
        vector<LocalPotentialType*>::iterator lit(lpp_m.begin()),lit_end(lpp_m.end());
        while(lit != lit_end) {
	  esum += (*lit)->evaluate(d_table->r(nn),d_table->rinv(nn)); ++lit;
        }
	//for(int ig=0; ig<grid_m.size(); ig++) 
	//  grid_m[ig]->index(d_table->r(nn));
	//for(int ip=0;ip<lpp_m.size(); ip++) 
	//  esum += lpp_m[ip]->evaluate(d_table->r(nn),d_table->rinv(nn));
      }
      return esum;
    }
  };

    ///the distance table containing electron-nuclei distances  
    DistanceTableData* d_table;
    ///the set of local-potentials (one for each ion)
    vector<RadialPotentialSet*> PP;
    ///unique index for each ion
    const ParticleSet::ParticleIndex_t& Centers;

    LocalPPotential(ParticleSet& ions, ParticleSet& els);
    
    ~LocalPPotential();

    inline Return_t evaluate(ParticleSet& P) {
      Value=0.0;
      //loop over all the ions
      for(int iat=0; iat<Centers.size(); iat++) {
	//evaluate the corresponding potential for the ion
        Value += PP[ Centers[iat] ]->evaluate(d_table,iat);
      }
      return Value;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

