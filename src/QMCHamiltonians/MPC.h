//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim and Ken Esler
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
#ifndef QMCPLUSPLUS_MPC_H
#define QMCPLUSPLUS_MPC_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   *\brief Calculates the Model Periodic Coulomb potential using PBCs
   */

  class UBspline_3d_d;

  class MPC: public QMCHamiltonianBase {
  private:
    UBspline_3d_d *VlongSpline;
    void compute_g_G(double &g_0_N, vector<double> &g_G_N, int N);
    void init_f_G();
    void init_spline();
    double Ecut;
    vector<PosType> Gvecs;

  public:
    ParticleSet& PtclRef;
    DistanceTableData* d_aa;

    // Store the average electron charge density in reciprocal space
    vector<ComplexType> RhoAvg_G;
    vector<RealType> f_G;
    // The G=0 component
    double f_0;

    bool FirstTime;
    int NumSpecies;
    int ChargeAttribIndx;
    int MemberAttribIndx;
    int NParticles;
    RealType myConst;
    RealType myRcut;
    vector<RealType> Zat,Zspec; 
    vector<int> NofSpecies;

    MPC(ParticleSet& ref, double cutoff);

    /// copy constructor
    MPC(const MPC& c);
    
    ~MPC();

    void resetTargetParticleSet(ParticleSet& P) { }

    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    /** Do nothing */
    bool put(xmlNodePtr cur);

    bool get(std::ostream& os) const {
      os << "MPC potential: " << PtclRef.getName();
      return true;
    }

    QMCHamiltonianBase* clone(ParticleSet& qp, TrialWaveFunction& psi);

    void initBreakup();

    Return_t evalSR();
    Return_t evalLR();
    Return_t evalConsts();

  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: esler $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: MPC.h 1581 2007-01-04 16:02:14Z esler $ 
 ***************************************************************************/

