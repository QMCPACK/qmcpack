//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
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
#ifndef OHMMS_QMC_GE_COREPOLPOTENTIAL_H
#define OHMMS_QMC_GE_COREPOLPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace ohmmsqmc {

  /**
     The effective core polarization
     \f[
     V_{CPP} = -\frac{1}{2}\sum_{C}\alpha_C {\bf f}_C \cdot {\bf f}_C
     \f]
     Electric field which acts on core \f$C\f$ due to the charges of 
     valence electrons \f$i\f$ and the other cores \f$C'\f$
     \f[
     {\bf f}_C = \sum_i \frac{{\bf r}_{Ci}}{r_{Ci}^3}C(r_{Ci},\rho_C)
     -\sum_{C' \ne C} \frac{{\bf R}_{CC'}}{R_{CC'}^3}Z_{C'} = 
     {\bf f}_C^e + {\bf f}_C^n
     \f]
     \f[
     {\bf r}_{Ci} = {\bf r}_i - {\bf r}_C
     \f]
     \f[
     {\bf R}_{CC'} = {\bf R}_{C'} - {\bf R}_{C}
     \f]

     \f$ C(r_{Ci},\rho_C) \f$ is a cut-off function for \f$ {\bf f}_C^e \f$ 
     with an adjustable parameter \f$ \rho_C. \f$

     \f[
     V_{CPP} = -\frac{1}{2}\sum_C \left\{ \:\:\: \sum_i 
     \frac{1}{r_{Ci}^4}C^2(r_{Ci},\rho_C) + 
     \sum_{i \ne j} \frac{{\bf r}_{Ci} \cdot {\bf r}_{Ci}}{r^3_{Ci}r^3_{Ci}}
     C(r_{Ci},\rho_C) C^2(r_{Cj},\rho_C) \right. 
     \f]
     \f[ 
     -2 \left. \sum_i \sum_{C' \ne C} \frac{{\bf r}_{Ci} \cdot 
     {\bf R}_{CC'}}{r^3_{Ci}R^3_{CC'}} Z_{C'}C(r_{Ci},\rho_C) 
     + \left| \sum_{C' \ne C}  \frac{{\bf R}_{CC'}}{R^3_{CC'}} Z_{C'} 
     \right|^2 \:\:\: \right\}
     \f]
  */
  struct GeCorePolPotential: public QMCHamiltonianBase {
    ///the number of ions 
    int nCenters;
    ///the number of electrons
    int nParticles;
    RealType alpha, r_b, r_binv;
    RealType eCoreCore;
    RealType C;
    ///the ion-electron DistanceTable
    DistanceTableData* d_ie;
    ///the ion-ion DistanceTable
    DistanceTableData* d_ii;
    ///CoreCoef(C) = 1.0 if C=Ge, 0.0 for all other ions
    vector<bool> CoreCoef;
    ///CoreCoreDipole(C,C') \f$= \frac{Z_{C'} {\bf R_{CC'}}}{R_{CC'}^3}\f$
    vector<PosType> CoreCoreDipole;
    ///ElCoreDipole(C,i) \f$= \frac{{\bf r_{Ci}}f({\bar{r_{bCi}}}{r_{Ci}^3}\f$
    Matrix<PosType> CoreElDipole;

    ///constructor
    GeCorePolPotential(ParticleSet& ions, ParticleSet& els);

    ~GeCorePolPotential();

    ValueType evaluate(ParticleSet& P);

    inline RealType fcpp(RealType z) {
      return pow((1.0-exp(-1.0*z*z)),2);
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

