//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/** @file LocalizedBasisSet.h
 * @author Jeongnim Kim
 * @brief A derived class from BasisSetBase
 *
 * This is intended as a replacement for MolecularOrbitalBase and
 * any other localized basis set.
 */
#ifndef QMCPLUSPLUS_LOCALIZEDBASISSET_H
#define QMCPLUSPLUS_LOCALIZEDBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus {

  /** Class for a molecular orbital basis
   *
   *The molecular orbital \f$ \psi_i \f$ can be written as a linear
   *combination of basis functions \f$ \{\phi\} \f$ such that
   \f[
   \psi_i ({\bf r}_j) = \sum_I \sum_k C_{ikI} \phi_{ikI}({\bf r}_j-{\bf R}_I).
   \f]
   *This class performs the evaluation of the basis functions and their
   *derivatives for each of the N-particles in a configuration.  All that 
   *is required to generate the actual molecular orbitals is to multiply
   *by the coefficient matrix.
   *
   *The template (C)entered(O)rbital(T)ype should provide the fuctions
   <ul>
   <li> evaluate(int source, int first, int nptcl, int offset, 
   VM& y, GM& dy, VM& d2y) {
   </ul>
   *An example being SphericalOrbitalSet
   */
  template<class COT>
  struct LocalizedBasisSet: public BasisSetBase {

    ///Reference to the center
    const ParticleSet& CenterSys;
    ///number of centers, e.g., ions
    int NumCenters;
    ///number of quantum particles
    int NumTargets;

    /** container to store the offsets of the basis functions
     *
     * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
     */
    vector<int>  BasisOffset;

    /** container of the pointers to the Atomic Orbitals
     *
     * AO[i] returns a Centered Orbital for an ion i
     */
    vector<COT*> LOBasis;

    /** container for the pointers to the Atomic Orbitals 
     *
     * the size of this container being determined by the number 
     * of unique centers
     */
    vector<COT*> LOBasisSet;

    ///the distance table (ion-electron)
    const DistanceTableData* myTable;

    /** constructor
     * @param ions ionic system
     * @param els electronic system
     */
    LocalizedBasisSet(ParticleSet& ions, ParticleSet& els): CenterSys(ions), myTable(0){ 
      myTable = DistanceTable::add(ions,els);
      NumCenters=CenterSys.getTotalNum();
      NumTargets=els.getTotalNum();
      LOBasis.resize(NumCenters,0);
      LOBasisSet.resize(CenterSys.getSpeciesSet().getTotalNum(),0);
      BasisOffset.resize(NumCenters+1);
    }

    /**
     @param atable the distance table (ion-electron)
     @brief Assign the distance table (ion-electron) and
     *determine the total number of basis states.
    */
    void setBasisSetSize(int nbs) { 
      if(nbs == BasisSetSize) return;

      if(myTable ==0) {
        app_error() << "LocalizedBasisSet cannot function without a distance table. Abort" << endl;
      }

      //reset the distance table for the atomic orbitals
      for(int i=0; i<LOBasisSet.size(); i++) LOBasisSet[i]->setTable(myTable);
      //evaluate the total basis dimension and offset for each center
      BasisOffset[0] = 0;
      for(int c=0; c<NumCenters; c++){
	BasisOffset[c+1] = BasisOffset[c]+LOBasis[c]->getBasisSetSize();
      }
      BasisSetSize = BasisOffset[NumCenters];

      resize(NumTargets);
    }

    void resetParameters(VarRegistry<RealType>& optVariables) 
    {
      //reset each unique basis functions
      for(int i=0; i<LOBasisSet.size(); i++) LOBasisSet[i]->resetParameters(optVariables);
    }
    
    /** reset the distance table with a new target P
     */
    void resetTargetParticleSet(ParticleSet& P) {
      LOGMSG("LocalizedBasisSet::resetTargetParticleSet")
      myTable = DistanceTable::add(CenterSys,P);
      for(int i=0; i<LOBasisSet.size(); i++) LOBasisSet[i]->setTable(myTable);
    }

    inline void 
    evaluateForWalkerMove(const ParticleSet& P) {
      for(int c=0; c<NumCenters;c++) {
        LOBasis[c]->evaluateForWalkerMove(c,0,P.getTotalNum(),BasisOffset[c],Y,dY,d2Y);
      }
    }

    inline void 
    evaluateForWalkerMove(const ParticleSet& P, int iat) {
      for(int c=0; c<NumCenters;c++) {
	LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,d2Phi);
        //int nn = myTable->M[c]+iat;
	//LOBasis[c]->evaluate(myTable->r(nn),myTable->rinv(nn), myTable->dr(nn), 
        //    BasisOffset[c],Phi,dPhi,d2Phi);
      }
    }

    inline void 
    evaluateForPtclMove(const ParticleSet& P, int iat)  {
      for(int c=0; c<NumCenters;c++) {
	//LOBasis[c]->evaluate(myTable->Temp[c].r1,myTable->Temp[c].rinv1, myTable->Temp[c].dr1, 
        //    BasisOffset[c],Phi);
	LOBasis[c]->evaluateForPtclMove(c,iat,BasisOffset[c],Phi);
      }
    }

    inline void 
    evaluateAllForPtclMove(const ParticleSet& P, int iat)  {
      for(int c=0; c<NumCenters;c++) {
	//LOBasis[c]->evaluate(myTable->Temp[c].r1,myTable->Temp[c].rinv1, myTable->Temp[c].dr1, 
        //    BasisOffset[c],Phi,dPhi,d2Phi);
	LOBasis[c]->evaluateAllForPtclMove(c,iat,BasisOffset[c],Phi,dPhi,d2Phi);
      }
    }

    /** add a new set of Centered Atomic Orbitals
     * @param icenter the index of the center
     * @param aos a set of Centered Atomic Orbitals
     */
    void add(int icenter, COT* aos) {
      aos->setTable(myTable);
      LOBasisSet[icenter]=aos;
      for(int i=0; i<NumCenters; i++) {
        if(CenterSys.GroupID[i] == icenter) LOBasis[i]=aos;
      }
    }

  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

