//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
/** @file SphericalBasisSet.h
 * @brief A basis set of spherical symmetry associated with a center
 */
#ifndef QMCPLUSPLUS_SPHERICALORBITAL_BASISSET_H
#define QMCPLUSPLUS_SPHERICALORBITAL_BASISSET_H

#include "Particle/DistanceTableData.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Numerics/SphericalTensor.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus {

  template<class ROT, class GT=DummyGrid>
  struct SphericalBasisSet 
  {
    typedef ROT                                                  RadialOrbital_t;
    typedef typename ROT::value_type                             value_type;
    typedef typename OrbitalSetTraits<value_type>::RealType      RealType;
    typedef typename OrbitalSetTraits<value_type>::ValueType     ValueType;
    typedef typename OrbitalSetTraits<value_type>::IndexType     IndexType;
    typedef typename OrbitalSetTraits<value_type>::PosType       PosType;
    typedef typename OrbitalSetTraits<value_type>::GradType      GradType;
    typedef typename OrbitalSetTraits<value_type>::HessType      HessType;
    typedef typename OrbitalSetTraits<value_type>::ValueVector_t ValueVector_t;
    typedef typename OrbitalSetTraits<value_type>::ValueMatrix_t ValueMatrix_t;
    typedef typename OrbitalSetTraits<value_type>::GradVector_t  GradVector_t;
    typedef typename OrbitalSetTraits<value_type>::GradMatrix_t  GradMatrix_t;
    typedef typename OrbitalSetTraits<value_type>::HessVector_t  HessVector_t;
    typedef typename OrbitalSetTraits<value_type>::HessMatrix_t  HessMatrix_t;
    typedef SphericalTensor<RealType,PosType>                    SphericalHarmonics_t;

    ///size of the basis set
    IndexType BasisSetSize;
    ///index of the CurrentCenter
    IndexType CurrentCenter;
    ///offset
    IndexType CurrentOffset;
    ///reference to a DistanceTableData (ion-electron)
    const DistanceTableData* myTable;
    ///spherical tensor unique to this set of SphericalOrbitals
    SphericalHarmonics_t Ylm;
    ///index of the corresponding real Spherical Harmonic with quantum 
    ///numbers \f$ (l,m) \f$
    vector<int> LM;
    /**index of the corresponding radial orbital with quantum 
      numbers \f$ (n,l) \f$ */
    vector<int> NL;
    ///container for the radial grid
    vector<GT*> Grids;
    ///container for the radial orbitals
    vector<ROT*> Rnl;
    ///container for the quantum-numbers
    vector<QuantumNumberType> RnlID;

    ///the constructor
    explicit SphericalBasisSet(int lmax, bool addsignforM=false): Ylm(lmax,addsignforM){ }

    ~SphericalBasisSet() { }

    SphericalBasisSet<ROT,GT>* makeClone() const
    {
      SphericalBasisSet<ROT,GT>* myclone=new SphericalBasisSet<ROT,GT>(*this);
      //for(int i=0; i<Grids.size(); ++i) myclone->Grids[i]=Grids[i]->makeClone();
      //for(int i=0; i<Rnl.size(); ++i) myclone->Rnl[i]=new ROT(*Rnl[i]);
      for(int i=0; i<Rnl.size(); ++i) 
        myclone->Rnl[i]=dynamic_cast<ROT*>(Rnl[i]->makeClone());
      return myclone;
    }

    void checkInVariables(opt_variables_type& active) 
    {
      for(int nl=0; nl<Rnl.size(); nl++) Rnl[nl]->checkInVariables(active);
    }

    void checkOutVariables(const opt_variables_type& active) 
    {
      for(int nl=0; nl<Rnl.size(); nl++) Rnl[nl]->checkOutVariables(active);
    }

    void resetParameters(const opt_variables_type& active) 
    { 
      for(int nl=0; nl<Rnl.size(); nl++) Rnl[nl]->resetParameters(active);
    }

    /** return the number of basis functions
     */
    inline int getBasisSetSize() const { return BasisSetSize;}//=NL.size(); }


    /** implement a BasisSetBase virutal function
     *
     * Use the size of LM to set BasisSetSize
     */
    inline void setBasisSetSize(int n) {
      BasisSetSize=LM.size();
    }

    /** reset the target ParticleSet
     *
     * Do nothing. Leave it to a composite object which owns this
     */
    void resetTargetParticleSet(ParticleSet& P) { }

    ///reset the DistanceTableData (ion-electron)
    inline void setTable(const DistanceTableData* atable) { 
      myTable = atable;
      BasisSetSize=LM.size();
    }

    ///set the current offset
    inline void setCenter(int c, int offset) {
      CurrentCenter=c;
      CurrentOffset=offset;
    }

    inline void 
    evaluateForWalkerMove(int c, int iat, int offset, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) {
      int nn = myTable->M[c]+iat;
      RealType r(myTable->r(nn));
      RealType rinv(myTable->rinv(nn));
      PosType  dr(myTable->dr(nn));

      Ylm.evaluateAll(dr);

      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit)->evaluateAll(r,rinv); ++rit;}
      
      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
	int nl(*nlit);//NL[ib];
	int lm(*lmit);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
	RealType drnloverr(rinv*rnl.dY);
	ValueType ang(Ylm.getYlm(lm));
	PosType gr_rad(drnloverr*dr);
	PosType gr_ang(Ylm.getGradYlm(lm));
        psi[offset]  = ang*rnl.Y;
        dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
	d2psi[offset] = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        ++nlit; ++lmit;++offset;
      }
    }

    inline void
    evaluateForWalkerMove(int c, int iat, int offset, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_Phi) {
      int nn = myTable->M[c]+iat;
      RealType r(myTable->r(nn));
      RealType rinv(myTable->rinv(nn));
      PosType  dr(myTable->dr(nn));

      Ylm.evaluateWithHessian(dr);

      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit)->evaluateAll(r,rinv); ++rit;}

      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
        int nl(*nlit);//NL[ib];
        int lm(*lmit);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
        RealType drnloverr(rinv*rnl.dY);
        ValueType ang(Ylm.getYlm(lm));
        PosType gr_rad(drnloverr*dr);
        PosType gr_ang(Ylm.getGradYlm(lm));
        HessType hess(Ylm.getHessYlm(lm));
        psi[offset]  = ang*rnl.Y;
        dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
// sloppy for now
        RealType temp1=rnl.d2Y*ang*rinv*rinv;
        RealType temp2=drnloverr*ang*rinv*rinv;
        for(int i=0; i<3; i++) {
          grad_grad_Phi[offset](i,i) = (temp1-temp2)*dr(i)*dr(i) 
                    + drnloverr*ang + rnl.Y*hess(i,i)
                    + 2*drnloverr*dr(i)*gr_ang(i);
          for(int j=i+1; j<3; j++) {
            grad_grad_Phi[offset](i,j) = (temp1-temp2)*dr(i)*dr(j)
                    + rnl.Y*hess(i,j)
                    + drnloverr*(dr(i)*gr_ang(j) + dr(j)*gr_ang(i));
            grad_grad_Phi[offset](j,i) = grad_grad_Phi[offset](i,j);
          }
        }
	//d2psi[offset] = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        ++nlit; ++lmit;++offset;
      }
    }

    inline void 
    evaluateForWalkerMove(int source, int first, int nptcl, int offset, ValueMatrix_t& y, GradMatrix_t& dy, ValueMatrix_t& d2y) {

      int nn = myTable->M[source]+first;//first pair of the particle subset
      for(int i=0, iat=first; i<nptcl; i++, iat++, nn++) {
	RealType r(myTable->r(nn));
	RealType rinv(myTable->rinv(nn));
	PosType dr(myTable->dr(nn));

	Ylm.evaluateAll(dr);

        typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
        while(rit != rit_end) {(*rit)->evaluateAll(r,rinv); ++rit;}
	
	int bindex(offset);
        vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 

        ValueType* restrict yptr = y[iat]+offset;
        GradType* restrict dyptr = dy[iat]+offset;
        ValueType* restrict d2yptr = d2y[iat]+offset;
        while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
          int nl(*nlit);//NL[ib];
          int lm(*lmit);//LM[ib];
          const ROT& rnl(*Rnl[nl]);
	  RealType drnloverr = rinv*rnl.dY;
	  ValueType ang = Ylm.getYlm(lm);
	  PosType gr_rad(drnloverr*dr);
	  PosType gr_ang(Ylm.getGradYlm(lm));
	  *yptr++ = ang*rnl.Y;
	  *dyptr++ = ang*gr_rad+rnl.Y*gr_ang;
	  *d2yptr++= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
	  //y(iat,bindex)= ang*rnl.Y;
	  //dy(iat,bindex) = ang*gr_rad+rnl.Y*gr_ang;
	  //d2y(iat,bindex)= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
          ++nlit; ++lmit;++bindex;
        }
      }
    }


    ///evaluate the value, gradient and laplacian of basis functions for the iath-particle
    inline void
    evaluateForWalkerMove(int c, int iat, int offset, Matrix<ValueType>& temp) {
      //RealType r(myTable->Temp[c].r1);
      //RealType rinv(myTable->Temp[c].rinv1);
      //PosType  dr(myTable->Temp[c].dr1);
      //
      int nn = myTable->M[c]+iat;
      RealType r(myTable->r(nn));
      RealType rinv(myTable->rinv(nn));
      PosType  dr(myTable->dr(nn));
      Ylm.evaluateAll(dr);
      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit++)->evaluateAll(r,rinv);}
      
      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 
      ValueType* restrict tptr=temp[offset];
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
        int nl(*nlit++);//NL[ib];
        int lm(*lmit++);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
        RealType drnloverr(rinv*rnl.dY);
        RealType ang(Ylm.getYlm(lm));
        PosType gr_rad(drnloverr*dr);
        PosType gr_ang(Ylm.getGradYlm(lm));
        //PosType g(ang*gr_rad+rnl.Y*gr_ang);
        *tptr++ = ang*rnl.Y;
        *tptr++ = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        *tptr++ = ang*gr_rad[0]+rnl.Y*gr_ang[0];
        *tptr++ = ang*gr_rad[1]+rnl.Y*gr_ang[1];
        *tptr++ = ang*gr_rad[2]+rnl.Y*gr_ang[2];
        //PAOps<RealType,DIM>::copy(ang*rnl.Y, ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang), ang*gr_rad+rnl.Y*gr_ang, temp[offset]); 
        //++nlit; ++lmit; 
      }
    }

    inline void
    evaluateForPtclMove(int source, int iat,  int offset, ValueVector_t& y) {
      RealType r(myTable->Temp[source].r1);
      //RealType rinv(myTable->Temp[source].rinv1);
      RealType rinv(1/r);
      PosType  dr(myTable->Temp[source].dr1);
      Ylm.evaluate(dr);
      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit)->evaluate(r,rinv); ++rit;}
      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
        y[offset++]=Ylm.getYlm(*lmit++)*Rnl[*nlit++]->Y;
      }
    }

    inline void
    evaluateAllForPtclMove(int source, int iat,  int offset, ValueVector_t& y,
        GradVector_t& dy, ValueVector_t& d2y) {
      RealType r(myTable->Temp[source].r1);
      RealType rinv(myTable->Temp[source].rinv1);
      PosType  dr(myTable->Temp[source].dr1);

      Ylm.evaluateAll(dr);
	
      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit)->evaluateAll(r,rinv); ++rit;}
      
      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
	int nl(*nlit);//NL[ib];
	int lm(*lmit);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
	RealType drnloverr(rinv*rnl.dY);
	ValueType ang(Ylm.getYlm(lm));
	PosType gr_rad(drnloverr*dr);
	PosType gr_ang(Ylm.getGradYlm(lm));
	y[offset]= ang*rnl.Y;
	dy[offset] = ang*gr_rad+rnl.Y*gr_ang;
	d2y[offset]= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        ++nlit; ++lmit;++offset;
      }
    }


    inline void
    evaluate(RealType r, RealType rinv, const PosType& dr, int offset, ValueVector_t& psi) {
      Ylm.evaluate(dr);
      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit)->evaluate(r,rinv); ++rit;}
      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
        psi[offset++]=Ylm.getYlm(*lmit++)*Rnl[*nlit++]->Y;
      }
    }

    inline void 
    evaluate(RealType r, RealType rinv, const PosType& dr, int offset, ValueVector_t& y, 
        GradVector_t& dy, ValueVector_t& d2y) {

      Ylm.evaluateAll(dr);
	
      typename vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end) {(*rit)->evaluateAll(r,rinv); ++rit;}
      
      vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin()); 
      while(nlit != nlit_end) { //for(int ib=0; ib<NL.size(); ib++, offset++) {
	int nl(*nlit);//NL[ib];
	int lm(*lmit);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
	RealType drnloverr(rinv*rnl.dY);
	ValueType ang(Ylm.getYlm(lm));
	PosType gr_rad(drnloverr*dr);
	PosType gr_ang(Ylm.getGradYlm(lm));
	y[offset]= ang*rnl.Y;
	dy[offset] = ang*gr_rad+rnl.Y*gr_ang;
	d2y[offset]= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        ++nlit; ++lmit;++offset;
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
