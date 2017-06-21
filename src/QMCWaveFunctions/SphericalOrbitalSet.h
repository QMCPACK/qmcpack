//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SPHERICALORBITALSET_H
#define QMCPLUSPLUS_SPHERICALORBITALSET_H

#include "Particle/DistanceTableData.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Numerics/SphericalTensor.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus
{

/**Class to represent a set of spherical orbitals centered at a common origin origin.
 *
 *Each basis orbital is represented by
 \f[
  \phi_i({\bf r_{jC}}) =  R_{n_i l_i}(r_{jC})
  \Re (Y_{l_i}^{m_i}({\bf \hat{r_{jC}}}))
  \f]
  Here \f${\bf r_{jC}} = {\bf r_j - R_C}\f$ denotes the position of
  particle \f$ j \f$ relative to the center \f$ C \f$.
  *
  The template (R)adial(O)rbital(T)ype should provide
  <ul>
  <li> evaluate(r,rinv)
  <li> Y
  <li> dY
  <li> d2Y
  </ul>
  The template (G)rid(T)ype should provide
  <ul>
  <li> index(r)
  </ul>
  Examples of ROT being GenericSTO and OneDimGridFunctor
  *
  Examples of GT being LogGrid, LinearGrid, LogGridZero
  and NumericalGrid
*/
template<class ROT, class GT=DummyGrid>
struct SphericalOrbitalSet
{

  typedef ROT                                                  RadialOrbital_t;
  typedef typename ROT::value_type                             value_type;
  typedef typename OrbitalSetTraits<value_type>::RealType      RealType;
  typedef typename OrbitalSetTraits<value_type>::ValueType     ValueType;
  typedef typename OrbitalSetTraits<value_type>::IndexType     IndexType;
  typedef typename OrbitalSetTraits<value_type>::PosType       PosType;
  typedef typename OrbitalSetTraits<value_type>::GradType      GradType;
  typedef typename OrbitalSetTraits<value_type>::ValueVector_t ValueVector_t;
  typedef typename OrbitalSetTraits<value_type>::ValueMatrix_t ValueMatrix_t;
  typedef typename OrbitalSetTraits<value_type>::GradVector_t  GradVector_t;
  typedef typename OrbitalSetTraits<value_type>::GradMatrix_t  GradMatrix_t;
  typedef SphericalTensor<RealType,PosType>           SphericalHarmonics_t;
  typedef VarRegistry<RealType>                       OptimizableSetType;

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
  std::vector<int> LM;
  /**index of the corresponding radial orbital with quantum
    numbers \f$ (n,l) \f$ */
  std::vector<int> NL;
  ///container for the radial grid
  std::vector<GT*> Grids;
  ///container for the radial orbitals
  std::vector<ROT*> Rnl;
  ///container for the quantum-numbers
  std::vector<QuantumNumberType> RnlID;

  ///the constructor
  explicit SphericalOrbitalSet(int lmax, bool addsignforM=false): Ylm(lmax,addsignforM) { }

  ~SphericalOrbitalSet() { }

  ///reset the internal values
  void resetParameters(OptimizableSetType& optVariables)
  {
    for(int nl=0; nl<Rnl.size(); nl++)
      Rnl[nl]->resetParameters(optVariables);
  }

  /** return the number of basis functions
   */
  inline int getBasisSetSize() const
  {
    return BasisSetSize;
  }//=NL.size(); }


  /** implement a BasisSetBase virutal function
   *
   * Use the size of LM to set BasisSetSize
   */
  inline void setBasisSetSize(int n)
  {
    BasisSetSize=LM.size();
  }

  /** reset the target ParticleSet
   *
   * Do nothing. Leave it to a composite object which owns this
   */
  void resetTargetParticleSet(ParticleSet& P) { }

  ///reset the DistanceTableData (ion-electron)
  inline void setTable(DistanceTableData* atable)
  {
    myTable = atable;
    BasisSetSize=LM.size();
  }

  ///set the current offset
  inline void setCenter(int c, int offset)
  {
    CurrentCenter=c;
    CurrentOffset=offset;
  }

  template<class VV, class GV>
  inline void
  evaluateBasis(int c, int iat, int offset, VV& psi, GV& dpsi, VV& d2psi)
  {
    int nn = myTable->M[c]+iat;
    RealType r(myTable->r(nn));
    RealType rinv(myTable->rinv(nn));
    PosType  dr(myTable->dr(nn));
    Ylm.evaluateAll(dr);
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit)->evaluateAll(r,rinv);
      ++rit;
    }
    std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
    while(nlit != nlit_end)
      //for(int ib=0; ib<NL.size(); ib++, offset++) {
    {
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
      ++nlit;
      ++lmit;
      ++offset;
    }
  }


  template<class VM>
  inline void
  evaluate(int source, int iat,  int offset, VM& y)
  {
    RealType r(myTable->Temp[source].r1);
    RealType rinv(myTable->Temp[source].rinv1);
    PosType  dr(myTable->Temp[source].dr1);
    Ylm.evaluate(dr);
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit)->evaluate(r,rinv);
      ++rit;
    }
    std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
    while(nlit != nlit_end)
      //for(int ib=0; ib<NL.size(); ib++, offset++) {
    {
      y(0,offset++)=Ylm.getYlm(*lmit++)*Rnl[*nlit++]->Y;
    }
  }

  template<class VM, class GM>
  inline void
  evaluate(int source, int iat,  int offset, VM& y, GM& dy, VM& d2y)
  {
    RealType r(myTable->Temp[source].r1);
    RealType rinv(myTable->Temp[source].rinv1);
    PosType  dr(myTable->Temp[source].dr1);
    Ylm.evaluateAll(dr);
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit)->evaluateAll(r,rinv);
      ++rit;
    }
    std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
    while(nlit != nlit_end)
      //for(int ib=0; ib<NL.size(); ib++, offset++) {
    {
      int nl(*nlit);//NL[ib];
      int lm(*lmit);//LM[ib];
      const ROT& rnl(*Rnl[nl]);
      RealType drnloverr(rinv*rnl.dY);
      ValueType ang(Ylm.getYlm(lm));
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(Ylm.getGradYlm(lm));
      y(0,offset)= ang*rnl.Y;
      dy(0,offset) = ang*gr_rad+rnl.Y*gr_ang;
      d2y(0,offset)= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

  /** For the center with index (source), evaluate all the basis functions beginning with index (offset).
   * @param source index of the center \f$ I \f$
   * @param first index of the first particle
   * @param nptcl number of particles
   * @param offset index of the first basis function
   * @param y return vector \f$ y[i,j] = \phi_j({\bf r_i-R_I}) \f$
   * @param dy return vector \f$ dy[i,j] =  {\bf \nabla}_i \phi_j({\bf r_i-R_I}) \f$
   * @param d2y return vector \f$ d2y[i,j] = \nabla^2_i \phi_j({\bf r_i-R_I}) \f$
   *
     The results are stored in the matrices:
     \f[ y[i,k] = \phi_k(r_{ic}) \f]
     \f[ dy[i,k] = {\bf \nabla}_i \phi_k(r_{ic}) \f]
     \f[ d2y[i,k] = \nabla_i^2 \phi_k(r_{ic}) \f]
     *
     The basis functions can be written in the form
     \f[ \phi({\bf R}) = F(r)G(r,\theta,\phi), \f]
     where \f$ F(r) = \frac{R(r)}{r^l} \f$ is related to the radial orbital
     and \f$ G(r,\theta,\phi) = r^lS_l^m(\theta,\phi) \f$ is related to
     the real spherical harmonic.  Evaluating the gradient and Laplacian
     leads to
     \f[  {\bf \nabla} \phi({\bf R}) =
     ({\bf \nabla } F(r)) G(r,\theta,\phi) +
     F(r) ({\bf \nabla} G(r,\theta,\phi))  \f]
     \f[
     {\bf \nabla} \phi({\bf R}) =
     \frac{dF(r)}{dr} G(r,\theta,\phi)\:{\bf \hat{r}} +
     F(r) {\bf \nabla} G(r,\theta,\phi)
     \f]
     and
     \f[
     \nabla^2 \phi({\bf R}) = (\nabla^2 F) G +2 \nabla F \cdot \nabla G
     + F (\nabla^2 G) \f]
     where
     \f[ \nabla^2 F(r) = \frac{2}{r}\frac{dF(r)}{dr} + \frac{d^2F(r)}{dr^2}
     \mbox{,     } \nabla F(r) = \frac{dF(r)}{dr}\:{\bf \hat{r}}
     \f]
     and
     \f[  \nabla^2 G(r,\theta,\phi) = \frac{1}{r^2}\frac{\partial}
     {\partial r} \left( r^2 \frac{\partial G}{\partial r}\right)
     - \frac{l(l+1)G}{r^2} = 0.
     \f]
  */
  template<class VM, class GM>
  inline void
  evaluate(int source, int first, int nptcl, int offset, VM& y, GM& dy, VM& d2y)
  {
    int nn = myTable->M[source]+first;//first pair of the particle subset
    for(int i=0, iat=first; i<nptcl; i++, iat++, nn++)
    {
      RealType r(myTable->r(nn));
      RealType rinv(myTable->rinv(nn));
      PosType dr(myTable->dr(nn));
      Ylm.evaluateAll(dr);
      typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end)
      {
        (*rit)->evaluateAll(r,rinv);
        ++rit;
      }
      int bindex(offset);
      std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
      while(nlit != nlit_end)
        //for(int ib=0; ib<NL.size(); ib++, offset++) {
      {
        int nl(*nlit);//NL[ib];
        int lm(*lmit);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
        RealType drnloverr = rinv*rnl.dY;
        ValueType ang = Ylm.getYlm(lm);
        PosType gr_rad(drnloverr*dr);
        PosType gr_ang(Ylm.getGradYlm(lm));
        y(iat,bindex)= ang*rnl.Y;
        dy(iat,bindex) = ang*gr_rad+rnl.Y*gr_ang;
        d2y(iat,bindex)= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        ++nlit;
        ++lmit;
        ++bindex;
      }
      //for(int ib=0; ib<NL.size(); ib++, bindex++) {
      //  int nl = NL[ib];
      //  int lm = LM[ib];
      //  RealType drnloverr = rinv*Rnl[nl]->dY;
      //  ValueType ang = Ylm.getYlm(lm);
      //  PosType gr_rad = drnloverr*dr;
      //  PosType gr_ang = Ylm.getGradYlm(lm);
      //  y(iat,bindex)= ang*Rnl[nl]->Y;
      //  dy(iat,bindex) = ang*gr_rad+Rnl[nl]->Y*gr_ang;
      //  d2y(iat,bindex)= ang*(2.0*drnloverr+Rnl[nl]->d2Y) + 2.0*dot(gr_rad,gr_ang);
      //}
    }
  }


};

}
#endif

