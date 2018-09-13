//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SphericalBasisSet.h
 * @brief A basis set of spherical symmetry associated with a center
 */
#ifndef QMCPLUSPLUS_SPHERICALORBITAL_BASISSET_H
#define QMCPLUSPLUS_SPHERICALORBITAL_BASISSET_H

#include "Particle/DistanceTableData.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Numerics/SphericalTensor.h"
#include "Numerics/CartesianTensor.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus
{

template<class ROT, class GT=DummyGrid>
struct SphericalBasisSet
{
  typedef ROT                                                  RadialOrbital_t;
  typedef GT                                                   GridType_t;
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
  typedef TinyVector<HessType, 3>                    GGGType;
  typedef Vector<GGGType>                            GGGVector_t;
  typedef Matrix<GGGType>                            GGGMatrix_t;
  typedef CartesianTensor<RealType,PosType>          CartesianHarmonics_t;
  typedef SphericalTensor<RealType,PosType>          SphericalHarmonics_t;

  ///size of the basis set
  IndexType BasisSetSize;
  ///index of the CurrentCenter
  IndexType CurrentCenter;
  ///offset
  IndexType CurrentOffset;
  ///reference to a DistanceTableData (ion-electron)
  const DistanceTableData* myTable;
  ///bool to chosse Spehrical/Cartesian should be templated
  bool useCartesian;
  ///spherical tensor unique to this set of SphericalOrbitals
  //AngularFunction_t* Ylm;
  SphericalHarmonics_t Ylm;
// mmorales: HACK HACK HACK!!!
//  to avoid having to template this function, which will change
//  all the builders
// remember to add laplacian of angular piece with cartesian gaussian
  CartesianHarmonics_t XYZ;
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
  explicit SphericalBasisSet(int lmax, bool addsignforM=false, bool useXYZ=false):Ylm(lmax,addsignforM),XYZ(useXYZ?lmax:0),useCartesian(useXYZ) {}

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
    for(int nl=0; nl<Rnl.size(); nl++)
      Rnl[nl]->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    for(int nl=0; nl<Rnl.size(); nl++)
      Rnl[nl]->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    for(int nl=0; nl<Rnl.size(); nl++)
      Rnl[nl]->resetParameters(active);
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
  inline void setTable(const DistanceTableData* atable)
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

  inline void
  evaluateForWalkerMove(int c, int iat, int offset, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    int nn = myTable->M[c]+iat;
    RealType r(myTable->r(nn));
    RealType rinv(myTable->rinv(nn));
    PosType  dr(myTable->dr(nn));
    if(useCartesian)
    {
      XYZ.evaluateAll(dr);
    }
    else
    {
      Ylm.evaluateAll(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<RealType>& laplYlm = useCartesian?XYZ.laplXYZ:Ylm.laplYlm;
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      psi[offset]  = ang*rnl.Y;
      dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
      d2psi[offset] = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang) + rnl.Y*laplYlm[lm];
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

  inline void
  evaluateForWalkerMove(int c, int iat, int offset, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_Phi)
  {
    int nn = myTable->M[c]+iat;
    RealType r(myTable->r(nn));
    RealType rinv(myTable->rinv(nn));
    PosType  dr(myTable->dr(nn));
    if(useCartesian)
    {
      XYZ.evaluateWithHessian(dr);
    }
    else
    {
      Ylm.evaluateWithHessian(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<Tensor<RealType,3> >& hessYlm = useCartesian?XYZ.hessXYZ:Ylm.hessYlm;
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      HessType hess(hessYlm[lm]);
      psi[offset]  = ang*rnl.Y;
      dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
// sloppy for now
      RealType temp1=rnl.d2Y*ang*rinv*rinv;
      RealType temp2=drnloverr*ang*rinv*rinv;
      for(int i=0; i<3; i++)
      {
        grad_grad_Phi[offset](i,i) = (temp1-temp2)*dr[i]*dr[i]
                                     + drnloverr*ang + rnl.Y*hess(i,i)
                                     + 2*drnloverr*dr[i]*gr_ang[i];
        for(int j=i+1; j<3; j++)
        {
          grad_grad_Phi[offset](i,j) = (temp1-temp2)*dr[i]*dr[j]
                                       + rnl.Y*hess(i,j)
                                       + drnloverr*(dr[i]*gr_ang[j] + dr[j]*gr_ang[i]);
          grad_grad_Phi[offset](j,i) = grad_grad_Phi[offset](i,j);
        }
      }
      //d2psi[offset] = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

  inline void
  evaluateForWalkerMove(int c, int iat, int offset, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_Phi, GGGVector_t& grad_grad_grad_Phi)
  {
    /*********************************************************************
          int la=1;
          PosType  dr0=1,dr1=1,dr2=1;
          RealType r0,r1,r2;
          RealType dh=0.00001;
          dr0[0] = 0.5;dr1[0] = dr0[0];dr2[0] = dr0[0];
          dr0[1] = 1.0;dr1[1] = dr0[1];dr2[1] = dr0[1];
          dr0[2] = 1.5;dr1[2] = dr0[2];dr2[2] = dr0[2];
          dr1[la] = dr0[la]+dh;dr2[la] = dr0[la]-dh;
          r0 = std::sqrt(dr0[0]*dr0[0] + dr0[1]*dr0[1] + dr0[2]*dr0[2]);
          r1 = std::sqrt(dr1[0]*dr1[0] + dr1[1]*dr1[1] + dr1[2]*dr1[2]);
          r2 = std::sqrt(dr2[0]*dr2[0] + dr2[1]*dr2[1] + dr2[2]*dr2[2]);
          int nb = NL.size();
          ValueVector_t psi0,psi1,psi2;
          GradVector_t dpsi0,dpsi1,dpsi2;
          HessVector_t hpsi0,hpsi1,hpsi2;
          GGGVector_t gpsi0,gpsi1,gpsi2;

           psi0.resize(nb);
           psi1.resize(nb);
           psi2.resize(nb);
           dpsi0.resize(nb);
           dpsi1.resize(nb);
           dpsi2.resize(nb);
           hpsi0.resize(nb);
           hpsi1.resize(nb);
           hpsi2.resize(nb);
           gpsi0.resize(nb);
           gpsi1.resize(nb);
           gpsi2.resize(nb);

           dummyEval(r0,dr0,psi0,dpsi0,hpsi0,gpsi0);
           dummyEval(r1,dr1,psi1,dpsi1,hpsi1,gpsi1);
           dummyEval(r2,dr2,psi2,dpsi2,hpsi2,gpsi2);

           for(int i=0; i<nb; i++) {
              std::cout <<"i: " <<i << std::endl
                  <<"dPhi_x: " <<dpsi0[i][la]-(psi1[i]-psi2[i])/(2*dh) << std::endl
                  <<"hPhi_0: " <<hpsi0[i](la,0)-(dpsi1[i][0]-dpsi2[i][0])/(2*dh) << std::endl
                  <<"hPhi_1: " <<hpsi0[i](la,1)-(dpsi1[i][1]-dpsi2[i][1])/(2*dh) << std::endl
                  <<"hPhi_2: " <<hpsi0[i](la,2)-(dpsi1[i][2]-dpsi2[i][2])/(2*dh) << std::endl
                  <<"gPhi_00: " <<gpsi0[i][la](0,0)-(hpsi1[i](0,0)-hpsi2[i](0,0))/(2*dh)  <<"  " <<gpsi0[i][la](0,0)  << std::endl
                  <<"gPhi_11: " <<gpsi0[i][la](1,1)-(hpsi1[i](1,1)-hpsi2[i](1,1))/(2*dh) <<"  " <<gpsi0[i][la](1,1) << std::endl
                  <<"gPhi_22: " <<gpsi0[i][la](2,2)-(hpsi1[i](2,2)-hpsi2[i](2,2))/(2*dh) <<"  " <<gpsi0[i][la](2,2) << std::endl
                  <<"gPhi_01: " <<gpsi0[i][la](0,1)-(hpsi1[i](0,1)-hpsi2[i](0,1))/(2*dh) <<"  " <<gpsi0[i][la](0,1) << std::endl
                  <<"gPhi_02: " <<gpsi0[i][la](0,2)-(hpsi1[i](0,2)-hpsi2[i](0,2))/(2*dh) <<"  " <<gpsi0[i][la](0,2) << std::endl
                  <<"gPhi_12: " <<gpsi0[i][la](1,2)-(hpsi1[i](1,2)-hpsi2[i](1,2))/(2*dh) <<"  " <<gpsi0[i][la](1,2) << std::endl;

           }

      //cout<<"psi: " <<psi0[0] << std::endl
      //    <<"dpsi: " <<dpsi0[0] << std::endl
      //    <<"hpsi: " <<hpsi0[0] << std::endl
      //    <<"gpsi: " <<gpsi0[0] << std::endl;

          APP_ABORT("Aborting after testing. \n");
    *********************************************************************/
    int nn = myTable->M[c]+iat;
    RealType r(myTable->r(nn));
    RealType rinv(myTable->rinv(nn));
    PosType  dr(myTable->dr(nn));
    PosType  drr(dr*rinv);
    GGGType ggg1,ggg2;
    for(int i=0; i<3; i++)
    {
      ggg1[i] = drr[i]*outerProduct(drr,drr);
      ggg2[i]=0.0;
      for(int j=0; j<3; j++)
      {
        for(int k=0; k<3; k++)
        {
          if(i==j)
            (ggg2[i])(j,k) += rinv*drr[k];
          if(i==k)
            (ggg2[i])(j,k) += rinv*drr[j];
          if(k==j)
            (ggg2[i])(j,k) += rinv*drr[i];
        }
      }
    }
    if(useCartesian)
    {
      XYZ.evaluateWithThirdDeriv(dr);
    }
    else
    {
      // FIX FIX FIX: not implemented for spherical basis
      //Ylm.evaluateWithThirdDeriv(dr);
      Ylm.evaluateWithHessian(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<Tensor<RealType,3> >& hessYlm = useCartesian?XYZ.hessXYZ:Ylm.hessYlm;
    std::vector<GGGType>& gggYlm = useCartesian?XYZ.gggXYZ:Ylm.gggYlm;
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit)->evaluateWithThirdDeriv(r,rinv);
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      HessType hess(hessYlm[lm]);
      GGGType Slm_ijk(gggYlm[lm]);
      HessType Rhs;
      GGGType& ggg=grad_grad_grad_Phi[offset];
      psi[offset]  = ang*rnl.Y;
      dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
      RealType temp1=rnl.d2Y-drnloverr;
      RealType temp2=rnl.d3Y-3.0*rinv*temp1;
      // hessian of radial piece
      for(int i=0; i<3; i++)
      {
        Rhs(i,i) = temp1*drr[i]*drr[i] + rnl.dY*rinv;
        for(int j=i+1; j<3; j++)
        {
          Rhs(i,j) = temp1*drr[i]*drr[j];
          Rhs(j,i) = Rhs(i,j);
        }
      }
      grad_grad_Phi[offset] = ang*Rhs + outerProduct(gr_rad,gr_ang)
                              + outerProduct(gr_ang,gr_rad) + rnl.Y*hess;
      for(int i=0; i<3; i++)
      {
        ggg[i] = ang*(ggg1[i]*temp2 + ggg2[i]*temp1)
                 + rnl.Y*Slm_ijk[i];
      }
      // I don't know how to make this part compact, so sloppy for now
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
          for(int k=0; k<3; k++)
          {
            ggg[i](j,k) += Rhs(j,k)*gr_ang[i] + Rhs(i,k)*gr_ang[j] + Rhs(i,j)*gr_ang[k] + gr_rad[k]*hess(i,j) + gr_rad[j]*hess(i,k) + gr_rad[i]*hess(j,k);
          }
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

  inline void
  evaluateThirdDerivOnly(int c, int iat, int offset, GGGVector_t& grad_grad_grad_Phi)
  {
    int nn = myTable->M[c]+iat;
    RealType r(myTable->r(nn));
    RealType rinv(myTable->rinv(nn));
    PosType  dr(myTable->dr(nn));
    PosType  drr(dr*rinv);
    GGGType ggg1,ggg2;
    for(int i=0; i<3; i++)
    {
      ggg1[i] = drr[i]*outerProduct(drr,drr);
      ggg2[i]=0.0;
      for(int j=0; j<3; j++)
      {
        for(int k=0; k<3; k++)
        {
          if(i==j)
            (ggg2[i])(j,k) += rinv*drr[k];
          if(i==k)
            (ggg2[i])(j,k) += rinv*drr[j];
          if(k==j)
            (ggg2[i])(j,k) += rinv*drr[i];
        }
      }
    }
    if(useCartesian)
    {
      XYZ.evaluateThirdDerivOnly(dr);
    }
    else
    {
      // not implemented, so only works for lmax <= 2
      // where ggg is zero anyway
      if(Ylm.lmax() > 2)
      {
        APP_ABORT("Need to implement Third derivatives in Spherical basis. \n");
      }
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<Tensor<RealType,3> >& hessYlm = useCartesian?XYZ.hessXYZ:Ylm.hessYlm;
    std::vector<GGGType>& gggYlm = useCartesian?XYZ.gggXYZ:Ylm.gggYlm;
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit)->evaluateWithThirdDeriv(r,rinv);
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      HessType hess(hessYlm[lm]);
      GGGType Slm_ijk(gggYlm[lm]);
      HessType Rhs;
      GGGType& ggg=grad_grad_grad_Phi[offset];
      RealType temp1=rnl.d2Y-drnloverr;
      RealType temp2=rnl.d3Y-3.0*rinv*temp1;
      // hessian of radial piece
      for(int i=0; i<3; i++)
      {
        Rhs(i,i) = temp1*drr[i]*drr[i] + rnl.dY*rinv;
        for(int j=i+1; j<3; j++)
        {
          Rhs(i,j) = temp1*drr[i]*drr[j];
          Rhs(j,i) = Rhs(i,j);
        }
      }
      for(int i=0; i<3; i++)
      {
        ggg[i] = ang*(ggg1[i]*temp2 + ggg2[i]*temp1)
                 + rnl.Y*Slm_ijk[i];
      }
      // I don't know how to make this part compact, so sloppy for now
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
          for(int k=0; k<3; k++)
          {
            ggg[i](j,k) += Rhs(j,k)*gr_ang[i] + Rhs(i,k)*gr_ang[j] + Rhs(i,j)*gr_ang[k] + gr_rad[k]*hess(i,j) + gr_rad[j]*hess(i,k) + gr_rad[i]*hess(j,k);
          }
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

  void evaluateValues(const DistanceTableData* dt, int c, int offset, ValueMatrix_t& phiM)
  {
    int nn = dt->M[c];
    for(int iat=0; iat<dt->targets(); iat++, nn++)
    {
      RealType r(dt->r(nn));
      RealType rinv(dt->rinv(nn));
      for(int ri=0;ri<Rnl.size(); ++ri) Rnl[ri]->evaluate(r,rinv);
      ValueType* restrict y = phiM[iat]+offset;
      if(useCartesian)
      {
        XYZ.evaluate(dt->dr(nn));
        for(int nl=0; nl<NL.size(); ++nl) (*y++)=XYZ.XYZ[LM[nl]]*Rnl[NL[nl]]->Y;
      }
      else
      {
        Ylm.evaluate(dt->dr(nn));
        for(int nl=0; nl<NL.size(); ++nl) (*y++)=Ylm.Ylm[LM[nl]]*Rnl[NL[nl]]->Y;
      }
    }
  }

  void
  dummyEval(RealType r, PosType  dr, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_Phi, GGGVector_t& grad_grad_grad_Phi)
  {
    int offset=0;
    RealType rinv(1.0/r);
    PosType  drr(dr*rinv);
    GGGType ggg1,ggg2;
    for(int i=0; i<3; i++)
    {
      ggg1[i] = drr[i]*outerProduct(drr,drr);
      ggg2[i]=0.0;
      for(int j=0; j<3; j++)
      {
        for(int k=0; k<3; k++)
        {
          if(i==j)
            (ggg2[i])(j,k) += rinv*drr[k];
          if(i==k)
            (ggg2[i])(j,k) += rinv*drr[j];
          if(k==j)
            (ggg2[i])(j,k) += rinv*drr[i];
        }
      }
    }
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit)->evaluateWithThirdDeriv(r,rinv);
      ++rit;
    }
    if(useCartesian)
    {
      XYZ.evaluateWithThirdDeriv(dr);
    }
    else
    {
      // FIX FIX FIX: not implemented for spherical basis
      //Ylm.evaluateWithThirdDeriv(dr);
      Ylm.evaluateWithHessian(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<Tensor<RealType,3> >& hessYlm = useCartesian?XYZ.hessXYZ:Ylm.hessYlm;
    std::vector<GGGType>& gggYlm = useCartesian?XYZ.gggXYZ:Ylm.gggYlm;
    std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
    while(nlit != nlit_end)
      //for(int ib=0; ib<NL.size(); ib++, offset++) {
    {
      int nl(*nlit);//NL[ib];
      int lm(*lmit);//LM[ib];
      const ROT& rnl(*Rnl[nl]);
      RealType drnloverr(rinv*rnl.dY);
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      HessType hess(hessYlm[lm]);
      GGGType Slm_ijk(gggYlm[lm]);
      HessType Rhs;
      GGGType& ggg=grad_grad_grad_Phi[offset];
      psi[offset]  = ang*rnl.Y;
      dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
      RealType temp1=rnl.d2Y-drnloverr;
      RealType temp2=rnl.d3Y-3.0*rinv*temp1;
      // hessian of radial piece
      for(int i=0; i<3; i++)
      {
        Rhs(i,i) = temp1*drr[i]*drr[i] + rnl.dY*rinv;
        for(int j=i+1; j<3; j++)
        {
          Rhs(i,j) = temp1*drr[i]*drr[j];
          Rhs(j,i) = Rhs(i,j);
        }
      }
      grad_grad_Phi[offset] = ang*Rhs + outerProduct(gr_rad,gr_ang)
                              + outerProduct(gr_ang,gr_rad) + rnl.Y*hess;
      // assuming third derivatives of angular piece are zero
      for(int i=0; i<3; i++)
      {
        ggg[i] = ang*(ggg1[i]*temp2 + ggg2[i]*temp1)
                 + rnl.Y*Slm_ijk[i];
      }
      // I don't know how to make this part compact, so sloppy for now
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
          for(int k=0; k<3; k++)
          {
            ggg[i](j,k) += Rhs(j,k)*gr_ang[i] + Rhs(i,k)*gr_ang[j] + Rhs(i,j)*gr_ang[k] + gr_rad[k]*hess(i,j) + gr_rad[j]*hess(i,k) + gr_rad[i]*hess(j,k);
          }
      ++nlit;
      ++lmit;
      ++offset;
    }
  }


  inline void
  evaluateForWalkerMove(int source, int first, int nptcl, int offset, ValueMatrix_t& y, GradMatrix_t& dy, ValueMatrix_t& d2y)
  {
    int nn = myTable->M[source]+first;//first pair of the particle subset
    for(int i=0, iat=first; i<nptcl; i++, iat++, nn++)
    {
      RealType r(myTable->r(nn));
      RealType rinv(myTable->rinv(nn));
      PosType dr(myTable->dr(nn));
      if(useCartesian)
      {
        XYZ.evaluateAll(dr);
      }
      else
      {
        Ylm.evaluateAll(dr);
      }
      std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
      std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
      std::vector<RealType>& laplYlm = useCartesian?XYZ.laplXYZ:Ylm.laplYlm;
      typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
      while(rit != rit_end)
      {
        (*rit)->evaluateAll(r,rinv);
        ++rit;
      }
      int bindex(offset);
      std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
      ValueType* restrict yptr = y[iat]+offset;
      GradType* restrict dyptr = dy[iat]+offset;
      ValueType* restrict d2yptr = d2y[iat]+offset;
      while(nlit != nlit_end)
        //for(int ib=0; ib<NL.size(); ib++, offset++) {
      {
        int nl(*nlit);//NL[ib];
        int lm(*lmit);//LM[ib];
        const ROT& rnl(*Rnl[nl]);
        RealType drnloverr = rinv*rnl.dY;
        ValueType ang = valueYlm[lm];
        PosType gr_rad(drnloverr*dr);
        PosType gr_ang(gradYlm[lm]);
        *yptr++ = ang*rnl.Y;
        *dyptr++ = ang*gr_rad+rnl.Y*gr_ang;
        *d2yptr++= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang)  + rnl.Y*laplYlm[lm];
        //y(iat,bindex)= ang*rnl.Y;
        //dy(iat,bindex) = ang*gr_rad+rnl.Y*gr_ang;
        //d2y(iat,bindex)= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
        ++nlit;
        ++lmit;
        ++bindex;
      }
    }
  }


  ///evaluate the value, gradient and laplacian of basis functions for the iath-particle
  inline void
  evaluateForWalkerMove(int c, int iat, int offset, Matrix<ValueType>& temp)
  {
    //RealType r(myTable->Temp[c].r1);
    //RealType rinv(myTable->Temp[c].rinv1);
    //PosType  dr(myTable->Temp[c].dr1);
    //
    int nn = myTable->M[c]+iat;
    RealType r(myTable->r(nn));
    RealType rinv(myTable->rinv(nn));
    PosType  dr(myTable->dr(nn));
    if(useCartesian)
    {
      XYZ.evaluateAll(dr);
    }
    else
    {
      Ylm.evaluateAll(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<RealType>& laplYlm = useCartesian?XYZ.laplXYZ:Ylm.laplYlm;
    typename std::vector<ROT*>::iterator rit(Rnl.begin()), rit_end(Rnl.end());
    while(rit != rit_end)
    {
      (*rit++)->evaluateAll(r,rinv);
    }
    std::vector<int>::iterator nlit(NL.begin()),nlit_end(NL.end()),lmit(LM.begin());
    ValueType* restrict tptr=temp[offset];
    while(nlit != nlit_end)
      //for(int ib=0; ib<NL.size(); ib++, offset++) {
    {
      int nl(*nlit++);//NL[ib];
      int lm(*lmit++);//LM[ib];
      const ROT& rnl(*Rnl[nl]);
      RealType drnloverr(rinv*rnl.dY);
      RealType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      //PosType g(ang*gr_rad+rnl.Y*gr_ang);
      *tptr++ = ang*rnl.Y;
      *tptr++ = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang) + rnl.Y*laplYlm[lm];
      *tptr++ = ang*gr_rad[0]+rnl.Y*gr_ang[0];
      *tptr++ = ang*gr_rad[1]+rnl.Y*gr_ang[1];
      *tptr++ = ang*gr_rad[2]+rnl.Y*gr_ang[2];
      //PAOps<RealType,DIM>::copy(ang*rnl.Y, ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang), ang*gr_rad+rnl.Y*gr_ang, temp[offset]);
      //++nlit; ++lmit;
    }
  }

  inline void
  evaluateForPtclMove(int source, int iat,  int offset, ValueVector_t& y)
  {
    RealType r(myTable->Temp[source].r1);
    //RealType rinv(myTable->Temp[source].rinv1);
    RealType rinv(1/r);
    PosType  dr(myTable->Temp[source].dr1);
    if(useCartesian)
    {
      XYZ.evaluate(dr);
    }
    else
    {
      Ylm.evaluate(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
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
      //y[offset++]=Ylm.getYlm(*lmit++)*Rnl[*nlit++]->Y;
      y[offset++]=valueYlm[*lmit++]*Rnl[*nlit++]->Y;
    }
  }

  inline void
  evaluateAllForPtclMove(int source, int iat,  int offset, ValueVector_t& y,
                         GradVector_t& dy, ValueVector_t& d2y)
  {
    RealType r(myTable->Temp[source].r1);
    RealType rinv(myTable->Temp[source].rinv1);
    PosType  dr(myTable->Temp[source].dr1);
    if(useCartesian)
    {
      XYZ.evaluateAll(dr);
    }
    else
    {
      Ylm.evaluateAll(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<RealType>& laplYlm = useCartesian?XYZ.laplXYZ:Ylm.laplYlm;
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      y[offset]= ang*rnl.Y;
      dy[offset] = ang*gr_rad+rnl.Y*gr_ang;
      d2y[offset]= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang) + rnl.Y*laplYlm[lm];
      ++nlit;
      ++lmit;
      ++offset;
    }
  }


  inline void
  evaluate(RealType r, RealType rinv, const PosType& dr, int offset, ValueVector_t& psi)
  {
    if(useCartesian)
    {
      XYZ.evaluate(dr);
    }
    else
    {
      Ylm.evaluate(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
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
      //psi[offset++]=Ylm.getYlm(*lmit++)*Rnl[*nlit++]->Y;
      psi[offset++]=valueYlm[*lmit++]*Rnl[*nlit++]->Y;
    }
  }

  inline void
  evaluate(RealType r, RealType rinv, const PosType& dr, int offset, ValueVector_t& y,
           GradVector_t& dy, ValueVector_t& d2y)
  {
    if(useCartesian)
    {
      XYZ.evaluateAll(dr);
    }
    else
    {
      Ylm.evaluateAll(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<RealType>& laplYlm = useCartesian?XYZ.laplXYZ:Ylm.laplYlm;
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      y[offset]= ang*rnl.Y;
      dy[offset] = ang*gr_rad+rnl.Y*gr_ang;
      d2y[offset]= ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang) + rnl.Y*laplYlm[lm];
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

  inline void
  evaluateAllForPtclMove(int source, int iat, int offset, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_Phi)
  {
    RealType r(myTable->Temp[source].r1);
    RealType rinv(myTable->Temp[source].rinv1);
    PosType  dr(myTable->Temp[source].dr1);
    if(useCartesian)
    {
      XYZ.evaluateWithHessian(dr);
    }
    else
    {
      Ylm.evaluateWithHessian(dr);
    }
    std::vector<RealType>& valueYlm = useCartesian?XYZ.XYZ:Ylm.Ylm;
    std::vector<PosType>& gradYlm = useCartesian?XYZ.gradXYZ:Ylm.gradYlm;
    std::vector<Tensor<RealType,3> >& hessYlm = useCartesian?XYZ.hessXYZ:Ylm.hessYlm;
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
      ValueType ang(valueYlm[lm]);
      PosType gr_rad(drnloverr*dr);
      PosType gr_ang(gradYlm[lm]);
      HessType hess(hessYlm[lm]);
      psi[offset]  = ang*rnl.Y;
      dpsi[offset] = ang*gr_rad+rnl.Y*gr_ang;
// sloppy for now
      RealType temp1=rnl.d2Y*ang*rinv*rinv;
      RealType temp2=drnloverr*ang*rinv*rinv;
      for(int i=0; i<3; i++)
      {
        grad_grad_Phi[offset](i,i) = (temp1-temp2)*dr[i]*dr[i]
                                     + drnloverr*ang + rnl.Y*hess(i,i)
                                     + 2*drnloverr*dr[i]*gr_ang[i];
        for(int j=i+1; j<3; j++)
        {
          grad_grad_Phi[offset](i,j) = (temp1-temp2)*dr[i]*dr[j]
                                       + rnl.Y*hess(i,j)
                                       + drnloverr*(dr[i]*gr_ang[j] + dr[j]*gr_ang[i]);
          grad_grad_Phi[offset](j,i) = grad_grad_Phi[offset](i,j);
        }
      }
      //d2psi[offset] = ang*(2.0*drnloverr+rnl.d2Y) + 2.0*dot(gr_rad,gr_ang);
      ++nlit;
      ++lmit;
      ++offset;
    }
  }

};

}
#endif

