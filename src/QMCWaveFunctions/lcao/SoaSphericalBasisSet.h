//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SphericalBasisSet.h
 * @brief A basis set of spherical symmetry associated with a center
 */
#ifndef QMCPLUSPLUS_SOA_SPHERICALORBITAL_BASISSET_H
#define QMCPLUSPLUS_SOA_SPHERICALORBITAL_BASISSET_H

namespace qmcplusplus
{

  template<typename ROT, typename SH>
    struct SoaSphericalBasisSet
    {
      typedef ROT                      RadialOrbital_t;
      typedef typename ROT::value_type value_type;

      ///size of the basis set
      size_t BasisSetSize;
      ///spherical harmonics
      SH Ylm;
      ///index of the corresponding real Spherical Harmonic with quantum numbers \f$ (l,m) \f$
      aligned_vector<int> LM;
      /**index of the corresponding radial orbital with quantum numbers \f$ (n,l) \f$ */
      aligned_vector<int> NL;
      ///container for the radial orbitals
      aligned_vector<ROT*> Rnl;
      ///container for the quantum-numbers
      std::vector<QuantumNumberType> RnlID;

      ///the constructor
      explicit SoaSphericalBasisSet(int lmax, bool addsignforM=false)
        :Ylm(lmax,addsignforM){}

      ~SoaSphericalBasisSet(){ } //cleanup

      SoaSphericalBasisSet<ROT,SH>* makeClone() const
      {
        SoaSphericalBasisSet<ROT,SH>* myclone=new SoaSphericalBasisSet<ROT,SH>(*this);
        for(int i=0; i<Rnl.size(); ++i)
          myclone->Rnl[i]=new ROT(Rnl[i]->makeClone());
        return myclone;
      }

      void checkInVariables(opt_variables_type& active)
      {
        for(size_t nl=0; nl<Rnl.size(); nl++)
          Rnl[nl]->checkInVariables(active);
      }

      void checkOutVariables(const opt_variables_type& active)
      {
        for(size_t nl=0; nl<Rnl.size(); nl++)
          Rnl[nl]->checkOutVariables(active);
      }

      void resetParameters(const opt_variables_type& active)
      {
        for(size_t nl=0; nl<Rnl.size(); nl++)
          Rnl[nl]->resetParameters(active);
      }

      /** return the number of basis functions
      */
      inline int getBasisSetSize() const
      {
        return BasisSetSize;
      }//=NL.size(); 


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

      ///set the current offset
      inline void setCenter(int c, int offset)
      { }


      template<typename VGL>
        inline void
        evaluateVGL(const T r, const PosType& dr, const size_t offset,  VGL& vgl)
        {
          CONSTEXPR T cone(1);
          CONSTEXPR T ctwo(2);
          const T x=dr[0], y=dr[1], z=dr[2];
          Ylm.evaluateVGL(x,y,z);

          const size_t nl_max=Rnl.size();
          T phi[nl_max]; 
          T dphi[nl_max];
          T d2phi[nl_max];

          for(size_t nl=0; nl<nl_max; ++nl)
            phi[nl]=Rnl[nl]->evaluate(r,dphi[nl],d2phi[nl]);

          //V,Gx,Gy,Gz,L
          T* restrict psi   =vgl[0]+offset; const T* restrict ylm_v=Ylm[0]; //value
          T* restrict dpsi_x=vgl[1]+offset; const T* restrict ylm_x=Ylm[1]; //gradX
          T* restrict dpsi_y=vgl[2]+offset; const T* restrict ylm_y=Ylm[2]; //gradY
          T* restrict dpsi_z=vgl[3]+offset; const T* restrict ylm_z=Ylm[3]; //gradZ
          T* restrict d2psi =vgl[4]+offset; const T* restrict ylm_l=Ylm[4]; //lap
          const T rinv=cone/r;
          const size_t ib_max=NL.size();
          for(size_t ib=0; ib<ib_max; ++ib)
          {
            const int nl(NL[ib]);
            const int lm(LM[ib]);
            const T drnloverr=rinv*dphi[nl];
            const T ang=ylm_v[lm];
            const T gr_x=drnloverr*x;
            const T gr_y=drnloverr*y;
            const T gr_z=drnloverr*z;
            const T ang_x=ylm_x[lm];
            const T ang_y=ylm_y[lm];
            const T ang_z=ylm_z[lm];
            const T vr=phi[nl];
            psi[ib]    = ang*vr;
            dpsi_x[ib] = ang*gr_x+vr*ang_x;
            dpsi_y[ib] = ang*gr_y+vr*ang_y;
            dpsi_z[ib] = ang*gr_z+vr*ang_z;
            d2psi[ib]  = ang*(ctwo*drnloverr+d2phi[nl]) + ctwo*(gr_x*ang_x+gr_y*ang_y+gr_z*ang_z)+vr*ylm_l[lm];
          }
        }

      inline void
        evaluateV(const T r, const PosType& dr, T* restrict psi)
        {
          CONSTEXPR T cone(1);
          CONSTEXPR T ctwo(2);
          const nl_max=Rnl.size();
          Ylm.evaluateV(dr[0],dr[1],dr[2]);
          T phi_r[nl_max];
          for(size_t nl=0; nl<nl_max; ++nl)
            phi_r[nl]=Rnl[nl]->evaluate(r);

          const T* restrict ylm_v=Ylm[0];
          for(size_t ib=0; ib<BasisSetSize; ++ib)
          {
            psi[ib]  = ylm_v[ LM[ib] ]*phi_r[ NL[ib] ];
          }
        }

    };

}
#endif
