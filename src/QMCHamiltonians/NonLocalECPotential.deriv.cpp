//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
//#include "Particle/VirtualParticleSet.h"
#include <simd/simd.hpp>
#include <Utilities/Timer.h>

namespace qmcplusplus
{
  NonLocalECPotential::Return_t
    NonLocalECPotential::evaluateValueAndDerivatives(ParticleSet& P,
        const opt_variables_type& optvars,
        const std::vector<RealType>& dlogpsi,
        std::vector<RealType>& dhpsioverpsi)
    {
      Value=0.0;
      for(int iat=0; iat<NumIons; iat++)
        if(PP[iat])
        {
          PP[iat]->randomize_grid(*(P.Sphere[iat]),UpdateMode[PRIMARY]);
          Value+=PP[iat]->evaluateValueAndDerivatives(P,iat,Psi,optvars,dlogpsi,dhpsioverpsi);
        }
      return Value;

      //int Nvars=optvars.size();
      //vector<RealType> dh_ana(Nvars,0.0);
      //Return_t e=0.0;
      //for(int iat=0; iat<NumIons; iat++)
      //  if(PP[iat])
      //  {
      //    PP[iat]->randomize_grid(*(P.Sphere[iat]),false);
      //    e += PP[iat]->evaluate(P,iat,Psi);
      //  }

      //vector<RealType> dh_diff(Nvars,0.0);
      //opt_variables_type wfVars(optvars),wfvar_prime(optvars);
      //RealType FiniteDiff = 1e-6;
      //QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
      //for (int i=0; i<Nvars ; i++)
      //{
      //  wfvar_prime[i] = wfVars[i]+ FiniteDiff;
      //  //     Psi.checkOutVariables(wfvar_prime);
      //  Psi.resetParameters(wfvar_prime);
      //  Psi.reset();
      //  RealType logpsiPlus = Psi.evaluateLog(P);
      //  RealType elocPlus=0.0;
      //  for(int iat=0; iat<NumIons; iat++) 
      //    if(PP[iat])                       
      //    {
      //      PP[iat]->randomize_grid(*(P.Sphere[iat]),false);
      //      elocPlus += PP[iat]->evaluate(P,iat,Psi);
      //    }

      //  wfvar_prime[i] = wfVars[i]- FiniteDiff;
      //  //     Psi.checkOutVariables(wfvar_prime);
      //  Psi.resetParameters(wfvar_prime);
      //  Psi.reset();
      //  RealType logpsiMinus = Psi.evaluateLog(P);
      //  RealType elocMinus=0.0;
      //  for(int iat=0; iat<NumIons; iat++) 
      //    if(PP[iat])                       
      //    {
      //      PP[iat]->randomize_grid(*(P.Sphere[iat]),false);
      //      elocMinus += PP[iat]->evaluate(P,iat,Psi);
      //    }
      //  dh_diff[i]= (elocPlus-elocMinus)*dh;
      //}
      //Psi.resetParameters(optvars);
      //Psi.reset();

      //app_log() << "DERIV " << std::endl;
      //for(int v=0;v<Nvars; ++v)
      //  app_log() << v << " " << dh_ana[v] << " " << dh_diff[v] << std::endl;

      //APP_ABORT("DONE");
      return Value;
    }

  /** evaluate the non-local potential of the iat-th ionic center
   * @param W electron configuration
   * @param iat ionic index
   * @param psi trial wavefunction
   * @param optvars optimizables 
   * @param dlogpsi derivatives of the wavefunction at W.R 
   * @param hdpsioverpsi derivatives of Vpp 
   * @param return the non-local component
   *
   * This is a temporary solution which uses TrialWaveFunction::evaluateDerivatives
   * assuming that the distance tables are fully updated for each ratio computation.
   */
  NonLocalECPComponent::RealType
    NonLocalECPComponent::evaluateValueAndDerivatives(ParticleSet& W,
        int iat, TrialWaveFunction& psi,
        const opt_variables_type& optvars,
        const std::vector<RealType>& dlogpsi,
        std::vector<RealType>& dhpsioverpsi)
    {
      Matrix<RealType> dratio(optvars.num_active_vars,nknot);
      std::vector<RealType> dlogpsi_t(dlogpsi.size(),0.0);
      std::vector<RealType> dhlogpsi_t(dlogpsi.size(),0.0);

      DistanceTableData* myTable = W.DistTables[myTableIndex];
      RealType esum=0.0;
      RealType pairpot;
      ParticleSet::ParticlePos_t deltarV(nknot);
      for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++)
      {
        register RealType r(myTable->r(nn));
        if(r>Rmax)
          continue;

        register RealType rinv(myTable->rinv(nn));
        register PosType  dr(myTable->dr(nn));

        //displacements wrt W.R[iel]
        for (int j=0; j < nknot ; j++) deltarV[j]=r*rrotsgrid_m[j]-dr;

        for (int j=0; j < nknot ; j++)
        {
          PosType pos_now=W.R[iel];
          W.makeMoveAndCheck(iel,deltarV[j]);
#if defined(QMC_COMPLEX)
          psiratio[j]=psi.ratio(W,iel)*std::cos(psi.getPhaseDiff());
#else
          psiratio[j]=psi.ratio(W,iel);
#endif
          psi.resetPhaseDiff();

          //use existing methods
          W.acceptMove(iel);

          std::fill(dlogpsi_t.begin(),dlogpsi_t.end(),0.0);
          psi.evaluateDerivatives(W, optvars, dlogpsi_t,dhlogpsi_t);
          for(int v=0; v<dlogpsi_t.size(); ++v) dratio(v,j)=dlogpsi_t[v];

          PosType md=-1.0*deltarV[j];
          W.makeMoveAndCheck(iel,md);
          W.acceptMove(iel);
        }

        for(int j=0; j<nknot; ++j) psiratio[j]*=sgridweight_m[j];

        for(int ip=0; ip< nchannel; ip++)
          vrad[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];
        // Compute spherical harmonics on grid
        for (int j=0, jl=0; j<nknot ; j++)
        {
          RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
          // Forming the Legendre polynomials
          lpol[0]=1.0;
          RealType lpolprev=0.0;
          for (int l=0 ; l< lmax ; l++)
          {
            lpol[l+1]=Lfactor1[l]*zz*lpol[l]-l*lpolprev;
            lpol[l+1]*=Lfactor2[l];
            lpolprev=lpol[l];
          }
          for(int l=0; l <nchannel; l++,jl++)
            Amat[jl]=lpol[ angpp_m[l] ];
        }
        if(nchannel==1)
        {
          pairpot = vrad[0]*BLAS::dot(nknot, &Amat[0],&psiratio[0]);
          for(int v=0; v<dhpsioverpsi.size(); ++v)
          {
            for(int j=0; j<nknot; ++j)
              dratio(v,j)=psiratio[j]*(dratio(v,j)-dlogpsi[v]);
            dhpsioverpsi[v]+=vrad[0]*BLAS::dot(nknot, &Amat[0],dratio[v]);
          }
        }
        else
        {
          BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
          pairpot = BLAS::dot(nchannel, &vrad[0], &wvec[0]);
          for(int v=0; v<dhpsioverpsi.size(); ++v)
          {
            for(int j=0; j<nknot; ++j)
              dratio(v,j)=psiratio[j]*(dratio(v,j)-dlogpsi[v]);
            BLAS::gemv(nknot, nchannel, &Amat[0], dratio[v], &wvec[0]);
            dhpsioverpsi[v]+=BLAS::dot(nchannel, &vrad[0], &wvec[0]);
          }
        }


        esum += pairpot;
      }   /* end loop over electron */

      return esum;
    }

}
