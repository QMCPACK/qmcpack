//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BACKFLOW_ELEC_ION_H
#define QMCPLUSPLUS_BACKFLOW_ELEC_ION_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include <cmath>
#include <vector>

namespace qmcplusplus
{
template<class FT>
class Backflow_eI : public BackflowFunctionBase
{
private:
  /// distance table index
  const int myTableIndex_;

public:
  std::vector<FT*> RadFun;
  std::vector<FT*> uniqueRadFun;
  std::vector<int> offsetPrms;

  Backflow_eI(ParticleSet& ions, ParticleSet& els)
    : BackflowFunctionBase(ions, els), myTableIndex_(els.addTable(ions, DT_SOA_PREFERRED))
  {
    resize(NumTargets, NumCenters);
  }

  //  build RadFun manually from builder class
  Backflow_eI(ParticleSet& ions, ParticleSet& els, FT* RF)
    : BackflowFunctionBase(ions, els), myTableIndex_(els.addTable(ions, DT_SOA_PREFERRED))
  {
    // same radial function for all centers by default
    uniqueRadFun.push_back(RF);
    for (int i = 0; i < NumCenters; i++)
      RadFun.push_back(RF);
  }

  ~Backflow_eI(){};

  BackflowFunctionBase* makeClone(ParticleSet& tqp)
  {
    Backflow_eI<FT>* clone = new Backflow_eI<FT>(CenterSys, tqp);
    clone->resize(NumTargets, NumCenters);
    clone->indexOfFirstParam = indexOfFirstParam;
    clone->offsetPrms        = offsetPrms;
    clone->numParams         = numParams;
    clone->derivs            = derivs;
    clone->uniqueRadFun.resize(uniqueRadFun.size());
    clone->RadFun.resize(RadFun.size());
    for (int i = 0; i < uniqueRadFun.size(); i++)
      clone->uniqueRadFun[i] = new FT(*(uniqueRadFun[i]));
    for (int i = 0; i < RadFun.size(); i++)
    {
      bool done = false;
      for (int k = 0; k < uniqueRadFun.size(); k++)
        if (RadFun[i] == uniqueRadFun[k])
        {
          done             = true;
          clone->RadFun[i] = clone->uniqueRadFun[k];
          break;
        }
      if (!done)
      {
        APP_ABORT("Error cloning Backflow_eI object. \n");
      }
    }
    return clone;
  }

  void reportStatus(std::ostream& os)
  {
    std::cerr << RadFun.size() << std::endl;
    std::cerr << isOptimizable() << std::endl;
    std::cerr << myTableIndex_ << std::endl;
    for (int i = 0; i < offsetPrms.size(); i++)
      std::cerr << offsetPrms[i] << std::endl;
    for (int i = 0; i < uniqueRadFun.size(); i++)
      uniqueRadFun[i]->reportStatus(os);
  }

  void resetParameters(const opt_variables_type& active)
  {
    for (int i = 0; i < uniqueRadFun.size(); i++)
      uniqueRadFun[i]->resetParameters(active);
  }

  void checkInVariables(opt_variables_type& active)
  {
    for (int i = 0; i < uniqueRadFun.size(); i++)
      uniqueRadFun[i]->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    for (int i = 0; i < uniqueRadFun.size(); i++)
      uniqueRadFun[i]->checkOutVariables(active);
  }

  inline bool isOptimizable()
  {
    for (int i = 0; i < uniqueRadFun.size(); i++)
      if (uniqueRadFun[i]->isOptimizable())
        return true;
    return false;
  }

  inline int indexOffset() { return RadFun[0]->myVars.where(0); }

  // only elements of qp coordinates that changed should be copied
  inline void acceptMove(int iat, int UpdateMode)
  {
    int num;
    switch (UpdateMode)
    {
    case ORB_PBYP_RATIO:
      num = UIJ.cols();
      for (int i = 0; i < num; i++)
        UIJ(iat, i) = UIJ_temp[i];
      break;
    case ORB_PBYP_PARTIAL:
      num = UIJ.cols();
      for (int i = 0; i < num; i++)
        UIJ(iat, i) = UIJ_temp[i];
      num = AIJ.cols();
      for (int i = 0; i < num; i++)
        AIJ(iat, i) = AIJ_temp[i];
      break;
    case ORB_PBYP_ALL:
      num = UIJ.cols();
      for (int i = 0; i < num; i++)
        UIJ(iat, i) = UIJ_temp[i];
      num = AIJ.cols();
      for (int i = 0; i < num; i++)
        AIJ(iat, i) = AIJ_temp[i];
      num = BIJ.cols();
      for (int i = 0; i < num; i++)
        BIJ(iat, i) = BIJ_temp[i];
      break;
    default:
      num = UIJ.cols();
      for (int i = 0; i < num; i++)
        UIJ(iat, i) = UIJ_temp[i];
      num = AIJ.cols();
      for (int i = 0; i < num; i++)
        AIJ(iat, i) = AIJ_temp[i];
      num = BIJ.cols();
      for (int i = 0; i < num; i++)
        BIJ(iat, i) = BIJ_temp[i];
      break;
    }
    UIJ_temp = 0.0;
    AIJ_temp = 0.0;
    BIJ_temp = 0.0;
  }

  inline void restore(int iat, int UpdateType)
  {
    UIJ_temp = 0.0;
    AIJ_temp = 0.0;
    BIJ_temp = 0.0;
  }

  void registerData(WFBufferType& buf)
  {
    FirstOfU = &(UIJ(0, 0)[0]);
    LastOfU  = FirstOfU + OHMMS_DIM * NumTargets * NumCenters;
    FirstOfA = &(AIJ(0, 0)[0]);
    LastOfA  = FirstOfA + OHMMS_DIM * OHMMS_DIM * NumTargets * NumCenters;
    FirstOfB = &(BIJ(0, 0)[0]);
    LastOfB  = FirstOfB + OHMMS_DIM * NumTargets * NumCenters;
    buf.add(FirstOfU, LastOfU);
    buf.add(FirstOfA, LastOfA);
    buf.add(FirstOfB, LastOfB);
  }

  /** calculate quasi-particle coordinates only
   */
  inline void evaluate(const ParticleSet& P, ParticleSet& QP)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluate(P,QP) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[i]->evaluate(myTable.r(nn), du, d2u);
        QP.R[j] += (UIJ(j, i) = uij * myTable.dr(nn)); // dr(ij) = r_j-r_i
      }
    }
#endif
  }


  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradVector_t& Bmat, HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluate(P,QP,Bmat_vec,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[i]->evaluate(myTable.r(nn), du, d2u);
        //PosType u = uij*myTable.dr(nn);
        du *= myTable.rinv(nn);
        QP.R[j] += (UIJ(j, i) = uij * myTable.dr(nn));
        HessType& hess = AIJ(j, i);
        hess           = du * outerProduct(myTable.dr(nn), myTable.dr(nn));
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
        Amat(j, j) += hess;
        //u = (d2u+4.0*du)*myTable.dr(nn);
        Bmat[j] += (BIJ(j, i) = (d2u + 4.0 * du) * myTable.dr(nn));
      }
    }
#endif
  }


  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat_full, HessMatrix_t& Amat)
  {
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
#ifdef ENABLE_SOA
    //     APP_ABORT("Backflow_eI.h::evaluate(P,QP,Bmat_full,Amat) not implemented for SoA\n");
    for (int jel = 0; jel < P.getTotalNum(); jel++)
    {
      const auto& dist  = myTable.Distances[jel];
      const auto& displ = myTable.Displacements[jel];
      for (int iat = 0; iat < NumCenters; iat++)
      {
        if (dist[iat] > 0)
        {
          RealType uij = RadFun[iat]->evaluate(dist[iat], du, d2u);
          du /= dist[iat];
          QP.R[jel] += (UIJ(jel, iat) = -uij * displ[iat]);
          HessType& hess = AIJ(jel, iat);
          hess           = du * outerProduct(displ[iat], displ[iat]);
          hess[0] += uij;
          hess[4] += uij;
          hess[8] += uij;
          Amat(jel, jel) += hess;
          // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
          //u = (d2u+4.0*du)*myTable.dr(nn);
          Bmat_full(jel, jel) += (BIJ(jel, iat) = -(d2u + 4.0 * du) * displ[iat]);
        }
      }
    }

#else
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[i]->evaluate(myTable.r(nn), du, d2u);
        du *= myTable.rinv(nn);
        //PosType u = uij*myTable.dr(nn);
        QP.R[j] += (UIJ(j, i) = uij * myTable.dr(nn));
        HessType& hess = AIJ(j, i);
        hess           = du * outerProduct(myTable.dr(nn), myTable.dr(nn));
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
        Amat(j, j) += hess;
        // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
        //u = (d2u+4.0*du)*myTable.dr(nn);
        Bmat_full(j, j) += (BIJ(j, i) = (d2u + 4.0 * du) * myTable.dr(nn));
      }
    }
#endif
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP, const std::vector<int>& index)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluatePbyP(P,QP,index_vec) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI = myTable.sources();
    int iat  = index[0];
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
    }
#endif
  }


  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP)
  {
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
#ifdef ENABLE_SOA
    // APP_ABORT("Backflow_eI.h::evaluatePbyP(P,iat,QP) not implemented for SoA\n");
    int maxI = myTable.sources();
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp_r[j], du, d2u);
      PosType u    = (UIJ_temp[j] = -uij * myTable.Temp_dr[j]) - UIJ(iat, j);
      newQP[iat] += u;
    }
#else
    int maxI = myTable.sources();
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
    }
#endif
  }

  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos_t& newQP,
                           const std::vector<int>& index,
                           HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluatePbyP(P,QP,index_vec,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI = myTable.sources();
    int iat  = index[0];
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      HessType& hess = AIJ_temp[j];
      hess           = (du * myTable.Temp[j].rinv1) * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
      // should I expand this??? Is the compiler smart enough???
      Amat(iat, iat) += (hess - AIJ(iat, j));
    }
#endif
  }

  inline void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP, HessMatrix_t& Amat)
  {
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
#ifdef ENABLE_SOA
    //     APP_ABORT("Backflow_eI.h::evaluatePbyP(P,iat,QP,Amat) not implemented for SoA\n");
    int maxI = myTable.sources();
    for (int j = 0; j < maxI; j++)
    {
      if (myTable.Temp_r[j] > 0)
      {
        RealType uij = RadFun[j]->evaluate(myTable.Temp_r[j], du, d2u);
        PosType u    = (UIJ_temp[j] = -uij * myTable.Temp_dr[j]) - UIJ(iat, j);
        newQP[iat] += u;
        HessType& hess = AIJ_temp[j];
        hess           = (du / myTable.Temp_r[j]) * outerProduct(myTable.Temp_dr[j], myTable.Temp_dr[j]);
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
        // should I expand this??? Is the compiler smart enough???
        Amat(iat, iat) += (hess - AIJ(iat, j));
      }
    }
#else
    int maxI = myTable.sources();
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      HessType& hess = AIJ_temp[j];
      hess           = (du * myTable.Temp[j].rinv1) * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
      // should I expand this??? Is the compiler smart enough???
      Amat(iat, iat) += (hess - AIJ(iat, j));
    }
#endif
  }

  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos_t& newQP,
                           const std::vector<int>& index,
                           GradMatrix_t& Bmat_full,
                           HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluatePbyP(P,QP,index_vec,Bmat,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI = myTable.sources();
    int iat  = index[0];
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      du *= myTable.Temp[j].rinv1;
      HessType& hess = AIJ_temp[j];
      hess           = du * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
      Amat(iat, iat) += (hess - AIJ(iat, j));
      BIJ_temp[j] = (d2u + 4.0 * du) * myTable.Temp[j].dr1;
      Bmat_full(iat, iat) += (BIJ_temp[j] - BIJ(iat, j));
    }
#endif
  }

  inline void evaluatePbyP(const ParticleSet& P,
                           int iat,
                           ParticleSet::ParticlePos_t& newQP,
                           GradMatrix_t& Bmat_full,
                           HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluatePbyP(P,iat,QP,Bmat,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI = myTable.sources();
    for (int j = 0; j < maxI; j++)
    {
      RealType uij = RadFun[j]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      du *= myTable.Temp[j].rinv1;
      HessType& hess = AIJ_temp[j];
      hess           = du * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
      Amat(iat, iat) += (hess - AIJ(iat, j));
      BIJ_temp[j] = (d2u + 4.0 * du) * myTable.Temp[j].dr1;
      Bmat_full(iat, iat) += (BIJ_temp[j] - BIJ(iat, j));
    }
#endif
  }

  /** calculate only Bmat
   *  This is used in pbyp moves, in updateBuffer()
   */
  inline void evaluateBmatOnly(const ParticleSet& P, GradMatrix_t& Bmat_full)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_eI.h::evaluateBmatOnly(P,QP,Bmat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[i]->evaluate(myTable.r(nn), du, d2u);
        Bmat_full(j, j) += (BIJ(j, i) = (d2u + 4.0 * du * myTable.rinv(nn)) * myTable.dr(nn));
      }
    }
#endif
  }

  /** calculate quasi-particle coordinates, Bmat and Amat
   *  calculate derivatives wrt to variational parameters
   */
  inline void evaluateWithDerivatives(const ParticleSet& P,
                                      ParticleSet& QP,
                                      GradMatrix_t& Bmat_full,
                                      HessMatrix_t& Amat,
                                      GradMatrix_t& Cmat,
                                      GradMatrix_t& Ymat,
                                      HessArray_t& Xmat)
  {
#ifdef ENABLE_SOA
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int jel = 0; jel < P.getTotalNum(); jel++)
    {
      const auto& dist  = myTable.Distances[jel];
      const auto& displ = myTable.Displacements[jel];
      for (int iat = 0; iat < NumCenters; iat++)
      {
        if (dist[iat] > 0)
        {
          RealType uij = RadFun[iat]->evaluate(dist[iat], du, d2u);
          int NPrms    = RadFun[iat]->NumParams;
          std::vector<TinyVector<RealType, 3>> derivsju(NPrms);
          RadFun[iat]->evaluateDerivatives(dist[iat], derivsju);
          du /= dist[iat];
          QP.R[jel] += (UIJ(jel, iat) = -uij * displ[iat]);
          HessType op    = outerProduct(displ[iat], displ[iat]);
          HessType& hess = AIJ(jel, iat);
          hess           = du * op;
          hess[0] += uij;
          hess[4] += uij;
          hess[8] += uij;
          Amat(jel, jel) += hess;
          Bmat_full(jel, jel) += (BIJ(jel, iat) = -(d2u + 4.0 * du) * displ[iat]);
          for (int prm = 0, la = indexOfFirstParam + offsetPrms[iat]; prm < NPrms; prm++, la++)
          {
            Cmat(la, jel) -= displ[iat] * derivsju[prm][0];
            Xmat(la, jel, jel) += (derivsju[prm][1] / dist[iat]) * op;
            Xmat(la, jel, jel)[0] += derivsju[prm][0];
            Xmat(la, jel, jel)[4] += derivsju[prm][0];
            Xmat(la, jel, jel)[8] += derivsju[prm][0];
            Ymat(la, jel) -= (derivsju[prm][2] + 4.0 * derivsju[prm][1] / dist[iat]) * displ[iat];
          }
        }
      }
    }

#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[i]->evaluate(myTable.r(nn), du, d2u);
        //           std::fill(derivs.begin(),derivs.end(),0.0);
        int NPrms = RadFun[i]->NumParams;
        std::vector<TinyVector<RealType, 3>> derivsju(NPrms);
        RadFun[i]->evaluateDerivatives(myTable.r(nn), derivsju);
        du *= myTable.rinv(nn);
        //PosType u = uij*myTable.dr(nn);
        QP.R[j] += (UIJ(j, i) = uij * myTable.dr(nn));
        HessType op    = outerProduct(myTable.dr(nn), myTable.dr(nn));
        HessType& hess = AIJ(j, i);
        hess           = du * op;
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
        Amat(j, j) += hess;
        // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
        //u = (d2u+4.0*du)*myTable.dr(nn);
        Bmat_full(j, j) += (BIJ(j, i) = (d2u + 4.0 * du) * myTable.dr(nn));
        for (int prm = 0, la = indexOfFirstParam + offsetPrms[i]; prm < NPrms; prm++, la++)
        {
          Cmat(la, j) += myTable.dr(nn) * derivsju[prm][0];
          Xmat(la, j, j) += (derivsju[prm][1] * myTable.rinv(nn)) * op;
          Xmat(la, j, j)[0] += derivsju[prm][0];
          Xmat(la, j, j)[4] += derivsju[prm][0];
          Xmat(la, j, j)[8] += derivsju[prm][0];
          Ymat(la, j) += (derivsju[prm][2] + 4.0 * derivsju[prm][1] * myTable.rinv(nn)) * myTable.dr(nn);
        }
      }
    }
#endif
  }
};
} // namespace qmcplusplus

#endif
