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


#ifndef QMCPLUSPLUS_BACKFLOW_E_E_H
#define QMCPLUSPLUS_BACKFLOW_E_E_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "Message/Communicate.h"
#include <cmath>

namespace qmcplusplus
{
template<class FT>
class Backflow_ee : public BackflowFunctionBase
{
private:
  /// distance table index
  const int myTableIndex_;

public:
  //number of groups of the target particleset
  std::vector<FT*> RadFun;
  std::vector<FT*> uniqueRadFun;
  std::vector<int> offsetPrms;
  int NumGroups;
  Matrix<int> PairID;
  bool first;

  Backflow_ee(ParticleSet& ions, ParticleSet& els)
    : BackflowFunctionBase(ions, els), myTableIndex_(els.addTable(els, DT_SOA_PREFERRED)), first(true)
  {
    resize(NumTargets, NumTargets);
    NumGroups = els.groups();
    PairID.resize(NumTargets, NumTargets);
    for (int i = 0; i < NumTargets; ++i)
      for (int j = 0; j < NumTargets; ++j)
        PairID(i, j) = els.GroupID[i] * NumGroups + els.GroupID[j];
    RadFun.resize(NumGroups * NumGroups, 0);
    offsetPrms.resize(NumGroups * NumGroups, 0);
  }

  ~Backflow_ee(){};

  BackflowFunctionBase* makeClone(ParticleSet& tqp)
  {
    Backflow_ee<FT>* clone = new Backflow_ee<FT>(tqp, tqp);
    clone->first           = false;
    clone->resize(NumTargets, NumTargets);
    clone->offsetPrms = offsetPrms;
    clone->numParams  = numParams;
    clone->derivs     = derivs;
    clone->offsetPrms = offsetPrms;
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
        APP_ABORT("Error cloning Backflow_ee object. \n");
      }
    }
    return clone;
  }

  void addFunc(int ia, int ib, FT* rf)
  {
    uniqueRadFun.push_back(rf);
    if (first)
    {
      // initialize all with rf the first time
      for (int i = 0; i < RadFun.size(); i++)
        RadFun[i] = rf;
      first = false;
    }
    else
    {
      RadFun[ia * NumGroups + ib] = rf;
      RadFun[ib * NumGroups + ia] = rf;
    }
  }

  void registerData(WFBufferType& buf)
  {
    FirstOfU = &(UIJ(0, 0)[0]);
    LastOfU  = FirstOfU + OHMMS_DIM * NumTargets * NumTargets;
    FirstOfA = &(AIJ(0, 0)[0]);
    LastOfA  = FirstOfA + OHMMS_DIM * OHMMS_DIM * NumTargets * NumTargets;
    FirstOfB = &(BIJ(0, 0)[0]);
    LastOfB  = FirstOfB + OHMMS_DIM * NumTargets * NumTargets;
    buf.add(FirstOfU, LastOfU);
    buf.add(FirstOfA, LastOfA);
    buf.add(FirstOfB, LastOfB);
  }

  void reportStatus(std::ostream& os)
  {
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

  inline void acceptMove(int iat, int UpdateMode)
  {
    int num;
    switch (UpdateMode)
    {
    case ORB_PBYP_RATIO:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        UIJ(iat, i) = UIJ_temp[i];
        UIJ(i, iat) = -1.0 * UIJ_temp[i];
      }
      break;
    case ORB_PBYP_PARTIAL:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        UIJ(iat, i) = UIJ_temp[i];
        UIJ(i, iat) = -1.0 * UIJ_temp[i];
      }
      num = AIJ.rows();
      for (int i = 0; i < num; i++)
      {
        AIJ(iat, i) = AIJ_temp[i];
        AIJ(i, iat) = AIJ_temp[i];
      }
      break;
    case ORB_PBYP_ALL:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        UIJ(iat, i) = UIJ_temp[i];
        UIJ(i, iat) = -1.0 * UIJ_temp[i];
      }
      num = AIJ.rows();
      for (int i = 0; i < num; i++)
      {
        AIJ(iat, i) = AIJ_temp[i];
        AIJ(i, iat) = AIJ_temp[i];
      }
      num = BIJ.rows();
      for (int i = 0; i < num; i++)
      {
        BIJ(iat, i) = BIJ_temp[i];
        BIJ(i, iat) = -1.0 * BIJ_temp[i];
      }
      break;
    default:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        UIJ(iat, i) = UIJ_temp[i];
        UIJ(i, iat) = -1.0 * UIJ_temp[i];
      }
      num = AIJ.rows();
      for (int i = 0; i < num; i++)
      {
        AIJ(iat, i) = AIJ_temp[i];
        AIJ(i, iat) = AIJ_temp[i];
      }
      num = BIJ.rows();
      for (int i = 0; i < num; i++)
      {
        BIJ(iat, i) = BIJ_temp[i];
        BIJ(i, iat) = -1.0 * BIJ_temp[i];
      }
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

  /** calculate quasi-particle coordinates only
   */
  inline void evaluate(const ParticleSet& P, ParticleSet& QP)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluate(P,QP) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[PairID(i, j)]->evaluate(myTable.r(nn), du, d2u);
        PosType u    = uij * myTable.dr(nn);
        QP.R[i] -= u; // dr(ij) = r_j-r_i
        QP.R[j] += u;
        UIJ(j, i) = u;
        UIJ(i, j) = -1.0 * u;
      }
    }
#endif
  }

  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradVector_t& Bmat, HessMatrix_t& Amat)
  {
    APP_ABORT("This shouldn't be called: Backflow_ee::evaluate(Bmat)");
    PosType du, d2u, temp;
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluate(P,QP,Bmat_vec,Amat) not implemented for SoA\n");
#else
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[PairID(i, j)]->evaluate(myTable.r(nn), du, d2u);
        PosType u    = uij * myTable.dr(nn);
        // UIJ = eta(r) * (r_i - r_j)
        UIJ(j, i) = u;
        UIJ(i, j) = -1.0 * u;
        du *= myTable.rinv(nn);
        QP.R[i] -= u;
        QP.R[j] += u;
//          Amat(i,j) -= du*(outerProduct(myTable.dr(nn),myTable.dr(nn)));
#if OHMMS_DIM == 3
        Amat(i, j)[0] -= uij;
        Amat(i, j)[4] -= uij;
        Amat(i, j)[8] -= uij;
#elif OHMMS_DIM == 2
        Amat(i, j)[0] -= uij;
        Amat(i, j)[3] -= uij;
#endif
        Amat(j, i) += Amat(i, j);
        Amat(i, i) -= Amat(i, j);
        Amat(j, j) -= Amat(i, j);
        // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
        u = 2.0 * (d2u + (OHMMS_DIM + 1.0) * du) * myTable.dr(nn);
        Bmat[i] -= u;
        Bmat[j] += u;
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
    for (int ig = 0; ig < NumGroups; ++ig)
    {
      for (int iat = P.first(ig), last = P.last(ig); iat < last; ++iat)
      {
        const auto& dist  = myTable.Distances[iat];
        const auto& displ = myTable.Displacements[iat];
        for (int jat = 0; jat < iat; ++jat)
        {
          if (dist[jat] > 0)
          {
            RealType uij = RadFun[PairID(iat, jat)]->evaluate(dist[jat], du, d2u);
            du /= dist[jat];
            PosType u     = uij * displ[jat];
            UIJ(jat, iat) = u;
            UIJ(iat, jat) = -1.0 * u;
            QP.R[iat] -= u;
            QP.R[jat] += u;
            HessType& hess = AIJ(iat, jat);
            hess           = du * outerProduct(displ[jat], displ[jat]);
#if OHMMS_DIM == 3
            hess[0] += uij;
            hess[4] += uij;
            hess[8] += uij;
#elif OHMMS_DIM == 2
            hess[0] += uij;
            hess[3] += uij;
#endif
            AIJ(jat, iat) = hess;
            Amat(iat, iat) += hess;
            Amat(jat, jat) += hess;
            Amat(iat, jat) -= hess;
            Amat(jat, iat) -= hess;
            GradType& grad = BIJ(jat, iat); // dr = r_j - r_i
            grad           = (d2u + (OHMMS_DIM + 1) * du) * displ[jat];
            BIJ(iat, jat)  = -1.0 * grad;
            Bmat_full(iat, iat) -= grad;
            Bmat_full(jat, jat) += grad;
            Bmat_full(iat, jat) += grad;
            Bmat_full(jat, iat) -= grad;
          }
        }
      }
    }
#else
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[PairID(i, j)]->evaluate(myTable.r(nn), du, d2u);
        du *= myTable.rinv(nn);
        PosType u = uij * myTable.dr(nn);
        UIJ(j, i) = u;
        UIJ(i, j) = -1.0 * u;
        QP.R[i] -= u;
        QP.R[j] += u;
        HessType& hess = AIJ(i, j);
        hess           = du * outerProduct(myTable.dr(nn), myTable.dr(nn));
#if OHMMS_DIM == 3
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
#elif OHMMS_DIM == 2
        hess[0] += uij;
        hess[3] += uij;
#endif
        AIJ(j, i) = hess;
        Amat(i, i) += hess;
        Amat(j, j) += hess;
        Amat(i, j) -= hess;
        Amat(j, i) -= hess;
        GradType& grad = BIJ(j, i); // dr = r_j - r_i
        grad           = (d2u + (OHMMS_DIM + 1) * du) * myTable.dr(nn);
        BIJ(i, j)      = -1.0 * grad;
        Bmat_full(i, i) -= grad;
        Bmat_full(j, j) += grad;
        Bmat_full(i, j) += grad;
        Bmat_full(j, i) -= grad;
      }
    }
#endif
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP, const std::vector<int>& index)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluatePbyP(P,QP,index_vec) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI = index.size();
    int iat  = index[0];
    for (int i = 1; i < maxI; i++)
    {
      int j = index[i];
      RealType uij = RadFun[PairID(iat, j)]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
    }
#endif
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP)
  {
#ifdef ENABLE_SOA
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < iat; i++)
    {
      // Temp[j].dr1 = (ri - rj)
      RealType uij = RadFun[PairID(iat, i)]->evaluate(myTable.Temp_r[i], du, d2u);
      PosType u    = (UIJ_temp[i] = -uij * myTable.Temp_dr[i]) - UIJ(iat, i);
      newQP[iat] += u;
      newQP[i] -= u;
    }
    for (int i = iat + 1; i < NumTargets; i++)
    {
      // Temp[j].dr1 = (ri - rj)
      RealType uij = RadFun[PairID(iat, i)]->evaluate(myTable.Temp_r[i], du, d2u);
      PosType u    = (UIJ_temp[i] = -uij * myTable.Temp_dr[i]) - UIJ(iat, i);
      newQP[iat] += u;
      newQP[i] -= u;
    }
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < iat; i++)
    {
      // Temp[j].dr1 = (ri - rj)
      RealType uij = RadFun[PairID(iat, i)]->evaluate(myTable.Temp[i].r1, du, d2u);
      PosType u    = (UIJ_temp[i] = uij * myTable.Temp[i].dr1) - UIJ(iat, i);
      newQP[iat] += u;
      newQP[i] -= u;
    }
    for (int i = iat + 1; i < NumTargets; i++)
    {
      // Temp[j].dr1 = (ri - rj)
      RealType uij = RadFun[PairID(iat, i)]->evaluate(myTable.Temp[i].r1, du, d2u);
      PosType u    = (UIJ_temp[i] = uij * myTable.Temp[i].dr1) - UIJ(iat, i);
      newQP[iat] += u;
      newQP[i] -= u;
    }
#endif
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos_t& newQP,
                           const std::vector<int>& index,
                           HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluatePbyP(P,QP,index_vec,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI = index.size();
    int iat  = index[0];
    for (int i = 1; i < maxI; i++)
    {
      int j        = index[i];
      RealType uij = RadFun[PairID(iat, j)]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
      HessType& hess = AIJ_temp[j];
      hess           = (du * myTable.Temp[j].rinv1) * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
#if OHMMS_DIM == 3
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
#elif OHMMS_DIM == 2
      hess[0] += uij;
      hess[3] += uij;
#endif
      HessType dA = hess - AIJ(iat, j);
      Amat(iat, iat) += dA;
      Amat(j, j) += dA;
      Amat(iat, j) -= dA;
      Amat(j, iat) -= dA;
    }
#endif
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP, HessMatrix_t& Amat)
  {
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
#ifdef ENABLE_SOA
    for (int j = 0; j < iat; j++)
    {
      if (myTable.Temp_r[j] > 0)
      {
        RealType uij = RadFun[PairID(iat, j)]->evaluate(myTable.Temp_r[j], du, d2u);
        PosType u    = (UIJ_temp[j] = -uij * myTable.Temp_dr[j]) - UIJ(iat, j);
        newQP[iat] += u;
        newQP[j] -= u;
        HessType& hess = AIJ_temp[j];
        hess           = (du / myTable.Temp_r[j]) * outerProduct(myTable.Temp_dr[j], myTable.Temp_dr[j]);
#if OHMMS_DIM == 3
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
#elif OHMMS_DIM == 2
        hess[0] += uij;
        hess[3] += uij;
#endif
        HessType dA = hess - AIJ(iat, j);
        Amat(iat, iat) += dA;
        Amat(j, j) += dA;
        Amat(iat, j) -= dA;
        Amat(j, iat) -= dA;
      }
    }
    for (int j = iat + 1; j < NumTargets; j++)
    {
      if (myTable.Temp_r[j] > 0)
      {
        RealType uij = RadFun[PairID(iat, j)]->evaluate(myTable.Temp_r[j], du, d2u);
        PosType u    = (UIJ_temp[j] = -uij * myTable.Temp_dr[j]) - UIJ(iat, j);
        newQP[iat] += u;
        newQP[j] -= u;
        HessType& hess = AIJ_temp[j];
        hess           = (du / myTable.Temp_r[j]) * outerProduct(myTable.Temp_dr[j], myTable.Temp_dr[j]);
#if OHMMS_DIM == 3
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
#elif OHMMS_DIM == 2
        hess[0] += uij;
        hess[3] += uij;
#endif
        HessType dA = hess - AIJ(iat, j);
        Amat(iat, iat) += dA;
        Amat(j, j) += dA;
        Amat(iat, j) -= dA;
        Amat(j, iat) -= dA;
      }
    }
#else
    for (int j = 0; j < iat; j++)
    {
      RealType uij = RadFun[PairID(iat, j)]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
      HessType& hess = AIJ_temp[j];
      hess           = (du * myTable.Temp[j].rinv1) * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
#if OHMMS_DIM == 3
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
#elif OHMMS_DIM == 2
      hess[0] += uij;
      hess[3] += uij;
#endif
      HessType dA = hess - AIJ(iat, j);
      Amat(iat, iat) += dA;
      Amat(j, j) += dA;
      Amat(iat, j) -= dA;
      Amat(j, iat) -= dA;
    }
    for (int j = iat + 1; j < NumTargets; j++)
    {
      RealType uij = RadFun[PairID(iat, j)]->evaluate(myTable.Temp[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * myTable.Temp[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
      HessType& hess = AIJ_temp[j];
      hess           = (du * myTable.Temp[j].rinv1) * outerProduct(myTable.Temp[j].dr1, myTable.Temp[j].dr1);
#if OHMMS_DIM == 3
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
#elif OHMMS_DIM == 2
      hess[0] += uij;
      hess[3] += uij;
#endif
      HessType dA = hess - AIJ(iat, j);
      Amat(iat, iat) += dA;
      Amat(j, j) += dA;
      Amat(iat, j) -= dA;
      Amat(j, iat) -= dA;
    }
#endif
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos_t& newQP,
                           const std::vector<int>& index,
                           GradMatrix_t& Bmat,
                           HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluatePbyP(P,QP,index_vec,Bmat,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    int maxI                                                = index.size();
    int iat                                                 = index[0];
    const std::vector<DistanceTableData::TempDistType>& TMP = myTable.Temp;
    for (int i = 1; i < maxI; i++)
    {
      int j        = index[i];
      RealType uij = RadFun[PairID(iat, j)]->evaluate(TMP[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * TMP[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
      du *= TMP[j].rinv1;
      HessType& hess = AIJ_temp[j];
      hess           = du * outerProduct(TMP[j].dr1, TMP[j].dr1);
#if OHMMS_DIM == 3
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
#elif OHMMS_DIM == 2
      hess[0] += uij;
      hess[3] += uij;
#endif
      HessType dA = hess - AIJ(iat, j);
      Amat(iat, iat) += dA;
      Amat(j, j) += dA;
      Amat(iat, j) -= dA;
      Amat(j, iat) -= dA;
      GradType& grad = BIJ_temp[j]; // dr = r_iat - r_j
      grad           = (d2u + (OHMMS_DIM + 1) * du) * TMP[j].dr1;
      GradType dg    = grad - BIJ(iat, j);
      Bmat(iat, iat) += dg;
      Bmat(j, j) -= dg;
      Bmat(iat, j) -= dg;
      Bmat(j, iat) += dg;
    }
#endif
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           int iat,
                           ParticleSet::ParticlePos_t& newQP,
                           GradMatrix_t& Bmat,
                           HessMatrix_t& Amat)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluatePbyP(P,iat,QP,Bmat,Amat) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    const std::vector<DistanceTableData::TempDistType>& TMP = myTable.Temp;
    for (int j = 0; j < iat; j++)
    {
      RealType uij = RadFun[PairID(iat, j)]->evaluate(TMP[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * TMP[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
      du *= TMP[j].rinv1;
      HessType& hess = AIJ_temp[j];
      hess           = du * outerProduct(TMP[j].dr1, TMP[j].dr1);
#if OHMMS_DIM == 3
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
#elif OHMMS_DIM == 2
      hess[0] += uij;
      hess[3] += uij;
#endif
      HessType dA = hess - AIJ(iat, j);
      Amat(iat, iat) += dA;
      Amat(j, j) += dA;
      Amat(iat, j) -= dA;
      Amat(j, iat) -= dA;
      GradType& grad = BIJ_temp[j]; // dr = r_iat - r_j
      grad           = (d2u + (OHMMS_DIM + 1) * du) * TMP[j].dr1;
      GradType dg    = grad - BIJ(iat, j);
      Bmat(iat, iat) += dg;
      Bmat(j, j) -= dg;
      Bmat(iat, j) -= dg;
      Bmat(j, iat) += dg;
    }
    for (int j = iat + 1; j < NumTargets; j++)
    {
      RealType uij = RadFun[PairID(iat, j)]->evaluate(TMP[j].r1, du, d2u);
      PosType u    = (UIJ_temp[j] = uij * TMP[j].dr1) - UIJ(iat, j);
      newQP[iat] += u;
      newQP[j] -= u;
      du *= TMP[j].rinv1;
      HessType& hess = AIJ_temp[j];
      hess           = du * outerProduct(TMP[j].dr1, TMP[j].dr1);
#if OHMMS_DIM == 3
      hess[0] += uij;
      hess[4] += uij;
      hess[8] += uij;
#elif OHMMS_DIM == 2
      hess[0] += uij;
      hess[3] += uij;
#endif
      HessType dA = hess - AIJ(iat, j);
      Amat(iat, iat) += dA;
      Amat(j, j) += dA;
      Amat(iat, j) -= dA;
      Amat(j, iat) -= dA;
      GradType& grad = BIJ_temp[j]; // dr = r_iat - r_j
      grad           = (d2u + (OHMMS_DIM + 1) * du) * TMP[j].dr1;
      GradType dg    = grad - BIJ(iat, j);
      Bmat(iat, iat) += dg;
      Bmat(j, j) -= dg;
      Bmat(iat, j) -= dg;
      Bmat(j, iat) += dg;
    }
#endif
  }

  /** calculate only Bmat
   *  This is used in pbyp moves, in updateBuffer()
   */
  inline void evaluateBmatOnly(const ParticleSet& P, GradMatrix_t& Bmat_full)
  {
#ifdef ENABLE_SOA
    APP_ABORT("Backflow_ee.h::evaluateBmatOnly(P,QP,Bmat_full) not implemented for SoA\n");
#else
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j        = myTable.J[nn];
        RealType uij = RadFun[PairID(i, j)]->evaluate(myTable.r(nn), du, d2u);
        PosType u    = (d2u + (OHMMS_DIM + 1) * du * myTable.rinv(nn)) * myTable.dr(nn);
        Bmat_full(i, i) -= u;
        Bmat_full(j, j) += u;
        Bmat_full(i, j) += u;
        Bmat_full(j, i) -= u;
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
    RealType du, d2u;
    const auto& myTable = P.getDistTable(myTableIndex_);
#ifdef ENABLE_SOA
    for (int ig = 0; ig < NumGroups; ++ig)
    {
      for (int iat = P.first(ig), last = P.last(ig); iat < last; ++iat)
      {
        const auto& dist  = myTable.Distances[iat];
        const auto& displ = myTable.Displacements[iat];
        for (int jat = 0; jat < iat; ++jat)
        {
          if (dist[jat] > 0)
          {
            RealType uij = RadFun[PairID(iat, jat)]->evaluate(dist[jat], du, d2u);
            //for(int q=0; q<derivs.size(); q++) derivs[q]=0.0; // I believe this is necessary
            //           std::fill(derivs.begin(),derivs.end(),0.0);
            int numParamJU = RadFun[PairID(iat, jat)]->NumParams;
            std::vector<TinyVector<RealType, 3>> derivsju(numParamJU);
            RadFun[PairID(iat, jat)]->evaluateDerivatives(dist[jat], derivsju);
            du /= dist[jat];
            PosType u     = uij * displ[jat];
            UIJ(jat, iat) = u;
            UIJ(iat, jat) = -1.0 * u;
            QP.R[iat] -= u;
            QP.R[jat] += u;
            HessType op    = outerProduct(displ[jat], displ[jat]);
            HessType& hess = AIJ(iat, jat);
            hess           = du * op;
#if OHMMS_DIM == 3
            hess[0] += uij;
            hess[4] += uij;
            hess[8] += uij;
#elif OHMMS_DIM == 2
            hess[0] += uij;
            hess[3] += uij;
#endif
            Amat(iat, iat) += hess;
            Amat(jat, jat) += hess;
            Amat(iat, jat) -= hess;
            Amat(jat, iat) -= hess;
            // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
            // d2u + (ndim+1)*du
            GradType& grad = BIJ(jat, iat); // dr = r_j - r_i
            grad           = (d2u + (OHMMS_DIM + 1) * du) * displ[jat];
            BIJ(iat, jat)  = -1.0 * grad;
            Bmat_full(iat, iat) -= grad;
            Bmat_full(jat, jat) += grad;
            Bmat_full(iat, jat) += grad;
            Bmat_full(jat, iat) -= grad;
            for (int prm = 0, la = indexOfFirstParam + offsetPrms[PairID(iat, jat)]; prm < numParamJU; prm++, la++)
            {
              PosType uk = displ[jat] * derivsju[prm][0];
              Cmat(la, iat) -= uk;
              Cmat(la, jat) += uk;
              Xmat(la, iat, jat) -= (derivsju[prm][1] / dist[jat]) * op;
#if OHMMS_DIM == 3
              Xmat(la, iat, jat)[0] -= derivsju[prm][0];
              Xmat(la, iat, jat)[4] -= derivsju[prm][0];
              Xmat(la, iat, jat)[8] -= derivsju[prm][0];
#elif OHMMS_DIM == 2
              Xmat(la, iat, jat)[0] -= derivsju[prm][0];
              Xmat(la, iat, jat)[3] -= derivsju[prm][0];
#endif
              Xmat(la, jat, iat) += Xmat(la, iat, jat);
              Xmat(la, iat, iat) -= Xmat(la, iat, jat);
              Xmat(la, jat, jat) -= Xmat(la, iat, jat);
              uk = 2.0 * (derivsju[prm][2] + (OHMMS_DIM + 1) * derivsju[prm][1] / dist[jat]) * displ[jat];
              Ymat(la, iat) -= uk;
              Ymat(la, jat) += uk;
            }
          }
        }
      }
    }
#else
    for (int i = 0; i < myTable.sources(); i++)
    {
      for (int nn = myTable.M[i]; nn < myTable.M[i + 1]; nn++)
      {
        int j          = myTable.J[nn];
        RealType uij   = RadFun[PairID(i, j)]->evaluate(myTable.r(nn), du, d2u);
        int numParamJU = RadFun[PairID(i, j)]->NumParams;
        std::vector<TinyVector<RealType, 3>> derivsju(numParamJU);
        RadFun[PairID(i, j)]->evaluateDerivatives(myTable.r(nn), derivsju);
        du *= myTable.rinv(nn);
        PosType u = uij * myTable.dr(nn);
        UIJ(j, i) = u;
        UIJ(i, j) = -1.0 * u;
        QP.R[i] -= u;
        QP.R[j] += u;
        HessType op    = outerProduct(myTable.dr(nn), myTable.dr(nn));
        HessType& hess = AIJ(i, j);
        hess           = du * op;
#if OHMMS_DIM == 3
        hess[0] += uij;
        hess[4] += uij;
        hess[8] += uij;
#elif OHMMS_DIM == 2
        hess[0] += uij;
        hess[3] += uij;
#endif
        Amat(i, i) += hess;
        Amat(j, j) += hess;
        Amat(i, j) -= hess;
        Amat(j, i) -= hess;
        GradType& grad = BIJ(j, i); // dr = r_j - r_i
        grad           = (d2u + (OHMMS_DIM + 1) * du) * myTable.dr(nn);
        BIJ(i, j)      = -1.0 * grad;
        Bmat_full(i, i) -= grad;
        Bmat_full(j, j) += grad;
        Bmat_full(i, j) += grad;
        Bmat_full(j, i) -= grad;
        for (int prm = 0, la = indexOfFirstParam + offsetPrms[PairID(i, j)]; prm < numParamJU; prm++, la++)
        {
          PosType uk = myTable.dr(nn) * derivsju[prm][0];
          Cmat(la, i) -= uk;
          Cmat(la, j) += uk;
          Xmat(la, i, j) -= (derivsju[prm][1] * myTable.rinv(nn)) * op;
#if OHMMS_DIM == 3
          Xmat(la, i, j)[0] -= derivsju[prm][0];
          Xmat(la, i, j)[4] -= derivsju[prm][0];
          Xmat(la, i, j)[8] -= derivsju[prm][0];
#elif OHMMS_DIM == 2
          Xmat(la, i, j)[0] -= derivsju[prm][0];
          Xmat(la, i, j)[3] -= derivsju[prm][0];
#endif
          Xmat(la, j, i) += Xmat(la, i, j);
          Xmat(la, i, i) -= Xmat(la, i, j);
          Xmat(la, j, j) -= Xmat(la, i, j);
          uk = 2.0 * (derivsju[prm][2] + (OHMMS_DIM + 1) * derivsju[prm][1] * myTable.rinv(nn)) * myTable.dr(nn);
          Ymat(la, i) -= uk;
          Ymat(la, j) += uk;
        }
      }
    }
#endif
  }
};

} // namespace qmcplusplus

#endif
