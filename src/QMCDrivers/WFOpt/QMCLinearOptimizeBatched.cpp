//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCLinearOptimizeBatched.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBatched.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#include "Numerics/LinearFit.h"
#include <iostream>
#include <fstream>


namespace qmcplusplus
{
QMCLinearOptimizeBatched::QMCLinearOptimizeBatched(const ProjectData& project_data,
                                                   MCWalkerConfiguration& w,
                                                   QMCDriverInput&& qmcdriver_input,
                                                   VMCDriverInput&& vmcdriver_input,
                                                   MCPopulation&& population,
                                                   SampleStack& samples,
                                                   Communicate* comm,
                                                   const std::string& QMC_driver_type)
    : QMCDriverNew(project_data,
                   std::move(qmcdriver_input),
                   std::move(population),
                   "QMCLinearOptimizeBatched::",
                   comm,
                   "QMCLinearOptimizeBatched"),
      param_tol(1e-4),
      wfNode(NULL),
      optNode(NULL),
      vmcdriver_input_(vmcdriver_input),
      samples_(samples),
      generate_samples_timer_(
          *timer_manager.createTimer("QMCLinearOptimizeBatched::GenerateSamples", timer_level_medium)),
      initialize_timer_(*timer_manager.createTimer("QMCLinearOptimizeBatched::Initialize", timer_level_medium)),
      eigenvalue_timer_(*timer_manager.createTimer("QMCLinearOptimizeBatched::Eigenvalue", timer_level_medium)),
      line_min_timer_(*timer_manager.createTimer("QMCLinearOptimizeBatched::Line_Minimization", timer_level_medium)),
      cost_function_timer_(*timer_manager.createTimer("QMCLinearOptimizeBatched::CostFunction", timer_level_medium)),
      W(w)
{
  //     //set the optimization flag
  qmc_driver_mode_.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  m_param.add(param_tol, "alloweddifference");
  //Set parameters for line minimization:
}

QMCLinearOptimizeBatched::RealType QMCLinearOptimizeBatched::getLowestEigenvector(Matrix<RealType>& A,
                                                                                  Matrix<RealType>& B,
                                                                                  std::vector<RealType>& ev)
{
  int Nl(ev.size());
  //   Getting the optimal worksize
  char jl('N');
  char jr('V');
  std::vector<RealType> alphar(Nl), alphai(Nl), beta(Nl);
  Matrix<RealType> eigenT(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<RealType> work(1);
  RealType tt(0);
  int t(1);
  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], &tt, &t, eigenT.data(),
               &Nl, &work[0], &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], &tt, &t, eigenT.data(),
               &Nl, &work[0], &lwork, &info);
  if (info != 0)
  {
    APP_ABORT("Invalid Matrix Diagonalization Function!");
  }
  std::vector<std::pair<RealType, int>> mappedEigenvalues(Nl);
  for (int i = 0; i < Nl; i++)
  {
    RealType evi(alphar[i] / beta[i]);
    if (std::abs(evi) < 1e10)
    {
      mappedEigenvalues[i].first  = evi;
      mappedEigenvalues[i].second = i;
    }
    else
    {
      mappedEigenvalues[i].first  = std::numeric_limits<RealType>::max();
      mappedEigenvalues[i].second = i;
    }
  }
  std::sort(mappedEigenvalues.begin(), mappedEigenvalues.end());
  for (int i = 0; i < Nl; i++)
    ev[i] = eigenT(mappedEigenvalues[0].second, i) / eigenT(mappedEigenvalues[0].second, 0);
  return mappedEigenvalues[0].first;
}


QMCLinearOptimizeBatched::RealType QMCLinearOptimizeBatched::getLowestEigenvector(Matrix<RealType>& A,
                                                                                  std::vector<RealType>& ev)
{
  int Nl(ev.size());
  //   Getting the optimal worksize
  RealType zerozero = A(0, 0);
  char jl('N');
  char jr('V');
  std::vector<RealType> alphar(Nl), alphai(Nl), beta(Nl);
  Matrix<RealType> eigenT(Nl, Nl);
  Matrix<RealType> eigenD(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<RealType> work(1);
  LAPACK::geev(&jl, &jr, &Nl, A.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenT.data(), &Nl, &work[0],
               &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);

  LAPACK::geev(&jl, &jr, &Nl, A.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenT.data(), &Nl, &work[0],
               &lwork, &info);
  if (info != 0)
  {
    APP_ABORT("Invalid Matrix Diagonalization Function!");
  }
  std::vector<std::pair<RealType, int>> mappedEigenvalues(Nl);
  for (int i = 0; i < Nl; i++)
  {
    RealType evi(alphar[i]);
    if ((evi < zerozero) && (evi > (zerozero - 1e2)))
    {
      mappedEigenvalues[i].first  = (evi - zerozero + 2.0) * (evi - zerozero + 2.0);
      mappedEigenvalues[i].second = i;
    }
    else
    {
      mappedEigenvalues[i].first  = std::numeric_limits<RealType>::max();
      mappedEigenvalues[i].second = i;
    }
  }
  std::sort(mappedEigenvalues.begin(), mappedEigenvalues.end());
  //         for (int i=0; i<4; i++) app_log()<<i<<": "<<alphar[mappedEigenvalues[i].second]<< std::endl;
  for (int i = 0; i < Nl; i++)
    ev[i] = eigenT(mappedEigenvalues[0].second, i) / eigenT(mappedEigenvalues[0].second, 0);
  return alphar[mappedEigenvalues[0].second];
  //     }
}

void QMCLinearOptimizeBatched::getNonLinearRange(int& first, int& last)
{
  std::vector<int> types;
  optTarget->getParameterTypes(types);
  first = 0;
  last  = types.size();
  //assume all non-linear coeffs are together.
  if (types[0] == optimize::LINEAR_P)
  {
    int i(0);
    while (i < types.size())
    {
      if (types[i] == optimize::LINEAR_P)
        first = i;
      i++;
    }
    first++;
  }
  else
  {
    int i(types.size() - 1);
    while (i >= 0)
    {
      if (types[i] == optimize::LINEAR_P)
        last = i;
      i--;
    }
  }
  //     returns the number of non-linear parameters.
  //    app_log()<<"line params: "<<first<<" "<<last<< std::endl;
}

QMCLinearOptimizeBatched::RealType QMCLinearOptimizeBatched::getNonLinearRescale(std::vector<RealType>& dP,
                                                                                 Matrix<RealType>& S)
{
  int first(0), last(0);
  getNonLinearRange(first, last);
  if (first == last)
    return 1.0;
  RealType rescale(1.0);
  RealType xi(0.5);
  RealType D(0.0);
  for (int i = first; i < last; i++)
    for (int j = first; j < last; j++)
      D += S(i + 1, j + 1) * dP[i + 1] * dP[j + 1];
  rescale = (1 - xi) * D / ((1 - xi) + xi * std::sqrt(1 + D));
  rescale = 1.0 / (1.0 - rescale);
  //     app_log()<<"rescale: "<<rescale<< std::endl;
  return rescale;
}

} // namespace qmcplusplus
