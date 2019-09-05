//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "Estimators/EstimatorManagerCrowd.h"
#include "Estimators/CollectablesEstimator.h"

namespace qmcplusplus
{

EstimatorManagerCrowd::EstimatorManagerCrowd(EstimatorManagerBase& em)
    : RecordCount(0),
      FieldWidth(20),
      MainEstimatorName(em.MainEstimatorName),
      Options(em.Options),
      Archive(0),
      DebugArchive(0),
      MainEstimator(0),
      Collectables(0),
      EstimatorMap(em.EstimatorMap),
      max4ascii(em.max4ascii)
{
  // For now I'm going to try to refactor away the clone pattern only at the manager level.
  // i.e. not continue into the scalar_estimators and collectables
  for (int i = 0; i < em.Estimators.size(); i++)
    scalar_estimators_.push_back(em.Estimators[i]->clone());
  MainEstimator = scalar_estimators_[EstimatorMap[MainEstimatorName]];
  if (em.Collectables)
    Collectables = em.Collectables->clone();
}


}
