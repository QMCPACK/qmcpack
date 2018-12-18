#ifndef QMCPLUSPLUS_AFQMC_WALKERDMESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_WALKERDMESTIMATOR_H

#include <Message/MPIObjectBase.h>
#include "AFQMC/config.h"
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/NewTimer.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.hpp"
#include "AFQMC/Wavefunctions/WavefunctionBase.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Wavefunctions/SlaterDetOperations.hpp"

namespace qmcplusplus
{

class WalkerDMEstimator : public EstimatorBase
{

  public:

  using WfnPtr = std::shared_ptr<WavefunctionBase>;

  WalkerDMEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo info,
        std::string title, xmlNodePtr cur) : 
                        EstimatorBase(info), TG(tg_), SDet(info) 
  {

    ncores_per_TG = TG.getNCoresPerTG();
    core_rank = TG.getLocalTGRank();
    writer = (TG.getGlobalRank()==0);

    for (int i = 0; i < 2*NMO; i++) {
      for (int j = 0; j < NMO; j++) {
        int ix = twoDimensionalIndexToOneDimension(i, j, NMO);
        elements.push_back(ix);
      }
    }

    // (weight, GF)
    data.resize(1+elements.size());
    tmp_data.resize(data.size());
    std::fill(data.begin(), data.end(), 0.0);
    walker_density_matrix.resize(2*NMO, NMO);

  }

  ~WalkerDMEstimator() {}

  void accumulate_step(WalkerSetWSetPtr wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSetWSetPtr wset)
  {
    int nW0 = wset->numWalkers(true);
    ComplexType eloc, weight, oa, ob, w, ooa, oob, tmp, factor;
    ComplexType *SM;
    for (int i=0; i<nW0; i++) {
      if (i%ncores_per_TG == core_rank) {
        if (!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6
            || std::isnan(wset->getWeight(i).real()) || !wset->isPMBufferFull(i)) continue;
        // 1. Construct Green's function / 1RDM.
        SDet.green_function(wset->getSM(i), wset->getSM(i), tmp, walker_density_matrix, true);
        // 2. Evaluate full Green's function.
        // 2.a Get walker weight at timestep M+N
        SM = wset->getWalker(i, w, tmp, ooa, oob);  // "impsampl"
        SM = wset->getWalker2(i, eloc, oa, ob);     // "estimator"
        weight = w * ((oa*ob)/(ooa*oob));
        if ((!std::isfinite(weight.real())) ||
            (!std::isfinite((eloc*weight).real()))) continue;
        data[0] += weight;
        // 3.b Update Green's function estimator.
        SPComplexMatrix::iterator itG = walker_density_matrix.begin();
        std::vector<ComplexType>::iterator itd = data.begin() + 1;
        // Update all the elements of the Green's function estimator.
        for (std::vector<int>::iterator it = elements.begin(); it != elements.end(); ++it) {
          *itd += weight * (*(itG + *it));
          ++itd;
        }
      } // walker loop
    }
    std::fill(tmp_data.begin(), tmp_data.end(), 0.0);
    MPI_Reduce(data.data(), tmp_data.data(), data.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    std::swap(data, tmp_data);
  }

  void tags(std::ofstream& out)
  {
    if (TG.getGlobalRank() == 0) {
        std::vector<int> ij;
        ij.resize(2);
        for (std::vector<int>::iterator it=elements.begin(); it!=elements.end(); ++it) {
          oneDimensionalIndexToTwoDimensions(*it, NMO, ij);
          if (*it < NMO*NMO) {
            out << "GW_up_" << ij[0] << "_"  << ij[1] << " ";
          } else {
            out << "GW_down_" << ij[0] << "_"  << ij[1] << " ";
          }
        }
    }
  }

  void print(std::ofstream& out, WalkerSetWSetPtr wset)
  {
    if (writer) {
      for (int i = 1; i < data.size(); i++) {
        out << std::setprecision(16) << " " << (data[i].real()/data[0].real()) << " ";
      }
    }
    // Zero our estimator array.
    std::fill(data.begin(), data.end(), 0);
  }

  private:

  inline void oneDimensionalIndexToTwoDimensions(int ix, int ncols, std::vector<int>& coordinate)
  {
    coordinate[0] = ix / ncols;
    coordinate[1] = ix % ncols;
  }
  inline int twoDimensionalIndexToOneDimension(int x, int y, int ncols)
  {
    return ncols*x + y;
  }

  afqmc::TaskGroup_& TG;

  bool writer;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  std::vector<ComplexType> data;
  // for MPI communication.
  std::vector<ComplexType> tmp_data;
  ComplexType gdenom = ComplexType(0.0, 0.0);

  RealType weight, weight_sub;
  RealType targetW = 1;
  int core_rank;
  int ncores_per_TG;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one = ComplexType(1.0, 0.0);

  SlaterDetOperations SDet;
  // Temporary storage for a given walker's density matrix.
  SPComplexMatrix walker_density_matrix;
  std::vector<int> elements;

  int num_heads_tg;
  MPI_Comm MPI_COMM_HEAD_OF_NODES;
  MPI_Comm MPI_COMM_NODE_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL_HEADS;

};
}

#endif
