#ifndef QMCPLUSPLUS_AFQMC_ONERDMESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_ONERDMESTIMATOR_H

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

class OneRdmEstimator: public EstimatorBase
{

  public:

  using WfnPtr = std::shared_ptr<WavefunctionBase>;

  OneRdmEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo info,
        std::string title, xmlNodePtr cur, WfnPtr wfn) : 
                                    EstimatorBase(info),TG(tg_), full(false), diagonal(false),
                                    row(false), select(false), SDet(info),wfn0(wfn) 
  {
    writer = (TG.getGlobalRank()==0);
    parse(cur);
    setup(); 
  }

  ~OneRdmEstimator() {}

  void accumulate_step(WalkerSetWSetPtr wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSetWSetPtr wset)
  {
    int nW0 = wset->numWalkers(true);
    ComplexType *sm, weight, oa, ob, w, ooa, oob, tmp;
    SPComplexType* dm;
    for (int i=0; i<nW0; i++) {
      if (!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6
          || std::isnan(wset->getWeight(i).real())) continue;
      if (i%ncores_per_TG == core_rank) {
        // Calculate this walkers green's function.
        wfn0->evaluateOneBodyMixedDensityMatrix(wset->getSM(i), walker_density_matrix,oa,ob);
        sm = wset->getWalker(i, w, tmp, ooa, oob);  // "impsampl"
        //sm = wset->getWalker2(i, eloc, oa, ob);     // "estimator"
        weight = w * ((oa*ob)/(ooa*oob));
        if (!std::isfinite(weight.real())) continue; 
        data[0] += weight;
        SPComplexMatrix::iterator itG = walker_density_matrix.begin();
        std::vector<ComplexType>::iterator itd = data.begin() + 1;
        // Update all the elements of the Green's function estimator.
        for (std::vector<int>::iterator it = elements.begin(); it != elements.end(); ++it) {
          *itd += weight * (*(itG + *it));
          ++itd;
        }
      } // walker loop
    }
//    std::fill(tmp_data.begin(), tmp_data.end(), 0.0);
//    MPI_Reduce(data.data(), tmp_data.data(), data.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
//    std::swap(data, tmp_data);
    TG.Global().reduce_in_place_n(data.begin(),data.size(),std::plus<>());
  }

  void tags(std::ofstream& out)
  {
    if (TG.getGlobalRank() == 0) {
        std::vector<int> ij;
        ij.resize(2);
        for (std::vector<int>::iterator it=elements.begin(); it!=elements.end(); ++it) {
          oneDimensionalIndexToTwoDimensions(*it, NMO, ij);
          if (*it < NMO*NMO) {
            out << "Gup_" << ij[0] << "_"  << ij[1] << " ";
          } else {
            out << "Gdown_" << ij[0] << "_"  << ij[1] << " ";
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

  void parse(xmlNodePtr cur)
  {
    if(cur==NULL) return;
    OhmmsAttributeSet oAttrib;
    std::string mode;
    oAttrib.add(mode, "mode");
    oAttrib.add(row_ix, "index");
    oAttrib.put(cur);
    if (mode == "full") {
      full = true;
    } else if (mode == "diagonal") {
      diagonal = true;
    } else if (mode == "row") {
      row = true;
    } else if (mode == "select") {
      select = true;
      std::vector<int> tmpUp, tmpDown;
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while (tcur != NULL) {
        // Populate elements with input indices of the green's function the user want's to
        // estimate.
        // Currently the user inputs a one-dimensional array which we interpret as
        // an array of tuples {(i,j)..} which are associated with G_{ij}.
        OhmmsAttributeSet t_param;
        std::string id;
        t_param.add(id, "id");
        t_param.put(tcur);
        if (id == "up_indices") {
          putContent(tmpUp, tcur);
        } else if (id == "down_indices") {
          putContent(tmpDown, tcur);
        }
        tcur = tcur->next;
      }
      // Convert to a one dimensional indexing scheme.
      // up spins
      for (std::vector<int>::iterator it=tmpUp.begin(); it != tmpUp.end(); it+=2) {
        int i = *it;
        int j = *(it+1);
        int index = twoDimensionalIndexToOneDimension(i, j, NMO);
        elements.push_back(index);
      }
      // down spins
      for (std::vector<int>::iterator it=tmpDown.begin(); it != tmpDown.end(); it+=2) {
        int i = *it;
        int j = *(it+1);
        int index = twoDimensionalIndexToOneDimension(i, j, NMO) + NMO*NMO;
        elements.push_back(index);
      }
    }
  }

  bool setup()
  {
    ncores_per_TG = TG.getNCoresPerTG();
    core_rank = TG.getLocalTGRank();

    int nrows = 2*NMO;

    // Manually construct indexing for some predetermined options.
    if (full) {
      for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < NMO; j++) {
          int ix = twoDimensionalIndexToOneDimension(i, j, NMO);
          elements.push_back(ix);
        }
      }
    } else if (diagonal) {
      for (int i = 0; i < NMO; i++) {
        int ix = twoDimensionalIndexToOneDimension(i, i, NMO);
        elements.push_back(ix);
      }
      for (int i = 0; i < NMO; i++) {
        int ix = twoDimensionalIndexToOneDimension(i+NMO, i, NMO);
        elements.push_back(ix);
      }
    } else if (row) {
      // Only pick elements above the diagonal.
      for (int j = row_ix; j < NMO; j++) {
        int ix = twoDimensionalIndexToOneDimension(row_ix, j, NMO);
        elements.push_back(ix);
      }
      for (int j = row_ix; j < NMO; j++) {
        int ix = twoDimensionalIndexToOneDimension(row_ix+NMO, j, NMO);
        elements.push_back(ix);
      }
    } else {
      // select option is dealt with during parsing of xml.
    }
    // +1 is for walker weight.
    data.resize(elements.size()+1);
    tmp_data.resize(data.size());
    std::fill(data.begin(), data.end(), 0.0);
    walker_density_matrix.resize(nrows, NMO);

    return true;
  }


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

  WfnPtr wfn0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  std::vector<ComplexType> data;
  // for MPI communication.
  std::vector<ComplexType> tmp_data;
  ComplexType gdenom = 0.0;

  RealType weight, weight_sub;
  RealType targetW = 1;
  int core_rank;
  int ncores_per_TG;

  // Print the full mixed estimator for the one-particle reduced density matrix.
  // This is the default behaviour.
  bool full = false;
  // Only print the diagonal.
  bool diagonal = false;
  // Only print a specific row specified on input.
  bool row = false;
  int row_ix;
  // Print a specific selection of elements contained in `elements' array.
  bool select = false;
  // Note that the user will input a two dimensional list of indices (tuples) {(i,j)}
  // which we then convirt to the one dimensional array `elements'.
  std::vector<int> elements;

  SlaterDetOperations SDet;
  // Temporary storage for a given walker's density matrix.
  SPComplexMatrix walker_density_matrix;

};
}

#endif
