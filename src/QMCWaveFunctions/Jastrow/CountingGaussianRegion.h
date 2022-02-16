//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//
// File created by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_NORMALIZED_GAUSSIAN_REGION_H
#define QMCPLUSPLUS_NORMALIZED_GAUSSIAN_REGION_H

#include "Configuration.h"
#include "QMCWaveFunctions/Jastrow/CountingGaussian.h"
#include "VariableSet.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

class CountingGaussianRegion
{
public:

  using RealType = QMCTraits::RealType;
  //using ValueType = QMCTraits::ValueType;
  using PosType = QMCTraits::PosType;
  using TensorType = QMCTraits::TensorType;

  using real_type = optimize::VariableSet::real_type;
  using opt_variables_type = optimize::VariableSet;

  // counting function pointers
  std::vector<std::unique_ptr<CountingGaussian>> C;

  RealType Nval_t;

  // reference gaussian id, pointer
  std::string Cref_id;
  std::unique_ptr<CountingGaussian> Cref;

  // number of electrons
  int num_els;
  // number of regions
  int num_regions;
  // counting function id
  std::vector<std::string> C_id;

  // value arrays
  Matrix<RealType> val;
  std::vector<RealType> sum;
  Matrix<PosType> grad;
  Matrix<RealType> lap;
  std::vector<RealType> Nval;

  // log values arrays
  Matrix<RealType> Lval;
  Matrix<PosType> Lgrad;
  Matrix<RealType> Llap;
  std::vector<RealType> Lmax;

  // temporary value arrays
  std::vector<RealType> val_t;
  std::vector<RealType> sum_t;
  std::vector<PosType> grad_t;
  std::vector<RealType> lap_t;

  //memory for temporary log value arrays
  std::vector<RealType> Lval_t;
  std::vector<PosType> Lgrad_t;
  std::vector<RealType> Llap_t;
  RealType Lmax_t;

  std::vector<RealType> _dLval_saved;


public:
  CountingGaussianRegion(ParticleSet& P) { num_els = P.getTotalNum(); }
  CountingGaussianRegion(int n) : num_els(n) { }

  int size() const { return num_regions; }

  const opt_variables_type& getVars(int I) { return C[I]->myVars; }

  int total_num_derivs()
  {
    int tnd = 0;
    for(int I = 0; I < C.size(); ++I)
      tnd += getVars(I).size_of_active();
    return tnd;
  }

  int max_num_derivs() const
  {
    auto comp  = [](const auto& a, const auto& b) { return a->myVars.size_of_active() < b->myVars.size_of_active(); };
    auto& Cmax = *(std::max_element(C.begin(), C.end(), comp));
    return Cmax->myVars.size();
  }

  inline RealType& dLval_saved(int I, int p, int i)
  {
    return _dLval_saved[I * max_num_derivs() * num_els + p * num_els + i];
  }

  void addFunc(std::unique_ptr<CountingGaussian> func, std::string fid)
  {
    C_id.push_back(fid);
    C.push_back(std::move(func));
  }

  void initialize()
  {
    num_regions = C.size();

    // resize arrays
    val.resize(num_regions, num_els);
    sum.resize(num_regions);
    grad.resize(num_regions, num_els);
    lap.resize(num_regions, num_els);
    Nval.resize(num_els);

    val_t.resize(num_regions);
    sum_t.resize(num_regions);
    grad_t.resize(num_regions);
    lap_t.resize(num_regions);

    Lval.resize(num_regions, num_els);
    Lgrad.resize(num_regions, num_els);
    Llap.resize(num_regions, num_els);
    Lmax.resize(num_els);

    Lval_t.resize(num_regions);
    Lgrad_t.resize(num_regions);
    Llap_t.resize(num_regions);

    // store log derivative values for single particle moves
    _dLval_saved.resize(max_num_derivs() * num_regions * num_els);
    for(int I = 0; I < C.size(); ++I)
      C[I]->myVars.resetIndex();
  }

  void checkInVariables(opt_variables_type& active)
  {
    for (auto it = C.begin(); it != C.end(); ++it)
      (*it)->checkInVariables(active);
  }
  void checkOutVariables(const opt_variables_type& active)
  {
    for (auto it = C.begin(); it != C.end(); ++it)
      (*it)->checkOutVariables(active);
  }
  void resetParameters(const opt_variables_type& active)
  {
    for (auto it = C.begin(); it != C.end(); ++it)
      (*it)->resetParameters(active);
  }

  void acceptMove(ParticleSet& P, int iat)
  {
    Nval[iat] = Nval_t;
    Lmax[iat] = Lmax_t;
    for (int I = 0; I < num_regions; ++I)
    {
      sum[I]       = sum_t[I];
      val(I, iat)  = val_t[I];
      grad(I, iat) = grad_t[I];
      lap(I, iat)  = lap_t[I];

      Lval(I, iat)  = Lval_t[I];
      Lgrad(I, iat) = Lgrad_t[I];
      Llap(I, iat)  = Llap_t[I];
    }
  }

  void restore(int iat)
  {
    for (int I = 0; I < C.size(); ++I)
      C[I]->restore(iat);
  }

  void reportStatus(std::ostream& os) 
  {
    // print some class variables:
    os << "    Region type: CountingGaussianRegion" << std::endl;
    os << "    Region optimizable parameters: " << total_num_derivs()  << std::endl << std::endl;
    os << "    Counting Functions: " << std::endl;
    for (int I = 0; I < C.size(); ++I)
      C[I]->reportStatus(os);
  }

  std::unique_ptr<CountingGaussianRegion> makeClone()
  {
    // create a new object and clone counting functions
    auto cr = std::make_unique<CountingGaussianRegion>(num_els);
    for (int i = 0; i < C.size(); ++i)
    {
      auto Ci = C[i]->makeClone(C_id[i]);
      cr->addFunc(std::move(Ci), C_id[i]);
    }
    // get the index of the reference gaussian
    cr->Cref_id = Cref_id;
    auto C_id_it  = std::find(C_id.begin(), C_id.end(), Cref_id);
    int ref_index = std::distance(C_id.begin(), C_id_it);
    // initialize each of the counting functions using the reference gaussian
    cr->Cref = cr->C[ref_index]->makeClone(Cref_id + "_ref");
    for(auto& ptr : cr->C)
      ptr->initialize(cr->Cref.get());
    cr->initialize();
    return cr;
  }

  bool put(xmlNodePtr cur)
  {
    // get the reference function
    OhmmsAttributeSet rAttrib;
    Cref_id = "g0";
    rAttrib.add(Cref_id, "reference_id");
    rAttrib.put(cur);
    // loop through array, find where Cref is
    auto C_id_it  = std::find(C_id.begin(), C_id.end(), Cref_id);
    int ref_index = std::distance(C_id.begin(), C_id_it);
    if (C_id_it == C_id.end())
      APP_ABORT("CountingGaussianRegion::put: reference function not found:" +
                (Cref_id == "none" ? " Cref not specified" : "\"" + Cref_id + "\""));
    // make a copy of the reference gaussian
    Cref = C[ref_index]->makeClone(Cref_id + "_ref");
    // initialize with reference gaussian
    for(auto& ptr : C)
      ptr->initialize(Cref.get());
    initialize();
    return true;
  }

  // evaluate using the log of the counting basis
  void evaluate(const ParticleSet& P)
  {
    // clear arrays
    std::fill(val.begin(), val.end(), 0);
    std::fill(sum.begin(), sum.end(), 0);
    std::fill(grad.begin(), grad.end(), 0);
    std::fill(lap.begin(), lap.end(), 0);
    std::fill(Lval.begin(), Lval.end(), 0);
    std::fill(Lgrad.begin(), Lgrad.end(), 0);
    std::fill(Llap.begin(), Llap.end(), 0);
    std::fill(Nval.begin(), Nval.end(), 0);
    std::fill(Lmax.begin(), Lmax.end(), 0);
    // temporary variables: Lval = ln(C), Lgrad = \nabla ln(C), Llap = \nabla^2 ln(C)
    for (int i = 0; i < num_els; ++i)
    {
      for (int I = 0; I < num_regions; ++I)
      {
        C[I]->evaluateLog(P.R[i], Lval(I, i), Lgrad(I, i), Llap(I, i));
        if (Lval(I, i) > Lmax[i])
          Lmax[i] = Lval(I, i);
      }
      // build counting function values; subtract off largest log value
      for (int I = 0; I < num_regions; ++I)
      {
        val(I, i) = std::exp(Lval(I, i) - Lmax[i]);
        Nval[i] += val(I, i);
      }
      PosType gLN_sum  = 0; // \sum\limits_I \nabla L_{Ii} N_{Ii}
      RealType lLN_sum = 0; // \sum\limits_I \nabla^2 L_{Ii} N_{Ii}
      // build normalized counting function value, intermediate values for gradient
      for (int I = 0; I < num_regions; ++I)
      {
        val(I, i) = val(I, i) / Nval[i];
        sum[I]  += val(I, i);
        gLN_sum += Lgrad(I, i) * val(I, i);
        lLN_sum += Llap(I, i) * val(I, i);
      }
      RealType gLgN_sum = 0; // \sum\limits_{I} \nabla L_{Ii} \cdot \nabla N_{Ii}
      // build gradient, intermediate values for laplacian
      for (int I = 0; I < num_regions; ++I)
      {
        grad(I, i) = (Lgrad(I, i) - gLN_sum) * val(I, i);
        gLgN_sum += dot(Lgrad(I, i), grad(I, i));
      }
      //build laplacian
      for (int I = 0; I < num_regions; ++I)
      {
        lap(I, i) = (Llap(I, i) - lLN_sum - gLgN_sum) * val(I, i) + dot(grad(I, i), Lgrad(I, i) - gLN_sum);
      }
    }
  }

  void evaluate_print(std::ostream& os, const ParticleSet& P)
  {
    for (auto it = C.begin(); it != C.end(); ++it)
      (*it)->evaluate_print(os, P);
    os << "CountingGaussianRegions::evaluate_print" << std::endl;
    os << "val: ";
    std::copy(val.begin(), val.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "sum: ";
    std::copy(sum.begin(), sum.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "grad: ";
    std::copy(grad.begin(), grad.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "lap: ";
    std::copy(lap.begin(), lap.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "Nval: ";
    std::copy(Nval.begin(), Nval.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "Lmax: ";
    std::copy(Lmax.begin(), Lmax.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl;
  }


  void evaluateTemp(const ParticleSet& P, int iat)
  {
    // clear arrays
    std::fill(val_t.begin(), val_t.end(), 0);
    std::fill(sum_t.begin(), sum_t.end(), 0);
    std::fill(grad_t.begin(), grad_t.end(), 0);
    std::fill(lap_t.begin(), lap_t.end(), 0);
    std::fill(Lval_t.begin(), Lval_t.end(), 0);
    std::fill(Lgrad_t.begin(), Lgrad_t.end(), 0);
    std::fill(Llap_t.begin(), Llap_t.end(), 0);

    Lmax_t = Lmax[iat];
    Nval_t = 0;
    // temporary variables
    for (int I = 0; I < num_regions; ++I)
    {
      C[I]->evaluateLog(P.getActivePos(), Lval_t[I], Lgrad_t[I], Llap_t[I]);
      if (Lval_t[I] > Lmax_t)
        Lmax_t = Lval_t[I];
    }
    // build counting function values; subtract off largest log value
    for (int I = 0; I < num_regions; ++I)
    {
      val_t[I] = std::exp(Lval_t[I] - Lmax_t);
      Nval_t += val_t[I];
    }
    PosType gLN_sum_t  = 0; // \sum\limits_I \nabla L_{Ii} N_{Ii}
    RealType lLN_sum_t = 0; // \sum\limits_I \nabla^2 L_{Ii} N_{Ii}
    // build normalized counting function value, intermediate values for gradient
    for (int I = 0; I < num_regions; ++I)
    {
      val_t[I] = val_t[I] / Nval_t;
      sum_t[I] = sum[I] + val_t[I] - val(I, iat);
      gLN_sum_t += Lgrad_t[I] * val_t[I];
      lLN_sum_t += Llap_t[I] * val_t[I];
    }
    RealType gLgN_sum_t = 0; // \sum\limits_{I} \nabla L_{Ii} \cdot \nabla N_{Ii}
    // build gradient, intermediate values for laplacian
    for (int I = 0; I < num_regions; ++I)
    {
      grad_t[I] = (Lgrad_t[I] - gLN_sum_t) * val_t[I];
      gLgN_sum_t += dot(Lgrad_t[I], grad_t[I]);
    }
    //build laplacian
    for (int I = 0; I < num_regions; ++I)
    {
      lap_t[I] = (Llap_t[I] - lLN_sum_t - gLgN_sum_t) * val_t[I] + dot(grad_t[I], Lgrad_t[I] - gLN_sum_t);
    }
  }

  void evaluateTemp_print(std::ostream& os, const ParticleSet& P)
  {
    os << "CountingGaussianRegion::evaluateTemp_print" << std::endl;
    os << "val_t: ";
    std::copy(val_t.begin(), val_t.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "sum_t: ";
    std::copy(sum_t.begin(), sum_t.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "grad_t: ";
    std::copy(grad_t.begin(), grad_t.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "lap_t: ";
    std::copy(lap_t.begin(), lap_t.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "Nval_t: " << Nval_t;
    os << std::endl << "Lmax_t: " << Lmax_t;
    os << std::endl;
  }


  // calculates derivatives of single particle move of particle with index iat
  void evaluateTempDerivatives(ParticleSet& P,
                               const int I, // index of the counting function parameter derivatives are associated with
                               int iat,
                               Matrix<RealType>& dNdiff)
  {
    // may assume evaluate and evaluateTemp has already been called
    int num_derivs = getVars(I).size();
    // get log derivatives
    static std::vector<RealType> dLval_t;
    static int mnd = max_num_derivs();
    dLval_t.resize(mnd);

    C[I]->evaluateLogTempDerivatives(P.getActivePos(), dLval_t);
    for (int J = 0; J < num_regions; ++J)
    {
      for (int p = 0; p < num_derivs; ++p)
      {
        RealType val_Ii = (I == J) - val(I, iat);
        RealType val_It = (I == J) - val_t[I];
        dNdiff(J, p)    = val_t[J] * dLval_t[p] * val_It - val(J, iat) * dLval_saved(I, p, iat) * val_Ii;
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P,
                           int I,
                           Matrix<PosType>& FCgrad,
                           Matrix<RealType>& dNsum,
                           Matrix<RealType>& dNggsum,
                           Matrix<RealType>& dNlapsum,
                           std::vector<RealType>& dNFNggsum)
  {
    evaluate(P);
    static std::vector<RealType> dLval;
    static std::vector<PosType> dLgrad;
    static std::vector<RealType> dLlap;
    static int mnd = max_num_derivs();
    dLval.resize(mnd);
    dLgrad.resize(mnd);
    dLlap.resize(mnd);

    int num_derivs = getVars(I).size();
    for (int i = 0; i < num_els; ++i)
    {
      // get log derivatives
      C[I]->evaluateLogDerivatives(P.R[i], dLval, dLgrad, dLlap);
      for (int J = 0; J < num_regions; ++J)
      {
        for (int p = 0; p < num_derivs; ++p)
        {
          RealType val_Ii = (I == J) - val(I, i);

          RealType dNval = val(J, i) * dLval[p] * val_Ii;
          PosType dNgrad =
              grad(J, i) * dLval[p] * val_Ii + val(J, i) * dLgrad[p] * val_Ii - val(J, i) * dLval[p] * grad(I, i);

          RealType dNlap = lap(J, i) * dLval[p] * val_Ii + 2 * dot(grad(J, i), dLgrad[p]) * val_Ii -
              2 * dot(grad(J, i), grad(I, i)) * dLval[p] + val(J, i) * dLlap[p] * val_Ii -
              2 * val(J, i) * dot(dLgrad[p], grad(I, i)) - val(J, i) * dLval[p] * lap(I, i);
          // accumulate
          dLval_saved(I, p, i) = dLval[p];
          PosType grad_i( std::real(P.G[i][0]), std::real(P.G[i][1]), std::real(P.G[i][2]) );
          dNsum(J, p)    += dNval;
          dNggsum(J, p)  += dot(dNgrad, grad_i);
          dNlapsum(J, p) += dNlap;
          dNFNggsum[p]   += dot(dNgrad, FCgrad(J, i));
        }
      }
    }
  } // end evaluateDerivatives
};

} // namespace qmcplusplus
#endif
