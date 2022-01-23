////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//
// File created by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H
#define QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/CountingGaussianRegion.h"

namespace qmcplusplus
{
template<class RegionType>
class CountingJastrow : public WaveFunctionComponent
{
protected:
  // number of electrons
  int num_els;
  // number of counting regions
  int num_regions;

  // debug options
  bool debug = false;
  int debug_seqlen;
  int debug_period;

  // Jastrow linear coefficients
  Matrix<RealType> F;

  // Counting Regions
  std::unique_ptr<RegionType> C;

  // Optimization Flags
  bool opt_F;
  //bool opt_G;
  bool opt_C;

  // Jastrow intermediate Matrix-vector products
  std::vector<RealType> FCsum;
  Matrix<PosType> FCgrad;
  Matrix<RealType> FClap;

  // grad dot grad and laplacian sums for evaluateDerivatives
  std::vector<RealType> FCggsum;
  std::vector<RealType> FClapsum;

  // Jastrow intermediate Matrix-vector products at proposed position
  std::vector<RealType> FCsum_t;
  std::vector<PosType> FCgrad_t;
  std::vector<RealType> FClap_t;

  // Jastrow exponent values and gradients (by particle index)
  RealType Jval;
  std::vector<PosType> Jgrad;
  std::vector<RealType> Jlap;

  // Jastrow exponent values and gradients at proposed position
  RealType Jval_t;
  std::vector<PosType> Jgrad_t;
  std::vector<RealType> Jlap_t;

  // containers for counting function derivative quantities
  Matrix<RealType> dCsum;
  Matrix<RealType> dCggsum;
  Matrix<RealType> dClapsum;
  std::vector<RealType> dCsum_t;
  std::vector<RealType> dCFCggsum;
  std::vector<int> dCindex;

  // first array index for opt_index, opt_id
  enum opt_var
  {
    OPT_F,
    //OPT_G,
    NUM_OPT_VAR
  };
  // vectors to store indices and names of active optimizable parameters
  std::array<std::vector<int>, NUM_OPT_VAR> opt_index;
  std::array<std::vector<std::string>, NUM_OPT_VAR> opt_id;

  //================================================================================

public:
  // constructor
  CountingJastrow(ParticleSet& P, std::unique_ptr<RegionType> c, const Matrix<RealType>& f)
      : WaveFunctionComponent("CountingJastrow"), F(f), C(std::move(c))
  {
    num_els = P.getTotalNum();
  }

  void checkInVariables(opt_variables_type& active) override
  {
    active.insertFrom(myVars);
    C->checkInVariables(active);
  }


  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.getIndex(active);
    C->checkOutVariables(active);
  }


  void resetParameters(const opt_variables_type& active) override
  {
    int ia, IJ, JI;
    std::string id;
    myVars.getIndex(active);
    for (int i = 0; i < myVars.size(); ++i)
    {
      ia = myVars.where(i);
      if (ia != -1)
        myVars[i] = active[ia];
    }
    // set F parameters from myVars
    for (int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
    {
      IJ         = opt_index[OPT_F][oi];
      JI         = num_regions * (IJ % num_regions) + IJ / num_regions;
      id         = opt_id[OPT_F][oi];
      ia         = active.getLoc(id);
      myVars[id] = active[ia];
      F(IJ)      = std::real(myVars[id]);
      F(JI)      = std::real(myVars[id]);
    }
    // reset parameters for counting regions
    C->resetParameters(active);
  }


  void initialize()
  {
    // allocate memory and assign variables
    num_regions = C->size();

    FCsum.resize(num_regions);
    FCgrad.resize(num_regions, num_els);
    FClap.resize(num_regions, num_els);

    FCggsum.resize(num_regions);
    FClapsum.resize(num_regions);

    FCsum_t.resize(num_regions);
    FCgrad_t.resize(num_regions);
    FClap_t.resize(num_regions);

    Jgrad.resize(num_els);
    Jlap.resize(num_els);

    Jgrad_t.resize(num_els);
    Jlap_t.resize(num_els);

    // check that F, C dimensions match
    if (F.size() != num_regions * num_regions)
    {
      std::ostringstream err;
      err << "CountingJastrow::initialize: F, C dimension mismatch: F: " << F.size() << ", C: " << num_regions
          << std::endl;
      APP_ABORT(err.str());
    }

    // for CountingRegion optimization: don't allocate every evalDeriv call
    int max_num_derivs = C->max_num_derivs();
    dCsum.resize(max_num_derivs, num_regions);
    dCggsum.resize(max_num_derivs, num_regions);
    dClapsum.resize(max_num_derivs, num_regions);
    dCFCggsum.resize(max_num_derivs);
    // register optimizable parameters
    std::ostringstream os;
    std::string id_F;
    for (int I = 0; I < num_regions; ++I)
      for (int J = I, IJ = I * num_regions + I; J < num_regions; ++J, ++IJ)
      {
        os.str("");
        os << "F_" << I << "_" << J;
        id_F = os.str();
        myVars.insert(id_F, F(IJ), (opt_F && I < (num_regions - 1)));
        opt_index[OPT_F].push_back(IJ);
        opt_id[OPT_F].push_back(id_F);
      }
    myVars.resetIndex();
  }


  void reportStatus(std::ostream& os) override
  {
    os << "    Number of counting regions: " << num_regions << std::endl;
    os << "    Total optimizable parameters: " << C->total_num_derivs() + myVars.size_of_active() << std::endl;
    os << "    F matrix optimizable parameters: " << myVars.size_of_active() << std::endl;
    if (debug)
    {
      os << "    Debug sample sequence length: " << debug_seqlen << std::endl;
      os << "    Debug sample periodicity: " << debug_period << std::endl;
    }
    os << std::endl;
    myVars.print(os, 6, true);
    os << std::endl;
    C->reportStatus(os);
  }


  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    evaluateExponents(P);
    for (int i = 0; i < num_els; ++i)
    {
      G[i] += Jgrad[i];
      L[i] += Jlap[i];
    }
    log_value_ = Jval;
    return log_value_;
  }


  void recompute(const ParticleSet& P) override
  {
    evaluateExponents(P);
    log_value_ = Jval;
  }

  void evaluateExponents(const ParticleSet& P)
  {
    // evaluate counting regions
    C->evaluate(P);
    std::fill(FCsum.begin(), FCsum.end(), 0);
    std::fill(FCgrad.begin(), FCgrad.end(), 0);
    std::fill(FClap.begin(), FClap.end(), 0);
    Jval = 0;
    std::fill(Jgrad.begin(), Jgrad.end(), 0);
    std::fill(Jlap.begin(), Jlap.end(), 0);

    // evaluate FC products
    for (int I = 0; I < num_regions; ++I)
    {
      for (int J = 0; J < num_regions; ++J)
      {
        FCsum[I] += F(I, J) * C->sum[J]; // MV
        for (int i = 0; i < num_els; ++i)
        {
          FCgrad(I, i) += F(I, J) * C->grad(J, i); // 3*nels*MV
          FClap(I, i) += F(I, J) * C->lap(J, i);   // nels*MV
        }
      }
    }
    // evaluate components of J
    for (int I = 0; I < num_regions; ++I)
    {
      Jval += FCsum[I] * C->sum[I]; // VV
      for (int i = 0; i < num_els; ++i)
      {
        Jgrad[i] += 2 * FCsum[I] * C->grad(I, i);                                      // 3*nels*VV
        Jlap[i] += 2 * FCsum[I] * C->lap(I, i) + 2 * dot(FCgrad(I, i), C->grad(I, i)); // nels*VV
      }
    }
    // print out results every so often
    if (debug)
    {
      static int exp_print_index = 0;
      if (exp_print_index < debug_seqlen)
        evaluateExponents_print(app_log(), P);
      ++exp_print_index;
      exp_print_index = exp_print_index % debug_period;
    }
  }


  void evaluateExponents_print(std::ostream& os, const ParticleSet& P)
  {
    // print counting regions
    C->evaluate_print(app_log(), P);
    // FCsum, FCgrad, FClap
    os << "CountingJastrow::evaluateExponents_print: ";
    os << std::endl << "FCsum: ";
    std::copy(FCsum.begin(), FCsum.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "FCgrad: ";
    std::copy(FCgrad.begin(), FCgrad.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "FClap: ";
    std::copy(FClap.begin(), FClap.end(), std::ostream_iterator<RealType>(os, ", "));
    // Jval, Jgrad, Jlap
    os << std::endl << "Jval: " << Jval;
    os << std::endl << "Jgrad: ";
    std::copy(Jgrad.begin(), Jgrad.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "Jlap:  ";
    std::copy(Jlap.begin(), Jlap.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << std::endl;
  }


  void evaluateTempExponents(ParticleSet& P, int iat)
  {
    // evaluate temporary counting regions
    C->evaluateTemp(P, iat);
    Jval_t = 0;
    std::fill(Jgrad_t.begin(), Jgrad_t.end(), 0);
    std::fill(Jlap_t.begin(), Jlap_t.end(), 0);
    std::fill(FCsum_t.begin(), FCsum_t.end(), 0);
    std::fill(FCgrad_t.begin(), FCgrad_t.end(), 0);
    std::fill(FClap_t.begin(), FClap_t.end(), 0);

    // evaluate temp FC arrays
    for (int I = 0; I < num_regions; ++I)
    {
      for (int J = 0; J < num_regions; ++J)
      {
        FCsum_t[I] += F(I, J) * C->sum_t[J];
        FCgrad_t[I] += F(I, J) * C->grad_t[J];
        FClap_t[I] += F(I, J) * C->lap_t[J];
      }
    }
    // evaluate components of the exponent
    for (int I = 0; I < num_regions; ++I)
    {
      Jval_t += C->sum_t[I] * FCsum_t[I];
      for (int i = 0; i < num_els; ++i)
      {
        if (i == iat)
        {
          Jgrad_t[i] += C->grad_t[I] * 2 * FCsum_t[I];
          Jlap_t[i] += C->lap_t[I] * 2 * FCsum_t[I] + 2 * dot(C->grad_t[I], FCgrad_t[I]);
        }
        else
        {
          Jgrad_t[i] += C->grad(I, i) * 2 * FCsum_t[I];
          Jlap_t[i] += C->lap(I, i) * 2 * FCsum_t[I] + 2 * dot(C->grad(I, i), FCgrad(I, i));
        }
      }
    }
    // print out results every so often
    if (debug)
    {
      static int expt_print_index = 0;
      if (expt_print_index < debug_seqlen)
        evaluateTempExponents_print(app_log(), P, iat);
      ++expt_print_index;
      expt_print_index = expt_print_index % debug_period;
    }
  }

  void evaluateTempExponents_print(std::ostream& os, ParticleSet& P, int iat)
  {
    // print counting regions
    C->evaluateTemp_print(app_log(), P);
    // FCsum, FCgrad, FClap
    os << "CountingJastrow::evaluateTempExponents_print: iat: " << iat;
    os << std::endl << "FCsum_t: ";
    std::copy(FCsum_t.begin(), FCsum_t.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "FCgrad_t: ";
    std::copy(FCgrad_t.begin(), FCgrad_t.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "FClap_t: ";
    std::copy(FClap_t.begin(), FClap_t.end(), std::ostream_iterator<RealType>(os, ", "));
    // Jval, Jgrad, Jlap
    os << std::endl << "Jval_t: " << Jval_t;
    os << std::endl << "Jgrad_t: ";
    std::copy(Jgrad_t.begin(), Jgrad_t.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "Jlap_t:  ";
    std::copy(Jlap_t.begin(), Jlap_t.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << std::endl;
  }

  GradType evalGrad(ParticleSet& P, int iat) override
  {
    evaluateExponents(P);
    log_value_ = Jval;
    return Jgrad[iat];
  }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    evaluateTempExponents(P, iat);
    grad_iat += Jgrad_t[iat];
    return std::exp(static_cast<PsiValueType>(Jval_t - Jval));
  }

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    C->acceptMove(P, iat);
    // update values for C, FC to those at proposed position
    // copy over temporary values
    for (int I = 0; I < num_regions; ++I)
    {
      FCsum[I]       = FCsum_t[I];
      FCgrad(I, iat) = FCgrad_t[I];
      FClap(I, iat)  = FClap_t[I];
    }
    // update exponent values to that at proposed position
    Jval       = Jval_t;
    log_value_ = Jval;
    for (int i = 0; i < num_els; ++i)
    {
      Jgrad[i] = Jgrad_t[i];
      Jlap[i]  = Jlap_t[i];
    }
  }

  void restore(int iat) override { C->restore(iat); }

  PsiValueType ratio(ParticleSet& P, int iat) override
  {
    evaluateTempExponents(P, iat);
    return std::exp(static_cast<PsiValueType>(Jval_t - Jval));
  }

  void registerData(ParticleSet& P, WFBufferType& buf) override
  {
    LogValueType logValue = evaluateLog(P, P.G, P.L);
    RealType* Jlap_begin  = &Jlap[0];
    RealType* Jlap_end    = Jlap_begin + Jlap.size();
    RealType* Jgrad_begin = &Jgrad[0][0];
    RealType* Jgrad_end   = Jgrad_begin + Jgrad.size() * DIM;
    DEBUG_PSIBUFFER(" CountingJastrow::registerData", buf.current());
    buf.add(&Jval, &Jval);
    buf.add(Jlap_begin, Jlap_end);
    buf.add(Jgrad_begin, Jgrad_end);
    DEBUG_PSIBUFFER(" CountingJastrow::registerData", buf.current());
  }

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override
  {
    LogValueType logValue = evaluateLog(P, P.G, P.L);
    RealType* Jlap_begin  = &Jlap[0];
    RealType* Jlap_end    = Jlap_begin + Jlap.size();
    RealType* Jgrad_begin = &Jgrad[0][0];
    RealType* Jgrad_end   = Jgrad_begin + Jgrad.size() * DIM;
    DEBUG_PSIBUFFER(" CountingJastrow::updateBuffer ", buf.current());
    buf.put(&Jval, &Jval);
    buf.put(Jlap_begin, Jlap_end);
    buf.put(Jgrad_begin, Jgrad_end);
    DEBUG_PSIBUFFER(" CountingJastrow::updateBuffer ", buf.current());
    return Jval;
  }

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override
  {
    RealType* Jlap_begin  = &Jlap[0];
    RealType* Jlap_end    = Jlap_begin + Jlap.size();
    RealType* Jgrad_begin = &Jgrad[0][0];
    RealType* Jgrad_end   = Jgrad_begin + Jgrad.size() * 3;
    DEBUG_PSIBUFFER(" CountingJastrow::copyFromBuffer ", buf.current());
    buf.get(&Jval, &Jval);
    buf.get(Jlap_begin, Jlap_end);
    buf.get(Jgrad_begin, Jgrad_end);
    DEBUG_PSIBUFFER(" CountingJastrow::copyFromBuffer ", buf.current());
    return;
  }

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override
  {
    auto cjc = std::make_unique<CountingJastrow>(tqp, C->makeClone(), F);
    cjc->setOptimizable(opt_C || opt_F);
    cjc->addOpt(opt_C, opt_F);
    cjc->addDebug(debug, debug_seqlen, debug_period);
    cjc->initialize();
    return cjc;
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override
  {
#ifdef QMC_COMPLEX
    APP_ABORT("CountingJastrow::evaluateDerivatives is not available on complex builds.");
#else
    evaluateExponents(P);
    // evaluate derivatives of F
    static int deriv_print_index = 0;
    if (opt_F)
    {
      for (int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
      {
        std::string id = opt_id[OPT_F][oi];
        int ia         = myVars.getIndex(id);
        if (ia == -1)
          continue; // ignore inactive parameters
        int IJ = opt_index[OPT_F][oi];
        int I  = IJ / num_regions;
        int J  = IJ % num_regions;
        // coefficient due to symmetry of F: \sum\limits_{I} F_{II} C_I^2 + \sum\limits_{J > I} 2 F_{IJ}*C_I*C_J
        RealType x       = (I == J) ? 1 : 2;
        RealType dJF_val = C->sum[I] * C->sum[J] * x;
        RealType dJF_gg = 0, dJF_lap = 0;
        for (int i = 0; i < num_els; ++i)
        {
          PosType grad_i(P.G[i]);
          dJF_gg += x * (dot(C->grad(I, i), grad_i) * C->sum[J] + C->sum[I] * dot(C->grad(J, i), grad_i));
          dJF_lap += x * (C->lap(I, i) * C->sum[J] + 2 * dot(C->grad(I, i), C->grad(J, i)) + C->lap(J, i) * C->sum[I]);
        }
        dlogpsi[ia] += dJF_val;
        dhpsioverpsi[ia] += -0.5 * dJF_lap - dJF_gg;
        if (debug && deriv_print_index < debug_seqlen)
        {
          app_log() << "  dJ/dF[" << I << "][" << J << "]; ia: " << ia << ",  dlogpsi: " << dlogpsi[ia]
                    << ", dhpsioverpsi: " << dhpsioverpsi[ia] << std::endl;
        }
      }
    }

    //  // evaluate partial derivatives of C
    if (opt_C)
    {
      // containers for CountingRegions' evaluateDerivatives calculations
      // blocks of dimension n_p x n_C
      // ex: dNsum = \sum\limits_k [ [dC1 / dp1] [dC2 / dp1] .. ]
      // Where each [dCi / dpj] block is n_p x 1 vector of derivatives of parameters
      //   for counting region j as summed over electron coordinate k

      // exception: dNFN ggsum is an n_p x 1 vector since it is an evaluation of a quadratic form:
      // \sum\limits_{kI} [\nabla_k dC_I/dpj] dot [ (F \nabla_k C)_I ]
      // since we have the premultiplied (F\nabla_k C) vector on hand.
      // make a lambda function FCgrad(I,i) which gives the appropriate element of FCgrad[iI]

      // clear some vectors
      std::fill(FCggsum.begin(), FCggsum.end(), 0);
      std::fill(FClapsum.begin(), FClapsum.end(), 0);

      // evaluate FCggsum
      for (int I = 0; I < num_regions; ++I)
      {
        for (int i = 0; i < num_els; ++i)
        {
          PosType grad_i(P.G[i]);
          FCggsum[I] += dot(FCgrad(I, i), grad_i);
          FClapsum[I] += FClap(I, i);
        }
      }
      // pointer to C->C[I]->myVars.Index
      // for pI in { 0 .. C->num_derivs(I) }
      //   dCindex->[pI]  is the index that corresponds to this parameter in active.
      //   i.e., active[dCindex->[pI]] <=> C->C[I]->myVars.Index[pI]

      // external print block
      if (debug && deriv_print_index < debug_seqlen)
      {
        app_log() << std::endl << "=== evaluateDerivatives ===" << std::endl;
        app_log() << "== print current exponent values ==" << std::endl;
        evaluateExponents_print(app_log(), P);
        app_log() << "== additional counting function terms ==" << std::endl;
        app_log() << "P.G: ";
        std::copy(P.G.begin(), P.G.end(), std::ostream_iterator<PosType>(app_log(), ", "));
        app_log() << std::endl << "FCgrad: ";
        std::copy(FCgrad.begin(), FCgrad.end(), std::ostream_iterator<PosType>(app_log(), ", "));
        app_log() << std::endl << "FClap: ";
        std::copy(FClap.begin(), FClap.end(), std::ostream_iterator<RealType>(app_log(), ", "));
        app_log() << std::endl << "FCggsum: ";
        std::copy(FCggsum.begin(), FCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
        app_log() << std::endl << "FClapsum: ";
        std::copy(FClapsum.begin(), FClapsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
        app_log() << std::endl;
      }

      for (int I = 0; I < num_regions; ++I)
      {
        // get the number of active parameters for the Ith counting region
        opt_variables_type I_vars = C->getVars(I);
        int I_num_derivs          = I_vars.size();
        // clear arrays before each evaluate
        std::fill(dCsum.begin(), dCsum.end(), 0);
        std::fill(dCggsum.begin(), dCggsum.end(), 0);
        std::fill(dClapsum.begin(), dClapsum.end(), 0);
        std::fill(dCFCggsum.begin(), dCFCggsum.end(), 0);
        // evaluate all derivatives for the Ith counting function
        C->evaluateDerivatives(P, I, FCgrad, dCsum, dCggsum, dClapsum, dCFCggsum);
        if (debug && deriv_print_index < debug_seqlen)
        {
          // print out current index information
          app_log() << std::endl;
          app_log() << "  == evaluateDerivatives for counting region " << I << ", num_derivs: " << I_num_derivs
                    << " ==" << std::endl;
          app_log() << "  Indices: ";
          std::copy(I_vars.Index.begin(), I_vars.Index.end(), std::ostream_iterator<int>(app_log(), ", "));
          app_log() << std::endl << "  Names: ";
          for (auto it = I_vars.NameAndValue.begin(); it != I_vars.NameAndValue.end(); ++it)
            app_log() << (*it).first << ", ";
          app_log() << std::endl << "  Values: ";
          for (auto it = I_vars.NameAndValue.begin(); it != I_vars.NameAndValue.end(); ++it)
            app_log() << (*it).second << ", ";
          // print out values from evaluate derivatives
          app_log() << std::endl << "  dCsum: ";
          std::copy(dCsum.begin(), dCsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
          app_log() << std::endl << "  dCggsum: ";
          std::copy(dCggsum.begin(), dCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
          app_log() << std::endl << "  dClapsum: ";
          std::copy(dClapsum.begin(), dClapsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
          app_log() << std::endl << "  dCFCggsum: ";
          std::copy(dCFCggsum.begin(), dCFCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
          app_log() << std::endl;
        }
        // loop over parameters for the Ith counting function
        for (int pI = 0; pI < I_num_derivs; ++pI)
        {
          // index for active optimizable variables
          int ia = I_vars.Index[pI];
          if (ia == -1)
            continue; // ignore inactive
          // middle laplacian term:
          dhpsioverpsi[ia] += -0.5 * (4.0 * dCFCggsum[pI]);
          if (debug && deriv_print_index < debug_seqlen)
          {
            app_log() << "    == evaluateDerivatives calculations ==" << std::endl;
            app_log() << "    pI: " << pI << ", name: " << I_vars.name(pI) << ", ia: " << ia << std::endl;
            app_log() << "    dCFCggsum: " << dCFCggsum[pI] << std::endl;
          }
          for (int J = 0; J < num_regions; ++J)
          {
            dlogpsi[ia] += dCsum(J, pI) * (2 * FCsum[J]);
            // grad dot grad terms
            dhpsioverpsi[ia] += -1.0 * (dCggsum(J, pI) * (2.0 * FCsum[J]) + dCsum(J, pI) * 2.0 * FCggsum[J]);
            // outer laplacian terms
            dhpsioverpsi[ia] += -0.5 * (2.0 * dCsum(J, pI) * FClapsum[J] + dClapsum(J, pI) * (2.0 * FCsum[J]));
            if (debug && deriv_print_index < debug_seqlen)
            {
              app_log() << "      J: " << J << std::endl;
              app_log() << "      dlogpsi term          : " << dCsum(J, pI) * (2 * FCsum[J]) << std::endl;
              app_log() << "      dhpsi/psi, graddotgrad: "
                        << -1.0 * (dCggsum(J, pI) * (2.0 * FCsum[J]) + dCsum(J, pI) * 2.0 * FCggsum[J]) << std::endl;
              app_log() << "      dhpsi/psi, laplacian  : "
                        << -0.5 * (2.0 * dCsum(J, pI) * FClapsum[J] + dClapsum(J, pI) * (2.0 * FCsum[J])) << std::endl;
            }
          }
        }
      }

    } // end opt_C
    // increment and modulo deriv_print_index
    if (debug)
    {
      app_log() << "Final derivatives: " << std::endl;
      app_log() << "  F derivatives: " << std::endl;
      for (int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
      {
        std::string id = opt_id[OPT_F][oi];
        int ia         = myVars.getIndex(id);
        if (ia == -1)
          continue; // ignore inactive parameters
        app_log() << "    ia: " << ia << ",  dlogpsi: " << dlogpsi[ia] << ", dhpsioverpsi: " << dhpsioverpsi[ia]
                  << std::endl;
      }
      app_log() << "  C derivatives: " << std::endl;
      for (int I = 0; I < num_regions; ++I)
      {
        app_log() << "    C[" << I << "] derivs: " << std::endl;
        // get the number of active parameters for the Ith counting region
        opt_variables_type I_vars = C->getVars(I);
        int I_num_derivs          = I_vars.size();
        for (int pI = 0; pI < I_num_derivs; ++pI)
        {
          // index for active optimizable variables
          int ia = I_vars.Index[pI];
          if (ia == -1)
            continue; // ignore inactive
          app_log() << "      ia: " << ia << ",  dlogpsi: " << dlogpsi[ia] << ", dhpsioverpsi: " << dhpsioverpsi[ia]
                    << std::endl;
        }
      }
      deriv_print_index = deriv_print_index % debug_period;
      deriv_print_index++;
    }
#endif
  }

  void addOpt(bool opt_C_flag, bool opt_F_flag)
  {
    opt_F = opt_F_flag;
    opt_C = opt_C_flag;
  }

  void addDebug(bool debug_flag, int seqlen, int period)
  {
    debug        = debug_flag;
    debug_seqlen = seqlen;
    debug_period = period;
  }
};

} // namespace qmcplusplus

#endif
