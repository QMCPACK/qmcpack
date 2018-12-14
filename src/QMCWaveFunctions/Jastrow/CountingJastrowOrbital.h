#ifndef QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H
#define QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/CountingRegion.h"

#include "boost/format.hpp"

namespace qmcplusplus
{

//#template <class RegionType> class CountingJastrowOrbital : public WaveFunctionComponent 
template <class RegionType>
class CountingJastrowOrbital: public WaveFunctionComponent 
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
  std::vector<RealType> G; 
  // Counting Regions
  RegionType* C;

  // Optimization Flags
  bool opt_F;
  bool opt_G;
  bool opt_C;

  // Jastrow intermediate Matrix-vector products 
  std::vector<RealType> FCsum;
  std::vector<GradType> FCgrad;
  std::vector<RealType> FClap;
  
  // grad dot grad and laplacian sums for evaluateDerivatives
  std::vector<RealType> FCggsum;
  std::vector<RealType> FClapsum;

  // Jastrow intermediate Matrix-vector products at proposed position
  std::vector<RealType> FCsum_t;
  std::vector<GradType> FCgrad_t;
  std::vector<RealType> FClap_t;
  
  // Jastrow exponent values and gradients (by particle index)
  RealType Jval; 
  std::vector<GradType> Jgrad;
  std::vector<RealType> Jlap;

  // Jastrow exponent values and gradients at proposed position
  RealType Jval_t;
  std::vector<GradType> Jgrad_t;
  std::vector<RealType> Jlap_t;

  // containers for counting function derivative quantities
  std::vector<RealType> dCsum;
  std::vector<RealType> dCsum_t;
  std::vector<RealType> dCggsum;
  std::vector<RealType> dClapsum;
  std::vector<RealType> dCFCggsum;
  std::vector<int> dCindex;
  
  // first array index for opt_index, opt_id 
  enum opt_var { OPT_F, OPT_G, NUM_OPT_VAR };
  // vectors to store indices and names of active optimizable parameters  
  std::array<std::vector<int>,NUM_OPT_VAR> opt_index; 
  std::array<std::vector<std::string>,NUM_OPT_VAR> opt_id;

//================================================================================

public:
  // constructor
  CountingJastrowOrbital(ParticleSet& P, RegionType* c, const Matrix<RealType>& f, const std::vector<RealType>& g):
    F(f), G(g), C(c)
  {
    num_els = P.getTotalNum();
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
    C->checkInVariables(active);
  }


  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    C->checkOutVariables(active);
  }


  void resetParameters(const opt_variables_type& active)
  {
    int ia, I, IJ, JI;
    std::string id;
    for(int i = 0; i < myVars.size(); ++i)
    {
      ia = myVars.where(i);
      if(ia != -1)
        myVars[i] = active[ia];
    }
    // set F parameters from myVars
    for(int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
    {
      IJ = opt_index[OPT_F][oi];
      JI = num_regions*(IJ%num_regions) + IJ/num_regions;
      id = opt_id[OPT_F][oi];
      F(IJ) = myVars[id];
      F(JI) = myVars[id];
    }
    // set G parameters from myVars
    for(int oi = 0; oi < opt_index[OPT_G].size(); ++oi)
    {
      I = opt_index[OPT_G][oi];
      id = opt_id[OPT_G][oi];
      G[I] = myVars[id];
    }
    // reset parameters for counting regions
    C->resetParameters(active);
  }


  void initialize()
  {
    // allocate memory and assign variables
    num_regions = C->size();
  
    FCsum.resize(num_regions);
    FCgrad.resize(num_regions*num_els);
    FClap.resize(num_regions*num_els);
  
    FCggsum.resize(num_regions);
    FClapsum.resize(num_regions);
    
    FCsum_t.resize(num_regions);
    FCgrad_t.resize(num_regions);
    FClap_t.resize(num_regions);
  
    Jgrad.resize(num_els);
    Jlap.resize(num_els);
  
    Jgrad_t.resize(num_els);
    Jlap_t.resize(num_els);
  
    // set G = 0 if using normalized counting functions
    if(C->normalized)
    {
      G.resize(num_regions);
      //std::fill(G.begin(),G.end(),0);
    }
    // check that F, C dimensions match
    if(F.size() != num_regions*num_regions)
    {
      std::ostringstream err;
      err << "CountingJastrowOrbital::initialize: F, C dimension mismatch: F: " << F.size() << ", C: " << num_regions << std::endl;
      APP_ABORT(err.str());
    }
    // check that G, C dimensions match
    if(G.size() != num_regions)
    {
      std::ostringstream err;
      err << "CountingJastrowOrbital::initialize: G, C dimension mismatch: G: " << G.size() << ", C: " << num_regions << std::endl;
      APP_ABORT(err.str());
    }
  
    // for CountingRegion optimization: don't allocate every evalDeriv call
    int max_num_derivs = C->max_num_derivs();
    dCsum.resize(max_num_derivs*num_regions);
    dCggsum.resize(max_num_derivs*num_regions);
    dClapsum.resize(max_num_derivs*num_regions);
    dCFCggsum.resize(max_num_derivs);
    // register optimizable parameters
    std::ostringstream os;
    std::string id_F, id_G;
    if(opt_F)
    {
      for(int I = 0; I < num_regions; ++I)
        for(int J = I, IJ = I*num_regions + I; J < num_regions; ++J, ++IJ)
        {
          // don't optimize bottom-right corner if regions are normalized
          if(!C->normalized || I < (num_regions - 1))
          {
            os.str("");
            os << "F_" << I << "_" << J;
            id_F = os.str();
            myVars.insert(id_F, F(IJ) ,opt_F);
            opt_index[OPT_F].push_back(IJ);
            opt_id[OPT_F].push_back(id_F);
          }
        }
    }

    // only use G when regions aren't normalized
    if(opt_G && !C->normalized)
    {
      for(int I = 0; I < num_regions; ++I)
      {
        os.str("");
        os << "G_" << I;
        id_G = os.str();
        myVars.insert(id_G,G[I],opt_G);
        opt_index[OPT_G].push_back(I);
        opt_id[OPT_G].push_back(id_G);
      }
    }
    reportStatus(app_log());
  }


  void reportStatus(std::ostream& os)
  {
    os << std::endl << "CountingJastrowOrbital::reportStatus begin" << std::endl;
    // print F matrix
    os << "  F matrix:" << std::endl;
    for(int I = 0; I < num_regions; ++I)
    {
      os << "  ";
      for(int J = 0, IJ = num_regions*I; J < num_regions; ++J, ++IJ)
      {
        os << boost::format("  %10.5f") % F(I,J);
      }
      os << std::endl;
    }
    os << "  opt_F: " << (opt_F?"true":"false");
    // print G vector
    if(!C->normalized)
    {
      os << "  G vector:" << std::endl << "  ";
      for(int I = 0; I < num_regions; ++I)
      {
        os << boost::format("  %10.5f") % G[I];
      }
      os << std::endl;
      os << ", opt_G: " << (opt_G?"true":"false");
    }
    // print additional information
    os << "  num_regions: " << num_regions << ", num_els: " << num_els << std::endl;
    if(debug)
    {
      os << "  debug_seqlen: " << debug_seqlen << std::endl;
      os << "  debug_period: " << debug_period << std::endl;
    }
    os << "  Optimizable variables:" << std::endl;
    myVars.print(os);
    os << std::endl;
    // print counting region status 
    C->reportStatus(os);
    app_log() << "CountingJastrowOrbital::reportStatus end" << std::endl;
  }


  void resetTargetParticleSet(ParticleSet& P)
  {
  }


  RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    evaluateExponents(P);
    for(int i = 0; i < num_els; ++i)
    {
      G[i] += Jgrad[i];
      L[i] += Jlap[i];
    }
    LogValue = Jval;
    return LogValue;
  }


  void recompute(ParticleSet& P)
  {
    evaluateExponents(P);
  }

  void evaluateExponents(ParticleSet& P)
  {
    // evaluate counting regions
    C->evaluate(P);
    std::fill(FCsum.begin(),FCsum.end(),0);
    std::fill(FCgrad.begin(),FCgrad.end(),0);
    std::fill(FClap.begin(),FClap.end(),0);
    Jval = 0;
    std::fill(Jgrad.begin(),Jgrad.end(),0);
    std::fill(Jlap.begin(),Jlap.end(),0);

    std::function<RealType&(int,int)> _F      = [&](int I, int J)->RealType& { return F(I*num_regions +J); };
    std::function<GradType&(int,int)> _FCgrad = [&](int I, int i)->GradType& { return FCgrad[I*num_els + i]; };
    std::function<RealType&(int,int)> _FClap = [&](int I, int i)->RealType& { return FClap[I*num_els + i]; };
    // evaluate FC products
    for(int I = 0; I < num_regions; ++I)
    {
      for(int J = 0; J < num_regions; ++J)
      {
        FCsum[I] += _F(I,J)*C->sum(J); // MV
        for(int i = 0; i < num_els; ++i)
        {
          _FCgrad(I,i) += _F(I,J)*C->grad(J,i); // 3*nels*MV
          _FClap(I,i) += _F(I,J)*C->lap(J,i); // nels*MV
        }
      }
    }
    // evaluate components of J
    for(int I = 0; I < num_regions; ++I)
    {
      Jval += (FCsum[I] + G[I])*C->sum(I); // VV
      for(int i = 0; i < num_els; ++i)
      {
        Jgrad[i] += (2*FCsum[I] + G[I])*C->grad(I,i); // 3*nels*VV
        Jlap[i] += (2*FCsum[I] + G[I])*C->lap(I,i) + 2*dot(_FCgrad(I,i),C->grad(I,i)); // nels*VV
      }
    }
    // print out results every so often
    if(debug)
    {
      static int exp_print_index = 0;
      if(exp_print_index < debug_seqlen)
        evaluateExponents_print(app_log(),P);
      ++exp_print_index;
      exp_print_index = exp_print_index % debug_period;
    }
  }


  void evaluateExponents_print(std::ostream& os, ParticleSet& P)
  {
    // print counting regions
    C->evaluate_print(app_log(),P);
    // FCsum, FCgrad, FClap
    os << "CountingJastrowOrbital::evaluateExponents_print: ";
    os << std::endl << "FCsum: ";
    std::copy(FCsum.begin(),FCsum.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "FCgrad: ";
    std::copy(FCgrad.begin(),FCgrad.end(), std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "FClap: ";
    std::copy(FClap.begin(),FClap.end(), std::ostream_iterator<RealType>(os,", "));
    // Jval, Jgrad, Jlap
    os << std::endl << "Jval: " << Jval;
    os << std::endl << "Jgrad: ";
    std::copy(Jgrad.begin(),Jgrad.end(), std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "Jlap:  ";
    std::copy(Jlap.begin(),Jlap.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl << std::endl;
  }


  void evaluateTempExponents(ParticleSet& P, int iat)
  {
    // evaluate temporary counting regions  
    C->evaluateTemp(P,iat);
    Jval_t = 0;
    std::fill(Jgrad_t.begin(),Jgrad_t.end(),0);
    std::fill(Jlap_t.begin(),Jlap_t.end(),0);
    std::fill(FCsum_t.begin(),FCsum_t.end(),0);
    std::fill(FCgrad_t.begin(),FCgrad_t.end(),0);
    std::fill(FClap_t.begin(),FClap_t.end(),0);
  
    std::function<RealType&(int,int)> _F      = [&](int I, int J)->RealType& { return F(I*num_regions +J); };
    std::function<const GradType&(int,int)> _FCgrad = [&](int I, int i)->const GradType&{ return FCgrad[I*num_els + i] ; };
    // evaluate temp FC arrays
    for(int I = 0; I < num_regions; ++I)
    {
      for(int J = 0; J < num_regions; ++J)
      {
        FCsum_t[I] += _F(I,J)*C->sum_t(J); 
        FCgrad_t[I] += _F(I,J)*C->grad_t(J);
        FClap_t[I] += _F(I,J)*C->lap_t(J);
      }
    }
    // evaluate components of the exponent
    for(int I = 0; I < num_regions; ++I)
    {
      Jval_t += C->sum_t(I)*(FCsum_t[I] + G[I]);
      for(int i = 0; i < num_els; ++i)
      {
        if(i == iat)
        {
          Jgrad_t[i] += C->grad_t(I)*(2*FCsum_t[I] + G[I]);
          Jlap_t[i]  += C->lap_t(I)*(2*FCsum_t[I] + G[I]) + 2*dot(C->grad_t(I),FCgrad_t[I]);
        }
        else
        {
          Jgrad_t[i] += C->grad(I,i)*(2*FCsum_t[I] + G[I]);
          Jlap_t[i]  += C->lap(I,i)*(2*FCsum_t[I] + G[I])  + 2*dot(C->grad(I,i),_FCgrad(I,i));
        }
      }
    }
    // print out results every so often
    if(debug)
    {
      static int expt_print_index = 0;
      if(expt_print_index < debug_seqlen)
        evaluateTempExponents_print(app_log(),P,iat);
      ++expt_print_index;
      expt_print_index = expt_print_index % debug_period;
    }
  }

  void evaluateTempExponents_print(std::ostream& os, ParticleSet& P, int iat)
  {
    // print counting regions
    C->evaluateTemp_print(app_log(),P);
    // FCsum, FCgrad, FClap
    os << "CountingJastrowOrbital::evaluateTempExponents_print: iat: " << iat;
    os << std::endl << "FCsum_t: ";
    std::copy(FCsum_t.begin(),FCsum_t.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "FCgrad_t: ";
    std::copy(FCgrad_t.begin(),FCgrad_t.end(), std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "FClap_t: ";
    std::copy(FClap_t.begin(),FClap_t.end(), std::ostream_iterator<RealType>(os,", "));
    // Jval, Jgrad, Jlap
    os << std::endl << "Jval_t: " << Jval_t;
    os << std::endl << "Jgrad_t: ";
    std::copy(Jgrad_t.begin(),Jgrad_t.end(), std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "Jlap_t:  ";
    std::copy(Jlap_t.begin(),Jlap_t.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl << std::endl;
  }
  
  GradType evalGrad(ParticleSet& P, int iat)
  {
    evaluateExponents(P);
    return Jgrad[iat];
  }

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    evaluateTempExponents(P,iat);
    grad_iat += Jgrad_t[iat];
    return std::exp(Jval_t - Jval);
  }

  void acceptMove(ParticleSet& P, int iat)
  {
    C->acceptMove(P,iat);
    // update values for C, FC to those at proposed position
    std::function<GradType&(int,int)> _FCgrad = [&](int I, int i)->GradType& { return FCgrad[I*num_els + i]; };
    std::function<RealType&(int,int)> _FClap = [&](int I, int i)->RealType& { return FClap[I*num_els + i]; };
    // copy over temporary values
    for(int I = 0; I < num_regions; ++I)
    {
      FCsum[I] = FCsum_t[I];
      _FCgrad(I,iat) = FCgrad_t[I];
      _FClap(I,iat) = FClap_t[I];
    }
    // update exponent values to that at proposed position
    Jval = Jval_t;
    for(int i = 0; i < num_els; ++i)
    {
      Jgrad[i] = Jgrad_t[i];
      Jlap[i] = Jlap_t[i];
    }
  }

  void restore(int iat)
  {
    C->restore(iat);
  }

  ValueType ratio(ParticleSet& P, int iat)
  {
    evaluateTempExponents(P,iat);
    return std::exp(Jval_t - Jval);
  }

  void registerData(ParticleSet& P, WFBufferType& buf)
  {
    RealType logValue = evaluateLog(P,P.G,P.L);
    RealType *Jlap_begin = &Jlap[0];
    RealType *Jlap_end = Jlap_begin + Jlap.size();
    RealType *Jgrad_begin = &Jgrad[0][0];
    RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*DIM;
    DEBUG_PSIBUFFER(" CountingJastrow::registerData",buf.current());
    buf.add(&Jval,&Jval);
    buf.add(Jlap_begin, Jlap_end);
    buf.add(Jgrad_begin,Jgrad_end);
    DEBUG_PSIBUFFER(" CountingJastrow::registerData",buf.current());
  }

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false)
  {
    RealType logValue = evaluateLog(P,P.G,P.L);
    RealType *Jlap_begin = &Jlap[0];
    RealType *Jlap_end = Jlap_begin + Jlap.size();
    RealType *Jgrad_begin = &Jgrad[0][0];
    RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*DIM;
    DEBUG_PSIBUFFER(" CountingJastrow::updateBuffer ",buf.current());
    buf.put(&Jval,&Jval);
    buf.put(Jlap_begin, Jlap_end);
    buf.put(Jgrad_begin,Jgrad_end);
    DEBUG_PSIBUFFER(" CountingJastrow::updateBuffer ",buf.current());
    return Jval;
  }

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    RealType *Jlap_begin = &Jlap[0];
    RealType *Jlap_end = Jlap_begin + Jlap.size();
    RealType *Jgrad_begin = &Jgrad[0][0];
    RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*3;
    DEBUG_PSIBUFFER(" CountingJastrow::copyFromBuffer ",buf.current());
    buf.get(&Jval,&Jval);
    buf.get(Jlap_begin, Jlap_end);
    buf.get(Jgrad_begin,Jgrad_end);
    DEBUG_PSIBUFFER(" CountingJastrow::copyFromBuffer ",buf.current());
    return;
  }

  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const 
  {
    CountingJastrowOrbital* cjo = new CountingJastrowOrbital(tqp, C, F, G);
    cjo->setOptimizable(opt_C || opt_G || opt_F);
    cjo->addOpt(opt_C, opt_G, opt_F);
    cjo->addDebug(debug, debug_seqlen, debug_period);
    cjo->initialize();
    return cjo;
  }

  void evaluateDerivatives(ParticleSet& P, const opt_variables_type& active, 
    std::vector<RealType>& dlogpsi, std::vector<RealType>& dhpsioverpsi)
  {
    evaluateExponents(P);
    // evaluate derivatives of F
    if(opt_F)
    {
      for(int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
      {
  
        std::string id = opt_id[OPT_F][oi];
        int ia = myVars.getIndex(id);
        if(ia == -1)
          continue; // ignore inactive parameters
        int IJ = opt_index[OPT_F][oi];
        int I = IJ/num_regions;
        int J = IJ%num_regions;
        // coefficient due to symmetry of F: \sum\limits_{I} F_{II} C_I^2 + \sum\limits_{J > I} 2 F_{IJ}*C_I*C_J
        RealType x = (I==J)?1:2;
        RealType dJF_val = C->sum(I)*C->sum(J)*x;
        RealType dJF_gg = 0, dJF_lap = 0;
        for(int i = 0; i < num_els; ++i)
        {
           dJF_gg += x*(dot(C->grad(I,i),P.G[i])*C->sum(J) + C->sum(I)*dot(C->grad(J,i),P.G[i]));
           dJF_lap += x*(C->lap(I,i)*C->sum(J) + 2*dot(C->grad(I,i),C->grad(J,i)) + C->lap(J,i)*C->sum(I));
        }
        dlogpsi[ia] += dJF_val;
        dhpsioverpsi[ia] += -0.5*dJF_lap - dJF_gg;
        
      }
    }
  
    // evaluate Derivatives of G
    if(opt_G && !C->normalized)
    {
      RealType dJG_val, dJG_gg, dJG_lap;
      for(int oi = 0; oi < opt_index[OPT_G].size(); ++oi)
      {
        std::string id = opt_id[OPT_G][oi];
        int ia = myVars.getIndex(id);
        if(ia == -1)
          continue; // ignore inactive params
        int I = opt_index[OPT_G][oi];
        RealType dJG_val = C->sum(I);
        RealType dJG_gg = dJG_lap = 0;
        for(int i = 0; i < num_els; ++i)
        {
           dJG_gg += dot(C->grad(I,i),P.G[i]);
           dJG_lap += C->lap(I,i);
        }
        dlogpsi[ia] += dJG_val;
        dhpsioverpsi[ia] += -0.5*dJG_lap - dJG_gg;
      }
    }
  //  // evaluate partial derivatives of C
    static int deriv_print_index = 0;
    if(opt_C)
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
      std::fill(FCggsum.begin(),FCggsum.end(),0);
      std::fill(FClapsum.begin(),FClapsum.end(),0);
  
      // easy-index functions for evaluateDerivatives calls
      std::function<const GradType&(int,int)> _FCgrad = [&](int I, int i)->const GradType&{ return FCgrad[I*num_els + i] ; };
      std::function<const RealType&(int,int)> _FClap  = [&](int I, int i)->const RealType&{ return FClap[I*num_els + i] ; };
      std::function<RealType&(int,int)> _dCsum        = [&](int I, int p)->RealType&{ return dCsum[p*num_regions + I]; }; 
      std::function<RealType&(int,int)> _dCggsum      = [&](int I, int p)->RealType&{ return dCggsum[p*num_regions + I] ; };
      std::function<RealType&(int,int)> _dClapsum     = [&](int I, int p)->RealType&{ return dClapsum[p*num_regions + I] ; };
      // evaluate FCggsum
      for(int I = 0; I < num_regions; ++I)
      {
        for(int i = 0; i < num_els; ++i)
        {
          FCggsum[I] += dot(_FCgrad(I,i),P.G[i]);
          FClapsum[I] += _FClap(I,i);
        }
      }
      // pointer to C->C[I]->myVars.Index
      // for pI in { 0 .. C->num_derivs(I) }
      //   dCindex->[pI]  is the index that corresponds to this parameter in active.
      //   i.e., active[dCindex->[pI]] <=> C->C[I]->myVars.Index[pI]
  
      // external print block
      if(debug && deriv_print_index < debug_seqlen)
      {
        app_log() << std::endl << "=== evaluateDerivatives ===" << std::endl;
        app_log() << "== print current exponent values ==" << std::endl;
        evaluateExponents_print(app_log(),P);
        app_log() << "== additional counting function terms ==" << std::endl;
        app_log() << "P.G: ";
        std::copy(P.G.begin(), P.G.end(), std::ostream_iterator<GradType>(app_log(), ", "));
        app_log() << std::endl << "FCgrad: ";
        std::copy(FCgrad.begin(), FCgrad.end(), std::ostream_iterator<GradType>(app_log(), ", "));
        app_log() << std::endl << "FClap: ";
        std::copy(FClap.begin(), FClap.end(), std::ostream_iterator<RealType>(app_log(), ", "));
        app_log() << std::endl << "FCggsum: ";
        std::copy(FCggsum.begin(), FCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
        app_log() << std::endl << "FClapsum: ";
        std::copy(FClapsum.begin(), FClapsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
        app_log() << std::endl;
      }
  
      for(int I = 0; I < num_regions; ++I)
      {
        // get the number of active parameters for the Ith counting region
        opt_variables_type I_vars = C->getVars(I); 
        int I_num_derivs = I_vars.size();
        // clear arrays before each evaluate
        std::fill(dCsum.begin(),dCsum.end(),0);
        std::fill(dCggsum.begin(),dCggsum.end(),0);
        std::fill(dClapsum.begin(),dClapsum.end(),0);
        std::fill(dCFCggsum.begin(),dCFCggsum.end(),0);
        // evaluate all derivatives for the Ith counting function
        C->evaluateDerivatives(P, I, _FCgrad, _dCsum, _dCggsum, _dClapsum, dCFCggsum);
        if(debug && deriv_print_index < debug_seqlen)
        {
          // print out current index information
          app_log() << std::endl;
          app_log() << "  == evaluateDerivatives for counting region " << I << ", num_derivs: " << I_num_derivs << " ==" << std::endl;
          app_log() << "  Indices: ";
          std::copy(I_vars.Index.begin(), I_vars.Index.end(), std::ostream_iterator<int>(app_log(),", "));
          app_log() << std::endl << "  Names: ";
          for(auto it = I_vars.NameAndValue.begin(); it != I_vars.NameAndValue.end(); ++it)
            app_log() << (*it).first << ", ";
          app_log() << std::endl << "  Values: ";
          for(auto it = I_vars.NameAndValue.begin(); it != I_vars.NameAndValue.end(); ++it)
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
        for(int pI = 0; pI < I_num_derivs; ++pI)
        {
          // index for active optimizable variables
          int ia = I_vars.Index[pI];
          if(ia == -1)
            continue; // ignore inactive
          // middle laplacian term: 
          dhpsioverpsi[ia] += -0.5*(4.0*dCFCggsum[pI]);
          if(debug && deriv_print_index < debug_seqlen)
          {
            app_log() << "    == evaluateDerivatives calculations ==" << std::endl;
            app_log() << "    pI: " << pI << ", name: " <<  I_vars.name(pI) <<  ", ia: " << ia << std::endl;
            app_log() << "    dCFCggsum: " << dCFCggsum[pI] << std::endl;
          }
          for(int J = 0; J < num_regions; ++J)
          {
            dlogpsi[ia] += _dCsum(J,pI)*(2*FCsum[J] + G[J]);
            // grad dot grad terms
            dhpsioverpsi[ia] += -1.0*( _dCggsum(J,pI)*(2.0*FCsum[J] + G[J]) + _dCsum(J,pI)*2.0*FCggsum[J]  );
            // outer laplacian terms
            dhpsioverpsi[ia] += -0.5*( 2.0*_dCsum(J,pI)*FClapsum[J] + _dClapsum(J,pI)*(2.0*FCsum[J] + G[J]) ) ;
            if(debug && deriv_print_index < debug_seqlen)
            {
              app_log() << "      J: " << J << std::endl;
              app_log() << "      dlogpsi term          : " << _dCsum(J,pI)*(2*FCsum[J] + G[J]) << std::endl;
              app_log() << "      dhpsi/psi, graddotgrad: " << -1.0*( _dCggsum(J,pI)*(2.0*FCsum[J] + G[J]) + _dCsum(J,pI)*2.0*FCggsum[J]  ) << std::endl;
              app_log() << "      dhpsi/psi, laplacian  : " << -0.5*( 2.0*_dCsum(J,pI)*FClapsum[J] + _dClapsum(J,pI)*(2.0*FCsum[J] + G[J]) ) << std::endl;
            }
  
  
          }
        }
      }
  
    } // end opt_C
    // increment and modulo deriv_print_index
    if(debug)
    {
      deriv_print_index = deriv_print_index % debug_period;
      deriv_print_index++;
    }
  }

  void addOpt(bool opt_C_flag, bool opt_G_flag, bool opt_F_flag)
  {
    opt_F = opt_F_flag;
    opt_G = opt_G_flag;
    opt_C = opt_C_flag;
  }

  void addDebug(bool debug_flag, int seqlen, int period)
  {
    debug = debug_flag;
    debug_seqlen = seqlen;
    debug_period = period;
  }

};

}

#endif
