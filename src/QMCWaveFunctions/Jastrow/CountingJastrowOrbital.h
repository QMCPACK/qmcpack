#ifndef QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H
#define QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"

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

//  // flag for using normalized counting regions
  bool C_norm;

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
  CountingJastrowOrbital(ParticleSet& P)
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
    if(C_norm)
    {
      G.resize(num_regions);
      std::fill(G.begin(),G.end(),0);
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
          if(!C_norm || I < (num_regions - 1))
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
    if(opt_G && !C_norm)
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
    if(!C_norm)
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
    //num_els = P.getTotalNum();
    //initialize();
    if(dPsi) dPsi->resetTargetParticleSet(P);
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

//================================================================================

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

    std::function<RealType&(int,int)> _F      = [&](int I, int J)->RealType& { return F[I*num_regions +J]; };
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
  
    std::function<RealType&(int,int)> _F      = [&](int I, int J)->RealType& { return F[I*num_regions +J]; };
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



  //void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi_all)
  //{
  //  APP_ABORT("WaveFunctionComponent::evaluateHessian is not implemented in "+ClassName+" class.");
  //}

  
  GradType evalGrad(ParticleSet& P, int iat)
  {
    evaluateExponents(P);
    return Jgrad[iat];
  }


  /////** return the logarithmic gradient for the iat-th particle
  //// * of the source particleset
  //// * @param Pquantum particle set
  //// * @param iat particle index
  //// * @return the gradient of the iat-th particle
  //// */
  //GradType evalGradSource(ParticleSet& P,
  //                                ParticleSet& source,
  //                                int iat);

  /////** Adds the gradient w.r.t. the iat-th particle of the
  //// *  source particleset (ions) of the logarithmic gradient
  //// *  and laplacian w.r.t. the target paritlceset (electrons).
  //// * @param P quantum particle set (electrons)
  //// * @param source classical particle set (ions)
  //// * @param iat particle index of source (ion)
  //// * @param the ion gradient of the elctron gradient
  //// * @param the ion gradient of the elctron laplacian.
  //// * @return the log gradient of psi w.r.t. the source particle iat
  //// */
  //GradType evalGradSource
  //(ParticleSet& P, ParticleSet& source, int iat,
  // TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
  // TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  /** evaluate the ratio of the new to old orbital value
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @param grad_iat Gradient for the active particle
   */
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    evaluateTempExponents(P,iat);
    grad_iat += Jgrad_t[iat];
    return std::exp(Jval_t - Jval);
  }

  /** a move for iat-th particle is accepted. Update the content for the next moves
   * @param P target ParticleSet
   * @param iat index of the particle whose new position was proposed
   */
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



  /** a move for iat-th particle is reject. Restore to the content.
   * @param iat index of the particle whose new position was proposed
   */
  void restore(int iat)
  {
    C->restore(iat);
  }

  /** evalaute the ratio of the new to old orbital value
   *@param P the active ParticleSet
   *@param iat the index of a particle
   *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'}\})\f$
   *
   *Specialized for particle-by-particle move.
   */
  ValueType ratio(ParticleSet& P, int iat)
  {
    evaluateTempExponents(P,iat);
    return std::exp(Jval_t - Jval);
  }

  /** For particle-by-particle move. Requests space in the buffer
   *  based on the data type sizes of the objects in this class.
   * @param P particle set
   * @param buf Anonymous storage
   */
  void registerData(ParticleSet& P, WFBufferType& buf)
  {
    // calculates logPsi and registers data with Pooled data ( .add)
    // underlying PooledData is a vector which is traversed sequentially
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
    return logValue;
  }

  /** For particle-by-particle move. Put the objects of this class
   *  in the walker buffer or forward the memory cursor.
   * @param P particle set
   * @param buf Anonymous storage
   * @param fromscratch request recomputing the precision critical
   *        pieces of wavefunction from scratch
   * @return log value of the wavefunction.
   */
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false)
  {
    // need to recompute?
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

  /** For particle-by-particle move. Copy data or attach memory
   *  from a walker buffer to the objects of this class.
   *  The log value, P.G and P.L contribution from the objects
   *  of this class are also added.
   * @param P particle set
   * @param buf Anonymous storage
   */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    // copy buffer to wavefunction variables (.get)
    // buf.get(...), same order as before
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

  /** make clone
   * @param tqp target Quantum ParticleSet
   * @param deepcopy if true, make a decopy
   *
   * If not true, return a proxy class
   */
  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const { return 0; };


//  virtual void multiplyDerivsByOrbR(std::vector<RealType>& dlogpsi)
//  {
//    RealType myrat = std::exp(LogValue)*std::cos(PhaseValue);
//    for(int j=0; j<myVars.size(); j++)
//    {
//      int loc=myVars.where(j);
//      dlogpsi[loc] *= myrat;
//    }
//  };


  void finalizeOptimization() { }

  ///** evaluate the ratios of one virtual move with respect to all the particles
  // * @param P reference particleset
  // * @param ratios \f$ ratios[i]=\{{\bf R}\}\rightarrow {r_0,\cdots,r_i^p=pos,\cdots,r_{N-1}}\f$
  // */
  //void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.refPtcl
   */
  void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
  }

  // function calls to pass to differential component dPsi
  void evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
    std::vector<ValueType>& ratios, Matrix<ValueType>& dratios)
  {
    evaluateRatios(VP, ratios);
    dPsi->evaluateDerivRatios(VP, optvars, dratios);
  }
  void evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars,
    std::vector<RealType>& dlogpsi, std::vector<RealType>& dhpsioverpsi)
  {
    dPsi->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }

//  void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
//    std::vector<RealType>& dgradlogpsi) 
//  {
//    dPsi->evaluateGradDerivatives(const ParticleSet::ParticleGradient_t&G_in,
//      std::vector<RealType>& dgradlogpsi);
//  }

  bool addRegion(RegionType* CR, Matrix<RealType>* F, std::vector<RealType>* G, bool opt_CR, bool opt_G, bool opt_F)
  {
  }

  bool addDebug(int seqlen, int period)
  {
  }

};

}

#endif
