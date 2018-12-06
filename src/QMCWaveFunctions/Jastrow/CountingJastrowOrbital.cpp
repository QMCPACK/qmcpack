#include "QMCWaveFunctions/Jastrow/CountingJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/CountingRegion.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Numerics/Blasf.h"

#include "OhmmsData/AttributeSet.h"
#include "Configuration.h"
#include <cmath>

#include "boost/format.hpp"

namespace qmcplusplus
{

typedef OrbitalBase::ValueType ValueType;
typedef OrbitalBase::RealType RealType;
typedef OrbitalBase::PosType PosType;
typedef OrbitalBase::GradType GradType;
typedef optimize::VariableSet::real_type real_type;
typedef optimize::VariableSet opt_variables_type;

CountingJastrowOrbital::CountingJastrowOrbital(ParticleSet& els) : num_els(els.getTotalNum()) 
{
  Optimizable = true;
}


CountingJastrowOrbital::~CountingJastrowOrbital() {}

// construct from xml file
bool CountingJastrowOrbital::put(xmlNodePtr cur)
{
  // give a report of where we are
  app_log() << "CountingJastrowOrbital::put" << std::endl;
  // get debug option
  OhmmsAttributeSet rAttrib;
  std::string debug_str = "no";
  std::string debug_seqstr = "5", debug_perstr = "10000";
  rAttrib.add(debug_str,"debug");
  rAttrib.add(debug_seqstr,"debug_sequence");
  rAttrib.add(debug_perstr,"debug_period");
  rAttrib.put(cur);
  std::transform(debug_str.begin(), debug_str.end(), debug_str.begin(), (int (*)(int))std::tolower);
  debug = (debug_str == "true");
  debug_seqlen = std::stoi(debug_seqstr);
  debug_period = std::stoi(debug_perstr);

  cur = cur->xmlChildrenNode;
  bool put_F = false, put_G = false, put_C = false;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "var")
    {
      OhmmsAttributeSet rAttrib;
      std::string name;
      rAttrib.add(name,"name");
      std::string opt = "no";
      rAttrib.add(opt,"opt");
      rAttrib.add(opt,"optimize");
      rAttrib.put(cur);
      // make options lowercase
      std::transform(name.begin(), name.end(), name.begin(), (int (*)(int))std::tolower);
      std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))std::tolower);
      if(name == "f")
      {
        // class variable
        opt_F = (opt == "yes" || opt == "true");
        // check if we take in symmetric component or entire matrix
        std::string form = "upper_triang";
        rAttrib.add(form,"form");
        std::transform(form.begin(), form.end(), form.begin(), (int (*)(int))std::tolower);
        if(form == "upper_triang")
        {
          std::vector<RealType> Fsym;
          put_F = putContent(Fsym,cur);
          // get dimension of full matrix, check that this works.
          int Fdim = (std::sqrt(8*Fsym.size() + 1) - 1)/2;
          if(Fsym.size() == Fdim*(Fdim+1)/2)
          {
            // set F from symmetric elements
            F.resize(Fdim*Fdim);
            std::function<RealType&(int,int)> _F      = [&](int I, int J)->RealType& { return F[I*Fdim +J]; };
            auto it = Fsym.begin();
            for(int I = 0; I < Fdim; ++I)
            {
              for(int J = I; J < Fdim; ++J)
              {
                _F(I,J) = _F(J,I) = (*it);
                ++it;
              }
            }
          }
          else
          {
            std::ostringstream err;
            err << "CountingJastrowOrbital::put: F cannot be the upper-triangular component of a square matrix: " << Fsym.size() << " != " << Fdim*(Fdim+1)/2 << std::endl;
            APP_ABORT(err.str());
          }
        }
        else if (form == "full_matrix")
          put_F = putContent(F, cur);
      }
      if(name == "g")
      {
        // class variable
        opt_G = (opt == "yes" || opt == "true");
        put_G = putContent(G,cur);
      }
    }
    if(cname == "region")
    {
      OhmmsAttributeSet rAttrib;
      std::string opt = "no";
      rAttrib.add(opt,"opt");
      rAttrib.add(opt,"optimize");
      std::string norm = "no";
      rAttrib.add(norm,"norm");
      rAttrib.add(norm,"normalize");
      rAttrib.put(cur);
      // make options lowercase
      std::transform(norm.begin(), norm.end(), norm.begin(), (int (*)(int))std::tolower);
      std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))std::tolower);
      // check type
      typedef CountingRegionBase CR;
      if( norm == "yes" || norm == "true")
      {
        typedef NormalizedCountingRegion CR1;
        C = new CR1(num_els);
        put_C = C->put(cur);
        C_norm = true;
      }
      else // not normalized
      {
        typedef CountingRegion CR2;
        C = new CR2(num_els);
        put_C = C->put(cur);
        C_norm = false;
      }
      opt_C = (opt == "yes" || opt == "true");
    }
    cur = cur->next;
  }
  // ensure that we collected all the info we need
  if( !put_F || !put_C || !(put_G || C_norm) )
  {
    std::ostringstream err;
    err << "CountingJastrowOrbital::put:  required variable unspecified: put_F: " << put_F << ", put_C: " << put_C << ", put_G || C_norm: " << (put_G||C_norm) << std::endl;
    APP_ABORT(err.str());
  }
  initialize();
  return true;
}

// called once at the end of put
void CountingJastrowOrbital::initialize()
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
          myVars.insert(id_F,F[IJ],opt_F,optimize::LOGLINEAR_P);
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
      myVars.insert(id_G,G[I],opt_G,optimize::LOGLINEAR_P);
      opt_index[OPT_G].push_back(I);
      opt_id[OPT_G].push_back(id_G);
    }
  }
  reportStatus(app_log());
}

// print stored variables
void CountingJastrowOrbital::reportStatus(std::ostream& os)
{
  os << std::endl << "CountingJastrowOrbital::reportStatus begin" << std::endl;
  // print F matrix
  os << "  F matrix:" << std::endl;
  for(int I = 0; I < num_regions; ++I)
  {
    os << "  ";
    for(int J = 0, IJ = num_regions*I; J < num_regions; ++J, ++IJ)
    {
      os << boost::format("  %10.5f") % F[IJ];
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
  os << "debug: " << debug << std::endl;
  os << "debug_seqlen: " << debug_seqlen << std::endl;
  os << "debug_period: " << debug_period << std::endl;
  os << "  Optimizable variables:" << std::endl;
  myVars.print(os);
  os << std::endl;
  // print counting region status 
  C->reportStatus(os);
  app_log() << "CountingJastrowOrbital::reportStatus end" << std::endl;
}



// === VMC calls ===  
// called at the beginning of every VMC run
void CountingJastrowOrbital::resetTargetParticleSet(ParticleSet& P) {}

inline RealType CountingJastrowOrbital::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  // calculates logPsi and registers data with Pooled data ( .add)
  // underlying PooledData is a vector which is traversed sequentially
  RealType logValue = evaluateLog(P,P.G,P.L);
  RealType *Jgrad_begin = &Jgrad[0][0];
  RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*3;
  DEBUG_PSIBUFFER(" CountingJastrow::registerData",buf.current());
  buf.add(&Jval,&Jval);
  buf.add(Jlap.begin(),Jlap.end());
  buf.add(Jgrad_begin,Jgrad_end);
  DEBUG_PSIBUFFER(" CountingJastrow::registerData",buf.current());
  return logValue;
}
// called in VMC runs before each step
// copy 
inline void CountingJastrowOrbital::copyFromBuffer(ParticleSet& P,PooledData<RealType>& buf)
{
  // copy buffer to wavefunction variables (.get)
  // buf.get(...), same order as before
  RealType *Jgrad_begin = &Jgrad[0][0];
  RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*3;
  DEBUG_PSIBUFFER(" CountingJastrow::copyFromBuffer ",buf.current());
  buf.get(&Jval,&Jval);
  buf.get(Jlap.begin(),Jlap.end());
  buf.get(Jgrad_begin,Jgrad_end);
  DEBUG_PSIBUFFER(" CountingJastrow::copyFromBuffer ",buf.current());
  return;
}
// called in VMC runs in advanceWalkers: return Gradient of particle at iat
GradType CountingJastrowOrbital::evalGrad(ParticleSet& P, int iat)
{
  evaluateExponents(P);
  return Jgrad[iat];
}
// called in VMC runs in advanceWalkers: calculate wfn and ratio at particle position iat
ValueType CountingJastrowOrbital::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  evaluateTempExponents(P,iat);
  grad_iat += Jgrad_t[iat];
  return std::exp(Jval_t - Jval);
}
// called in VMC runs after each substep if step is accepted
void CountingJastrowOrbital::acceptMove(ParticleSet& P, int iat)
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


// called in VMC runs after each substep if step is rejected
void CountingJastrowOrbital::restore(int iat) { C->restore(iat); }

// called in VMC runs after each set of substeps (after each step)
inline RealType CountingJastrowOrbital::updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch)
{
  // need to recompute?
  RealType logValue = evaluateLog(P,P.G,P.L);
  RealType *Jgrad_begin = &Jgrad[0][0];
  RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*3;
  DEBUG_PSIBUFFER(" CountingJastrow::updateBuffer ",buf.current());
  buf.put(&Jval,&Jval);
  buf.put(Jlap.begin(),Jlap.end());
  buf.put(Jgrad_begin,Jgrad_end);
  DEBUG_PSIBUFFER(" CountingJastrow::updateBuffer ",buf.current());
  return Jval;
}
// called every nBlocksBetweenRecompute: default 0
void CountingJastrowOrbital::recompute(ParticleSet& P)
{
  app_log() << "CountingJastrowOrbital::recompute" << std::endl;
  app_log() << "initial values:" << std::endl;
  app_log() << "Jval: " << Jval << std::endl;
  app_log() << "Jgrad: " << std::endl;
  std::copy(Jgrad.begin(),Jgrad.end(), std::ostream_iterator<GradType>(app_log(),", "));
  app_log() << std::endl << "Jlap: ";
  std::copy(Jlap.begin(),Jlap.end(), std::ostream_iterator<RealType>(app_log(),", "));
  // only recompute internal values
  evaluateExponents(P);
  app_log() << "final values:" << std::endl;
  app_log() << "Jval: " << Jval << std::endl;
  app_log() << "Jgrad: " << std::endl;
  std::copy(Jgrad.begin(),Jgrad.end(), std::ostream_iterator<GradType>(app_log(),", "));
  app_log() << std::endl << "Jlap: ";
  std::copy(Jlap.begin(),Jlap.end(), std::ostream_iterator<RealType>(app_log(),", "));
}

// called after recompute
inline RealType CountingJastrowOrbital::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  RealType logValue = evaluateLog(P,P.G,P.L);
  RealType *Jgrad_begin = &Jgrad[0][0];
  RealType *Jgrad_end = Jgrad_begin + Jgrad.size()*3;
  DEBUG_PSIBUFFER(" CountingJastrow::evaluateLog ",buf.current());
  buf.put(&Jval,&Jval);
  buf.put(Jlap.begin(),Jlap.end());
  buf.put(Jgrad_begin,Jgrad_end);
  DEBUG_PSIBUFFER(" CountingJastrow::evaluateLog",buf.current());
  return Jval;
}

// === linear method calls ===
// called from QMCCostFunctionBase.put: once at the beginning of each qmcDriver

// places local optimizable parameters in the global set
void CountingJastrowOrbital::checkInVariables(opt_variables_type& active)
{
  active.insertFrom(myVars);
  C->checkInVariables(active);
}

// gets the index of local parameters in the global set:
// after myVars.getIndex(active), for p = myVars.Index[i]:
// myVars[i] <=> active[p] correspond to the same parameter
void CountingJastrowOrbital::checkOutVariables(const opt_variables_type& active)
{
  myVars.getIndex(active);
  C->checkOutVariables(active);
}

// called from LMYEngineCost, four times (initial position, three shifts)
// copies parameter values from active to myVars and updates any associated wavefunction data
void CountingJastrowOrbital::resetParameters(const opt_variables_type& active)
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
    F[IJ] = myVars[id];
    F[JI] = myVars[id];
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
  app_log() << "CountingJastrowOrbital::resetParameters" << std::endl;
  reportStatus(app_log());
  app_log() << std::endl;
}




void CountingJastrowOrbital::evaluateExponents(ParticleSet& P)
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

void CountingJastrowOrbital::evaluateExponents_print(std::ostream& os, ParticleSet& P)
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

void CountingJastrowOrbital::evaluateTempExponents(ParticleSet& P, int iat)
{
  // evaluate temporary counting regions  
  C->evaluateTemp(P,iat);
  Jval_t = 0;
  std::fill(Jgrad_t.begin(),Jgrad_t.end(),0);
  std::fill(Jlap_t.begin(),Jlap_t.end(),0);
  std::fill(FCsum_t.begin(),FCsum_t.end(),0);
  std::fill(FCgrad_t.begin(),FCgrad_t.end(),0);
  std::fill(FClap_t.begin(),FClap_t.end(),0);

  // test using blas calls - do this later  
//  std::vector<RealType> FCsum_t_blas;
//  std::vector<GradType> FCgrad_t_blas;
//  std::vector<RealType> FClap_t_blas;
//  FCsum_t_blas.resize(FCsum_t.size(), 0);
//  FCgrad_t_blas.resize(FCgrad_t.size(), 0);
//  FClap_t_blas.resize(FClap_t.size(), 0);
//  dgemv('N', num_regions, num_regions, 1.0, &F[0], num_regions, &C->_sum_t[0], 1, 0.0, &FCsum_t_blas[0], 1 );
//  dgemm('N', 'N', num_regions, num_regions, 3, 1.0, &F[0], num_regions, &C->_grad_t[0][0], 3, 0.0, &FCgrad_t_blas[0][0], 3);
//  dgemv('N', num_regions, num_regions, 1.0, &F[0], num_regions, &C->_lap_t[0], 1, 0.0, &FClap_t_blas[0], 1 );
    
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
//  RealType Jval_t_blas = BLAS::dot(num_regions, C->_sum_t.begin(), )
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

void CountingJastrowOrbital::evaluateTempExponents_print(std::ostream& os, ParticleSet& P, int iat)
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

// called in evaluateDeltaLog when no derivative buffer is available.
RealType CountingJastrowOrbital::evaluateLog(ParticleSet& P, 
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL)
{
  evaluateExponents(P);
  for(int i = 0; i < num_els; ++i)
  {
    dG[i] += Jgrad[i];
    dL[i] += Jlap[i];
  }
  LogValue = Jval;
  return LogValue;
}

// called in evaluateValueAndDerivatives for pseudopotential contributions to local energy derivatives
// returns the wavefunction ratio according to moving the particle at iat.
ValueType CountingJastrowOrbital::ratio(ParticleSet& P, int iat)
{
  evaluateTempExponents(P,iat);
  return std::exp(Jval_t - Jval);
//  APP_ABORT("CountingJastrowOrbital::ratio(P,iat) unimplemented");
}


// derivatives of single-particle update for ecp nonlocal quadrature
void CountingJastrowOrbital::evaluateTempDerivatives(ParticleSet& P, 
                         const opt_variables_type& active, 
                         RealType& ratioval,
                         std::vector<RealType>& dlogpsi_t,
                         int iat,
                         PosType dr)
{
  // assume that current state is determined by having called
  // evaluateDerivatives(P,...) immediately before this function
  P.makeMoveAndCheck(iat,dr);
  ratioval = ratio(P,iat);
  // all non-temp variables are set to values associated with position P
  // all temp (_t) variables are set to values for moved position 
  // evaluate log of F parameter derivs at moved position

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
      RealType dJF_val = x*(C->sum_t(I)*C->sum_t(J));
      dlogpsi_t[ia] += dJF_val;
    }
  }
  // evaluate partial derivatives of G at moved position
  if(opt_G)
  {
    for(int oi = 0; oi < opt_index[OPT_G].size(); ++oi)
    {
      std::string id = opt_id[OPT_G][oi];
      int ia = myVars.getIndex(id);
      if(ia == -1)
        continue; // ignore inactive params
      int I = opt_index[OPT_G][oi];
      RealType dJG_val = C->sum_t(I);
      dlogpsi_t[ia] += dJG_val;
    }
  }

  if(opt_C)
  {
    // difference; easier to calculate than absolute values
    static std::vector<RealType> dCdiff;
    static int max_num_derivs = C->max_num_derivs();
    dCdiff.resize(max_num_derivs*num_regions);
    // easy-index functions for evaluateDerivatives calls
    std::function<RealType&(int,int)> _dCsum  = [&](int I, int p)->RealType&{ return dCsum[p*num_regions + I]; }; 
    std::function<RealType&(int,int)> _dCdiff = [&](int I, int p)->RealType&{ return dCdiff[p*num_regions + I]; }; 
    // pointer to C->C[I]->myVars.Index
    // for pI in { 0 .. C->num_derivs(I) }
    //   dCindex->[pI]  is the index that corresponds to this parameter in active.
    //   i.e., active[dCindex->[pI]] <=> C->C[I]->myVars.Index[pI]
    std::fill(dCdiff.begin(), dCdiff.end(), 0);
    for(int I = 0; I < num_regions; ++I)
    {
      // get the number of active parameters for the Ith counting region
      opt_variables_type I_vars = C->getVars(I); 
      int I_num_derivs = I_vars.size();
      // evaluateTempDerivatives increments difference of derivative to dCdiff 
      C->evaluateTempDerivatives(P, I, iat, _dCdiff);
      // loop over parameters for the Ith counting function
      for(int pI = 0; pI < I_num_derivs; ++pI)
      {
        // index for active optimizable variables
        int ia = I_vars.Index[pI];
        if(ia == -1)
          continue; // ignore inactive
        for(int J = 0; J < num_regions; ++J)
        {
          dlogpsi_t[ia] += (_dCsum(J,pI) + _dCdiff(J,pI))*(2*FCsum_t[J] + G[J]);
        }
      }
    }

  } // end opt_C

  // move particle back to the original position
  P.makeMoveAndCheck(iat,-1.0*dr);
}


// Calculates derivatives for use in the linear method.
void CountingJastrowOrbital::evaluateDerivatives(ParticleSet& P, 
                         const opt_variables_type& active, 
                         std::vector<RealType>& dlogpsi, 
                         std::vector<RealType>& dhpsioverpsi)
{
  evaluateExponents(P);
  // indices, strings
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

  // evaluate partial derivatives of G
  if(opt_G)
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

  // evaluate partial derivatives of C
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
  deriv_print_index = deriv_print_index % debug_period;
  deriv_print_index++;
}


// === DMC method calls === 
// called in DMC runs in advanceWalkers
ValueType CountingJastrowOrbital::ratio(ParticleSet& P, 
                                        int iat,
                                        ParticleSet::ParticleGradient_t& dG,
                                        ParticleSet::ParticleLaplacian_t& dL)
{
  evaluateTempExponents(P,iat);
  RealType logRatio = Jval_t - Jval;
  for(int i = 0; i < num_els; ++i)
  {
    dG[i] += Jgrad_t[i] - Jgrad[i];
    dL[i] += Jlap_t[i] - Jlap[i];
  }
  return std::exp(logRatio);
}

// storageType = 0 ; store everything. 
// checkConfigurations in 
void registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf, int storageType) {}

// TrialWaveFunction::evaluate(P), doesn't seem to actually be used: 
// evaluateLog is used instead in most drivers, but implementation required by OrbitalBase
ValueType CountingJastrowOrbital::evaluate(ParticleSet& P, 
                                           ParticleSet::ParticleGradient_t& G,
                                           ParticleSet::ParticleLaplacian_t& L) 
{
  return std::exp(evaluateLog(P,G,L));
}

// called in VMCCUDA advanceWalkers.
void CountingJastrowOrbital::update(ParticleSet& P, 
                                    ParticleSet::ParticleGradient_t& dG,
                                    ParticleSet::ParticleLaplacian_t& dL, int iat) {}

OrbitalBasePtr CountingJastrowOrbital::makeClone(ParticleSet& P) const
{
  app_log() << "  CountingJastrowOrbital::makeClone" << std::endl;
  CountingJastrowOrbital* cjo = new CountingJastrowOrbital(P);
  // flags
  cjo->debug = debug;
  cjo->debug_seqlen = debug_seqlen;
  cjo->debug_period = debug_period;
  cjo->C_norm = C_norm;
  cjo->opt_F = opt_F;
  cjo->opt_G = opt_G;
  cjo->opt_C = opt_C;
  //resize and copy over cjo parameters
  cjo->F.resize(F.size());
  cjo->G.resize(G.size());
  for(int IJ = 0; IJ < F.size(); IJ++)
    cjo->F[IJ] = F[IJ];
  for(int I = 0; I < G.size(); I++)
    cjo->G[I] = G[I];
  //std::copy(F.begin(), F.end(), cjo->F.begin());
  //std::copy(G.begin(), G.end(), cjo->G.begin());
  // make copy of counting regions
  cjo->C = C->makeClone();
  cjo->initialize();
  cjo->evaluateExponents(P);
  cjo->myVars = myVars;
  return cjo;
}

}

