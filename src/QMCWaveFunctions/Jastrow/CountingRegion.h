#ifndef QMCPLUSPLUS_COUNTINGREGION_H
#define QMCPLUSPLUS_COUNTINGREGION_H
#include "CountingFunctors.h"
#include "Optimize/VariableSet.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "OhmmsData/AttributeSet.h"

#include <functional>
#include "boost/format.hpp"

namespace qmcplusplus
{

typedef OrbitalBase::ValueType ValueType;
typedef OrbitalBase::RealType RealType;
typedef OrbitalBase::PosType PosType;
typedef OrbitalBase::GradType GradType;
typedef optimize::VariableSet::real_type real_type;
typedef optimize::VariableSet opt_variables_type;

// CountingRegion interface
//   Primarily serves as an intermediary between the full counting jastrow orbital
// and the jastrow basis functions. Handles Counting Function storage, evaluation, and normalization.
//   Contains no optimizable parameters itself (...yet), but manages the parameters of the underlying
// basis functions by handing information (index maps) up to the counting jastrow orbital.
struct CountingRegionBase
{
  // number of electrons
  int num_els;
  // C.size()
  int num_regions;
  // counting function pointers
  std::vector<CountingFunctorBase*> C;
  // counting function id
  std::vector<std::string> C_id;

  // value arrays
  std::vector<RealType> _val;
  std::vector<RealType> _sum;
  std::vector<GradType> _grad;
  std::vector<RealType> _lap;

  // log values arrays
  std::vector<RealType> _Lval;
  std::vector<GradType> _Lgrad;
  std::vector<RealType> _Llap;

  // temporary value arrays
  std::vector<RealType> _val_t;
  std::vector<RealType> _sum_t;
  std::vector<GradType> _grad_t;
  std::vector<RealType> _lap_t;
  
  //memory for temporary log value arrays
  std::vector<RealType> _Lval_t;
  std::vector<GradType> _Lgrad_t;
  std::vector<RealType> _Llap_t;

  CountingRegionBase(int ne): num_els(ne) {}

  int size() const                  { return num_regions; }

  const opt_variables_type& getVars(int I) {return C[I]->myVars;}
  
  int max_num_derivs() const   
    { 
      auto comp = [](CountingFunctorBase* a, CountingFunctorBase* b){ return a->myVars.size() < b->myVars.size(); };
      CountingFunctorBase* Cmax = *(std::max_element(C.begin(),C.end(),comp));
      return Cmax->myVars.size();
    }

  // sum/grad/lap getters: index convention defined in base
  inline RealType& val(int I, int i)  { return _val[I*num_els + i]; }
  inline RealType& sum(int I)         { return _sum[I]; }
  inline GradType& grad(int I, int i) { return _grad[I*num_els + i]; }
  inline RealType& lap(int I, int i)  { return _lap[I*num_els + i]; }
  inline RealType& val_t(int I)       { return _val_t[I]; }
  inline RealType& sum_t(int I)       { return _sum_t[I]; }
  inline GradType& grad_t(int I)      { return _grad_t[I]; }
  inline RealType& lap_t(int I)       { return _lap_t[I]; }

  inline RealType& Lval(int I, int i)  { return _Lval[I*num_els + i]; }
  inline GradType& Lgrad(int I, int i) { return _Lgrad[I*num_els + i]; }
  inline RealType& Llap(int I, int i)  { return _Llap[I*num_els + i]; }

  inline RealType& Lval_t(int I)  { return _Lval_t[I]; }
  inline GradType& Lgrad_t(int I) { return _Lgrad_t[I]; }
  inline RealType& Llap_t(int I)  { return _Llap_t[I]; }

  virtual bool put(xmlNodePtr cur) = 0;
  
  virtual void initialize() 
  {
    app_log() << "CountingRegionBase::initialize" << std::endl;
    num_regions = C.size();
    // resize arrays
    _val.resize(num_regions*num_els);
    _sum.resize(num_regions);
    _grad.resize(num_regions*num_els);
    _lap.resize(num_regions*num_els);

    _val_t.resize(num_regions);
    _sum_t.resize(num_regions);
    _grad_t.resize(num_regions);
    _lap_t.resize(num_regions);

    _Lval.resize(num_regions*num_els);
    _Lgrad.resize(num_regions*num_els);
    _Llap.resize(num_regions*num_els);

    _Lval_t.resize(num_regions);
    _Lgrad_t.resize(num_regions);
    _Llap_t.resize(num_regions);
  }

  virtual void reportStatus(std::ostream& os) { }

  virtual CountingRegionBase* makeClone() = 0;

  virtual void restore(int iat) {}
  void checkInVariables(opt_variables_type& active)
  {
    for(auto it = C.begin(); it != C.end(); ++it)
      (*it)->checkInVariables(active);
  }
  void checkOutVariables(const opt_variables_type& active) 
  {
    for(auto it = C.begin(); it != C.end(); ++it)
      (*it)->checkOutVariables(active);
  }
  void resetParameters(const opt_variables_type& active)
  {
    for(auto it = C.begin(); it != C.end(); ++it)
      (*it)->resetParameters(active);
  }

  void acceptMove(ParticleSet& P, int iat)
  {
    for(int I = 0; I < num_regions; ++I)
    {
      sum(I) += val_t(I) - val(I,iat);
      val(I,iat) = val_t(I);
      grad(I,iat) = grad_t(I);
      lap(I,iat) = lap_t(I);
    }
  }
  virtual void evaluate(ParticleSet& P) = 0;
  virtual void evaluateTemp(ParticleSet& P, int iat) = 0;

  virtual void evaluate_print(std::ostream& os, ParticleSet& P)
  { 
    os << "CountingRegionBase::evaluate_print: not implemented" << std::endl;
  }
  virtual void evaluateTemp_print(std::ostream& os, ParticleSet& P) 
  {
    os << "CountingRegionBase::evaluateTemp_print: not implemented." << std::endl;
  }

  // using lambda functions for unambiguous multi-indexing
  virtual void evaluateTempDerivatives(ParticleSet& P, 
                                       int I, 
                                       int iat,
                                       std::function<RealType&(int,int)> dCsum) = 0;

  virtual void evaluateDerivatives(ParticleSet& P, 
                                   int I, 
                                   std::function<const GradType&(int,int)> FCgrad, 
                                   std::function<RealType&(int,int)> dCsum, 
                                   std::function<RealType&(int,int)> dCggsum,
                                   std::function<RealType&(int,int)> dClapsum, 
                                   std::vector<RealType>& dCFCggsum) = 0;

};

// Set of elementary counting functions
struct CountingRegion : public CountingRegionBase
{
  CountingRegion(int ne) : CountingRegionBase(ne) {}

  void reportStatus(std::ostream& os)
  {
    // print some class variables:
    os << "CountingRegion::reportStatus begin" << std::endl;
    os << "num_els: " << num_els << ", num_regions: " << num_regions << std::endl;
    os << "Counting Functions: " << std::endl;
    for(int I = 0; I < C.size(); ++I)
      C[I]->reportStatus(os);
    os << "CountingRegion::reportStatus end" << std::endl;
  }

  bool put(xmlNodePtr cur)
  {
    app_log() << "CountingRegion::put" << std::endl;
    // get the function type
    OhmmsAttributeSet rAttrib;
    std::string ftype = "none";
    rAttrib.add(ftype,"function");
    rAttrib.add(ftype,"ftype");
    rAttrib.put(cur);
    //std::transform(ftype.begin(),ftype.end(),ftype.begin(), static_cast<std::function<int(int)> >(std::tolower) );
    std::transform(ftype.begin(),ftype.end(),ftype.begin(), (int(*)(int))std::tolower );
    cur = cur->xmlChildrenNode;
    // name indices
    int i = 0;
    while(cur != NULL)
    {
      bool put_f = false;
      CountingFunctorBase* f;
      std::string cname((const char*)(cur->name));
      if(cname == "function")
      {
        // default id
        std::ostringstream s_id;
        s_id << boost::format("C_%d") %i;
        std::string f_id = s_id.str();
        // get function id
        OhmmsAttributeSet rAttrib2;
        rAttrib2.add(f_id,"f_id");
        rAttrib2.add(f_id,"id");
        rAttrib2.add(f_id,"name");
        rAttrib2.put(cur);
        // construct function object
        if(ftype == "gaussian")
        {
          f = new GaussianCountingFunctor(f_id);
          put_f = f->put(cur);
        }
        else if(ftype == "fermidirac" || ftype == "fermi-dirac")
        {
          f = new FermiDiracCountingFunctor(f_id);
          put_f = f->put(cur);
        }
        else
          APP_ABORT("CountingRegion::put:  unrecognized function name: \"" + ftype + "\"");
        ++i; // increment id index counter
        if(!put_f)
          APP_ABORT("CountingRegion::put: bad put: \"" + f_id + "\"");
        // add to counting region array
        C_id.push_back(f_id);
        C.push_back(f);
      }
      cur = cur->next;
    }
    CountingRegionBase::initialize();
    return true;
  }

  CountingRegionBase* makeClone()
  {
    app_log() << "  CountingRegion::makeClone" << std::endl;
    CountingRegion* cr = new CountingRegion(num_els);
    // copy class variables set in put()
    cr->C_id.resize( C_id.size() );
    cr->C.resize( C.size() );
    for(int i = 0; i < C.size(); ++i)
    {
      cr->C_id[i] = std::string(C_id[i]);
      cr->C[i] = C[i]->makeClone(C_id[i]);
    }
//    for(int i = 0; i < C.size(); ++i)
//    {
//      cr->C_id.push_back( std::string(C_id[i]) );
//      cr->C.push_back( C[i]->makeClone(C_id[i]) );
//    }
    // initialize 
    cr->initialize();
    return cr;
  } 

  void evaluate(ParticleSet& P) 
  {
    // clear arrays
    std::fill(_val.begin(),_val.end(),0);
    std::fill(_sum.begin(),_sum.end(),0);
    std::fill(_grad.begin(),_grad.end(),0);
    std::fill(_lap.begin(),_lap.end(),0);
    // temporary variables
    RealType Cval, Clap;
    GradType Cgrad;
    for(int I = 0; I < num_regions; ++I)
      for(int i = 0; i < num_els; ++i)
      {
        C[I]->evaluate(P.R[i],Cval,Cgrad,Clap);
        val(I,i)  = Cval;
        sum(I) += Cval;
        grad(I,i) = Cgrad;
        lap(I,i)  = Clap;
      }
  }


  void evaluateTemp(ParticleSet& P, int iat) 
  {
    // clear arrays
    std::fill(_val_t.begin(),_val_t.end(),0);
    std::fill(_sum_t.begin(),_sum_t.end(),0);
    std::fill(_grad_t.begin(),_grad_t.end(),0);
    std::fill(_lap_t.begin(),_lap_t.end(),0);
    // temporary variables
    RealType Cval_t, Clap_t;
    GradType Cgrad_t;
    for(int I = 0; I < num_regions; ++I)
    {
      C[I]->evaluate(P.R[iat], Cval_t, Cgrad_t, Clap_t);
      grad_t(I) = Cgrad_t;
      lap_t(I) = Clap_t;
      val_t(I) = Cval_t;
      sum_t(I) = sum(I) + val_t(I) - val(I,iat);
    }
  }

  // inputs:
  // - ParticleSet P: container that specifies particle positions (P.G) and total wavefunction gradient (P.G)
  // - int I: index of counting function
  // - index function FC grad ( FCgrad(I,i) provides a reference to (F\nabla_i C(r_i))_I 
  // 
  // note that x_{I_p} corresponds to the p^th optimizable parameter of the I^th function

  // outputs: 
  // - dCsum[p] = \sum\limits_i \frac{\partial C_I(r_i}}{\partial x_{I_p}} 
  // - dCggsum[p] = \sum\limits_i \nabla_i \ln(\Psi(r_i)) \cdot \nabla_i \frac{\partial C_I(r_i)}{\partial x_{I_p}}
  // - dClapsum[p] = \sum\limits_i \nabla_i^2 \frac{\partial C_I(r_i)}{\partial x_{I_p}} 
  // - dCFCggsum[p] = \sum\limits_i \nabla_i \frac{\partial C_I(r_i)}{\partial x_{I_p}} \cdot (F\nabla_i C_(r_i))
  // - active_index: pointer to vector: p = (active_index*)[i] ==> myVars[p] <=> active[i] 
  void evaluateDerivatives(ParticleSet& P, 
                          const int I, // index of the counting function parameter derivatives are associated with
                          std::function<const GradType&(int,int)> FCgrad, 
                          std::function<RealType&(int,int)> dCsum, 
                          std::function<RealType&(int,int)> dCggsum,
                          std::function<RealType&(int,int)> dClapsum, 
                          std::vector<RealType>& dCFCggsum)
  {
    static std::vector<RealType> dCval;
    static std::vector<GradType> dCgrad;
    static std::vector<RealType> dClap;
    static int mnd = max_num_derivs();

    dCval.resize(mnd);
    dCgrad.resize(mnd);
    dClap.resize(mnd);
    for(int i = 0; i < num_els; ++i)
    {
      // derivatives of C[J] wrt pI are zero: only need C[I] derivatives
      C[I]->evaluateDerivatives(P.R[i],dCval,dCgrad,dClap);
      int num_derivs = getVars(I).size(); 
      for(int p = 0; p < num_derivs; ++p)
      {
        dCsum(I,p) += dCval[p];
        dCggsum(I,p) += dot(dCgrad[p],P.G[i]);
        dClapsum(I,p) += dClap[p];
        dCFCggsum[p] += dot(dCgrad[p],FCgrad(I,i));
      }
    }
  }

  void evaluateTempDerivatives(ParticleSet& P, 
                               const int I, // index of the counting function parameter derivatives are associated with
                               int iat,
                               std::function<RealType&(int,int)> dCsum)
  {
  }


};

// Set of elementary counting functions, normalized
struct NormalizedCountingRegion : public CountingRegion
{

  // normalization value, by counting function: 
  // Nval[i] = \sum\limits_{I} Cval(I,i)
  std::vector<RealType> Nval;
  RealType Nval_t;

  // maximum of log of unnormalized counting functions, by electron
  std::vector<RealType> Lmax; 
  RealType Lmax_t;

  std::vector<RealType> _dLval_saved;

  // reference orbital id
  std::string Cref_id;
  // pointer to a copy of the reference orbital
  CountingFunctorBase* Cref;

  // data saved for single-particle-move ratio deriv evaluation
  inline RealType& dLval_saved(int I, int p, int i)  { return _dLval_saved[I*max_num_derivs()*num_els + p*num_els + i]; }

  NormalizedCountingRegion(int ne) : CountingRegion(ne) {}

  bool put(xmlNodePtr cur)
  {
    // do CountingRegion put
    CountingRegion::put(cur);
    // get the reference function
    OhmmsAttributeSet rAttrib;
    Cref_id = "none"; // class variable
    rAttrib.add(Cref_id,"ref_id");
    rAttrib.add(Cref_id,"reference_id");
    rAttrib.add(Cref_id,"Cref_id");
    rAttrib.put(cur);
    // loop through array, find where Cref is
    auto C_id_it = std::find(C_id.begin(),C_id.end(),Cref_id);
    // get index of the reference
    int ref_index = std::distance(C_id.begin(),C_id_it);
    if(Cref_id == "none" || C_id_it == C_id.end())
      APP_ABORT("NormalizedCountingRegion::put: reference function not found:"+ (Cref_id == "none"?" Cref not specified":"\"" + Cref_id + "\"")); 
    // make a copy of the reference gaussian
    Cref = C[ref_index]->makeClone(Cref_id + "_ref");
    // divide all gaussians by the reference
    for(auto it = C.begin(); it != C.end(); ++it)
    {
      (*it)->divide_eq(Cref);
    }
    initialize();
    return true;
  }

  CountingRegionBase* makeClone()
  {
    app_log() << "  NormalizedCountingRegion::makeClone" << std::endl;
    NormalizedCountingRegion* ncr = new NormalizedCountingRegion(num_els);
    // copy class variables set in put()
    for(int i = 0; i < C.size(); ++i)
    {
      ncr->C_id.push_back( std::string(C_id[i]) );
      ncr->C.push_back( C[i]->makeClone(C_id[i]) );
    }
    ncr->Cref_id = std::string(Cref_id);
    ncr->Cref = Cref->makeClone(Cref_id + "_ref");
    // initialize 
    ncr->initialize();
    return ncr;
  }
  
  void initialize()
  {
    app_log() << "NormalizedCountingRegion::initialize" << std::endl;
    CountingRegionBase::initialize();
    // allocate memory for normalization 
    Nval.resize(num_els);
    Lmax.resize(num_els);
    // store log derivative values for single particle moves
    _dLval_saved.resize(max_num_derivs()*num_regions*num_els);
  }

  void reportStatus(std::ostream& os)
  {
    os << "NormalizedCountingRegion::reportStatus begin" << std::endl;
    // print some class variables:
    os << "num_els: " << num_els << ", num_regions: " << num_regions << std::endl;
    os << "Counting Functions: " << std::endl;
    Cref->reportStatus(os);
    for(int I = 0; I < C.size(); ++I)
      C[I]->reportStatus(os);
    os << "NormalizedCountingRegion::reportStatus end" << std::endl;
  }

  // evaluate using the log of the counting basis
  void evaluate(ParticleSet& P)
  {
    // clear arrays
    std::fill(_val.begin(),_val.end(),0);
    std::fill(_sum.begin(),_sum.end(),0);
    std::fill(_grad.begin(),_grad.end(),0);
    std::fill(_lap.begin(),_lap.end(),0);
    std::fill(_Lval.begin(),_Lval.end(),0);
    std::fill(_Lgrad.begin(),_Lgrad.end(),0);
    std::fill(_Llap.begin(),_Llap.end(),0);
    std::fill(Nval.begin(),Nval.end(),0);
    std::fill(Lmax.begin(),Lmax.end(),0);
    // temporary variables: Lval = ln(C), Lgrad = \nabla ln(C), Llap = \nabla^2 ln(C)
    for(int i = 0; i < num_els; ++i)
    {
      for(int I = 0; I < num_regions; ++I)
      {
        C[I]->evaluateLog(P.R[i],Lval(I,i),Lgrad(I,i),Llap(I,i));
        if(Lval(I,i) > Lmax[i])
          Lmax[i] = Lval(I,i);

      }
      // build counting function values; subtract off largest log value
      for(int I = 0; I < num_regions; ++I)
      {
        val(I,i) = std::exp(Lval(I,i) - Lmax[i]);
        Nval[i] += val(I,i);
      }
      GradType gLN_sum = 0; // \sum\limits_I \nabla L_{Ii} N_{Ii}
      RealType lLN_sum = 0; // \sum\limits_I \nabla^2 L_{Ii} N_{Ii}
      // build normalized counting function value, intermediate values for gradient
      for(int I = 0; I < num_regions; ++I)
      {
        val(I,i) = val(I,i) / Nval[i];
        sum(I) += val(I,i);
        gLN_sum += Lgrad(I,i) * val(I,i);
        lLN_sum += Llap(I,i) * val(I,i);
      }
      RealType gLgN_sum = 0; // \sum\limits_{I} \nabla L_{Ii} \cdot \nabla N_{Ii}
      // build gradient, intermediate values for laplacian
      for(int I = 0; I < num_regions; ++I)
      {
        grad(I,i) = (Lgrad(I,i) - gLN_sum)*val(I,i);
        gLgN_sum += dot(Lgrad(I,i), grad(I,i));
      }
      //build laplacian
      for(int I = 0; I < num_regions; ++I)
      {
        lap(I,i) = (Llap(I,i) - lLN_sum - gLgN_sum)*val(I,i) + dot(grad(I,i),Lgrad(I,i) - gLN_sum);
      }
    }
  }

  void evaluate_print(std::ostream& os, ParticleSet& P)
  {
    for(auto it = C.begin(); it != C.end(); ++it)
      (*it)->evaluate_print(os,P);
    os << "NormalizedCountingRegions::evaluate_print" << std::endl;
    os << "val: ";
    std::copy(_val.begin(),_val.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "sum: ";
    std::copy(_sum.begin(),_sum.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "grad: ";
    std::copy(_grad.begin(),_grad.end(),std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "lap: ";
    std::copy(_lap.begin(),_lap.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "Nval: ";
    std::copy(Nval.begin(),Nval.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "Lmax: "; 
    std::copy(Lmax.begin(), Lmax.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl;
  }


  void evaluateTemp(ParticleSet& P, int iat)
  {
    // clear arrays
    std::fill(_val_t.begin(),_val_t.end(),0);
    std::fill(_sum_t.begin(),_sum_t.end(),0);
    std::fill(_grad_t.begin(),_grad_t.end(),0);
    std::fill(_lap_t.begin(),_lap_t.end(),0);
    std::fill(_Lval_t.begin(),_Lval_t.end(),0);
    std::fill(_Lgrad_t.begin(),_Lgrad_t.end(),0);
    std::fill(_Llap_t.begin(),_Llap_t.end(),0);

    Lmax_t = Lmax[iat];
    Nval_t = 0;
    // temporary variables
    for(int I = 0; I < num_regions; ++I)
    {
      C[I]->evaluateLog(P.R[iat],Lval_t(I),Lgrad_t(I),Llap_t(I));
      if(Lval_t(I) > Lmax_t)
        Lmax_t = Lval_t(I);
    }
    // build counting function values; subtract off largest log value
    for(int I = 0; I < num_regions; ++I)
    {
      val_t(I) = std::exp(Lval_t(I) - Lmax_t);
      Nval_t += val_t(I);
    }
    GradType gLN_sum_t = 0; // \sum\limits_I \nabla L_{Ii} N_{Ii}
    RealType lLN_sum_t = 0; // \sum\limits_I \nabla^2 L_{Ii} N_{Ii}
    // build normalized counting function value, intermediate values for gradient
    for(int I = 0; I < num_regions; ++I)
    {
      val_t(I) = val_t(I) / Nval_t;
      sum_t(I) = sum(I) + val_t(I) - val(I,iat);
      gLN_sum_t += Lgrad_t(I) * val_t(I);
      lLN_sum_t += Llap_t(I) * val_t(I);
    }
    RealType gLgN_sum_t = 0; // \sum\limits_{I} \nabla L_{Ii} \cdot \nabla N_{Ii}
    // build gradient, intermediate values for laplacian
    for(int I = 0; I < num_regions; ++I)
    {
      grad_t(I) = (Lgrad_t(I) - gLN_sum_t)*val_t(I);
      gLgN_sum_t += dot(Lgrad_t(I), grad_t(I));
    }
    //build laplacian
    for(int I = 0; I < num_regions; ++I)
    {
      lap_t(I) = (Llap_t(I) - lLN_sum_t - gLgN_sum_t)*val_t(I) + dot(grad_t(I), Lgrad_t(I) - gLN_sum_t);
    }
    
  }

  void evaluateTemp_print(std::ostream& os, ParticleSet& P)
  {
    for(auto it = C.begin(); it != C.end(); ++it)
      (*it)->evaluate_print(os,P);
    os << "NormalizedCountingRegions::evaluateTemp_print" << std::endl;
    os << "val_t: ";
    std::copy(_val_t.begin(),_val_t.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "sum_t: ";
    std::copy(_sum_t.begin(),_sum_t.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "grad_t: ";
    std::copy(_grad_t.begin(),_grad_t.end(),std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "lap_t: ";
    std::copy(_lap_t.begin(),_lap_t.end(),std::ostream_iterator<RealType>(os,", "));
    os << std::endl << "Nval_t: " << Nval_t;
    os << std::endl << "Lmax_t: " << Lmax_t;
    os << std::endl;
  }


  void acceptMove(ParticleSet& P, int iat)
  {
    Nval[iat] = Nval_t;
    Lmax[iat] = Lmax_t;
    for(int I = 0; I < num_regions; ++I)
    {
      sum(I) = sum_t(I); 
      val(I,iat) = val_t(I);
      grad(I,iat) = grad_t(I);
      lap(I,iat) = lap_t(I);

      Lval(I,iat) = Lval_t(I);
      Lgrad(I,iat) = Lgrad_t(I);
      Llap(I,iat) = Llap_t(I);
    }
  }

  // calculates derivatives of single particle move of particle with index iat
  void evaluateTempDerivatives(ParticleSet& P, 
                               const int I, // index of the counting function parameter derivatives are associated with
                               int iat,
                               std::function<RealType&(int,int)> dNdiff)
  {
    // may assume evaluate and evaluateTemp has already been called
    int num_derivs = getVars(I).size();
    // get log derivatives
    static std::vector<RealType> dLval_t;
    static int mnd = max_num_derivs();
    dLval_t.resize(mnd);

    C[I]->evaluateLogTempDerivatives(P.R[iat], dLval_t);
    for(int J = 0; J < num_regions; ++J)
    {
      for(int p = 0; p < num_derivs; ++p)
      {
        RealType val_Ii = (I==J) - val(I,iat);
        RealType val_It = (I==J) - val_t(I);
        dNdiff(J,p) = val_t(J)*dLval_t[p]*val_It - val(J,iat)*dLval_saved(I,p,iat)*val_Ii;
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P, 
                          int I, 
                          std::function<const GradType&(int,int)> FCgrad, 
                          std::function<RealType&(int,int)> dNsum, 
                          std::function<RealType&(int,int)> dNggsum,
                          std::function<RealType&(int,int)> dNlapsum, 
                          std::vector<RealType>& dNFNggsum)
  {
    evaluate(P);
    static std::vector<RealType> dLval;
    static std::vector<GradType> dLgrad;
    static std::vector<RealType> dLlap;
    static int mnd = max_num_derivs();
    dLval.resize(mnd);
    dLgrad.resize(mnd);
    dLlap.resize(mnd);

    int num_derivs = getVars(I).size();
    for(int i = 0; i < num_els; ++i)
    {
      // get log derivatives
      //C[I]->evaluateLogDerivatives(P.R[i], dCval, dLgrad, dLlap);
      C[I]->evaluateLogDerivatives(P.R[i], dLval, dLgrad, dLlap);
      for(int J = 0; J < num_regions; ++J)
      {
        for(int p = 0; p < num_derivs; ++p)
        {
          RealType val_Ii = (I==J) - val(I,i);

          RealType dNval = val(J,i)*dLval[p]*val_Ii;
          GradType dNgrad = grad(J,i)*dLval[p]*val_Ii + val(J,i)*dLgrad[p]*val_Ii  - val(J,i)*dLval[p]*grad(I,i);

          RealType dNlap = lap(J,i)*dLval[p]*val_Ii + 2*dot(grad(J,i),dLgrad[p])*val_Ii - 2*dot(grad(J,i),grad(I,i))*dLval[p] 
                                                      + val(J,i)*dLlap[p]*val_Ii          - 2*val(J,i)*dot(dLgrad[p],grad(I,i)) 
                                                                                          - val(J,i)*dLval[p]*lap(I,i);         
          // accumulate
          dLval_saved(I,p,i) = dLval[p];
          dNsum(J,p) += dNval;
          dNggsum(J,p) += dot(dNgrad, P.G[i]);
          dNlapsum(J,p) += dNlap;
          dNFNggsum[p] += dot(dNgrad, FCgrad(J,i));
        }
      }
    }
  } // end evaluateDerivatives

};

}


#endif
