#ifndef QMCPLUSPLUS_COUNTINGFUNCTORS_H
#define QMCPLUSPLUS_COUNTINGFUNCTORS_H

#include "QMCWaveFunctions/OrbitalBase.h"
#include "Optimize/VariableSet.h"
#include "OhmmsData/AttributeSet.h"
#include "boost/format.hpp"

#include <array>
#include <formic/utils/lapack_interface.h>
#include "Numerics/Blasf.h"

namespace qmcplusplus
{


typedef OrbitalBase::ValueType ValueType;
typedef OrbitalBase::PosType PosType;
typedef OrbitalBase::GradType GradType;
typedef OrbitalBase::RealType RealType;
typedef OrbitalBase::TensorType TensorType;
typedef optimize::VariableSet::real_type real_type;
typedef optimize::VariableSet opt_variables_type;

// base class for optimizable counting functors
struct CountingFunctorBase
{
 
  // container for optimizable variables
  opt_variables_type myVars;

  // identifier string (mostly for labelling opimizable variables)
  std::string id;

  // most recent evaluations
  RealType Fval;
  GradType Fgrad;
  RealType Flap; 

  inline CountingFunctorBase(std::string f_id) : id(f_id) {}
  virtual bool put(xmlNodePtr cur) = 0;
  virtual void reportStatus(std::ostream& os) = 0;

  // parameter handling
  void checkInVariables(opt_variables_type& active) { active.insertFrom(myVars); }
  void checkOutVariables(const opt_variables_type& active) { myVars.getIndex(active); }
  virtual void resetParameters(const opt_variables_type& active) = 0;
  // return a pointer to the myVars index array by reference
  // return the name of the i^th parameter
  const std::vector<std::pair<std::string,real_type> >* getNameAndValue() { return &myVars.NameAndValue; }

  virtual void initialize() = 0;

  virtual void multiply_eq(const CountingFunctorBase* rhs) = 0;
  virtual void divide_eq(const CountingFunctorBase* rhs) = 0;

  virtual CountingFunctorBase* makeClone(std::string clone_id) = 0;

  virtual void evaluate(PosType r, RealType& fval, GradType& fgrad, RealType& flap) = 0;
  virtual void evaluateLog(PosType r, RealType& lval, GradType& lgrad, RealType& llap) = 0;
  virtual void evaluateDerivatives(PosType r, std::vector<RealType>& dfval, std::vector<GradType>& dfgrad, std::vector<RealType>& dflap) = 0;
  virtual void evaluateLogDerivatives(PosType r, std::vector<RealType>& dlval, std::vector<GradType>& dlgrad, std::vector<RealType>& dllap) = 0;
  virtual void evaluateLogTempDerivatives(PosType r, std::vector<RealType>& dlval) = 0;

  // debug helper functions
  virtual void evaluate_print(std::ostream& os, ParticleSet& P)
  { 
    os << "CountingFunctorBase::evaluate_print: not implemented." << std::endl;
  }

};

//struct NormalizableBase
//{
//  // functions used for normalization
//  void multiply_eq(const CountingFunctorBase* rhs);
//  void divide_eq(const CountingFunctorBase* rhs);
//};

// optimizable three-dimensional gaussian
struct GaussianCountingFunctor : public CountingFunctorBase
{
  // three-dimensional gaussian of the form: exp(x^T A x -2*B^t x + C)
  //   designed for correspondence with the factored form:
  //   (x - mu)^T A (x - mu) + K: B <=> Tmu
  //   The factored form is nice because we can easily think of each of the parameters as
  //   ellipsoid-point-weight (T,mu,K, respectively)
  // in the case of planar functions, were we to optimize the factored form directly
  //   we would find that mu shoots off to infinity as T gets really small if a planar form is 
  //   optimal which is pathological behavior.
  // We take either either factored or unfactored form as inputs and work with the factored
  //   form during optimization.

  // enumerations for axes and T parameters
  enum A_vars { XX, XY, XZ, YY, YZ, ZZ, NUM_A };
  enum B_vars { X, Y, Z, DIM};
  // working parameter variables
  // quadratic coefficients
  TensorType A;
  // linear coefficients
  PosType B;
  // constant coeff
  RealType C;

  // opt variables: vector of bools: one for each parameter
  std::vector<bool> opt_A;
  std::vector<bool> opt_B;
  bool opt_C;

  // constructor
  GaussianCountingFunctor(std::string f_id) : CountingFunctorBase(f_id) {}

  void reportStatus(std::ostream& os)
  {
    os << "GaussianCountingFunctor::reportStatus begin" << std::endl;
    os << "id: " << id << std::endl;
    os << "  A: ";
    for(int I = 0; I < 3; ++I)
    {
      os << std::endl;
      for(int J = 0, IJ = I*3; J < 3; ++J, ++IJ)
        os << boost::format("  %10.5f") % A[IJ];
    }
    os << std::endl << "  opt_A: ";
    std::copy(opt_A.begin(), opt_A.end(), std::ostream_iterator<bool>(os,", "));
    os << std::endl << "  B: "; 
    for(auto it = B.begin(); it != B.end(); ++it)
      os << boost::format("  %10.5f") % (*it);
    os << std::endl << "  opt_B: "; 
    std::copy(opt_B.begin(), opt_B.end(), std::ostream_iterator<bool>(os,", "));
    os << std::endl << "  C: " << C << std::endl;
    os << "  opt_C: " << opt_C << std::endl; 
    os << "  registered optimizable variables:" << std::endl;
    myVars.print(os);
    os << "GaussianCountingFunctor::reportStatus end" << std::endl;
  }

  // converts eulerform to matrix
  void euler_to_matrix(const std::array<RealType,6>& euler, TensorType& matrix)
  {
    enum euler_vars { SIGMA_X, SIGMA_Y, SIGMA_Z, THETA_X, THETA_Y, THETA_Z, NUM_EULER};
    // trig functions of euler angles
    RealType cx = std::cos(euler[THETA_X]);
    RealType sx = std::sin(euler[THETA_X]);
    RealType cy = std::cos(euler[THETA_Y]);
    RealType sy = std::sin(euler[THETA_Y]);
    RealType cz = std::cos(euler[THETA_Z]);
    RealType sz = std::sin(euler[THETA_Z]);
    // euler rotation matrices
    TensorType Rx(1,0,0, 0,cx,sx, 0,-sx,cx);
    TensorType Ry(cy,0,-sy, 0,1,0, sy,0,cy);
    TensorType Rz(cz,sz,0, -sz,cz,0, 0,0,1);
    // build matrix of eigenvalues
    RealType dx = 1/(euler[SIGMA_X]*euler[SIGMA_X]);
    RealType dy = 1/(euler[SIGMA_Y]*euler[SIGMA_Y]);
    RealType dz = 1/(euler[SIGMA_Z]*euler[SIGMA_Z]);
    TensorType Sigma(dx,0,0, 0,dy,0, 0,0,dz); 
    // order of multiplication matters
    TensorType Ut = dot(dot(Rz,Ry),Rx);
    TensorType U = transpose(Ut);
    matrix = dot(dot(U,Sigma),Ut);
  }

  // builds the full matrix from its upper triangle
  void triang_to_matrix(const std::array<RealType,6>& triang, TensorType& matrix)
  {
    auto it = triang.begin();
    for(int i = 0; i < 3; ++i)
      for(int j = i; j < 3; ++j)
      {
        matrix(i,j) = matrix(j,i) = *it;
        ++it;
      }
  }

  void d_to_b(const PosType& d, const TensorType& a, PosType& b)
  {
    b = dot(a,d); 
  }

  void k_to_c(const RealType& k, const TensorType& a, const PosType& d, RealType& c)
  {
    c = k + dot(d,dot(a,d));
  }

  bool put(xmlNodePtr cur)
  {
    app_log() << "GaussianCountingFunctor::put" << std::endl;
    bool put_A = false, put_B = false, put_C = false, put_D = false, put_K = false;
    // alternate inputs
    std::array<RealType,6> A_euler;
    std::array<RealType,6> A_triang;
    PosType D = 0;
    RealType K = 0;
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "var")
      {
        std::string opt = "false";
        std::string name = "none";
        OhmmsAttributeSet rAttrib;
        rAttrib.add(name,"name");
        rAttrib.add(opt,"opt");
        rAttrib.add(opt,"opt_bits");
        rAttrib.put(cur);
        // check if opt is a binary string
        bool is_bitstr = std::all_of(opt.begin(), opt.end(), [&](char c){ return c == '0' || c == '1';} );
        std::transform(opt.begin(),opt.end(),opt.begin(), (int (*)(int)) std::tolower);
        std::transform(name.begin(),name.end(),name.begin(), (int (*)(int)) std::tolower);
        if(name == "aeuler") // input is the eulerform of A
        {
          // copy opt bits to class 
          opt_A.resize(6);
          if(opt.size() == 6 && is_bitstr) 
            std::transform(opt.begin(),opt.end(),opt_A.begin(),[&](char c){return (c == '1');} );
          else
            // default opt = true
            std::fill(opt_A.begin(),opt_A.end(), (opt == "true") );
          put_A = putContent(A_euler.begin(),A_euler.end(),cur);
          // perform conversion
          euler_to_matrix(A_euler,A);
        }
        if(name == "a" || name == "atriang") // input is the upper-triangle of A
        {
          // bits still correspond to sigma, theta optimization 
          opt_A.resize(6);
          if(opt.size() == 6 && is_bitstr)
            std::transform(opt.begin(),opt.end(),opt_A.begin(),[&](char c){return (c == '1');} );
          // default opt = true
          else
            std::fill(opt_A.begin(),opt_A.end(), (opt == "true") );
          put_A = putContent(A_triang.begin(),A_triang.end(),cur);
          // perform conversion
          triang_to_matrix(A_triang,A);
        }
        if(name == "aelem") // input is A, by element
        {
          opt_A.resize(6);
          if(opt.size() == 6 && is_bitstr)
            std::transform(opt.begin(),opt.end(),opt_A.begin(),[&](char c){return (c == '1');} );
          // default opt = true
          else
            std::fill(opt_A.begin(),opt_A.end(), (opt == "true") );
          put_A = putContent(A.begin(),A.end(),cur);
        }
        if(name == "b")
        {
          opt_B.resize(3);
          if(opt.size() == 3 && is_bitstr)
            std::transform(opt.begin(),opt.end(),opt_B.begin(),[&](char c){return (c == '1');} );
          // default opt = true
          else
            std::fill(opt_B.begin(),opt_B.end(), (opt == "true") );
          put_B = putContent(B.begin(),B.end(),cur);
        }
        if(name == "d")
        {
          opt_B.resize(3);
          if(opt.size() == 3 && is_bitstr)
            std::transform(opt.begin(),opt.end(),opt_B.begin(),[&](char c){return (c == '1');} );
          else
            std::fill(opt_B.begin(),opt_B.end(), (opt == "true") );
          put_D = putContent(D.begin(),D.end(),cur);
        }
        if(name == "c")
        {
          opt_C = (opt == "true");
          put_C = putContent(C,cur);
        }
        if(name == "k")
        {
          opt_C = (opt == "true");
          put_K = putContent(K,cur);
        }
      }
      cur = cur->next;
    }
    // convert B <=> D 
    if(put_B && put_D)
      APP_ABORT("GaussianCountingFunctor::put: overdetermined: both B and D are specified");
    if(put_C && put_K)
      APP_ABORT("GaussianCountingFunctor::put: overdetermined: both C and K are specified");
    // convert D,B using A
    if(put_D)
      d_to_b(D,A,B);
    if(put_K)
      k_to_c(K,A,D,C);
    // check that we got everything
    if(!(put_A && (put_B || put_D) && (put_C || put_K)) )
    {
      std::ostringstream err;
      err << "GaussianCountingFunctor::put:  required variable unspecified: ";
      err << "put_A: " << put_A << ", put_B: " << put_B << ", put_C: " << put_C << ", put_D: " << put_D << ", put_K: " << put_K << std::endl;
      APP_ABORT(err.str());
    }
    initialize();
    return true;
  }

  void initialize()
  {
    app_log() << "GaussianCountingFunctor::initialize: " << id << std::endl;
    // register and update optimizable variables to current values
    if(opt_A[XX]) { myVars.insert(id+"_A_xx", A(X,X), opt_A[XX], optimize::OTHER_P); myVars[id+"_A_xx"] = A(X,X); }
    if(opt_A[XY]) { myVars.insert(id+"_A_xy", A(X,Y), opt_A[XY], optimize::OTHER_P); myVars[id+"_A_xy"] = A(X,Y); }
    if(opt_A[XZ]) { myVars.insert(id+"_A_xz", A(X,Z), opt_A[XZ], optimize::OTHER_P); myVars[id+"_A_xz"] = A(X,Z); }
    if(opt_A[YY]) { myVars.insert(id+"_A_yy", A(Y,Y), opt_A[YY], optimize::OTHER_P); myVars[id+"_A_yy"] = A(Y,Y); }
    if(opt_A[YZ]) { myVars.insert(id+"_A_yz", A(Y,Z), opt_A[YZ], optimize::OTHER_P); myVars[id+"_A_yz"] = A(Y,Z); }
    if(opt_A[ZZ]) { myVars.insert(id+"_A_zz", A(Z,Z), opt_A[ZZ], optimize::OTHER_P); myVars[id+"_A_zz"] = A(Z,Z); }
    if(opt_B[X])  { myVars.insert(id+"_B_x", B[X], opt_B[X], optimize::OTHER_P);     myVars[id+"_B_x"] = B[X]; }
    if(opt_B[Y])  { myVars.insert(id+"_B_y", B[Y], opt_B[Y], optimize::OTHER_P);     myVars[id+"_B_y"] = B[Y]; } 
    if(opt_B[Z])  { myVars.insert(id+"_B_z", B[Z], opt_B[Z], optimize::OTHER_P);     myVars[id+"_B_z"] = B[Z]; }
    if(opt_C)     { myVars.insert(id+"_C", C, opt_C, optimize::OTHER_P);             myVars[id+"_C"] = C; }
  }


  void resetParameters(const opt_variables_type& active)
  {
    // set myVars from active
    for(int i = 0; i < myVars.size(); ++i)
    {
      int ia = myVars.where(i);
      if(ia != -1)
        myVars[i] = active[ia];
    }
    app_log() << "GaussianCountingFunctor::resetParameters" << std::endl;
    // set local variables from myVars
    if(opt_A[XX]) A(0,0) = myVars[id+"_A_xx"];
    if(opt_A[YY]) A(1,1) = myVars[id+"_A_yy"];
    if(opt_A[ZZ]) A(2,2) = myVars[id+"_A_zz"];
    if(opt_A[XY]) A(1,0) = A(0,1) = myVars[id+"_A_xy"];
    if(opt_A[XZ]) A(2,0) = A(0,2) = myVars[id+"_A_xz"];
    if(opt_A[YZ]) A(2,1) = A(1,2) = myVars[id+"_A_yz"];
    if(opt_B[X]) B[X] = myVars[id+"_B_x"];
    if(opt_B[Y]) B[Y] = myVars[id+"_B_y"];
    if(opt_B[Z]) B[Z] = myVars[id+"_B_z"];
    if(opt_C) C = myVars[id+"_C"];
  }

  // f = std::exp(x^T A x - 2B^T x + C)
  void evaluate(PosType r, RealType& fval, GradType& fgrad, RealType& flap)
  {
    PosType Ar = dot(A,r);
    RealType x = dot(Ar-2*B,r) + C;
    fval = std::exp(x);
    fgrad = 2*(Ar-B)*fval;
    flap = 4*dot(Ar-B,Ar-B)*fval + 2*trace(A)*fval;
  }

  // f = x^T A x - 2B^T x + C
  void evaluateLog(PosType r, RealType& lval, GradType& lgrad, RealType& llap)
  {
    PosType Ar = dot(A,r);
    lval = dot(Ar-2*B,r) + C;
    lgrad = 2*(Ar-B);
    llap = 2*trace(A);
  }

  // evaluate parameter derivatives of log of the wavefunction
  void evaluateLogDerivatives(PosType r,
                              std::vector<RealType>& dlval,
                              std::vector<GradType>& dlgrad,
                              std::vector<RealType>& dllap)
  {
    int p = 0;
    if(opt_A[XX]) { dlval[p] =   r[X]*r[X]; dlgrad[p] = 2*GradType(r[X],0,0);    dllap[p] = 2; ++p; }
    if(opt_A[XY]) { dlval[p] = 2*r[X]*r[Y]; dlgrad[p] = 2*GradType(r[Y],r[X],0); dllap[p] = 0; ++p; }
    if(opt_A[XZ]) { dlval[p] = 2*r[X]*r[Z]; dlgrad[p] = 2*GradType(r[Z],0,r[X]); dllap[p] = 0; ++p; }
    if(opt_A[YY]) { dlval[p] =   r[Y]*r[Y]; dlgrad[p] = 2*GradType(0,r[Y],0);    dllap[p] = 2; ++p; }
    if(opt_A[YZ]) { dlval[p] = 2*r[Y]*r[Z]; dlgrad[p] = 2*GradType(0,r[Z],r[Y]); dllap[p] = 0; ++p; }
    if(opt_A[ZZ]) { dlval[p] =   r[Z]*r[Z]; dlgrad[p] = 2*GradType(0,0,r[Z]);    dllap[p] = 2; ++p; }
    if(opt_B[X])  { dlval[p] = -2*r[X];     dlgrad[p] = -2*GradType(1,0,0);      dllap[p] = 0; ++p; }
    if(opt_B[Y])  { dlval[p] = -2*r[Y];     dlgrad[p] = -2*GradType(0,1,0);      dllap[p] = 0; ++p; }
    if(opt_B[Z])  { dlval[p] = -2*r[Z];     dlgrad[p] = -2*GradType(0,0,1);      dllap[p] = 0; ++p; }
    if(opt_C)     { dlval[p] = 1;           dlgrad[p] = 0;                       dllap[p] = 0; ++p; }
  }

  void evaluateLogTempDerivatives(PosType r,
                                   std::vector<RealType>& dlval)
  {
    int p = 0;
    if(opt_A[XX]) { dlval[p] =   r[X]*r[X]; ++p; }
    if(opt_A[XY]) { dlval[p] = 2*r[X]*r[Y]; ++p; }
    if(opt_A[XZ]) { dlval[p] = 2*r[X]*r[Z]; ++p; }
    if(opt_A[YY]) { dlval[p] =   r[Y]*r[Y]; ++p; }
    if(opt_A[YZ]) { dlval[p] = 2*r[Y]*r[Z]; ++p; }
    if(opt_A[ZZ]) { dlval[p] =   r[Z]*r[Z]; ++p; }
    if(opt_B[X])  { dlval[p] = -2*r[X];     ++p; }
    if(opt_B[Y])  { dlval[p] = -2*r[Y];     ++p; }
    if(opt_B[Z])  { dlval[p] = -2*r[Z];     ++p; }
    if(opt_C)     { dlval[p] = 1;           ++p; }
  }

  void evaluate_print(std::ostream& os, ParticleSet& P)
  {
    // calculate all intermediates and values for electrons in a particle set
    std::vector<PosType> r_vec;
    std::vector<PosType> Ar_vec;
    std::vector<RealType> x_vec;
    std::vector<RealType> fval_vec;
    std::vector<GradType> fgrad_vec;
    std::vector<RealType> flap_vec;
    RealType x, fval, flap;
    PosType Ar, r;
    GradType fgrad;
    for(auto it = P.R.begin(); it != P.R.end(); ++it)
    {
      r = *it;
      Ar = dot(A,r);
      x = dot(Ar-2*B,r) + C;
      fval = std::exp(x);
      fgrad = 2*(Ar-B)*fval;
      flap = 4*dot(Ar-B,Ar-B)*fval + 2*trace(A)*fval;

      r_vec.push_back(r);
      Ar_vec.push_back(Ar);
      x_vec.push_back(x);
      fval_vec.push_back(fval);
      fgrad_vec.push_back(fgrad);
      flap_vec.push_back(flap);
    }
    os << "CountingFunctor::evaluate_print: id: " << id << std::endl;
    os << "r: ";
    std::copy(r_vec.begin(),r_vec.end(), std::ostream_iterator<PosType>(os,", "));
    os << std::endl <<  "Ar: ";
    std::copy(Ar_vec.begin(),Ar_vec.end(), std::ostream_iterator<PosType>(os,", "));
    os << std::endl <<  "x: ";
    std::copy(x_vec.begin(),x_vec.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl <<  "fval: ";
    std::copy(fval_vec.begin(),fval_vec.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl <<  "fgrad: ";
    std::copy(fgrad_vec.begin(),fgrad_vec.end(), std::ostream_iterator<GradType>(os,", "));
    os << std::endl << "flap: ";
    std::copy(flap_vec.begin(),flap_vec.end(), std::ostream_iterator<RealType>(os,", "));
    os << std::endl;
  }

  void evaluateDerivative_A(A_vars q, PosType r, RealType& dfval, GradType& dfgrad, RealType& dflap)
  {
    // Fval, Fgrad, Flap are up-to-date function values
    // x correponds to the exponent value: x = ln f; dx are param derivs
    RealType dxval = 0, dxlap = 0;
    GradType dxgrad = 0;
    if(q == XX) { dxval = r[X]*r[X]; dxgrad = GradType(2*r[X],0,0);  dxlap = 2; }
    if(q == YY) { dxval = r[Y]*r[Y]; dxgrad = GradType(0,2*r[Y],0);  dxlap = 2; }
    if(q == ZZ) { dxval = r[Z]*r[Z]; dxgrad = GradType(0,0,2*r[Z]);  dxlap = 2; }
    // off-diagonal terms: factor of two is since A is symmetric
    if(q == XY) { dxval = 2*r[X]*r[Y]; dxgrad = 2*GradType(r[Y],r[X],0); dxlap = 0; }
    if(q == XZ) { dxval = 2*r[X]*r[Z]; dxgrad = 2*GradType(r[Z],0,r[X]); dxlap = 0; }
    if(q == YZ) { dxval = 2*r[Y]*r[Z]; dxgrad = 2*GradType(0,r[Z],r[Y]); dxlap = 0; }
    dfval  = Fval*dxval;
    dfgrad = Fgrad*dxval + Fval*dxgrad;
    dflap  = Flap*dxval + 2*dot(Fgrad,dxgrad) + dxlap*Fval; 
  }


  void evaluateDerivative_B(B_vars q, PosType r, RealType& dfval, GradType& dfgrad, RealType& dflap)
  {
    RealType dxval = 0, dxlap = 0;
    GradType dxgrad = 0;
    if(q == X) { dxval = -2*r[X]; dxgrad = -2*GradType(1,0,0); dxlap = 0; }
    if(q == Y) { dxval = -2*r[Y]; dxgrad = -2*GradType(0,1,0); dxlap = 0; }
    if(q == Z) { dxval = -2*r[Z]; dxgrad = -2*GradType(0,0,1); dxlap = 0; }
    dfval  = Fval*dxval;
    dfgrad = Fgrad*dxval + Fval*dxgrad;
    dflap  = Flap*dxval + 2*dot(Fgrad,dxgrad) + dxlap*Fval; 
  }

  void evaluateDerivatives(PosType r, 
                           std::vector<RealType>& dfval, 
                           std::vector<GradType>& dfgrad, 
                           std::vector<RealType>& dflap)
  {
    evaluate(r, Fval, Fgrad, Flap);
    int p = 0;
    if(opt_A[XX]) { evaluateDerivative_A(XX, r, dfval[p], dfgrad[p], dflap[p]); ++p; }
    if(opt_A[XY]) { evaluateDerivative_A(XY, r, dfval[p], dfgrad[p], dflap[p]); ++p; }
    if(opt_A[XZ]) { evaluateDerivative_A(XZ, r, dfval[p], dfgrad[p], dflap[p]); ++p; } 
    if(opt_A[YY]) { evaluateDerivative_A(YY, r, dfval[p], dfgrad[p], dflap[p]); ++p; } 
    if(opt_A[YZ]) { evaluateDerivative_A(YZ, r, dfval[p], dfgrad[p], dflap[p]); ++p; } 
    if(opt_A[ZZ]) { evaluateDerivative_A(ZZ, r, dfval[p], dfgrad[p], dflap[p]); ++p; }
    if(opt_B[X])  { evaluateDerivative_B(X,  r, dfval[p], dfgrad[p], dflap[p]); ++p; }
    if(opt_B[Y])  { evaluateDerivative_B(Y,  r, dfval[p], dfgrad[p], dflap[p]); ++p; }
    if(opt_B[Z])  { evaluateDerivative_B(Z,  r, dfval[p], dfgrad[p], dflap[p]); ++p; }
    if(opt_C)          
    {
      dfval[p] = Fval;
      dfgrad[p] = Fgrad;
      dflap[p] = Flap;
      ++p;
    }
  }

  GaussianCountingFunctor* makeClone(std::string clone_id)
  {
    app_log() << "  GaussianCountingFunctor::makeClone" << std::endl;
    GaussianCountingFunctor* rptr = new GaussianCountingFunctor(clone_id);
    // copy class variables defined in put
    //std::copy(A.begin(),A.end(), rptr->A.begin());
    //std::copy(B.begin(),B.end(), rptr->B.begin());
    for(int i = 0; i < A.size(); ++i)
      rptr->A[i] = A[i];
    for(int i = 0; i < B.size(); ++i)
      rptr->B[i] = B[i];
    rptr->C = C;
    rptr->opt_A.resize(opt_A.size());
    rptr->opt_B.resize(opt_B.size());
    for(int i = 0; i < opt_A.size(); ++i)
      rptr->opt_A[i] = opt_A[i];
    for(int i = 0; i < opt_B.size(); ++i)
      rptr->opt_B[i] = opt_B[i];
    //std::copy(opt_A.begin(), opt_A.end(), rptr->opt_A.begin());
    //std::copy(opt_B.begin(), opt_B.end(), rptr->opt_B.begin());
    rptr->opt_C = opt_C;
    //rptr->initialize();
    rptr->myVars = myVars;
    return rptr;
  }

  // this <= this * rhs
  void multiply_eq(const CountingFunctorBase* rhs) 
  {
    const GaussianCountingFunctor* rhs_gcf = dynamic_cast<const GaussianCountingFunctor*>(rhs);
    // check that conversion succeeded
    if(rhs_gcf == NULL)
      APP_ABORT("GaussianCountingFunctor::multiply_eq: can't multiply by non-gaussian");
    app_log() << "  GaussianCountingFunctor::multiply_eq" << std::endl;
    // calculate product
    A = A + rhs_gcf->A;
    B = B + rhs_gcf->B;
    C = C + rhs_gcf->C;
    initialize();
  }

  // this <- this / rhs
  void divide_eq(const CountingFunctorBase* rhs) 
  {
    const GaussianCountingFunctor* rhs_gcf = dynamic_cast<const GaussianCountingFunctor*>(rhs);
    // check that conversion succeeded
    if(rhs_gcf == NULL)
      APP_ABORT("GaussianCountingFunctor::divide_eq: can't divide by non-gaussian");
    app_log() << "GaussianCountingFunctor::divide_eq" << std::endl;
    // calculate product
    A = A - rhs_gcf->A;
    B = B - rhs_gcf->B;
    C = C - rhs_gcf->C;
    initialize();
    app_log() << "optvar after divide: " << std::endl;
    myVars.print(app_log());
  }

};


// unimplemented so far: perhaps make these using a set of GaussianCountingFunctors?
// this would primarily be to provide the closure for products and quotients.
struct FermiDiracCountingFunctor : public CountingFunctorBase 
{
  // same as pairwise-normalized gaussians: same parameters as a single gaussian (as above)

  // covariance matrix
  // displacement
  // constant
  // opt variables: std::vector bool (one for each parameter individually?)
  // id-strings for each optimizable variable
  // See GaussianFunctor

  // constructor
  FermiDiracCountingFunctor(std::string f_id) : CountingFunctorBase(f_id)
  {
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  void reportStatus(std::ostream& os) {}

  void resetParameters(const opt_variables_type& active) 
  {
  }

  void evaluate(PosType r, RealType& fval, GradType& fgrad, RealType& flap)
  {
  }
  
  void evaluateLog(PosType r, RealType& lval, GradType& lgrad, RealType& llap)
  {
  }

  void evaluateDerivatives(PosType r, std::vector<RealType>& df, std::vector<GradType>& dfgrad, std::vector<RealType>& dflap)
  {
  }

  void evaluateLogDerivatives(PosType r,
                              std::vector<RealType>& dlval,
                              std::vector<GradType>& dlgrad,
                              std::vector<RealType>& dllap)
  {
  }

  void evaluateLogTempDerivatives(PosType r,
                                  std::vector<RealType>& dlval)
  {
  }

  FermiDiracCountingFunctor* makeClone(std::string clone_id) {}

  void initialize() {}

  void multiply_eq(const CountingFunctorBase* rhs)
  {
  }

  void divide_eq(const CountingFunctorBase* rhs)
  {
  }

};

}

#endif
