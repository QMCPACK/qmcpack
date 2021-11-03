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

#ifndef QMCPLUSPLUS_GAUSSIAN_FUNCTOR_H
#define QMCPLUSPLUS_GAUSSIAN_FUNCTOR_H

#include "OhmmsData/AttributeSet.h"
#include "VariableSet.h"
#include <array>

namespace qmcplusplus
{

// doesn't inherit from OptimizableFunctorBase since this is a function of the entire position vector
class CountingGaussian 
{
  using RealType = QMCTraits::RealType;
  using PosType = QMCTraits::PosType;
  using GradType = QMCTraits::PosType; // complex not implemented
  using TensorType = QMCTraits::TensorType;

  using real_type = optimize::VariableSet::real_type;
  using opt_variables_type = optimize::VariableSet;

  // enumerations for axis parameters
  enum A_vars
  {
    XX,
    XY,
    XZ,
    YY,
    YZ,
    ZZ,
    NUM_A
  };
  enum B_vars
  {
    X,
    Y,
    Z,
    DIM
  };


  // opt variables: vector of bools: one for each parameter
  std::vector<bool> opt_A;
  std::vector<bool> opt_B;
  bool opt_C;

  // id string
  std::string id;
  
  // reference gaussian
  CountingGaussian* gref = NULL;

  // most recent evaluations
  RealType Fval;
  GradType Fgrad;
  RealType Flap;

public:
  TensorType A;
  PosType B;
  RealType C;

  // optimizable variables
  opt_variables_type myVars;

  CountingGaussian(std::string fid) { id = fid; }

  // constructor for voronoi region
  CountingGaussian(std::string fid, RealType alpha, PosType r, bool opt)
  { 
    id = fid;
    // set opt flags
    opt_A.resize(6);
    opt_B.resize(3);
    for(auto it = opt_A.begin(); it != opt_A.end(); ++it)
      *it = opt;
    for(auto it = opt_B.begin(); it != opt_B.end(); ++it)
      *it = opt;
    opt_C = opt;
    // set A
    A(0,0) = A(1,1) = A(2,2) = -1.0*alpha;
    A(0,1) = A(1,0) = A(0,2) = A(2,0) = A(1,2) = A(2,1) = 0;
    // set B
    d_to_b(r, A, B);
    // set C
    k_to_c(0, A, r, C);
  }
  
  // initialize with another gaussian reference
  void initialize(CountingGaussian* ref)
  {
    // register and update optimizable variables to current values
    myVars.insert(id + "_A_xx", A(X, X), opt_A[XX], optimize::OTHER_P);
    myVars[id + "_A_xx"] = A(X, X);
    myVars.insert(id + "_A_xy", A(X, Y), opt_A[XY], optimize::OTHER_P);
    myVars[id + "_A_xy"] = A(X, Y);
    myVars.insert(id + "_A_xz", A(X, Z), opt_A[XZ], optimize::OTHER_P);
    myVars[id + "_A_xz"] = A(X, Z);
    myVars.insert(id + "_A_yy", A(Y, Y), opt_A[YY], optimize::OTHER_P);
    myVars[id + "_A_yy"] = A(Y, Y);
    myVars.insert(id + "_A_yz", A(Y, Z), opt_A[YZ], optimize::OTHER_P);
    myVars[id + "_A_yz"] = A(Y, Z);
    myVars.insert(id + "_A_zz", A(Z, Z), opt_A[ZZ], optimize::OTHER_P);
    myVars[id + "_A_zz"] = A(Z, Z);
    myVars.insert(id + "_B_x", B[X], opt_B[X], optimize::OTHER_P);
    myVars[id + "_B_x"] = B[X];
    myVars.insert(id + "_B_y", B[Y], opt_B[Y], optimize::OTHER_P);
    myVars[id + "_B_y"] = B[Y];
    myVars.insert(id + "_B_z", B[Z], opt_B[Z], optimize::OTHER_P);
    myVars[id + "_B_z"] = B[Z];
    myVars.insert(id + "_C", C, opt_C, optimize::OTHER_P);
    myVars[id + "_C"] = C;
    // transform according to reference
    if(gref != NULL)
      this->divide_eq(gref);
  }

  void restore(int iat) {}

  void checkInVariables(opt_variables_type& active) { active.insertFrom(myVars); }


  void checkOutVariables(const opt_variables_type& active) { myVars.getIndex(active); }


  void resetParameters(const opt_variables_type& active)
  {
    int ia;
    std::string vid;
    if (opt_A[XX])
    {
      vid = id + "_A_xx";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      A(0, 0) = std::real(myVars[vid]);
    }
    if (opt_A[YY])
    {
      vid = id + "_A_yy";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      A(1, 1) = std::real(myVars[vid]);
    }
    if (opt_A[ZZ])
    {
      vid = id + "_A_zz";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      A(2, 2) = std::real(myVars[vid]);
    }
    if (opt_A[XY])
    {
      vid = id + "_A_xy";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      A(1, 0) = A(0, 1) = std::real(myVars[vid]);
    }
    if (opt_A[XZ])
    {
      vid = id + "_A_xz";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      A(2, 0) = A(0, 2) = std::real(myVars[vid]);
    }
    if (opt_A[YZ])
    {
      vid = id + "_A_yz";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      A(2, 1) = A(1, 2) = std::real(myVars[vid]);
    }
    if (opt_B[X])
    {
      vid = id + "_B_x";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      B[X] = std::real(myVars[vid]);
    }
    if (opt_B[Y])
    {
      vid = id + "_B_y";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      B[Y] = std::real(myVars[vid]);
    }
    if (opt_B[Z])
    {
      vid = id + "_B_z";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      B[Z] = std::real(myVars[vid]);
    }
    if (opt_C)
    {
      vid = id + "_C";
      ia = active.getLoc(vid);
      myVars[vid] = active[ia];
      C = std::real(myVars[vid]);
    }
    // transform according to reference
    if(gref != NULL)
      this->divide_eq(gref);
  }

  void reportStatus(std::ostream& os)
  {
    os << "    type: CountingGaussian" << std::endl;
    os << "    id: " << id << std::endl;
    app_log() << std::endl;
    myVars.print(os, 6, true);
    app_log() << std::endl;
  }

  std::unique_ptr<CountingGaussian> makeClone(std::string fid) const
  {
    auto rptr = std::make_unique<CountingGaussian>(fid);
    for (int i = 0; i < A.size(); ++i)
      rptr->A[i] = A[i];
    for (int i = 0; i < B.size(); ++i)
      rptr->B[i] = B[i];
    rptr->C = C;
    rptr->opt_A.resize(opt_A.size());
    rptr->opt_B.resize(opt_B.size());
    for (int i = 0; i < opt_A.size(); ++i)
      rptr->opt_A[i] = opt_A[i];
    for (int i = 0; i < opt_B.size(); ++i)
      rptr->opt_B[i] = opt_B[i];
    rptr->opt_C  = opt_C;
    rptr->myVars = myVars;
    return rptr;
  }

  // builds the full matrix from its upper triangle
  void triang_to_matrix(const std::array<RealType, 6>& triang, TensorType& matrix)
  {
    auto it = triang.begin();
    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j)
      {
        matrix(i, j) = matrix(j, i) = *it;
        ++it;
      }
  }

  void d_to_b(const PosType& d, const TensorType& a, PosType& b) { b = dot(a, d); }

  void k_to_c(const RealType& k, const TensorType& a, const PosType& d, RealType& c) { c = k + dot(d, dot(a, d)); }

  bool put(xmlNodePtr cur)
  {
    //app_log() << "CountingGaussian::put" << std::endl;
    bool put_A = false, put_B = false, put_C = false, put_D = false, put_K = false;
    // alternate inputs
    std::array<RealType, 6> A_triang;
    PosType D  = 0;
    RealType K = 0;
    cur        = cur->xmlChildrenNode;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if (cname == "var")
      {
        std::string opt  = "false";
        std::string name = "none";
        OhmmsAttributeSet rAttrib;
        rAttrib.add(name, "name");
        rAttrib.add(opt, "opt");
        rAttrib.add(opt, "opt_bits");
        rAttrib.put(cur);
        // check if opt is a binary string
        bool is_bitstr = std::all_of(opt.begin(), opt.end(), [&](char c) { return c == '0' || c == '1'; });
        std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))std::tolower);
        std::transform(name.begin(), name.end(), name.begin(), (int (*)(int))std::tolower);
        if (name == "a") // input is the upper-triangle of A
        {
          opt_A.resize(6);
          if (opt.size() == 6 && is_bitstr)
            std::transform(opt.begin(), opt.end(), opt_A.begin(), [&](char c) { return (c == '1'); });
          // default opt = true
          else
            std::fill(opt_A.begin(), opt_A.end(), (opt == "true"));
          put_A = putContent(A_triang.begin(), A_triang.end(), cur);
          // perform conversion
          triang_to_matrix(A_triang, A);
        }
        if (name == "b")
        {
          opt_B.resize(3);
          if (opt.size() == 3 && is_bitstr)
            std::transform(opt.begin(), opt.end(), opt_B.begin(), [&](char c) { return (c == '1'); });
          // default opt = true
          else
            std::fill(opt_B.begin(), opt_B.end(), (opt == "true"));
          put_B = putContent(B.begin(), B.end(), cur);
        }
        if (name == "d")
        {
          opt_B.resize(3);
          if (opt.size() == 3 && is_bitstr)
            std::transform(opt.begin(), opt.end(), opt_B.begin(), [&](char c) { return (c == '1'); });
          else
            std::fill(opt_B.begin(), opt_B.end(), (opt == "true"));
          put_D = putContent(D.begin(), D.end(), cur);
        }
        if (name == "c")
        {
          opt_C = (opt == "true");
          put_C = putContent(C, cur);
        }
        if (name == "k")
        {
          opt_C = (opt == "true");
          put_K = putContent(K, cur);
        }
      }
      cur = cur->next;
    }
    // convert B <=> D
    if (put_B && put_D)
      APP_ABORT("CountingGaussian::put: overdetermined: both B and D are specified");
    if (put_C && put_K)
      APP_ABORT("CountingGaussian::put: overdetermined: both C and K are specified");
    // convert D,B using A
    if (put_D)
      d_to_b(D, A, B);
    if (put_K)
      k_to_c(K, A, D, C);
    // check that we got everything
    if (!(put_A && (put_B || put_D) && (put_C || put_K)))
    {
      std::ostringstream err;
      err << "CountingGaussian::put:  required variable unspecified: ";
      err << "put_A: " << put_A << ", put_B: " << put_B << ", put_C: " << put_C << ", put_D: " << put_D
          << ", put_K: " << put_K << std::endl;
      //      APP_ABORT(err.str());
      app_log() << err.str();
    }
    return true;
  }

  void divide_eq(const CountingGaussian* rhs)
  {
    A = A - rhs->A;
    B = B - rhs->B;
    C = C - rhs->C;
  }


  // f = std::exp(x^T A x - 2B^T x + C)
  void evaluate(PosType r, RealType& fval, GradType& fgrad, RealType& flap)
  {
    PosType Ar = dot(A, r);
    RealType x = dot(Ar - 2 * B, r) + C;
    fval       = std::exp(x);
    fgrad      = 2 * (Ar - B) * fval;
    flap       = 4 * dot(Ar - B, Ar - B) * fval + 2 * trace(A) * fval;
  }

  void evaluateLog(PosType r, RealType& lval, GradType& lgrad, RealType& llap)
  {
    PosType Ar = dot(A, r);
    lval       = dot(Ar - 2 * B, r) + C;
    lgrad      = 2 * (Ar - B);
    llap       = 2 * trace(A);
  }

  void evaluateDerivative_A(A_vars q, PosType r, RealType& dfval, GradType& dfgrad, RealType& dflap)
  {
    // Fval, Fgrad, Flap are up-to-date function values
    // x correponds to the exponent value: x = ln f; dx are param derivs
    RealType dxval = 0, dxlap = 0;
    GradType dxgrad = 0;
    if (q == XX)
    {
      dxval  = r[X] * r[X];
      dxgrad = GradType(2 * r[X], 0, 0);
      dxlap  = 2;
    }
    if (q == YY)
    {
      dxval  = r[Y] * r[Y];
      dxgrad = GradType(0, 2 * r[Y], 0);
      dxlap  = 2;
    }
    if (q == ZZ)
    {
      dxval  = r[Z] * r[Z];
      dxgrad = GradType(0, 0, 2 * r[Z]);
      dxlap  = 2;
    }
    // off-diagonal terms: factor of two is since A is symmetric
    if (q == XY)
    {
      dxval  = 2 * r[X] * r[Y];
      dxgrad = 2 * GradType(r[Y], r[X], 0);
      dxlap  = 0;
    }
    if (q == XZ)
    {
      dxval  = 2 * r[X] * r[Z];
      dxgrad = 2 * GradType(r[Z], 0, r[X]);
      dxlap  = 0;
    }
    if (q == YZ)
    {
      dxval  = 2 * r[Y] * r[Z];
      dxgrad = 2 * GradType(0, r[Z], r[Y]);
      dxlap  = 0;
    }
    dfval  = Fval * dxval;
    dfgrad = Fgrad * dxval + Fval * dxgrad;
    dflap  = Flap * dxval + 2 * dot(Fgrad, dxgrad) + dxlap * Fval;
  }

  void evaluateDerivative_B(B_vars q, PosType r, RealType& dfval, GradType& dfgrad, RealType& dflap)
  {
    RealType dxval = 0, dxlap = 0;
    GradType dxgrad = 0;
    if (q == X)
    {
      dxval  = -2 * r[X];
      dxgrad = -2 * GradType(1, 0, 0);
      dxlap  = 0;
    }
    if (q == Y)
    {
      dxval  = -2 * r[Y];
      dxgrad = -2 * GradType(0, 1, 0);
      dxlap  = 0;
    }
    if (q == Z)
    {
      dxval  = -2 * r[Z];
      dxgrad = -2 * GradType(0, 0, 1);
      dxlap  = 0;
    }
    dfval  = Fval * dxval;
    dfgrad = Fgrad * dxval + Fval * dxgrad;
    dflap  = Flap * dxval + 2 * dot(Fgrad, dxgrad) + dxlap * Fval;
  }

  void evaluateDerivatives(PosType r,
                           std::vector<RealType>& dfval,
                           std::vector<GradType>& dfgrad,
                           std::vector<RealType>& dflap)
  {
    evaluate(r, Fval, Fgrad, Flap);
    int p = 0;
    if (opt_A[XX])
    {
      evaluateDerivative_A(XX, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_A[XY])
    {
      evaluateDerivative_A(XY, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_A[XZ])
    {
      evaluateDerivative_A(XZ, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_A[YY])
    {
      evaluateDerivative_A(YY, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_A[YZ])
    {
      evaluateDerivative_A(YZ, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_A[ZZ])
    {
      evaluateDerivative_A(ZZ, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_B[X])
    {
      evaluateDerivative_B(X, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_B[Y])
    {
      evaluateDerivative_B(Y, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_B[Z])
    {
      evaluateDerivative_B(Z, r, dfval[p], dfgrad[p], dflap[p]);
      ++p;
    }
    if (opt_C)
    {
      dfval[p]  = Fval;
      dfgrad[p] = Fgrad;
      dflap[p]  = Flap;
      ++p;
    }
  }

  // evaluate parameter derivatives of log of the wavefunction
  void evaluateLogDerivatives(PosType r,
                              std::vector<RealType>& dlval,
                              std::vector<GradType>& dlgrad,
                              std::vector<RealType>& dllap)
  {
    int p = 0;
    if (opt_A[XX])
    {
      dlval[p]  = r[X] * r[X];
      dlgrad[p] = 2 * GradType(r[X], 0, 0);
      dllap[p]  = 2;
      ++p;
    }
    if (opt_A[XY])
    {
      dlval[p]  = 2 * r[X] * r[Y];
      dlgrad[p] = 2 * GradType(r[Y], r[X], 0);
      dllap[p]  = 0;
      ++p;
    }
    if (opt_A[XZ])
    {
      dlval[p]  = 2 * r[X] * r[Z];
      dlgrad[p] = 2 * GradType(r[Z], 0, r[X]);
      dllap[p]  = 0;
      ++p;
    }
    if (opt_A[YY])
    {
      dlval[p]  = r[Y] * r[Y];
      dlgrad[p] = 2 * GradType(0, r[Y], 0);
      dllap[p]  = 2;
      ++p;
    }
    if (opt_A[YZ])
    {
      dlval[p]  = 2 * r[Y] * r[Z];
      dlgrad[p] = 2 * GradType(0, r[Z], r[Y]);
      dllap[p]  = 0;
      ++p;
    }
    if (opt_A[ZZ])
    {
      dlval[p]  = r[Z] * r[Z];
      dlgrad[p] = 2 * GradType(0, 0, r[Z]);
      dllap[p]  = 2;
      ++p;
    }
    if (opt_B[X])
    {
      dlval[p]  = -2 * r[X];
      dlgrad[p] = -2 * GradType(1, 0, 0);
      dllap[p]  = 0;
      ++p;
    }
    if (opt_B[Y])
    {
      dlval[p]  = -2 * r[Y];
      dlgrad[p] = -2 * GradType(0, 1, 0);
      dllap[p]  = 0;
      ++p;
    }
    if (opt_B[Z])
    {
      dlval[p]  = -2 * r[Z];
      dlgrad[p] = -2 * GradType(0, 0, 1);
      dllap[p]  = 0;
      ++p;
    }
    if (opt_C)
    {
      dlval[p]  = 1;
      dlgrad[p] = 0;
      dllap[p]  = 0;
      ++p;
    }
  }

  void evaluateLogTempDerivatives(PosType r, std::vector<RealType>& dlval)
  {
    int p = 0;
    if (opt_A[XX])
    {
      dlval[p] = r[X] * r[X];
      ++p;
    }
    if (opt_A[XY])
    {
      dlval[p] = 2 * r[X] * r[Y];
      ++p;
    }
    if (opt_A[XZ])
    {
      dlval[p] = 2 * r[X] * r[Z];
      ++p;
    }
    if (opt_A[YY])
    {
      dlval[p] = r[Y] * r[Y];
      ++p;
    }
    if (opt_A[YZ])
    {
      dlval[p] = 2 * r[Y] * r[Z];
      ++p;
    }
    if (opt_A[ZZ])
    {
      dlval[p] = r[Z] * r[Z];
      ++p;
    }
    if (opt_B[X])
    {
      dlval[p] = -2 * r[X];
      ++p;
    }
    if (opt_B[Y])
    {
      dlval[p] = -2 * r[Y];
      ++p;
    }
    if (opt_B[Z])
    {
      dlval[p] = -2 * r[Z];
      ++p;
    }
    if (opt_C)
    {
      dlval[p] = 1;
      ++p;
    }
  }

  void evaluate_print(std::ostream& os, const ParticleSet& P)
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
    for (auto it = P.R.begin(); it != P.R.end(); ++it)
    {
      r     = *it;
      Ar    = dot(A, r);
      x     = dot(Ar - 2 * B, r) + C;
      fval  = std::exp(x);
      fgrad = 2 * (Ar - B) * fval;
      flap  = 4 * dot(Ar - B, Ar - B) * fval + 2 * trace(A) * fval;

      r_vec.push_back(r);
      Ar_vec.push_back(Ar);
      x_vec.push_back(x);
      fval_vec.push_back(fval);
      fgrad_vec.push_back(fgrad);
      flap_vec.push_back(flap);
    }
    os << "CountingGaussian::evaluate_print: id: " << id << std::endl;
    os << "r: ";
    std::copy(r_vec.begin(), r_vec.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "Ar: ";
    std::copy(Ar_vec.begin(), Ar_vec.end(), std::ostream_iterator<PosType>(os, ", "));
    os << std::endl << "x: ";
    std::copy(x_vec.begin(), x_vec.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "fval: ";
    std::copy(fval_vec.begin(), fval_vec.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl << "fgrad: ";
    std::copy(fgrad_vec.begin(), fgrad_vec.end(), std::ostream_iterator<GradType>(os, ", "));
    os << std::endl << "flap: ";
    std::copy(flap_vec.begin(), flap_vec.end(), std::ostream_iterator<RealType>(os, ", "));
    os << std::endl;
  }
};

} // namespace qmcplusplus
#endif
