#include<cmath>

#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{

enum
{
  DIM = OHMMS_DIM
};

using IntType     = int;
using RealType    = QMCTraits::RealType;
using IntVec      = TinyVector<IntType, DIM>;
using RealVec     = TinyVector<RealType, DIM>;
using RealMat     = Tensor<RealType, DIM>;
using PosArray    = std::vector<RealVec>;
using ChargeArray = std::vector<RealType>;



class RspaceMadelungTerm
{
  public:

  RealMat a;
  RealType rconst;

  RspaceMadelungTerm(RealMat a_in,RealType rconst_in)
  {
    a      = a_in;
    rconst = rconst_in;
  }

  RealType operator()(IntVec i)
  {
    RealVec Rv  = dot(i,a);
    RealType R  = std::sqrt(dot(Rv,Rv));
    RealType rm = std::erfc(rconst*R)/R;
    return rm;
  }
};



class KspaceMadelungTerm
{
  public:

  RealMat b;
  RealType kconst;
  RealType kfactor;

  KspaceMadelungTerm(RealMat b_in,RealType kconst_in,RealType kfactor_in)
  {
    b       = b_in;
    kconst  = kconst_in;
    kfactor = kfactor_in;
  }

  RealType operator()(IntVec i)
  {
    RealVec  Kv = dot(i,b);
    RealType K2 = dot(Kv,Kv);
    RealType km = kfactor*std::exp(kconst*K2)/K2;
    return km;
  }
};


class RspaceEwaldTerm
{
  public:

  RealVec r;
  RealMat a;
  RealType rconst;

  RspaceEwaldTerm(RealVec r_in,RealMat a_in,RealType rconst_in)
  {
    r      = r_in;
    a      = a_in;
    rconst = rconst_in;
  }

  RealType operator()(IntVec i)
  {
    RealVec Rv = dot(i,a);
    for(IntType d: {0,1,2})
      Rv[d] -= r[d];
    RealType R  = std::sqrt(dot(Rv,Rv));
    RealType rm = std::erfc(rconst*R)/R;
    return rm;
  }
};


class KspaceEwaldTerm
{
  public:

  RealVec r;
  RealMat b;
  RealType kconst;
  RealType kfactor;

  KspaceEwaldTerm(RealVec r_in,RealMat b_in,RealType kconst_in,RealType kfactor_in)
  {
    r       = r_in;
    b       = b_in;
    kconst  = kconst_in;
    kfactor = kfactor_in;
  }

  RealType operator()(IntVec i)
  {
     RealVec  Kv = dot(i,b);
     RealType K2 = dot(Kv,Kv);
     RealType Kr = dot(Kv,r);
     RealType km = kfactor*std::exp(kconst*K2)*std::cos(Kr)/K2;
     return km;
  }
};


template<typename T>
RealType gridSum(T& function,bool zero=true,RealType tol=1e-11)
{
  RealType dv = 1e99;
  RealType v  = 0.0;
  IntType im  = 0;
  IntType jm  = 0;
  IntType km  = 0;
  IntVec  iv;
  if(zero)
  {
    iv = 0;
    v += function(iv);
  }
  while(std::abs(dv)>tol)
  {
    dv = 0.0;
    im += 1;
    for(IntType i: {-im,im})
      for(IntType j=-jm; j<jm+1; ++j)
        for(IntType k=-km; k<km+1; ++k)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          dv += function(iv);
        }
    jm += 1;
    for(IntType j: {-jm,jm})
      for(IntType k=-km; k<km+1; ++k)
        for(IntType i=-im; i<im+1; ++i)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          dv += function(iv);
        }
    km += 1;
    for(IntType k: {-km,km})
      for(IntType i=-im; i<im+1; ++i)
        for(IntType j=-jm; j<jm+1; ++j)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          dv += function(iv);
        }
    v += dv;
  }

  return v;
}


RealType madelungSum(RealMat a,RealType tol=1e-10)
{
  RealType volume = std::abs(det(a));
  RealMat b       = 2*M_PI*transpose(inverse(a));
  RealType rconv  = 8*std::pow(3.*volume/(4*M_PI),1./3);
  RealType kconv  = 2*M_PI/rconv;

  RealType rconst  = kconv;
  RealType kconst  = -1./(4*std::pow(kconv,2));
  RealType kfactor = 4*M_PI/volume;

  RspaceMadelungTerm rfunc(a,rconst);
  KspaceMadelungTerm kfunc(b,kconst,kfactor);

  RealType cval = -M_PI/(std::pow(kconv,2)*volume)-2*kconv/std::sqrt(M_PI);
  RealType rval = gridSum(rfunc,false,tol);
  RealType kval = gridSum(kfunc,false,tol);

  RealType ms = cval + rval + kval;

  return ms;
}


RealType ewaldSum(RealVec r,RealMat a,RealType tol=1e-10)
{
  RealType volume = std::abs(det(a));
  RealMat b       = 2*M_PI*transpose(inverse(a));
  RealType rconv  = 8*std::pow(3.*volume/(4*M_PI),1./3);
  RealType kconv  = 2*M_PI/rconv;

  RealType rconst  = 1./(std::sqrt(2.)*kconv);
  RealType kconst  = -std::pow(kconv,2)/2;
  RealType kfactor = 4*M_PI/volume;

  RspaceEwaldTerm rfunc(r,a,rconst);
  KspaceEwaldTerm kfunc(r,b,kconst,kfactor);

  RealType cval = -2*M_PI*std::pow(kconv,2)/volume;
  RealType rval = gridSum(rfunc,true,tol);
  RealType kval = gridSum(kfunc,false,tol);

  RealType es = cval + rval + kval;

  return es;
}


RealType ewaldEnergy(RealMat a,PosArray R,ChargeArray Q,RealType tol=1e-10)
{
  IntType N = R.size();

  RealType qqmax=0.0;
  for(size_t i=0; i<N; ++i)
    qqmax = std::max(std::abs(Q[i]*Q[i]),qqmax);

  RealType vm = madelungSum(a,tol*2./qqmax);

  RealType ve = 0.0;
  // Madelung terms
  for(size_t i=0; i<N; ++i)
    ve += Q[i]*Q[i]*vm/2;
  // Pair interaction terms
  for(size_t i=0; i<N; ++i)
    for(size_t j=0; j<i; ++j)
    {
      RealType qq = Q[i]*Q[j];
      ve += qq*ewaldSum(R[i]-R[j],a,tol/qq);
    }
  return ve;
}

}
