//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_FFTABLEVECTOR_H
#define QMCPLUSPLUS_FFTABLEVECTOR_H
#include "config.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"

// In the future should have a constructor that takes a Vector as an
// argument and then the FFTAbleVector's data is just the data from
// that vector (copy the pointer so that FFTing the data causes the
// data in the original Vector to be transformed.

// also in the future move everything to the base class and just have
// an FFTEngine handed to this by some builder

namespace qmcplusplus
{

// dummy base
template<typename precision>
class FFTAbleVectorBase : public Vector<precision>
{
protected:
  int NumPts;
  // would like to use precision here, but precision is a complex type
  // and this needs to be a real type commensurate with precision
  double ForwardNorm;
  double BackwardNorm;
public:
  FFTAbleVectorBase() : NumPts(1), ForwardNorm(1.0), BackwardNorm(1.0)
  {
    ;
  }
  virtual ~FFTAbleVectorBase()
  {
    ;
  }
  FFTAbleVectorBase(const FFTAbleVectorBase& rhs) : Vector<precision>(rhs), NumPts(rhs.NumPts),
    ForwardNorm(rhs.ForwardNorm), BackwardNorm(rhs.BackwardNorm)
  {
    ;
  }
  virtual FFTAbleVectorBase* clone() const = 0;
  inline void setForwardNorm(double fn)
  {
    ForwardNorm = fn;
  }
  inline void setBackwardNorm(double bn)
  {
    BackwardNorm = bn;
  }
  virtual void transformForward() = 0;
  virtual void transformBackward() = 0;

};

template<unsigned dimensions, typename precision, template<unsigned, class> class FFTEngine>
class FFTAbleVector : public FFTAbleVectorBase<precision>
{
private:
  typedef FFTAbleVectorBase<precision> base;
  using base::NumPts;
  using base::ForwardNorm;
  using base::BackwardNorm;

  TinyVector<int, dimensions> SizeDims;
  FFTEngine<dimensions, precision> MyEngine;
  void initialize(const TinyVector<int, dimensions>& DimSizes)
  {
    for (int i = 0; i < dimensions; i++)
    {
      SizeDims[i] = DimSizes[i];
    }
    for (unsigned i = 0; i < dimensions; i++)
      NumPts *= DimSizes[i];
    this->resize(NumPts);
    MyEngine.initialize(SizeDims.begin(), this->data());
  }
  /*
  typedef T      Type_t
  typedef C      Contianer_t;
  typedef Vector<T,C> This_t;
  typedef typename Container_t::iterator iterator;
  typedef typename Container_t::const_iterator const_iterator;
   */
public:
  FFTAbleVector(int sizeDim1, int sizeDim2 = 0, int sizeDim3 = 0, int sizeDim4 = 0,
                int sizeDim5 = 0, int sizeDim6 = 0, int sizeDim7 = 0, int sizeDim8 = 0)
  {
    TinyVector<int, 8> DimSizes;
    DimSizes[0] = sizeDim1;
    DimSizes[1] = sizeDim2;
    DimSizes[2] = sizeDim3;
    DimSizes[3] = sizeDim4;
    DimSizes[4] = sizeDim5;
    DimSizes[5] = sizeDim6;
    DimSizes[6] = sizeDim7;
    DimSizes[7] = sizeDim8;
    initialize(DimSizes);
  }

  FFTAbleVector(const TinyVector<int, dimensions>& DimSizes)
  {
    initialize(DimSizes);
  }

  FFTAbleVector(const FFTAbleVector& rhs) : FFTAbleVectorBase<precision>(rhs)
  {
    SizeDims = rhs.SizeDims;
    MyEngine.initialize(SizeDims.begin(), this->data());
  }

  FFTAbleVector* clone() const
  {
    return new FFTAbleVector(*this);
  }

  inline void transformForward()
  {
    MyEngine.transformForward(this->data());
    (*this) *= ForwardNorm;
  }

  inline void transformBackward()
  {
    MyEngine.transformBackward(this->data());
    (*this) *= BackwardNorm;
  }
};
}
#endif
