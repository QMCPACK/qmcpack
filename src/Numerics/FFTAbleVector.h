#ifndef QMCPLUSPLUS_FFTABLEVECTOR_H
#define QMCPLUSPLUS_FFTABLEVECTOR_H
#include "ohmms-config.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"

// In the future should have a constructor that takes a Vector as an
// argument and then the FFTAbleVector's data is just the data from
// that vector (copy the pointer so that FFTing the data causes the 
// data in the original Vector to be transformed.

namespace APPNAMESPACE {
template<unsigned dimensions, typename precision, template<unsigned, class> class FFTEngine>
class FFTAbleVector : public Vector<precision> {
private:
  int NumPts;
  // would like to use precision here, but precision is a complex type
  // and this needs to be a real type commensurate with precision
  double ForwardNorm; 
  double BackwardNorm;
  TinyVector<int, dimensions> SizeDims;
  FFTEngine<dimensions, precision> MyEngine;
  void initialize(const TinyVector<int, dimensions>& DimSizes) {
    for (int i = 0; i < dimensions; i++) {
      SizeDims[i] = DimSizes[i];
    }
    for (unsigned i = 0; i < dimensions; i++) NumPts *= DimSizes[i];
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
                int sizeDim5 = 0, int sizeDim6 = 0, int sizeDim7 = 0, int sizeDim8 = 0) : 
  NumPts(1), ForwardNorm(1.0), BackwardNorm(1.0) {
    TinyVector<int, 8> DimSizes;
    DimSizes[0] = sizeDim1; DimSizes[1] = sizeDim2; DimSizes[2] = sizeDim3;
    DimSizes[3] = sizeDim4; DimSizes[4] = sizeDim5; DimSizes[5] = sizeDim6;
    DimSizes[6] = sizeDim7; DimSizes[7] = sizeDim8;
    initialize(DimSizes);
  }
  
  FFTAbleVector(const TinyVector<int, dimensions>& DimSizes) : NumPts(1), ForwardNorm(1.0), BackwardNorm(1.0) {
    initialize(DimSizes);
  }
  
  FFTAbleVector(const FFTAbleVector& rhs) : NumPts(rhs.NumPts), ForwardNorm(rhs.ForwardNorm),
  BackwardNorm(rhs.NumPts), SizeDims(rhs.SizeDims), Vector<precision>(rhs) {
    MyEngine.initialize(SizeDims.begin(), this->data());
  }
  
  inline void transformForward() {
    MyEngine.transformForward(this->data());
    (*this) *= ForwardNorm;
  }
  
  inline void transformBackward() {
    MyEngine.transformBackward(this->data());
    (*this) *= BackwardNorm;
  }
  
  inline void setForwardNorm(double fn) { ForwardNorm = fn; }
  inline void setBackwardNorm(double bn) { BackwardNorm = bn; }
};
}
#endif
