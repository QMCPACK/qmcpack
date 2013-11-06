//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VECTOR_VIEWER_H
#define QMCPLUSPLUS_VECTOR_VIEWER_H

namespace qmcplusplus
{

/** a class to map a memory sequence to a vector
 * @tparam T datatype
 */
template<typename T>
struct VectorViewer
{
  T* data_ptr;
  int data_size;

  inline VectorViewer(T* a, int n):data_ptr(a),data_size(n) {}

  inline T* data()
  {
    return data_ptr;
  }

  inline int size() const
  {
    return data_size;
  }

  inline T& operator[](int i)
  {
    return data_ptr[i];
  }

  inline T operator[](int i) const
  {
    return data_ptr[i];
  }
};

}

#endif
