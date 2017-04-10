//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_VECTOR_ESTIMATO_IMPL_H
#define QMCPLUSPLUS_VECTOR_ESTIMATO_IMPL_H
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{

/** Class to manage observables of vector types
 *
 * Container and utility class for observables like gofr, sk.
 * Do not inherit from VectorEstimatorImpl but use it.
 */
template<typename T>
struct VectorEstimatorImpl
{

  ///name of the observables
  std::string Name;
  /** current accumulative data
   *
   * d_data.size() = 2 * size of vector elements
   * For the i-th element
   * - d_data[2*i] += \f$= \sum_{s} O_i^s\f$
   * - d_data[2*i+1] += \f$= \sum_{s} O_i^s O_i^s\f$
   */
  Vector<T> d_data;
  ///d_sum[i] block sum
  Vector<T> d_sum;
  /////d_sum2[i] block sum of squared
  //Vector<T> d_sum2;

  ///default constructor
  inline VectorEstimatorImpl() {}

  /** constructor
   * @param n size of the data
   */
  explicit inline VectorEstimatorImpl(int n)
  {
    resize(n);
  }

  /** constructor
   * @param a name of this data
   * @param n size of this data
   */
  VectorEstimatorImpl(const std::string& a, int n):Name(a)
  {
    resize(n);
  }

  ///copy constructor
  VectorEstimatorImpl(const VectorEstimatorImpl& est):
    d_data(est.d_data), d_sum(est.d_sum)
  {}

  ///destructo
  ~VectorEstimatorImpl() {}

  /** resize the data containers
   * @param n number of elements of the vector observable
   */
  inline void resize(int n)
  {
    d_data.resize(n);
    d_sum.resize(n);
    //d_sum2.resize(n);
  }

  inline size_t size() const
  {
    return d_sum.size();
  }
  inline void init()
  {
    d_data=T();
    d_sum=T();
    //d_sum2=T();
  }

  /// zero the active data
  inline void reset()
  {
    d_data=T();
  }

  /** accumulate expectation values
   * @param first start of vector data
   * @param last end of vector data
   * @param wit start of weight data
   */
  template<typename InputIterator, typename WeightIterator>
  inline void accumulate(InputIterator first, InputIterator last, WeightIterator wit, T norm)
  {
    typename std::vector<T>::iterator it(d_data.begin());
    while(first != last)
    {
      *it++ += norm*(*first++)*(*wit++);
      //T v=(*first1)*(*first2++);//w[i]*v[i]
      //(*it++)+=v;
      //(*it++)+=v*(*first1++);//w[i]*v[i]*v[i]
    }
  }

  /** accumulate expectation values
   * @param first start of vector data
   * @param last end of vector data
   * @param wm normalization factor
   */
  template<typename InputIterator, typename T1>
  inline void accumulate(InputIterator first, InputIterator last,  T1 wm)
  {
    typename std::vector<T>::iterator it(d_data.begin());
    while(first != last)
    {
      *it++  += wm*(*first++);
      //T v=wm*(*first1)*(*first2++);
      //(*it++)+=v;
      //(*it++)+=v*(*first1++);
    }
  }

  /** accumulate expectation values
   *\param awalker a single walker
   *\param wgt the weight
   */
  template<typename InputIterator, typename T1>
  inline void accumulate(InputIterator first, T1 wgt)
  {
    typename std::vector<T>::iterator it(d_data.begin());
    typename std::vector<T>::iterator it_end(d_data.end());
    while(it != it_end)
    {
      (*it++) += wgt*(*first++);
      //(*it++)+=wgt*(*first);
      //(*it++)+=wgt*(*first)*(*first);
      //++first;
    }
  }

  /** accumulate expectation values
   *\param awalker a single walker
   *\param wgt the weight
   */
  template<typename InputIterator>
  inline void accumulate(InputIterator first)
  {
    typename std::vector<T>::iterator it(d_data.begin());
    typename std::vector<T>::iterator it_end(d_data.end());
    while(it != it_end)
    {
      *it++ += *first++;
      //(*it++)+=(*first)*(*first);
      //++first;
    }
  }


  template<typename T1>
  inline void takeBlockAverage(T1 wgtnorm)
  {
    typename std::vector<T>::const_iterator it(d_data.begin());
    typename std::vector<T>::const_iterator it_end(d_data.end());
    typename std::vector<T>::iterator sit(d_sum.begin());
    //typename std::vector<T>::iterator s2it(d_sum2.begin());
    while(it != it_end)
    {
      (*sit++)=wgtnorm*(*it++);
      //(*s2it++)=wgtnorm*(*it++);
    }
  }

  template<typename ForwardIterator>
  ForwardIterator putMessage(ForwardIterator cur) const
  {
    copy(d_sum.begin(), d_sum.end(),cur);
    return cur+d_sum.size();
  }

  template<typename ForwardIterator>
  ForwardIterator getMessage(ForwardIterator cur)
  {
    ForwardIterator last=cur+d_sum.size();
    copy(cur,last,d_sum.begin());
    return last;
  }

  void print(std::ostream& os)
  {
    for(int i=0; i<d_sum.size(); i++)
      os << std::setw(20) << d_sum[i];
    os << std::endl;
  }
};

/** specialization for VectorTempImpl<double> */
template<>
struct HDFAttribIO<VectorEstimatorImpl<double> >
{
  typedef VectorEstimatorImpl<double> DataType_t;
  hid_t groupID;
  hid_t vsetID, vspaceID;
  //hid_t v2setID, v2spaceID;
  hsize_t maxdims[2];
  hsize_t curdims[2];
  hsize_t dims[2];
  hsize_t offset[2];
  DataType_t&  ref;

  HDFAttribIO(DataType_t& a): vsetID(-1),ref(a)
  {
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = a.size();
    curdims[0] = 1;
    curdims[1] = a.size();
    dims[0] = 1;
    dims[1] = a.size();
    offset[0]=0;
    offset[1]=0;
  }

  ~HDFAttribIO()
  {
    if(vsetID>-1)
    {
      //H5Sclose(v2spaceID);
      //H5Dclose(v2setID);
      H5Sclose(vspaceID);
      H5Dclose(vsetID);
    }
  }

  inline void reserve(hid_t grp)
  {
    const hsize_t RANK=2;
    vspaceID=H5Screate_simple(RANK, dims, maxdims);
    hid_t p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,RANK,dims);
    //vsetID= H5Dcreate(grp,name,H5T_NATIVE_DOUBLE,vspaceID,p);
    vsetID= H5Dcreate(grp,"value",H5T_NATIVE_DOUBLE,vspaceID,p);
    hid_t memspace = H5Screate_simple(RANK, dims, NULL);
    hid_t ret = H5Dwrite(vsetID, H5T_NATIVE_DOUBLE, memspace, vspaceID, H5P_DEFAULT, ref.d_sum.data());
    H5Sclose(memspace);
    H5Pclose(p);
    //v2spaceID=H5Screate_simple(RANK, dims, maxdims);
    //p = H5Pcreate (H5P_DATASET_CREATE);
    //H5Pset_chunk(p,RANK,dims);
    //v2setID= H5Dcreate(grp,"value2",H5T_NATIVE_DOUBLE,v2spaceID,p);
    //memspace = H5Screate_simple(RANK, dims, NULL);
    //ret = H5Dwrite(v2setID, H5T_NATIVE_DOUBLE, memspace, v2spaceID, H5P_DEFAULT, ref.d_sum2.data());
    //H5Sclose(memspace);
    //H5Pclose(p);
  }

  inline void write(hid_t grp, const char* name)
  {
    const hsize_t RANK=2;
    H5Dextend(vsetID,curdims);
    H5Sset_extent_simple(vspaceID,RANK,curdims,maxdims);
    H5Sselect_hyperslab(vspaceID, H5S_SELECT_SET, offset, NULL, dims, NULL);
    hid_t memspace = H5Screate_simple(RANK, dims, NULL);
    hid_t ret = H5Dwrite(vsetID, H5T_NATIVE_DOUBLE, memspace, vspaceID, H5P_DEFAULT, ref.d_sum.data());
    H5Sclose(memspace);
    //H5Dextend(v2setID,curdims);
    //H5Sset_extent_simple(v2spaceID,RANK,curdims,maxdims);
    //H5Sselect_hyperslab(v2spaceID, H5S_SELECT_SET, offset, NULL, dims, NULL);
    //memspace = H5Screate_simple(RANK, dims, NULL);
    //ret = H5Dwrite(v2setID, H5T_NATIVE_DOUBLE, memspace, v2spaceID, H5P_DEFAULT, ref.d_sum2.data());
    //H5Sclose(memspace);
    curdims[0]++;
    offset[0]++;
  }

  inline void read(hid_t grp, const char* name)
  {
  }

};
}

#endif
