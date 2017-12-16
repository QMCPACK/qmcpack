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
    
    


#ifndef OHMMS_PTCLPAIRLISTBASE_H
#define OHMMS_PTCLPAIRLISTBASE_H

#include <vector>
#include <string>

#ifndef OHMMS_PETE_VECTOR_H
#include "OhmmsPETE/OhmmsVector.h"
#endif
#ifndef OHMMS_TINYVECTOR_H
#include "OhmmsPETE/TinyVector.h"
#endif
#ifndef	OHMMS_TENSOR_H
#include "OhmmsPETE/Tensor.h"
#endif
/*  Container for nearest-neighbor data.
 *
 *  <ul>
 *  \li M[i] = locator for the first nearest-neighbor ptcl of the i-th ptcl
 *  \li M[i+1] - M[i] = number of nn ptcls of the i-thc ptcl
 *  \li Attributes of integer, scalar, vector and tenor types for a particle
 *      pair are accessed via index \f$l\f$ where \f$ M[i]\le l < M[i+1]\f$.
 *   <ul>
 *  \li J[l] = index of the neighbor atom
 *  \li R[l] = distance
 *  \li dR[l] = directional consine
 *  </ul>
 *  \li A new named attribute can be added.
 *  \li Estimating the mixmum bound for the number of nearest neighbors is important.
 *  </ul>
 */
template<class T, unsigned D>
class PtclPairListBase
{

public:

  typedef T               Scalar_t;
  typedef TinyVector<T,D> Vector_t;
  typedef Tensor<T,D>     Tensor_t;

  typedef std::vector<int>      IndexList_t;
  typedef std::vector<T>        ScalarList_t;
  typedef std::vector<Vector_t> VectorList_t;
  typedef std::vector<Tensor_t> TensorList_t;

  //typedef Vector<int, std::vector<int> >           IndexList_t;
  //typedef Vector<T, std::vector<T> >               ScalarList_t;
  //typedef Vector<Vector_t, std::vector<Vector_t> > VectorList_t;
  //typedef Vector<Tensor_t, std::vector<Tensor_t> > TensorList_t;
  //typedef std::map<std::string,int>  PairFieldMap_t;

  IndexList_t  I;//!< ID
  IndexList_t  M;//!< Locator
  IndexList_t  J;//!< Particle ID of a neighbor
  ScalarList_t R;//!< Distance between a pair ptcls
  VectorList_t dR;//!< Directional cosine between a pair ptcls

  PtclPairListBase();//!< Constructor
  ~PtclPairListBase();//!< Destructor

  //!< Number of particles
  inline int getLocalNum() const
  {
    return LocalNum;
  }

  //!< Number of pairs
  inline int getTotNadj() const
  {
    return M[LocalNum];
  }

  //!< Maximum number of neighbors per atom for memory management
  inline int getMaxNadj() const
  {
    return MaxNN;
  }

  //!< Returns a number of neighbors of the i-th ptcl.
  inline int nadj(int i) const
  {
    return M[i+1]-M[i];
  }

  //!< Returns the id of j-th neighbor for i-th ptcl
  inline int iadj(int i, int j) const
  {
    return J[M[i] +j];
  }

  //!< Location to insert a i-j pair
  inline int loc(int i, int j) const
  {
    return M[i] + j;
  }

  //!< Sets the maximum number of neighbors
  inline void setMaxNadj(int nnmax)
  {
    //if(nnmax > MaxNN)  MaxNN = nnmax;
    MaxNN = nnmax;
  }

  void resize(unsigned m); // Adds m pairs for each attribute
  void resize(unsigned m, unsigned maxnn)
  {
    //if(M.size() < m+1) {
    //  M.resize(m+1);
    //  I.resize(m);
    //}
    M.resize(m+1);
    I.resize(m);
    LocalNum = m;
    resize(maxnn*(m+1));
    MaxNN = maxnn;
    //if(maxnn > MaxNN) { MaxNN = maxnn; resize(MaxNN*M.size()); }
  }//!< Resizes neighbor lists. Removes everything.

  void print(std::ostream& os) const
  {
    os << "Total Number of Pairs" << M[LocalNum] << std::endl;
    os << "Total Number of Pair Attributes "  << NumAttrib << std::endl;
  }

  inline void reset()
  {
    //!< Reset the internal counter to accumulate nnlist
    CurI = -1;
    TotalNumAdj = 0; // total number of nn pairs
    CurNN = 0; // number of nn pairs for an atom i
    M[0] = 0;
  }

  inline void add(int nn, Scalar_t rsq, const Vector_t& dr)
  {
    R[nn] = rsq;
    dR[nn] = dr;
  }

  inline void add(int i, int j, Scalar_t rsq, const Vector_t& dr)
  {
    if(CurI != i)
    {
      CurI = i;
      if(CurI > 0)
        M[CurI] = M[CurI-1]+CurNN;
      CurNN = 0;
    }
    if(j>=0)
    {
      J[TotalNumAdj] = j;
      R[TotalNumAdj] = rsq;
      dR[TotalNumAdj] = dr;
      TotalNumAdj++;
      CurNN++;
    }
  }


protected:
  bool FullList;
  int NumAttrib, LocalNum, TotalNumAdj, MaxNN;
  int CurI, CurNN;
};

template<class T, unsigned D>
PtclPairListBase<T,D>::PtclPairListBase():LocalNum(0), MaxNN(0)
{
  TotalNumAdj = 0;  // total number of pairs
  NumAttrib = 3;
}

template<class T, unsigned D>
PtclPairListBase<T,D>::~PtclPairListBase()
{
}

template<class T, unsigned D>
void PtclPairListBase<T,D>::resize(unsigned m)
{
  if(m > J.size())
    // delete everything
  {
    J.resize(m); // index
    R.resize(m); // distance
    dR.resize(m);// directional cosine
  }//!< Creates m pairs for each attribute
}

#endif

