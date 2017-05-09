//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_UNIFORMCARTESIANGRID_H
#define OHMMS_UNIFORMCARTESIANGRID_H

/**@file UniformCartesianGrid.h
 *@brief Class declarations for UniformCartesianGrid
 */
#include "Utilities/DistributedIndex.h"
#include "OhmmsPETE/TinyVector.h"
#include <iostream>

namespace qmcplusplus
{
/* \class UniformCartesianGrid
 * \brief Class to manage a uniform grid.
 *
 * Does nothing and needs to be specialized.
 */
template<class T, unsigned D>
struct UniformCartesianGrid {};

/** Specialization of UniformCartesianGrid<T,D> to 3-Dim domains.
 *
 *Regardless of the shape of the domain, the view of UniformCartesianGrid
 *is a cubic cell. The 3-D position vector is valid only if \f$v = ([0,1),[0,1),[0,1))\f$.
 *
 *The main function of this class is to partition a domain defined by
 *3 orthogonal axes, \f${\hat x}, {\hat y}\f$ and \f${\hat z}\f$.
 *A region or domain is specified by three indices (i,j,k)
 *and they satisfy the conditions \f[ i \in [0,NP[0])\f] \f[ j \in
 *[0,NP[1])\f] \f[ k \in [0,NP[2])\f] The subdomains are stored in an
 *one-dimensional array and key is used to map (i,j,k) to the
 *corresponding domain.  Member function key is provided in
 *anticipation of the need to use a sparse storage of the subdomains
 *to handle inhomogeneous systems.
 */
template<class T>
class UniformCartesianGrid<T,3>
{

public:

  typedef DistributedIndex::iterator iterator;
  typedef DistributedIndex::const_iterator const_iterator;

  ///default constructor
  inline UniformCartesianGrid()
  {
    NumGrids = 1;
    NP[0] = 1;
    NP[1] = 1;
    NP[2] = 1;
    Delta[0] = 1;
    Delta[1] = 1;
    Delta[2] =1;
    InvDelta[0] = 1;
    InvDelta[1] = 1;
    InvDelta[2] =1;
  }

  /////copy constructor
  //inline UniformCartesianGrid(const UniformCartesianGrid<T,3>& gr) {
  //  makeCopy(gr);
  //}

  /**constructor
   *@param ng a 3-Dim index vector that sets the partition
   */
  template<class IV>
  inline UniformCartesianGrid(const IV& ng)
  {
    setGrid(ng);
  }


  ///destructor
  virtual inline ~UniformCartesianGrid() { }

  ///copy operator
  UniformCartesianGrid<T,3>& operator=(const  UniformCartesianGrid<T,3>& gr)
  {
    makeCopy(gr);
    return *this;
  }

  ///copy function
  inline void makeCopy(const  UniformCartesianGrid<T,3>& gr)
  {
    setGrid(gr.NP);
  }

  ///return the number of grid in the i-th direction
  inline int size(int i) const
  {
    return NP[i];
  }

  /**set the parameters to partition a 3-Dim domain
   *@param ng a 3-Dim index vector that sets the partition
   */
  template<class IV>
  void setGrid(const IV& ng)
  {
    NumGrids = ng[0]*ng[1]*ng[2];
    NP[0] = ng[0];
    NP[1] = ng[1];
    NP[2] = ng[2];
    InvDelta[0] =static_cast<T>(ng[0]);
    InvDelta[1] =static_cast<T>(ng[1]);
    InvDelta[2] =static_cast<T>(ng[2]);
    Delta[0] =1./InvDelta[0];
    Delta[1] =1./InvDelta[1];
    Delta[2] =1./InvDelta[2];
//     CellKey.resize(ng[0]*ng[1]*ng[2]);
//     int ikey=0;
//     for(int ig=0; ig<ng[0]; ig++)
//       for(int jg=0; jg<ng[1]; jg++)
// 	for(int kg=0; kg<ng[2]; kg++,ikey++) CellKey(ikey)=ikey;
  }

  ///return the total number of sub domains/grids
  inline int getTotalNum() const
  {
    return NumGrids;
  }

  /**get a unique key for a domain(i,j,k)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the key for a domain(i,j,k)
     *
     *@warning broken. Should fail.
  */
  inline int key(int i, int j, int k) const
  {
    return CellKey[k+NP[2]*(j+NP[1]*i)];
  }

  /**get a index of a domain(i,j,k)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the storage index of a domain(i,j,k)
  */
  inline int loc(int i, int j, int k) const
  {
    if(i<0)
      i+=NP[0];
    if(i>=NP[0])
      i-= NP[0];
    if(j<0)
      j+=NP[1];
    if(j>=NP[1])
      j-= NP[1];
    if(k<0)
      i+=NP[2];
    if(k>=NP[2])
      k-= NP[2];
    return k+NP[2]*(j+NP[1]*i);
  }

  /**get a domain index of a 3-D vector v(x,y,z)
   @param x the position in the first dimension
   @param y the position in the second dimension
   @param z the position in the third dimension
   @return the storage index of a domain whose position is p(x,y,z)
   */
  inline int loc(T x, T y, T z) const
  {
    return int(z*InvDelta[2])+NP[2]*(int(y*InvDelta[1])+
                                     NP[1]*int(x*InvDelta[0]));
  }

  /**get a domain index of a 3-D vector p
   @param p a 3-D vector with operator[]
   @return the storage index of a domain whose position is p
   */
  template<class Pos_t>
  inline
  int loc(const Pos_t& p) const
  {
//     int i=int(p[0]*InvDelta[0]);
//     int j=int(p[1]*InvDelta[1]);
//     int k=int(p[2]*InvDelta[2]);
//     if(i<0) i+=NP[0];
//     else if(i>=NP[0]) i-= NP[0];
//     if(j<0) j+=NP[1];
//     else if(j>=NP[1]) j-= NP[1];
//     if(k<0) i+=NP[2];
//     else if(k>=NP[2]) k-= NP[2];
//     return k+NP[2]*(j+NP[1]*i);
    return
      int(p[2]*InvDelta[2])
      +NP[2]*(int(p[1]*InvDelta[1])+NP[1]*int(p[0]*InvDelta[0]));
  }

  /**get a 3-D index vector for a position r
   *@param r the position
   *@return a index vector containing the cell indices in the three directions
   */
  inline TinyVector<int,3> index(const TinyVector<T,3>& r)
  {
    return TinyVector<int,3>(static_cast<int>(r[0]*InvDelta[0]),
                             static_cast<int>(r[1]*InvDelta[1]),
                             static_cast<int>(r[2]*InvDelta[2]));
  }

  /**get a center position of a domain(i,j,k)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the position of the domain(i,j,k)
   */
  inline TinyVector<T,3> center(int i, int j, int k) const
  {
    return TinyVector<T,3>((static_cast<T>(i)+0.5)*Delta[0],
                           (static_cast<T>(j)+0.5)*Delta[1],
                           (static_cast<T>(k)+0.5)*Delta[2]);
  }

  /**get a center position of a domain(i,j,k)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the position of the domain(i,j,k)
   */
  inline TinyVector<T,3> vertex(int i, int j, int k) const
  {
    return TinyVector<T,3>((static_cast<T>(i))*Delta[0],
                           (static_cast<T>(j))*Delta[1],
                           (static_cast<T>(k))*Delta[2]);
  }

  /**distribute the domains over the processors
   *@param ntot the total number of processors
   *
   *The processor can be a MPI node or an OpenMP thread.
   *The domains are distributed over the processors as evenly as possible.
   *However, the computataionl nodes are not factored in for load balancing.
   */
  inline void distribute(int ntot)
  {
    I.distribute(ntot,NumGrids);
  }

  /**distribute the domains over the processors [first, last)
   *@param first iterator for the starting node
   *@param last iterator for the ending node
   *
   *A similar function to the void distribute(int ntot).
   */
  template <class _InputIterator>
  void distribute(_InputIterator first, _InputIterator last)
  {
    I.distribute(first,last);
  }

  /**print the domain partition information
   *@param os std::ostream to write to
   *
   *This is for debug/test only.
   */
  void printGrid(std::ostream& os) const
  {
    const int nw = 15;
    os << "number of sub regions " << NumGrids << std::endl;
    os << "Attribute distribution  = " << std::endl;
    os << "grid spacing = " << std::setw(nw) << Delta[0] << std::setw(nw) << Delta[1]
       << std::setw(nw) << Delta[2] << std::endl;
    TinyVector<T,3> origin;
    for(int ig=0; ig<NP[0]; ig++)
    {
      origin[0] = static_cast<T>(ig)*Delta[0];
      for(int jg=0; jg<NP[1]; jg++)
      {
        origin[1] = static_cast<T>(jg)*Delta[1];
        for(int kg=0; kg<NP[2]; kg++)
        {
          origin[2] = static_cast<T>(kg)*Delta[2];
          os << std::setw(nw) << origin[0] << " | " << std::setw(nw) << origin[0]+Delta[0]
             << std::setw(nw) << origin[1] << " | " << std::setw(nw) << origin[1]+Delta[1]
             << std::setw(nw) << origin[2] << " | " << std::setw(nw) << origin[2]+Delta[2]
             << std::endl;
        }
      }
    }
  }

  /** resize the grouping of domains
   *@param ng the number of groups
   */
  inline void resizeGroup(int ng)
  {
    PID.resize(ng);
  }

  ///return the processor ID of the group ig
  inline int group(int ig) const
  {
    return PID[ig];
  }
  ///assign the processor ID of the group ig
  inline int& group(int ig)
  {
    return PID[ig];
  }

  inline void createData()
  {
    I.create(NumGrids);
  }
  inline void clearData()
  {
    I.clear();
  }
  inline void addData(int ig, int iat)
  {
    I.add(ig,iat);
  }

  inline DistributedIndex& getDataSets()
  {
    return I;
  }
  inline const DistributedIndex& getDataSets() const
  {
    return I;
  }

  inline int getNumData() const
  {
    return I.getNumData();
  }
  inline int getMaxDataPerGrid() const
  {
    return I.capacity();
  }
  inline int getNumDataSets() const
  {
    return I.getNumDataSets();
  }
  inline int getNumData(int i) const
  {
    return I.size(i);
  }
  inline int firstData(int i) const
  {
    return I.M[i];
  }
  inline int lastData(int i) const
  {
    return I.M[i+1];
  }
  inline iterator beginData(int i)
  {
    return I.begin(i);
  }
  inline iterator endData(int i)
  {
    return I.end(i);
  }
  inline const_iterator beginData(int i) const
  {
    return I.begin(i);
  }
  inline const_iterator endData(int i) const
  {
    return I.end(i);
  }

  void printData(std::ostream& os) const
  {
    I.print(os);
  }

  int NumGrids, NP[3];
  T Delta[3],InvDelta[3];
  DistributedIndex I;
  std::vector<int> PID;
  std::vector<int> CellKey;
};

template<class T>
class UniformCartesianGrid<T,2>
{

public:

  typedef DistributedIndex::iterator iterator;
  typedef DistributedIndex::const_iterator const_iterator;

  ///default constructor
  inline UniformCartesianGrid()
  {
    NumGrids = 1;
    NP[0] = 1;
    NP[1] = 1;
    Delta[0] = 1;
    Delta[1] = 1;
    InvDelta[0] = 1;
    InvDelta[1] = 1;
  }

  ///copy constructor
  inline UniformCartesianGrid(const UniformCartesianGrid<T,2>& gr)
  {
    makeCopy(gr);
  }

  /**constructor
   *@param ng a 3-Dim index vector that sets the partition
   */
  template<class IV>
  inline UniformCartesianGrid(const IV& ng)
  {
    setGrid(ng);
  }


  ///destructor
  virtual inline ~UniformCartesianGrid() { }

  ///copy operator
  UniformCartesianGrid<T,2>& operator=(const  UniformCartesianGrid<T,2>& gr)
  {
    makeCopy(gr);
    return *this;
  }

  ///copy function
  inline void makeCopy(const  UniformCartesianGrid<T,2>& gr)
  {
    setGrid(gr.NP);
  }

  ///return the number of grid in the i-th direction
  inline int size(int i) const
  {
    return NP[i];
  }

  /**set the parameters to partition a 3-Dim domain
   *@param ng a 3-Dim index vector that sets the partition
   */
  template<class IV>
  void setGrid(const IV& ng)
  {
    NumGrids = ng[0]*ng[1];
    NP[0] = ng[0];
    NP[1] = ng[1];
    InvDelta[0] =static_cast<T>(ng[0]);
    InvDelta[1] =static_cast<T>(ng[1]);
    Delta[0] =1./InvDelta[0];
    Delta[1] =1./InvDelta[1];
  }

  ///return the total number of sub domains/grids
  inline int getTotalNum() const
  {
    return NumGrids;
  }

  /**get a unique key for a domain(i,j,k)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the key for a domain(i,j,k)
     *
     *@warning broken. Should fail.
  */
  inline int key(int i, int j) const
  {
    return CellKey[j+NP[1]*i];
  }

  /**get a index of a domain(i,j)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the storage index of a domain(i,j,k)
  */
  inline int loc(int i, int j) const
  {
    if(i<0)
      i+=NP[0];
    if(i>=NP[0])
      i-= NP[0];
    if(j<0)
      j+=NP[1];
    if(j>=NP[1])
      j-= NP[1];
    return j+NP[1]*i;
  }

  /**get a domain index of a 3-D vector v(x,y,z)
   @param x the position in the first dimension
   @param y the position in the second dimension
   @param z the position in the third dimension
   @return the storage index of a domain whose position is p(x,y,z)
   */
  inline int loc(T x, T y) const
  {
    return (int(y*InvDelta[1])+NP[1]*int(x*InvDelta[0]));
  }

  /**get a domain index of a 3-D vector p
   @param p a 3-D vector with operator[]
   @return the storage index of a domain whose position is p
   */
  template<class Pos_t>
  inline
  int loc(const Pos_t& p) const
  {
    return (int(p[1]*InvDelta[1])+NP[1]*int(p[0]*InvDelta[0]));
  }

  /**get a 3-D index vector for a position r
   *@param r the position
   *@return a index vector containing the cell indices in the three directions
   */
  inline TinyVector<int,2> index(const TinyVector<T,2>& r)
  {
    return TinyVector<int,2>(static_cast<int>(r[0]*InvDelta[0]),
                             static_cast<int>(r[1]*InvDelta[1]));
  }

  /**get a center position of a domain(i,j)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the position of the domain(i,j,k)
   */
  inline TinyVector<T,2> center(int i, int j) const
  {
    return TinyVector<T,2>((static_cast<T>(i)+0.5)*Delta[0],
                           (static_cast<T>(j)+0.5)*Delta[1]);
  }

  /**get a center position of a domain(i,j)
     @param i the index of the first dimension
     @param j the index of the second dimension
     @param k the index of the third dimension
     @return the position of the domain(i,j,k)
   */
  inline TinyVector<T,2> vertex(int i, int j) const
  {
    return TinyVector<T,2>((static_cast<T>(i))*Delta[0],
                           (static_cast<T>(j))*Delta[1]);
  }

  /**distribute the domains over the processors
   *@param ntot the total number of processors
   *
   *The processor can be a MPI node or an OpenMP thread.
   *The domains are distributed over the processors as evenly as possible.
   *However, the computataionl nodes are not factored in for load balancing.
   */
  inline void distribute(int ntot)
  {
    I.distribute(ntot,NumGrids);
  }

  /**distribute the domains over the processors [first, last)
   *@param first iterator for the starting node
   *@param last iterator for the ending node
   *
   *A similar function to the void distribute(int ntot).
   */
  template <class _InputIterator>
  void distribute(_InputIterator first, _InputIterator last)
  {
    I.distribute(first,last);
  }

  /**print the domain partition information
   *@param os std::ostream to write to
   *
   *This is for debug/test only.
   */
  void printGrid(std::ostream& os) const
  {
    const int nw = 15;
    os << "number of sub regions " << NumGrids << std::endl;
    os << "Attribute distribution  = " << std::endl;
    os << "grid spacing = " << std::setw(nw) << Delta[0] << std::setw(nw) << Delta[1]
       << std::endl;
    TinyVector<T,2> origin;
    for(int ig=0; ig<NP[0]; ig++)
    {
      origin[0] = static_cast<T>(ig)*Delta[0];
      for(int jg=0; jg<NP[1]; jg++)
      {
        origin[1] = static_cast<T>(jg)*Delta[1];
        os << std::setw(nw) << origin[0] << " | " << std::setw(nw) << origin[0]+Delta[0]
           << std::setw(nw) << origin[1] << " | " << std::setw(nw) << origin[1]+Delta[1]
           << std::endl;
      }
    }
  }

  /** resize the grouping of domains
   *@param ng the number of groups
   */
  inline void resizeGroup(int ng)
  {
    PID.resize(ng);
  }

  ///return the processor ID of the group ig
  inline int group(int ig) const
  {
    return PID[ig];
  }
  ///assign the processor ID of the group ig
  inline int& group(int ig)
  {
    return PID[ig];
  }

  inline void createData()
  {
    I.create(NumGrids);
  }
  inline void clearData()
  {
    I.clear();
  }
  inline void addData(int ig, int iat)
  {
    I.add(ig,iat);
  }

  inline DistributedIndex& getDataSets()
  {
    return I;
  }
  inline const DistributedIndex& getDataSets() const
  {
    return I;
  }

  inline int getNumData() const
  {
    return I.getNumData();
  }
  inline int getMaxDataPerGrid() const
  {
    return I.capacity();
  }
  inline int getNumDataSets() const
  {
    return I.getNumDataSets();
  }
  inline int getNumData(int i) const
  {
    return I.size(i);
  }
  inline int firstData(int i) const
  {
    return I.M[i];
  }
  inline int lastData(int i) const
  {
    return I.M[i+1];
  }
  inline iterator beginData(int i)
  {
    return I.begin(i);
  }
  inline iterator endData(int i)
  {
    return I.end(i);
  }
  inline const_iterator beginData(int i) const
  {
    return I.begin(i);
  }
  inline const_iterator endData(int i) const
  {
    return I.end(i);
  }

  void printData(std::ostream& os) const
  {
    I.print(os);
  }

  int NumGrids, NP[2];
  T Delta[2],InvDelta[2];
  DistributedIndex I;
  std::vector<int> PID;
  std::vector<int> CellKey;
};

}
#endif

