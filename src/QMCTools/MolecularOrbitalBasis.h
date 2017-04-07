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
    
    



#ifndef QMCPLUSPLUS_MOLECULARORBITALBASIS_H
#define QMCPLUSPLUS_MOLECULARORBITALBASIS_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

class DistanceTableData;

/** Class for a molecular orbital basis
 *
 *The molecular orbital \f$ \psi_i \f$ can be written as a linear
 *combination of basis functions \f$ \{\phi\} \f$ such that
 \f[
 \psi_i ({\bf r}_j) = \sum_I \sum_k C_{ikI} \phi_{ikI}({\bf r}_j-{\bf R}_I).
 \f]
 *This class performs the evaluation of the basis functions and their
 *derivatives for each of the N-particles in a configuration.  All that
 *is required to generate the actual molecular orbitals is to multiply
 *by the coefficient matrix.
 *
 *The template (C)entered(O)rbital(T)ype should provide the fuctions
 <ul>
 <li> evaluate(int source, int first, int nptcl, int offset,
 VM& y, GM& dy, VM& d2y) {
 </ul>
 *An example being SphericalOrbitalSet
 */
template<class COT>
class MolecularOrbitalBasis: public QMCTraits
{
public:

  const ParticleSet* IonConfig;
  ///constructor
  MolecularOrbitalBasis(): NumWalkers(1),myTable(0), IonConfig(0) { }

  ///constructor with ncenters
  explicit MolecularOrbitalBasis(int ncenters):NumWalkers(1), myTable(0), IonConfig(0)
  {
    for(int ic=0; ic<ncenters; ic++)
      AOs.push_back(0);
  }

  /**
   @param atable the distance table (ion-electron)
   @brief Assign the distance table (ion-electron) and
   *determine the total number of basis states.
  */
  void setTable(DistanceTableData* atable)
  {
    IonConfig=&(atable->origin());
    myTable = atable;
    //first set up the center list
    const ParticleSet& ptcl_ref = myTable->origin();
    I.resize(myTable->centers());
    for(int ic=0; ic<I.size(); ic++)
    {
      I[ic] = ptcl_ref.GroupID[ic];
    }
    //reset the distance table for the atomic orbitals
    for(int i=0; i<AOs.size(); i++)
      AOs[i]->setTable(myTable);
    //evaluate the total basis dimension and offset for each center
    Basis.resize(I.size()+1);
    Basis[0] = 0;
    for(int c=0; c<I.size(); c++)
    {
      Basis[c+1] = Basis[c]+AOs[I[c]]->getBasisSetSize();
    }
    TotalBasis = Basis[I.size()];
  }

  void reset()
  {
    // DO NOTHING
  }

  void resetParameters(VarRegistry<RealType>& optVariables)
  {
    for(int i=0; i<AOs.size(); i++)
      AOs[i]->resetParameters(optVariables);
  }


  /** return the size of the basis set
   */
  inline int size() const
  {
    return TotalBasis;
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P)
  {
    LOGMSG("MolecularOrbitalBasis::resetTargetParticleSet")
    myTable = DistanceTable::add(*IonConfig,P,DT_AOS);
    for(int i=0; i<AOs.size(); i++)
      AOs[i]->setTable(myTable);
  }

  /**
     @param nptcl number of particles
     @brief resize the containers for data
     for nptcl particles and TotalBasis basis functions
  */
  inline void resize(int nptcl)
  {
    NumPtcls = nptcl;
    Y.resize(nptcl,TotalBasis);
    dY.resize(nptcl,TotalBasis);
    d2Y.resize(nptcl,TotalBasis);
  }


  /** Evalaute the values of  basis functions
     @param P input configuration containing N particles
     @brief For each center, evaluate all the Atomic Orbitals belonging to that center.
    *
    The data members contain
    - \f[ Y[i,j] =  \phi_j(r_i-R_I) \f]
    - \f[ dY[i,j] = {\bf \nabla}_i \phi_j(r_i-R_I) \f]
    - \f[ d2Y[i,j] = \nabla^2_i \phi_j(r_i-R_I), \f]
    where \f$ {\bf R_I} \f$ is the correct center for the basis function \f$ \phi_j \f$
  */
  inline void
  evaluate(const ParticleSet& P)
  {
    for(int c=0; c<I.size(); c++)
    {
      AOs[I[c]]->evaluate(c,0,P.getTotalNum(),Basis[c],Y,dY,d2Y);
    }
  }

  inline void
  evaluate(const ParticleSet& P, int iat)
  {
    for(int c=0; c<I.size(); c++)
    {
      AOs[I[c]]->evaluate(c,iat,Basis[c],Y);
    }
  }

  inline void
  evaluateAll(const ParticleSet& P, int iat)
  {
    for(int c=0; c<I.size(); c++)
    {
      AOs[I[c]]->evaluate(c,iat,Basis[c],Y,dY,d2Y);
    }
  }

  /** add a new set of Centered Atomic Orbitals
   *@param aos a set of Centered Atomic Orbitals
   */
  inline void add(COT* aos)
  {
    ERRORMSG("!!!!!! this should not be used !!!!!! ")
    AOs.push_back(aos);
  }

  /** add a new set of Centered Atomic Orbitals
   * @param aos a set of Centered Atomic Orbitals
   * @param icenter the index of the center
   */
  inline void add(COT* aos, int icenter)
  {
    AOs[icenter]=aos;
  }

  void print(std::ostream& os)
  {
    for(int i=0; i<Y.rows(); i++)
    {
      for(int j=0; j<Y.cols(); j++)
      {
        os << Y(i,j) << " " << dY(i,j) << " " << d2Y(i,j) << std::endl;
      }
    }
  }

  ///row i of matrix Y
  inline const ValueType* restrict y(int i)
  {
    return &Y(i,0);
  }
  ///row i of vector matrix dY
  inline const GradType* restrict dy(int i)
  {
    return &dY(i,0);
  }
  ///row i of matrix d2Y
  inline const ValueType* restrict d2y(int i)
  {
    return &d2Y(i,0);
  }

#ifdef USE_FASTWALKER
  inline const ValueType* restrict y(int iw, int ia)
  {
    return &Y(iw+NumWalkers*ia,0);
  }
  inline const GradType* restrict dy(int iw, int ia)
  {
    return &dY(iw+NumWalkers*ia,0);
  }
  inline const ValueType* restrict d2y(int iw, int ia)
  {
    return &d2Y(iw+NumWalkers*ia,0);
  }
#else
  inline const ValueType* restrict y(int iw, int ia)
  {
    return &Y(iw*NumPtcls+ia,0);
  }
  inline const GradType* restrict dy(int iw, int ia)
  {
    return &dY(iw*NumPtcls+ia,0);
  }
  inline const ValueType* restrict d2y(int iw, int ia)
  {
    return &d2Y(iw*NumPtcls+ia,0);
  }
#endif

  ///the number of particles
  int NumPtcls;
  ///the number of walkers
  int NumWalkers;
  ///total number of basis functions
  int TotalBasis;
  ///container for the id's of the centers (ions),
  ///several centers may share the same id
  std::vector<int>  I;
  ///container to store the offsets of the basis functions,
  ///the number of basis states for center J is Basis[J+1]-Basis[J]
  std::vector<int>  Basis;
  ///container for the pointers to the Atomic Orbitals,
  ///the size of this container being determined by the number
  ///of unique centers
  std::vector<COT*> AOs;
  ///matrix to store values \f$ Y[i,j] = \phi_j(r_i) \f$
  Matrix<ValueType> Y;
  ///matrix to store gradients \f$ dY[i,j] = {\bf \nabla}_i \phi_j(r_i) \f$
  Matrix<GradType>  dY;
  ///matrix to store laplacians \f$ d2Y[i,j] = \nabla^2_i \phi_j(r_i) \f$
  Matrix<ValueType> d2Y;
  ///the distance table (ion-electron)
  DistanceTableData* myTable;
};

}
#endif


