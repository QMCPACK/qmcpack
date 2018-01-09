//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SoaCuspCorrectionBasisSet.h
 *
 * Convert CorrectingOrbitalBasisSet using MultiQuinticSpline1D<T>
 */
#ifndef QMCPLUSPLUS_SOA_CUSPCORRECTION_BASISSET_H
#define QMCPLUSPLUS_SOA_CUSPCORRECTION_BASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

  struct CuspCorrectionBase
  {
    typedef QMCTraits::Value_type ValueType;
    typedef SoaVectorContainer<value_type,OHMMS_PRECISION+2> VGLVector_t;

    virtual CuspCorrectionBase* makeClone() const = 0;
    virtual void addVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)=0;
    virtual void addV(const ParticleSet& P, int iat, ValueType* restrict vals)=0;

  };

  /** Handles a set of correction orbitals per atom 
   *
   * Reduction over the orbitals - beware of the reduction problem
   */ 
  template<typename T>
  typedef CuspCorrectionAtomicBasis
  {
    tyepdef MultiQuinticSpline1D<T> RadialSetType;

    RadialSetType AOs;
    aligned_vector<size_t> ID;

    inline void evaluate(const T r, T* restrict vals) const
    {
      const size_t nr=AOs.size();
      T phi[nr];
      AOs.evaluate(r,phi);
      for(size_t i=0; i<nr; ++i)
        vals[ID[i]]+=phi[i];
    }

    inline void evaluate(const T r, T* restrict u, T* restrict du, T* restrict d2u) const
    {
      const size_t nr=AOs.size();
      T phi[nr];
      T dphi[nr];
      T d2phi[nr];
      AOs.evaluate(r,phi,phi,d2phi);
      for(size_t i=0; i<nr; ++i)
      {
        const size_t j=ID[i];
        u[j]  += phi[i];
        du[j] += dphi[i];
        d2u[j]+= d2phi[i];
      }
    }
  }

/** A localized basis set derived from BasisSetBase<typename COT::value_type>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 */
struct SoaCuspCorrection: public CuspCorrectionBase  //: public BasisSetBase<typename COT::value_type>
{
  //typedef typename OrbitalSetTraits<value_type>::VGLVector_t VGLVector_t;

  ///number of centers, e.g., ions
  size_t NumCenters;
  ///number of quantum particles
  size_t NumTargets;
  ///number of quantum particles
  int myTableIndex;
  ///size of the basis set
  int BasisSetSize;
  ///Reference to the center
  const ParticleSet::ParticleIndex_t& IonID;

  ///COMPLEX WON'T WORK
  typedef CuspCorrectionAtomicBasis<RealType> COT;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number  of unique centers
   */
  aligned_vector<COT*> LOBasisSet;

  Matrix<RealType> myVGL;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  SoaCuspCorrection(ParticleSet& ions, ParticleSet& els): IonID(ions.GroupID)
  {
    myTableIndex=els.addTable(ions,DT_SOA);
    NumCenters=ions.getTotalNum();
    NumTargets=els.getTotalNum();
    LOBasisSet.resize(ions.getSpeciesSet().getTotalNum(),0);
  }

  /** copy constructor */
  SoaCuspCorrection(const SoaCuspCorrection& a)=default;

  /** makeClone */
  SoaCuspCorrection<COT>* makeClone() const
  {
    SoaCuspCorrection<COT>* myclone=new SoaCuspCorrection<COT>(*this);
    for(int i=0; i<LOBasisSet.size(); ++i)
      myclone->LOBasisSet[i]=LOBasisSet[i]->makeClone();
    return myclone;
  }


  /** set BasisSetSize and allocate mVGL container
   */
  void setBasisSetSize(int nbs)
  {
    BasisSetSize=nbs;
    //THIS NEEDS TO BE FIXE for OpenMP
    myVGL.resize(3,getAlignedSize(BasisSetSize));
  }

  /** compute VGL 
   * @param P quantum particleset
   * @param iat active particle
   * @param vgl Matrix(5,BasisSetSize)
   * @param trialMove if true, use Temp_r/Temp_dr
   */
  inline void evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
  {

    constexpr RealType czero(0);

    myVGL=czero;

    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const auto restrict dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];
    for(int c=0; c<NumCenters; c++)
    {
      LOBasisSet[IonID[c]]->evaluate(dist[c],myVGL[0],myVGL[1],myVGL[2]);
    }

    {
      const auto restrict v_in=myVGL[0];
      const auto restrict g_in=myVGL[1];
      const auto restrict l_in=myVGL[2];
      auto restrict v_out =vgl.data(0);
      auto restrict gx_out=vgl.data(1);
      auto restrict gy_out=vgl.data(2);
      auto restrict gz_out=vgl.data(3);
      auto restrict l_out =vgl.data(4);
      for(size_t i=0; i<BasisSetSize; ++i)
      {
        v_out[i] +=v_in[i];
        gx_out[i]+=g_in[i];
        gy_out[i]+=g_in[i];
        gz_out[i]+=g_in[i];
        l_out[i] +=l_in[i];
      }
    }
  }

  /** compute values for the iat-paricle move
   *
   * Always uses Temp_r and Temp_dr
   */
  inline void evaluateV(const ParticleSet& P, int iat, ValueType* restrict vals)
  {
    constexpr RealType czero(0);
    ValueType* tmp_vals=myVGL[0];

    std::fill_n(tmp_vals,myVGL.capacity(),czero);

    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    const auto restrict dist = (P.activePtcl==iat)? d_table->Temp_r.data(): d_table->Distances[iat];

    //THIS IS SERIAL, only way to avoid this is to use myVGL 
    for(int c=0; c<NumCenters; c++)
      LOBasisSet[IonID[c]]->evaluate(dist[c],tmp_vals);

    {//collect
      const auto restrict v_in=myVGL[0];
      for(size_t i=0; i<BasisSetSize; ++i)
      {
        vals[i] +=v_in[i];
      }
    }

  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos)
  {
    LOBasisSet[icenter]=aos;
  }
};
}
#endif
