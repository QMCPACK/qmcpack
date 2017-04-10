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
    
    
#ifndef QMCPLUSPLUS_SLATERDETERMINANT_ORBITAL_H
#define QMCPLUSPLUS_SLATERDETERMINANT_ORBITAL_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalTraits.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/DiracDeterminant.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief An AntiSymmetric OrbitalBase composed of DiracDeterminants.
 *
 A SlaterDeterminant is a product of DiracDeterminants
 \f[
 S({\bf R}) = \prod_i D_i({\bf r_{first,i},\ldots,r_{last,i}}).
 \f]
 *Typically \f$ S({\bf R}) \f$ is a product of an up and down
 *DiracDeterminant such that
 \f[
 S({\bf R}) = D_{\uparrow}({\bf r_1,r_2,\ldots,r_{N_{\uparrow}}})
 D_{\downarrow}({\bf r_{N_{\uparrow +1}},\ldots,r_{N}})
 \f]
 *
 *SlaterDeterminant is a composite class which collects each of
 *the compotents in a simple manner: multiplications of determinants
 *and addition of gradients and laplacian terms in a linear pattern.
 *
 *The (S)ingle(P)article(O)rbitalSet template parameter is an
 *engine which fills in single-particle orbital terms.
 *
 *@note MultiSlaterDeterminant is a linear combination of SlaterDeterminants.
 *@todo Use BasisSet for particle-by-particle update
 */
template<class SPOSet>
class SlaterDeterminant: public OrbitalBase
{

public:

  typedef DiracDeterminant<SPOSet> Determinant_t;
  typedef typename SPOSet::BasisSet_t BasisSet_t;

  /// constructor
  SlaterDeterminant():BasisSet(0)
  {
    M.resize(3,0);
    Optimizable=false;
  }

  ///destructor
  ~SlaterDeterminant() { }

  ///add a new DiracDeterminant to the list of determinants
  void add(Determinant_t* det)
  {
    int last=Dets.size();
    Dets.push_back(det);
    M[last+1]=M[last]+Dets[last]->rows();
    DetID.insert(DetID.end(),det->rows(),last);
  }

  ///reset all the Dirac determinants, Optimizable is true
  ///optimizations  are disabled
  inline void checkInVariables(opt_variables_type& active)
  { }

  inline void checkOutVariables(const opt_variables_type& active)
  { }

  inline void resetParameters(const opt_variables_type& active)
  { }

  inline void reportStatus(std::ostream& os)
  { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    //BasisSet->resetTargetParticleSet(P);
    LOGMSG("\nSlaterDeterminant::resetTargetParticleSet")
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(P);
  }

  void setBasisSet(BasisSet_t* bs)
  {
    BasisSet=bs;
  }
  /** Calculate the value of the Slater determinant for the input configuration.
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@return SlaterDeterminant value
   *
   aAdd the gradient and laplacian contribution of the Slater determinant to G(radient) and L(aplacian)
   *for local energy calculations.
   */
  inline ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L)
  {
    ValueType psi = 1.0;
    for(int i=0; i<Dets.size(); i++)
      psi *= Dets[i]->evaluate(P,G,L);
    return psi;
  }


  inline RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L)
  {
    //@attention BasisSet::evaluate is to be called but the due to the bugs, it is commented out.
    //if(BasisSet == 0)
    //{
    //  ERRORMSG("SlaterDeterminant::BasisSet is not assigned")
    //  OHMMS::Controller->abort();
    //}
    //BasisSet->evaluate(P);
    ValueType psi = 1.0;
    for(int i=0; i<Dets.size(); i++)
      psi *= Dets[i]->evaluate(P,G,L);
    return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  }

  ///return the total number of Dirac determinants
  inline int size() const
  {
    return Dets.size();
  }

  ///return the dimension of the i-th Dirac determinant
  inline int size(int i) const
  {
    return Dets[i]->cols();
  }

  /** similar to evaluateLog
   */
  RealType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    //BasisSet->evaluate(P);
    ValueType psi = 1.0;
    for(int i=0; i<Dets.size(); i++)
      psi *= Dets[i]->registerData(P,buf);
    return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  }

  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
                        bool fromscratch=false)
  {
    ValueType psi = 1.0;
    for(int i=0; i<Dets.size(); i++)
      psi *= Dets[i]->updateBuffer(P,buf);
    return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  }

  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->copyFromBuffer(P,buf);
  }

  /** reimplements the virtual function
   *
   * The DiractDeterminants of SlaterDeterminant need to save the inverse
   * of the determinant matrix to evaluate ratio
   */
  void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->dumpToBuffer(P,buf);
  }

  /** reimplements the virtual function
   *
   * Matching function to dumpToBuffer.
   */
  void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->dumpFromBuffer(P,buf);
  }

  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
  {
    //LogValue=0.0;
    //int sign_tot = 1;
    //for(int i=0; i<Dets.size(); i++)
    //{
    //  LogValue += Dets[i]->evaluateLog(P,buf);
    //  sign_tot *=Dets[i]->DetSign;
    //}
    //PhaseValue=evaluatePhase(sign_tot);
    //return LogValue;
    ValueType r=1.0;
    for(int i=0; i<Dets.size(); i++)
      r *= Dets[i]->evaluateLog(P,buf);
    return evaluateLogAndPhase(r,PhaseValue);
    //return r;
  }

  inline ValueType ratio(ParticleSet& P, int iat,
                         ParticleSet::ParticleGradient_t& dG,
                         ParticleSet::ParticleLaplacian_t& dL)
  {
    return Dets[DetID[iat]]->ratio(P,iat,dG,dL);
  }

  inline ValueType logRatio(ParticleSet& P, int iat,
                            ParticleSet::ParticleGradient_t& dG,
                            ParticleSet::ParticleLaplacian_t& dL)
  {
    ValueType r = Dets[DetID[iat]]->ratio(P,iat,dG,dL);
    return evaluateLogAndPhase(r,PhaseValue);
  }

  inline void restore(int iat)
  {
    return Dets[DetID[iat]]->restore(iat);
  }

  inline void acceptMove(ParticleSet& P, int iat)
  {
    Dets[DetID[iat]]->acceptMove(P,iat);
  }

  ValueType
  ratio(ParticleSet& P, int iat)
  {
    return Dets[DetID[iat]]->ratio(P,iat);
  }

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat)
  {
    return Dets[DetID[iat]]->update(P,dG,dL,iat);
  }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    SlaterDeterminant<SPOSet>* myclone= new SlaterDeterminant<SPOSet>(*this);
    //myclone->M=M;
    //myclone->DetID=DetID;
    //myclone->Dets.resize(Dets.size(),0);
    for(int i=0; i<Dets.size(); i++)
    {
      Determinant_t* adet=new Determinant_t(*Dets[i]);
      adet->Phi=Dets[i]->clonePhi();//assign a new SPOSet
      adet->resetTargetParticleSet(tqp);
      myclone->Dets[i]=adet;
    }
    return myclone;
  }

private:
  std::vector<int> M;
  std::vector<int> DetID;
  ///container for the DiracDeterminants
  std::vector<Determinant_t*>  Dets;
  BasisSet_t* BasisSet;
};
}
#endif
