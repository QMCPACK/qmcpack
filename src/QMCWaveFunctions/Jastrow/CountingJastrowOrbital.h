#ifndef QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H
#define QMC_PLUS_PLUS_COUNTING_JASTROW_ORBITAL_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"

namespace qmcplusplus
{

//#template <class RegionType> class CountingJastrowOrbital : public WaveFunctionComponent 
template <class RegionType>
class CountingJastrowOrbital: public WaveFunctionComponent 
{
  ///** flag to set the optimization mode */
  //bool IsOptimizing;

  ///** boolean to set optimization
  // *
  // * If true, this object is actively modified during optimization
  // */
  //bool Optimizable;

  ///** true, if FermionWF */
  //bool IsFermionWF;

  ///** true, if it is done with derivatives */
  //bool derivsDone;

  ///** true, if compute for the ratio instead of buffering */
  //bool Need2Compute4PbyP;

  ///** define the level of storage in derivative buffer **/
  //int DerivStorageType;

  ///** flag to calculate and return ionic derivatives */
  //bool ionDerivs;

  //int parameterType;

  ///** current update mode */
  //int UpdateMode;

  ///** current \f$\log\phi \f$
  // */
  //RealType LogValue;

  ///** current phase
  // */

  //RealType PhaseValue;
  ///** Pointer to the differential orbital of this object
  // *
  // * If dPsi=0, this orbital is constant with respect to the optimizable variables
  // */
  //DiffWaveFunctionComponentPtr dPsi;

  ///** A vector for \f$ \frac{\partial \nabla \log\phi}{\partial \alpha} \f$
  // */
  //GradVectorType dLogPsi;

  ///** A vector for \f$ \frac{\partial \nabla^2 \log\phi}{\partial \alpha} \f$
  // */
  //ValueVectorType d2LogPsi;

  ///** Name of the class derived from WaveFunctionComponent
  // */
  //std::string ClassName;

  /////list of variables this orbital handles
  //opt_variables_type myVars;

  /////Bytes in WFBuffer
  //size_t Bytes_in_WFBuffer;

  /////assign a differential orbital
  //virtual void setDiffOrbital(DiffWaveFunctionComponentPtr d);

  /////assembles the full value from LogValue and PhaseValue
  //ValueType getValue() const

public:
  // constructor
  CountingJastrowOrbital(ParticleSet& P)
  {

  }

  void checkInVariables(opt_variables_type& active)
  {
  }

  void checkOutVariables(const opt_variables_type& active)
  {
  }

  void resetParameters(const opt_variables_type& active)
  {
  }

  void reportStatus(std::ostream& os)
  {
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) { return RealType(0); };

  void recompute(ParticleSet& P) {};

  //void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi_all)
  //{
  //  APP_ABORT("WaveFunctionComponent::evaluateHessian is not implemented in "+ClassName+" class.");
  //}

  /** return the current gradient for the iat-th particle
   * @param Pquantum particle set
   * @param iat particle index
   * @return the gradient of the iat-th particle
   */
  GradType evalGrad(ParticleSet& P, int iat) { return GradType(0); };

  /////** return the logarithmic gradient for the iat-th particle
  //// * of the source particleset
  //// * @param Pquantum particle set
  //// * @param iat particle index
  //// * @return the gradient of the iat-th particle
  //// */
  //GradType evalGradSource(ParticleSet& P,
  //                                ParticleSet& source,
  //                                int iat);

  /////** Adds the gradient w.r.t. the iat-th particle of the
  //// *  source particleset (ions) of the logarithmic gradient
  //// *  and laplacian w.r.t. the target paritlceset (electrons).
  //// * @param P quantum particle set (electrons)
  //// * @param source classical particle set (ions)
  //// * @param iat particle index of source (ion)
  //// * @param the ion gradient of the elctron gradient
  //// * @param the ion gradient of the elctron laplacian.
  //// * @return the log gradient of psi w.r.t. the source particle iat
  //// */
  //GradType evalGradSource
  //(ParticleSet& P, ParticleSet& source, int iat,
  // TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
  // TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  /** evaluate the ratio of the new to old orbital value
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @param grad_iat Gradient for the active particle
   */
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) { return ValueType(0); };

  /** a move for iat-th particle is accepted. Update the content for the next moves
   * @param P target ParticleSet
   * @param iat index of the particle whose new position was proposed
   */
  void acceptMove(ParticleSet& P, int iat) { return; } ;

  /** a move for iat-th particle is reject. Restore to the content.
   * @param iat index of the particle whose new position was proposed
   */
  void restore(int iat) { return; };

  /** evalaute the ratio of the new to old orbital value
   *@param P the active ParticleSet
   *@param iat the index of a particle
   *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'}\})\f$
   *
   *Specialized for particle-by-particle move.
   */
  ValueType ratio(ParticleSet& P, int iat) { return ValueType(0); };

  /** For particle-by-particle move. Requests space in the buffer
   *  based on the data type sizes of the objects in this class.
   * @param P particle set
   * @param buf Anonymous storage
   */
  void registerData(ParticleSet& P, WFBufferType& buf) {};

  /** For particle-by-particle move. Put the objects of this class
   *  in the walker buffer or forward the memory cursor.
   * @param P particle set
   * @param buf Anonymous storage
   * @param fromscratch request recomputing the precision critical
   *        pieces of wavefunction from scratch
   * @return log value of the wavefunction.
   */
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false) {};

  /** For particle-by-particle move. Copy data or attach memory
   *  from a walker buffer to the objects of this class.
   *  The log value, P.G and P.L contribution from the objects
   *  of this class are also added.
   * @param P particle set
   * @param buf Anonymous storage
   */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) {};

  /** make clone
   * @param tqp target Quantum ParticleSet
   * @param deepcopy if true, make a decopy
   *
   * If not true, return a proxy class
   */
  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const { return 0; };


//  virtual void multiplyDerivsByOrbR(std::vector<RealType>& dlogpsi)
//  {
//    RealType myrat = std::exp(LogValue)*std::cos(PhaseValue);
//    for(int j=0; j<myVars.size(); j++)
//    {
//      int loc=myVars.where(j);
//      dlogpsi[loc] *= myrat;
//    }
//  };


  void finalizeOptimization() { }

  ///** evaluate the ratios of one virtual move with respect to all the particles
  // * @param P reference particleset
  // * @param ratios \f$ ratios[i]=\{{\bf R}\}\rightarrow {r_0,\cdots,r_i^p=pos,\cdots,r_{N-1}}\f$
  // */
  //void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.refPtcl
   */
  void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
  }

  // function calls to pass to differential component dPsi
  void evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
    std::vector<ValueType>& ratios, Matrix<ValueType>& dratios)
  {
    evaluateRatios(VP, ratios);
    dPsi->evaluateDerivRatios(VP, optvars, dratios);
  }
  void evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars,
    std::vector<RealType>& dlogpsi, std::vector<RealType>& dhpsioverpsi)
  {
    dPsi->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }

//  void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
//    std::vector<RealType>& dgradlogpsi) 
//  {
//    dPsi->evaluateGradDerivatives(const ParticleSet::ParticleGradient_t&G_in,
//      std::vector<RealType>& dgradlogpsi);
//  }

  bool addRegion(RegionType* CR, Matrix<RealType>* F, std::vector<RealType>* G, bool opt_CR, bool opt_G, bool opt_F);

  bool addDebug(int seqlen, int period);

};

}

#endif
