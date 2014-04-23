#ifndef QMCPLUSPLUS_LRBASIS_H
#define QMCPLUSPLUS_LRBASIS_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus
{

/** @ingroup longrange
 *\brief Base-class for long-range breakups.
 *
 * Contains 3 important functions: c(n,k), h(n,r), hintr2(r)
 *  which evaluate the n'th basis function at k or r.
 *  \f[ hintr2_n = \int dr h_n(r) r^2, \f]
 */

class LRBasis: public QMCTraits
{
protected:
  ///size of the basis elements
  int BasisSize;
  ///Real-space cutoff for short-range part
  RealType m_rc;

  ///Typedef for the lattice-type. We don't need the full particle-set.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;

  //A copy of the lattice, so that we have the cell-vectors + volume.
  ParticleLayout_t& Lattice;

public:

  //Constructor: copy lattice and init rc
  LRBasis(ParticleLayout_t& ref): Lattice(ref), m_rc(0.0)
  {
    /*Do nothing*/
  }

  virtual ~LRBasis() { }

  inline int NumBasisElem() const
  {
    return BasisSize;
  }

  //Real-space basis function + integral: override these
  virtual RealType h(int n, RealType r) const = 0;
  virtual RealType hintr2(int n) = 0;
  virtual RealType dh_dr(int n, RealType r) const {return 0.0;};
  //k-space basis function: override this
  virtual RealType c(int m, RealType k) = 0;
  //k-space basis function k space derivative.
  virtual RealType dc_dk(int m, RealType k) {return 0.0;};
  
//
// df(m,r) is included for legacy reasons.  Please use dh_dr
//

  inline RealType df(int m, RealType r) const { return dh_dr(m,r); };
  
  ///$f(r,{tn})$ returns the value of $\sum_n t_n*h_{\alpha n}(r)$ 
  ///   r is radial position (scalar)
  ///   vector<RealType> coefs are the {tn} optimized breakup coefficients.
  
/*
 * 
 * name: f: function interpolated by optimized-breakup.
 * @param r The radial coordinate.
 * @param coefs A vector of optimized breakup coefficents.
 * @return The value of $f(r) = \sum_n t_n h_{n \alpha}(r)$
 * 
 */
  
  inline RealType f(RealType r, vector<RealType>& coefs)
  {
    RealType f = 0.0;
    //RealType df = myFunc.df(r, rinv);
    for(int n=0; n<coefs.size(); n++)
      f += coefs[n]*h(n,r);
    return f;	  
  }
  
/*
 * 
 * name: df_dr
 * @param r The radial coordinate.
 * @param coefs A vector of optimized breakup coefficents.
 * @return $\frac{df}{dr}$.
 * 
 */
  
  inline RealType df_dr(RealType r, vector<RealType>& coefs)
  {
    RealType df = 0.0;
    //RealType df = myFunc.df(r, rinv);
    for(int n=0; n<coefs.size(); n++)
      df += coefs[n]*dh_dr(n,r);
    return df;
  }
  
/*
 * 
 * name: fk
 * @param k |k|
 * @param coefs Optimized breakup coefficients
 * @return The fourier transform of $f(r)$
 * 
 */
  
  
  inline RealType fk(RealType k, vector<RealType> coefs)
  {
	RealType fk = 0.0;
	for (int n=0; n<coefs.size(); n++)
	  fk += coefs[n]*c(n,k);
	return fk;	  
  }
  
/*
 * 
 * name: dfk_dk
 * @param k |k|
 * @param coefs Optimized breakup coefficients
 * @return $\frac{df_k}{dk}$
 * 
 */
  
  inline RealType dfk_dk(RealType k, vector<RealType> coefs)
  {
	RealType dfk = 0.0;
	for (int n=0; n<coefs.size(); n++)
	  dfk += coefs[n]*dc_dk(n,k);
	return dfk;
  }
    
  //May need extra functionality when resetting rc. Override this.
  virtual void set_rc(RealType rc) = 0;
  inline RealType get_rc()
  {
    return m_rc;
  }
  inline RealType get_CellVolume()
  {
    return Lattice.Volume;
  }
  inline ParticleLayout_t& get_Lattice()
  {
    return Lattice;
  }
  inline void set_Lattice(ParticleLayout_t& ref)
  {
    Lattice = ref;
  }
};

}

#endif
