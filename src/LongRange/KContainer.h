#ifndef QMCPLUSPLUS_KCONTAINER_H
#define QMCPLUSPLUS_KCONTAINER_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/OhmmsInfo.h"
namespace qmcplusplus
{

/** Container for k-points
 *
 * It generates a set of k-points that are unit-translations of the reciprocal-space
 * cell. K-points are generated within a spherical cutoff set by the supercell
 */
class KContainer: public QMCTraits
{
private:

  //Function to return a unique number for each kVector
  inline long GetHashOfVec(const TinyVector<int,3>& inpv, int hashparam)
  {
    return inpv[2] + hashparam * (inpv[1] + hashparam * inpv[0]);
  }

  inline long GetHashOfVec(const TinyVector<int,2>& inpv, int hashparam)
  {
    return (inpv[1] + hashparam * inpv[0]);
  }

  inline long GetHashOfVec(const TinyVector<int,1>& inpv, int hashparam)
  {
    return inpv[0];
  }

  /// The cutoff up to which k-vectors are generated. 
  RealType kcutoff;
  /// kcutoff*kcutoff
  RealType kcut2; 

public:
  //Typedef for the lattice-type. We don't need the full particle-set.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  ///typedef of vector containers
  typedef vector<PosType>  VContainer_t;
  ///typedef of scalar containers
  typedef vector<RealType> SContainer_t;

  ///number of k-points
  int numk;

  /** maximum integer translations of reciprocal cell within kc.
   *
   * Last index is max. of first dimension+1
   */
  TinyVector<int,DIM+1> mmax;

  /** K-vector in reduced coordinates
   */
  vector<TinyVector<int,DIM> > kpts;
  /** K-vector in Cartesian coordinates
   */
  VContainer_t kpts_cart;
  /** squre of kpts in Cartesian coordniates
   */
  SContainer_t ksq;
  /** Given a k index, return index to -k
   */
  vector<int> minusk;
  /** kpts which belong to the ith-shell [kshell[i], kshell[i+1]) */
  vector<int> kshell;

  /** k points sorted by the |k|  excluding |k|=0
   *
   * The first for |k|
   * The second for a map to the full index. The size of the second is the degeneracy.
   */
  //std::map<int,std::vector<int>*>  kpts_sorted;

  /** update k-vectors 
   * @param sc supercell
   * @param kc cutoff radius in the K
   * @param useSphere if true, use the |K|
   */
  void UpdateKLists(ParticleLayout_t& lattice, RealType kc, bool useSphere=true);

private:
  /** compute approximate parallelpiped that surrounds kc
   * @param lattice supercell
   */
  void FindApproxMMax(ParticleLayout_t& lattice);
  /** construct the container for k-vectors */
  void BuildKLists(ParticleLayout_t& lattice, bool useSphere);
};

}

#endif
