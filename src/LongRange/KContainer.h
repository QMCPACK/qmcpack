#ifndef QMCPLUSPLUS_KCONTAINER_H
#define QMCPLUSPLUS_KCONTAINER_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  /** @ingroup longrange
   *\brief A handler class for generating and storing lists of
   * k-points that are unit-translations of the reciprocal-space
   * cell. K-points are generated within a spherical cutoff, kc
   */

  class KContainer: public QMCTraits {
  private:
    //Typedef for the lattice-type. We don't need the full particle-set.
    typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
    ///typedef of vector containers
    typedef vector<TinyVector<RealType,3> > VContainer_t;
    ///typedef of scalar containers
    typedef vector<RealType>                SContainer_t;

    //Function to return a unique number for each kVector
    long GetHashOfVec(const TinyVector<int,3>& inpv, int hashparam) {
      return inpv[2] + hashparam * (inpv[1] + hashparam * inpv[0]);
    }
    
    //The cutoff up to which k-vectors are generated.
    RealType kcutoff;
    RealType kcut2; //kcutoff*kcutoff

  public:
    //Maximum integer translations of reciprocal cell within kc.
    //Last index is max. of first 3.
    TinyVector<int,4> mmax;

    //K-vector list
    /** K-vector in reduced coordinates 
     */
    vector<TinyVector<int,3> > kpts; 
    /** K-vector in Cartesian coordinates 
     */
    VContainer_t kpts_cart; 
    /** squre of kpts in Cartesian coordniates 
     */
    SContainer_t ksq; 
    /** Given a k index, return index to -k
     */
    vector<int> minusk; 
    int numk;

    //A copy of the lattice, so that we have the cell-vectors
    ParticleLayout_t& Lattice;
    
  public:
    //Constructor
    KContainer(ParticleLayout_t& ref);
    //Copy Constructor
    //    KContainer(const KContainer& ref);
    //Destructor
    ~KContainer();
    //Overloaded assignment operator
    KContainer& operator=(const KContainer&);

    //Public Methods:
    // UpdateKLists() - call for new k or when lattice changed.
    void UpdateKLists(ParticleLayout_t& ref, RealType kc, bool useSphere=true);
    void UpdateKLists(RealType kc, bool useSphere=true);

  private:
    //Private Methods:
    // FindApproxMMax - compute approximate parallelpiped that surrounds kc
    // BuildKLists - Correct mmax and fill lists of k-vectors.
    void FindApproxMMax();
    void BuildKLists(bool useSphere);
  };

}

#endif
