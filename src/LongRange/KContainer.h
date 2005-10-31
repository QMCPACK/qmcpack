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

    //The cutoff up to which k-vectors are generated.
    RealType kcutoff;
    RealType kcut2; //kcutoff*kcutoff

  public:

    //Maximum integer translations of reciprocal cell within kc.
    //Last index is max. of first 3.
    TinyVector<int,4> mmax;

    //K-vector list
    vector<TinyVector<int,3> > kpts; //In reduced coordinates
    vector<TinyVector<RealType,3> > kpts_cart; //In Cartesian coordinates
    vector<int> minusk; //Given a k index, return index to -k.
    int numk;

    //A copy of the lattice, so that we have the cell-vectors
    ParticleLayout_t& Lattice;
    
  public:
    //Constructor
    KContainer(ParticleLayout_t& ref);
    //Destructor
    ~KContainer();

    //Public Methods:
    // UpdateKLists() - call for new k or when lattice changed.
    void UpdateKLists(ParticleLayout_t& ref, RealType kc);

  private:
    //Private Methods:
    // FindApproxMMax - compute approximate parallelpiped that surrounds kc
    // BuildKLists - Correct mmax and fill lists of k-vectors.
    void FindApproxMMax();
    void BuildKLists();
  };

}

#endif
