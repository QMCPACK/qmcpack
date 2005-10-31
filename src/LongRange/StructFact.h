#ifndef QMCPLUSPLUS_STRUCTFACT_H
#define QMCPLUSPLUS_STRUCTFACT_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "LongRange/KContainer.h"
#include <complex>

namespace qmcplusplus {

  /** @ingroup longrange
   *\brief Calculates the structure-factor for a particle set
   *
   * \f[ rho_{\bf k} = \sum_{i} e^{i{\bf k}.{\bf r_i}}, \f]
   *
   * The structure-factor is species-dependent - access with GroupID
   */

  class StructFact: public QMCTraits {
    //Type is hard-wired. Rhok must be complex
    //First index: k-point (ordered as per kpts and kptscart)
    //Second index: GroupID from PtclRef.
  private:
    typedef TinyVector<OHMMS_PRECISION,3> Position_t;

  public:
    Matrix<complex<RealType> > rhok;
    //Maximum reciprocal cell translations in kc.
    //Last index is max. of first 3.
    //    TinyVector<int,4> mmax; 
    ParticleSet& PtclRef;
    //K-Vector List.
    //    vector<TinyVector<int,3> > kpts;
    //    vector<TinyVector<RealType,3> > kpts_cart;
    KContainer KLists;
    
  public:
    //Constructor - copy ParticleSet and init. k-shells
    StructFact(ParticleSet& ref, RealType kc);
    //Destructor
    ~StructFact();
      
    //Public Methods:
    ///Recompute Rhok if lattice changed
    void UpdateNewCell(RealType kc);
    /// Update Rhok if 1 particle moved
    void Update1Part(Position_t rold,Position_t rnew,int GroupID);
    /// Update Rhok if all particles moved
    void UpdateAllPart();

  private:
    //Private Methods
    ///Compute all rhok elements from the start
    void FillRhok();
    ///Smart update of rhok for 1-particle move. Simply supply old+new position
    void UpdateRhok(Position_t rold,Position_t rnew,int GroupID);

  };
}

#endif
