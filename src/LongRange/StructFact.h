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
    typedef PooledData<RealType>  BufferType;

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
    //Copy Constructor
    StructFact(const StructFact& ref);
    //Destructor
    ~StructFact();
    //Need to overload assignment operator.
    //Default doesn't work because we have non-static reference members.
    StructFact& operator=(const StructFact& ref);
      
    //Public Methods:
    ///Recompute Rhok if lattice changed
    void UpdateNewCell(RealType kc);
    /// Update Rhok if 1 particle moved
    void Update1Part(Position_t rold,Position_t rnew,int GroupID);
    /// Update Rhok if all particles moved
    void UpdateAllPart();

    //Buffer methods. For PbyP MC where data must be stored in an anonymous
    //buffer between iterations
    /** @brief register rhok data to buf so that it can copyToBuffer and copyFromBuffer
     *
     * This function is used for particle-by-particle MC methods to register structure factor
     * to an anonymous buffer.
     */
    inline void registerData(BufferType& buf) {
      //Buffertype is presently always real. Rhok is complex. Trick to force saving
      RealType* first = static_cast<RealType*>(&(rhok(0,0).real()));
      RealType* last = first + rhok.size()*2;
      buf.add(first,last);
    };
    /** @brief register rhok data to buf so that it can copyToBuffer and copyFromBuffer
     *
     * This function is used for particle-by-particle MC methods to register structure factor
     * to an anonymous buffer.
     */
    inline void updateBuffer(BufferType& buf) {
      //Buffertype is presently always real. Rhok is complex. Trick to force saving
      RealType* first = static_cast<RealType*>(&(rhok(0,0).real()));
      RealType* last = first + rhok.size()*2;
      buf.put(first,last);
    };
    /** @brief copy the data to an anonymous buffer
     *
     * Any data that will be used by the next iteration will be copied to a buffer.
     */
    inline void copyToBuffer(BufferType& buf) {
      //Buffertype is presently always real. Rhok is complex. Trick to force saving
      RealType* first = static_cast<RealType*>(&(rhok(0,0).real()));
      RealType* last = first + rhok.size()*2;
      buf.put(first,last);
    };
    /** @brief copy the data from an anonymous buffer
     *
     * Any data that was used by the previous iteration will be copied from a buffer.
     */
    inline void copyFromBuffer(BufferType& buf) {
      //Buffertype is presently always real. Rhok is complex. Trick to force saving
      RealType* first = static_cast<RealType*>(&(rhok(0,0).real()));
      RealType* last = first + rhok.size()*2;
      buf.get(first,last);
    };

  private:
    //Private Methods
    ///Compute all rhok elements from the start
    void FillRhok();
    ///Smart update of rhok for 1-particle move. Simply supply old+new position
    void UpdateRhok(Position_t rold,Position_t rnew,int GroupID);

  };
}

#endif
