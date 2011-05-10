#ifndef QMCPLUSPLUS_STRUCTFACT_H
#define QMCPLUSPLUS_STRUCTFACT_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "LongRange/KContainer.h"
#include "OhmmsPETE/OhmmsVector.h"

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
    typedef PooledData<RealType>  BufferType;

  public:
    /** false, if the structure factor is not actively update
     *
     * Default is false. Particle-by-particle update functions, makeMove, acceptMove and rejectMove
     * are costly and do not need to be performed unless Hamiltonian uses pbyp.
     */
    bool DoUpdate;
    ///reference particle set
    ParticleSet& PtclRef;
    /** Rhok[alpha][k]
     *
     * Structure factor Rhok[alpha][k] \f$ \equiv \rho_{k}^{\alpha} = \sum_{i}e^{i{\bf k}\cdot{\bf r_i}}\f$
     */
    Matrix<ComplexType> rhok;
    ///eikr[particle-index][K]
    Matrix<ComplexType> eikr;
    /// K-Vector List.
    KContainer KLists;
    ///eikr[K] for a proposed move
    Vector<ComplexType> eikr_temp;
    /** Constructor - copy ParticleSet and init. k-shells
     * @param ref Reference particle set
     * @param kc cutoff for k
     */
    StructFact(ParticleSet& ref, RealType kc);
    ////Copy Constructor
    //StructFact(const StructFact& ref);
    //Destructor
    ~StructFact();
    //Need to overload assignment operator.
    //Default doesn't work because we have non-static reference members.
    StructFact& operator=(const StructFact& ref);
      
    /** Recompute Rhok if lattice changed
     * @param kc cut-off K
     */
    void UpdateNewCell(RealType kc);
    /**  Update Rhok if all particles moved
     */
    void UpdateAllPart();

    /** evaluate eikr_temp for eikr for the proposed move
     * @param active index of the moved particle
     * @param pos proposed position
     */
    void makeMove(int active, const PosType& pos);
    /** update eikr and rhok with eikr_temp
     * @param active index of the moved particle  
     */
    void acceptMove(int active);
    /** discard any temporary data
     * @param active index of the moved particle  
     *
     * Do nothing
     */
    void rejectMove(int active);
    /// Update Rhok if 1 particle moved
    //void Update1Part(const PosType& rold, const PosType& rnew,int iat,int GroupID);

    //Buffer methods. For PbyP MC where data must be stored in an anonymous
    //buffer between iterations
    /** @brief register rhok data to buf so that it can copyToBuffer and copyFromBuffer
     *
     * This function is used for particle-by-particle MC methods to register structure factor
     * to an anonymous buffer.
     */
    inline void registerData(BufferType& buf) {
      buf.add(rhok.first_address(),rhok.last_address());
      buf.add(eikr.first_address(),eikr.last_address());
    }

    /** @brief put rhok data to buf
     *
     * This function is used for particle-by-particle MC methods
     */
    inline void updateBuffer(BufferType& buf) {
      buf.put(rhok.first_address(),rhok.last_address());
      buf.put(eikr.first_address(),eikr.last_address());
    }
    /** @brief copy the data to an anonymous buffer
     *
     * Any data that will be used by the next iteration will be copied to a buffer.
     */
    inline void copyToBuffer(BufferType& buf) {
      buf.put(rhok.first_address(),rhok.last_address());
      buf.put(eikr.first_address(),eikr.last_address());
    }
    /** @brief copy the data from an anonymous buffer
     *
     * Any data that was used by the previous iteration will be copied from a buffer.
     */
    inline void copyFromBuffer(BufferType& buf) {
      buf.get(rhok.first_address(),rhok.last_address());
      buf.get(eikr.first_address(),eikr.last_address());
    }

  private:
    ///data for recursive evaluation for a given position
    Matrix<ComplexType> C;
    ///Compute all rhok elements from the start
    void FillRhok();
    ///Smart update of rhok for 1-particle move. Simply supply old+new position
    void UpdateRhok(const PosType& rold,
        const PosType& rnew,int iat,int GroupID);
    ///resize the internal data
    void resize();
  };
}

#endif
