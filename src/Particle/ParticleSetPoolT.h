//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of
// Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign Mark A. Berrill,
//                    berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PARTICLESETPOOLT_H
#define QMCPLUSPLUS_PARTICLESETPOOLT_H

#include "Message/MPIObjectBase.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/MCWalkerConfigurationT.h"
#include "ParticleSetT.h"
#include "SimulationCellT.h"

namespace qmcplusplus
{
/** @ingroup qmcapp
 * @brief Manage a collection of ParticleSet objects
 *
 * This object handles \<particleset\> elements and
 * functions as a builder class for ParticleSet objects.
 */
template <typename T>
class ParticleSetPoolT : public MPIObjectBase
{
public:
    using PoolType =
        std::map<std::string, const std::unique_ptr<ParticleSetT<T>>>;

    /** constructor
     * @param aname xml tag
     */
    ParticleSetPoolT(Communicate* c, const char* aname = "particleset");
    ~ParticleSetPoolT();

    ParticleSetPoolT(const ParticleSetPoolT&) = delete;
    ParticleSetPoolT&
    operator=(const ParticleSetPoolT&) = delete;
    ParticleSetPoolT(ParticleSetPoolT&& pset) noexcept;
    ParticleSetPoolT&
    operator=(ParticleSetPoolT&&) = default;

    bool
    put(xmlNodePtr cur);
    bool
    get(std::ostream& os) const;
    void
    reset();

    void
    output_particleset_info(Libxml2Document& doc, xmlNodePtr root);

    /** initialize the supercell shared by all the particle sets
     *
     *  return value is never checked anywhere
     *  side effect simulation_cell_ UPtr<ParticleLayout> is set
     *  to particle layout created on heap.
     *  This is later directly assigned to pset member variable Lattice.
     */
    bool
    readSimulationCellXML(xmlNodePtr cur);

    /// return true, if the pool is empty
    inline bool
    empty() const
    {
        return myPool.empty();
    }

    /** add a ParticleSet* to the pool with its ownership transferred
     * ParticleSet built outside the ParticleSetPool must be constructed with
     * the simulation cell from this->simulation_cell_.
     */
    void
    addParticleSet(std::unique_ptr<ParticleSetT<T>>&& p);

    /** get a named ParticleSet
     * @param pname name of the ParticleSet
     * @return a MCWalkerConfiguration object with pname
     *
     * When the named ParticleSet is not in this object, return 0.
     */
    ParticleSetT<T>*
    getParticleSet(const std::string& pname);

    /** get a named MCWalkerConfiguration
     * @param pname name of the MCWalkerConfiguration
     * @return a MCWalkerConfiguration object with pname
     *
     * When the named MCWalkerConfiguration is not in this object, return 0.
     */
    MCWalkerConfigurationT<T>*
    getWalkerSet(const std::string& pname);

    /** get the Pool object
     */
    inline const PoolType&
    getPool() const
    {
        return myPool;
    }

    /// get simulation cell
    const auto&
    getSimulationCell() const
    {
        return *simulation_cell_;
    }

    /// set simulation cell
    void
    setSimulationCell(const SimulationCellT<T>& simulation_cell)
    {
        *simulation_cell_ = simulation_cell;
    }

    /** randomize a particleset particleset/@random='yes' &&
     * particleset@random_source exists
     */
    void
    randomize();

private:
    /** global simulation cell
     *
     * updated by
     * - readSimulationCellXML() parsing <simulationcell> element
     * - setSimulationCell()
     */
    std::unique_ptr<SimulationCellT<T>> simulation_cell_;
    /** List of ParticleSet owned
     *
     * Each ParticleSet has to have a unique name which is used as a key for the
     * map.
     */
    PoolType myPool;
    /** xml node for random initialization.
     *
     * randomize() process initializations just before starting qmc sections
     */
    std::vector<xmlNodePtr> randomize_nodes;
};
} // namespace qmcplusplus
#endif
