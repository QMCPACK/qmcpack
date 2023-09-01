//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of
// Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National
//                    Laboratory Jeongnim Kim, jeongnim.kim@gmail.com,
//                    University of Illinois at Urbana-Champaign Mark A.
//                    Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_BASISSETFACTORYT_H
#define QMCPLUSPLUS_BASISSETFACTORYT_H

#include "QMCWaveFunctions/SPOSetBuilderT.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
template <typename T>
class SPOSetBuilderFactoryT : public MPIObjectBase
{
public:
    using SPOMap = typename SPOSetT<T>::SPOMap;
    using PSetMap =
        std::map<std::string, const std::unique_ptr<ParticleSetT<T>>>;

    /** constructor
     * \param comm communicator
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
    SPOSetBuilderFactoryT(
        Communicate* comm, ParticleSetT<T>& els, const PSetMap& psets);

    ~SPOSetBuilderFactoryT();

    std::unique_ptr<SPOSetBuilderT<T>>
    createSPOSetBuilder(xmlNodePtr rootNode);

    /** returns a named sposet from the pool
     *  only use in serial portion of execution
     *  ie during initialization prior to threaded code
     */
    const SPOSetT<T>*
    getSPOSet(const std::string& name) const;

    void
    buildSPOSetCollection(xmlNodePtr cur);

    bool
    empty() const
    {
        return sposets.empty();
    }

    /** add an SPOSet to sposets map.
     * This is only used to handle legacy SPOSet input styles without using
     * sposet_collection
     */
    void addSPOSet(std::unique_ptr<SPOSetT<T>>);

    SPOMap&&
    exportSPOSets()
    {
        return std::move(sposets);
    }

private:
    /// reference to the target particle
    ParticleSetT<T>& targetPtcl;

    /// reference to the particle pool
    const PSetMap& ptclPool;

    /// list of all sposets created by the builders of this factory
    SPOMap sposets;

    static std::string basisset_tag;
};
} // namespace qmcplusplus
#endif
