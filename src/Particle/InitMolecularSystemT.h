//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of
// Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_INITMOLECULARSYSTEMT_H
#define QMCPLUSPLUS_INITMOLECULARSYSTEMT_H

#include "OhmmsData/OhmmsElementBase.h"
#include "ParticleSetTraits.h"

#include <map>

namespace qmcplusplus
{
template <typename T>
class ParticleSetT;
template <typename T>
class ParticleSetPoolT;

/* Engine to initialize the initial electronic structure for a molecular system
 */
template <typename T>
class InitMolecularSystemT : public OhmmsElementBase
{
public:
    using RealType = typename ParticleSetTraits<T>::RealType;

    InitMolecularSystemT(ParticleSetPoolT<T>& pset, const char* aname = "mosystem");

    bool
    get(std::ostream& os) const override;
    bool
    put(std::istream& is) override;
    bool
    put(xmlNodePtr cur) override;
    void
    reset() override;

    /** initialize els for an atom
     */
    void
    initAtom(ParticleSetT<T>* ions, ParticleSetT<T>* els);
    /** initialize els position for a molecule
     *
     * Use the valence of each ionic species on a sphere
     */
    void
    initMolecule(ParticleSetT<T>* ions, ParticleSetT<T>* els);
    /** initialize els for the systems with a mixed boundary
     *
     * Use the bound of the ionic systems and uniform random positions within a
     * reduced box
     */
    void
    initWithVolume(ParticleSetT<T>* ions, ParticleSetT<T>* els);

private:
    /** pointer to ParticleSetPool
     *
     * QMCHamiltonian needs to know which ParticleSet object
     * is used as an input object for the evaluations.
     * Any number of ParticleSet can be used to describe
     * a QMCHamiltonian.
     */
    ParticleSetPoolT<T>& ptclPool;
};
} // namespace qmcplusplus
#endif
