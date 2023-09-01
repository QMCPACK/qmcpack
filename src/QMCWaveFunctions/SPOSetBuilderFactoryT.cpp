//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at
//                    Urbana-Champaign Miguel Morales, moralessilva2@llnl.gov,
//                    Lawrence Livermore National Laboratory Jeremy McMinnis,
//                    jmcminis@gmail.com, University of Illinois at
//                    Urbana-Champaign Jaron T. Krogel, krogeljt@ornl.gov, Oak
//                    Ridge National Laboratory Jeongnim Kim,
//                    jeongnim.kim@gmail.com, University of Illinois at
//                    Urbana-Champaign Mark A. Berrill, berrillma@ornl.gov, Oak
//                    Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "SPOSetBuilderFactoryT.h"

#include "ModernStringUtils.hpp"
#include "QMCWaveFunctions/ElectronGas/FreeOrbitalBuilderT.h"
#include "QMCWaveFunctions/HarmonicOscillator/SHOSetBuilderT.h"
#include "QMCWaveFunctions/SPOSetScannerT.h"
#if OHMMS_DIM == 3
#include "QMCWaveFunctions/LCAO/LCAOSpinorBuilderT.h"
#include "QMCWaveFunctions/LCAO/LCAOrbitalBuilderT.h"
#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#endif

#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#endif
#endif
#include "Message/MPIObjectBase.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/CompositeSPOSetT.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
template <typename T>
struct LCAOSpinorBuilderMaker
{
    template <typename... TArgs>
    std::unique_ptr<LCAOSpinorBuilderT<T>>
    operator()(TArgs&&...) const
    {
        throw std::runtime_error(
            "lcao spinors not compatible with non-complex value types");
    }
};

template <typename T>
struct LCAOSpinorBuilderMaker<std::complex<T>>
{
    template <typename... TArgs>
    std::unique_ptr<LCAOSpinorBuilderT<std::complex<T>>>
    operator()(TArgs&&... args) const
    {
        return std::make_unique<LCAOSpinorBuilderT<std::complex<T>>>(
            std::forward<TArgs>(args)...);
    }
};

template <typename T>
const SPOSetT<T>*
SPOSetBuilderFactoryT<T>::getSPOSet(const std::string& name) const
{
    if (auto spoit = sposets.find(name); spoit == sposets.end()) {
        // keep this commented until legacy input styles are moved.
        // In legacy input styles, this look up may fail and need to build
        // SPOSetT on the fly.
        return nullptr;
    }
    else
        return spoit->second.get();
}

/** constructor
 * \param els reference to the electrons
 * \param psi reference to the wavefunction
 * \param ions reference to the ions
 */
template <typename T>
SPOSetBuilderFactoryT<T>::SPOSetBuilderFactoryT(
    Communicate* comm, ParticleSetT<T>& els, const PSetMap& psets) :
    MPIObjectBase(comm),
    targetPtcl(els),
    ptclPool(psets)
{
    ClassName = "SPOSetBuilderFactoryT";
}

template <typename T>
SPOSetBuilderFactoryT<T>::~SPOSetBuilderFactoryT()
{
    DEBUG_MEMORY("SPOSetBuilderFactoryT::~SPOSetBuilderFactoryT");
}

template <typename T>
std::unique_ptr<SPOSetBuilderT<T>>
SPOSetBuilderFactoryT<T>::createSPOSetBuilder(xmlNodePtr rootNode)
{
    ReportEngine PRE(ClassName, "createSPOSetBuilder");
    std::string sourceOpt("ion0");
    std::string type("");
    std::string name("");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(sourceOpt, "source");
    aAttrib.add(type, "type");
    aAttrib.add(name, "name");

    if (rootNode != NULL)
        aAttrib.put(rootNode);

    std::string type_in = type;
    type = lowerCase(type);

    // when name is missing, type becomes the input
    if (name.empty())
        name = type_in;

    std::unique_ptr<SPOSetBuilderT<T>> bb;

    if (type == "composite") {
        app_log() << "Composite SPO set with existing SPOSets." << std::endl;
        bb = std::make_unique<CompositeSPOSetBuilderT<T>>(myComm, *this);
    }
    else if (type == "jellium" || type == "heg" || type == "free") {
        app_log() << "Free-particle SPO set" << std::endl;
        bb = std::make_unique<FreeOrbitalBuilderT<T>>(
            targetPtcl, myComm, rootNode);
    }
    else if (type == "sho") {
        app_log() << "Harmonic Oscillator SPO set" << std::endl;
        bb = std::make_unique<SHOSetBuilderT<T>>(targetPtcl, myComm);
    }
#if OHMMS_DIM == 3
    else if (type.find("spline") < type.size()) {
        if (targetPtcl.isSpinor()) {
#ifdef QMC_COMPLEX
            app_log() << "Einspline Spinor Set\n";
            // FIXME
            // bb = std::make_unique<EinsplineSpinorSetBuilder>(targetPtcl,
            // ptclPool, myComm, rootNode);
#else
            PRE.error("Use of einspline spinors requires QMC_COMPLEX=1.  "
                      "Rebuild with this option");
#endif
        }
        else {
#if defined(HAVE_EINSPLINE)
            PRE << "EinsplineSetBuilder:  using libeinspline for B-spline "
                   "orbitals.\n";
            // FIXME
            // bb = std::make_unique<EinsplineSetBuilder>(targetPtcl, ptclPool,
            // myComm, rootNode);
#else
            PRE.error("Einspline is missing for B-spline orbitals", true);
#endif
        }
    }
    else if (type == "molecularorbital" || type == "mo") {
        ParticleSetT<T>* ions = nullptr;
        // initialize with the source tag
        auto pit(ptclPool.find(sourceOpt));
        if (pit == ptclPool.end())
            PRE.error("Missing basisset/@source.", true);
        else
            ions = pit->second.get();
        if (targetPtcl.isSpinor()) {
            try {
                bb = LCAOSpinorBuilderMaker<T>{}(
                    targetPtcl, *ions, myComm, rootNode);
            }
            catch (const std::exception& e) {
                PRE.error(e.what());
            }
        }
        else
            bb = std::make_unique<LCAOrbitalBuilderT<T>>(
                targetPtcl, *ions, myComm, rootNode);
    }
#endif // OHMMS_DIM==3
    PRE.flush();

    if (!bb)
        myComm->barrier_and_abort("SPOSetBuilderFactoryT::createSPOSetBuilder "
                                  "SPOSetBuilderT creation failed.");

    app_log() << "  Created SPOSetT builder named '" << name << "' of type "
              << type << std::endl;
    return bb;
}

template <typename T>
void
SPOSetBuilderFactoryT<T>::buildSPOSetCollection(xmlNodePtr cur)
{
    std::string collection_name;
    std::string collection_type;
    OhmmsAttributeSet attrib;
    attrib.add(collection_name, "name");
    attrib.add(collection_type, "type");
    attrib.put(cur);

    // use collection_type as collection_name if collection_name is not given
    if (collection_name.empty())
        collection_name = collection_type;

    app_summary() << std::endl;
    app_summary() << "   Single particle orbitals (SPO) collection"
                  << std::endl;
    app_summary() << "   -----------------------------------------"
                  << std::endl;
    app_summary() << "    Name: " << collection_name
                  << "   Type input: " << collection_type << std::endl;
    app_summary() << std::endl;

    // create the SPOSetT builder
    auto bb = createSPOSetBuilder(cur);

    // going through a list of sposet entries
    int nsposets = 0;
    processChildren(
        cur, [&](const std::string& cname, const xmlNodePtr element) {
            if (cname == "sposet") {
                addSPOSet(
                    std::unique_ptr<SPOSetT<T>>(bb->createSPOSet(element)));
                nsposets++;
            }
            if (cname == "rotated_sposet") {
                addSPOSet(std::unique_ptr<SPOSetT<T>>(
                    bb->createRotatedSPOSet(element)));
                nsposets++;
            }
        });

    if (nsposets == 0)
        myComm->barrier_and_abort(
            "SPOSetBuilderFactoryT::buildSPOSetCollection  no <sposet/> "
            "elements found");

    // going through a list of spo_scanner entries
    processChildren(
        cur, [&](const std::string& cname, const xmlNodePtr element) {
            if (cname == "spo_scanner")
                if (myComm->rank() == 0) {
                    SPOSetScannerT<T> ascanner(sposets, targetPtcl, ptclPool);
                    ascanner.put(element);
                }
        });
}

template <typename T>
void
SPOSetBuilderFactoryT<T>::addSPOSet(std::unique_ptr<SPOSetT<T>> spo)
{
    if (spo->getName().empty())
        myComm->barrier_and_abort(
            "sposet created in sposet_collection must have a name!");

    if (sposets.find(spo->getName()) != sposets.end())
        myComm->barrier_and_abort("The name of each sposet must be unique! '" +
            spo->getName() + "' exists.");
    else
        sposets.emplace(spo->getName(), std::move(spo));
}

template <typename T>
std::string SPOSetBuilderFactoryT<T>::basisset_tag = "basisset";

template class SPOSetBuilderFactoryT<std::complex<double>>;
template class SPOSetBuilderFactoryT<std::complex<float>>;
template class SPOSetBuilderFactoryT<double>;
template class SPOSetBuilderFactoryT<float>;
} // namespace qmcplusplus
