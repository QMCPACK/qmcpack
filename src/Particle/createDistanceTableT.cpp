//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at
// Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of
//                    Illinois at Urbana-Champaign Jeongnim Kim,
//                    jeongnim.kim@gmail.com, University of Illinois at
//                    Urbana-Champaign Jaron T. Krogel, krogeljt@ornl.gov, Oak
//                    Ridge National Laboratory Mark A. Berrill,
//                    berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/createDistanceTableT.h"

#include "CPU/SIMD/algorithm.hpp"
#include "Particle/DistanceTableT.h"
#include "Particle/SoaDistanceTableAAT.h"
#include "Particle/SoaDistanceTableAATOMPTarget.h"
#include "Particle/SoaDistanceTableABT.h"
#include "Particle/SoaDistanceTableABTOMPTarget.h"

namespace qmcplusplus
{
/** Adding SymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return index of the distance table with the name
 */
template <typename T>
std::unique_ptr<DistanceTableT<T>>
createDistanceTableAAT(ParticleSetT<T>& s, std::ostream& description)
{
    using RealType = typename ParticleSetT<T>::RealType;
    enum
    {
        DIM = OHMMS_DIM
    };
    const int sc = s.getLattice().SuperCellEnum;
    std::unique_ptr<DistanceTableT<T>> dt;
    std::ostringstream o;
    o << "  Distance table for similar particles (A-A):" << std::endl;
    o << "    source/target: " << s.getName() << std::endl;
    o << "    Using structure-of-arrays (SoA) data layout" << std::endl;

    if (sc == SUPERCELL_BULK) {
        if (s.getLattice().DiagonalOnly) {
            o << "    Distance computations use orthorhombic periodic cell in "
                 "3D."
              << std::endl;
            dt = std::make_unique<
                SoaDistanceTableAAT<T, DIM, PPPO + SOA_OFFSET>>(s);
        }
        else {
            if (s.getLattice().WignerSeitzRadius >
                s.getLattice().SimulationCellRadius) {
                o << "    Distance computations use general periodic cell in "
                     "3D with corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableAAT<T, DIM, PPPG + SOA_OFFSET>>(s);
            }
            else {
                o << "    Distance computations use general periodic cell in "
                     "3D without corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableAAT<T, DIM, PPPS + SOA_OFFSET>>(s);
            }
        }
    }
    else if (sc == SUPERCELL_SLAB) {
        if (s.getLattice().DiagonalOnly) {
            o << "    Distance computations use orthorhombic code for periodic "
                 "cell in 2D."
              << std::endl;
            dt = std::make_unique<
                SoaDistanceTableAAT<T, DIM, PPNO + SOA_OFFSET>>(s);
        }
        else {
            if (s.getLattice().WignerSeitzRadius >
                s.getLattice().SimulationCellRadius) {
                o << "    Distance computations use general periodic cell in "
                     "2D with corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableAAT<T, DIM, PPNG + SOA_OFFSET>>(s);
            }
            else {
                o << "    Distance computations use general periodic cell in "
                     "2D without corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableAAT<T, DIM, PPNS + SOA_OFFSET>>(s);
            }
        }
    }
    else if (sc == SUPERCELL_WIRE) {
        o << "    Distance computations use periodic cell in one dimension."
          << std::endl;
        dt = std::make_unique<
            SoaDistanceTableAAT<T, DIM, SUPERCELL_WIRE + SOA_OFFSET>>(s);
    }
    else // open boundary condition
    {
        o << "    Distance computations use open boundary conditions in 3D."
          << std::endl;
        dt = std::make_unique<
            SoaDistanceTableAAT<T, DIM, SUPERCELL_OPEN + SOA_OFFSET>>(s);
    }

    description << o.str() << std::endl;
    return dt;
}

template std::unique_ptr<DistanceTableT<double>>
createDistanceTableAAT<double>(
    ParticleSetT<double>& t, std::ostream& description);
template std::unique_ptr<DistanceTableT<float>>
createDistanceTableAAT<float>(
    ParticleSetT<float>& t, std::ostream& description);
template std::unique_ptr<DistanceTableT<std::complex<double>>>
createDistanceTableAAT<std::complex<double>>(
    ParticleSetT<std::complex<double>>& t, std::ostream& description);
template std::unique_ptr<DistanceTableT<std::complex<float>>>
createDistanceTableAAT<std::complex<float>>(
    ParticleSetT<std::complex<float>>& t, std::ostream& description);

/** Adding AsymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return index of the distance table with the name
 */
template <typename T>
std::unique_ptr<DistanceTableT<T>>
createDistanceTableABT(
    const ParticleSetT<T>& s, ParticleSetT<T>& t, std::ostream& description)
{
    using RealType = typename ParticleSetT<T>::RealType;
    enum
    {
        DIM = OHMMS_DIM
    };
    const int sc = t.getLattice().SuperCellEnum;
    std::unique_ptr<DistanceTableT<T>> dt;
    std::ostringstream o;
    o << "  Distance table for dissimilar particles (A-B):" << std::endl;
    o << "    source: " << s.getName() << "  target: " << t.getName()
      << std::endl;
    o << "    Using structure-of-arrays (SoA) data layout" << std::endl;

    if (sc == SUPERCELL_BULK) {
        if (s.getLattice().DiagonalOnly) {
            o << "    Distance computations use orthorhombic periodic cell in "
                 "3D."
              << std::endl;
            dt = std::make_unique<
                SoaDistanceTableABT<T, DIM, PPPO + SOA_OFFSET>>(s, t);
        }
        else {
            if (s.getLattice().WignerSeitzRadius >
                s.getLattice().SimulationCellRadius) {
                o << "    Distance computations use general periodic cell in "
                     "3D with corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableABT<T, DIM, PPPG + SOA_OFFSET>>(s, t);
            }
            else {
                o << "    Distance computations use general periodic cell in "
                     "3D without corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableABT<T, DIM, PPPS + SOA_OFFSET>>(s, t);
            }
        }
    }
    else if (sc == SUPERCELL_SLAB) {
        if (s.getLattice().DiagonalOnly) {
            o << "    Distance computations use orthorhombic code for periodic "
                 "cell in 2D."
              << std::endl;
            dt = std::make_unique<
                SoaDistanceTableABT<T, DIM, PPNO + SOA_OFFSET>>(s, t);
        }
        else {
            if (s.getLattice().WignerSeitzRadius >
                s.getLattice().SimulationCellRadius) {
                o << "    Distance computations use general periodic cell in "
                     "2D with corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableABT<T, DIM, PPNG + SOA_OFFSET>>(s, t);
            }
            else {
                o << "    Distance computations use general periodic cell in "
                     "2D without corner image checks."
                  << std::endl;
                dt = std::make_unique<
                    SoaDistanceTableABT<T, DIM, PPNS + SOA_OFFSET>>(s, t);
            }
        }
    }
    else if (sc == SUPERCELL_WIRE) {
        o << "    Distance computations use periodic cell in one dimension."
          << std::endl;
        dt = std::make_unique<
            SoaDistanceTableABT<T, DIM, SUPERCELL_WIRE + SOA_OFFSET>>(s, t);
    }
    else // open boundary condition
    {
        o << "    Distance computations use open boundary conditions in 3D."
          << std::endl;
        dt = std::make_unique<
            SoaDistanceTableABT<T, DIM, SUPERCELL_OPEN + SOA_OFFSET>>(s, t);
    }

    description << o.str() << std::endl;
    return dt;
}

template std::unique_ptr<DistanceTableT<double>>
createDistanceTableABT<double>(const ParticleSetT<double>& s,
    ParticleSetT<double>& t, std::ostream& description);
template std::unique_ptr<DistanceTableT<float>>
createDistanceTableABT<float>(const ParticleSetT<float>& s,
    ParticleSetT<float>& t, std::ostream& description);
template std::unique_ptr<DistanceTableT<std::complex<double>>>
createDistanceTableABT<std::complex<double>>(
    const ParticleSetT<std::complex<double>>& s,
    ParticleSetT<std::complex<double>>& t, std::ostream& description);
template std::unique_ptr<DistanceTableT<std::complex<float>>>
createDistanceTableABT<std::complex<float>>(
    const ParticleSetT<std::complex<float>>& s,
    ParticleSetT<std::complex<float>>& t, std::ostream& description);
} // namespace qmcplusplus
