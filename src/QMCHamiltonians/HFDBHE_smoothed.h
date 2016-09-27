//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D.C. Yang, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HFDBHE_SMOOTHED_H
#define QMCPLUSPLUS_HFDBHE_SMOOTHED_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#if defined(HAVE_MKL)
#include <mkl_vml_functions.h>
#endif

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief HFDBHE_smoothed for the indentical source and target particle sets.
 */
struct HFDBHE_smoothed_phy: public QMCHamiltonianBase
{
  Return_t rc, mirrorshift, smooth;
  // remember that the default units are Kelvins and Bohrs
  DistanceTableData* d_table;
  ParticleSet* PtclRef;

  HFDBHE_smoothed_phy(ParticleSet& P);

  ~HFDBHE_smoothed_phy() { }

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t dampF(Return_t r)
  {
    const Return_t D = 8.301460704;
    if (r < D)
      return std::exp(-(D/r - 1.0)*(D/r - 1.0));
    else
      return 1.0;
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HFD-B(He)_smoothed_phy: " << PtclRef->getName();
    return true;
  }

  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

  void addCorrection(QMCHamiltonian& targetH);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HFDBHE_smoothed_phy(qp);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex = plist.add(myName);
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex] = Value;
  }
};



/** @ingroup hamiltonian
 *@brief HFDBHE_smoothed_np for the indentical source and target particle sets.
 */
struct HFDBHE_smoothed_aux: public QMCHamiltonianBase
{
  const HFDBHE_smoothed_phy* phyH;
  Return_t tailcorr;

  // epsilon = 3.42016039e-5, rm = 5.607384357
  // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
  HFDBHE_smoothed_aux(const HFDBHE_smoothed_phy* orig): phyH(orig)
  {
    Return_t N0 = phyH->PtclRef->getTotalNum(), rho = N0/phyH->PtclRef->Lattice.Volume, rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    //      tailcorr = -0.000274387876825;	// N = 32, RHO = 0.022 /Angstrom^3 : 11x
    // Note the 2 powers of N
    // Specific value for:
    // Equilibrium density 0.022 Angstroms^-3
    // N = 32
    // Will be different for other densities.
    // Generalize later.
    //      tailcorr = -0.0002704899426;	// N = 32, RHO = 0.3648 sigma^2 : 01x
    //      tailcorr = -0.0002650846591;	// N = 64, RHO = 0.3648 sigma^3 : 02x
    //      tailcorr = -0.0002618259317;	// N = 128, RHO = 0.3648 sigma^2 : 03x
    //      tailcorr = -0.0002598311735;	// N = 256, RHO = 0.3648 sigma^2 : 04x
#if defined(HAVE_MKL)
    tailcorr = 2.0*M_PI*rho*N0*(-26.21440794815934*std::pow(rc,-7.0) - 2.822014763523686*std::pow(rc,-5.0) - 0.4870028588757336*std::pow(rc,-3.0)
                                + 6.394277070147666*std::exp(-std::pow(3.455073733699311 + 0.26965201272692674*rc,2.0))
                                * (-1.3471677840430893e7 + 1.0514001507549516e6*rc + 1.3142529690505246e13*std::exp((1.8633351728239138 + 0.07271220796768266*rc)*rc)
                                   * (erfc(3.455073733699311 + 0.26965201272692674*rc)) ) );
    app_log() << "  HFDBHE_smoothed_aux tail correction is " << tailcorr << std::endl;
#else
    tailcorr = 0.0;
    app_log() << "  Tail correction without Intel MKL is not supported.  Do this manually." << std::endl
              << "  tailcorr value set to 0." << std::endl;
#endif
  }

  ~HFDBHE_smoothed_aux() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    Return_t N0 = phyH->PtclRef->getTotalNum(), rho = N0/phyH->PtclRef->Lattice.Volume, rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    //      tailcorr = -0.000274387876825;	// N = 32, RHO = 0.022 /Angstroms^3 : 11x
    //      tailcorr = -0.0002704899426;	// N = 32, RHO = 0.3648 sigma^2 : 01x
    //      tailcorr = -0.0002650846591;	// N = 64, RHO = 0.3648 sigma^3 : 02x
    //      tailcorr = -0.0002618259317;	// N = 128, RHO = 0.3648 sigma^2 : 03x
    //      tailcorr = -0.0002598311735;	// N = 256, RHO = 0.3648 sigma^2 : 04x
#if defined(HAVE_MKL)
    tailcorr = 2.0*M_PI*rho*N0*(-26.21440794815934*std::pow(rc,-7.0) - 2.822014763523686*std::pow(rc,-5.0) - 0.4870028588757336*std::pow(rc,-3.0)
                                + 6.394277070147666*std::exp(-std::pow(3.455073733699311 + 0.26965201272692674*rc,2.0))
                                * (-1.3471677840430893e7 + 1.0514001507549516e6*rc + 1.3142529690505246e13*std::exp((1.8633351728239138 + 0.07271220796768266*rc)*rc)
                                   * (erfc(3.455073733699311 + 0.26965201272692674*rc)) ) );
    app_log() << "  HFDBHE_smoothed_aux tail correction is " << tailcorr << std::endl;
#else
    tailcorr = 0.0;
    app_log() << "  Tail correction without Intel MKL is not supported.  Do this manually." << std::endl
              << "  tailcorr value set to 0." << std::endl;
#endif
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value = -(phyH->smooth) + tailcorr;
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

#if defined(HAVE_MKL)
  double erfc(const double x)
  {
    const double *a;
    a = &x;
    double y[1];
    vdErfc(1, a, y);
    return y[0];
  }
#endif

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HFD-B(He)_smoothed_aux: " << phyH->PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return false;
  }
};
}
#endif
