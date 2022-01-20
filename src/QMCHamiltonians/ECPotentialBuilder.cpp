//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ECPotentialBuilder.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "QMCHamiltonians/L2Potential.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/OneDimNumGridFunctor.h"
#ifdef QMC_CUDA
#include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
#include "QMCHamiltonians/LocalECPotential_CUDA.h"
#include "QMCHamiltonians/NonLocalECPotential_CUDA.h"
#endif

namespace qmcplusplus
{
/** constructor
 *\param ions the positions of the ions
 *\param els the positions of the electrons
 *\param psi trial wavefunction
 */
ECPotentialBuilder::ECPotentialBuilder(QMCHamiltonian& h,
                                       ParticleSet& ions,
                                       ParticleSet& els,
                                       TrialWaveFunction& psi,
                                       Communicate* c)
    : MPIObjectBase(c),
      hasLocalPot(false),
      hasNonLocalPot(false),
      hasSOPot(false),
      hasL2Pot(false),
      targetH(h),
      IonConfig(ions),
      targetPtcl(els),
      targetPsi(psi)
{}

bool ECPotentialBuilder::put(xmlNodePtr cur)
{
  if (localPot.empty())
  {
    int ng = IonConfig.getSpeciesSet().getTotalNum();
    localZeff.resize(ng, 1);
    localPot.resize(ng);
    nonLocalPot.resize(ng);
    soPot.resize(ng);
    L2Pot.resize(ng);
  }
  std::string ecpFormat;
  std::string NLPP_algo;
  std::string use_DLA;
  std::string pbc;
  std::string forces;
  std::string physicalSO;

  OhmmsAttributeSet pAttrib;
  pAttrib.add(ecpFormat, "format", {"table", "xml"});
  pAttrib.add(NLPP_algo, "algorithm", {"", "batched", "non-batched"});
  pAttrib.add(use_DLA, "DLA", {"no", "yes"});
  pAttrib.add(pbc, "pbc", {"yes", "no"});
  pAttrib.add(forces, "forces", {"no", "yes"});
  pAttrib.add(physicalSO, "physicalSO", {"yes", "no"});
  pAttrib.put(cur);

  if (NLPP_algo.empty())
#ifdef ENABLE_OFFLOAD
    NLPP_algo = "batched";
#else
    NLPP_algo = "non-batched";
#endif

  bool doForces = (forces == "yes") || (forces == "true");
  if (use_DLA == "yes")
    app_log() << "    Using determinant localization approximation (DLA)" << std::endl;
  if (ecpFormat == "xml")
  {
    useXmlFormat(cur);
  }
  else
  {
    useSimpleTableFormat();
  }

  ///create LocalECPotential
  bool usePBC = !(IonConfig.getLattice().SuperCellEnum == SUPERCELL_OPEN || pbc == "no");


  if (hasLocalPot)
  {
    if (IonConfig.getLattice().SuperCellEnum == SUPERCELL_OPEN || pbc == "no")
    {
#ifdef QMC_CUDA
      std::unique_ptr<LocalECPotential_CUDA> apot = std::make_unique<LocalECPotential_CUDA>(IonConfig, targetPtcl);
#else
      std::unique_ptr<LocalECPotential> apot = std::make_unique<LocalECPotential>(IonConfig, targetPtcl);
#endif
      for (int i = 0; i < localPot.size(); i++)
        if (localPot[i])
          apot->add(i, std::move(localPot[i]), localZeff[i]);
      targetH.addOperator(std::move(apot), "LocalECP");
    }
    else
    {
      if (doForces)
        app_log() << "  Will compute forces in CoulombPBCAB.\n" << std::endl;
#ifdef QMC_CUDA
      std::unique_ptr<CoulombPBCAB_CUDA> apot = std::make_unique<CoulombPBCAB_CUDA>(IonConfig, targetPtcl, doForces);
#else
      std::unique_ptr<CoulombPBCAB> apot     = std::make_unique<CoulombPBCAB>(IonConfig, targetPtcl, doForces);
#endif
      for (int i = 0; i < localPot.size(); i++)
      {
        if (localPot[i])
          apot->add(i, std::move(localPot[i]));
      }
      targetH.addOperator(std::move(apot), "LocalECP");
    }
  }
  if (hasNonLocalPot)
  {
#ifdef QMC_CUDA
    std::unique_ptr<NonLocalECPotential_CUDA> apot =
        std::make_unique<NonLocalECPotential_CUDA>(IonConfig, targetPtcl, targetPsi, usePBC, doForces,
                                                   use_DLA == "yes");
#else
    std::unique_ptr<NonLocalECPotential> apot =
        std::make_unique<NonLocalECPotential>(IonConfig, targetPtcl, targetPsi, doForces, use_DLA == "yes");
#endif
    int nknot_max = 0;
    for (int i = 0; i < nonLocalPot.size(); i++)
    {
      if (nonLocalPot[i])
      {
        nknot_max = std::max(nknot_max, nonLocalPot[i]->getNknot());
        if (NLPP_algo == "batched")
          nonLocalPot[i]->initVirtualParticle(targetPtcl);
        apot->addComponent(i, std::move(nonLocalPot[i]));
      }
    }
    app_log() << "\n  Using NonLocalECP potential \n"
              << "    Maximum grid on a sphere for NonLocalECPotential: " << nknot_max << std::endl;
    if (NLPP_algo == "batched")
      app_log() << "    Using batched ratio computing in NonLocalECP" << std::endl;

    targetH.addOperator(std::move(apot), "NonLocalECP");
  }
  if (hasSOPot)
  {
#ifndef QMC_COMPLEX
    APP_ABORT("SOECPotential evaluations require complex build. Rebuild with -D QMC_COMPLEX=1\n");
#endif
    if (physicalSO == "yes")
      app_log() << "    Spin-Orbit potential included in local energy" << std::endl;
    else if (physicalSO == "no")
      app_log() << "    Spin-Orbit potential is not included in local energy" << std::endl;
    else
      APP_ABORT("physicalSO must be set to yes/no. Unknown option given\n");

    std::unique_ptr<SOECPotential> apot = std::make_unique<SOECPotential>(IonConfig, targetPtcl, targetPsi);
    int nknot_max                       = 0;
    int sknot_max                       = 0;
    for (int i = 0; i < soPot.size(); i++)
    {
      if (soPot[i])
      {
        nknot_max = std::max(nknot_max, soPot[i]->getNknot());
        sknot_max = std::max(sknot_max, soPot[i]->getSknot());
        apot->addComponent(i, std::move(soPot[i]));
      }
    }
    app_log() << "\n  Using SOECP potential \n"
              << "    Maximum grid on a sphere for SOECPotential: " << nknot_max << std::endl;
    app_log() << "    Maximum grid for Simpson's rule for spin integral: " << sknot_max << std::endl;

    if (physicalSO == "yes")
      targetH.addOperator(std::move(apot), "SOECP"); //default is physical operator
    else
      targetH.addOperator(std::move(apot), "SOECP", false);
  }
  if (hasL2Pot)
  {
    std::unique_ptr<L2Potential> apot = std::make_unique<L2Potential>(IonConfig, targetPtcl, targetPsi);
    for (int i = 0; i < L2Pot.size(); i++)
      if (L2Pot[i])
        apot->add(i, std::move(L2Pot[i]));
    app_log() << "\n  Using L2 potential" << std::endl;
    targetH.addOperator(std::move(apot), "L2");
  }

  app_log().flush();
  return true;
}

void ECPotentialBuilder::useXmlFormat(xmlNodePtr cur)
{
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "pseudo")
    {
      std::string href("none");
      std::string ionName("none");
      std::string format("xml");
      int nrule = -1;
      //RealType rc(2.0);//use 2 Bohr
      OhmmsAttributeSet hAttrib;
      hAttrib.add(href, "href");
      hAttrib.add(ionName, "elementType");
      hAttrib.add(ionName, "symbol");
      hAttrib.add(format, "format");
      hAttrib.add(nrule, "nrule");
      //hAttrib.add(rc,"cutoff");
      hAttrib.put(cur);
      SpeciesSet& ion_species(IonConfig.getSpeciesSet());
      int speciesIndex      = ion_species.findSpecies(ionName);
      int chargeIndex       = ion_species.findAttribute("charge");
      int AtomicNumberIndex = ion_species.findAttribute("atomicnumber");
      if (AtomicNumberIndex == -1)
        AtomicNumberIndex = ion_species.findAttribute("atomic_number");
      bool success = false;
      if (speciesIndex < ion_species.getTotalNum())
      {
        app_log() << std::endl << "  Adding pseudopotential for " << ionName << std::endl;

        ECPComponentBuilder ecp(ionName, myComm, nrule);
        if (format == "xml")
        {
          if (href == "none")
          {
            success = ecp.put(cur);
          }
          else
          {
            success = ecp.parse(href, cur);
          }
        }
        else if (format == "casino")
        {
          //success=ecp.parseCasino(href,rc);
          success = ecp.parseCasino(href, cur);
        }
        if (success)
        {
          if (OHMMS::Controller->rank() == 0)
            ecp.printECPTable();
          if (ecp.pp_loc)
          {
            localPot[speciesIndex]  = std::move(ecp.pp_loc);
            localZeff[speciesIndex] = ecp.Zeff;
            hasLocalPot             = true;
          }
          if (ecp.pp_nonloc)
          {
            hasNonLocalPot            = true;
            nonLocalPot[speciesIndex] = std::move(ecp.pp_nonloc);
          }
          if (ecp.pp_so)
          {
            hasSOPot            = true;
            soPot[speciesIndex] = std::move(ecp.pp_so);
          }
          if (ecp.pp_L2)
          {
            hasL2Pot            = true;
            L2Pot[speciesIndex] = std::move(ecp.pp_L2);
          }
          if (chargeIndex == -1)
          {
            app_error() << "  Ion species " << ionName << " needs parameter \'charge\'" << std::endl;
            success = false;
          }
          else
          {
            RealType ion_charge = ion_species(chargeIndex, speciesIndex);
            if (std::fabs(ion_charge - ecp.Zeff) > 1e-4)
            {
              app_error() << "  Ion species " << ionName << " charge " << ion_charge << " pseudopotential charge "
                          << ecp.Zeff << " mismatch!" << std::endl;
              success = false;
            }
          }
          if (AtomicNumberIndex == -1)
          {
            app_error() << "  Ion species " << ionName << " needs parameter \'atomicnumber\'" << std::endl;
            success = false;
          }
          else
          {
            int atomic_number = ion_species(AtomicNumberIndex, speciesIndex);
            if (atomic_number != ecp.AtomicNumber)
            {
              app_error() << "  Ion species " << ionName << " atomicnumber " << atomic_number
                          << " pseudopotential atomic-number " << ecp.AtomicNumber << " mismatch!" << std::endl;
              success = false;
            }
          }
        }
      }
      else
      {
        app_error() << "  Ion species " << ionName << " is not found." << std::endl;
      }
      if (!success)
      {
        app_error() << "  Failed to add pseudopotential for element " << ionName << std::endl;
        myComm->barrier_and_abort("ECPotentialBuilder::useXmlFormat failed!");
      }
    }
    cur = cur->next;
  }
}

/** reimplement simple table format used by NonLocalPPotential
 */
void ECPotentialBuilder::useSimpleTableFormat()
{
  SpeciesSet& Species(IonConfig.getSpeciesSet());
  int ng(Species.getTotalNum());
  int icharge(Species.addAttribute("charge"));
  for (int ig = 0; ig < ng; ig++)
  {
    std::vector<RealType> grid_temp, pp_temp;
    std::string species(Species.speciesName[ig]);
    std::string fname = species + ".psf";
    std::ifstream fin(fname.c_str(), std::ios_base::in);
    if (!fin)
    {
      ERRORMSG("Could not open file " << fname)
      exit(-1);
    }
    // Read Number of potentials (local and non) for this atom
    int npotentials;
    fin >> npotentials;
    int lmax      = -1;
    int numnonloc = 0;
    RealType rmax(0.0);
    app_log() << "  ECPotential for " << species << std::endl;
    std::unique_ptr<NonLocalECPComponent> mynnloc;
    typedef OneDimCubicSpline<RealType> CubicSplineFuncType;
    for (int ij = 0; ij < npotentials; ij++)
    {
      int angmom, npoints;
      fin >> angmom >> npoints;
      OneDimNumGridFunctor<RealType> inFunc;
      inFunc.put(npoints, fin);
      if (angmom < 0)
      //local potential, input is rescale by -r/z
      {
        RealType zinv = -1.0 / Species(icharge, ig);
        int ng        = npoints - 1;
        RealType rf   = 5.0;
        ng            = static_cast<int>(rf * 100) + 1; //use 1e-2 resolution
        auto agrid    = std::make_unique<LinearGrid<RealType>>();
        agrid->set(0, rf, ng);
        std::vector<RealType> pp_temp(ng);
        pp_temp[0] = 0.0;
        for (int j = 1; j < ng; j++)
        {
          RealType r((*agrid)[j]);
          pp_temp[j] = r * zinv * inFunc.splint(r);
        }
        pp_temp[ng - 1] = 1.0;
        auto app        = std::make_unique<RadialPotentialType>(std::move(agrid), pp_temp);
        app->spline();
        localPot[ig] = std::move(app);
        app_log() << "    LocalECP l=" << angmom << std::endl;
        app_log() << "      Linear grid=[0," << rf << "] npts=" << ng << std::endl;
        hasLocalPot = true; //will create LocalECPotential
      }
      else
      {
        hasNonLocalPot = true; //will create NonLocalECPotential
        if (!mynnloc)
          mynnloc = std::make_unique<NonLocalECPComponent>();
        RealType rf = inFunc.rmax();
        auto agrid  = std::make_unique<LinearGrid<RealType>>();
        int ng      = static_cast<int>(rf * 100) + 1;
        agrid->set(0.0, rf, ng);
        app_log() << "    NonLocalECP l=" << angmom << " rmax = " << rf << std::endl;
        app_log() << "      Linear grid=[0," << rf << "] npts=" << ng << std::endl;
        std::vector<RealType> pp_temp(ng);
        //get the value
        pp_temp[0] = inFunc(0);
        for (int j = 1; j < ng; j++)
        {
          pp_temp[j] = inFunc.splint((*agrid)[j]);
        }
        auto app = new RadialPotentialType(std::move(agrid), pp_temp);
        app->spline();
        mynnloc->add(angmom, app);
        lmax = std::max(lmax, angmom);
        rmax = std::max(rmax, rf);
        numnonloc++;
      }
      if (mynnloc)
      {
        mynnloc->setLmax(lmax);
        mynnloc->setRmax(rmax);
        app_log() << "    Maximum cutoff of NonLocalECP " << rmax << std::endl;
      }
    }
    fin.close();
    if (mynnloc)
    {
      int numsgridpts   = 0;
      std::string fname = species + ".sgr";
      std::ifstream fin(fname.c_str(), std::ios_base::in);
      if (!fin)
      {
        app_error() << "Could not open file " << fname << std::endl;
        exit(-1);
      }
      PosType xyz;
      RealType weight;
      while (fin >> xyz >> weight)
      {
        mynnloc->addknot(xyz, weight);
        numsgridpts++;
      }
      //cout << "Spherical grid : " << numsgridpts << " points" << std::endl;
      mynnloc->resize_warrays(numsgridpts, numnonloc, lmax);
      nonLocalPot[ig] = std::move(mynnloc);
    }
  } //species
}
} // namespace qmcplusplus
