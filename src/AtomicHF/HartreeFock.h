//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_ATOMICHF_HARTREEFOCK_H
#define OHMMS_ATOMICHF_HARTREEFOCK_H
#include "AtomicHF/HFConfiguration.h"
#include "AtomicHF/HFAtomicOrbitals.h"
#include "AtomicHF/RadialPotentialSet.h"
#include "Numerics/RnlRelTransform.h"
#include "Numerics/RnlTransform.h"
#include "Numerics/HarmTransform.h"
#include "Numerics/Numerov.h"
#include "Numerics/RadialFunctorUtility.h"
#include "OhmmsData/libxmldefs.h"
#include <libxml++/libxml++.h>
#include <ofstream>

namespace ohmmshf
{

struct HartreeFock
{

  typedef HFAtomicOrbitals::value_type value_type;
  typedef HFAtomicOrbitals::RadialOrbital_t RadialOrbital_t;

  int maxiter;
  value_type eig_tol, scf_tol, ratio;
  RadialPotentialSet& Pot;
  HFAtomicOrbitals& Psi;

  HartreeFock(RadialPotentialSet& pot, HFAtomicOrbitals& psi,
              const xmlpp::Node* root):
    Pot(pot), Psi(psi), maxiter(1000), eig_tol(1e-12),
    scf_tol(1e-8), ratio(0.35)
  {
    put(root);
  }


  bool put(const xmlpp::Node* q)
  {
    using namespace xmlpp;
    NodeSet eset = q->find("./EigenSolve");
    if(eset.empty())
    {
      WARNMSG("Using default values for eigensolver.")
      XMLReport("maximum iterations = " << maxiter)
      XMLReport("eigentolerance = " << eig_tol)
      XMLReport("scftolerance = " << scf_tol)
      XMLReport("ratio = " << ratio)
    }
    else
    {
      NodeSet pset = eset[0]->find("./Parameter");
      for(int i=0; i<pset.size(); i++)
      {
        Element* cur = dynamic_cast<Element*>(pset[i]);
        Attribute* att = cur->get_attribute("name");
        if(att)
        {
          const std::string& aname = att->get_value();
          xmlNode* ccur = cur->cobj();
          if(aname == "max_iter")
          {
            putContent(maxiter,ccur);
            XMLReport("maximum iterations = " << maxiter)
          }
          else
            if(aname == "eig_tol")
            {
              putContent(eig_tol,ccur);
              XMLReport("eigentolerance = " << eig_tol)
            }
            else
              if(aname == "en_tol")
              {
                putContent(scf_tol,ccur);
                XMLReport("scftolerance = " << scf_tol)
              }
              else
                if(aname == "mix_ratio")
                {
                  putContent(ratio,ccur);
                  XMLReport("ratio = " << ratio)
                }
        }
      }
    }
  }


  template<class Transform_t>
  void run_hf(Transform_t* fake, int norb)
  {
    typedef Numerov<Transform_t, RadialOrbital_t> Numerov_t;
    value_type Vtotal,KEnew, KEold,E;
    std::vector<value_type> energy(Pot.size()), epsilon(norb);
    value_type lowerbound, upperbound;
    int iter = 0;
    Vtotal = Pot.evaluate(Psi,energy,norb);
    Pot.mix(0.0);
    KEnew = Pot.calcKE(Psi,0,norb);
    std::string label("spdf");
    std::ofstream log_stream("atomicHF.log");
    log_stream.precision(8);
    do
    {
      KEold = KEnew;
      value_type eigsum = 0.0;
      ///loop over the orbitals
      for(int ob=0; ob < norb; ob++)
      {
        ///set up the transformer
        Transform_t es(Pot.V[ob], Psi.N[ob], Psi.L[ob],
                       Psi.CuspParam, Psi.EigenBoundParam);
        ///initialize the numerov solver
        Numerov_t numerov(es,Psi(ob));
        ///calculate the lower and upper bounds for the eigenvalue
        es.EnergyBound(lowerbound,upperbound);
        ///perform the numerov algorithm
        ///calculate the eigenvalue and the corresponding orbital
        eigsum += (epsilon[ob] =
                     numerov.solve(lowerbound, upperbound, eig_tol));
        log_stream << Psi.N[ob]<< label[Psi.L[ob]] << '\t' << epsilon[ob] << std::endl;
      }
      log_stream << std::endl;
      ///normalize the orbitals
      Psi.normalize(norb);
      ///restrict the orbitals
      Psi.applyRestriction(norb);
      ///calculate the new kinetic energy
      KEnew = Pot.calcKE(Psi,eigsum,norb);
      ///the total energy
      E = KEnew + Vtotal;
      ///for the new orbitals Psi, calculate the new SCF potentials
      Vtotal = Pot.evaluate(Psi,energy,norb);
      Pot.applyRestriction(Psi);
      Pot.mix(ratio);
      log_stream.precision(10);
      log_stream << "Iteration #" << iter+1 << std::endl;
      log_stream << "KE    = " << std::setw(15) << KEnew
                 << "  PE     = " << std::setw(15) << Vtotal << std::endl;
      log_stream << "PE/KE = " << std::setw(15) << Vtotal/KEnew
                 << "  Energy = " << std::setw(15) << E << std::endl;
      log_stream << std::endl;
      iter++;
      ///continue the loop until the kinetic energy converges
    }
    while(std::abs(KEnew-KEold)>scf_tol && iter<maxiter);
  }

  inline void solve( std::string pottype, std::string gridtype, int norb)
  {
    if(pottype == "Harmonic")
    {
      if(gridtype == "linear")
      {
        HarmLinearTransform<RadialOrbital_t> *afake=NULL;
        run_hf(afake,norb);
      }
      else
        if(gridtype == "log")
        {
          HarmLogTransform<RadialOrbital_t> *afake=NULL;
          run_hf(afake,norb);
        }
    }
    else
      if(pottype == "Nuclear_Scalar_Rel")
      {
        if(gridtype == "linear")
        {
          RnlRelLinearTransform<RadialOrbital_t> *afake=NULL;
          run_hf(afake,norb);
        }
        else
          if(gridtype == "log")
          {
            RnlRelLogTransform<RadialOrbital_t> *afake=NULL;
            run_hf(afake,norb);
          }
      }
      else
        if(pottype == "Nuclear")
        {
          if(gridtype == "linear")
          {
            RnlLinearTransform<RadialOrbital_t> *afake=NULL;
            run_hf(afake,norb);
          }
          else
            if(gridtype == "log")
            {
              RnlLogTransform<RadialOrbital_t> *afake=NULL;
              run_hf(afake,norb);
            }
        }
        else
          if(pottype == "Pseudo")
          {
            if(gridtype == "linear")
            {
              RnlLinearTransform<RadialOrbital_t> *afake=NULL;
              run_hf(afake,norb);
            }
            else
              if(gridtype == "log")
              {
                RnlLogTransform<RadialOrbital_t> *afake=NULL;
                run_hf(afake,norb);
              }
          }
          else
            return;
  }
};
}

#endif
