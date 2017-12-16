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
    
    


#include "AtomicHF/HFConfiguration.h"
#include "Numerics/Clebsch_Gordan.h"
#include "AtomicHF/fillshells.h"
#include "AtomicHF/SphericalPotential/ZOverRPotential.h"
#include "AtomicHF/SphericalPotential/HarmonicPotential.h"
#include "AtomicHF/SphericalPotential/StepPotential.h"
#include "AtomicHF/SphericalPotential/RegularLinearTransform.h"
#include "AtomicHF/SphericalPotential/RegularLogTransform.h"
#include "AtomicHF/SphericalPotential/NuclearLinearTransform.h"
#include "AtomicHF/SphericalPotential/NuclearLogTransform.h"
#include "Numerics/Numerov.h"
#include "Numerics/RadialFunctorUtility.h"
#include "SQD/HartreeFock.h"

namespace ohmmshf
{

HartreeFock::HartreeFock(RadialPotentialSet& pot,
                         SphericalOrbitalTraits::BasisSetType& psi,
                         const xmlpp::Node* root):
  Pot(pot), Psi(psi), maxiter(1000), eig_tol(1e-12),
  scf_tol(1e-8), ratio(0.35)
{
  LogFileName = "atomiHF.log";
  put(root);
}

void HartreeFock::setRoot(const std::string& aroot)
{
  LogFileName = aroot;
  LogFileName.append(".log");
}

/**
   @param fake Transformation object
   @param norb the number of eigen vectors to be obtained
   @brief perform self-consistent Hartree-Fock calculations
   *
   *The first argument is used to tell the compiler which transform object is used but
   *not meaningful for the calculations. This is to force the compiler inline
   *everything possible. A better solution will be implement using traits.
   */
template<class Transform_t>
inline void
HartreeFock::run(Transform_t* fake, int norb)
{
  typedef Numerov<Transform_t, RadialOrbital_t> Numerov_t;
  value_type Vtotal,KEnew, KEold,E;
  value_type lowerbound, upperbound;
  std::vector<value_type> energy(Pot.size());
  eigVal.resize(norb);
  int iter = 0;
  Vtotal = Pot.evaluate(Psi,energy,norb);
  Pot.mix(0.0);
  KEnew = Pot.calcKE(Psi,0,norb);
  std::string label("spdf");
  std::ofstream log_stream(LogFileName.c_str());
  log_stream.precision(8);
  do
  {
    KEold = KEnew;
    value_type eigsum = 0.0;
    //loop over the orbitals
    for(int ob=0; ob < norb; ob++)
    {
      //set the number of nodes of the eigen vector
      Pot.V[ob].setNumOfNodes(Pot.getNumOfNodes(Psi.N[ob],Psi.L[ob]));
      //calculate the lower and upper bounds for the eigenvalue
      Pot.EnergyBound(lowerbound,upperbound);
      //set up the transformer
      Transform_t es(Pot.V[ob], Psi.N[ob], Psi.L[ob],Psi.CuspParam, Pot.getMass());
      //initialize the numerov solver
      Numerov_t numerov(es,Psi(ob));
      //calculate the eigenvalue and the corresponding orbital
      eigsum += (eigVal[ob] =
                   numerov.solve(lowerbound, upperbound, eig_tol));
      log_stream << Psi.N[ob]<< label[Psi.L[ob]] << '\t' << eigVal[ob] << std::endl;
    }
    log_stream << std::endl;
    //normalize the orbitals
    Psi.normalize(norb);
    //restrict the orbitals
    Psi.applyRestriction(norb);
    //calculate the new kinetic energy
    KEnew = Pot.calcKE(Psi,eigsum,norb);
    //the total energy
    E = KEnew + Vtotal;
    //for the new orbitals Psi, calculate the new SCF potentials
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

/**
   @param pottype the type of potential
   @param gridtype the type of the radial grid
   @param norb the number of target eigen vectors
   @brief Instantiate a transformation function based on the potential and grid type and call run.
   *
   */
void HartreeFock::solve( std::string pottype, std::string gridtype, int norb)
{
  if(pottype == "Harmonic" || pottype == "Step")
  {
    if(gridtype == "linear")
    {
      RegularLinearTransform<RadialOrbital_t> *afake=NULL;
      run(afake,norb);
    }
    else
      if(gridtype == "log")
      {
        RegularLogTransform<RadialOrbital_t> *afake=NULL;
        run(afake,norb);
      }
  }
  else
    if(pottype == "Nuclear")
    {
      if(gridtype == "linear")
      {
        NuclearLinearTransform<RadialOrbital_t> *afake=NULL;
        run(afake,norb);
      }
      else
        if(gridtype == "log")
        {
          NuclearLogTransform<RadialOrbital_t> *afake=NULL;
          run(afake,norb);
        }
    }
    else
      return;
}

/**
   @param q the current xml node which contains the parameter definitions for the eigen solver
   @return true if successful
   @brief Set the parameters for the eigen solver
   *
   *Aailable parameters
   * - max_iter, the maximum self-consistent iterations, default=1000
   * - eig_tol, the eigen-value tolerance, default = \f$1 e^{-12}\f$
   * - en_tol, the tolerance of the self-consistent loops, default = \f$1 e^{-8}\f$
   * - mix_ratio, the mixing ratio of the charge density, default = 0.35
   */
bool HartreeFock::put(const xmlpp::Node* q)
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

bool parseXMLFile(RadialPotentialSet& Pot,
                  SphericalOrbitalTraits::BasisSetType& Psi,
                  std::string& name,
                  std::string& pottype,
                  std::string& gridtype,
                  const xmlpp::Node* root)
{
  using namespace xmlpp;
  int nshells = 0;
  OneDimGridBase<double>* grid;
  NodeSet aset = root->find("//Atom");
  if(aset.empty())
  {
    ERRORMSG("Must include atomic information.")
    return false;
  }
  else
  {
    Element* cur = dynamic_cast<Element*>(aset[0]);
    Attribute* att1 = cur->get_attribute("name");
    Attribute* att2 = cur->get_attribute("num_closed_shells");
    if(att1)
    {
      name = att1->get_value();
      XMLReport("Name = " << name)
    }
    else
    {
      ERRORMSG("Must specify name.")
      return false;
    }
    if(att2)
    {
      nshells = atoi(att2->get_value().c_str());
      XMLReport("Number of Closed Shells = " << nshells)
    }
    else
    {
      WARNMSG("Number of Closed Shells = " << nshells)
    }
    NodeSet gset = aset[0]->find("./Grid");
    if(gset.empty())
    {
      ERRORMSG("No grid information")
      return false;
    }
    else
    {
      Element* p = dynamic_cast<Element*>(gset[0]);
      Element::AttributeList atts = p->get_attributes();
      Element::AttributeList::iterator it = atts.begin();
      double scale = 1.0;
      double min = 0.001;
      double max = 1000.0;
      int npts = 2001;
      while(it != atts.end())
      {
        const std::string& aname = (*it)->get_name();
        if(aname == "type")
        {
          gridtype = (*it)->get_value();
        }
        else
          if(aname == "scale")
          {
            scale = atof((*it)->get_value().c_str());
          }
        it++;
      }
      XMLReport("Grid type = " << gridtype)
      NodeSet gpset = gset[0]->find("./Parameter");
      if(gpset.empty())
      {
        WARNMSG("Using default grid values")
      }
      else
      {
        for(int i=0; i<gpset.size(); i++)
        {
          Element* cur = dynamic_cast<Element*>(gpset[i]);
          Attribute* att = cur->get_attribute("name");
          if(att)
          {
            const std::string& aname = att->get_value();
            xmlNode* curc = cur->cobj();
            if(aname == "min")
              putContent(min,curc);
            else
              if(aname == "max")
                putContent(max,curc);
              else
                if(aname == "npts")
                  putContent(npts,curc);
          }
        }
      }
      if(gridtype == "log")
      {
        grid = new LogGrid<double>;
        min/=scale;
        max/=sqrt(scale);
        grid->set(min,max,npts);
        XMLReport("rmin = " << min << ", rmax = " << max
                  //	<< ", dh = " << grid->dh() << ", npts = " << npts)
                  << ", npts = " << npts)
      }
      else
        if(gridtype == "linear")
        {
          grid = new LinearGrid<double>;
          min/=scale;
          max/=sqrt(scale);
          grid->set(min,max,npts);
          XMLReport("rmin = " << min << ", rmax = " << max
                    //	   << ", dh = " << grid->dh() << ", npts = " << npts)
                    << ", npts = " << npts)
        }
        else
        {
          ERRORMSG("Grid Type Options: Log or Linear.")
          return false;
        }
    }
    // pass grid to wavefunction
    Psi.m_grid = grid;
    // class to generate Clebsh Gordan coeff.
    Clebsch_Gordan* CG_coeff = NULL;
    NodeSet potset = aset[0]->find("./Potential");
    if(potset.empty())
    {
      ERRORMSG("Must provide potential information.")
      return false;
    }
    else
    {
      Element* cur = dynamic_cast<Element*>(potset[0]);
      xmlNodePtr cur_sub = cur->cobj()->xmlChildrenNode;
      while(cur_sub != NULL)
      {
        std::string pname((const char*)(cur_sub->name));
        if(pname == "Parameter")
        {
          xmlAttrPtr att=cur_sub->properties;
          while(att != NULL)
          {
            std::string vname((const char*)(att->children->content));
            if(vname == "mass")
            {
              double m=1.0;
              putContent(m,cur_sub);
              Pot.setMass(m);
              XMLReport("The effective mass is set to " << Pot.getMass())
            }
            att = att->next;
          }
        }
        cur_sub = cur_sub->next;
      }
      Attribute* att = cur->get_attribute("type");
      if(att)
      {
        pottype = att->get_value();
        NodeSet oset = aset[0]->find("./OrbitalSet");
        if(oset.empty())
        {
          ERRORMSG("Must specify OrbitalSet")
          return false;
        }
        else
        {
          Element* cur = dynamic_cast<Element*>(oset[0]);
          Attribute* att = cur->get_attribute("rest_type");
          if(att)
            Psi.Restriction = att->get_value();
          XMLReport("Orbital restriction type = " << Psi.Restriction)
          xmlNode* curc = cur->cobj();
          if(pottype == "Harmonic" || pottype == "Step")
          {
            if(nshells)
              FillShellsHarmPot(Psi,nshells);
          }
          else
          {
            if(nshells)
              FillShellsNucPot(Psi,nshells);
          }
          Psi.put(curc);
          Pot.initialize(Psi);
          LOGMSG("No of Orbs = " << Psi.size());
          LOGMSG("Orbital | Orbital ID");
          for(int j=0; j < Psi.size(); j++)
          {
            LOGMSG(j << " " << Psi.ID[j]);
          }
          LOGMSG("ID | IDcount");
          for(int j=0; j < Psi.NumUniqueOrb; j++)
          {
            LOGMSG(j << " " << Psi.IDcount[j]);
          }
        }
        XMLReport("Potential Type = " << pottype)
      }
      else
      {
        ERRORMSG("Must provide Potential Type.")
        return false;
      }
      if(pottype == "Nuclear")
      {
        XMLReport("Creating a Nuclear Potential.")
        double Z;
        NodeSet gpset = potset[0]->find("./Parameter");
        if(gpset.empty())
        {
          Z = Psi.NumOrb;
          WARNMSG("Default: Z = " << Z)
        }
        else
        {
          for(int i=0; i<gpset.size(); i++)
          {
            Element* cur = dynamic_cast<Element*>(gpset[i]);
            Attribute* att = cur->get_attribute("name");
            if(att)
            {
              const std::string& aname = att->get_value();
              xmlNode* curc = cur->cobj();
              if(aname == "Z")
                putContent(Z,curc);
              XMLReport("Z = " << Z)
            }
          }
        }
        // determine maximum angular momentum
        int lmax = 0;
        for(int i=0; i<Psi.L.size(); i++)
          if(Psi.L[i] > lmax)
            lmax = Psi.L[i];
        XMLReport("Maximum Angular Momentum = " << lmax)
        CG_coeff = new Clebsch_Gordan(lmax);
        Pot.add(new ZOverRPotential(Z), true);
        Pot.add(new HartreePotential(CG_coeff));
        Pot.add(new ExchangePotential(CG_coeff));
        Psi.CuspParam = Z;
      } // if Nuclear
      else
        if(pottype == "Harmonic")
        {
          XMLReport("Creating a Harmonic Potential.")
          double Omega=0.5;
          NodeSet gpset = potset[0]->find("./Parameter");
          if(gpset.empty())
          {
            Omega = Psi.NumOrb;
            WARNMSG("Default: Omega = " << Omega)
          }
          else
          {
            for(int i=0; i<gpset.size(); i++)
            {
              Element* cur = dynamic_cast<Element*>(gpset[i]);
              Attribute* att = cur->get_attribute("name");
              if(att)
              {
                const std::string& aname = att->get_value();
                xmlNode* curc = cur->cobj();
                if(aname == "Omega")
                  putContent(Omega,curc);
                else
                  if(aname == "Z")
                    putContent(Omega,curc);
                XMLReport("Omega = " << Omega)
              }
            }
          }
          // determine maximum angular momentum
          int lmax = 0;
          for(int i=0; i<Psi.L.size(); i++)
            if(Psi.L[i] > lmax)
              lmax = Psi.L[i];
          lmax++; // increment by 1
          XMLReport("Maximum Angular Momentum = " << lmax)
          CG_coeff = new Clebsch_Gordan(lmax);
          Pot.add(new HarmonicPotential(Omega),true);
          Pot.add(new HartreePotential(CG_coeff));
          Pot.add(new ExchangePotential(CG_coeff));
          Psi.CuspParam = 0.0;
        } // if Harmonic
        else
          if(pottype == "Step")
          {
            XMLReport("Creating a step-function Potential.")
            int lmax = 0;
            for(int i=0; i<Psi.L.size(); i++)
              if(Psi.L[i] > lmax)
                lmax = Psi.L[i];
            lmax++; // increment by 1
            XMLReport("Maximum Angular Momentum = " << lmax)
            CG_coeff = new Clebsch_Gordan(lmax);
            StepPotential* apot = new StepPotential;
            apot->put(dynamic_cast<Element*>(potset[0])->cobj());
            Pot.add(apot,true);
            Pot.add(new HartreePotential(CG_coeff));
            Pot.add(new ExchangePotential(CG_coeff));
            Psi.CuspParam = 0.0;
          } // steps
    }
  }
  return true;
}
}

// int main(int argc, char **argv) {
//  using namespace ohmmshf;
//  using namespace xmlpp;
//  DomParser parser;
//  parser.parse_file(argv[1]);
//  Node* root = parser.get_document()->get_root_node(); //deleted by DomParser.
//  std::string element;
//  std::string pottype;
//  std::string gridtype;
//  RadialPotentialSet Pot;
//  SphericalOrbitalTraits::BasisSetType Psi;
//  parseXMLFile(Pot,Psi,element,pottype,gridtype,root);
//  HartreeFock HFSolver(Pot,Psi,root);
//  HFSolver.solve(pottype,gridtype,Psi.size());
//  Psi.print(element);
//  return 1;
// }





