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
#include "AtomicHF/Clebsch_Gordan.h"
#include "AtomicHF/HartreeFock.h"
#include "AtomicHF/fillshells.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/OhmmsSpecies.h"
#include <fstream>
#include <libxml/xpath.h>

namespace ohmmshf
{

bool parseXMLFile(RadialPotentialSet& Pot, HFAtomicOrbitals& Psi,
                  std::string& name,std::string& pottype,std::string& gridtype,
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
                  << ", dh = " << grid->dh() << ", npts = " << npts)
      }
      else
        if(gridtype == "linear")
        {
          grid = new LinearGrid<double>;
          min/=scale;
          max/=sqrt(scale);
          grid->set(min,max,npts);
          XMLReport("rmin = " << min << ", rmax = " << max
                    << ", dh = " << grid->dh() << ", npts = " << npts)
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
    Clebsch_Gordan* CG_coeff;
    NodeSet potset = aset[0]->find("./Potential");
    if(potset.empty())
    {
      ERRORMSG("Must provide potential information.")
      return false;
    }
    else
    {
      Element* cur = dynamic_cast<Element*>(potset[0]);
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
          if(pottype == "Harmonic")
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
        Pot.add(new ZOverRFunctor(Z));
        Pot.add(new HartreePotential(CG_coeff));
        Pot.add(new ExchangePotential(CG_coeff));
        Psi.CuspParam = Z;
        Psi.EigenBoundParam = Z;
      } // if Nuclear
      else
        if(pottype == "Harmonic")
        {
          XMLReport("Creating a Harmonic Potential.")
          double Omega;
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
          Pot.add(new HarmonicFunctor(Omega));
          Pot.add(new HartreePotential(CG_coeff));
          Pot.add(new ExchangePotential(CG_coeff));
          Psi.CuspParam = 0.0;
          Psi.EigenBoundParam = Omega;
        } // if Harmonic
        else
          if(pottype == "Nuclear_Scalar_Rel")
          {
            XMLReport("Creating a Nuclear Scalar Relativistic Potential.")
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
            Pot.add(new ZOverRFunctor(Z));
            Pot.add(new HartreePotential(CG_coeff));
            Pot.add(new ExchangePotential(CG_coeff));
            Psi.CuspParam = Z;
            Psi.EigenBoundParam = Z;
          } // if Nuclear_Scalar_Rel
          else
            if(pottype == "Pseudo")
            {
              XMLReport("Creating a Starkloff-Joannopoulos Pseudo-Potential.")
              double rc, lambda, Zeff;
              NodeSet gpset = potset[0]->find("./Parameter");
              if(gpset.empty())
              {
                ERRORMSG("Must Provide core radius and decay for PseudoPotential.")
                return false;
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
                    if(aname == "rc")
                    {
                      putContent(rc,curc);
                      XMLReport("Core radius = " << rc)
                    }
                    else
                      if(aname == "lambda")
                      {
                        putContent(lambda,curc);
                        XMLReport("Decay = " << lambda)
                      }
                      else
                        if(aname == "Zeff")
                        {
                          putContent(Zeff,curc);
                          XMLReport("Effective core charge = " << Zeff)
                        }
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
              Pot.add(new PseudoPotential(Zeff,rc,lambda));
              Pot.add(new HartreePotential(CG_coeff));
              Pot.add(new ExchangePotential(CG_coeff));
              Psi.CuspParam = 0.0;
              Psi.EigenBoundParam = 2.0*Zeff/rc;
            } // if Pseudo
    }
  }
  return true;
}
}


int main(int argc, char **argv)
{
  using namespace ohmmshf;
  using namespace xmlpp;
  DomParser parser;
  parser.parse_file(argv[1]);
  Node* root
  = parser.get_document()->get_root_node(); //deleted by DomParser.
  std::string element;
  std::string pottype;
  std::string gridtype;
  RadialPotentialSet Pot;
  HFAtomicOrbitals Psi;
  parseXMLFile(Pot,Psi,element,pottype,gridtype,root);
  HartreeFock HFSolver(Pot,Psi,root);
  HFSolver.solve(pottype,gridtype,Psi.size());
  Psi.print(element);
  return 1;
}





