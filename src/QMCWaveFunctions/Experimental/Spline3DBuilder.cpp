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
    
    

#include "OhmmsData/libxmldefs.h"
#include "QMCWaveFunctions/Spline3D/Spline3DBuilder.h"
#include "QuantumSystem/AnalyticOrbitalSet.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

#include <iostream>

namespace qmcplusplus
{

//Spline3DSet Spline3DBuilder::orbitals;

bool
Spline3DBuilder::put(xmlNodePtr cur)
{
  if(!d_orbitals)
    d_orbitals = new Spline3DSet;
  //typedef DiracDeterminant<Spline3DSet> DDet_t;
  //typedef SlaterDeterminant<Spline3DSet> SDet_t;
  typedef DiracDeterminant<SPOSet_t> DDet_t;
  typedef SlaterDeterminant<SPOSet_t> SDet_t;
  std::map<std::string,int> SplineID;
  int nbasis=0;
  SDet_t* sdet = new SDet_t;
  ///go thru the tree
  ///cur = DeterminantSet
  cur = cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "BasisSet")
    {
      d_orbitals->put(cur);
    }
    else
      if(cname =="SlaterDeterminant")
      {
        int first = 0;
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL)
        {
          std::string cname2((const char*)(tcur->name));
          if(cname2 == "Determinant")
          {
            SPOSet_t* swfs = new SPOSet_t;
            int norb
            =atoi((const char*)(xmlGetProp(tcur, (const xmlChar *)"orbitals"))) ;
            XMLReport("The number of orbitals for a Dirac Determinant " << norb)
            XMLReport("The number of basis functions " << swfs->size())
            xmlNodePtr t = tcur->xmlChildrenNode;
            while(t != NULL)
            {
              std::string oname((const char*)(t->name));
              if(oname == "Orbital")
              {
                Spline3D* neworb
                = d_orbitals->getOrbital((const char*)xmlGetProp(t,(xmlChar*)"name"));
                if(neworb)
                  swfs->add(neworb);
                //int iorb = atoi((const char*)xmlGetProp(t,(xmlChar*)"index"));
                //                swfs->add(orbitals.getOrbital(iorb));
              }
              t = t->next;
            }
            // for(xmlNodePtr t=cur1->xmlChildrenNode; t != NULL; t=t->next) {
            //  std::string awfs((const char*)xmlGetProp(t,(xmlChar*)"name"));
            //  int iorb= SplineID[awfs];
            DDet_t* ddet = new DDet_t(*swfs,first);
            ddet->set(first,swfs->size());
            XMLReport("Adding a determinant to the SlaterDeterminant of size "
                      << swfs->size() << " begining with orbital " << first)
            first += swfs->size();
            sdet->add(ddet);
          } //dirac-determinant
          tcur = tcur->next;
        }
      }//slater-determinant
    cur = cur->next;
  }
  wfs_ref.add(sdet);
  grid_ref = d_orbitals->getFullGrid();
  return true;
}
}
