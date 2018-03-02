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
    
    
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/GaussianBasisSet.h"
#include "Numerics/SlaterBasisSet.h"
#include "QMCWaveFunctions/Jastrow/CBSOBuilder.h"
#include "QMCWaveFunctions/Jastrow/WMFunctor.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

CBSOBuilder::CBSOBuilder(xmlNodePtr cur):
  Normalized(true),m_rcut(5.0),NumGridPoints(51)
{
  if(cur != NULL)
  {
    putCommon(cur);
  }
}

bool
CBSOBuilder::putCommon(xmlNodePtr cur)
{
  return true;
}

void
CBSOBuilder::setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter)
{
  m_orbitals = oset;
  m_species = acenter;
}

bool
CBSOBuilder::addGrid(xmlNodePtr cur)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  RealType delta=-1.0;
  OhmmsAttributeSet attrib;
  attrib.add(m_rcut,"rf");
  attrib.add(NumGridPoints,"npts");
  attrib.add(delta,"step");
  attrib.put(cur);
  if(delta>0)
    NumGridPoints = static_cast<int>(m_rcut/delta)+1;
  app_log() << "   CBSOBuilder::addGrid Rcut = " << m_rcut << "  NumGridPoints = " <<  NumGridPoints << std::endl;
  return true;
}

/** Add a new Slater Type Orbital with quantum numbers \f$(n,l,m,s)\f$
 * \param cur  the current xmlNode to be processed
 * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
 * \return true is succeeds
 *
 This function puts the STO on a logarithmic grid and calculates the boundary
 conditions for the 1D Cubic Spline.  The derivates at the endpoint
 are assumed to be all zero.  Note: for the radial orbital we use
 \f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
 radial orbital and \f$ l \f$ is the angular momentum.
*/
bool
CBSOBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  std::string radtype("Gaussian");
  OhmmsAttributeSet attrib;
  attrib.add(radtype,"type");
  attrib.put(cur);
  int lastRnl = m_orbitals->Rnl.size();
  m_nlms = nlms;
  int L= m_nlms[1];
  typedef RadialOrbitalType::FNIN InputFuncType;
  InputFuncType* infunc=0;
  //create one
  if(radtype == "Gaussian" || radtype == "GTO")
  {
    GaussianCombo<RealType> *gset = new GaussianCombo<RealType>(L,Normalized);
    gset->putBasisGroup(cur);
    app_log() << "   " << L << " basisGroup  contains " << gset->size() << " radial functors." << std::endl;
    infunc=gset;
  }
  else
    if(radtype == "Slater")
    {
      SlaterCombo<RealType> *sto_set=new SlaterCombo<RealType>(L,Normalized);
      sto_set->putBasisGroup(cur);
      infunc=sto_set;
    }
    else
      if(radtype == "WM")
      {
        xmlNodePtr tcur=cur->children;
        int nr=0;
        while(tcur != NULL)
        {
          std::string cname((const char*)(tcur->name));
          if(cname == "radfunc")
          {
            OhmmsAttributeSet rAttrib;
            RealType rcut=m_rcut;//use the global cutoff but can overwrite
            rAttrib.add(rcut,"rcut");
            rAttrib.put(tcur);
            infunc = new WMFunctor<RealType>(1.0,rcut);
            infunc->put(tcur);
            nr++;
          }
          tcur=tcur->next;
        }
      }
  if(infunc)
  {
    app_log()
        << "   CBSOBuilder::addRadialOrbital  input " << radtype
        << " output = CubicBspline " << std::endl;
    RadialOrbitalType *radorb=new RadialOrbitalType;
    radorb->initialize(infunc,m_rcut,NumGridPoints);
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(m_nlms);
    //if(lastRnl && m_orbitals->Rnl.size()> lastRnl) {
    //  //LOGMSG("\tSetting GridManager of " << lastRnl << " radial orbital to false")
    //  m_orbitals->Rnl[lastRnl]->setGridManager(false);
    //}
  }
  return true;
}

}
