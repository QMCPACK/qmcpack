//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/GaussianBasisSet.h"
#include "QMCTools/GridMolecularOrbitals.h"
#include "QMCTools/GTO2GridBuilder.h"
#include "QMCFactory/OneDimGridFactory.h"

namespace qmcplusplus
{
bool
GTO2GridBuilder::putCommon(xmlNodePtr cur)
{
  const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
  if(a)
  {
    if(xmlStrEqual(a,(const xmlChar*)"no"))
      Normalized=false;
  }
  return true;
}

/** Process basisGroup node
 * @param cur current xml node
 * @param nlms Quantum Numbers
 * @return true when successful
 *
 * cur is the basisGroup node, basisGroup := (radfunc)+
 * For basisGroup with l=0 and multiple radfunc (contracted S),
 * the exponent of a Slater-type orbital is provided as an alternative.
 */
bool
GTO2GridBuilder::addRadialOrbital(xmlNodePtr cur,
                                  const QuantumNumberType& nlms)
{
  int n=nlms[0];
  int l=nlms[1];
  std::string b_name((const char*)xmlGetProp(cur,(const xmlChar*)"rid"));
  //Using default <radfunc exponent="alpha" contraction="c"/>
  GaussianCombo<RealType> gaussian(l,Normalized);
  gaussian.putBasisGroup(cur);
  //pointer to the grid
  GridType* agrid = m_orbitals->Grids[0];
  RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);
  //spline the slater type orbital
  Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gaussian, *radorb);
  transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());
  //add the radial orbital to the list
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(nlms);
  //guess the exponent of a corresponding STO of the contracted S basisGroup
  if(l ==0 && gaussian.size()>1)
  {
    RealType peak(0.0);
    int gMax=0;
    for(int ig=1; ig<agrid->size()-1; ig++)
    {
      RealType r=(*agrid)(ig);
      RealType y=(*radorb)(ig);
      RealType mag=r*r*y*y;
      if(mag>peak)
      {
        peak=mag;
        gMax=ig;
      }
    }
    std::ostringstream slater;
    slater<<" Possible substitution "<<b_name <<" by a Slater-type orbital\n";
    slater<<"  <basisGroup rid=\""<<b_name<<"\" n=\""<<n<<"\" l=\""<<l<<"\" type=\"Slater\">\n";
    slater<<"    <radfunc exponent=\""<<1.0/(*agrid)(gMax)<<"\" contraction=\"1.0\"/>\n";
    slater<<"  </basisGroup>\n";
    xmlAddPrevSibling(cur,xmlNewComment((const xmlChar*)slater.str().c_str()));
  }
  return true;
}

/** Default function to add a radial grid to the list of radial grids.
 * \param cur the current xmlNode to be processed
 * \return true if succeeds
 */
bool
GTO2GridBuilder::addGrid(xmlNodePtr cur)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  XMLReport("Converting analytic orbitals to radial grid functions. Modify to use zero-based grid.")
  m_orbitals->Grids.push_back(OneDimGridFactory::createGrid(cur));
  return true;
}
}
