//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#include "Utilities/OhmmsInfo.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/HDFNumericAttrib.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/RadialGridFunctorBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.h"

namespace ohmmsqmc {

  int DFactorial(int l) {
    if(l == 1) 
      return 1;
    else 
      return DFactorial(l-2);
  }

  template<class T>
  GaussianCombo<T>::GaussianCombo(int l){
    L = static_cast<T>(l);
    //Everything related to L goes to NormL and NormPow
    NormL = pow(2,L+1)*sqrt(2.0/static_cast<T>(DFactorial(2*l+1)))*pow(2.0/M_PI,0.25);
    NormPow = 0.5*(L+1.0)+0.25;
  }

  template<class T>
  bool GaussianCombo<T>::put(xmlNodePtr cur) {
    InParam.push_back(cur);
  }

  template<class T>
  void GaussianCombo<T>::reset() {
    int n=gset.size();
    while(n<InParam.size()) {
      gset.push_back(new BasicGaussian<T>());
      n++;
    }
    for(int i=0; i<InParam.size(); i++) {
      const xmlChar* aptr = xmlGetProp(InParam[i],(const xmlChar *)"alpha");
      const xmlChar* cptr = xmlGetProp(InParam[i],(const xmlChar *)"c");
      if(aptr == 0 || cptr == 0) continue;
      T alpha = atof((const char*)aptr);
      T c = atof((const char*)cptr);
      T norm = c*NormL*pow(alpha,NormPow); //get the normalization factor
      gset[i]->reset(alpha,norm);
    }
  }


  bool
  GTO2GridBuilder::addRadialOrbital(xmlNodePtr cur, 
				    const QuantumNumberType& nlms) {

    RadialOrbitalType *radorb =  NULL;

    int n=nlms[0];
    int l=nlms[1];

    GaussianCombo<ValueType> gaussian(l);
 
    xmlNodePtr s = cur->xmlChildrenNode;
    while(s != NULL) {
      string cname((const char*)(s->name));
      if(cname == "Rnl") {
        gaussian.put(s);
      }
      s = s->next;
    }

    gaussian.reset();

    //pointer to the grid
    GridType* agrid = m_orbitals->Grids[0];
    radorb = new OneDimCubicSpline<ValueType>(agrid);

    //spline the slater type orbital
    Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gaussian, *radorb);
    transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());

    //add the radial orbital to the list
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(nlms);

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
