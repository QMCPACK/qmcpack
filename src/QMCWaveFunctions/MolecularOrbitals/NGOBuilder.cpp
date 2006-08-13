//////////////////////////////////////////////////////////////////
// (c) Copyright 2003 by Jeongnim Kim and Jordan Vincent
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
#include "Numerics/GaussianBasisSet.h"
#include "Numerics/SlaterBasisSet.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCFactory/OneDimGridFactory.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
namespace qmcplusplus {

  NGOBuilder::NGOBuilder(xmlNodePtr cur): 
    Normalized(true),m_rcut(-1.0){
      if(cur != NULL) {
        putCommon(cur);
      }
  }

  bool 
    NGOBuilder::putCommon(xmlNodePtr cur) {
    const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
    if(a) {
      if(xmlStrEqual(a,(const xmlChar*)"no")) Normalized=false;
    }
    return true;
  }

  void 
    NGOBuilder::setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter) { 
    m_orbitals = oset;
    m_species = acenter;
  }

  bool 
    NGOBuilder::addGrid(xmlNodePtr cur) {
    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }
    GridType *agrid = OneDimGridFactory::createGrid(cur);
    m_orbitals->Grids.push_back(agrid);
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
    NGOBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms) {

    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    const xmlChar *tptr = xmlGetProp(cur,(const xmlChar*)"type");
    string radtype("Gaussian");
    if(tptr) radtype = (const char*)tptr;

    tptr = xmlGetProp(cur,(const xmlChar*)"rmax");
    if(tptr) m_rcut = atof((const char*)tptr);

    int lastRnl = m_orbitals->Rnl.size();

    m_nlms = nlms;
    if(radtype == "Gaussian") {
      addGaussian(cur);
    } else if(radtype == "Slater") {
      addSlater(cur);
    } else if(radtype == "Pade") {
      addPade(cur);
    }

    if(lastRnl && m_orbitals->Rnl.size()> lastRnl) {
      //LOGMSG("\tSetting GridManager of " << lastRnl << " radial orbital to false")
      m_orbitals->Rnl[lastRnl]->setGridManager(false);
    }

    return true;
  }

  void NGOBuilder::addGaussian(xmlNodePtr cur) {

    int L= m_nlms[1];
    GaussianCombo<RealType> gset(L,Normalized);
    gset.putBasisGroup(cur);

    GridType* agrid = m_orbitals->Grids[0];
    RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);

    if(m_rcut<0) m_rcut = agrid->rmax();
    Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gset, *radorb);
    transform.generate(agrid->rmin(),m_rcut,agrid->size());

    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(m_nlms);
  }

  void NGOBuilder::addSlater(xmlNodePtr cur) {

    ////pointer to the grid
    GridType* agrid = m_orbitals->Grids[0];
    RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);

    SlaterCombo<RealType> sto(m_nlms[1],Normalized);
    sto.putBasisGroup(cur);

    //spline the slater type orbital
    Transform2GridFunctor<SlaterCombo<RealType>,RadialOrbitalType> transform(sto, *radorb);
    transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());
    //transform.generate(agrid->rmax());
    //add the radial orbital to the list
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(m_nlms);

  }

  template<class T>
  struct PadeOrbital: public RadialOrbitalBase<T> {
  
    typedef T value_type;
    T a0,a1,a2,a3,rcut;
    std::string  nodeName;
  
    explicit 
      PadeOrbital(const char* node_name="radfunc"):
        a0(1.0),a1(-0.5),a2(0.0),a3(-0.2),rcut(4.0),nodeName(node_name){}
  
    ~PadeOrbital(){ }
  
    void reset() {}
  
    inline value_type f(value_type r) const {
      return a0*std::exp((a1+a2*r)*r/(1.0e0+a3*r));
    }
  
    inline value_type df(value_type r) const {
      value_type t = 1.0/(1.0e0+a3*r);
      value_type z=(a1+a2*r)*r*t;
      value_type res = a0*std::exp(z);
      return res*(a1+2.0*a2*r-z*a3)*t;
    }
  
    bool putBasisGroup(xmlNodePtr cur) {
      cur = cur->children;
      while(cur != NULL) {
        string cname((const char*)cur->name);
        if(cname == nodeName) {
          OhmmsAttributeSet radAttrib;
          radAttrib.add(a0,"a0"); radAttrib.add(a1,"a1"); 
          radAttrib.add(a2,"a2"); radAttrib.add(a3,"a3");
          radAttrib.put(cur);
          rcut = -1.0/a3-std::numeric_limits<T>::epsilon();
          LOGMSG("Pade compoenent (a0,a1,a2,a3,rcut) = " << a0 << " " << a1 << " " << a2 << " " << a3 << " " << rcut)
        }
        cur=cur->next;
      }
      return true;
    }
  };

  void NGOBuilder::addPade(xmlNodePtr cur) {

    GridType* agrid = m_orbitals->Grids[0];
    RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);

    PadeOrbital<RealType> pade;
    pade.putBasisGroup(cur);

    //spline the slater type orbital
    Transform2GridFunctor<PadeOrbital<RealType>,RadialOrbitalType> transform(pade, *radorb);
    if(pade.rcut>0) 
      transform.generate(pade.rcut);
    else 
      transform.generate(agrid->rmax());
    //add the radial orbital to the list
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(m_nlms);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
