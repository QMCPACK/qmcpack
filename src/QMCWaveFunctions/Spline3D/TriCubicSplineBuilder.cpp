#include "OhmmsData/libxmldefs.h"
#include "QMCWaveFunctions/TriCubicSplineBuilder.h"
#include "QuantumSystem/AnalyticOrbitalSet.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus {

  bool TriCubicSplineBuilder::put(xmlNodePtr cur) {

    
    if(!m_orbitals) m_orbitals = new TriCubicSplineSet;

    typedef DiracDeterminant<SPOSet_t> DDet_t;
    typedef SlaterDeterminant<SPOSet_t> SDet_t;

    map<string,int> SplineID;
    int nbasis=0;

    SDet_t* sdet = new SDet_t;
    /// go through the tree cur = DeterminantSet
    cur = cur->xmlChildrenNode;
    while(cur!=NULL){
      string cname((const char*)(cur->name));
      if(cname == "BasisSet") {
	m_orbitals->put(cur,grid_ref);  /// read in qGridi and qGridf,
			       /// initialise Schr_Grid         
      } else if(cname =="SlaterDeterminant") {
	int first = 0;
	xmlNodePtr tcur = cur->xmlChildrenNode;
	while(tcur != NULL) {
	  string cname2((const char*)(tcur->name));
	  if(cname2 == "Determinant") {
	    SPOSet_t* swfs = new SPOSet_t;
	    int norb = atoi
	      ((const char*)(xmlGetProp(tcur, (const xmlChar *)"orbitals"))) ;
	    XMLReport
	      ("The number of orbitals for a Dirac Determinant " << norb)
	      XMLReport("The number of basis functions " << swfs->size())
	    xmlNodePtr t = tcur->xmlChildrenNode;
	    while(t != NULL) {
	      string oname((const char*)(t->name));
	      if(oname == "Orbital") {
		/// ??? what is this doing ???
		TriCubicSpline* neworb = m_orbitals->getOrbital
		  ((const char*)xmlGetProp(t,(xmlChar*)"name"));
		if(neworb) swfs->add(neworb);
	      }
	      t = t->next;
	    }
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
    //    grid_ref = m_orbitals->getFullGrid();   /// ?? why needed??

    return true;	
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
