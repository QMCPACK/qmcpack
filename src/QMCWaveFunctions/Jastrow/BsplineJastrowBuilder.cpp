#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "BsplineJastrowBuilder.h"
#include "BsplineFunctor.h"
#include "OneBodyJastrowFunction.h"

namespace qmcplusplus {

  bool
  BsplineJastrowBuilder::put(xmlNodePtr cur)
  {
    xmlNodePtr kids = cur->xmlChildrenNode;

    //set this jastrow function to be not optimizable
    //  if(targetPsi.VarList.size() == cur_var) {
    //    J1->setOptimizable(false);
    //  }


    OneBodyJastrow<BsplineFunctor<double> >* J1 = new OneBodyJastrow<BsplineFunctor<double> >(sourcePtcl, targetPtcl);
    J1->setOptimizable(true);

    while (kids != NULL) {
      std::string kidsname = (char*)kids->name;
      if (kidsname == "function") {
	BsplineFunctor<double> *functor = new BsplineFunctor<double>();
	functor->put (kids);
	functor->addOptimizables(targetPsi.VarList);
	
	// Find the number of the source species
	SpeciesSet &sSet = sourcePtcl.getSpeciesSet();
	int numSpecies = sSet.getTotalNum();
	int num = sSet.findSpecies (functor->elementType);
	if (num >= numSpecies) {
	  app_error() << "Unknown source species in BsplineJastrowBuilder: " 
		      << functor->elementType <<endl;
	  abort();
	}
	J1->addFunc (num, functor);
      }
      kids = kids->next;
    }

    targetPsi.addOrbital(J1);

    return true;
  }

}
