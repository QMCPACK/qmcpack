#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "BsplineJastrowBuilder.h"
#include "BsplineFunctor.h"
#include "OneBodyJastrowOrbital.h"
#include "TwoBodyJastrowOrbital.h"


namespace qmcplusplus {

  bool
  BsplineJastrowBuilder::put(xmlNodePtr cur)
  {
    xmlNodePtr kids = cur->xmlChildrenNode;

    //set this jastrow function to be not optimizable
    //  if(targetPsi.VarList.size() == cur_var) {
    //    J1->setOptimizable(false);
    //  }

    // Create a one-body Jastrow
    if (sourcePtcl != NULL) {
      OneBodyJastrowOrbital<BsplineFunctor<double> > *J1 = 
	new OneBodyJastrowOrbital<BsplineFunctor<double> >(*sourcePtcl, targetPtcl);
      while (kids != NULL) {
	std::string kidsname = (char*)kids->name;
	if (kidsname == "correlation") {
	  BsplineFunctor<double> *functor = new BsplineFunctor<double>();
	  functor->put (kids);
	  functor->addOptimizables(targetPsi.VarList);
	  
	  // Find the number of the source species
	  SpeciesSet &sSet = sourcePtcl->getSpeciesSet();
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
      J1->setOptimizable(true);
    } 
    // Create a two-body Jastrow
    else {
      TwoBodyJastrowOrbital<BsplineFunctor<double> > *J2 = 
	new TwoBodyJastrowOrbital<BsplineFunctor<double> >(targetPtcl);
      
      std::map<std::string,BsplineFunctor<double>*> functorMap;
      while (kids != NULL) {
	std::string kidsname = (char*)kids->name;
	if (kidsname == "correlation") {
	  BsplineFunctor<double> *functor = new BsplineFunctor<double>();
	  functor->put (kids);
	  functor->addOptimizables(targetPsi.VarList);
	  string pairType = functor->pairType;
	  if ((pairType != "uu") && (pairType != "ud") &&
	      (pairType != "dd") && (pairType != "dd")) {
	    app_error() << "Unrecognized pair type " << pairType 
			<< " in BsplineJastrowBuilder.\n";
	    abort();
	  }
	  functorMap[pairType] = functor;
	  J2->insert (pairType, functor);
	}
	kids = kids->next;
      }
      BsplineFunctor<double> *Juu, *Jud, *Jdu, *Jdd;
      std::map<std::string,BsplineFunctor<double>*>::iterator iter;

      // Find up-up section
      iter = functorMap.find("uu");
      if (iter == functorMap.end()) iter = functorMap.find("dd");
      if (iter == functorMap.end()) {
	app_error() << "You must specify a correlation section for pairType=\"uu\" or \"dd\".\n";
	abort();
      }
      Juu = iter->second;
      
      // Find up-down section
      iter = functorMap.find("ud");
      if (iter == functorMap.end()) iter = functorMap.find("du");
      if (iter == functorMap.end()) {
	app_error() << "You must specify a correlation section for pairType=\"ud\" or \"du\".\n";
	abort();
      }
      Jud = iter->second;

      // Find down-up section
      iter = functorMap.find("du");
      if (iter == functorMap.end()) iter = functorMap.find("ud");
      if (iter == functorMap.end()) {
	app_error() << "You must specify a correlation section for pairType=\"du\" or \"ud\".\n";
	abort();
      }
      Jdu = iter->second;
      
      // Find down-down section
      iter = functorMap.find("dd");
      if (iter == functorMap.end()) iter = functorMap.find("uu");
      if (iter == functorMap.end()) {
	app_error() << "You must specify a correlation section for pairType=\"uu\" or \"dd\".\n";
	abort();
      }
      Jdd = iter->second;
      


      J2->addFunc(Juu);
      J2->addFunc(Jud);
      J2->addFunc(Jdu);
      J2->addFunc(Jdd);

      targetPsi.addOrbital(J2);
      J2->setOptimizable(true);
    }

    return true;
  }

}
