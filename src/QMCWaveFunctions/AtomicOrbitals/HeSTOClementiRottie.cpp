//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
// -*- C++ -*-
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.h"
#include "Utilities/OhmmsInfo.h"

namespace ohmmsqmc {


  HePresetHFBuilder::HePresetHFBuilder(TrialWaveFunction& wfs,
				       ParticleSet& ions,
				       ParticleSet& els):
    OrbitalBuilderBase(wfs) 
  {   

    DistanceTableData* d_ei = DistanceTable::getTable(DistanceTable::add(ions,els));

    typedef DiracDeterminant<HePresetHF> Det_t;
    HePresetHF* heorb = new HePresetHF;
    heorb->setTable(d_ei);
    
    Det_t *DetU = new Det_t(*heorb,0);
    DetU->set(0,1);
    Det_t* DetD = new Det_t(*heorb,1);
    DetD->set(1,1);
    
    SlaterDeterminant<HePresetHF> *asymmpsi=new SlaterDeterminant<HePresetHF>;
    
    asymmpsi->add(DetU);
    asymmpsi->add(DetD);

    wfs_ref.add(asymmpsi);
  }

  bool 
  HePresetHFBuilder::put(xmlNodePtr cur){

    DistanceTableData* d_ei = NULL;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == dtable_tag) {
	string dtable((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
	dtable.append((const char*)(xmlGetProp(cur,(const xmlChar *)"target")));
	d_ei= DistanceTable::getTable(dtable.c_str());
	LOGMSG("DistanceTable is selected = " << dtable)
      }
      cur = cur->next;
    }
    typedef DiracDeterminant<HePresetHF> Det_t;
    HePresetHF* heorb = new HePresetHF;
    heorb->setTable(d_ei);

    Det_t *DetU = new Det_t(*heorb,0);
    DetU->set(0,1);
    Det_t* DetD = new Det_t(*heorb,1);
    DetD->set(1,1);

    SlaterDeterminant<HePresetHF> *asymmpsi=new SlaterDeterminant<HePresetHF>;

    asymmpsi->add(DetU);
    asymmpsi->add(DetD);

    wfs_ref.add(asymmpsi);
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
