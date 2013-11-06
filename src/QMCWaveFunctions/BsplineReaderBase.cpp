//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file BsplineReaderBase.cpp
 *
 * Implement super function
 */
#include <QMCWaveFunctions/EinsplineSetBuilder.h>
#include <QMCWaveFunctions/BsplineReaderBase.h>
#include <QMCWaveFunctions/SPOSetComboNoCopy.h>
#include "Message/CommOperators.h"
namespace qmcplusplus
{
  void BsplineReaderBase::get_psi_g(int ti, int spin, int ib, Vector<complex<double> >& cG)
  {
    int ncg=0;
    if(myComm->rank()==0)
    {
      string path=psi_g_path(ti,spin,ib);
      HDFAttribIO<Vector<complex<double> > >  h_cG(cG);
      h_cG.read (mybuilder->H5FileID, path.c_str());
      ncg=cG.size();
    }
    myComm->bcast(ncg);
    if(ncg != mybuilder->MaxNumGvecs)
    {
      APP_ABORT("Failed : ncg != MaxNumGvecs");
    }
    myComm->bcast(cG);
  }

  SPOSetBase* BsplineReaderBase::create_spline_set(int spin, EinsplineSet* orbitalset)
  {
    vector<BandInfo>& SortBands=mybuilder->SortBands;
    BandInfoGroup vals;
    vals.GroupID=0;
    vals.selectBands(mybuilder->SortBands,0,mybuilder->NumDistinctOrbitals);
    return create_spline_set(spin,orbitalset,vals);

    //Test SPOSetComboNoCopy that can have multiple SPOSets
    //SPOSetComboNoCopy* bb=new SPOSetComboNoCopy;
    ////create SPOSet for the semi core states e=(-1000,-3.0)
    //BandInfoGroup cores0;
    //cores0.selectBands(mybuilder->SortBands,-1000.0,-3.0);
    //cores0.GroupID=0;
    //SPOSetBase* bandone=create_spline_set(spin,orbitalset,cores0);
    //
    ////create SPOSet for the rest with a coarse grid
    //TinyVector<int,3> mesh_saved=MeshSize;
    //for(int i=0; i<3; ++i) MeshSize[i] -= mesh_saved[i]/4;
    //BandInfoGroup cores1;
    //cores1.selectBands(mybuilder->SortBands,cores0.getNumDistinctOrbitals(),mybuilder->NumDistinctOrbitals);
    //cores1.GroupID=1;
    //SPOSetBase* bandtwo=create_spline_set(spin,orbitalset,cores1);
    //
    ////add them to bb
    //bb->add(bandone);
    //bb->add(bandtwo);
    //bb->setOrbitalSetSize(orbitalset->getOrbitalSetSize());
    //return bb;
  }
}

