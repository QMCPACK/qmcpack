//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   e-mail: jeongnim.kim@gmail.com
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file BandInfo.cpp
 */
#include <QMCWaveFunctions/BandInfo.h>
#include <QMCWaveFunctions/SPOSetInfo.h>
namespace qmcplusplus
{
  BandInfoGroup::BandInfoGroup(): FirstSPO(0), NumSPOs(0), FirstBand(0) {}

  void BandInfoGroup::selectBands(const vector<BandInfo>& bigspace, double emin, double emax)
  {
    myBands.clear();
    if(emin>emax) return;

    int iorb=0;
    int N=bigspace.size();
    bool skipit=true;
    int n_lower=0; 
    do
    {
      if(bigspace[iorb].Energy>=emin) break;
      n_lower += (bigspace[iorb].MakeTwoCopies)?2:1;
      ++iorb;
    } while(iorb<N);

    if(iorb>=N) 
    {
      APP_ABORT("BandInfoGroup::selectBands failed due to iorb>=N");
    }

    FirstSPO=n_lower;
    FirstBand=iorb;
    NumSPOs=0;
    while(iorb<N)
    {
      if(bigspace[iorb].Energy>=emax) break;
      myBands.push_back(bigspace[iorb]);
      NumSPOs += (bigspace[iorb].MakeTwoCopies)?2:1;
      ++iorb;
    }

    app_log() << "BandInfoGroup::selectBands using energy ["<< emin << "," << emax << ")"  << endl;
    app_log() << "  Number of distinct bands " << myBands.size() << endl;
    app_log() << "  First Band index " << FirstBand << endl;
    app_log() << "  First SPO index " << FirstSPO << endl;
    app_log() << "  Size of SPOs " << NumSPOs << endl;

    //for(int i=0; i<myBands.size(); ++i)
    //  app_log() << myBands[i].TwistIndex << " " << myBands[i].Energy << endl;
    app_log().flush();
  }

  void BandInfoGroup::selectBands(const vector<BandInfo>& bigspace, 
      int first_orb, int num_spos, bool relative)
  {

    app_log() << "BandInfoGroup::selectBands bigspace has " << bigspace.size() 
      << " distinct orbitals " << endl;
    myBands.clear();

    int iorb=0;
    int N=bigspace.size();
    bool skipit=true;
    int n_lower=0; 
    do
    {
      if(iorb>=first_orb) break;
      n_lower += (bigspace[iorb].MakeTwoCopies)?2:1;
      ++iorb;
    } while(iorb<N);

    if(iorb>=N) 
    {
      APP_ABORT("BandInfoGroup::selectBands failed due to iorb>=N");
    }

    FirstSPO=(relative)?n_lower:0;
    FirstBand=iorb;
    NumSPOs=0;
    int ns_max=num_spos-1;
    while(iorb<N && NumSPOs<num_spos)
    {
      //if(myBands.size()>=num_orbs) break;
      myBands.push_back(bigspace[iorb]);
      NumSPOs += (NumSPOs<ns_max && bigspace[iorb].MakeTwoCopies)?2:1;
      ++iorb;
    }

    app_log() << "BandInfoGroup::selectBands using distinct orbitals ["<< first_orb << "," << iorb << ")"  << endl;
    app_log() << "  Number of distinct bands " << myBands.size() << endl;
    app_log() << "  First Band index " << FirstBand << endl;
    app_log() << "  First SPO index " << FirstSPO << endl;
    app_log() << "  Size of SPOs " << NumSPOs << endl;
  }

//  void BandInfoGroup::selectBands(const vector<BandInfo>& bigspace, int first_orb, int last_orb)
//  {
//    app_log() << "BandInfoGroup::selectBands bigspace has " << bigspace.size() << " distinct orbitals " << endl;
//    myBands.clear();
//
//    int iorb=0;
//    int N=bigspace.size();
//    bool skipit=true;
//    int n_lower=0; 
//    do
//    {
//      if(iorb>=first_orb) break;
//      n_lower += (bigspace[iorb].MakeTwoCopies)?2:1;
//      ++iorb;
//    } while(iorb<N);
//
//    if(iorb>=N) 
//    {
//      APP_ABORT("BandInfoGroup::selectBands failed due to iorb>=N");
//    }
//
//    if(first_orb != iorb)
//    {
//      APP_ABORT("Cannot locate the first SPO ");
//    }
//
//    int num_orbs=last_orb-first_orb;
//    if(num_orbs<=0)
//    {
//      APP_ABORT("BandInfoGroup::selectBands(bigspace,first_orb,last_orb) Illegal range ");
//    }
//
//    FirstSPO=n_lower;
//    FirstBand=iorb;
//    NumSPOs=0;
//    while(iorb<N)
//    {
//      if(myBands.size()>=num_orbs) break;
//      myBands.push_back(bigspace[iorb]);
//      NumSPOs += (bigspace[iorb].MakeTwoCopies)?2:1;
//      ++iorb;
//    }
//
//    //for(int i=0; i<myBands.size(); ++i)
//    //  app_log() << myBands[i].TwistIndex << " " << myBands[i].Energy << endl;
//
//    app_log().flush();
//  }
}

