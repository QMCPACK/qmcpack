//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file BsplineReaderBase.cpp
 *
 * Implement super function
 */
#include <QMCWaveFunctions/EinsplineSetBuilder.h>
#include <QMCWaveFunctions/BsplineReaderBase.h>
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "qmc_common.h"
namespace qmcplusplus
{
  BsplineReaderBase::BsplineReaderBase(EinsplineSetBuilder* e)
    : mybuilder(e), MeshSize(0), myFirstSPO(0),myNumOrbs(0),GridFactor(1),Rcut(0)
  {
    myComm=mybuilder->getCommunicator();
  }

  void BsplineReaderBase::get_psi_g(int ti, int spin, int ib, Vector<std::complex<double> >& cG)
  {
    int ncg=0;
    if(myComm->rank()==0)
    {
      std::string path=psi_g_path(ti,spin,ib);
      HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
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

  BsplineReaderBase::~BsplineReaderBase() 
  {
  }

  inline std::string make_bandinfo_filename(const std::string& root, int spin, int twist, const Tensor<int,3>& tilematrix, int gid)
  {
    std::ostringstream oo;
    oo<<root 
      << ".tile_"
      << tilematrix(0,0) <<tilematrix(0,1) <<tilematrix(0,2)
      << tilematrix(1,0) <<tilematrix(1,1) <<tilematrix(1,2)
      << tilematrix(2,0) <<tilematrix(2,1) <<tilematrix(2,2)
      << ".spin_"<< spin << ".tw_" << twist
      ;
    if(gid>=0) oo << ".g"<<gid;
    return oo.str();
  }


  inline std::string make_bandgroup_name(const std::string& root, int spin, int twist, const Tensor<int,3>& tilematrix, int first, int last)
  {
    std::ostringstream oo;
    oo<<root 
      << ".tile_"
      << tilematrix(0,0) <<tilematrix(0,1) <<tilematrix(0,2)
      << tilematrix(1,0) <<tilematrix(1,1) <<tilematrix(1,2)
      << tilematrix(2,0) <<tilematrix(2,1) <<tilematrix(2,2)
      << ".spin_"<< spin << ".tw_" << twist
      <<".l"<<first<<"u"<<last;
    return oo.str();
  }

  void BsplineReaderBase::setCommon(xmlNodePtr cur)
  {
    OhmmsAttributeSet a;
    a.add(Rcut,"rmax_core");
    a.add(GridFactor,"dilation");
    a.put(cur);

    app_log() << "Rcut = " << Rcut << std::endl;
    app_log() << "dilation = " << GridFactor << std::endl;

  }

  SPOSetBase* BsplineReaderBase::create_spline_set(int spin, xmlNodePtr cur)
  {
    int ns(0);
    OhmmsAttributeSet a;
    a.add(ns,"size");
    a.put(cur);

    if(ns==0) 
      APP_ABORT_TRACE(__FILE__,__LINE__, "parameter/@size missing");

    if(spo2band.empty()) 
      spo2band.resize(mybuilder->states.size());

    std::vector<BandInfo>& fullband=(*(mybuilder->FullBands[spin]));

    if(spo2band[spin].empty())
    {
      spo2band[spin].reserve(fullband.size());
      if(mybuilder->states[spin]==0) mybuilder->states[spin]=new SPOSetInfo;
      mybuilder->clear_states(spin);
      initialize_spo2band(spin,fullband,*mybuilder->states[spin],spo2band[spin]);
    }

    BandInfoGroup vals;
    vals.TwistIndex=fullband[0].TwistIndex;
    vals.GroupID=0;
    vals.myName=make_bandgroup_name(mybuilder->getName(),spin,mybuilder->TwistNum,mybuilder->TileMatrix,0,ns);
    vals.selectBands(fullband,0, ns, false);

    size_t mem_now=qmc_common.memory_allocated;
    SPOSetBase* newspo=create_spline_set(spin,vals);       
    qmc_common.print_memory_change("BsplineSetReader", mem_now);
    return newspo;

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

  SPOSetBase* BsplineReaderBase::create_spline_set(int spin, xmlNodePtr cur, SPOSetInputInfo& input_info)
  {
    if(spo2band.empty()) 
      spo2band.resize(mybuilder->states.size());

    std::vector<BandInfo>& fullband=(*(mybuilder->FullBands[spin]));

    if(spo2band[spin].empty())
    {
      spo2band[spin].reserve(fullband.size());
      if(mybuilder->states[spin]==0) mybuilder->states[spin]=new SPOSetInfo;
      mybuilder->clear_states(spin);
      initialize_spo2band(spin,fullband,*mybuilder->states[spin],spo2band[spin]);
    }

    BandInfoGroup vals;
    vals.TwistIndex=fullband[0].TwistIndex;
    vals.GroupID=0;
    vals.myName=make_bandgroup_name(mybuilder->getName(),spin,mybuilder->TwistNum,mybuilder->TileMatrix
        ,input_info.min_index(),input_info.max_index());
    vals.selectBands(fullband,
        spo2band[spin][input_info.min_index()], 
        input_info.max_index()-input_info.min_index(),false);
    //vals.FirstSPO=0;
    //vals.NumSPOs=input_info.max_index()-input_info.min_index();
    size_t mem_now=qmc_common.memory_allocated;
    SPOSetBase* newspo=create_spline_set(spin,vals);       
    qmc_common.print_memory_change("BsplineSetReader", mem_now);
    return newspo;
  }

  /** build index tables to map a state to band with k-point folidng
   * @param bigspace full BandInfo constructed by EinsplineSetBuilder
   * @param sposet SPOSetInfo owned by someone, most likely EinsplinseSetBuilder
   * @param spo2band spo2band[i] is the index in bigspace
   *
   * At gamma or arbitrary kpoints with complex wavefunctions, spo2band[i]==i
   */
  void BsplineReaderBase::initialize_spo2band(int spin, const std::vector<BandInfo>& bigspace, SPOSetInfo& sposet, std::vector<int>& spo2band)
  {
    spo2band.reserve(bigspace.size());
    int ns=0;
    for(int i=0; i<bigspace.size(); ++i)
    {
      spo2band.push_back(i);
      SPOInfo a(ns,bigspace[i].Energy);
      sposet.add(a);
      ns++;
      if(bigspace[i].MakeTwoCopies)
      {
        spo2band.push_back(i);
        SPOInfo b(ns,bigspace[i].Energy);
        sposet.add(b);
        ns++;
      }
    }

    //write to a file
    const Communicate* comm=myComm;
    if(comm->rank()) return;

    std::string aname= make_bandinfo_filename(mybuilder->getName(),spin
        , mybuilder->TwistNum, mybuilder->TileMatrix,comm->getGroupID());
    aname+=".bandinfo.dat";

    std::ofstream o(aname.c_str());
    char s[1024];
    ns=0;
    typedef QMCTraits::PosType PosType;
    o << "#  Band    State   TwistIndex BandIndex Energy      Kx      Ky      Kz      K1      K2      K3    KmK " << std::endl;   
    for(int i=0; i<bigspace.size(); ++i)
    {
      int ti   = bigspace[i].TwistIndex;
      int bi   = bigspace[i].BandIndex;
      double e = bigspace[i].Energy;
      int nd = (bigspace[i].MakeTwoCopies)?2:1;
      PosType k= mybuilder->PrimCell.k_cart(mybuilder->TwistAngles[ti]);
      sprintf (s, "%8d %8d %8d %8d %12.6f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6d\n", 
          i, ns, ti, bi, e, k[0], k[1], k[2], 
          mybuilder->TwistAngles[ti][0], mybuilder->TwistAngles[ti][1], mybuilder->TwistAngles[ti][2],nd);
      o<<s;
      ns+=nd;
    }
  }
}

