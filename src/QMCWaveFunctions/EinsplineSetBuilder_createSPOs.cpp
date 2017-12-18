//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Numerics/e2iphi.h"
#include "simd/vmath.hpp"
#include "qmc_common.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
#include "Numerics/HDFSTLAttrib.h"
#include "ParticleIO/ESHDFParticleParser.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Particle/DistanceTableData.h"
#include <fftw3.h>
#include <Utilities/ProgressReportEngine.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineReaderBase.h"
#include "QMCWaveFunctions/EinsplineAdoptor.h"

namespace qmcplusplus
{

  ///create R2R, real wavefunction in double
  BsplineReaderBase* createBsplineRealDouble(EinsplineSetBuilder* e, bool hybrid_rep);
  ///create R2R, real wavefunction in float
  BsplineReaderBase* createBsplineRealSingle(EinsplineSetBuilder* e, bool hybrid_rep);
  ///create C2C or C2R, complex wavefunction in double
  BsplineReaderBase* createBsplineComplexDouble(EinsplineSetBuilder* e, bool hybrid_rep);
  ///create C2C or C2R, complex wavefunction in single
  BsplineReaderBase* createBsplineComplexSingle(EinsplineSetBuilder* e, bool hybrid_rep);
  ///disable truncated orbitals for now
  BsplineReaderBase* createTruncatedSingle(EinsplineSetBuilder* e, int celltype)
  {
    return nullptr;
  }
  BsplineReaderBase* createTruncatedDouble(EinsplineSetBuilder* e, int celltype)
  {
    return nullptr;
  }

void EinsplineSetBuilder::set_metadata(int numOrbs, int TwistNum_inp)
{
  // 1. set a lot of internal parameters in the EinsplineSetBuilder class
  //  e.g. TileMatrix, UseRealOrbitals, DistinctTwists, MakeTwoCopies.
  // 2. this is also where metadata for the orbitals are read from the wavefunction hdf5 file
  //  and broacasted to MPI groups. Variables broadcasted are listed in 
  //  EinsplineSetBuilderCommon.cpp EinsplineSetBuilder::BroadcastOrbitalInfo()
  //   

  Timer orb_info_timer;
  // The tiling can be set by a simple vector, (e.g. 2x2x2), or by a
  // full 3x3 matrix of integers.  If the tilematrix was not set in
  // the input file...
  bool matrixNotSet = true;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      matrixNotSet = matrixNotSet && (TileMatrix(i,j) == 0);
  // then set the matrix to what may have been specified in the
  // tiling vector
  if (matrixNotSet)
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        TileMatrix(i,j) = (i==j) ? TileFactor[i] : 0;
  char buff[1000];
  if (myComm->rank() == 0)
  {
    snprintf (buff, 1000, "  TileMatrix = \n [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d ]\n",
             TileMatrix(0,0), TileMatrix(0,1), TileMatrix(0,2),
             TileMatrix(1,0), TileMatrix(1,1), TileMatrix(1,2),
             TileMatrix(2,0), TileMatrix(2,1), TileMatrix(2,2));
    app_log() << buff;
  }  
  if (numOrbs == 0)
  {
    app_error() << "You must specify the number of orbitals in the input file.\n";
    APP_ABORT("EinsplineSetBuilder::createSPOSet");
  }
  else
    app_log() << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";
  orb_info_timer.restart();
  /////////////////////////////////////////////////////////////////
  // Read the basic orbital information, without reading all the //
  // orbitals themselves.                                        //
  /////////////////////////////////////////////////////////////////
  if (myComm->rank() == 0)
    if (!ReadOrbitalInfo())
    {
      app_error() << "Error reading orbital info from HDF5 file.  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::createSPOSet");
    }
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  myComm->barrier();
  orb_info_timer.restart();
  BroadcastOrbitalInfo();

  app_log() <<  "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  app_log().flush();

  // setup primitive cell and supercell
  PrimCell.set(Lattice);
  SuperCell.set(SuperLattice);
  GGt=dot(transpose(PrimCell.G), PrimCell.G);
  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
    AtomicOrbitals[iat].Lattice = Lattice;

  // Now, analyze the k-point mesh to figure out the what k-points  are needed
  TwistNum = TwistNum_inp;
  AnalyzeTwists2();
}

SPOSetBase*
EinsplineSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  update_token(__FILE__,__LINE__,"createSPOSetFromXML");
  //use 2 bohr as the default when truncated orbitals are used based on the extend of the ions
  BufferLayer=2.0;
  SPOSetBase *OrbitalSet;
  int numOrbs = 0;
  qafm=0;
  int sortBands(1);
  int spinSet = 0;
  int TwistNum_inp=0;

  std::string sourceName;
  std::string spo_prec("double");
  std::string truncate("no");
  std::string hybrid_rep("no");
  std::string use_einspline_set_extended("no"); // use old spline library for high-order derivatives, e.g. needed for backflow optimization
#if defined(QMC_CUDA)
  std::string useGPU="yes";
#else
  std::string useGPU="no";
#endif
  NewTimer* spo_timer = new NewTimer("einspline::CreateSPOSetFromXML", timer_level_medium);
  TimerManager.addTimer(spo_timer);
  spo_timer->start();

  {
    OhmmsAttributeSet a;
    a.add (H5FileName, "href");
    a.add (TileFactor, "tile");
    a.add (sortBands,  "sort");
    a.add (qafm,  "afmshift");
    a.add (TileMatrix, "tilematrix");
    a.add (TwistNum_inp,   "twistnum");
    a.add (givenTwist,   "twist");
    a.add (sourceName, "source");
    a.add (MeshFactor, "meshfactor");
    a.add (hybrid_rep, "hybridrep");
    a.add (useGPU,     "gpu");
    a.add (spo_prec,   "precision");
    a.add (truncate,   "truncate");
    a.add (BufferLayer, "buffer");
    a.add (use_einspline_set_extended,"use_old_spline");
    a.add (myName, "tag");
#if defined(QMC_CUDA)
    a.add (gpu::MaxGPUSpineSizeMB, "Spline_Size_Limit_MB");
#endif

    a.put (XMLRoot);
    a.add (numOrbs,    "size");
    a.add (numOrbs,    "norbs");
    a.add(spinSet,"spindataset"); a.add(spinSet,"group");
    a.put (cur);

    if(myName.empty()) myName="einspline";

  }

  SourcePtcl=ParticleSets[sourceName];
  if(SourcePtcl==0)
  {
    APP_ABORT("Einspline needs the source particleset");
  }
  else
  { //keep the one-body distance table index 
#if defined(ENABLE_SOA)
    myTableIndex=TargetPtcl.addTable(*SourcePtcl,DT_SOA_PREFERRED);
#else
    myTableIndex=TargetPtcl.addTable(*SourcePtcl,DT_AOS);
#endif
  }

  ///////////////////////////////////////////////
  // Read occupation information from XML file //
  ///////////////////////////////////////////////
  std::vector<int> Occ_Old(0,0);
  Occ.resize(0,0);
  bool NewOcc(false);

  {
    OhmmsAttributeSet oAttrib;
    oAttrib.add(spinSet,"spindataset");
    oAttrib.put(cur);
  }

  xmlNodePtr spo_cur=cur;
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      std::string occ_mode("ground");
      occ_format="energy";
      particle_hole_pairs=0;
      OhmmsAttributeSet oAttrib;
      oAttrib.add(occ_mode,"mode");
      oAttrib.add(spinSet,"spindataset");
      oAttrib.add(occ_format,"format");
      oAttrib.add(particle_hole_pairs,"pairs");
      oAttrib.put(cur);
      if(occ_mode == "excited")
      {
        putContent(Occ,cur);
      }
      else if(occ_mode != "ground")
      {
        app_error() << "Only ground state occupation currently supported "
          << "in EinsplineSetBuilder.\n";
        APP_ABORT("EinsplineSetBuilder::createSPOSet");
      }
    }
    cur = cur->next;
  }
  if (Occ != Occ_Old)
  {
    NewOcc=true;
    Occ_Old = Occ;
  }
  else
    NewOcc=false;
#if defined(QMC_CUDA)
  app_log() << "\t  QMC_CUDA=1 Overwriting the einspline storage on the host to double precision.\n";
  spo_prec="double"; //overwrite
  truncate="no"; //overwrite
#endif
#if defined(MIXED_PRECISION)
  app_log() << "\t  MIXED_PRECISION=1 Overwriting the einspline storage to single precision.\n";
  spo_prec="single"; //overwrite
#endif
  H5OrbSet aset(H5FileName, spinSet, numOrbs);
  std::map<H5OrbSet,SPOSetBase*,H5OrbSet>::iterator iter;
  iter = SPOSetMap.find (aset);
  if ((iter != SPOSetMap.end() ) && (!NewOcc) && (qafm==0))
  {
    qafm=0;
    app_log() << "SPOSet parameters match in EinsplineSetBuilder:  "
              << "cloning EinsplineSet object.\n";
    return iter->second->makeClone();
  }

  if(FullBands[spinSet]==0) FullBands[spinSet]=new std::vector<BandInfo>;

  // Ensure the first SPO set must be spinSet==0
  // to correctly initialize key data of EinsplineSetBuilder
  if ( SPOSetMap.size()==0 && spinSet!=0 )
  {
    app_error() << "The first SPO set must have spindataset=\"0\"" << std::endl;
    abort();
  }

  // set the internal parameters
  if (spinSet == 0) set_metadata(numOrbs,TwistNum_inp);
  //if (use_complex_orb == "yes") UseRealOrbitals = false; // override given user input

  // look for <backflow>, would be a lot easier with xpath, but I cannot get it to work
  bool has_backflow = false;

  xmlNodePtr wf  = XMLRoot->parent; // <wavefuntion>
  xmlNodePtr kid = wf->children;
  while (kid != NULL)
  {
    std::string tag((const char*)(kid->name));
    if (tag=="determinantset" || tag=="sposet_builder")
    {
      xmlNodePtr kid1 = kid->children;
      while (kid1 != NULL)
      {
        std::string tag1((const char*)(kid1->name));
        if (tag1=="backflow")
        {
          has_backflow = true;
        }
        kid1 = kid1->next;
      } 
    }
    kid = kid->next; 
  }

  if (has_backflow && UseRealOrbitals) APP_ABORT("backflow optimization is broken with UseRealOrbitals");
  if (has_backflow && use_einspline_set_extended!="yes") APP_ABORT("backflow optimization does not yet function with EinsplinAdoptor, please add use_old_spline=\"yes\" to <determinantset> or <sposet_builder>.");

  //////////////////////////////////
  // Create the OrbitalSet object
  //////////////////////////////////
  Timer mytimer;
  mytimer.restart();
  OccupyBands(spinSet, sortBands, numOrbs);
  if(spinSet==0) TileIons();

  bool use_single= (spo_prec == "single" || spo_prec == "float");

#if !defined(QMC_COMPLEX)
  if (UseRealOrbitals)
  {
    //if(TargetPtcl.Lattice.SuperCellEnum != SUPERCELL_BULK && truncate=="yes")
    if(MixedSplineReader==0)
    {
      if(truncate=="yes")
      {
        if(use_single)
          MixedSplineReader=createTruncatedSingle(this,TargetPtcl.Lattice.SuperCellEnum);
        else
          MixedSplineReader=createTruncatedDouble(this,TargetPtcl.Lattice.SuperCellEnum);
      }
      else
      {
        if(use_single)
          MixedSplineReader= createBsplineRealSingle(this, hybrid_rep=="yes");
        else
          MixedSplineReader= createBsplineRealDouble(this, hybrid_rep=="yes");
      }
    }
  }
  else
#endif
  {
    if(MixedSplineReader==0)
    {
      if(truncate == "yes")
      {
        app_log() << "  Truncated orbitals with multiple kpoints are not supported yet!" << std::endl;
      }
      if(use_single)
        MixedSplineReader= createBsplineComplexSingle(this, hybrid_rep=="yes");
      else
        MixedSplineReader= createBsplineComplexDouble(this, hybrid_rep=="yes");
    }
  }

  MixedSplineReader->setCommon(XMLRoot);
  size_t delta_mem=qmc_common.memory_allocated;
  // temporary disable the following function call, Ye Luo
  // RotateBands_ESHDF(spinSet, dynamic_cast<EinsplineSetExtended<std::complex<double> >*>(OrbitalSet));
  HasCoreOrbs=bcastSortBands(spinSet,NumDistinctOrbitals,myComm->rank()==0);
  SPOSetBase* bspline_zd=MixedSplineReader->create_spline_set(spinSet,spo_cur);
  if(!bspline_zd)
    APP_ABORT_TRACE(__FILE__,__LINE__,"Failed to create SPOSetBase*");
  delta_mem=qmc_common.memory_allocated-delta_mem;
  app_log() <<"  MEMORY allocated SplineAdoptorReader " << (delta_mem>>20) << " MB" << std::endl;
  OrbitalSet = bspline_zd;
#if defined(MIXED_PRECISION)
  if(use_einspline_set_extended=="yes")
  {
    app_error() << "Option use_old_spline is not supported by the mixed precision build!" << std::endl;
    abort();
  }
#else
#ifndef QMC_CUDA
  if(use_einspline_set_extended=="yes")
#endif
  {
    EinsplineSet *new_OrbitalSet;
    if (UseRealOrbitals)
    {
      EinsplineSetExtended<double> *temp_OrbitalSet;
#if defined(QMC_CUDA)
      if (AtomicOrbitals.size() > 0)
        temp_OrbitalSet = new EinsplineSetHybrid<double>;
      else
#endif
        temp_OrbitalSet = new EinsplineSetExtended<double>;
      MixedSplineReader->export_MultiSpline(&(temp_OrbitalSet->MultiSpline));
      temp_OrbitalSet->MultiSpline->num_splines = NumDistinctOrbitals;
      temp_OrbitalSet->resizeStorage(NumDistinctOrbitals, NumValenceOrbs);
      //set the flags for anti periodic boundary conditions
      temp_OrbitalSet->HalfG = dynamic_cast<SplineAdoptorBase<double,3> *>(OrbitalSet)->HalfG;
      new_OrbitalSet = temp_OrbitalSet;
    }
    else
    {
      EinsplineSetExtended<std::complex<double> > *temp_OrbitalSet;
#if defined(QMC_CUDA)
      if (AtomicOrbitals.size() > 0)
        temp_OrbitalSet = new EinsplineSetHybrid<std::complex<double> >;
      else
#endif
        temp_OrbitalSet = new EinsplineSetExtended<std::complex<double> >;
      MixedSplineReader->export_MultiSpline(&(temp_OrbitalSet->MultiSpline));
      temp_OrbitalSet->MultiSpline->num_splines = NumDistinctOrbitals;
      temp_OrbitalSet->resizeStorage(NumDistinctOrbitals, NumValenceOrbs);
      for (int iorb=0, num=0; iorb<NumDistinctOrbitals; iorb++)
      {
        int ti = (*FullBands[spinSet])[iorb].TwistIndex;
        temp_OrbitalSet->kPoints[iorb] = PrimCell.k_cart(TwistAngles[ti]);
        temp_OrbitalSet->MakeTwoCopies[iorb] = (num < (numOrbs-1)) && (*FullBands[spinSet])[iorb].MakeTwoCopies;
        num += temp_OrbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
      }
      new_OrbitalSet = temp_OrbitalSet;
    }
    //set the internal parameters
    setTiling(new_OrbitalSet,numOrbs);
    OrbitalSet = new_OrbitalSet;
  }
#endif
#ifdef Ye_debug
#ifndef QMC_COMPLEX
  if (myComm->rank()==0 && OrbitalSet->MuffinTins.size() > 0)
  {
    FILE *fout  = fopen ("TestMuffins.dat", "w");
    Vector<double> phi(numOrbs), lapl(numOrbs);
    Vector<PosType> grad(numOrbs);
    ParticleSet P;
    P.R.resize(6);
    for (int i=0; i<P.R.size(); i++)
      P.R[i] = PosType (0.0, 0.0, 0.0);
    PosType N = 0.25*PrimCell.a(0) + 0.25*PrimCell.a(1) + 0.25*PrimCell.a(2);
    for (double x=-1.0; x<=1.0; x+=0.0000500113412)
    {
      // for (double x=-0.003; x<=0.003; x+=0.0000011329343481381) {
      P.R[0] = x * (PrimCell.a(0) + 0.914*PrimCell.a(1) +
                    0.781413*PrimCell.a(2));
      double r = std::sqrt(dot(P.R[0], P.R[0]));
      double rN = std::sqrt(dot(P.R[0]-N, P.R[0]-N));
      OrbitalSet->evaluate(P, 0, phi, grad, lapl);
      // OrbitalSet->evaluate(P, 0, phi);
      fprintf (fout, "%1.12e ", r*x/std::abs(x));
      for (int j=0; j<numOrbs; j++)
      {
        double gmag = std::sqrt(dot(grad[j],grad[j]));
        fprintf (fout, "%16.12e ",
                 /*phi[j]*phi[j]**/(-5.0/r  -0.5*lapl[j]/phi[j]));
        // double E = -5.0/r -0.5*lapl[j]/phi[j];
        fprintf (fout, "%16.12e ", phi[j]);
        fprintf (fout, "%16.12e ", gmag);
      }
      fprintf (fout, "\n");
    }
    fclose(fout);
  }
#endif
#endif
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadBands " << mytimer.elapsed() << std::endl;
  SPOSetMap[aset] = OrbitalSet;
  //if (sourceName.size() && (ParticleSets.find(sourceName) == ParticleSets.end()))
  //{
  //  app_log() << "  EinsplineSetBuilder creates a ParticleSet " << sourceName << std::endl;
  //  ParticleSet* ions=new ParticleSet;
  //  ions->Lattice=TargetPtcl.Lattice;
  //  ESHDFIonsParser ap(*ions,H5FileID,myComm);
  //  ap.put(XMLRoot);
  //  ap.expand(TileMatrix);
  //  ions->setName(sourceName);
  //  ParticleSets[sourceName]=ions;
  //  //overwrite the lattice and assign random
  //  if(TargetPtcl.Lattice.SuperCellEnum)
  //  {
  //    TargetPtcl.Lattice=ions->Lattice;
  //    makeUniformRandom(TargetPtcl.R);
  //    TargetPtcl.R.setUnit(PosUnit::LatticeUnit);
  //    TargetPtcl.convert2Cart(TargetPtcl.R);
  //    TargetPtcl.createSK();
  //  }
  //}
#ifdef QMC_CUDA
  if (useGPU == "yes" || useGPU == "1")
  {
    app_log() << "Initializing GPU data structures.\n";
    OrbitalSet->initGPU();
  }
#endif
  spo_timer->stop();
  return OrbitalSet;
}

SPOSetBase* EinsplineSetBuilder::createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input_info)
{
  update_token(__FILE__,__LINE__,"createSPOSet(cur,input_info)");

  if(MixedSplineReader==0)
  {
    APP_ABORT_TRACE(__FILE__,__LINE__,"EinsplineSetExtended<T> cannot create a SPOSet");
  }

  std::string aname;
  int spinSet(0);
  OhmmsAttributeSet a;
  a.add(aname,"name");
  a.add(spinSet,"spindataset");
  a.add(spinSet,"group");
  a.put(cur);

  if(aname.empty()) 
  {
    APP_ABORT_TRACE(__FILE__,__LINE__,"Missing sposet@name");
  }

  //allow only non-overlapping index sets and use the max index as the identifier
  int norb=input_info.max_index();
  H5OrbSet aset(H5FileName, spinSet, norb);

  SPOSetBase* bspline_zd=MixedSplineReader->create_spline_set(spinSet,cur,input_info);
  //APP_ABORT_TRACE(__FILE__,__LINE__,"DONE");
  if(bspline_zd)
    SPOSetMap[aset] = bspline_zd;
  return bspline_zd;
}

}


