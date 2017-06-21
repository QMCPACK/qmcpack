//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Utilities/ProgressReportEngine.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "ParticleIO/ESHDFParticleParser.h"
#include "ParticleBase/RandomSeqGenerator.h"

#define EINSPLINE_HELP_DEBUG
#if defined(EINSPLINE_HELP_DEBUG)
#include <QMCWaveFunctions/einspline_helper.hpp>
#endif

namespace qmcplusplus
{

#ifdef QMC_CUDA
inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_d *in,
    multi_UBspline_3d_s_cuda* &out)
{
  out = create_multi_UBspline_3d_s_cuda_conv (in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_d *in,
    multi_UBspline_3d_d_cuda * &out)
{
  out = create_multi_UBspline_3d_d_cuda(in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_c_cuda* &out)
{
  out = create_multi_UBspline_3d_c_cuda_conv (in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_z_cuda * &out)
{
  out = create_multi_UBspline_3d_z_cuda(in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_d_cuda * &out)
{
  app_error() << "Attempted to convert complex CPU spline into a real "
              << " GPU spline.\n";
  abort();
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_s_cuda * &out)
{
  app_error() << "Attempted to convert complex CPU spline into a real "
              << " GPU spline.\n";
  abort();
}
#endif

bool
EinsplineSetBuilder::ReadOrbitalInfo_ESHDF()
{
  app_log() << "  Reading orbital file in ESHDF format.\n";
//     TinyVector<int,3> version;
  HDFAttribIO<TinyVector<int,3> > h_version(Version);
  h_version.read (H5FileID, "/version");
  app_log() << "  ESHDF orbital file version "
            << Version[0] << "." << Version[1] << "." << Version[2] << std::endl;
  HDFAttribIO<Tensor<double,3> > h_Lattice(Lattice);
  h_Lattice.read      (H5FileID, "/supercell/primitive_vectors");
  RecipLattice = 2.0*M_PI*inverse(Lattice);
  SuperLattice = dot(TileMatrix, Lattice);
  char buff[1000];
  snprintf (buff, 1000,
            "  Lattice = \n    [ %9.6f %9.6f %9.6f\n"
            "      %9.6f %9.6f %9.6f\n"
            "      %9.6f %9.6f %9.6f ]\n",
            Lattice(0,0), Lattice(0,1), Lattice(0,2),
            Lattice(1,0), Lattice(1,1), Lattice(1,2),
            Lattice(2,0), Lattice(2,1), Lattice(2,2));
  app_log() << buff;
  snprintf (buff, 1000,
            "  SuperLattice = \n    [ %9.6f %9.6f %9.6f\n"
            "      %9.6f %9.6f %9.6f\n"
            "      %9.6f %9.6f %9.6f ]\n",
            SuperLattice(0,0), SuperLattice(0,1), SuperLattice(0,2),
            SuperLattice(1,0), SuperLattice(1,1), SuperLattice(1,2),
            SuperLattice(2,0), SuperLattice(2,1), SuperLattice(2,2));
  CheckLattice();
  app_log() << buff;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      LatticeInv(i,j) = RecipLattice(i,j)/(2.0*M_PI);
  int have_dpsi = 0;
  int NumAtomicOrbitals = 0;
  HDFAttribIO<int> h_NumBands(NumBands), h_NumElectrons(NumElectrons),
              h_NumSpins(NumSpins), h_NumTwists(NumTwists), h_NumCore(NumCoreStates),
              h_NumMuffinTins(NumMuffinTins), h_have_dpsi(have_dpsi),
              h_NumAtomicOrbitals(NumAtomicOrbitals);
  NumCoreStates = NumMuffinTins = NumTwists = NumSpins =
                                    NumBands = NumElectrons = NumAtomicOrbitals = 0;
  have_dpsi = false;
  h_NumBands.read      (H5FileID, "/electrons/kpoint_0/spin_0/number_of_states");
  h_NumCore.read       (H5FileID, "/electrons/kpoint_0/spin_0/number_of_core_states");
  h_NumElectrons.read  (H5FileID, "/electrons/number_of_electrons");
  h_NumSpins.read      (H5FileID, "/electrons/number_of_spins");
  h_NumTwists.read     (H5FileID, "/electrons/number_of_kpoints");
  h_NumMuffinTins.read (H5FileID, "/muffin_tins/number_of_tins");
  h_have_dpsi.read     (H5FileID, "/electrons/have_dpsi");
  h_NumAtomicOrbitals.read (H5FileID, "/electrons/number_of_atomic_orbitals");
  HaveOrbDerivs = have_dpsi;
  app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons
            << ", spins=" << NumSpins << ", twists=" << NumTwists
            << ", muffin tins=" << NumMuffinTins
            << ", core states=" << NumCoreStates << std::endl;
  app_log() << "atomic orbital=" << NumAtomicOrbitals << std::endl;
  if (TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1)
    app_log() << "  Using a " << TileFactor[0] << "x" << TileFactor[1]
              << "x" << TileFactor[2] << " tiling factor.\n";
  //////////////////////////////////
  // Read ion types and locations //
  //////////////////////////////////
  Vector<int> species_ids;
  HDFAttribIO<Vector<int> > h_species_ids(species_ids);
  h_species_ids.read (H5FileID, "/atoms/species_ids");
  int num_species;
  HDFAttribIO<int> h_num_species(num_species);
  h_num_species.read (H5FileID, "/atoms/number_of_species");
  std::vector<int> atomic_numbers(num_species);
  for (int isp=0; isp<num_species; isp++)
  {
    std::ostringstream name;
    name << "/atoms/species_" << isp << "/atomic_number";
    HDFAttribIO<int> h_atomic_number (atomic_numbers[isp]);
    h_atomic_number.read(H5FileID, name.str().c_str());
  }
  IonTypes.resize(species_ids.size());
  for (int i=0; i<species_ids.size(); i++)
    IonTypes[i] = atomic_numbers[species_ids[i]];
  HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
  h_IonPos.read   (H5FileID, "/atoms/positions");
  for (int i=0; i<IonTypes.size(); i++)
    app_log() << "Atom type(" << i << ") = " << IonTypes(i) << std::endl;
  /////////////////////////////////////
  // Read atomic orbital information //
  /////////////////////////////////////
  AtomicOrbitals.resize(NumAtomicOrbitals);
  for (int iat=0; iat<NumAtomicOrbitals; iat++)
  {
    AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
    int lmax, polynomial_order, spline_points;
    RealType cutoff_radius, polynomial_radius, spline_radius;
    PosType position;
    HDFAttribIO<int> h_lmax(lmax), h_polynomial_order(polynomial_order),
                h_spline_points(spline_points);
    HDFAttribIO<RealType> h_cutoff_radius(cutoff_radius),
                h_polynomial_radius(polynomial_radius),
                h_spline_radius(spline_radius);
    HDFAttribIO<PosType> h_position(position);
    std::ostringstream groupstream;
    groupstream << "/electrons/atomic_orbital_" << iat << "/";
    std::string groupname = groupstream.str();
    h_lmax.read              (H5FileID, (groupname + "lmax"             ).c_str());
    h_polynomial_order.read  (H5FileID, (groupname + "polynomial_order" ).c_str());
    h_spline_points.read     (H5FileID, (groupname + "spline_points"    ).c_str());
    h_cutoff_radius.read     (H5FileID, (groupname + "cutoff_radius"    ).c_str());
    h_polynomial_radius.read (H5FileID, (groupname + "polynomial_radius").c_str());
    h_spline_radius.read     (H5FileID, (groupname + "spline_radius"    ).c_str());
    h_position.read          (H5FileID, (groupname + "position"         ).c_str());
    orb.set_pos (position);
    orb.set_lmax (lmax);
    orb.set_cutoff (cutoff_radius);
    orb.set_spline (spline_radius, spline_points);
    orb.set_polynomial (polynomial_radius, polynomial_order);
  }
  ///////////////////////////
  // Read the twist angles //
  ///////////////////////////
  TwistAngles.resize(NumTwists);
  TwistSymmetry.resize(NumTwists);
  TwistWeight.resize(NumTwists);
  for (int ti=0; ti<NumTwists; ti++)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/reduced_k";
    HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
    h_Twist.read (H5FileID, path.str().c_str());
    std::ostringstream sym_path;
    sym_path << "/electrons/kpoint_" << ti << "/symgroup";
    HDFAttribIO<int> h_Sym(TwistSymmetry[ti]);
    h_Sym.read (H5FileID, sym_path.str().c_str());
    std::ostringstream nsym_path;
    nsym_path << "/electrons/kpoint_" << ti << "/numsym";
    HDFAttribIO<int> h_Nsym(TwistWeight[ti]);
    h_Nsym.read (H5FileID, nsym_path.str().c_str());
    app_log()<<"HERE I AM"<< TwistSymmetry[ti]<<"  "<<TwistWeight[ti] << std::endl;
    // Early versions from wfconv had wrong sign convention for
    // k-points.  EinsplineSet uses the opposite sign convention
    // from most DFT codes.
    if (version[0] >= 2)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        TwistAngles[ti][dim] *= -1.0;
    snprintf (buff, 1000, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n",
              TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    app_log() << buff;
  }
  //////////////////////////////////////////////////////////
  // If the density has not been set in TargetPtcl, and   //
  // the density is available, read it in and save it     //
  // in TargetPtcl.                                       //
  //////////////////////////////////////////////////////////
  // FIXME:  add support for more than one spin density
  if (!TargetPtcl.Density_G.size())
  {
    HDFAttribIO<std::vector<TinyVector<int,OHMMS_DIM> > >
    h_reduced_gvecs(TargetPtcl.DensityReducedGvecs);
    HDFAttribIO<Array<RealType,OHMMS_DIM> >
    h_density_r (TargetPtcl.Density_r);
    TinyVector<int,3> mesh;
    h_reduced_gvecs.read (H5FileID, "/electrons/density/gvectors");
    int numG = TargetPtcl.DensityReducedGvecs.size();
    // Convert primitive G-vectors to supercell G-vectors
    // Also, flip sign since ESHDF format uses opposite sign convention
    #pragma omp parallel for
    for (int iG=0; iG < numG; iG++)
      TargetPtcl.DensityReducedGvecs[iG] =
        -1 * dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
    app_log() << "  Read " << numG << " density G-vectors.\n";
    for (int ispin=0; ispin<NumSpins; ispin++)
    {
      std::ostringstream density_r_path, density_g_path;
      density_r_path << "/electrons/density/spin_" << ispin << "/density_r";
      density_g_path << "/electrons/density/spin_" << ispin << "/density_g";
      h_density_r.read (H5FileID, density_r_path.str().c_str());
      if (TargetPtcl.DensityReducedGvecs.size())
      {
        app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
        std::vector<ComplexType> density_G;
        HDFAttribIO<std::vector<ComplexType > > h_density_G (density_G);
        h_density_G.read (H5FileID, density_g_path.str().c_str());
        if (!density_G.size())
        {
          app_error() << "  Density reduced G-vectors defined, but not the"
                      << " density.\n";
          abort();
        }
        else
        {
          if (ispin == 0)
            TargetPtcl.Density_G = density_G;
          else
            for (int iG=0; iG<density_G.size(); iG++)
              TargetPtcl.Density_G[iG] += density_G[iG];
        }
      }
    }
  }
  //////////////////////////////////////////////////////////
  // If the density has not been set in TargetPtcl, and   //
  // the density is available, read it in and save it     //
  // in TargetPtcl.                                       //
  //////////////////////////////////////////////////////////
  // FIXME:  add support for more than one spin potential
  if (!TargetPtcl.VHXC_r[0].size())
  {
    HDFAttribIO<std::vector<TinyVector<int,OHMMS_DIM> > >
    h_reduced_gvecs(TargetPtcl.VHXCReducedGvecs);
    TinyVector<int,3> mesh;
    h_reduced_gvecs.read (H5FileID, "/electrons/VHXC/gvectors");
    int numG = TargetPtcl.VHXCReducedGvecs.size();
    // Convert primitive G-vectors to supercell G-vectors
    // Also, flip sign since ESHDF format uses opposite sign convention
    #pragma omp parallel for
    for (int iG=0; iG < numG; iG++)
      TargetPtcl.VHXCReducedGvecs[iG] =
        -1 * dot(TileMatrix, TargetPtcl.VHXCReducedGvecs[iG]);
    app_log() << "  Read " << numG << " VHXC G-vectors.\n";
    for (int ispin=0; ispin<NumSpins; ispin++)
    {
      HDFAttribIO<Array<RealType,OHMMS_DIM> >
      h_VHXC_r (TargetPtcl.VHXC_r[ispin]);
      std::ostringstream VHXC_r_path, VHXC_g_path;
      VHXC_r_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_r";
      VHXC_g_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_g";
      h_VHXC_r.read (H5FileID, VHXC_r_path.str().c_str());
      if (TargetPtcl.VHXCReducedGvecs.size())
      {
        app_log() << "  EinsplineSetBuilder found VHXC in the HDF5 file.\n";
        std::vector<ComplexType> VHXC_G;
        HDFAttribIO<std::vector<ComplexType > > h_VHXC_G (VHXC_G);
        h_VHXC_G.read (H5FileID, VHXC_g_path.str().c_str());
        if (!VHXC_G.size())
        {
          app_error() << "  VHXC reduced G-vectors defined, but not the"
                      << " VHXC.\n";
          abort();
        }
        else
          TargetPtcl.VHXC_G[ispin] = VHXC_G;
      }
    }
  }
  HaveLocalizedOrbs = false;
  return true;
}


void
EinsplineSetBuilder::OccupyBands_ESHDF(int spin, int sortBands)
{
  //testing
  if (myComm->rank() != 0)
    return;
  SortBands.clear();
  int maxOrbs(0);
  for (int ti=0; ti<DistinctTwists.size(); ti++)
  {
    int tindex = DistinctTwists[ti];
    // First, read valence states
    std::ostringstream ePath;
    ePath << "/electrons/kpoint_" << tindex << "/spin_"
          << spin << "/eigenvalues";
    std::vector<double> eigvals;
    HDFAttribIO<std::vector<double> > h_eigvals(eigvals);
    h_eigvals.read(H5FileID, ePath.str().c_str());
    for (int bi=0; bi<NumBands; bi++)
    {
      BandInfo band;
      band.IsCoreState = false;
      band.TwistIndex = tindex;
      band.BandIndex  = bi;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      band.Energy = eigvals[bi];
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
      if (MakeTwoCopies[ti])
        maxOrbs+=2;
      else
        maxOrbs++;
    }
    // Now, read core states
    for (int cs=0; cs<NumCoreStates; cs++)
    {
      BandInfo band;
      band.IsCoreState = true;
      band.TwistIndex = tindex;
      band.BandIndex  = cs;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      HDFAttribIO<double> h_energy(band.Energy);
      h_energy.read   (H5FileID, (CoreStatePath(ti,cs)+"eigenvalue").c_str());
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
      if (MakeTwoCopies[ti])
        maxOrbs+=2;
      else
        maxOrbs++;
    }
  }
  // Now sort the bands by energy
  if (sortBands==2)
  {
    app_log() << "Sorting the bands by index now:\n";
    sort (SortBands.begin(), SortBands.end(), sortByIndex);
  }
  else
    if (sortBands==1)
    {
      app_log() << "Sorting the bands now:\n";
      sort (SortBands.begin(), SortBands.end());
    }
  //occupy the ground state first
  std::vector<int> gsOcc(maxOrbs);
  int N_gs_orbs=OrbitalSet->getOrbitalSetSize();
  int nocced(0);
  for (int ti=0; ti<SortBands.size(); ti++)
  {
    if (nocced<N_gs_orbs)
    {
      if (SortBands[ti].MakeTwoCopies && (N_gs_orbs-nocced>1))
      {
        nocced+=2;
        gsOcc[ti]=2;
      }
      else
        if ( (SortBands[ti].MakeTwoCopies && (N_gs_orbs-nocced==1)) || !SortBands[ti].MakeTwoCopies )
        {
          nocced+=1;
          gsOcc[ti]=1;
        }
    }
  }
  if (occ_format=="energy")
  {
    // To get the occupations right.
    std::vector<int> Removed(0,0);
    std::vector<int> Added(0,0);
    for(int ien=0; ien<Occ.size(); ien++)
    {
      if (Occ[ien]<0)
        Removed.push_back(-Occ[ien]);
      else
        if (Occ[ien]>0)
          Added.push_back(Occ[ien]);
    }
    if(Added.size()-Removed.size() != 0)
    {
      app_log()<<"need to add and remove same number of orbitals. "<< Added.size()<<" "<<Removed.size()<< std::endl;
      APP_ABORT("ChangedOccupations");
    }
    std::vector<int> DiffOcc(maxOrbs,0);
    //Probably a cleaner way to do this.
    for(int i=0; i<Removed.size(); i++)
      DiffOcc[Removed[i]-1]-=1;
    for(int i=0; i<Added.size(); i++)
      DiffOcc[Added[i]-1]+=1;
    std::vector<int> SumOrb(SortBands.size(),0);
    int doi(0);
    for(int i=0; i<SumOrb.size(); i++)
    {
      if(SortBands[i].MakeTwoCopies)
      {
        SumOrb[i]=  gsOcc[i]+DiffOcc[doi++];
        SumOrb[i]+= DiffOcc[doi++];
      }
      else
        SumOrb[i]=gsOcc[i]+DiffOcc[doi++];
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for(int i=0; i<SumOrb.size(); i++)
    {
      if(SumOrb[i]==2)
      {
        SortBands[i].MakeTwoCopies=true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else
        if (SumOrb[i]==1)
        {
          SortBands[i].MakeTwoCopies=false;
          ReOrderedBands.push_back(SortBands[i]);
        }
        else
          if (SumOrb[i]==0)
          {
            SortBands[i].MakeTwoCopies=false;
            RejectedBands.push_back(SortBands[i]);
          }
          else
          {
            app_log()<<" Trying to add the same orbital ("<<i<<") less than zero or more than 2 times."<< std::endl;
            APP_ABORT("Sorting Excitation");
          }
    }
    ReOrderedBands.insert(ReOrderedBands.end(),RejectedBands.begin(),RejectedBands.end());
    SortBands=ReOrderedBands;
  }
  else
    if (occ_format=="band")
    {
      app_log()<<"  Occupying bands based on (bi,ti) data."<< std::endl;
      if(Occ.size() != particle_hole_pairs*4)
      {
        app_log()<<" Need Occ = pairs*4. Occ is (ti,bi) of removed, then added."<< std::endl;
        app_log()<<Occ.size()<<" "<<particle_hole_pairs<< std::endl;
        APP_ABORT("ChangedOccupations");
      }
      int cnt(0);
      for(int ien=0; ien<SortBands.size(); ien++)
      {
        if((Occ[cnt] == SortBands[ien].TwistIndex)&&(Occ[cnt+1] == SortBands[ien].BandIndex))
          if(cnt<particle_hole_pairs*2)
          {
            gsOcc[ien]-=1;
            cnt+=2;
            app_log()<<"removing orbital "<<ien<< std::endl;
          }
          else
          {
            gsOcc[ien]+=1;
            app_log()<<"adding orbital "<<ien<< std::endl;
            cnt+=2;
          }
      }
      std::vector<BandInfo> ReOrderedBands;
      std::vector<BandInfo> RejectedBands;
      for(int i=0; i<SortBands.size(); i++)
      {
        if(gsOcc[i]==2)
        {
          SortBands[i].MakeTwoCopies=true;
          ReOrderedBands.push_back(SortBands[i]);
        }
        else
          if (gsOcc[i]==1)
          {
            SortBands[i].MakeTwoCopies=false;
            ReOrderedBands.push_back(SortBands[i]);
          }
          else
            if (gsOcc[i]==0)
            {
              SortBands[i].MakeTwoCopies=false;
              RejectedBands.push_back(SortBands[i]);
            }
            else
            {
              app_log()<<" Trying to add the same orbital ("<<i<<") less than zero or more than 2 times."<< std::endl;
              APP_ABORT("Sorting Excitation");
            }
      }
      ReOrderedBands.insert(ReOrderedBands.end(),RejectedBands.begin(),RejectedBands.end());
      SortBands=ReOrderedBands;
    }
  //for(int sw=0;sw<Removed.size();sw++){
  //  app_log()<<" Swapping two orbitals "<<Removed[sw]<<" and "<<Added[sw]<< std::endl;
  //  BandInfo tempband(SortBands[Removed[sw]-1]);
  //  SortBands[Removed[sw]-1] = SortBands[Added[sw]-1];
  //  SortBands[Added[sw]-1] = tempband;
  //}
  int orbIndex = 0;
  int numOrbs = 0;
  NumValenceOrbs=0;
  NumCoreOrbs=0;
  while (numOrbs < OrbitalSet->getOrbitalSetSize())
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs += 2;
    else
      numOrbs++;
    if (SortBands[orbIndex].IsCoreState)
      NumCoreOrbs++;
    else
      NumValenceOrbs++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct orbitals.\n";
  app_log() << "There are " << NumCoreOrbs << " core states and "
            << NumValenceOrbs << " valence states.\n";
}


void
EinsplineSetBuilder::RotateBands_ESHDF
(int spin, EinsplineSetExtended<std::complex<double > >* orbitalSet)
{
  bool root = (myComm->rank()==0);
  if (root)
  {
    rotationMatrix.resize(0);
    rotatedOrbitals.resize(0);
    xmlNodePtr kids=XMLRoot->children;
    while(kids != NULL)
    {
      std::string cname((const char*)(kids->name));
      if(cname == "rotationmatrix")
        putContent(rotationMatrix,kids);
      else
        if(cname=="rotatedorbitals")
          putContent(rotatedOrbitals,kids);
      kids=kids->next;
    }
    if ((rotatedOrbitals.size()*rotatedOrbitals.size() != rotationMatrix.size()) && (rotationMatrix.size()!=0))
    {
      app_log()<<" Rotation Matrix is wrong dimension. "<<rotationMatrix.size()<<" should be "<<rotatedOrbitals.size()*rotatedOrbitals.size()<< std::endl;
    }
    else
      if (rotationMatrix.size()>0)
      {
        app_log()<<" Rotating between: ";
        for (int i=0; i<rotatedOrbitals.size(); i++)
          app_log()<<rotatedOrbitals[i]<<" ";
        app_log()<< std::endl;
        app_log()<<" Using the following rotation"<< std::endl;
        for (int i=0; i<rotatedOrbitals.size(); i++)
        {
          for (int j=0; j<rotatedOrbitals.size(); j++)
            app_log()<<rotationMatrix[rotatedOrbitals.size()*i+j]<<" ";
          app_log()<< std::endl;
        }
      }
    if ((rotationMatrix.size()>0) && (rotatedOrbitals.size()>0) )
    {
      int N = NumDistinctOrbitals;
      int num(0);
      for (int iorb=0, indx=0; iorb<N; iorb++)
      {
        num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
        if (num==rotatedOrbitals[indx])
        {
          rotatedOrbitals[indx]=iorb;
          indx++;
        }
      }
      //simple copy file function. make backup.
      std::string backupName = H5FileName+"_bkup";
      std::ifstream fin(H5FileName.c_str(), std::ios::in | std::ios::binary);
      std::ofstream fout(backupName.c_str() , std::ios::in); // open with this mode to check whether file exists
      //       std::ofstream fout(backupName.c_str(), std::ios::out | std::ios::binary);
      if (fin.fail())
      {
        // reset status flags
        fin.clear();
        std::cout << " source file does not exist, try it again"<< std::endl;
        exit( 0 );
      }
      if (!fout.fail())
      {
        fout.close();
        std::cout << " destination file already exists, backup completed"<< std::endl;
      }
      else
      {
        fout.close();
        fout.open(backupName.c_str() , std::ios::out | std::ios::binary); // change to writting mode
        int BUFFER_SIZE = 128;
        char buffer[BUFFER_SIZE];
        while (!fin.eof() )
        {
          fin.read( buffer, BUFFER_SIZE);
          if (fin.bad())
          {
            std::cout << "Error reading data" << std::endl;
            exit( 0 );
          }
          else
            fout.write(buffer, fin.gcount());
        }
      }
      fin.close();
      fout.close();
      int nx, ny, nz, bi, ti;
      std::vector<Array<std::complex<double>,3> > allRotatedSplines;
      Array<std::complex<double>,3> splineData;
      TinyVector<int,3> mesh;
      // Find the orbital mesh size
      HDFAttribIO<TinyVector<int,3> > h_mesh(mesh);
      h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
      h_mesh.read (H5FileID, "/electrons/mesh");
      //     myComm->bcast(mesh);
      nx=mesh[0];
      ny=mesh[1];
      nz=mesh[2];
      splineData.resize(nx,ny,nz);
      for (int i=0; i<rotatedOrbitals.size(); i++)
      {
        int iorb = rotatedOrbitals[i];
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        double e = SortBands[iorb].Energy;
        PosType k;
        PosType twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        fprintf (stderr, "  Rotating state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d \n",
                 ti, bi, e, k[0], k[1], k[2], myComm->rank() );
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        std::string psirName = path.str() + "psi_r";
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(splineData);
        h_splineData.read(H5FileID, psirName.c_str());
        if ((splineData.size(0) != nx) ||
            (splineData.size(1) != ny) ||
            (splineData.size(2) != nz))
        {
          fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
          fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
          abort();
        }
        allRotatedSplines.push_back(splineData);
      }
      app_log()<< std::endl;
      std::vector<Array<std::complex<double>,3> > allOriginalSplines(allRotatedSplines);
      #pragma omp parallel for
      for (int i=0; i<rotatedOrbitals.size(); i++)
        for (int ix=0; ix<nx; ix++)
          for (int iy=0; iy<ny; iy++)
            for (int iz=0; iz<nz; iz++)
              allRotatedSplines[i](ix,iy,iz)=0.0;
      #pragma omp parallel for
      for (int i=0; i<rotatedOrbitals.size(); i++)
      {
        for(int j=0; j<rotatedOrbitals.size(); j++)
        {
          for (int ix=0; ix<nx; ix++)
            for (int iy=0; iy<ny; iy++)
              for (int iz=0; iz<nz; iz++)
                allRotatedSplines[i](ix,iy,iz) += rotationMatrix[i*rotatedOrbitals.size()+j] * allOriginalSplines[j](ix,iy,iz);
        }
      }
      for (int i=0; i<rotatedOrbitals.size(); i++)
      {
        int iorb = rotatedOrbitals[i];
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        std::string psirName = path.str() + "psi_r";
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(allRotatedSplines[i],true);
        h_splineData.write(H5FileID, psirName.c_str());
      }
//      for (int i=0;i<rotatedOrbitals.size();i++){
//           int iorb = rotatedOrbitals[i];
//           int ti   = SortBands[iorb].TwistIndex;
//           int bi   = SortBands[iorb].BandIndex;
//
//           std::ostringstream path;
//           path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
//           std::string psirName = path.str() + "psi_r";
//
//           HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(allOriginalSplines[i]);
//           h_splineData.read(H5FileID, psirName.c_str());
//         }
    }
    else
      app_log()<<" No rotations defined"<< std::endl;
  }
}

void
EinsplineSetBuilder::RotateBands_ESHDF
(int spin, EinsplineSetExtended<double>* orbitalSet)
{
  bool root = (myComm->rank()==0);
  if (root)
  {
    rotationMatrix.resize(0);
    rotatedOrbitals.resize(0);
    xmlNodePtr kids=XMLRoot->children;
    while(kids != NULL)
    {
      std::string cname((const char*)(kids->name));
      if(cname == "rotationmatrix")
        putContent(rotationMatrix,kids);
      else
        if(cname=="rotatedorbitals")
          putContent(rotatedOrbitals,kids);
      kids=kids->next;
    }
    if ((rotatedOrbitals.size()*rotatedOrbitals.size() != rotationMatrix.size()) && (rotationMatrix.size()!=0))
    {
      app_log()<<" Rotation Matrix is wrong dimension. "<<rotationMatrix.size()<<" should be "<<rotatedOrbitals.size()*rotatedOrbitals.size()<< std::endl;
    }
    else
      if (rotationMatrix.size()>0)
      {
        app_log()<<" Rotating between: ";
        for (int i=0; i<rotatedOrbitals.size(); i++)
          app_log()<<rotatedOrbitals[i]<<" ";
        app_log()<< std::endl;
        app_log()<<" Using the following rotation"<< std::endl;
        for (int i=0; i<rotatedOrbitals.size(); i++)
        {
          for (int j=0; j<rotatedOrbitals.size(); j++)
            app_log()<<rotationMatrix[rotatedOrbitals.size()*i+j]<<" ";
          app_log()<< std::endl;
        }
      }
    if ((rotationMatrix.size()>0) && (rotatedOrbitals.size()>0) )
    {
      int N = NumDistinctOrbitals;
      int num(0);
      for (int iorb=0, indx=0; iorb<N; iorb++)
      {
        num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
        if (num==rotatedOrbitals[indx])
        {
          rotatedOrbitals[indx]=iorb;
          indx++;
        }
      }
      //simple copy file function. make backup.
      std::string backupName = H5FileName+"_bkup";
      std::ifstream fin(H5FileName.c_str(), std::ios::in | std::ios::binary);
      std::ofstream fout(backupName.c_str() , std::ios::in); // open with this mode to check whether file exists
      //       std::ofstream fout(backupName.c_str(), std::ios::out | std::ios::binary);
      if (fin.fail())
      {
        // reset status flags
        fin.clear();
        std::cout << " source file does not exist, try it again"<< std::endl;
        exit( 0 );
      }
      if (!fout.fail())
      {
        fout.close();
        std::cout << " destination file already exists, backup completed"<< std::endl;
      }
      else
      {
        fout.close();
        fout.open(backupName.c_str() , std::ios::out | std::ios::binary); // change to writting mode
        int BUFFER_SIZE = 128;
        char buffer[BUFFER_SIZE];
        while (!fin.eof() )
        {
          fin.read( buffer, BUFFER_SIZE);
          if (fin.bad())
          {
            std::cout << "Error reading data" << std::endl;
            exit( 0 );
          }
          else
            fout.write(buffer, fin.gcount());
        }
      }
      fin.close();
      fout.close();
      int nx, ny, nz, bi, ti;
      std::vector<Array<std::complex<double>,3> > allRotatedSplines;
      Array<std::complex<double>,3> splineData;
      TinyVector<int,3> mesh;
      // Find the orbital mesh size
      HDFAttribIO<TinyVector<int,3> > h_mesh(mesh);
      h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
      h_mesh.read (H5FileID, "/electrons/mesh");
      //     myComm->bcast(mesh);
      nx=mesh[0];
      ny=mesh[1];
      nz=mesh[2];
      splineData.resize(nx,ny,nz);
      for (int i=0; i<rotatedOrbitals.size(); i++)
      {
        int iorb = rotatedOrbitals[i];
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        double e = SortBands[iorb].Energy;
        PosType k;
        PosType twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        fprintf (stderr, "  Rotating state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d \n",
                 ti, bi, e, k[0], k[1], k[2], myComm->rank() );
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        std::string psirName = path.str() + "psi_r";
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(splineData);
        h_splineData.read(H5FileID, psirName.c_str());
        if ((splineData.size(0) != nx) ||
            (splineData.size(1) != ny) ||
            (splineData.size(2) != nz))
        {
          fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
          fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
          abort();
        }
        allRotatedSplines.push_back(splineData);
      }
      app_log()<< std::endl;
      std::vector<Array<std::complex<double>,3> > allOriginalSplines(allRotatedSplines);
      #pragma omp parallel for
      for (int i=0; i<rotatedOrbitals.size(); i++)
        for (int ix=0; ix<nx; ix++)
          for (int iy=0; iy<ny; iy++)
            for (int iz=0; iz<nz; iz++)
              allRotatedSplines[i](ix,iy,iz)=0.0;
      #pragma omp parallel for
      for (int i=0; i<rotatedOrbitals.size(); i++)
      {
        for(int j=0; j<rotatedOrbitals.size(); j++)
        {
          for (int ix=0; ix<nx; ix++)
            for (int iy=0; iy<ny; iy++)
              for (int iz=0; iz<nz; iz++)
                allRotatedSplines[i](ix,iy,iz) += rotationMatrix[i*rotatedOrbitals.size()+j] * allOriginalSplines[j](ix,iy,iz);
        }
      }
      for (int i=0; i<rotatedOrbitals.size(); i++)
      {
        int iorb = rotatedOrbitals[i];
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        std::string psirName = path.str() + "psi_r";
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(allRotatedSplines[i],true);
        h_splineData.write(H5FileID, psirName.c_str());
      }
//      for (int i=0;i<rotatedOrbitals.size();i++){
//           int iorb = rotatedOrbitals[i];
//           int ti   = SortBands[iorb].TwistIndex;
//           int bi   = SortBands[iorb].BandIndex;
//
//           std::ostringstream path;
//           path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
//           std::string psirName = path.str() + "psi_r";
//
//           HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(allOriginalSplines[i]);
//           h_splineData.read(H5FileID, psirName.c_str());
//         }
    }
    else
      app_log()<<" No rotations defined"<< std::endl;
  }
}

void
EinsplineSetBuilder::ReadBands_ESHDF
(int spin, EinsplineSetExtended<std::complex<double > >* orbitalSet)
{
  ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF(EinsplineSetExtended<std::complex<double > >*");
  bool root = myComm->rank()==0;
  // bcast other stuff
  myComm->bcast (NumDistinctOrbitals);
  myComm->bcast (NumValenceOrbs);
  myComm->bcast (NumCoreOrbs);
  int N = NumDistinctOrbitals;
  orbitalSet->kPoints.resize(N);
  orbitalSet->MakeTwoCopies.resize(N);
  orbitalSet->StorageValueVector.resize(N);
  orbitalSet->BlendValueVector.resize(N);
  orbitalSet->StorageLaplVector.resize(N);
  orbitalSet->BlendLaplVector.resize(N);
  orbitalSet->StorageGradVector.resize(N);
  orbitalSet->BlendGradVector.resize(N);
  orbitalSet->StorageHessVector.resize(N);
  orbitalSet->StorageGradHessVector.resize(N);
  orbitalSet->phase.resize(N);
  orbitalSet->eikr.resize(N);
  orbitalSet->NumValenceOrbs = NumValenceOrbs;
  orbitalSet->NumCoreOrbs    = NumCoreOrbs;
  // Read in k-points
  int numOrbs = orbitalSet->getOrbitalSetSize();
  int num = 0;
  if (root)
  {
    for (int iorb=0; iorb<N; iorb++)
    {
      int ti = SortBands[iorb].TwistIndex;
      PosType twist  = TwistAngles[ti];
      orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
      orbitalSet->MakeTwoCopies[iorb] =
        (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
      num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
    }
  }
  myComm->bcast(orbitalSet->kPoints);
  myComm->bcast(orbitalSet->MakeTwoCopies);
  // First, check to see if we have already read this in
  H5OrbSet set(H5FileName, spin, N);
  bool havePsig=ReadGvectors_ESHDF();
  ///// @todo IMPROVE BELOW
  //MeshSize = 0; //TinyVector<int,3>(0,0,0);
  //// Find the orbital mesh size
  //if (root) {
  //  HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize);
  //  h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
  //  h_mesh.read (H5FileID, "/electrons/mesh");
  //}
  //// HACK HACK HACK
  ////bool havePsir = false;//MeshSize[0] != 0;
  //bool havePsir = MeshSize[0] != 0;
  //myComm->bcast(havePsir);
  //if(havePsir)
  //  myComm->bcast(MeshSize);
  //else
  //  ReadGvectors_ESHDF();
  app_log() << "MeshSize = (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  /// IMPROVE ABOVE
  int nx, ny, nz, bi, ti;
  nx=MeshSize[0];
  ny=MeshSize[1];
  nz=MeshSize[2];
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = PERIODIC;
  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;
  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;
  zBC.rCode = PERIODIC;
  x_grid.start = 0.0;
  x_grid.end = 1.0;
  x_grid.num = nx;
  y_grid.start = 0.0;
  y_grid.end = 1.0;
  y_grid.num = ny;
  z_grid.start = 0.0;
  z_grid.end = 1.0;
  z_grid.num = nz;
  // Create the multiUBspline object
  orbitalSet->MultiSpline =
    create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  //////////////////////////////////////
  // Create the MuffinTin APW splines //
  //////////////////////////////////////
  orbitalSet->MuffinTins.resize(NumMuffinTins);
  for (int tin=0; tin<NumMuffinTins; tin++)
  {
    orbitalSet->MuffinTins[tin].Atom = tin;
    orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
    orbitalSet->MuffinTins[tin].set_lattice(Lattice);
    orbitalSet->MuffinTins[tin].init_APW
    (MT_APW_rgrids[tin], MT_APW_lmax[tin],
     NumValenceOrbs);
  }
  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
  {
    AtomicOrbitals[iat].set_num_bands(NumValenceOrbs);
    AtomicOrbitals[iat].allocate();
  }
  int isComplex=1;
  if (root)
  {
    HDFAttribIO<int> h_isComplex(isComplex);
    h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");
  }
  myComm->bcast(isComplex);
  if (!isComplex)
  {
    APP_ABORT("Expected complex orbitals in ES-HDF file, but found real ones.");
  }
  EinsplineSetBuilder::RotateBands_ESHDF(spin, orbitalSet);
  bool isCore = bcastSortBands(N,root);
  if(isCore)
  {
    APP_ABORT("Core states not supported by ES-HDF yet.");
  }
  Array<std::complex<double>,3> splineData(nx,ny,nz);
  if(havePsig)
  {
    FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
    FFTplan = fftw_plan_dft_3d
              (MeshSize[0], MeshSize[1], MeshSize[2],
               reinterpret_cast<fftw_complex*>(FFTbox.data()),
               reinterpret_cast<fftw_complex*>(FFTbox.data()),
               +1, FFTW_ESTIMATE);
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      Vector<std::complex<double> > cG;
      int ncg=0;
      int ti=SortBands[iorb].TwistIndex;
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti    //SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
        HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
        h_cG.read (H5FileID, path.str().c_str());
        ncg=cG.size();
      }
      myComm->bcast(ncg);
      if(ncg != Gvecs[ti].size())
      {
        APP_ABORT("Failed : ncg != Gvecs[ti].size()");
      }
      if(!root)
        cG.resize(ncg);
      myComm->bcast(cG);
      unpack4fftw(cG,Gvecs[ti],MeshSize,FFTbox);
      fftw_execute (FFTplan);
      fix_phase_rotate(FFTbox,splineData,TwistAngles[ti]);
      set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
    }
    fftw_destroy_plan(FFTplan);
  }
  else
  {
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      //check dimension
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
        HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(splineData);
        h_splineData.read(H5FileID, path.str().c_str());
      }
      myComm->bcast(splineData);
      set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
    }
  }
//    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
//    {
//      int ti   = SortBands[iorb].TwistIndex;
//      int bi   = SortBands[iorb].BandIndex;
//      double e = SortBands[iorb].Energy;
//
//      PosType twist=TwistAngles[ti];
//      PosType k= orbitalSet->PrimLattice.k_cart(twist);
//
//      if(root)
//      {
//        fprintf (stderr, "  Valence state:  ti=%3d  bi=%3d energy=%8.5f k=(%7.4f, %7.4f, %7.4f) rank=%d\n",
//            ti, bi, e, k[0], k[1], k[2], myComm->rank());
//        std::ostringstream path;
//        path << "/electrons/kpoint_" << ti << "/spin_" << spin
//          << "/state_" << bi << "/";
//        if (havePsir) {
//          std::string psirName = path.str() + "psi_r";
//
//          HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(splineData);
//          h_splineData.read(H5FileID, psirName.c_str());
//          if ((splineData.size(0) != nx) ||
//              (splineData.size(1) != ny) ||
//              (splineData.size(2) != nz)) {
//            fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
//            fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
//            abort();
//          }
//        }
//        else {
//          std::string psiGname = path.str() + "psi_g";
//          // Array<std::complex<double>,1> cG;
//          // HDFAttribIO<Array<std::complex<double>,1> >  h_cG(cG);
//          Vector<std::complex<double> > cG;
//          HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
//          h_cG.read (H5FileID, psiGname.c_str());
//          assert (cG.size() == Gvecs[ti].size());
//#if defined(EINSPLINE_HELP_DEBUG)
//          unpack4fftw(cG,Gvecs[ti],MeshSize,FFTbox);
//#else
//          FFTbox = std::complex<double>();
//#pragma omp parallel for
//          for (int iG=0; iG<cG.size(); iG++) {
//            TinyVector<int,3> index = Gvecs[ti][iG];
//            index[0] = ((index[0] + MeshSize[0])%MeshSize[0]);
//            index[1] = ((index[1] + MeshSize[1])%MeshSize[1]);
//            index[2] = ((index[2] + MeshSize[2])%MeshSize[2]);
//            FFTbox(index[0], index[1], index[2]) = cG[iG];
//          }
//#endif
//          fftw_execute (FFTplan);
//          // Now, rotate the phase of the orbitals so that neither
//          // the real part nor the imaginary part are very near
//          // zero.  This sometimes happens in crystals with high
//          // symmetry at special k-points.
//
//#if defined(EINSPLINE_HELP_DEBUG)
//          fix_phase_rotate(FFTbox,splineData,TwistAngles[ti]);
//#else
//          double rNorm=0.0, iNorm=0.0;
//#pragma omp parallel for reduction(+:rNorm,iNorm)
//          for (int ix=0; ix<nx; ix++) {
//            PosType ru;
//            //	      ru[0] = (RealType)ix / (RealType)(nx-1);
//            ru[0] = (RealType)ix / (RealType)nx;
//            for (int iy=0; iy<ny; iy++){
//              //		ru[1] = (RealType)iy / (RealType)(ny-1);
//              ru[1] = (RealType)iy / (RealType)ny;
//              for (int iz=0; iz<nz; iz++) {
//                //		  ru[2] = (RealType)iz / (RealType)(nz-1);
//                ru[2] = (RealType)iz / (RealType)nz;
//                double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
//                double s, c;
//                sincos(phi, &s, &c);
//                std::complex<double> phase(c,s);
//                std::complex<double> z = phase*FFTbox(ix,iy,iz);
//
//                rNorm += z.real()*z.real();
//                iNorm += z.imag()*z.imag();
//              }
//            }
//          }
//          double arg = std::atan2(iNorm, rNorm);
//          //cerr << "Phase = " << arg/M_PI << " pi.\n";
//          double s,c;
//          sincos(0.5*(0.25*M_PI-arg), &s, &c);
//          std::complex<double> phase(c,s);
//          //rNorm=0.0; iNorm=0.0;
//          //#pragma omp parallel for reduction(+:rNorm,iNorm)
//#pragma omp parallel for
//          for (int ix=0; ix<nx; ix++)
//            for (int iy=0; iy<ny; iy++)
//              for (int iz=0; iz<nz; iz++) {
//                std::complex<double> z =
//                  splineData(ix,iy,iz) = phase*FFTbox(ix,iy,iz);
//                //rNorm += z.real()*z.real();
//                //iNorm += z.imag()*z.imag();
//              }
//          //arg = std::atan2(iNorm, rNorm);
//#endif
//        }
//      }
//      myComm->bcast(splineData);
//      myComm->bcast(twist);
//      set_multi_UBspline_3d_z(orbitalSet->MultiSpline, ival, splineData.data());
//    }
//
  for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
  {
    PosType twist=TwistAngles[SortBands[iorb].TwistIndex];
    // Read atomic orbital information
    for (int iat=0; iat<AtomicOrbitals.size(); iat++)
    {
      app_log() << "Reading orbital " << iat << " for band " << ival << std::endl;
      AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
      Array<std::complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm),
            poly_coefs(orb.PolyOrder+1,orb.Numlm);
      if (root)
      {
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
        std::ostringstream spline_path, poly_path;
        spline_path << path.str() << "radial_spline_" << iat;
        poly_path   << path.str() << "poly_coefs_"    << iat;
        HDFAttribIO<Array<std::complex<double>,2> > h_radial_spline(radial_spline);
        HDFAttribIO<Array<std::complex<double>,2> > h_poly_coefs(poly_coefs);
        h_radial_spline.read(H5FileID, spline_path.str().c_str());
        h_poly_coefs.read   (H5FileID, poly_path.str().c_str());
        // std::cerr << "radial_spline.size = (" << radial_spline.size(0)
        // 	 << ", " << radial_spline.size(1) << ")\n";
        // std::cerr << "poly_coefs.size = (" << poly_coefs.size(0)
        // 	 << ", " << poly_coefs.size(1) << ")\n";
      }
      myComm->bcast(radial_spline);
      myComm->bcast(poly_coefs);
      AtomicOrbitals[iat].set_band (ival, radial_spline, poly_coefs, twist);
    }
    // Now read muffin tin data
    for (int tin=0; tin<NumMuffinTins; tin++)
    {
      // app_log() << "Reading data for muffin tin " << tin << std::endl;
      PosType twist, k;
      int lmax = MT_APW_lmax[tin];
      int numYlm = (lmax+1)*(lmax+1);
      Array<std::complex<double>,2>
      u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
      Array<std::complex<double>,1> du_lm_dr (numYlm);
      if (root)
      {
        int ti   = SortBands[iorb].TwistIndex;
        int bi   = SortBands[iorb].BandIndex;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        std::string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
        std::string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
        HDFAttribIO<Array<std::complex<double>,2> > h_u_lm_r(u_lm_r);
        HDFAttribIO<Array<std::complex<double>,1> > h_du_lm_dr(du_lm_dr);
        h_u_lm_r.read(H5FileID, uName.c_str());
        h_du_lm_dr.read(H5FileID, duName.c_str());
      }
      myComm->bcast(u_lm_r);
      myComm->bcast(du_lm_dr);
      myComm->bcast(k);
      double Z = (double)IonTypes(tin);
      OrbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
    }
  }
  orbitalSet->AtomicOrbitals = AtomicOrbitals;
  for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
    orbitalSet->AtomicOrbitals[i].registerTimers();
  ExtendedMap_z[set] = orbitalSet->MultiSpline;
}



void
EinsplineSetBuilder::ReadBands_ESHDF
(int spin, EinsplineSetExtended<double>* orbitalSet)
{
  ReportEngine PRE("EinsplineSetBuilder","ReadBands_ESHDF(EinsplineSetExtended<double>*");
  std::vector<AtomicOrbital<double> > realOrbs(AtomicOrbitals.size());
  for (int iat=0; iat<realOrbs.size(); iat++)
  {
    AtomicOrbital<std::complex<double> > &corb (AtomicOrbitals[iat]);
    realOrbs[iat].set_pos  (corb.Pos);
    realOrbs[iat].set_lmax (corb.lMax);
    realOrbs[iat].set_cutoff (corb.CutoffRadius);
    realOrbs[iat].set_spline (corb.SplineRadius, corb.SplinePoints);
    realOrbs[iat].set_polynomial (corb.PolyRadius, corb.PolyOrder);
    realOrbs[iat].Lattice = corb.Lattice;
  }
  bool root = myComm->rank()==0;
  // bcast other stuff
  myComm->bcast (NumDistinctOrbitals);
  myComm->bcast (NumValenceOrbs);
  myComm->bcast (NumCoreOrbs);
  int N = NumDistinctOrbitals;
  orbitalSet->kPoints.resize(N);
  orbitalSet->MakeTwoCopies.resize(N);
  orbitalSet->StorageValueVector.resize(N);
  orbitalSet->BlendValueVector.resize(N);
  orbitalSet->StorageLaplVector.resize(N);
  orbitalSet->BlendLaplVector.resize(N);
  orbitalSet->StorageGradVector.resize(N);
  orbitalSet->BlendGradVector.resize(N);
  orbitalSet->StorageHessVector.resize(N);
  orbitalSet->StorageGradHessVector.resize(N);
  orbitalSet->phase.resize(N);
  orbitalSet->eikr.resize(N);
  orbitalSet->NumValenceOrbs = NumValenceOrbs;
  orbitalSet->NumCoreOrbs    = NumCoreOrbs;
  orbitalSet->FirstOrderSplines.resize(IonPos.size());
  // Read in k-points
  int numOrbs = orbitalSet->getOrbitalSetSize();
  int num = 0;
  if (root)
  {
    for (int iorb=0; iorb<N; iorb++)
    {
      int ti = SortBands[iorb].TwistIndex;
      PosType twist  = TwistAngles[ti];
      orbitalSet->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(twist);
      orbitalSet->MakeTwoCopies[iorb] =
        (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
      num += orbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
    }
    PosType twist0 = TwistAngles[SortBands[0].TwistIndex];
    for (int i=0; i<OHMMS_DIM; i++)
      if (std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)
        orbitalSet->HalfG[i] = 1;
      else
        orbitalSet->HalfG[i] = 0;
    EinsplineSetBuilder::RotateBands_ESHDF(spin, orbitalSet);
  }
  myComm->bcast(orbitalSet->kPoints);
  myComm->bcast(orbitalSet->MakeTwoCopies);
  myComm->bcast(orbitalSet->HalfG);
  // First, check to see if we have already read this in
  H5OrbSet set(H5FileName, spin, N);
  bool havePsir=!ReadGvectors_ESHDF();
  app_log() << "MeshSize = (" << MeshSize[0] << ", "
            << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  int nx, ny, nz, bi, ti;
  nx=MeshSize[0];
  ny=MeshSize[1];
  nz=MeshSize[2];
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  if (orbitalSet->HalfG[0])
  {
    xBC.lCode = ANTIPERIODIC;
    xBC.rCode = ANTIPERIODIC;
  }
  else
  {
    xBC.lCode = PERIODIC;
    xBC.rCode = PERIODIC;
  }
  if (orbitalSet->HalfG[1])
  {
    yBC.lCode = ANTIPERIODIC;
    yBC.rCode = ANTIPERIODIC;
  }
  else
  {
    yBC.lCode = PERIODIC;
    yBC.rCode = PERIODIC;
  }
  if (orbitalSet->HalfG[2])
  {
    zBC.lCode = ANTIPERIODIC;
    zBC.rCode = ANTIPERIODIC;
  }
  else
  {
    zBC.lCode = PERIODIC;
    zBC.rCode = PERIODIC;
  }
  x_grid.start = 0.0;
  x_grid.end = 1.0;
  x_grid.num = nx;
  y_grid.start = 0.0;
  y_grid.end = 1.0;
  y_grid.num = ny;
  z_grid.start = 0.0;
  z_grid.end = 1.0;
  z_grid.num = nz;
  // Create the multiUBspline object
  orbitalSet->MultiSpline =
    create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  if (HaveOrbDerivs)
  {
    orbitalSet->FirstOrderSplines.resize(IonPos.size());
    for (int ion=0; ion<IonPos.size(); ion++)
      for (int dir=0; dir<OHMMS_DIM; dir++)
        orbitalSet->FirstOrderSplines[ion][dir] =
          create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, NumValenceOrbs);
  }
  //////////////////////////////////////
  // Create the MuffinTin APW splines //
  //////////////////////////////////////
  orbitalSet->MuffinTins.resize(NumMuffinTins);
  for (int tin=0; tin<NumMuffinTins; tin++)
  {
    orbitalSet->MuffinTins[tin].Atom = tin;
    orbitalSet->MuffinTins[tin].set_center (MT_centers[tin]);
    orbitalSet->MuffinTins[tin].set_lattice(Lattice);
    orbitalSet->MuffinTins[tin].init_APW
    (MT_APW_rgrids[tin], MT_APW_lmax[tin],
     NumValenceOrbs);
  }
  for (int iat=0; iat<realOrbs.size(); iat++)
  {
    realOrbs[iat].set_num_bands(NumValenceOrbs);
    realOrbs[iat].allocate();
  }
  int isComplex;
  if (root)
  {
    HDFAttribIO<int> h_isComplex(isComplex);
    h_isComplex.read(H5FileID, "/electrons/psi_r_is_complex");
  }
  myComm->bcast(isComplex);
  bool isCore = bcastSortBands(N,root);
  if(isCore)
  {
    APP_ABORT("Core states not supported by ES-HDF yet.");
  }
  //this is common
  Array<double,3> splineData(nx,ny,nz);
  if(havePsir)
  {
    if(isComplex)
    {
      app_log() << "   Reading complex psi_r and convert to real" << std::endl;
      Array<std::complex<double>,3> rawData;
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        if(root)
        {
          std::ostringstream path;
          path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
               << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
          HDFAttribIO<Array<std::complex<double>,3> >  h_splineData(rawData);
          h_splineData.read(H5FileID, path.str().c_str());
        }
        myComm->bcast(rawData);
        //multiply twist factor and project on the real
        fix_phase_c2r(rawData,splineData,TwistAngles[ti]);
        set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
      }
    }
    else
    {
      app_log() << "   Reading real psi_r" << std::endl;
      for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
      {
        if(root)
        {
          std::ostringstream path;
          path << "/electrons/kpoint_" << SortBands[iorb].TwistIndex
               << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_r";
          HDFAttribIO<Array<double,3> >  h_splineData(splineData);
          h_splineData.read(H5FileID, path.str().c_str());
        }
        myComm->bcast(splineData);
        set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
      }
    }
  }
  else
  {
    FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
    FFTplan = fftw_plan_dft_3d
              (MeshSize[0], MeshSize[1], MeshSize[2],
               reinterpret_cast<fftw_complex*>(FFTbox.data()),
               reinterpret_cast<fftw_complex*>(FFTbox.data()),
               +1, FFTW_ESTIMATE);
    for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
    {
      Vector<std::complex<double> > cG;
      int ncg=0;
      int ti=SortBands[iorb].TwistIndex;
      if(root)
      {
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti    //SortBands[iorb].TwistIndex
             << "/spin_" << spin << "/state_" << SortBands[iorb].BandIndex << "/psi_g";
        std::cout << "FFT on " <<  path.str() << std::endl;
        HDFAttribIO<Vector<std::complex<double> > >  h_cG(cG);
        h_cG.read (H5FileID, path.str().c_str());
        ncg=cG.size();
      }
      myComm->bcast(ncg);
      if(ncg != Gvecs[ti].size())
      {
        APP_ABORT("Failed : ncg != Gvecs[ti].size()");
      }
      if(!root)
        cG.resize(ncg);
      myComm->bcast(cG);
      unpack4fftw(cG,Gvecs[ti],MeshSize,FFTbox);
      fftw_execute (FFTplan);
      fix_phase_rotate(FFTbox,splineData,TwistAngles[ti]);
      set_multi_UBspline_3d_d (orbitalSet->MultiSpline, ival, splineData.data());
    }
    fftw_destroy_plan(FFTplan);
  }
  for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
  {
    // Read atomic orbital information
    for (int iat=0; iat<realOrbs.size(); iat++)
    {
      app_log() << "Reading orbital " << iat << " for band " << ival << std::endl;
      AtomicOrbital<double> &orb = realOrbs[iat];
      //AtomicOrbital<std::complex<double> > &orb = realOrbs[iat];
      Array<std::complex<double>,2> radial_spline(orb.SplinePoints,orb.Numlm),
            poly_coefs(orb.PolyOrder+1,orb.Numlm);
      int ti   = SortBands[iorb].TwistIndex;
      if (root)
      {
        int bi   = SortBands[iorb].BandIndex;
        std::ostringstream path;
        path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/";
        std::ostringstream spline_path, poly_path;
        spline_path << path.str() << "radial_spline_" << iat;
        poly_path   << path.str() << "poly_coefs_"    << iat;
        HDFAttribIO<Array<std::complex<double>,2> > h_radial_spline(radial_spline);
        HDFAttribIO<Array<std::complex<double>,2> > h_poly_coefs(poly_coefs);
        h_radial_spline.read(H5FileID, spline_path.str().c_str());
        h_poly_coefs.read   (H5FileID, poly_path.str().c_str());
      }
      myComm->bcast(radial_spline);
      myComm->bcast(poly_coefs);
      realOrbs[iat].set_band (ival, radial_spline, poly_coefs, TwistAngles[ti]);
    }
  }
  for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
  {
    // Now read muffin tin data
    for (int tin=0; tin<NumMuffinTins; tin++)
    {
      // app_log() << "Reading data for muffin tin " << tin << std::endl;
      PosType twist, k;
      int lmax = MT_APW_lmax[tin];
      int numYlm = (lmax+1)*(lmax+1);
      Array<std::complex<double>,2>
      u_lm_r(numYlm, MT_APW_num_radial_points[tin]);
      Array<std::complex<double>,1> du_lm_dr (numYlm);
      int ti   = SortBands[iorb].TwistIndex;
      if (root)
      {
        int bi   = SortBands[iorb].BandIndex;
        twist = TwistAngles[ti];
        k = orbitalSet->PrimLattice.k_cart(twist);
        std::string uName  = MuffinTinPath (ti, bi,tin) + "u_lm_r";
        std::string duName = MuffinTinPath (ti, bi,tin) + "du_lm_dr";
        HDFAttribIO<Array<std::complex<double>,2> > h_u_lm_r(u_lm_r);
        HDFAttribIO<Array<std::complex<double>,1> > h_du_lm_dr(du_lm_dr);
        h_u_lm_r.read(H5FileID, uName.c_str());
        h_du_lm_dr.read(H5FileID, duName.c_str());
      }
      myComm->bcast(u_lm_r);
      myComm->bcast(du_lm_dr);
      myComm->bcast(k);
      double Z = (double)IonTypes(tin);
      OrbitalSet->MuffinTins[tin].set_APW (ival, k, u_lm_r, du_lm_dr, Z);
    }
  }
  //FIX HaveOrbDerivs after debugging
//	// Now read orbital derivatives if we have them
//	if (HaveOrbDerivs) {
//	  for (int ion=0; ion<IonPos.size(); ion++)
//	    for (int dim=0; dim<OHMMS_DIM; dim++) {
//	      if (root) {
//		int ti   = SortBands[iorb].TwistIndex;
//		int bi   = SortBands[iorb].BandIndex;
//
//		app_log() << "Reading orbital derivative for ion " << ion
//			  << " dim " << dim << " spin " << spin << " band "
//			  << bi << " kpoint " << ti << std::endl;
//		ostringstream path;
//		path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << bi << "/"
//		     << "dpsi_" << ion << "_" << dim << "_r";
//		string psirName = path.str();
//		if (isComplex) {
//		  HDFAttribIO<Array<std::complex<double>,3> > h_rawData(rawData);
//		  h_rawData.read(H5FileID, psirName.c_str());
//		  if ((rawData.size(0) != nx) ||
//		      (rawData.size(1) != ny) ||
//		      (rawData.size(2) != nz)) {
//		    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
//		    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
//		    abort();
//		  }
//#pragma omp parallel for
//		  for (int ix=0; ix<nx; ix++) {
//		  PosType ru;
//		    ru[0] = (RealType)ix / (RealType)nx;
//		    for (int iy=0; iy<ny; iy++) {
//		      ru[1] = (RealType)iy / (RealType)ny;
//		      for (int iz=0; iz<nz; iz++) {
//			ru[2] = (RealType)iz / (RealType)nz;
//			double phi = -2.0*M_PI*dot (ru, TwistAngles[ti]);
//			double s, c;
//			sincos(phi, &s, &c);
//			complex<double> phase(c,s);
//			complex<double> z = phase*rawData(ix,iy,iz);
//			splineData(ix,iy,iz) = z.real();
//		      }
//		    }
//		  }
//		}
//		else {
//		  HDFAttribIO<Array<double,3> >  h_splineData(splineData);
//		  h_splineData.read(H5FileID, psirName.c_str());
//		  if ((splineData.size(0) != nx) ||
//		      (splineData.size(1) != ny) ||
//		      (splineData.size(2) != nz)) {
//		    fprintf (stderr, "Error in EinsplineSetBuilder::ReadBands.\n");
//		    fprintf (stderr, "Extended orbitals should all have the same dimensions\n");
//		    abort();
//		  }
//		}
//	      }
//	      myComm->bcast(splineData);
//	      set_multi_UBspline_3d_d
//		(orbitalSet->FirstOrderSplines[ion][dim], ival, splineData.data());
//	    }
//	}
//
//
//
  orbitalSet->AtomicOrbitals = realOrbs;
  for (int i=0; i<orbitalSet->AtomicOrbitals.size(); i++)
    orbitalSet->AtomicOrbitals[i].registerTimers();
  ExtendedMap_d[set] = orbitalSet->MultiSpline;
}

}

