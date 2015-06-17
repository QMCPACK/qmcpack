//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////
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
#include <fftw3.h>
#include <Utilities/ProgressReportEngine.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/EinsplineAdoptor.h"
#include "QMCWaveFunctions/SplineC2XAdoptor.h"
#include "QMCWaveFunctions/SplineR2RAdoptor.h"
#include "QMCWaveFunctions/SplineMixedAdoptor.h"

#include "QMCWaveFunctions/BsplineReaderBase.h"
#include "QMCWaveFunctions/SplineAdoptorReaderP.h"
#include "QMCWaveFunctions/SplineAdoptorReaderInterface.h"
#include "QMCWaveFunctions/SplineMixedAdoptorReaderP.h"

#include "Interfaces/ESInterfaceBase.h"
#include "Interfaces/ESHDF5/ESHDF5Interface.h"
namespace qmcplusplus
{

extern bool sortByIndex(BandInfo leftB, BandInfo rightB);

bool EinsplineSetBuilder::ReadOrbitalInfo_Interface(ESInterfaceBase* interface)
{
  update_token(__FILE__,__LINE__,"ReadOrbitalInfo_Interface");
  interface->getVersion();
//  app_log() << "  Reading orbital file in ESHDF format.\n";
//  HDFAttribIO<TinyVector<int,3> > h_version(Version);
//  h_version.read (H5FileID, "/version");
//  app_log() << "  ESHDF orbital file version "
//            << Version[0] << "." << Version[1] <<"." << Version[2] << endl;
//  HDFAttribIO<Tensor<double,3> > h_Lattice(Lattice);
//  h_Lattice.read      (H5FileID, "/supercell/primitive_vectors");
  interface->getPrimVecs(Lattice);
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
  int have_dpsi = false;
  int NumAtomicOrbitals = 0;
  NumCoreStates = NumMuffinTins = NumTwists = NumSpins = NumBands = NumAtomicOrbitals = 0;
  //vector<int> nels_spin(2);
  //nels_spin[0]=TargetPtcl.last(0)-TargetPtcl.first(0);
  //nels_spin[1]=TargetPtcl.getTotalNum()-nels_spin[0];
  NumElectrons=TargetPtcl.getTotalNum();
//  HDFAttribIO<int> h_NumBands(NumBands),
//              h_NumSpins(NumSpins), h_NumTwists(NumTwists), h_NumCore(NumCoreStates),
//              h_NumMuffinTins(NumMuffinTins), h_have_dpsi(have_dpsi),
//              h_NumAtomicOrbitals(NumAtomicOrbitals);
//  h_NumBands.read      (H5FileID, "/electrons/kpoint_0/spin_0/number_of_states");
//  h_NumCore.read       (H5FileID, "/electrons/kpoint_0/spin_0/number_of_core_states");
  //h_NumElectrons.read  (H5FileID, "/electrons/number_of_electrons");
//  h_NumSpins.read      (H5FileID, "/electrons/number_of_spins");
//  h_NumTwists.read     (H5FileID, "/electrons/number_of_kpoints");
//  h_NumMuffinTins.read (H5FileID, "/muffin_tins/number_of_tins");
//  h_have_dpsi.read     (H5FileID, "/electrons/have_dpsi");
//  h_NumAtomicOrbitals.read (H5FileID, "/electrons/number_of_atomic_orbitals");

  NumBands = interface->getNumBands();
  NumSpins = interface->getNumSpins();
  NumTwists = interface->getNumTwists();
  NumCoreStates = interface->getNumCoreStates();
  NumMuffinTins = interface->getNumMuffinTins();
  have_dpsi = interface->getHaveDPsi();
  NumAtomicOrbitals = interface->getNumAtomicOrbitals();

  HaveOrbDerivs = have_dpsi;
  app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons
            << ", spins=" << NumSpins << ", twists=" << NumTwists
            << ", muffin tins=" << NumMuffinTins
            << ", core states=" << NumCoreStates << endl;
  app_log() << "atomic orbital=" << NumAtomicOrbitals << endl;
  if (TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1)
    app_log() << "  Using a " << TileFactor[0] << "x" << TileFactor[1]
              << "x" << TileFactor[2] << " tiling factor.\n";
  //////////////////////////////////
  // Read ion types and locations //
  //////////////////////////////////
  Vector<int> species_ids;
  interface->getSpeciesIDs(species_ids);

  int num_species = species_ids.size();
  //HDFAttribIO<Vector<int> > h_species_ids(species_ids);
  //h_species_ids.read (H5FileID, "/atoms/species_ids");
//  int num_species(0);
//  num_species = interface->getNumSpecies();
  //HDFAttribIO<int> h_num_species(num_species);
  //h_num_species.read (H5FileID, "/atoms/number_of_species");
  vector<int> atomic_numbers(num_species);
//  for (int isp=0; isp<num_species; isp++)
//  {
//    ostringstream name;
//    name << "/atoms/species_" << isp << "/atomic_number";
//    HDFAttribIO<int> h_atomic_number (atomic_numbers[isp]);
//    h_atomic_number.read(H5FileID, name.str().c_str());
//  }

  interface->getAtomicNumbers(atomic_numbers);
 // int num_species = interface->getNumSpecies();
  IonTypes.resize(species_ids.size());
  for (int i=0; i<species_ids.size(); i++)
    IonTypes[i] = atomic_numbers[species_ids[i]];
 // HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
 // h_IonPos.read   (H5FileID, "/atoms/positions");
  interface->getIonPositions(IonPos);

  for (int i=0; i<IonTypes.size(); i++)
    app_log() << "Atom type(" << i << ") = " << IonTypes(i) << endl;
  /////////////////////////////////////
  // Read atomic orbital information //
  /////////////////////////////////////
 // AtomicOrbitals.resize(NumAtomicOrbitals);
  interface->getAtomicOrbitals(AtomicOrbitals);

/*  for (int iat=0; iat<NumAtomicOrbitals; iat++)
  {
    AtomicOrbital<complex<double> > &orb = AtomicOrbitals[iat];
    int lmax, polynomial_order, spline_points;
    RealType cutoff_radius, polynomial_radius, spline_radius;
    PosType position;
    HDFAttribIO<int> h_lmax(lmax), h_polynomial_order(polynomial_order),
                h_spline_points(spline_points);
    HDFAttribIO<RealType> h_cutoff_radius(cutoff_radius),
                h_polynomial_radius(polynomial_radius),
                h_spline_radius(spline_radius);
    HDFAttribIO<PosType> h_position(position);
    ostringstream groupstream;
    groupstream << "/electrons/atomic_orbital_" << iat << "/";
    string groupname = groupstream.str();
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
  }*/
  ///////////////////////////
  // Read the twist angles //
  ///////////////////////////
  interface->getTwistData(TwistAngles, TwistWeight, TwistSymmetry);

/*  TwistAngles.resize(NumTwists);
  TwistSymmetry.resize(NumTwists);
  TwistWeight.resize(NumTwists);
  for (int ti=0; ti<NumTwists; ti++)
  {
    ostringstream path;
    path << "/electrons/kpoint_" << ti << "/reduced_k";
    HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
    h_Twist.read (H5FileID, path.str().c_str());
    if ((Version[0] >= 2) and (Version[1] >= 1))
    {
      ostringstream sym_path;
      sym_path << "/electrons/kpoint_" << ti << "/symgroup";
      HDFAttribIO<int> h_Sym(TwistSymmetry[ti]);
      h_Sym.read (H5FileID, sym_path.str().c_str());
      ostringstream nsym_path;
      nsym_path << "/electrons/kpoint_" << ti << "/numsym";
      HDFAttribIO<int> h_Nsym(TwistWeight[ti]);
      h_Nsym.read (H5FileID, nsym_path.str().c_str());
    }
    // Early versions from wfconv had wrong sign convention for
    // k-points.  EinsplineSet uses the opposite sign convention
    // from most DFT codes.
    if (Version[0] >= 2)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        TwistAngles[ti][dim] *= -1.0;
    //       snprintf (buff, 1000, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n",
    //           TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    //       app_log() << buff;
  } */

  if(qmc_common.use_density)
  {
    APP_ABORT("So...  Density not implemented yet in interface.  Because I'm lazy")
    //////////////////////////////////////////////////////////
    // Only if it is bulk: If the density has not been set in TargetPtcl, and   //
    // the density is available, read it in and save it     //
    // in TargetPtcl.                                       //
    //////////////////////////////////////////////////////////
    if(TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_BULK)
    {
      // FIXME:  add support for more than one spin density
      if (!TargetPtcl.Density_G.size())
      {
        HDFAttribIO<vector<TinyVector<int,OHMMS_DIM> > >
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
          ostringstream density_r_path, density_g_path;
          density_r_path << "/electrons/density/spin_" << ispin << "/density_r";
          density_g_path << "/electrons/density/spin_" << ispin << "/density_g";
          h_density_r.read (H5FileID, density_r_path.str().c_str());
          if (TargetPtcl.DensityReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
            vector<ComplexType> density_G;
            HDFAttribIO<vector<ComplexType > > h_density_G (density_G);
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
        HDFAttribIO<vector<TinyVector<int,OHMMS_DIM> > >
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
          ostringstream VHXC_r_path, VHXC_g_path;
          VHXC_r_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_r";
          VHXC_g_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_g";
          h_VHXC_r.read (H5FileID, VHXC_r_path.str().c_str());
          if (TargetPtcl.VHXCReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found VHXC in the HDF5 file.\n";
            vector<ComplexType> VHXC_G;
            HDFAttribIO<vector<ComplexType > > h_VHXC_G (VHXC_G);
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
    }
  }
  else
  {
    app_log() << "   Skip initialization of the density" << endl;
  }
  return true;
}


void EinsplineSetBuilder::OccupyBands_Interface(ESInterfaceBase* interface, int spin, int sortBands, int numOrbs)
{
  update_token(__FILE__,__LINE__,"OccupyBands_Interface");
  if (myComm->rank() != 0)
    return;

  vector<BandInfo>& SortBands(*FullBands[spin]);
  SortBands.clear(); //??? can exit if SortBands is already made?
  int maxOrbs(0);
  for (int ti=0; ti<DistinctTwists.size(); ti++)
  {
    int tindex = DistinctTwists[ti];
//     First, read valence states
//    ostringstream ePath;
//    ePath << "/electrons/kpoint_" << tindex << "/spin_"
//          << spin << "/eigenvalues";
    vector<double> eigvals;
//    HDFAttribIO<vector<double> > h_eigvals(eigvals);
//    h_eigvals.read(H5FileID, ePath.str().c_str());
    interface->getOrbEigenvals(spin,ti,eigvals);
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
      APP_ABORT("Core states not supported with interface yet")
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
  else if (sortBands==1)
  {
    app_log() << "Sorting the bands now:\n";
    sort (SortBands.begin(), SortBands.end());
  }

  vector<int> gsOcc(maxOrbs);
  int N_gs_orbs=numOrbs;
  int nocced(0), ntoshift(0);
  for (int ti=0; ti<SortBands.size(); ti++)
  {
    if (nocced<N_gs_orbs)
    {
      if (SortBands[ti].MakeTwoCopies && (N_gs_orbs-nocced>1))
      {
        nocced+=2;
        ntoshift++;
        gsOcc[ti]=2;
      }
      else if ( (SortBands[ti].MakeTwoCopies && (N_gs_orbs-nocced==1)) || !SortBands[ti].MakeTwoCopies )
      {
        nocced+=1;
        ntoshift++;
        gsOcc[ti]=1;
      }
    }
  }
//    if(qafm!=0)
//    {
//      app_log()<<"Finding AFM pair for first "<<ntoshift<<" orbitals."<<endl;
//
//      for (int ti=0; ti<ntoshift; ti++)
//      {
//        bool found(false);
//        PosType ku = TwistAngles[SortBands[ti].TwistIndex];
//        PosType k1 = OrbitalSet->PrimLattice.k_cart(ku);
//        for (int tj=0; tj<TwistAngles.size(); tj++)
//        {
//          if(tj!=SortBands[ti].TwistIndex)
//          {
//            ku=TwistAngles[tj];
//            PosType k2 = OrbitalSet->PrimLattice.k_cart(ku);
//            double dkx = abs(k1[0] - k2[0]);
//            double dky = abs(k1[1] - k2[1]);
//            double dkz = abs(k1[2] - k2[2]);
//            bool rightK = ((dkx<qafm+0.0001)&&(dkx>qafm-0.0001)&&(dky<0.0001)&&(dkz<0.0001));
//            if(rightK)
//            {
//              SortBands[ti].TwistIndex = tj;
//              //               app_log()<<"swapping: "<<ti<<" "<<tj<<endl;
//              found=true;
//              break;
//            }
//          }
//        }
//        if(!found)
//        {
//          if((abs(k1[1])<qafm+0.0001)&&(abs(k1[1])>qafm-0.0001)) k1[1]*=-1;
//          else if((abs(k1[2])<qafm+0.0001)&&(abs(k1[2])>qafm-0.0001)) k1[2]*=-1;
//
//          for (int tj=0; tj<TwistAngles.size(); tj++)
//          {
//            if(tj!=SortBands[ti].TwistIndex)
//            {
//              ku=TwistAngles[tj];
//              PosType k2 = OrbitalSet->PrimLattice.k_cart(ku);
//              double dkx = abs(k1[0] - k2[0]);
//              double dky = abs(k1[1] - k2[1]);
//              double dkz = abs(k1[2] - k2[2]);
//              bool rightK = ((dkx<qafm+0.0001)&&(dkx>qafm-0.0001)&&(dky<0.0001)&&(dkz<0.0001));
//              if(rightK)
//              {
//                SortBands[ti].TwistIndex = tj;
//                //               app_log()<<"swapping: "<<ti<<" "<<tj<<endl;
//                found=true;
//                break;
//              }
//            }
//          }
//        }
//
//        if(!found)
//        {
//          app_log()<<"Need twist: ("<<k1[0]+qafm<<","<<k1[1]<<","<<k1[2]<<")"<<endl;
//          app_log()<<"Did not find afm pair for orbital: "<<ti<<", twist index: "<<SortBands[ti].TwistIndex<<endl;
//          APP_ABORT("EinsplineSetBuilder::OccupyBands_ESHDF");
//        }
//      }
//    }
  if (occ_format=="energy")
  {
    // To get the occupations right.
    vector<int> Removed(0,0);
    vector<int> Added(0,0);
    for(int ien=0; ien<Occ.size(); ien++)
    {
      if (Occ[ien]<0)
        Removed.push_back(-Occ[ien]);
      else if (Occ[ien]>0)
        Added.push_back(Occ[ien]);
    }
    if(Added.size()-Removed.size() != 0)
    {
      app_log()<<"need to add and remove same number of orbitals. "<< Added.size()<<" "<<Removed.size()<<endl;
      APP_ABORT("ChangedOccupations");
    }
    vector<int> DiffOcc(maxOrbs,0);
    //Probably a cleaner way to do this.
    for(int i=0; i<Removed.size(); i++)
      DiffOcc[Removed[i]-1]-=1;
    for(int i=0; i<Added.size(); i++)
      DiffOcc[Added[i]-1]+=1;
    vector<int> SumOrb(SortBands.size(),0);
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
    vector<BandInfo> ReOrderedBands;
    vector<BandInfo> RejectedBands;
    for(int i=0; i<SumOrb.size(); i++)
    {
      if(SumOrb[i]==2)
      {
        SortBands[i].MakeTwoCopies=true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i]==1)
      {
        SortBands[i].MakeTwoCopies=false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i]==0)
      {
        SortBands[i].MakeTwoCopies=false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log()<<" Trying to add the same orbital ("<<i<<") less than zero or more than 2 times."<<endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(),RejectedBands.begin(),RejectedBands.end());
    SortBands=ReOrderedBands;
  }
  else if (occ_format=="band")
  {
    app_log()<<"  Occupying bands based on (bi,ti) data."<<endl;
    if(Occ.size() != particle_hole_pairs*4)
    {
      app_log()<<" Need Occ = pairs*4. Occ is (ti,bi) of removed, then added."<<endl;
      app_log()<<Occ.size()<<" "<<particle_hole_pairs<<endl;
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
          app_log()<<"removing orbital "<<ien<<endl;
        }
        else
        {
          gsOcc[ien]+=1;
          app_log()<<"adding orbital "<<ien<<endl;
          cnt+=2;
        }
    }
    vector<BandInfo> ReOrderedBands;
    vector<BandInfo> RejectedBands;
    for(int i=0; i<SortBands.size(); i++)
    {
      if(gsOcc[i]==2)
      {
        SortBands[i].MakeTwoCopies=true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i]==1)
      {
        SortBands[i].MakeTwoCopies=false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i]==0)
      {
        SortBands[i].MakeTwoCopies=false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log()<<" Trying to add the same orbital ("<<i<<") less than zero or more than 2 times."<<endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(),RejectedBands.begin(),RejectedBands.end());
    SortBands=ReOrderedBands;
  }
  //for(int sw=0;sw<Removed.size();sw++){
  //  app_log()<<" Swapping two orbitals "<<Removed[sw]<<" and "<<Added[sw]<<endl;
  //  BandInfo tempband(SortBands[Removed[sw]-1]);
  //  SortBands[Removed[sw]-1] = SortBands[Added[sw]-1];
  //  SortBands[Added[sw]-1] = tempband;
  //}
  int orbIndex = 0;
  int numOrbs_counter = 0;
  NumValenceOrbs=0;
  NumCoreOrbs=0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
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

/*SPOSetBase*
EinsplineSetBuilder::createSPOSetFromInterface(xmlNodePtr cur)
{
  update_token(__FILE__,__LINE__,"createSPOSetFromInterface");
  //use 2 bohr as the default when truncated orbitals are used based on the extend of the ions
  BufferLayer=2.0;
  SPOSetBase *OrbitalSet;
  int numOrbs = 0;
  qafm=0;
  int sortBands(1);
  int spinSet = 0;

  string sourceName;
  string spo_prec("double");
  string truncate("no");
#if defined(QMC_CUDA)
  string useGPU="yes";
#else
  string useGPU="no";
#endif

  {
    OhmmsAttributeSet a;
    a.add (H5FileName, "href");
    a.add (TileFactor, "tile");
    a.add (sortBands,  "sort");
    a.add (qafm,  "afmshift");
    a.add (TileMatrix, "tilematrix");
    a.add (TwistNum,   "twistnum");
    a.add (givenTwist,   "twist");
    a.add (sourceName, "source");
    a.add (MeshFactor, "meshfactor");
    a.add (useGPU,     "gpu");
    a.add (spo_prec,   "precision");
    a.add (truncate,   "truncate");
    a.add (BufferLayer, "buffer");
    a.add (myName, "tag");

    a.put (XMLRoot);
    a.add (numOrbs,    "size");
    a.add (numOrbs,    "norbs");
    a.add(spinSet,"spindataset"); a.add(spinSet,"group");
    a.put (cur);

    if(myName.empty()) myName="einspline";

  }

   H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  ESHDF5Interface* interface = new ESHDF5Interface(myComm);
  interface->put(XMLRoot);
  interface->initialize();

  SourcePtcl=ParticleSets[sourceName];
  if(SourcePtcl==0)
  {
    APP_ABORT("Einspline needs the source particleset");
  }
  else
  { //keep the one-body distance table index 
    myTableIndex=TargetPtcl.addTable(*SourcePtcl);
  }

  ///////////////////////////////////////////////
  // Read occupation information from XML file //
  ///////////////////////////////////////////////
  vector<int> Occ_Old(0,0);
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
    string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      string occ_mode("ground");
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
  app_log() << "\t  QMC_CUDA=1 Overwriting the precision of the einspline storage on the host.\n";
  spo_prec="double"; //overwrite
  truncate="no"; //overwrite
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

  if(FullBands[spinSet]==0) FullBands[spinSet]=new vector<BandInfo>;

  Timer mytimer;
  if(spinSet==0)
  {
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
        TileMatrix(i,j) = (i==j) ? TileFactor(i) : 0;
  if (myComm->rank() == 0)
    fprintf (stderr, " [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d ]\n",
             TileMatrix(0,0), TileMatrix(0,1), TileMatrix(0,2),
             TileMatrix(1,0), TileMatrix(1,1), TileMatrix(1,2),
             TileMatrix(2,0), TileMatrix(2,1), TileMatrix(2,2));
  if (numOrbs == 0)
  {
    app_error() << "You must specify the number of orbitals in the input file.\n";
    APP_ABORT("EinsplineSetBuilder::createSPOSet");
  }
  else
    app_log() << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";
  mytimer.restart();
  /////////////////////////////////////////////////////////////////
  // Read the basic orbital information, without reading all the //
  // orbitals themselves.                                        //
  /////////////////////////////////////////////////////////////////
  if (myComm->rank() == 0)
    if (!ReadOrbitalInfo_ESHDF())
   // if (!ReadOrbitalInfo_Interface(interface))
    {
      app_error() << "Error reading orbital info from interface.  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::createSPOSet");
    }
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << mytimer.elapsed() << endl;
  myComm->barrier();
  mytimer.restart();
  BroadcastOrbitalInfo();

  app_log() <<  "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << mytimer.elapsed() << endl;
  app_log().flush();

  // Now, analyze the k-point mesh to figure out the what k-points  are needed                                                    //
  PrimCell.set(Lattice);
  SuperCell.set(SuperLattice);
  GGt=dot(transpose(PrimCell.G), PrimCell.G);

  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
    AtomicOrbitals[iat].Lattice = Lattice;
  // Copy supercell into the ParticleSets
//     app_log() << "Overwriting XML lattice with that from the ESHDF file.\n";
//     PtclPoolType::iterator piter;
//     for(piter = ParticleSets.begin(); piter != ParticleSets.end(); piter++)
//       piter->second->Lattice.copy(SuperCell);
  AnalyzeTwists2();

  } //use spinSet==0 to initialize shared properties of orbitals

  //////////////////////////////////
  // Create the OrbitalSet object
  //////////////////////////////////
  mytimer.restart();

  //OccupyBands_Interface(interface, spinSet, sortBands, numOrbs);
  OccupyBands_ESHDF(spinSet, sortBands, numOrbs);
//#ifdef QMC_CUDA.  RAY: No CUDA for now.

  // the CPU branch
  if(spinSet==0) TileIons();

  bool use_single= (spo_prec == "single" || spo_prec == "float");

  if (UseRealOrbitals)
  {
    //RAY:  No mixed spline reader, no truncated orbitals, and no single precision assumed.  Following commented lines remove that functionality.
    //See CreateSPOSetFromXML for removed block.  

    MixedSplineReader= new SplineAdoptorReader<SplineR2RAdoptor<double,double,3> >(this);  //RAY:  Moved from above block.  
//    MixedSplineReader= new SplineAdoptorReaderInterface<SplineR2RAdoptor<double,double,3> >(this,interface);  //RAY:  Moved from above block.  
    MixedSplineReader->setCommon(XMLRoot);
    HasCoreOrbs=bcastSortBands(spinSet,NumDistinctOrbitals,myComm->rank()==0);
    SPOSetBase* bspline_zd=MixedSplineReader->create_spline_set(spinSet,spo_cur);
    if(bspline_zd)
      SPOSetMap[aset] = bspline_zd;
    else
      APP_ABORT_TRACE(__FILE__,__LINE__,"Failed to create SPOSetBase*");

    OrbitalSet = bspline_zd;
  }
  else
  {

    MixedSplineReader= new SplineAdoptorReader<SplineC2RPackedAdoptor<double,double,3> >(this);
    //MixedSplineReader= new SplineAdoptorReaderInterface<SplineC2RPackedAdoptor<double,double,3> >(this,interface);

    MixedSplineReader->setCommon(XMLRoot);
    size_t delta_mem=qmc_common.memory_allocated;
    // temporary disable the following function call, Ye Luo
    // RotateBands_ESHDF(spinSet, dynamic_cast<EinsplineSetExtended<complex<double> >*>(OrbitalSet));
    HasCoreOrbs=bcastSortBands(spinSet,NumDistinctOrbitals,myComm->rank()==0);
    SPOSetBase* bspline_zd=MixedSplineReader->create_spline_set(spinSet,spo_cur);
    if(bspline_zd)
      SPOSetMap[aset] = bspline_zd;
    else
      APP_ABORT_TRACE(__FILE__,__LINE__,"Failed to create SPOSetBase*");

    delta_mem=qmc_common.memory_allocated-delta_mem;
    app_log() <<"  MEMORY allocated SplineAdoptorReader " << (delta_mem>>20) << " MB" << endl;
    OrbitalSet = bspline_zd;
  }
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadBands " << mytimer.elapsed() << endl;
  SPOSetMap[aset] = OrbitalSet;

  interface->finalize();
  delete interface;

  return OrbitalSet;
	
}*/


SPOSetBase*
EinsplineSetBuilder::createSPOSetFromInterface(xmlNodePtr cur)
{

  update_token(__FILE__,__LINE__,"createSPOSetFromInterface");
  //use 2 bohr as the default when truncated orbitals are used based on the extend of the ions
  BufferLayer=2.0;
  SPOSetBase *OrbitalSet;
  int numOrbs = 0;
  qafm=0;
  int sortBands(1);
  int spinSet = 0;

  string sourceName;
  string spo_prec("double");
  string truncate("no");
#if defined(QMC_CUDA)
  string useGPU="yes";
#else
  string useGPU="no";
#endif

  {
    OhmmsAttributeSet a;
    a.add (H5FileName, "href");
    a.add (TileFactor, "tile");
    a.add (sortBands,  "sort");
    a.add (qafm,  "afmshift");
    a.add (TileMatrix, "tilematrix");
    a.add (TwistNum,   "twistnum");
    a.add (givenTwist,   "twist");
    a.add (sourceName, "source");
    a.add (MeshFactor, "meshfactor");
    a.add (useGPU,     "gpu");
    a.add (spo_prec,   "precision");
    a.add (truncate,   "truncate");
    a.add (BufferLayer, "buffer");
    a.add (myName, "tag");

    a.put (XMLRoot);
    a.add (numOrbs,    "size");
    a.add (numOrbs,    "norbs");
    a.add(spinSet,"spindataset"); a.add(spinSet,"group");
    a.put (cur);

    if(myName.empty()) myName="einspline";

  }

   H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  ESHDF5Interface* interface = new ESHDF5Interface(myComm);
  interface->put(XMLRoot);
  interface->initialize();

  SourcePtcl=ParticleSets[sourceName];
  if(SourcePtcl==0)
  {
    APP_ABORT("Einspline needs the source particleset");
  }
  else
  { //keep the one-body distance table index 
    myTableIndex=TargetPtcl.addTable(*SourcePtcl);
  }

  ///////////////////////////////////////////////
  // Read occupation information from XML file //
  ///////////////////////////////////////////////
  vector<int> Occ_Old(0,0);
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
    string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      string occ_mode("ground");
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
  app_log() << "\t  QMC_CUDA=1 Overwriting the precision of the einspline storage on the host.\n";
  spo_prec="double"; //overwrite
  truncate="no"; //overwrite
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

  if(FullBands[spinSet]==0) FullBands[spinSet]=new vector<BandInfo>;

  Timer mytimer;
  if(spinSet==0)
  {
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
        TileMatrix(i,j) = (i==j) ? TileFactor(i) : 0;
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
  mytimer.restart();
  /////////////////////////////////////////////////////////////////
  // Read the basic orbital information, without reading all the //
  // orbitals themselves.                                        //
  /////////////////////////////////////////////////////////////////
  if (myComm->rank() == 0)
  //  if(!ReadOrbitalInfo())
    if (!ReadOrbitalInfo_Interface(interface))
    {
      app_error() << "Error reading orbital info from HDF5 file.  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::createSPOSet");
    }
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << mytimer.elapsed() << endl;
  myComm->barrier();
  mytimer.restart();
  BroadcastOrbitalInfo();

  app_log() <<  "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << mytimer.elapsed() << endl;
  app_log().flush();

  // Now, analyze the k-point mesh to figure out the what k-points  are needed                                                    //
  PrimCell.set(Lattice);
  SuperCell.set(SuperLattice);
  GGt=dot(transpose(PrimCell.G), PrimCell.G);

  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
    AtomicOrbitals[iat].Lattice = Lattice;
  // Copy supercell into the ParticleSets
//     app_log() << "Overwriting XML lattice with that from the ESHDF file.\n";
//     PtclPoolType::iterator piter;
//     for(piter = ParticleSets.begin(); piter != ParticleSets.end(); piter++)
//       piter->second->Lattice.copy(SuperCell);
  AnalyzeTwists2();

  } //use spinSet==0 to initialize shared properties of orbitals

  //////////////////////////////////
  // Create the OrbitalSet object
  //////////////////////////////////
  mytimer.restart();

  //OccupyBands_ESHDF(spinSet, sortBands, numOrbs);
  OccupyBands_Interface(interface,spinSet, sortBands, numOrbs);
#ifdef QMC_CUDA
  // the GPU branch
/*  EinsplineSet *new_OrbitalSet;

  if (AtomicOrbitals.size() > 0)
  {
    if (UseRealOrbitals)
      new_OrbitalSet = new EinsplineSetHybrid<double>;
    else
      new_OrbitalSet = new EinsplineSetHybrid<complex<double> >;
  }
  else
  {
    if (UseRealOrbitals)
      new_OrbitalSet = new EinsplineSetExtended<double>;
    else
      new_OrbitalSet = new EinsplineSetExtended<complex<double> >;
  }

  //set the internal parameters
  if(spinSet==0) { setTiling(new_OrbitalSet,numOrbs); TileIons(); }

  OrbitalSet = new_OrbitalSet;

  if (UseRealOrbitals)
  {
    app_log() << ">>>> Creating EinsplineSetExtended<double> <<<< " << endl;
    EinsplineSetExtended<double> *restrict orbitalSet =
      dynamic_cast<EinsplineSetExtended<double>* > (OrbitalSet);
    if (Format == ESHDF)
      ReadBands_ESHDF(spinSet,orbitalSet);
    else
      ReadBands(spinSet, orbitalSet);
  }
  else
  {
    app_log() << ">>>> Creating EinsplineSetExtended<complex<double> > <<<< " << endl;
    EinsplineSetExtended<complex<double> > *restrict orbitalSet =
      dynamic_cast<EinsplineSetExtended<complex<double> >*>(OrbitalSet);
    if (Format == ESHDF)
      ReadBands_ESHDF(spinSet,orbitalSet);
    else
      ReadBands(spinSet, orbitalSet);
  }*/
#else
  // the CPU branch
  if(spinSet==0) TileIons();

  bool use_single= (spo_prec == "single" || spo_prec == "float");

  if (UseRealOrbitals)
  {
    //if(TargetPtcl.Lattice.SuperCellEnum != SUPERCELL_BULK && truncate=="yes")
    if(MixedSplineReader==0)
    {
      if(truncate=="yes")
      {
        if(use_single)
        {
//          if(TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN)
//            MixedSplineReader= new SplineMixedAdoptorReaderInterface<SplineOpenAdoptor<float,double,3> >(this,interface);
//          else if(TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_SLAB)
//            MixedSplineReader= new SplineMixedAdoptorReaderInterface<SplineMixedAdoptor<float,double,3> >(this,interface);
//          else
            MixedSplineReader= new SplineAdoptorReaderInterface<SplineR2RAdoptor<float,double,3> >(this,interface);
        }
        else
        {
//          if(TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN)
//            MixedSplineReader= new SplineMixedAdoptorReaderInterface<SplineOpenAdoptor<double,double,3> >(this,interface);
//          else if(TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_SLAB)
//            MixedSplineReader= new SplineMixedAdoptorReaderInterface<SplineMixedAdoptor<double,double,3> >(this,interface);
//          else
            MixedSplineReader= new SplineAdoptorReaderInterface<SplineR2RAdoptor<double,double,3> >(this,interface);
        }
      }
      else
      {
        if(use_single)
          MixedSplineReader= new SplineAdoptorReaderInterface<SplineR2RAdoptor<float,double,3> >(this,interface);
        else
          MixedSplineReader= new SplineAdoptorReaderInterface<SplineR2RAdoptor<double,double,3> >(this,interface);
      }
    }

    MixedSplineReader->setCommon(XMLRoot);
    HasCoreOrbs=bcastSortBands(spinSet,NumDistinctOrbitals,myComm->rank()==0);
    SPOSetBase* bspline_zd=MixedSplineReader->create_spline_set(spinSet,spo_cur);
    if(bspline_zd)
      SPOSetMap[aset] = bspline_zd;
    else
      APP_ABORT_TRACE(__FILE__,__LINE__,"Failed to create SPOSetBase*");

    OrbitalSet = bspline_zd;
  }
  else
  {
    if(MixedSplineReader==0)
    {
      if(truncate == "yes")
      {
        app_log() << "  Truncated orbitals with multiple kpoints are not supported yet!" << endl;
      }
      if(use_single)
      {
#if defined(QMC_COMPLEX)
        MixedSplineReader= new SplineAdoptorReaderInterface<SplineC2CPackedAdoptor<float,double,3> >(this,interface);
#else
        MixedSplineReader= new SplineAdoptorReaderInterface<SplineC2RPackedAdoptor<float,double,3> >(this,interface);
#endif
      }
      else
      {
#if defined(QMC_COMPLEX)
        MixedSplineReader= new SplineAdoptorReaderInterface<SplineC2CPackedAdoptor<double,double,3> >(this,interface);
#else
        MixedSplineReader= new SplineAdoptorReaderInterface<SplineC2RPackedAdoptor<double,double,3> >(this,interface);
#endif
      }
    }

    MixedSplineReader->setCommon(XMLRoot);
    size_t delta_mem=qmc_common.memory_allocated;
    // temporary disable the following function call, Ye Luo
    // RotateBands_ESHDF(spinSet, dynamic_cast<EinsplineSetExtended<complex<double> >*>(OrbitalSet));
    HasCoreOrbs=bcastSortBands(spinSet,NumDistinctOrbitals,myComm->rank()==0);
    SPOSetBase* bspline_zd=MixedSplineReader->create_spline_set(spinSet,spo_cur);
    if(bspline_zd)
      SPOSetMap[aset] = bspline_zd;
    else
      APP_ABORT_TRACE(__FILE__,__LINE__,"Failed to create SPOSetBase*");

    delta_mem=qmc_common.memory_allocated-delta_mem;
    app_log() <<"  MEMORY allocated SplineAdoptorReader " << (delta_mem>>20) << " MB" << endl;
    OrbitalSet = bspline_zd;
  }
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
      fprintf (fout, "%1.12e ", r*x/std::fabs(x));
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
// the end of CPU branch
#endif
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadBands " << mytimer.elapsed() << endl;
  SPOSetMap[aset] = OrbitalSet;
  //if (sourceName.size() && (ParticleSets.find(sourceName) == ParticleSets.end()))
  //{
  //  app_log() << "  EinsplineSetBuilder creates a ParticleSet " << sourceName << endl;
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
  return OrbitalSet;
}

bool EinsplineSetBuilder::ReadGvectors_Interface(ESInterfaceBase* interface)
{
  update_token(__FILE__,__LINE__,"ReadGvectors_Interface");
  bool root=myComm->rank() ==0;
  //this is always ugly
  MeshSize = 0;
  int hasPsig=1;

/*
#if defined(__bgp__)||(__bgq__)
  //Ray:  Replace with interface function calls.  
  if(root)
  {
    hid_t gid=H5Dopen(H5FileID,"/electrons/kpoint_0/spin_0/state_0/psi_g"); //Ray:  HWI
    if(gid<0)   
      hasPsig=0;
    H5Dclose(gid);  //Ray: HWI
  }
  //Ray
  myComm->bcast(hasPsig);
#else
  if(root)
  {
    //Ray:  Replace with interface function calls.  
    HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize); //Ray:  HWI
    h_mesh.read (H5FileID, "/electrons/psi_r_mesh");  //Ray:  HWI
    h_mesh.read (H5FileID, "/electrons/mesh"); //Ray:  HWI
  }
  myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
#endif
*/
  hasPsig = interface->getHasPsiG();
  interface->getMeshSize(MeshSize);
  if(hasPsig)
  {
    int nallowed=257;
    int allowed[] =
    {
      72,75,80,81,90,96,100,108,120,125,
      128,135,144,150,160,162,180,192,200,216,
      225,240,243,250,256,270,288,300,320,324,
      360,375,384,400,405,432,450,480,486,500,
      512,540,576,600,625,640,648,675,720,729,
      750,768,800,810,864,900,960,972,1000,1024,
      1080,1125,1152,1200,1215,1250,1280,1296,1350,1440,
      1458,1500,1536,1600,1620,1728,1800,1875,1920,1944,
      2000,2025,2048,2160,2187,2250,2304,2400,2430,2500,
      2560,2592,2700,2880,2916,3000,3072,3125,3200,3240,
      3375,3456,3600,3645,3750,3840,3888,4000,4050,4096,
      4320,4374,4500,4608,4800,4860,5000,5120,5184,5400,
      5625,5760,5832,6000,6075,6144,6250,6400,6480,6561,
      6750,6912,7200,7290,7500,7680,7776,8000,8100,8192,
      8640,8748,9000,9216,9375,9600,9720,10000,10125,10240,
      10368,10800,10935,11250,11520,11664,12000,12150,12288,12500,
      12800,12960,13122,13500,13824,14400,14580,15000,15360,15552,
      15625,16000,16200,16384,16875,17280,17496,18000,18225,18432,
      18750,19200,19440,19683,20000,20250,20480,20736,21600,21870,
      22500,23040,23328,24000,24300,24576,25000,25600,25920,26244,
      27000,27648,28125,28800,29160,30000,30375,30720,31104,31250,
      32000,32400,32768,32805,33750,34560,34992,36000,36450,36864,
      37500,38400,38880,39366,40000,40500,40960,41472,43200,43740,
      45000,46080,46656,46875,48000,48600,49152,50000,50625,51200,
      51840,52488,54000,54675,55296,56250,57600,58320,59049,60000,
      60750,61440,62208,62500,64000,64800,65536
    };
    int numk=0;
    MaxNumGvecs=0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int,3> maxIndex(0,0,0);
    Gvecs.resize(NumTwists);
    {
      int numg=0;
      if(root)
      {
        //Ray:  Replace with interface
//        ostringstream Gpath;  //Ray:  HWI
//        Gpath    << "/electrons/kpoint_0/gvectors";  //Ray:  HWI 
//        HDFAttribIO<vector<TinyVector<int,3> > > h_Gvecs(Gvecs[0]); //Ray:  HWI
//        h_Gvecs.read (H5FileID, Gpath.str().c_str()); //Ray:  HWI
        interface->getReducedGVecs(Gvecs);
        numg=Gvecs[0].size(); //Ray:  Replace with interface. 
      }
      myComm->bcast(numg);
      if(!root)
        Gvecs[0].resize(numg);
      myComm->bcast(Gvecs[0]);
      MaxNumGvecs=Gvecs[0].size();
      for (int ig=0; ig<Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } //done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0*MeshFactor*maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0*MeshFactor*maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0*MeshFactor*maxIndex[2]);
    //only use 2^a 3^b 5^c where a>=2  up to 65536
    int *ix=lower_bound(allowed,allowed+nallowed,MeshSize[0]);
    int *iy=lower_bound(allowed,allowed+nallowed,MeshSize[1]);
    int *iz=lower_bound(allowed,allowed+nallowed,MeshSize[2]);
    MeshSize[0]=(MeshSize[0]>128)? *ix:(MeshSize[0]+MeshSize[0]%2);
    MeshSize[1]=(MeshSize[1]>128)? *iy:(MeshSize[1]+MeshSize[1]%2);
    MeshSize[2]=(MeshSize[2]>128)? *iz:(MeshSize[2]+MeshSize[2]%2);
    if(Version[0]<2) //Ray:  HWI
    {
      //get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  ESHDF::Version " << Version << endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. Regenerate orbital files using updated QE." << endl;
      for(int k=0; k<DistinctTwists.size(); ++k)
      {
        int ik=DistinctTwists[k];
        if(ik==0)
          continue; //already done
        int numg=0;
        if(root)
        {
          //Ray:  Replace with interface
//          ostringstream Gpath; //Ray:  HWI
//          Gpath    << "/electrons/kpoint_" << ik << "/gvectors"; //Ray:  HWI
//          HDFAttribIO<vector<TinyVector<int,3> > > h_Gvecs(Gvecs[ik]); //Ray:  HWI
//          h_Gvecs.read (H5FileID, Gpath.str().c_str()); //Ray:  HWI
          interface->getReducedGVecs(Gvecs,ik);
          numg=Gvecs[ik].size(); //Ray
        }
        myComm->bcast(numg);
        if(numg==0)
        {
          //copy kpoint_0, default
          Gvecs[ik]=Gvecs[0];
        }
        else
        {
          if(numg !=  MaxNumGvecs)
          {
            ostringstream o;
            o<< "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
             << " This is not supported anymore. Rerun pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if(!root)
            Gvecs[ik].resize(numg);
          myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << endl;
  app_log().flush();
  return hasPsig;
}

}
