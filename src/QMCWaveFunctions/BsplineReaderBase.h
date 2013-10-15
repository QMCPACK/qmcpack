//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file BsplineReaderBase.h
 *
 * base class to read data and manage spline tables
 */
#ifndef QMCPLUSPLUS_BSPLINE_READER_BASE_H
#define QMCPLUSPLUS_BSPLINE_READER_BASE_H
#include <mpi/collectives.h>
#include <mpi/point2point.h>
namespace qmcplusplus
{

/**
 * Each SplineAdoptor needs a reader derived from BsplineReaderBase.
 * This base class handles common chores
 * - check_twists : read gvectors, set twists for folded bands if needed, and set the phase for the special K
 * - set_grid : create the basic grid and boundary conditions for einspline
 * Note that template is abused but it works.
 */
struct BsplineReaderBase
{
  ///pointer to the EinsplineSetBuilder
  EinsplineSetBuilder* mybuilder;
  ///communicator
  Communicate* myComm;
  ///mesh size
  TinyVector<int,3> MeshSize;
  ///first index of the SPO set
  int myFirstSPO;
  ///number of orbitals to be created
  int myNumOrbs;

  BsplineReaderBase(EinsplineSetBuilder* e)
    : mybuilder(e), MeshSize(0), myFirstSPO(0),myNumOrbs(0)
  {
    myComm=mybuilder->getCommunicator();
  }

  ///** copy minimal informatino from EinsplineSet to manage SplineAdoptor
  // */
  //template<typename SPE>
  //void init(EinsplineSet* in, SPE* out)
  //{
  //  out->PrimLattice=in->PrimLattice;
  //  out->SuperLattice=in->SuperLattice;
  //  out->GGt=in->GGt;
  //  out->setOrbitalSetSize(in->getOrbitalSetSize());
  //}

  /** read gvectors and set the mesh, and prepare for einspline
   */
  template<typename GT, typename BCT>
  inline bool set_grid(const TinyVector<int,3>& halfg, GT* xyz_grid, BCT* xyz_bc)
  {
    //This sets MeshSize from the input file
    bool havePsig=mybuilder->ReadGvectors_ESHDF();

    //If this MeshSize is not initialized, use the meshsize set by the input based on FFT grid and meshfactor
    if(MeshSize[0]==0) MeshSize=mybuilder->MeshSize;

    app_log() << "  Using meshsize=" << MeshSize 
            << "\n  vs input meshsize=" << mybuilder->MeshSize << endl;

    xyz_grid[0].start = 0.0;
    xyz_grid[0].end = 1.0;
    xyz_grid[0].num = MeshSize[0];
    xyz_grid[1].start = 0.0;
    xyz_grid[1].end = 1.0;
    xyz_grid[1].num = MeshSize[1];
    xyz_grid[2].start = 0.0;
    xyz_grid[2].end = 1.0;
    xyz_grid[2].num = MeshSize[2];
    for(int j=0; j<3; ++j)
    {
      if(halfg[j])
      {
        xyz_bc[j].lCode=ANTIPERIODIC;
        xyz_bc[j].rCode=ANTIPERIODIC;
      }
      else
      {
        xyz_bc[j].lCode=PERIODIC;
        xyz_bc[j].rCode=PERIODIC;
      }
    }
    return havePsig;
  }

  /** initialize twist-related data for N orbitals
   */
  template<typename SPE>
  inline void check_twists(EinsplineSet* orbitalSet, SPE* bspline, const BandInfoGroup& bandgroup)
  {
    //init(orbitalSet,bspline);
    bspline->PrimLattice=orbitalSet->PrimLattice;
    bspline->SuperLattice=orbitalSet->SuperLattice;
    bspline->GGt=orbitalSet->GGt;

    int N = bandgroup.getNumDistinctOrbitals();
    int numOrbs=bandgroup.getNumSPOs();

    bspline->setOrbitalSetSize(numOrbs);
    bspline->resizeStorage(N,N);

    bspline->first_spo=bandgroup.getFirstSPO();
    bspline->last_spo =bandgroup.getLastSPO();

    int num = 0;
    const vector<BandInfo>& cur_bands=bandgroup.myBands;
    for (int iorb=0; iorb<N; iorb++)
    {
      int ti = cur_bands[iorb].TwistIndex;
      bspline->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(mybuilder->TwistAngles[ti]); //twist);
      bspline->MakeTwoCopies[iorb] = (num < (numOrbs-1)) && cur_bands[iorb].MakeTwoCopies;
      num += bspline->MakeTwoCopies[iorb] ? 2 : 1;
    }

    app_log() << "NumDistinctOrbitals " << N
              << " numOrbs = " << numOrbs << endl;

    bspline->HalfG=0;
    TinyVector<int,3> bconds=mybuilder->TargetPtcl.Lattice.BoxBConds;
    if(!bspline->is_complex)
    {
      //no k-point folding, single special k point (G, L ...)
      TinyVector<double,3> twist0 = mybuilder->TwistAngles[cur_bands[0].TwistIndex];
      for (int i=0; i<3; i++)
        if (bconds[i] && ((std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)))
          bspline->HalfG[i] = 1;
        else
          bspline->HalfG[i] = 0;
      app_log() << "  TwistIndex = " << cur_bands[0].TwistIndex << " TwistAngle " << twist0 << endl;
      app_log() <<"   HalfG = " << bspline->HalfG << endl;
    }
    app_log().flush();
  }

  /** return the path name in hdf5
   */
  inline string psi_g_path(int ti, int spin, int ib)
  {
    ostringstream path;
    path << "/electrons/kpoint_" << ti
         << "/spin_" << spin << "/state_" << ib << "/psi_g";
    return path.str();
  }

  /** return the path name in hdf5
   */
  inline string psi_r_path(int ti, int spin, int ib)
  {
    ostringstream path;
    path << "/electrons/kpoint_" << ti
         << "/spin_" << spin << "/state_" << ib << "/psi_r";
    return path.str();
  }

  /** read/bcast psi_g
   * @param ti twist index
   * @param spin spin index
   * @param ib band index
   * @param cG psi_g as stored in hdf5
   */
  void get_psi_g(int ti, int spin, int ib, Vector<complex<double> >& cG);

  virtual ~BsplineReaderBase() {}

  /** create the actual spline sets
   */
  virtual SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalset, const BandInfoGroup& bandgroup)=0;

  /** create the spline set */
  SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalset);
};

}
#endif
