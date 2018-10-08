//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPLINEBUILDER_H
#define QMCPLUSPLUS_SPLINEBUILDER_H

#include <map>

#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/SPOSetInfo.h"
#include "QMCWaveFunctions/SPOSetInputInfo.h"
#include "QMCWaveFunctions/SPOSet.h"

#include "QMCWaveFunctions/BandInfo.h"

#include "Configuration.h"

namespace qmcplusplus
{


  
class SplineBuilder
{
public:
  struct CentersInfo
  {
    std::vector<int> lmax, spline_npoints, GroupID;
    std::vector<double> spline_radius, cutoff, inner_cutoff, non_overlapping_radius;
    std::vector<TinyVector<double,OHMMS_DIM> > ion_pos;
    int Ncenters;

    CentersInfo(): Ncenters(0) {};

    void resize(int ncenters)
    {
      Ncenters=ncenters;
      lmax.resize(ncenters, -1);
      spline_npoints.resize(ncenters, -1);
      GroupID.resize(ncenters, 0);
      spline_radius.resize(ncenters, -1.0);
      inner_cutoff.resize(ncenters, -1.0);
      non_overlapping_radius.resize(ncenters, -1.0);
      cutoff.resize(ncenters, -1.0);
      ion_pos.resize(ncenters);
    }
  };
public:

  //typedef std::map<std::string,ParticleSet*> PtclPoolType;
  using UnitCellType = CrystalLattice<ParticleSet::Scalar_t,OHMMS_DIM>;

  /** initialize the Antisymmetric wave function for electrons
   * @param cur the current xml node
   */
  virtual SPOSet<>* createSPOSetFromXML(xmlNodePtr cur) = 0;

  /** a specific but clean code path in createSPOSetFromXML, for PBC, double, ESHDF
   * @param cur the current xml node
   */
  virtual void set_metadata(int numOrbs, int TwistNum_inp) = 0;

  /** initialize with the existing SPOSet */
  virtual SPOSet<>* createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info) = 0;

  //////////////////////////////////////
  // HDF5-related data  and functions //
  //////////////////////////////////////
  virtual bool ReadOrbitalInfo () = 0;
  virtual bool ReadOrbitalInfo_ESHDF () = 0;
  virtual void BroadcastOrbitalInfo() = 0;
  virtual bool CheckLattice() = 0;

  /** read gvectors for each twist
   * @return true, if psi_g is found
   */
  virtual bool ReadGvectors_ESHDF() = 0;

  virtual bool TwistPair (QMCTraits::PosType a, QMCTraits::PosType b) = 0;
  virtual void AnalyzeTwists2() = 0;
  virtual void TileIons() = 0;
  virtual void OccupyBands(int spin, int sortBands, int numOrbs) = 0;
  virtual void OccupyBands_ESHDF(int spin, int sortBands, int numOrbs) = 0;

  // virtual void CopyBands(int numOrbs) = 0;

  // This returns the path in the HDF5 file to the group for orbital
  // with twist ti and band bi
  virtual std::string OrbitalPath   (int ti, int bi) = 0;
  virtual std::string CoreStatePath (int ti, int bi) = 0;
  virtual std::string MuffinTinPath (int ti, int bi, int tin) = 0;

  virtual bool bcastSortBands(int splin, int N, bool root) = 0;

  //Accessors to implementation variables
  virtual TinyVector<int, 3> getMeshSize() = 0;
  virtual UnitCellType getPrimCell() = 0;
  virtual UnitCellType getSuperCell() = 0;
  virtual std::vector<int> getSuper2Prim() = 0;
  virtual Tensor<double, OHMMS_DIM> getGGt() = 0;
  virtual std::vector<TinyVector<double,OHMMS_DIM> >& getTwistAngles() = 0;
  virtual ParticleSet& getTargetPtcl() = 0;
  virtual ParticleSet& getSourcePtcl() = 0;
  virtual CentersInfo& getAtomicCentersInfo() = 0;
  virtual int getTwistNum() = 0;
  virtual std::vector<std::vector<TinyVector<int,3> > >& getGvecs() = 0;
  virtual std::string getH5FileName() = 0;
  virtual hid_t getH5FileID() = 0;
  virtual int getNumberSpinStates() = 0;
  virtual std::vector<BandInfo>& getFullBandsBySpin(int) = 0;
  virtual int getMaxNumGvecs() = 0;
  virtual std::vector<SPOSetInfo*>& getStates() = 0;
  virtual Tensor<int,OHMMS_DIM>& getTileMatrix() = 0;
  //inline void update_token(const char* f, int l, const char* msg) 
  //{}
  virtual std::string getName() = 0;
  // Why are we getting our mpi communicator from the spline set builder?
  virtual Communicate* getCommunicator() = 0;

};

}
#endif
