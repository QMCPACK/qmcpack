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

#ifndef QMCPLUSPLUS_MOCKEINSPLINESETBUILDER_H
#define QMCPLUSPLUS_MOCKEINSPLINESETBUILDER_H

#include "QMCWaveFunctions/SplineBuilder.h"

namespace qmcplusplus
{
  
class MockEinsplineSetBuilder : public SplineBuilder
{
public:
  MockEinsplineSetBuilder()
  {
    OHMMS::Controller->initialize(0, NULL);
    myComm = OHMMS::Controller;
  };
  
  MockEinsplineSetBuilder(const MockEinsplineSetBuilder&) = default;

  /** initialize the Antisymmetric wave function for electrons
   * @param cur the current xml node
   */
  SPOSet* createSPOSetFromXML(xmlNodePtr cur) { return nullptr; };

  /** a specific but clean code path in createSPOSetFromXML, for PBC, double, ESHDF
   * @param cur the current xml node
   */
  void set_metadata(int numOrbs, int TwistNum_inp) { return; };

  /** initialize with the existing SPOSet */
  SPOSet* createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info) { return nullptr; };

  Communicate* myComm;
  //////////////////////////////////////
  // HDF5-related data  and functions //
  //////////////////////////////////////
  bool ReadOrbitalInfo () { return false; };
  bool ReadOrbitalInfo_ESHDF () { return false; };
  void BroadcastOrbitalInfo() { return; };
  bool CheckLattice() { return false; };

  /** read gvectors for each twist
   * @return true, if psi_g is found
   */
  bool ReadGvectors_ESHDF() { return false; };

  bool TwistPair (QMCTraits::PosType a, QMCTraits::PosType b) { return false; };
  void AnalyzeTwists2() { return; };
  void TileIons() { return; };
  void OccupyBands(int spin, int sortBands, int numOrbs) { return; };
  void OccupyBands_ESHDF(int spin, int sortBands, int numOrbs) { return; };

  // virtual void CopyBands(int numOrbs) = 0;

  // This returns the path in the HDF5 file to the group for orbital
  // with twist ti and band bi
  std::string OrbitalPath   (int ti, int bi) { return "mock_String"; };
  std::string CoreStatePath (int ti, int bi) { return "mock_String"; };
  std::string MuffinTinPath (int ti, int bi, int tin) { return "mock_String"; };

  bool bcastSortBands(int splin, int N, bool root) { return false; };

  std::vector<TinyVector<double,OHMMS_DIM> > TwistAngles;
  ParticleSet TargetPtcl;
  ParticleSet SourcePtcl;
  CentersInfo AtomicCentersInfo;
  std::vector<std::vector<TinyVector<int,3> > > Gvecs;
  std::vector<std::vector<BandInfo>*> FullBands;
  std::vector<SPOSetInfo*> states;
  Tensor<int,OHMMS_DIM> TileMatrix;
  
  //Accessors to implementation variables
  TinyVector<int, 3> getMeshSize() { return TinyVector<int, 3>(); };
  UnitCellType getPrimCell() { return UnitCellType(); };
  UnitCellType getSuperCell() { return UnitCellType(); };
  std::vector<int> getSuper2Prim() { return std::vector<int>(); };
  Tensor<double, OHMMS_DIM> getGGt() { return Tensor<double, OHMMS_DIM>(); };
  std::vector<TinyVector<double,OHMMS_DIM> >& getTwistAngles() { return TwistAngles; };
  ParticleSet& getTargetPtcl() { return TargetPtcl; };
  ParticleSet& getSourcePtcl() { return SourcePtcl; };
  CentersInfo& getAtomicCentersInfo() { return AtomicCentersInfo; };
  int getTwistNum() { return int(); };
  std::vector<std::vector<TinyVector<int,3> > >& getGvecs() { return Gvecs; };
  std::string getH5FileName() { return "mock_file_name"; };
  hid_t getH5FileID() { return 0; };
  int getNumberSpinStates() { return 1; };
  std::vector<BandInfo>& getFullBandsBySpin(int spin) { return *(FullBands[spin]); };
  int getMaxNumGvecs() { return 0; };
  std::vector<SPOSetInfo*>& getStates() { return states; };
  Tensor<int,OHMMS_DIM>& getTileMatrix() { return TileMatrix; };
  
  //inline void update_token(const char* f, int l, const char* msg) 
  //{}
  std::string getName() { return "mock_String"; };
  // Why are we getting our mpi communicator from the spline set builder?
  Communicate* getCommunicator() { return myComm; };

};


}
#endif //QMCPLUSPLUS_MOCKEINSPLINESETBUILDER_H
