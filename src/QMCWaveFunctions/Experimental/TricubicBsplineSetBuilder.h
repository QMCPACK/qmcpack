//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_TRICUBIC_BSPLINESETBUILDER_H
#define QMCPLUSPLUS_TRICUBIC_BSPLINESETBUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/Bspline3DSetBase.h"
#include "Numerics/HDFNumericAttrib.h"

namespace qmcplusplus
{

class PWParameterSet;

/**@ingroup WFSBuilder
 * A builder class for a set of Spline functions
 */
class TricubicBsplineSetBuilder: public BasisSetBuilder
{

public:

  typedef Bspline3DSetBase              BsplineBasisType;
  typedef Bspline3DSetBase::StorageType StorageType;
  typedef std::map<std::string,ParticleSet*>      PtclPoolType;

  /** real-space orbital type
   */
  struct RSOType
  {
    PosType Center;
    PosType Origin;
    StorageType* Coeffs;
    RSOType():Coeffs(0) {}
    RSOType(const PosType& c, const PosType o, StorageType* d): Center(c),Origin(o),Coeffs(d) {}
  };

  /** constructor
   * @param p target ParticleSet
   * @param psets a set of ParticleSet objects
   */
  TricubicBsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur);

  ///destructor
  ~TricubicBsplineSetBuilder();

  ///implement put function to process an xml node
  bool put(xmlNodePtr cur);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   */
  SPOSetBase* createSPOSet(xmlNodePtr cur);

private:

  /** set of StorageType*
   *
   * Key is $1#$2 where $1 is the hdf5 file name and $2 is the band indenx
   */
  //static std::map<std::string,StorageType*> BigDataSet;
  static std::map<std::string,RSOType*> BigDataSet;
  ///if true, grid is open-ended [0,nx) x [0,ny) x [0, nz)
  bool OpenEndGrid;
  ///if true, a grid is translated so that the center of localized orbitals conincides with the cell center
  bool TranslateGrid;
  ///if true, the input grid is not fixed.
  bool FloatingGrid;
  ///boolean to enable debug with EG
  bool DebugWithEG;
  ///number of SPOs created by this buidlder
  int CurSPOSize;
  ///twist angle
  PosType TwistAngle;
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///reference to a ParticleSetPool
  PtclPoolType& ptclPool;
  ///base lattice for the numerical grids
  ParticleSet::ParticleLayout_t basisLattice;
  ///number of grid points for each direction
  TinyVector<IndexType,DIM> BoxGrid;
  ///number of copies of the lattice for tiling
  TinyVector<IndexType,DIM> BoxDup;
  ///coordiate of the box origin
  PosType LowerBox;
  ///coordiate of the box bound
  PosType UpperBox;
  ///parameter set for h5 tags
  PWParameterSet* myParam;
  ///save xml node
  xmlNodePtr rootNode;
  ///hdf5 handler to clean up
  hid_t hfileID;
  ///pointer to the BsplineBasisType use
  BsplineBasisType* activeBasis;
  ///grid of the orinal grid
  Bspline3DSetBase::GridType dataKnot;
  ///current hdf5 file name
  std::string curH5Fname;
  ///set of WFSetType*
  std::map<std::string,BsplineBasisType*> myBasis;
  ///single-particle orbital sets
  std::map<std::string,SPOSetBase*> mySPOSet;

  ///set Bspline functions
  void setBsplineBasisSet(xmlNodePtr cur);

  ///initialize the grid of the spline orbitals
  void initGrid();

  ///initialize the center and origin of orbitals
  void getCenterAndOrigin(const char* hroot, const std::vector<int>& occSet, int norb);

  ///read numerical orbitals and spline them
  template<typename Tin, typename Tout>
  void readData(const char* hroot, const std::vector<int>& occSet, int spinIndex, int degeneracy);

  ///read numerical orbitals and spline them in parallel
  template<typename Tin, typename Tout>
  void readDataOMP(const char* hroot, const std::vector<int>& occSet, int spinIndex, int degeneracy);

  void readComplex2RealDataWithTruncation(const char* hroot, const std::vector<int>& occSet,
                                          int spinIndex, int degeneracy);

  ///a function to test with EG
  //SPOSetBase* createSPOSetWithEG();
};
}
#endif
