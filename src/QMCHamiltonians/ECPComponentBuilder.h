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


/** @file
 * @brief Declaration of a builder class for an ECP component for an ionic type
 */
#ifndef QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#define QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "QMCHamiltonians/SOECPComponent.h"
#include "QMCHamiltonians/L2Potential.h"

namespace qmcplusplus
{
struct ECPComponentBuilder : public MPIObjectBase, public QMCTraits
{
  typedef LocalECPotential::GridType GridType;
  typedef ParticleSet::Scalar_t mRealType;
  typedef OneDimGridBase<mRealType> mGridType;
  typedef LocalECPotential::RadialPotentialType RadialPotentialType;

  int NumNonLocal;
  int Lmax, Llocal, Nrule, Srule;
  int NumSO;  //The number of spin-orbit channels.
  int LmaxSO; //The maximum angular momentum of spin-orbit channels.
  int AtomicNumber;
  RealType Zeff;
  RealType RcutMax;
  std::string Species;
  std::unique_ptr<mGridType> grid_global;
  std::map<std::string, std::unique_ptr<mGridType>> grid_inp;
  std::unique_ptr<RadialPotentialType> pp_loc;
  std::unique_ptr<NonLocalECPComponent> pp_nonloc;
  std::unique_ptr<SOECPComponent> pp_so; //Spin-orbit potential component.
  std::unique_ptr<L2RadialPotential> pp_L2;
  std::map<std::string, int> angMon;

  ECPComponentBuilder(const std::string& aname, Communicate* c, int nrule = -1);

  bool parse(const std::string& fname, xmlNodePtr cur);
  bool put(xmlNodePtr cur);
  void addSemiLocal(xmlNodePtr cur);
  void buildLocal(xmlNodePtr cur);
  void buildSemiLocalAndLocal(std::vector<xmlNodePtr>& semiPtr);
  void buildL2(xmlNodePtr cur);

  bool parseCasino(const std::string& fname, xmlNodePtr cur); //std::string& fname, RealType rc);
  //bool parseCasino(std::string& fname, RealType rc);
  // This sets the spherical quadrature rule used to apply the
  // projection operators.  rule can be 1 to 7.  See
  // J. Chem. Phys. 95 (3467) (1991)
  // Rule     # points     lexact
  //  1           1          0
  //  2           4          2
  //  3           6          3
  //  4          12          5
  //  5          18          5
  //  6          26          7
  //  7          50         11
  void SetQuadratureRule(int rule);

  std::unique_ptr<mGridType> createGrid(xmlNodePtr cur, bool useLinear = false);
  RadialPotentialType* createVrWithBasisGroup(xmlNodePtr cur, mGridType* agrid);
  RadialPotentialType* createVrWithData(xmlNodePtr cur, mGridType* agrid, int rCorrection = 0);

  void doBreakUp(const std::vector<int>& angList,
                 const Matrix<mRealType>& vnn,
                 RealType rmax,
                 mRealType Vprefactor = 1.0);

  /** brief buildSO - takes the previously parsed angular momenta and spin-orbit tabulated potentials and uses
  **     them to construct SOECPComponent* pp_so.  This is called in "doBreakUp". 
  **
  ** param std::vector<int>& angList  The angular momentum for each SO potential.
  ** param Matrix<mRealType>& vnnso (npot x ngrid) matrix storing all tabulated SO potentials.
  ** param RealType rmax  max r on the specified grid.
  ** param mRealType Vprefactor  optional scale factor.
  **
  ** return void
  **
  **/
  void buildSO(const std::vector<int>& angList,
               const Matrix<mRealType>& vnnso,
               RealType rmax,
               mRealType Vprefactor = 1.0);

  void printECPTable();
  bool read_pp_file(const std::string& fname);
};

// Read a file into a memory buffer.
// Under MPI, it reads the file with one node and broadcasts the contents to all the other nodes.

class ReadFileBuffer
{
  char* cbuffer;
  std::ifstream* fin;
  Communicate* myComm;
  int get_file_length(std::ifstream* f) const;

public:
  bool is_open;
  int length;
  ReadFileBuffer(Communicate* c) : cbuffer(NULL), fin(NULL), myComm(c), is_open(false), length(0) {}
  bool open_file(const std::string& fname);
  bool read_contents();
  char* contents() { return cbuffer; }
  void reset();

  ~ReadFileBuffer()
  {
    delete[] cbuffer;
    if (fin)
      delete fin;
  }
};


} // namespace qmcplusplus
#endif
