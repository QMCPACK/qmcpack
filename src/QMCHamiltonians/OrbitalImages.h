//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ORBITAL_IMAGES_H
#define QMCPLUSPLUS_ORBITAL_IMAGES_H

#include "QMCHamiltonians/OperatorBase.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"

namespace qmcplusplus
{
/** "Estimator" to produce files for orbital plotting.
 *  Orbitals are evaluated on a uniform grid and written to ascii files.
 *  Only format currently supported is xsf (XCrySDen).
 *  Can print real, imag, abs, and abs^2 of each orbital.
 *  This class should work with any SPOSet.
 *  All work is performed by omp thread 0 of mpi task 0. 
 *
 *  TO USE THIS YOU NEED TO KNOW THE NAME OF THE SPOSET(S)!!!
 *    For example, using sposet_builder the names are prominently displayed:
 * 
 *      <sposet_builder type="bspline" href="pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.10">
 *         <sposet type="bspline" name="spo_ud" size="2" spindataset="0"/>
 *      </sposet_builder>
 * 
 *    In this case a single sposet named "spo_ud" exists.
 *    
 *    If you are using Slater-Jastrow w/o sposet_builder 
 *    the sposets should be named updet and downdet.   
 * 
 * 
 *  To make xsf files, add xml similar to the following to <hamiltonian/>:
 *
 *    minimal working example: single sposet named spo_ud on a 20x20x20 grid
 *
 *      <estimator name="OrbitalImages" type="orbitalimages" ions="ion0">
 *        <parameter name="sposets"> spo_ud </parameter>
 *        <parameter name="grid"> 20 20 20 </parameter>
 *      </estimator>
 *      
 *    as above, but print real, imag, abs, and abs^2 of orbitals (default real)  
 *      
 *      <estimator name="OrbitalImages" type="orbitalimages" ions="ion0">
 *        <parameter name="sposets"> spo_ud </parameter>
 *        <parameter name="grid"> 20 20 20 </parameter>
 *        <parameter name="value"> real imag abs abs2 </parameter>
 *      </estimator>
 *      
 *    up and down sposets named spo_u and spo_d w/ individual orbitals selected  
 *      
 *      <estimator name="OrbitalImages" type="orbitalimages" ions="ion0">
 *        <parameter name="sposets"> spo_u spo_d </parameter>
 *        <parameter name="spo_u">   13 24 37    </parameter>
 *        <parameter name="spo_d">   10 18 29 41 </parameter>
 *        <parameter name="grid">    20 20 20    </parameter>
 *        <parameter name="value"> real imag abs abs2 </parameter>
 *      </estimator>
 *
 *    user defined cell by cell center and axes (openbc's, subset of pbc cell)
 *
 *      <estimator name="OrbitalImages" type="orbitalimages" ions="ion0">
 *        <parameter name="sposets"> spo_ud </parameter>
 *        <parameter name="grid"> 20 20 20 </parameter>
 *        <parameter name="center"> 2.5 2.5 2.5 </parameter>
 *        <parameter name="cell"> 
 *           5.0  0.0  0.0 
 *           0.0  5.0  0.0 
 *           0.0  0.0  5.0 
 *        </parameter>
 *      </estimator>
 *
 *    user defined cell by cell corner and axes (default corner is 0 0 0)
 *
 *      <estimator name="OrbitalImages" type="orbitalimages" ions="ion0">
 *        <parameter name="sposets"> spo_ud </parameter>
 *        <parameter name="grid"> 20 20 20 </parameter>
 *        <parameter name="corner"> 0.0 0.0 0.0 </parameter>
 *        <parameter name="cell"> 
 *           5.0  0.0  0.0 
 *           0.0  5.0  0.0 
 *           0.0  0.0  5.0 
 *        </parameter>
 *      </estimator>
 *
 *    store only a batch of orbital values in memory
 *    this can save on memory for very large grids
 *    but will evaluate all orbitals #orbitals/batch_size times
 *
 *      <estimator name="OrbitalImages" type="orbitalimages" ions="ion0">
 *        <parameter name="sposets"> spo_ud </parameter>
 *        <parameter name="grid"> 200 200 200 </parameter>
 *        <parameter name="batch_size"> 10    </parameter>
 *      </estimator>
 *
 */
class OrbitalImages : public OperatorBase
{
public:
  enum
  {
    DIM = OHMMS_DIM
  };

  using ValueVector = SPOSet::ValueVector;
  using GradVector  = SPOSet::GradVector;
  using Lattice_t   = ParticleSet::ParticleLayout;
  using PSPool      = std::map<std::string, ParticleSet*>;


  ///derivative types
  enum derivative_types_enum
  {
    value_d = 0,
    gradient_d,
    laplacian_d
  };


  ///options for orbital value to write
  enum value_types_enum
  {
    real_val = 0,
    imag_val,
    abs_val,
    abs2_val
  };

  ///options for orbital output file format
  enum formats_enum
  {
    xsf = 0
  };

  ///at put() ion particleset is obtained from ParticleSetPool
  PSPool& psetpool;

  ///electron particleset
  ParticleSet* Peln;

  ///ion particleset
  ParticleSet* Pion;

  ///mpi communicator
  Communicate* comm;


  ///file format selection
  formats_enum format;

  ///orbital value selections
  std::vector<value_types_enum> value_types;

  ///write out derivatives in addition to values
  bool derivatives;

  ///names of sposets to evaluate
  std::vector<std::string> sposet_names;

  ///indices of orbitals within each sposet to evaluate
  const std::shared_ptr<std::vector<std::vector<int>>> sposet_indices;

  ///sposets obtained by name from WaveFunctionFactory
  std::vector<std::unique_ptr<SPOSet>> sposets;

  ///evaluate points at grid cell centers instead of edges
  bool center_grid;

  ///cell bounding the evaluation grid, default is simulation cell
  Lattice_t cell;

  ///location of cell corner, positions in the cell are corner+uvec*cell
  PosType corner;

  ///number of grid points in each direction (cell axis)
  TinyVector<int, DIM> grid;

  ///stride to generate grid in arbitrary dimensions
  TinyVector<int, DIM> gdims;

  ///total number of grid points
  int npoints;

  ///number of orbitals to store in memory at a time (batch_size*npoints)
  int batch_size;

  ///temporary vector to hold values of all orbitals at a single point
  ValueVector spo_vtmp;

  ///temporary vector to hold gradients of all orbitals at a single point
  GradVector spo_gtmp;

  ///temporary vector to hold laplacians of all orbitals at a single point
  ValueVector spo_ltmp;

  ///temporary array to hold values of a batch of orbitals at all grid points
  Array<ValueType, 2> batch_values;

  ///temporary array to hold gradients of a batch of orbitals at all grid points
  Array<GradType, 2> batch_gradients;

  ///temporary array to hold laplacians of a batch of orbitals at all grid points
  Array<ValueType, 2> batch_laplacians;

  ///temporary array to hold values of a single orbital at all grid points
  std::vector<ValueType> orbital;

  //constructors
  OrbitalImages(ParticleSet& P, PSPool& PSP, Communicate* mpicomm, const WaveFunctionFactory& factory);
  OrbitalImages(const OrbitalImages& other);

  //standard interface
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) final;

  ///read xml input
  bool put(xmlNodePtr cur) override;

  ///hijack estimator evaluate to evaluate and write all orbitals
  Return_t evaluate(ParticleSet& P) override;

  //optional standard interface
  //void getRequiredTraces(TraceManager& tm);
  //void setRandomGenerator(RandomGenerator* rng);

  //required for Collectables interface
  void addObservables(PropertySetType& plist, BufferType& olist) override {}
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const override {}

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P) override {}
  void setObservables(PropertySetType& plist) override {}
  void setParticlePropertyList(PropertySetType& plist, int offset) override {}
#if !defined(REMOVE_TRACEMANAGER)
  void checkout_scalar_arrays(TraceManager& tm) {}
  void collect_scalar_samples() {}
  void delete_scalar_arrays() {}
#endif

  //obsolete?
  bool get(std::ostream& os) const override { return false; }

  //local functions
  ///write brief report of configuration data
  void report(const std::string& pad);

  ///write a single orbital to file
  void write_orbital(const std::string& sponame,
                     int index,
                     std::vector<ValueType>& orb,
                     value_types_enum value_type,
                     derivative_types_enum derivative_type = value_d,
                     int dimension                         = 0);

  ///write a single orbital to an xsf file
  void write_orbital_xsf(const std::string& sponame,
                         int index,
                         std::vector<ValueType>& orb,
                         value_types_enum value_type,
                         derivative_types_enum derivative_type = value_d,
                         int dimension                         = 0);

private:
  /// reference to the sposet_builder_factory
  const WaveFunctionFactory& wf_factory_;
};

} // namespace qmcplusplus

#endif
