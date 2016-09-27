//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_H
#define QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <QMCWaveFunctions/CompositeSPOSet.h>
#include <ParticleBase/RandomSeqGenerator.h>

namespace qmcplusplus
{

class DensityMatrices1B : public QMCHamiltonianBase
{

 public:

  enum{DIM=OHMMS_DIM};

  typedef ValueType Value_t;
  typedef GradType  Grad_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef ParticleSet::ParticleLayout_t Lattice_t;
  typedef Vector<Value_t> Vector_t;
  typedef Matrix<Value_t> Matrix_t;
  typedef std::vector<PosType>  pts_t;
  typedef std::vector<RealType> dens_t;

  enum integrators {uniform_grid=0,uniform,density,no_integrator};

  enum evaluators {loop=0,matrix,no_evaluator};

  enum samplings {volume_based=0,metropolis,no_sampling};

  //data members
  bool energy_mat;
  CompositeSPOSet basis_functions;
  ValueVector_t   basis_values;
  ValueVector_t   basis_norms;
  GradVector_t    basis_gradients;
  ValueVector_t   basis_laplacians;
  ValueVector_t   integrated_values;
  bool warmed_up;
  std::vector<PosType> rsamples;
  Vector<RealType> sample_weights;
  RealType dens;
  PosType  drift;
  int nindex;
  int eindex;
  TrialWaveFunction& Psi;
  Lattice_t& Lattice;
  ParticleSet& Pq;
  const ParticleSet* Pc;
  TraceSample<TraceReal>* w_trace;
  TraceSample<TraceComp>* T_trace;
  CombinedTraceSample<TraceReal>* Vq_trace;
  CombinedTraceSample<TraceReal>* Vc_trace;
  CombinedTraceSample<TraceReal>* Vqq_trace;
  CombinedTraceSample<TraceReal>* Vqc_trace;
  CombinedTraceSample<TraceReal>* Vcc_trace;
  std::vector<Value_t> E_samp;
  CombinedTraceSample<TraceReal>* E_trace;


  bool initialized;
  bool normalized;
  int basis_size;
  int samples;
  int nparticles;
  int nspecies;
  std::vector<int> species_size;
  std::vector<std::string> species_name;
  std::vector<Vector_t*> E_N;
  std::vector<Matrix_t*> Phi_NB,Psi_NM,Phi_Psi_NB,N_BB,E_BB;
  Matrix_t Phi_MB;
  bool check_overlap;
  bool check_derivatives;

//#define DMCHECK
#ifdef DMCHECK
  std::vector<Vector_t*> E_Ntmp;
  std::vector<Matrix_t*> Phi_NBtmp,Psi_NMtmp,Phi_Psi_NBtmp,N_BBtmp,E_BBtmp;
  Matrix_t Phi_MBtmp;
#endif

  integrators integrator;
  evaluators  evaluator;
  samplings   sampling;
  int         points;
  RealType    scale;
  PosType     center,rcorner;
  RealType    volume;
  bool        periodic;
  int         warmup;
  RealType    timestep;
  bool        use_drift;
  int         nmoves;
  int         naccepted;
  RealType    acceptance_ratio;
  bool write_acceptance_ratio;
  bool write_rstats;

  int ind_dims[DIM];
  RealType metric;

  PosType  rpcur;
  PosType  dpcur;
  RealType rhocur;

  RandomGenerator_t* uniform_random;

  //constructor/destructor
  DensityMatrices1B(ParticleSet& P, TrialWaveFunction& psi,ParticleSet* Pcl);
  DensityMatrices1B(DensityMatrices1B& master, ParticleSet& P, TrialWaveFunction& psi);
  ~DensityMatrices1B();

  //standard interface
  QMCHamiltonianBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur);
  Return_t evaluate(ParticleSet& P);
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P); 
  }

  //optional standard interface
  void get_required_traces(TraceManager& tm);
  void setRandomGenerator(RandomGenerator_t* rng);

  //required for Collectables interface
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P)                      { }
  void setObservables(PropertySetType& plist)                      { }
  void setParticlePropertyList(PropertySetType& plist, int offset) { }
  void contribute_scalar_quantities()                              { }
  void checkout_scalar_quantities(TraceManager& tm)                { }
  void collect_scalar_quantities()                                 { }
  void delete_scalar_quantities()                                  { }

  //obsolete?
  bool get(std::ostream& os) const { return false; }

  //local functions
  //  initialization/finalization
  void reset();
  void set_state(xmlNodePtr cur);
  void set_state(DensityMatrices1B& master);
  void initialize();
  void finalize();
  void normalize();
  //  printing
  void report(const std::string& pad="");
  //  sample generation
  void warmup_sampling();
  void generate_samples(RealType weight,int steps=0);
  void generate_uniform_grid(RandomGenerator_t& rng);
  void generate_uniform_samples(RandomGenerator_t& rng);
  void generate_density_samples(bool save,int steps,RandomGenerator_t& rng);
  void diffusion(RealType sqt,PosType& diff);
  void density_only(const PosType& r,RealType& dens);
  void density_drift(const PosType& r,RealType& dens,PosType& drift);
  //  basis & wavefunction ratio matrix construction
  void get_energies(std::vector<Vector_t*>& E_n);
  void generate_sample_basis(Matrix_t& Phi_mb);
  void generate_sample_ratios(std::vector<Matrix_t*> Psi_nm);
  void generate_particle_basis(ParticleSet& P,std::vector<Matrix_t*>& Phi_nb);
  //  basis set updates
  void update_basis(const PosType& r);
  void update_basis_d012(const PosType& r);
  //  testing
  void test_overlap();
  void test_derivatives();
  //  original loop implementation
  void integrate(ParticleSet& P,int n);
  Return_t evaluate_loop(ParticleSet& P);
  //  matrix implementation
  Return_t evaluate_check(ParticleSet& P);
  Return_t evaluate_matrix(ParticleSet& P);
  //  postprocessing
  void postprocess_density(const std::string& infile,const std::string& species,
                           pts_t& points,dens_t& density,dens_t& density_err);


  bool match(Value_t e1, Value_t e2, RealType tol=1e-12);
  bool same(Vector_t& v1, Vector_t& v2, RealType tol=1e-6);
  bool same(Matrix_t& m1, Matrix_t& m2, RealType tol=1e-6);
  void compare(const std::string& name, Vector_t& v1, Vector_t& v2,
               bool write=false,bool diff_only=true);
  void compare(const std::string& name, Matrix_t& m1, Matrix_t& m2,
                      bool write=false,bool diff_only=true);

  
};

}

#endif
