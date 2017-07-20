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
    
    

#include <QMCHamiltonians/DensityMatrices1B.h>
#include <OhmmsData/AttributeSet.h>
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <Numerics/MatrixOperators.h>
#include <Utilities/IteratorUtility.h>
#include <Utilities/string_utils.h>
#include <QMCWaveFunctions/BasisSetFactory.h>



namespace qmcplusplus
{


  using MatrixOperators::product;
  using MatrixOperators::product_AtB;
  using MatrixOperators::diag_product;


  DensityMatrices1B::DensityMatrices1B(ParticleSet& P, TrialWaveFunction& psi, ParticleSet* Pcl) 
    : Lattice(P.Lattice),Pq(P),Psi(psi),Pc(Pcl)
  {
    reset();
  }


  DensityMatrices1B::DensityMatrices1B(DensityMatrices1B& master, ParticleSet& P, TrialWaveFunction& psi) 
    : QMCHamiltonianBase(master), Lattice(P.Lattice),Pq(P),Psi(psi),Pc(master.Pc)
  {
    app_log()<<"dm1b deepcopy"<< std::endl;
    reset();
    set_state(master);
    basis_functions.clone_from(master.basis_functions);
    initialize();
    app_log()<<"dm1b end deepcopy"<< std::endl;
  }


  DensityMatrices1B::~DensityMatrices1B()
  {
    if(initialized)
      finalize();
  }


  QMCHamiltonianBase* DensityMatrices1B::makeClone(ParticleSet& P, TrialWaveFunction& psi)
  {
    return new DensityMatrices1B(*this,P,psi);
  }


  void DensityMatrices1B::reset()
  {
    // uninitialized data
    w_trace        = NULL;
    T_trace        = NULL;
    Vq_trace       = NULL;
    Vc_trace       = NULL;
    Vqq_trace      = NULL;
    Vqc_trace      = NULL;
    Vcc_trace      = NULL;
    basis_size     = -1;
    nindex         = -1;
    eindex         = -1;
    uniform_random = NULL;
    // basic HamiltonianBase info
    UpdateMode.set(COLLECTABLE,1);
    // default values
    energy_mat = false;
    integrator = uniform_grid;
    evaluator  = loop;
    sampling   = volume_based;
    scale      = 1.0;
    center     = 0.0;
    points     = 10;
    samples    = 10;
    warmup     = 30;
    timestep   = .5;
    use_drift  = false;
    warmed_up  = false;
    metric     = 1.0;
    write_acceptance_ratio = false;
    write_rstats  = false;
    normalized    = false;
    check_overlap = false;
    check_derivatives = false;
    // trace data is required
    request.request_scalar("weight");
    request.request_array("Kinetic_complex");
    request.request_array("Vq");
    request.request_array("Vc");
    request.request_array("Vqq");
    request.request_array("Vqc");
    request.request_array("Vcc");
    // has not been initialized
    initialized = false;
  }


  bool DensityMatrices1B::put(xmlNodePtr cur)
  {
    app_log()<<"dm1b put"<< std::endl;

    // read in parameters from the input xml
    set_state(cur);

    // resize local data and perform warmup sampling, if necessary
    initialize();

    // write a report to the log
    report();

    if(check_overlap)
      test_overlap();

    app_log()<<"dm1b end put"<< std::endl;
    return true;
  }


  void DensityMatrices1B::set_state(xmlNodePtr cur)
  {
    app_log()<<"dm1b set_state"<< std::endl;

    bool center_defined=false;
    std::string emstr="no";
    std::string igstr="uniform_grid";
    std::string evstr="loop";
    std::string costr="no";
    std::string cdstr="no";
    std::string arstr="no";
    std::string udstr="no";
    std::string wrstr="no";
    std::string nmstr="no";
    std::vector<std::string> sposets;

    xmlNodePtr element = cur->xmlChildrenNode;
    while(element!=NULL)
    {
      std::string ename((const char*)element->name);
      if(ename=="parameter")
      {
        std::string name((const char*)(xmlGetProp(element,(const xmlChar*)"name")));
        if(name== "basis")      
          putContent(sposets,element);   
        else if(name=="energy_matrix") 
          putContent(emstr,element);        
        else if(name=="integrator") 
          putContent(igstr,element);        
        else if(name=="evaluator") 
          putContent(evstr,element);        
        else if(name=="scale") 
          putContent(scale,element);        
        else if(name=="center")
        {
          putContent(center,element);        
          center_defined = true;
        }
        else if(name=="points") 
          putContent(points,element);       
        else if(name=="samples") 
          putContent(samples,element);      
        else if(name=="warmup") 
          putContent(warmup,element);       
        else if(name=="timestep") 
          putContent(timestep,element);     
        else if(name=="use_drift") 
          putContent(udstr,element);     
        else if(name=="check_overlap") 
          putContent(costr,element);     
        else if(name=="check_derivatives") 
          putContent(cdstr,element);     
        else if(name=="acceptance_ratio") 
          putContent(arstr,element);     
        else if(name=="rstats") 
          putContent(wrstr,element);     
        else if(name=="normalized") 
          putContent(nmstr,element);     
      }
      element = element->next;
    }

    if(scale>1.0+1e-10){
      APP_ABORT("DensityMatrices1B::put  scale must be less than one");
    }else if(scale<0.0-1e-10)
      APP_ABORT("DensityMatrices1B::put  scale must be greater than zero");

    // get volume and cell information
    Lattice.reset();
    if(!center_defined)
      center = Lattice.Center;
    volume     = Lattice.Volume*std::exp(DIM*std::log(scale));
    periodic   = Lattice.SuperCellEnum!=SUPERCELL_OPEN;
    rcorner    = center - scale*Lattice.Center;

    energy_mat = emstr=="yes";
    if(igstr=="uniform_grid")
    {
      integrator = uniform_grid;
      sampling   = volume_based;
      samples    = pow(points,DIM);
      metric     = volume/samples;
      ind_dims[0] = pow(points,DIM-1);
      for(int d=1;d<DIM;++d)
        ind_dims[d] = ind_dims[d-1]/points;
    }
    else if(igstr=="uniform")
    {
      integrator = uniform;
      sampling   = volume_based;
      metric     = volume/samples;
    }
    else if(igstr=="density")
    {
      integrator = density;
      sampling   = metropolis;
      metric     = 1.0/samples;
    }
    else
      APP_ABORT("DensityMatrices1B::set_state  invalid integrator\n  valid options are: uniform_grid, uniform, density");

    if(evstr=="loop")
      evaluator=loop;
    else if(evstr=="matrix")
      evaluator=matrix;
    else
      APP_ABORT("DensityMatrices1B::set_state  invalid evaluator\n  valid options are: loop, matrix");

    normalized        = nmstr=="yes";
    use_drift         = udstr=="yes";
    check_overlap     = costr=="yes";
    check_derivatives = cdstr=="yes";
    write_rstats      = wrstr=="yes";
    write_acceptance_ratio = arstr=="yes";


    // get the sposets that form the basis
    if(sposets.size()==0)
      APP_ABORT("DensityMatrices1B::put  basis must have at least one sposet");

    for(int i=0;i<sposets.size();++i)
    {
      app_log()<<"  sposets requested "<<sposets.size()<<" "<<i<<" "<<sposets[i]<< std::endl;
      app_log()<<"  size before "<<basis_functions.size()<< std::endl;
      basis_functions.add(get_sposet(sposets[i]));
      app_log()<<"  size after  "<<basis_functions.size()<< std::endl;
    }
    basis_size = basis_functions.size();

    if(basis_size<1)
      APP_ABORT("DensityMatrices1B::put  basis_size must be greater than one");


    app_log()<<"dm1b end set_state"<< std::endl;
  }


  void DensityMatrices1B::set_state(DensityMatrices1B& master)
  {
    app_log()<<"dm1b set_state master"<< std::endl;
    basis_size = master.basis_size;
    energy_mat = master.energy_mat;
    integrator = master.integrator;
    evaluator  = master.evaluator;
    sampling   = master.sampling;
    scale      = master.scale;
    points     = master.points;
    samples    = master.samples;
    warmup     = master.warmup;
    timestep   = master.timestep;
    use_drift  = master.use_drift;
    volume     = master.volume;
    periodic   = master.periodic;
    metric     = master.metric;
    rcorner    = master.rcorner;
    normalized = master.normalized;
    for(int d=0;d<DIM;++d)
      ind_dims[d] = master.ind_dims[d];
    app_log()<<"dm1b end set_state master"<< std::endl;
  }


  void DensityMatrices1B::initialize()
  {
    app_log()<<"dm1b initialize"<< std::endl;
    // get particle information
    SpeciesSet& species = Pq.getSpeciesSet();
    nparticles = Pq.getTotalNum();
    nspecies   = species.size();
    int natt   = species.numAttributes();
    int isize  = species.addAttribute("membersize");
    if(isize==natt)
      APP_ABORT("DensityMatrices1B::set_state  Species set does not have the required attribute 'membersize'");
    for(int s=0;s<nspecies;++s)
      species_size.push_back(species(isize,s));
    for(int s=0;s<nspecies;++s)
      species_name.push_back(species.speciesName[s]);

    // allocate space
    basis_values.resize(basis_size);
    integrated_values.resize(basis_size);
    basis_norms.resize(basis_size);
    for(int i=0;i<basis_size;++i)
      basis_norms[i] = 1.0;

    rsamples.resize(samples);
    sample_weights.resize(samples);

    if(evaluator==matrix)
    {
      Phi_MB.resize(samples,basis_size);
      for(int s=0;s<nspecies;++s)
      {
        int specs_size = species_size[s];
        Phi_NB.push_back(    new Matrix_t(specs_size,basis_size));
        Psi_NM.push_back(    new Matrix_t(specs_size,samples   ));
        Phi_Psi_NB.push_back(new Matrix_t(specs_size,basis_size));
        N_BB.push_back(      new Matrix_t(basis_size,basis_size));
        if(energy_mat)
        {
          E_N.push_back( new Vector_t(specs_size));
          E_BB.push_back(new Matrix_t(basis_size,basis_size));
        }
      }

#ifdef DMCHECK
      Phi_MBtmp.resize(samples,basis_size);
      for(int s=0;s<nspecies;++s)
      {
        int specs_size = species_size[s];
        Phi_NBtmp.push_back(    new Matrix_t(specs_size,basis_size));
        Psi_NMtmp.push_back(    new Matrix_t(specs_size,samples   ));
        Phi_Psi_NBtmp.push_back(new Matrix_t(specs_size,basis_size));
        N_BBtmp.push_back(      new Matrix_t(basis_size,basis_size));
        if(energy_mat)
        {
          E_Ntmp.push_back( new Vector_t(specs_size));
          E_BBtmp.push_back(new Matrix_t(basis_size,basis_size));
        }
      }
#endif
    }

    if(sampling==metropolis)
    {
      basis_gradients.resize(basis_size);
      basis_laplacians.resize(basis_size);
    }

    if(!normalized)
      normalize();

    initialized = true;

    app_log()<<"dm1b end initialize"<< std::endl;
  }


  void DensityMatrices1B::finalize()
  {
    delete_iter(     Phi_NB.begin(),     Phi_NB.end() );
    delete_iter(     Psi_NM.begin(),     Psi_NM.end() );
    delete_iter( Phi_Psi_NB.begin(), Phi_Psi_NB.end() );
    delete_iter(       N_BB.begin(),       N_BB.end() );
    if(energy_mat)
    {
      delete_iter(  E_N.begin(),  E_N.end() );
      delete_iter( E_BB.begin(), E_BB.end() );
    }

#ifdef DMCHECK
    delete_iter(     Phi_NBtmp.begin(),     Phi_NBtmp.end() );
    delete_iter(     Psi_NMtmp.begin(),     Psi_NMtmp.end() );
    delete_iter( Phi_Psi_NBtmp.begin(), Phi_Psi_NBtmp.end() );
    delete_iter(       N_BBtmp.begin(),       N_BBtmp.end() );
    if(energy_mat)
    {
      delete_iter(  E_Ntmp.begin(),  E_Ntmp.end() );
      delete_iter( E_BBtmp.begin(), E_BBtmp.end() );
    }
#endif

  }


  void DensityMatrices1B::report(const std::string& pad)
  {
    std::vector<std::string> integrator_list;
    integrator_list.push_back("uniform_grid");
    integrator_list.push_back("uniform");
    integrator_list.push_back("density");
    integrator_list.push_back("no_integrator");
    std::vector<std::string> evaluator_list;
    evaluator_list.push_back("loop");
    evaluator_list.push_back("matrix");
    evaluator_list.push_back("no_evaluator");
    std::vector<std::string> sampling_list;
    sampling_list.push_back("volume_based");
    sampling_list.push_back("metropolis");
    sampling_list.push_back("no_sampling");

    //ostream& out = app_log();
    std::ostream& out = std::cout;

    out<<pad<<"DensityMatrices1B"<< std::endl;

    out<<pad<<"  integrator  = "<< integrator_list[(int)integrator] << std::endl; 
    out<<pad<<"  sampling    = "<< sampling_list[  (int)sampling  ] << std::endl; 
    out<<pad<<"  evaluator   = "<< evaluator_list[ (int)evaluator ] << std::endl; 
    out<<pad<<"  periodic    = "<< periodic    << std::endl;
    if(sampling==volume_based)
    {
      PosType rmax = rcorner + 2*scale*Lattice.Center;
      out<<pad<<"  points      = "<< points      << std::endl;
      out<<pad<<"  scale       = "<< scale       << std::endl;
      out<<pad<<"  center      = "<< center      << std::endl;
      out<<pad<<"  rmin        = "<< rcorner     << std::endl;
      out<<pad<<"  rmax        = "<< rmax        << std::endl;
      out<<pad<<"  volume      = "<< volume      << std::endl;
    }
    else if(sampling==metropolis)
    {
      out<<pad<<"  warmup      = "<< warmup      << std::endl;
      out<<pad<<"  timestep    = "<< timestep    << std::endl;
    }
    out<<pad<<"  metric      = "<< metric      << std::endl;
    out<<pad<<"  nparticles  = "<< nparticles  << std::endl;
    out<<pad<<"  nspecies    = "<< nspecies    << std::endl;
    for(int s=0;s<nspecies;++s)
      out<<pad<<"    species "<<s<<" = "<< species_size[s] << std::endl;
    out<<pad<<"  basis_size  = "<< basis_size  << std::endl;
    out<<pad<<"  samples     = "<< samples     << std::endl;
    out<<pad<<"  energy_mat  = "<< energy_mat  << std::endl;
    out<<pad<<"  initialized = "<< initialized << std::endl;
    out<<pad<<"  rsamples : "<< rsamples.size()<< std::endl;
    if(evaluator==matrix)
    {
      out<<pad<<"  Phi_MB   : "<< Phi_MB.rows()<<" "<<Phi_MB.cols()<< std::endl;
      for(int s=0;s<nspecies;++s)
      {
        out<<pad<<"  matrices/vectors for species "<<s<< std::endl;
        out<<pad<<"    E_N        : "<< E_N[s]->size()<< std::endl;
        out<<pad<<"    Phi_NB     : "<< Phi_NB[s]->rows()<<" "<<Phi_NB[s]->cols()<< std::endl;
        out<<pad<<"    Psi_NM     : "<< Psi_NM[s]->rows()<<" "<<Psi_NM[s]->cols()<< std::endl;
        out<<pad<<"    Phi_Psi_NB : "<< Phi_Psi_NB[s]->rows()<<" "<<Phi_Psi_NB[s]->cols()<< std::endl;
        out<<pad<<"    N_BB       : "<< N_BB[s]->rows()<<" "<<N_BB[s]->cols()<< std::endl;
        out<<pad<<"    E_BB       : "<< E_BB[s]->rows()<<" "<<E_BB[s]->cols()<< std::endl;
      }
    }
    out<<pad<<"  basis_norms"<< std::endl;
    for(int i=0;i<basis_size;++i)
      out<<pad<<"    "<<i<<"  "<<basis_norms[i]<< std::endl;
    out<<pad<<"end DensityMatrices1B"<< std::endl;
  }


  void DensityMatrices1B::get_required_traces(TraceManager& tm)
  {
    app_log()<<"dm1b get_required_traces"<< std::endl;
    w_trace = tm.get_real_trace("weight");
    if(energy_mat)
    {
      T_trace   = tm.get_complex_trace(Pq,"Kinetic_complex");
      Vq_trace  = tm.get_real_combined_trace(Pq,"Vq");
      Vqq_trace = tm.get_real_combined_trace(Pq,"Vqq");
      Vqc_trace = tm.get_real_combined_trace(Pq,"Vqc");
      if(Pc)
      {
        Vc_trace  = tm.get_real_combined_trace(*Pc,"Vc");
        Vcc_trace = tm.get_real_combined_trace(*Pc,"Vcc");
      }

      E_samp.resize(nparticles);

    } 
    have_required_traces = true;
    app_log()<<"dm1b end get_required_traces"<< std::endl;
  }


  void DensityMatrices1B::setRandomGenerator(RandomGenerator_t* rng)
  {
    app_log()<<"dm1b setRandomGenerator"<< std::endl;
    uniform_random = rng;
    app_log()<<"dm1b end setRandomGenerator"<< std::endl;
  }


  void DensityMatrices1B::addObservables(PropertySetType& plist,BufferType& collectables)
  {
    app_log()<<"dm1b addObservables"<< std::endl;
#if defined(QMC_COMPLEX)
    int nentries = 2*basis_size*basis_size*nspecies;
#else
    int nentries = basis_size*basis_size*nspecies;
#endif
    myIndex=collectables.current();
    nindex=myIndex;
    std::vector<RealType> ntmp(nentries);
    collectables.add(ntmp.begin(),ntmp.end());
    if(energy_mat)
    {
      eindex = nindex+nentries;
      std::vector<RealType> etmp(nentries);
      collectables.add(etmp.begin(),etmp.end());
    }
    app_log()<<"dm1b end addObservables"<< std::endl;
  }


  void DensityMatrices1B::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const 
  {
    app_log()<<"dm1b registerCollectables"<< std::endl;
#if defined(QMC_COMPLEX)
    std::vector<int> ng(3);
    ng[0] = basis_size;
    ng[1] = basis_size;
    ng[2] = 2;
    int nentries = ng[0]*ng[1]*ng[2];
#else
    std::vector<int> ng(2);
    ng[0] = basis_size;
    ng[1] = basis_size;
    int nentries = ng[0]*ng[1];
#endif

    std::string dname = myName;
    hid_t dgid=H5Gcreate(gid,dname.c_str(),0);

    std::string nname = "number_matrix";
    hid_t ngid=H5Gcreate(dgid,nname.c_str(),0);
    for(int s=0;s<nspecies;++s)
    {
      observable_helper* oh;
      oh = new observable_helper(species_name[s]);
      oh->set_dimensions(ng,nindex+s*nentries);
      oh->open(ngid);
      h5desc.push_back(oh);
    }

    if(energy_mat)
    { 
      std::string ename = "energy_matrix";
      hid_t egid=H5Gcreate(dgid,ename.c_str(),0);
      for(int s=0;s<nspecies;++s)
      {
        observable_helper* oh;
        oh = new observable_helper(species_name[s]);
        oh->set_dimensions(ng,eindex+s*nentries);
        oh->open(egid);
        h5desc.push_back(oh);
      }
    }
    app_log()<<"dm1b end registerCollectables"<< std::endl;
  }


  void DensityMatrices1B::warmup_sampling()
  {
    if(sampling==metropolis)
    {
      if(!warmed_up)
      {
        diffusion(std::sqrt(timestep),rpcur);
        rpcur += center;
        if(integrator==density)
          density_drift(rpcur,rhocur,dpcur);
        else
          APP_ABORT("DensityMatrices1B::warmup_sampling  invalid integrator");
      }
      generate_samples(1.0,warmup);
      warmed_up = true;
    }
  }


  DensityMatrices1B::Return_t DensityMatrices1B::evaluate(ParticleSet& P)
  {
    if(have_required_traces)
    {
      if(check_derivatives)
        test_derivatives();
      if(evaluator==loop)
        evaluate_loop(P);
      else if(evaluator==matrix)
        evaluate_matrix(P);
      else
        APP_ABORT("DensityMatrices1B::evaluate  invalid evaluator");
    }
    return 0.0;
  }


  DensityMatrices1B::Return_t DensityMatrices1B::evaluate_matrix(ParticleSet& P)
  {
    //perform warmup sampling the first time
    if(!warmed_up)
      warmup_sampling();
    // get weight and single particle energy trace data
    RealType weight=w_trace->sample[0]*metric;
    if(energy_mat)
      get_energies(E_N);               // energies        : particles x 1
    // compute sample positions (monte carlo or deterministic)
    generate_samples(weight);
    // compute basis and wavefunction ratio values in matrix form 
    generate_sample_basis(Phi_MB);     // basis           : samples   x basis_size
    generate_sample_ratios(Psi_NM);    // qmcplusplus::conj(Psi ratio) : particles x samples
    generate_particle_basis(P,Phi_NB); // qmcplusplus::conj(basis)     : particles x basis_size
    // perform integration via matrix products
    for(int s=0;s<nspecies;++s)
    {
      Matrix_t& Psi_nm     = *Psi_NM[s];
      Matrix_t& Phi_Psi_nb = *Phi_Psi_NB[s];
      Matrix_t& Phi_nb     = *Phi_NB[s];
      diag_product(Psi_nm,sample_weights,Psi_nm);
      product(Psi_nm,Phi_MB,Phi_Psi_nb);       // ratio*basis : particles x basis_size
      product_AtB(Phi_nb,Phi_Psi_nb,*N_BB[s]);  // qmcplusplus::conj(basis)^T*ratio*basis : basis_size^2
      if(energy_mat)
      {
        Vector_t& E = *E_N[s];
        diag_product(E,Phi_nb,Phi_nb);         // diag(energies)*qmcplusplus::conj(basis)
        product_AtB(Phi_nb,Phi_Psi_nb,*E_BB[s]);// (energies*qmcplusplus::conj(basis))^T*ratio*basis
      }
    }
    // accumulate data into collectables
    const int basis_size2 = basis_size*basis_size;
    int ij=nindex;
    for(int s=0;s<nspecies;++s)
    {
      //int ij=nindex; // for testing
      const Matrix_t& NDM = *N_BB[s];
      for(int n=0;n<basis_size2;++n)
      {
        Value_t val = NDM(n);
        P.Collectables[ij] += real(val);
        ij++;
#if defined(QMC_COMPLEX)
        P.Collectables[ij] += imag(val);
        ij++;
#endif
      }
    }
    if(energy_mat)
    {
      int ij=eindex;
      for(int s=0;s<nspecies;++s)
      {
        //int ij=eindex; // for testing
        const Matrix_t& EDM = *E_BB[s];
        for(int n=0;n<basis_size2;++n)
        {
          Value_t val = EDM(n);
          P.Collectables[ij] += real(val);
          ij++;
#if defined(QMC_COMPLEX)
          P.Collectables[ij] += imag(val);
          ij++;
#endif
        }
      }
    }


    // jtk come back to this
    //   why are matrices so similar across species?
    //   O atom, einspline, points=20, blocks=1, steps=1
    //
    //app_log()<<" ntraces = ";
    //for(int s=0;s<nspecies;++s)
    //{
    //  const Matrix_t& NDM = *N_BB[s];
    //  Value_t ntrace = 0.0;
    //  for(int i=0;i<basis_size;++i)
    //    ntrace+=NDM(i,i);
    //  app_log()<<ntrace<<"  ";
    //}
    //app_log()<< std::endl;
    //
    //app_log()<<"nmatrices"<< std::endl;
    //for(int s=0;s<nspecies;++s)
    //{
    //  app_log()<<" species "<<s<< std::endl;
    //  const Matrix_t& NDM = *N_BB[s];
    //  for(int i=0;i<basis_size;++i)
    //  {
    //    for(int j=0;j<basis_size;++j)
    //      app_log()<<"  "<<real(NDM(i,j));
    //    app_log()<< std::endl;
    //  }
    //}
    //
    //app_log()<<"positions"<< std::endl;
    //int p=0;
    //for(int s=0;s<nspecies;++s)
    //{
    //  app_log()<<" species "<<s<< std::endl;
    //  for(int ps=0;ps<species_size[s];++ps,++p)
    //    app_log()<<"  "<<p<<"  "<<P.R[p]-P.Lattice.Center<< std::endl;
    //}
    //
    ////app_log()<<"basis_values"<< std::endl;
    ////int p=0;
    ////for(int s=0;s<nspecies;++s)
    ////{
    ////  app_log()<<" species "<<s<< std::endl;
    ////  for(int ps=0;ps<species_size[s];++ps,++p)
    ////    app_log()<<"  "<<p<<"  "<<P.R[p]<< std::endl;
    ////}

    
#ifdef DMCHECK
    report();
    app_log()<<"DM Check"<< std::endl;
    evaluate_check(P);
    compare("  Phi_MB",Phi_MB,Phi_MBtmp);
    for(int s=0;s<nspecies;++s)
    {
      app_log()<<"  species "<<s<< std::endl;
      compare("    E_N       ", *E_N[s]       , *E_Ntmp[s]);       
      compare("    Phi_NB    ", *Phi_NB[s]    , *Phi_NBtmp[s]);    
      compare("    Psi_NM    ", *Psi_NM[s]    , *Psi_NMtmp[s]);    
      compare("    Phi_Psi_NB", *Phi_Psi_NB[s], *Phi_Psi_NBtmp[s],true);
      compare("    N_BB      ", *N_BB[s]      , *N_BBtmp[s],true);      
      compare("    E_BB      ", *E_BB[s]      , *E_BBtmp[s],true);       
    }
    app_log()<<"end DM Check"<< std::endl;
    APP_ABORT("DM Check");
#endif


    return 0.0;
  }


  DensityMatrices1B::Return_t DensityMatrices1B::evaluate_check(ParticleSet& P)
  {
#ifdef DMCHECK
    APP_ABORT("DensityMatrices1B::evaluate_check  use of E_trace in this function needs to be replaces with get_energies() and E_samp");
    int n=0;
    for(int s=0;s<nspecies;++s)
    {
      Matrix_t& Phi_mb     =  Phi_MBtmp;
      Matrix_t& Psi_nm     = *Psi_NMtmp[s];
      Matrix_t& Phi_Psi_nb = *Phi_Psi_NBtmp[s];
      Matrix_t& Phi_nb     = *Phi_NBtmp[s];
      Vector_t& E_n        = *E_Ntmp[s];
      Matrix_t& N_bb       = *N_BBtmp[s];
      Matrix_t& E_bb       = *E_BBtmp[s];

      for(int ij=0;ij<basis_size*basis_size;++ij)
        N_bb(ij) = 0.0;
      for(int ij=0;ij<basis_size*basis_size;++ij)
        E_bb(ij) = 0.0;

      for(int ns=0;ns<species_size[s];++ns,++n)
      {
        std::fill(integrated_values.begin(),integrated_values.end(),0.0);
        for(int m=0;m<samples;++m)
        {
          PosType& rsamp = rsamples[m];
          update_basis(rsamp);
          PosType dr = rsamp-P.R[n];
          P.makeMove(n,dr);
          Value_t ratio = sample_weights[m]*qmcplusplus::conj(Psi.full_ratio(P,n));
          P.rejectMove(n);
          for(int i=0;i<basis_size;++i)
          {
            integrated_values[i]+=ratio*basis_values[i];
            Phi_mb(m,i) = basis_values[i];
          }
          Psi_nm(ns,m) = ratio;
        }
        update_basis(P.R[n]);
        for(int i=0;i<basis_size;++i)
          Phi_Psi_nb(ns,i) = integrated_values[i];
        for(int i=0;i<basis_size;++i)
          Phi_nb(ns,i) = qmcplusplus::conj(basis_values[i]);
        for(int i=0;i<basis_size;++i)
        {
          Value_t phi_i = qmcplusplus::conj(basis_values[i]);
          for(int j=0;j<basis_size;++j)
          {
            Value_t val = phi_i*integrated_values[j];
            N_bb(i,j) +=  val;
          }
        }
        if(energy_mat)
        {
          RealType e_n = E_trace->sample[n]; //replace this with traces access later
          E_n[ns] = e_n;
          for(int i=0;i<basis_size;++i)
          {
            Value_t ephi_i = e_n*qmcplusplus::conj(basis_values[i]);
            Phi_nb(ns,i) = ephi_i;
            for(int j=0;j<basis_size;++j)
            {
              Value_t val = ephi_i*integrated_values[j];
              E_bb(i,j) += val;
            }
          }
        }
      }
    }
#endif
    return 0.0;
  }


  DensityMatrices1B::Return_t DensityMatrices1B::evaluate_loop(ParticleSet& P)
  {
    const int basis_size2 = basis_size*basis_size;
    if(!warmed_up)
      warmup_sampling();
    RealType weight=w_trace->sample[0]*metric;
    int nparticles = P.getTotalNum();
    generate_samples(weight);
    int n=0;
    for(int s=0;s<nspecies;++s)
    {
      for(int ns=0;ns<species_size[s];++ns,++n)
      {
        integrate(P,n);
        update_basis(P.R[n]);
        int ij = nindex + s*basis_size2;
        for(int i=0;i<basis_size;++i)
        {
          Value_t phi_i = qmcplusplus::conj(basis_values[i]);
          for(int j=0;j<basis_size;++j)
          {
            Value_t val = phi_i*integrated_values[j];
            P.Collectables[ij] += real(val);
            ij++;
#if defined(QMC_COMPLEX)
            P.Collectables[ij] += imag(val);
            ij++;
#endif
          }
        }
        if(energy_mat)
        {
          RealType e_n = E_trace->sample[n]; //replace this with traces access later
          int ij = eindex + s*basis_size2;
          for(int i=0;i<basis_size;++i)
          {
            Value_t ephi_i = e_n*qmcplusplus::conj(basis_values[i]);
            for(int j=0;j<basis_size;++j)
            {
              Value_t val = ephi_i*integrated_values[j];
              P.Collectables[ij] += real(val);
              ij++;
#if defined(QMC_COMPLEX)
              P.Collectables[ij] += imag(val);
              ij++;
#endif
            }
          }
        }
      }
    }
    return 0.0;
  }


  inline void DensityMatrices1B::generate_samples(RealType weight,int steps)
  {
    RandomGenerator_t& rng = *uniform_random;
    bool save = false;
    if(steps==0)
    {
      save  = true;
      steps = samples;
    }
    if(integrator==uniform_grid)
      generate_uniform_grid(rng);
    else if(integrator==uniform)
      generate_uniform_samples(rng);
    else if(integrator==density)
      generate_density_samples(save,steps,rng);
    
    if(save)
    {
      if(sampling==metropolis)
      {
        sample_weights*=weight;
      }
      else
      {
        std::fill(sample_weights.begin(),sample_weights.end(),weight);
      }
    }


    // temporary check
    if(write_rstats && omp_get_thread_num()==0)
    {
      PosType rmin = std::numeric_limits<RealType>::max();
      PosType rmax = -std::numeric_limits<RealType>::max();
      PosType rmean = 0.0;
      PosType rstd = 0.0;
      for(int s=0;s<rsamples.size();++s)
        for(int d=0;d<DIM;++d)
        {
          RealType rd = rsamples[s][d];
          rmin[d] = std::min(rmin[d],rd);
          rmax[d] = std::max(rmax[d],rd);
          rmean[d] += rd;
          rstd[d] += rd*rd;
        }
      rmean/=rsamples.size();
      rstd/=rsamples.size();
      for(int d=0;d<DIM;++d)
        rstd[d] = std::sqrt(rstd[d]-rmean[d]*rmean[d]);
      app_log()<<"\nrsamples properties:"<< std::endl;
      app_log()<<"  rmin  = "<<rmin<< std::endl;
      app_log()<<"  rmax  = "<<rmax<< std::endl;
      app_log()<<"  rmean = "<<rmean<< std::endl;
      app_log()<<"  rstd  = "<<rstd<< std::endl;
    }
  }


  inline void DensityMatrices1B::generate_uniform_grid(RandomGenerator_t& rng)
  {
    PosType rp;
    PosType ushift = 0.0;
    RealType du = scale/points;
    for(int d=0;d<DIM;++d)
      ushift[d] += rng()*du;
    for(int s=0;s<samples;++s)
    {
      int nrem = s;
      for(int d=0;d<DIM-1;++d)
      {
        int ind = nrem/ind_dims[d];
        rp[d] = ind*du + ushift[d];
        nrem-=ind*ind_dims[d];
      }
      rp[DIM-1] = nrem*du + ushift[DIM-1];
      rsamples[s] = Lattice.toCart(rp) + rcorner;
    }
  }


  inline void DensityMatrices1B::generate_uniform_samples(RandomGenerator_t& rng)
  {
    PosType rp;
    for(int s=0;s<samples;++s)
    {
      for(int d=0;d<DIM;++d)
        rp[d] = scale*rng();
      rsamples[s] = Lattice.toCart(rp) + rcorner;
    }
  }


  inline void DensityMatrices1B::generate_density_samples(bool save,int steps,RandomGenerator_t& rng)
  {
    RealType sqt = std::sqrt(timestep);
    RealType ot  = 1.0/timestep;
    PosType  r   = rpcur;   //current position
    PosType  d   = dpcur;   //current drift
    RealType rho = rhocur;  //current density
    for(int s=0;s<steps;++s)
    {
      nmoves++;
      PosType  n,rp,dp,ds;        //diffusion, trial pos/drift, drift sum
      RealType rhop,ratio,Pacc;   //trial density, dens ratio, acc prob
      diffusion(sqt,n);           //get diffusion
      if(use_drift)
      {
        rp    = r + n + d;          //update trial position
        density_drift(rp,rhop,dp);  //get trial drift and density
        ratio = rhop/rho;           //density ratio
        ds    = dp+d;               //drift sum
        Pacc  = ratio*std::exp(-ot*(dot(n,ds)+.5*dot(ds,ds))); //acceptance probability
      }
      else
      {
        rp    = r + n;              //update trial position
        density_only(rp,rhop);      //get trial density
        ratio = rhop/rho;           //density ratio
        Pacc  = ratio;              //acceptance probability
      }
      if(rng()<Pacc)
      {//accept move
        r   = rp;
        d   = dp;
        rho = rhop;
        naccepted++;
      }
      if(save)
      {
        rsamples[s]       = r;
        sample_weights[s] = 1.0/rho;
      }
    }
    acceptance_ratio = RealType(naccepted)/nmoves;
    
    if(write_acceptance_ratio && omp_get_thread_num()==0)
      app_log()<<"dm1b  acceptance_ratio = "<<acceptance_ratio<< std::endl;

    rpcur  = r;
    dpcur  = d;
    rhocur = rho;
  }


  inline void DensityMatrices1B::diffusion(RealType sqt,PosType& diff)
  {
    assignGaussRand(&diff[0],DIM,*uniform_random);
    diff *= sqt;
  }

  
  inline void DensityMatrices1B::density_only(const PosType& r,RealType& dens)
  {
    update_basis(r);
    dens=0.0;
    for(int i=0;i<basis_size;++i)
    {
      Value_t b = basis_values[i];
      dens += std::abs(qmcplusplus::conj(b)*b);
    }
    dens /= basis_size;
  }

  
  inline void DensityMatrices1B::density_drift(const PosType& r,RealType& dens,PosType& drift)
  {
    update_basis_d012(r);
    dens=0.0;
    drift=0.0;
    for(int i=0;i<basis_size;++i)
    {
      const Grad_t& bg = basis_gradients[i];
      Value_t b  = basis_values[i];
      Value_t bc = qmcplusplus::conj(b);
      dens  += std::abs(bc*b);
      for(int d=0;d<DIM;++d)
        drift[d] += prod_real(bc,bg[d]);
    }
    drift *= timestep/dens;
    dens  /= basis_size;
  }


  typedef DensityMatrices1B::RealType RealType;
  typedef DensityMatrices1B::Value_t Value_t;


  inline RealType accum_constant(CombinedTraceSample<TraceReal>* etrace,RealType weight=1.0)
  {
    RealType E=0.0;
    if(etrace)
    {
      etrace->combine();
      for(int p=0;p<etrace->sample.size();++p)
        E += etrace->sample[p];
    }
    return weight*E;
  }

  template<typename T>
  inline void accum_sample(std::vector<Value_t>& E_samp,CombinedTraceSample<T>* etrace,RealType weight=1.0)
  {
    if(etrace)
    {
      etrace->combine();
      for(int p=0;p<etrace->sample.size();++p)
        E_samp[p] += weight*etrace->sample[p];
    }
  }

  template<typename T>
  inline void accum_sample(std::vector<Value_t>& E_samp,TraceSample<T>* etrace,RealType weight=1.0)
  {
#ifdef QMC_COMPLEX
    if(etrace)
      for(int p=0;p<etrace->sample.size();++p)
        E_samp[p] += weight*etrace->sample[p];
#else
    if(etrace)
      for(int p=0;p<etrace->sample.size();++p)
        E_samp[p] += weight*real(etrace->sample[p]);
#endif
  }


  void DensityMatrices1B::get_energies(std::vector<Vector_t*>& E_n)
  {
    Value_t Vc = 0;
    Vc += accum_constant(Vc_trace);
    Vc += accum_constant(Vcc_trace);
    Vc /= nparticles;
    std::fill(E_samp.begin(),E_samp.end(),Vc);
    accum_sample(E_samp,T_trace);
    accum_sample(E_samp,Vq_trace);
    accum_sample(E_samp,Vqq_trace);
    accum_sample(E_samp,Vqc_trace,2.0);

    int p=0;
    for(int s=0;s<nspecies;++s)
    {
      Vector_t& E = *E_n[s];
      for(int ps=0;ps<E.size();++ps,++p)
        E[ps] = E_samp[p];
    }

    //E_trace->combine();
    //RealType E = 0.0;
    //for(int p=0;p<nparticles;++p)
    //  E += E_samp[p];
    //app_log()<<"  E = "<<E<<"  "<<E_trace->sample[0]<< std::endl;
    //APP_ABORT("dm1b::get_energies  check sp traces");
  }


  void DensityMatrices1B::generate_sample_basis(Matrix_t& Phi_mb)
  {
    int mb=0;
    for(int m=0;m<samples;++m)
    {
      update_basis(rsamples[m]);
      for(int b=0;b<basis_size;++b,++mb)
        Phi_mb(mb) = basis_values[b];
    }
  }


  void DensityMatrices1B::generate_sample_ratios(std::vector<Matrix_t*> Psi_nm)
  {
    int p=0;
    for(int s=0;s<nspecies;++s)
    {
      int nm=0;
      Matrix_t& P_nm = *Psi_nm[s];
      for(int n=0;n<species_size[s];++n,++p)
      {
        PosType& Rp = Pq.R[p];
        for(int m=0;m<samples;++m,++nm)
        {
          Pq.makeMove(p,rsamples[m]-Rp);
          P_nm(nm) = qmcplusplus::conj(Psi.full_ratio(Pq,p));
          Pq.rejectMove(p);
        }
      }
    }
  }


  void DensityMatrices1B::generate_particle_basis(ParticleSet& P,std::vector<Matrix_t*>& Phi_nb)
  {
    int p=0;
    for(int s=0;s<nspecies;++s)
    {
      int nb=0;
      Matrix_t& P_nb = *Phi_nb[s];
      for(int n=0;n<species_size[s];++n,++p)
      {
        update_basis(P.R[p]);
        for(int b=0;b<basis_size;++b,++nb)
          P_nb(nb) = qmcplusplus::conj(basis_values[b]);
      }
    }
  }


  inline void DensityMatrices1B::integrate(ParticleSet& P,int n)
  {
    std::fill(integrated_values.begin(),integrated_values.end(),0.0);
    for(int s=0;s<samples;++s)
    {
      PosType& rsamp = rsamples[s];
      update_basis(rsamp);
      P.makeMove(n,rsamp-P.R[n]);
      Value_t ratio = sample_weights[s]*qmcplusplus::conj(Psi.full_ratio(P,n));
      P.rejectMove(n);
      for(int i=0;i<basis_size;++i)
        integrated_values[i]+=ratio*basis_values[i];
    }
  }


  inline void DensityMatrices1B::update_basis(const PosType& r)
  {
    Pq.makeMove(0,r-Pq.R[0]);
    basis_functions.evaluate(Pq,0,basis_values);
    Pq.rejectMove(0);
    for(int i=0;i<basis_size;++i)
      basis_values[i]*=basis_norms[i];
  }


  inline void DensityMatrices1B::update_basis_d012(const PosType& r)
  {
    Pq.makeMove(0,r-Pq.R[0]);
    basis_functions.evaluate(Pq,0,basis_values,basis_gradients,basis_laplacians);
    Pq.rejectMove(0);
    for(int i=0;i<basis_size;++i)
      basis_values[i]*=basis_norms[i];
    for(int i=0;i<basis_size;++i)
      basis_gradients[i]*=basis_norms[i];
    for(int i=0;i<basis_size;++i)
      basis_laplacians[i]*=basis_norms[i];
  }


  inline void DensityMatrices1B::normalize()
  {
    int   ngrid = std::max(200,points);
    int   ngtot = pow(ngrid,DIM);
    RealType du = scale/ngrid;
    RealType dV = volume/ngtot;
    PosType rp;
    ValueVector_t bnorms;
    int gdims[DIM];
    gdims[0] = pow(ngrid,DIM-1);
    for(int d=1;d<DIM;++d)
      gdims[d] = gdims[d-1]/ngrid;
    bnorms.resize(basis_size);
    for(int i=0;i<basis_size;++i)
      bnorms[i] = 0.0;
    std::fill(basis_norms.begin(),basis_norms.end(),1.0);
    for(int p=0;p<ngtot;++p)
    {
      int nrem = p;
      for(int d=0;d<DIM-1;++d)
      {
        int ind = nrem/gdims[d];
        rp[d] = ind*du+du/2;
        nrem-=ind*gdims[d];
      }
      rp[DIM-1] = nrem*du+du/2;
      rp = Lattice.toCart(rp) + rcorner;
      update_basis(rp);
      for(int i=0;i<basis_size;++i)
        bnorms[i] += qmcplusplus::conj(basis_values[i])*basis_values[i]*dV;
    }
    for(int i=0;i<basis_size;++i)
      basis_norms[i] = 1.0/std::sqrt(real(bnorms[i]));
    normalized = true;
  }


  inline void DensityMatrices1B::test_overlap()
  {
    int ngrid = std::max(200,points);
    int ngtot = pow(ngrid,DIM);

    PosType rp;
    RealType du = scale/ngrid;
    RealType dV = volume/ngtot;

    PosType rmin = std::numeric_limits<RealType>::max();
    PosType rmax = -std::numeric_limits<RealType>::max();
    int gdims[DIM];
    gdims[0] = pow(ngrid,DIM-1);
    for(int d=1;d<DIM;++d)
      gdims[d] = gdims[d-1]/ngrid;

    Array<Value_t,2> omat;
    omat.resize(basis_size,basis_size);
    for(int i=0;i<basis_size;++i)
      for(int j=0;j<basis_size;++j)
        omat(i,j) = 0.0;

    for(int p=0;p<ngtot;++p)
    {
      int nrem = p;
      for(int d=0;d<DIM-1;++d)
      {
        int ind = nrem/gdims[d];
        rp[d] = ind*du+du/2;
        nrem-=ind*gdims[d];
      }
      rp[DIM-1] = nrem*du+du/2;
      rp = Lattice.toCart(rp) + rcorner;
      update_basis(rp);
      for(int i=0;i<basis_size;++i)
        for(int j=0;j<basis_size;++j)
          omat(i,j) += qmcplusplus::conj(basis_values[i])*basis_values[j]*dV;
      for(int d=0;d<DIM;++d)
      {
        rmin[d] = std::min(rmin[d],rp[d]);
        rmax[d] = std::max(rmax[d],rp[d]);
      }
    }

    app_log()<<"DensityMatrices1B::test_overlap  checking overlap matrix"<< std::endl;
    app_log()<<"  rmin = "<<rmin<< std::endl;
    app_log()<<"  rmax = "<<rmax<< std::endl;
    app_log()<<"  overlap scale "<<std::abs(omat(0,0))<< std::endl;
    app_log()<<"  overlap matrix:"<< std::endl;
    for(int i=0;i<basis_size;++i){
      app_log()<< std::endl;
      for(int j=0;j<basis_size;++j)
        app_log()<<std::abs(omat(i,j))/std::abs(omat(0,0))<<" ";
    }
    app_log()<< std::endl;
    
    APP_ABORT("DensityMatrices1B::test_overlap");
  }


  void DensityMatrices1B::test_derivatives()
  {
    app_log()<<"DensityMatrices1B::test_derivatives  checking drift"<< std::endl;

    PosType r,rtmp;

    RealType delta = 1e-5;

    RealType dens,densp,densm;
    PosType  drift,driftn,drifttmp;

    app_log()<<"  warming up"<< std::endl;
    warmup_sampling();
    app_log()<<"  generating samples"<< std::endl;
    generate_samples(1.0);

    app_log()<<"  testing derivatives at sample points"<< std::endl;
    for(int s=0;s<rsamples.size();++s)
    {
      r = rsamples[s];

      density_drift(r,dens,drift);

      for(int d=0;d<DIM;++d)
      {
        rtmp = r;

        rtmp[d] = r[d]+delta;
        density_drift(rtmp,densp,drifttmp);

        rtmp[d] = r[d]-delta;
        density_drift(rtmp,densm,drifttmp);

        driftn[d] = (densp-densm)/(2*delta);
      }
      driftn *= .5*timestep/dens;

      app_log()<<s<< std::endl;
      app_log()<<"  "<<driftn<< std::endl;
      app_log()<<"  "<<drift<< std::endl;
    }

    APP_ABORT("DensityMatrices1B::test_derivatives");

  }


  void DensityMatrices1B::
  postprocess_density(const std::string& infile,const std::string& species,
                      pts_t& points,dens_t& density,dens_t& density_err)
  {
    std::ifstream datafile;
    datafile.open(infile.c_str());
    if(!datafile.is_open())
      APP_ABORT("DensityMatrices1B::postprocess_density\n  could not open file: "+infile);

    const int nelem = basis_size*basis_size;
    int n=0;
    std::vector<RealType> nmr;
    std::vector<RealType> nmi;
    std::vector<RealType> nmr_err;
    std::vector<RealType> nmi_err;
    RealType value;    
    while(datafile>>value)
    {
      if(n<nelem)
        nmr.push_back(value);
      else if(n<2*nelem)
        nmr_err.push_back(value);
      else if(n<3*nelem)
        nmi.push_back(value);
      else if(n<4*nelem)
        nmi_err.push_back(value);
      n++;
    }
    std::vector<Value_t> nm;
    std::vector<Value_t> nmv;
    nm.resize(nelem);
    nmv.resize(nelem);
#if defined(QMC_COMPLEX)
    const int nmat = 4;
    for(int ij=0;ij<nelem;++ij)
      nm[ij] = Value_t(nmr[ij],nmi[ij]);
    for(int ij=0;ij<nelem;++ij)
    {
      Value_t s = Value_t(nmr_err[ij],nmi_err[ij]);
      nmv[ij] = qmcplusplus::conj(s)*s; 
    }
#else
    const int nmat = 2;
    for(int ij=0;ij<nelem;++ij)
      nm[ij] = nmr[ij];
    for(int ij=0;ij<nelem;++ij)
      nmv[ij] = nmr_err[ij]*nmr_err[ij];
#endif
    if(n!=nmat*nelem)
    {
      app_log()<<"DensityMatrices1B::postprocess_density\n  file "<<infile<<"\n  contains "<< n <<" values\n  expected "<<nmat*nelem<<" values"<< std::endl;
      APP_ABORT("DensityMatrices1B::postprocess_density");
    }

    //evaluate the density
    //  this error analysis neglects spatial (i,j) correlations
    for(int p=0;p<points.size();++p)
    {
      RealType d  = 0.0;
      RealType de = 0.0;
      update_basis(points[p]);
      int ij=0;
      for(int i=0;i<basis_size;++i)
        for(int j=0;j<basis_size;++j,++ij)
        {
          Value_t aij = basis_values[i]*qmcplusplus::conj(basis_values[j]);
          d  += real(aij*nm[ij]);
          de += real(qmcplusplus::conj(aij)*aij*nmv[ij]);
        }
      de = std::sqrt(de);
      density[p]     = d;
      density_err[p] = de;
    }
    
  }





  bool DensityMatrices1B::match(Value_t e1, Value_t e2, RealType tol)
  {
    return std::abs(e1-e2) < tol;
    //return std::abs(2*(e1-e2)/(e1+e2)) < tol;
  }


  bool DensityMatrices1B::same(Vector_t& v1, Vector_t& v2, RealType tol)
  {
    if(v1.size()!=v2.size())
      APP_ABORT("DensityMatrices1B::same(vector)  vectors differ in size");
    bool sm = true;
    for(int i=0;i<v1.size();++i)
      sm &= match(v1[i],v2[i]);
    return sm;
  }

  bool DensityMatrices1B::same(Matrix_t& m1, Matrix_t& m2, RealType tol)
  {
    if(m1.rows()!=m2.rows() || m1.cols()!=m2.cols())
      APP_ABORT("DensityMatrices1B::same(matrix)  matrices differ in size");
    bool sm = true;
    int n = m1.rows()*m1.cols();
    for(int i=0;i<n;++i)
      sm &= match(m1(i),m2(i));
    return sm;
  }

  void DensityMatrices1B::compare(const std::string& name, Vector_t& v1, Vector_t& v2,bool write,bool diff_only)
  {
    bool sm = same(v1,v2);
    std::string result="differ";
    if(sm)
      result = "agree";
    app_log()<<name<<" "<<result<< std::endl;
    if(write && !sm)
      for(int i=0;i<v1.size();++i)
        app_log()<<"      "<<i<<" "<<real(v1[i])<<" "<<real(v2[i])
                 <<" "<<real(v1[i]/v2[i])<<" "<<real(v2[i]/v1[i])<< std::endl;
  }

  void DensityMatrices1B::compare(const std::string& name, Matrix_t& m1, Matrix_t& m2,bool write,bool diff_only)
  {
    bool sm = same(m1,m2);
    std::string result="differ";
    if(sm)
      result = "agree";
    app_log()<<name<<" "<<result<< std::endl;
    if(write && !sm)
      for(int i=0;i<m1.rows();++i)
        for(int j=0;j<m1.cols();++j)
          if(!diff_only || !match(m1(i,j),m2(i,j)) )
            app_log()<<"      "<<i<<" "<<j<<" "<<real(m1(i,j))<<" "<<real(m2(i,j))
                     <<" "<<real(m1(i,j)/m2(i,j))<< std::endl;
  }








}
