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
    
    

#include <QMCWaveFunctions/HarmonicOscillator/SHOSetBuilder.h>
#include <QMCWaveFunctions/SPOSetInputInfo.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/IteratorUtility.h>
#include <Utilities/string_utils.h>


namespace qmcplusplus
{

  SHOSetBuilder::SHOSetBuilder(ParticleSet& P) 
  : Ps(P)
  {
    legacy = false;
    app_log()<<"Constructing SHOSetBuilder"<< std::endl;
    reset();
  }


  SHOSetBuilder::~SHOSetBuilder()
  {
  }


  void SHOSetBuilder::reset()
  {
    nstates = 0;
    mass    = -1.0;
    energy  = -1.0;
    length  = -1.0;
    center  = 0.0;
  }
  

  SPOSetBase* SHOSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
  {
    APP_ABORT("SHOSetBuilder::createSPOSetFromXML  SHOSetBuilder should not use legacy interface");

    app_log() << "SHOSetBuilder::createSHOSet(xml) " << std::endl;

    SPOSetInputInfo input(cur);

    return createSPOSet(cur,input);
  }


  SPOSetBase* SHOSetBuilder::createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input)
  {
    app_log() << "SHOSetBuilder::createSHOSet(indices) " << std::endl;

    using std::sqrt;
    using std::ceil;

    reset();

    // read parameters
    std::string spo_name = "sho";
    OhmmsAttributeSet attrib;
    attrib.add(spo_name,"name"     );
    attrib.add(spo_name,"id"       );
    attrib.add(mass,    "mass"     );
    attrib.add(energy,  "energy"   );
    attrib.add(energy,  "frequency");
    attrib.add(length,  "length"   );
    attrib.add(center,  "center"   );
    attrib.add(nstates, "size"     );
    attrib.put(cur);
    if(energy<0.0)
      energy = 1.0;
    if(mass<0.0 && length<0.0)
      length = 1.0;
    if(mass<0.0)
      mass = 1.0/(energy*length*length);
    else if(length<0.0)
      length = 1.0/sqrt(mass*energy);

    // initialize states and/or adjust basis
    int smax = -1;
    if(input.has_index_info)
      smax = std::max(smax,input.max_index());
    if(input.has_energy_info)
    {
      smax = std::max(smax,(int)ceil(input.max_energy()/energy));
    }
    if(smax<0)
      APP_ABORT("SHOSetBuilder::Initialize\n  invalid basis size");
    update_basis_states(smax);

    // create sho state request
    indices_t& indices = input.get_indices(states);
    std::vector<SHOState*> sho_states;
    for(int i=0;i<indices.size();++i)
      sho_states.push_back(basis_states[indices[i]]);

    // make the sposet
    SHOSet* sho = new SHOSet(length,center,sho_states);

    sho->report("  ");
    //sho->test_derivatives();
    //sho->test_overlap();
    //APP_ABORT("SHOSetBuilder check");

    return sho;
  }


  void SHOSetBuilder::update_basis_states(int smax)
  {
    using std::ceil;
    using std::sqrt;
    using std::sort;
    using std::exp;
    using std::log;

    int states_required = smax-basis_states.size()+1;
    if(states_required>0)
    {
      RealType N = smax+1;
      if(DIM==1)
        nmax = smax;
      else if(DIM==2)
        nmax = ceil(.5*sqrt(8.*N+1.)-1.5);
      else if(DIM==3)
      {
        RealType f = exp( 1.0/3.0 * log(81.*N+3.*sqrt(729.*N*N-3.)) );
        nmax = ceil(f/3.+1./f-2.);
      }else
        APP_ABORT("SHOSetBuilder::update_basis_states  dimensions other than 1, 2, or 3 are not supported");
      int ndim = nmax+1;
      ind_dims[DIM-1] = 1;
      for(int d=DIM-2;d>-1;--d)
        ind_dims[d] = ind_dims[d+1]*ndim;
      int s=0;
      int ntot = pow(ndim,DIM);
      TinyVector<int,DIM> qnumber;
      for(int m=0;m<ntot;++m)
      {
        int n=0; // principal quantum number
        int nrem = m;
        for(int d=0;d<DIM;++d)
        {
          int i = nrem/ind_dims[d];
          nrem -= i*ind_dims[d];
          qnumber[d]=i;
          n+=i;
        }
        if(n<=nmax)
        {
          SHOState* st;
          if(s<basis_states.size())
            st = basis_states[s];
          else
          {
            st = new SHOState();
            basis_states.add(st);
          }
          RealType e = energy*(n+.5*DIM);
          st->set(qnumber,e);
          s++;
        }
      }
      basis_states.energy_sort(1e-6,true);
    }

    // reset energy scale even if no states need to be added
    for(int i=0;i<basis_states.size();++i)
    {
      SHOState& state = *basis_states[i];
      const TinyVector<int,DIM>& qnumber = state.quantum_number;
      int n = 0;
      for(int d=0;d<DIM;++d)
        n += qnumber[d];
      state.energy = energy*(n+.5*DIM);
    }

    //somewhat redundant, but necessary
    clear_states(0);
    states[0]->finish(basis_states.states);

    if(basis_states.size()<=smax)
      APP_ABORT("SHOSetBuilder::update_basis_states  failed to make enough states");
  }


  void SHOSetBuilder::report(const std::string& pad)
  {
    app_log()<<pad<<"SHOSetBuilder report"<< std::endl;
    app_log()<<pad<<"  dimension = "<< DIM << std::endl;
    app_log()<<pad<<"  mass      = "<< mass << std::endl;
    app_log()<<pad<<"  frequency = "<< energy << std::endl;
    app_log()<<pad<<"  energy    = "<< energy << std::endl;
    app_log()<<pad<<"  length    = "<< length << std::endl;
    app_log()<<pad<<"  center    = "<< center << std::endl;
    app_log()<<pad<<"  nstates   = "<< nstates << std::endl;
    app_log()<<pad<<"  nmax      = "<< nmax << std::endl;
    app_log()<<pad<<"  ind_dims  = "<< ind_dims << std::endl;
    app_log()<<pad<<"  # basis states = "<< basis_states.size() << std::endl;
    app_log()<<pad<<"  basis_states"<< std::endl;
    for(int s=0;s<basis_states.size();++s)
      basis_states[s]->report(pad+"  "+int2string(s)+" ");
    app_log()<<pad<<"end SHOSetBuilder report"<< std::endl;
    app_log().flush();
  }

}
