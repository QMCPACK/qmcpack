//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCHamiltonians/StaticStructureFactor.h>
#include <OhmmsData/AttributeSet.h>
#include <LongRange/KContainer.h>
#include <LongRange/StructFact.h>

namespace qmcplusplus
{

  StaticStructureFactor::StaticStructureFactor(ParticleSet& P)
      : Pinit(P)
  {
#ifndef USE_REAL_STRUCT_FACTOR
    APP_ABORT("StaticStructureFactor: please recompile with USE_REAL_STRUCT_FACTOR=1");
#endif
    if(P.Lattice.SuperCellEnum==SUPERCELL_OPEN)
      APP_ABORT("StaticStructureFactor is incompatible with open boundary conditions");

    // get particle information
    SpeciesSet& species = P.getSpeciesSet();
    nspecies = species.size();
    for(int s=0;s<nspecies;++s)
      species_name.push_back(species.speciesName[s]);
    reset();
  }


  void StaticStructureFactor::reset()
  {
    myName   = "StaticStructureFactor";
    UpdateMode.set(COLLECTABLE,1);
    ecut     = -1.0;
    nkpoints = -1;
  }

  
  QMCHamiltonianBase* StaticStructureFactor::makeClone(ParticleSet& P, TrialWaveFunction& Psi)
  {
    return new StaticStructureFactor(*this);
  }
 

  bool StaticStructureFactor::put(xmlNodePtr cur)
  {
    using std::sqrt;
    reset();
    const k2_t& k2_init = Pinit.SK->KLists.ksq;

    string write_report = "no";
    OhmmsAttributeSet attrib;
    attrib.add(myName,"name");
    //attrib.add(ecut,"ecut");
    attrib.add(write_report,"report");
    attrib.put(cur);

    if(ecut<0.0)
      nkpoints = k2_init.size();
    else
    {
      RealType k2cut = 2.0*ecut;
      nkpoints = 0;
      for(int i=0;i<k2_init.size();++i)
        if(k2_init[i]<k2cut)
          nkpoints++;
    }
    
    if(nkpoints==0)
      APP_ABORT("StaticStructureFactor::put  could not find any kpoints");
    
    ecut = .5*k2_init[nkpoints-1];

    if(write_report=="yes")
      report("  ");

    return true;
  }


  void StaticStructureFactor::report(const string& pad)
  {
    app_log()<<pad<<"StaticStructureFactor report"<<endl;
    app_log()<<pad<<"  name     = "<< myName   <<endl;
    app_log()<<pad<<"  ecut     = "<< ecut     <<endl;
    app_log()<<pad<<"  nkpoints = "<< nkpoints <<endl;
    app_log()<<pad<<"  nspecies = "<< nspecies <<endl;
    for(int s=0;s<nspecies;++s)
      app_log()<<pad<<"    species["<<s<<"] = "<< species_name[s] <<endl;
    app_log()<<pad<<"end StaticStructureFactor report"<<endl;
  }


  void StaticStructureFactor::addObservables(PropertySetType& plist,BufferType& collectables)
  {
    myIndex=collectables.current();
    vector<RealType> tmp(nspecies*2*nkpoints); // real & imag parts
    collectables.add(tmp.begin(),tmp.end());
  }


  void StaticStructureFactor::registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const 
  {
    hid_t sgid=H5Gcreate(gid,myName.c_str(),0);
    vector<int> ng(2);
    ng[0] = 2;
    ng[1] = nkpoints;
    for(int s=0;s<nspecies;++s)
    {
      observable_helper* oh = new observable_helper(species_name[s]);
      oh->set_dimensions(ng,myIndex+s*2*nkpoints);
      oh->open(sgid);
      h5desc.push_back(oh);
    }
  }


  StaticStructureFactor::Return_t StaticStructureFactor::evaluate(ParticleSet& P)
  {
    RealType w=tWalker->Weight;
    const Matrix<RealType>& rhok_r = P.SK->rhok_r;
    const Matrix<RealType>& rhok_i = P.SK->rhok_i;
    int nkptot = rhok_r.cols();
    for(int s=0;s<nspecies;++s)
    {
      int kc      = myIndex + s*2*nkpoints;
      //int kstart  = s*nkptot;
      //for(int k=kstart;k<kstart+nkpoints;++k,++kc)
      //  P.Collectables[kc] += w*rhok_r(k);
      //for(int k=kstart;k<kstart+nkpoints;++k,++kc)
      //  P.Collectables[kc] += w*rhok_i(k);
      for(int k=0;k<nkpoints;++k,++kc)
        P.Collectables[kc] += w*rhok_r(s,k);
      for(int k=0;k<nkpoints;++k,++kc)
        P.Collectables[kc] += w*rhok_i(s,k);
    }
    return 0.0;
  }


  void StaticStructureFactor::
  postprocess_density(const string& infile,const string& species,
                      pts_t& points,dens_t& density,dens_t& density_err)
  {
    ifstream datafile;
    datafile.open(infile.c_str());
    if(!datafile.is_open())
      APP_ABORT("StaticStructureFactor::postprocess_density\n  could not open file: "+infile);

    const int nk = nkpoints;
    int n=0;
    vector<RealType> skr;
    vector<RealType> ski;
    vector<RealType> skrv;
    vector<RealType> skiv;
    RealType value;    
    while(datafile>>value)
    {
      if(n<nk)
        skr.push_back(value);
      else if(n<2*nk)
        ski.push_back(value);
      else if(n<3*nk)
        skrv.push_back(value*value);
      else if(n<4*nk)
        skiv.push_back(value*value);
      n++;
    }
    if(skiv.size()!=nkpoints)
    {
      app_log()<<"StaticStructureFactor::postprocess_density\n  file "<<infile<<"\n  contains "<< n <<" values\n  expected "<<4*nkpoints<<" values"<<endl;
      APP_ABORT("StaticStructureFactor::postprocess_density");
    }

    //evaluate the density
    //  this error analysis neglects spatial (k,k') correlations
    using std::cos;
    using std::sin;
    using std::sqrt;
    const vector<PosType>& kpoints = Pinit.SK->KLists.kpts_cart;
    for(int p=0;p<points.size();++p)
    {
      RealType d  = 0.0;
      RealType de = 0.0;
      PosType& r  = points[p];
      for(int k=0;k<nkpoints;++k)
      {
        RealType kr = dot(kpoints[k],r);
        RealType cr = cos(kr);
        RealType sr = sin(kr);
        d  +=    cr*skr[k]  +    sr*ski[k];
        de += cr*cr*skrv[k] + sr*sr*skiv[k];
      }
      de = sqrt(de);
      density[p]     = d;
      density_err[p] = de;
    }
    
  }

}
