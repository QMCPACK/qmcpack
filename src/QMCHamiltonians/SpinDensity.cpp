//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include <QMCHamiltonians/SpinDensity.h>
#include <OhmmsData/AttributeSet.h>
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

  SpinDensity::SpinDensity(ParticleSet& P)
  {
    // get particle information
    SpeciesSet& species = P.getSpeciesSet();
    nspecies = species.size();
    int isize  = species.addAttribute("membersize");
    if(isize==species.numAttributes())
      APP_ABORT("SpinDensity(P)  Species set does not have the required attribute 'membersize'");
    for(int s=0;s<nspecies;++s)
      species_size.push_back(species(isize,s));
    for(int s=0;s<nspecies;++s)
      species_name.push_back(species.speciesName[s]);
    reset();

    //jtk: spin density only works for periodic bc's for now
    //     abort if using open boundary conditions
    bool open_bcs=(P.Lattice.SuperCellEnum == SUPERCELL_OPEN);
    if(open_bcs)
    {
      APP_ABORT("SpinDensity is not implemented for open boundary conditions at present\n  please contact the developers if you need this feature");
    }

    Ptmp = &P;
  }


  void SpinDensity::reset()
  {
    myName = "SpinDensity";
    UpdateMode.set(COLLECTABLE,1);
    corner = 0.0;
  }

  
  QMCHamiltonianBase* SpinDensity::makeClone(ParticleSet& P, TrialWaveFunction& Psi)
  {
    return new SpinDensity(*this);
  }
 

  bool SpinDensity::put(xmlNodePtr cur)
  {
    using std::ceil;
    using std::sqrt;
    reset();
    std::string write_report = "no";
    OhmmsAttributeSet attrib;
    attrib.add(myName,"name");
    attrib.add(write_report,"report");
    attrib.put(cur);

    bool have_dr = false;
    bool have_grid = false;
    bool have_center = false;
    bool have_corner = false;
    bool have_cell = false;

    PosType dr;
    PosType center;
    Tensor<RealType,DIM> axes;

    int test_moves = 0;

    xmlNodePtr element = cur->xmlChildrenNode;
    while(element!=NULL)
    {
      std::string ename((const char*)element->name);
      if(ename=="parameter")
      {
        std::string name((const char*)(xmlGetProp(element,(const xmlChar*)"name")));
        if(name=="dr")      
        {
          have_dr = true;
          putContent(dr,element);   
        }
        else if(name=="grid") 
        {
          have_grid = true;
          putContent(grid,element);        
        }
        else if(name=="corner") 
        {
          have_corner = true;
          putContent(corner,element);        
        }
        else if(name=="center") 
        {
          have_center = true;
          putContent(center,element);        
        }
        else if(name=="cell") 
        {
          have_cell = true;
          putContent(axes,element);        
        }
        else if(name=="test_moves") 
          putContent(test_moves,element);        
      }
      element = element->next;
    }

    if(have_dr && have_grid)
    {
      APP_ABORT("SpinDensity::put  dr and grid are provided, this is ambiguous");
    }
    else if(!have_dr && !have_grid)
      APP_ABORT("SpinDensity::put  must provide dr or grid");

    if(have_corner && have_center)
      APP_ABORT("SpinDensity::put  corner and center are provided, this is ambiguous");
    if(have_cell)
    {
      cell.set(axes);
      if(!have_corner && !have_center)
        APP_ABORT("SpinDensity::put  must provide corner or center");
    }
    else
      cell = Ptmp->Lattice;

    if(have_center)
      corner = center-cell.Center;

    if(have_dr)
      for(int d=0;d<DIM;++d)
        grid[d] = (int)ceil(sqrt(dot(cell.Rv[d],cell.Rv[d]))/dr[d]);

    npoints = 1;
    for(int d=0;d<DIM;++d)
      npoints *= grid[d];
    gdims[0] = npoints/grid[0];
    for(int d=1;d<DIM;++d)
      gdims[d] = gdims[d-1]/grid[d];

    if(write_report=="yes")
      report("  ");
    if(test_moves>0)
      test(test_moves,*Ptmp);

    return true;
  }


  void SpinDensity::report(const std::string& pad)
  {
    app_log()<<pad<<"SpinDensity report"<< std::endl;
    app_log()<<pad<<"  dim     = "<< DIM << std::endl;
    app_log()<<pad<<"  npoints = "<< npoints << std::endl;
    app_log()<<pad<<"  grid    = "<< grid << std::endl;
    app_log()<<pad<<"  gdims   = "<< gdims << std::endl;
    app_log()<<pad<<"  corner  = "<< corner << std::endl;
    app_log()<<pad<<"  center  = "<< corner+cell.Center << std::endl;
    app_log()<<pad<<"  cell " << std::endl;
    for(int d=0;d<DIM;++d)
      app_log()<<pad<<"    "<< d <<" "<< cell.Rv[d] << std::endl;
    app_log()<<pad<<"  end cell " << std::endl;
    app_log()<<pad<<"  nspecies = "<< nspecies << std::endl;
    for(int s=0;s<nspecies;++s)
      app_log()<<pad<<"    species["<<s<<"]"<<" = "<<species_name[s]<<" "<<species_size[s]<< std::endl;
    app_log()<<pad<<"end SpinDensity report"<< std::endl;
  }


  void SpinDensity::addObservables(PropertySetType& plist,BufferType& collectables)
  {
    myIndex=collectables.current();
    std::vector<RealType> tmp(nspecies*npoints);
    collectables.add(tmp.begin(),tmp.end());
  }


  void SpinDensity::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const 
  {
    hid_t sgid=H5Gcreate(gid,myName.c_str(),0);

    //vector<int> ng(DIM);
    //for(int d=0;d<DIM;++d)
    //  ng[d] = grid[d];

    std::vector<int> ng(1);
    ng[0] = npoints;

    for(int s=0;s<nspecies;++s)
    {
      observable_helper* oh = new observable_helper(species_name[s]);
      oh->set_dimensions(ng,myIndex+s*npoints);
      oh->open(sgid);
      h5desc.push_back(oh);
    }
  }


  SpinDensity::Return_t SpinDensity::evaluate(ParticleSet& P)
  {
    RealType w=tWalker->Weight;
    int p=0;
    int offset = myIndex;
    for(int s=0;s<nspecies;++s,offset+=npoints)
      for(int ps=0;ps<species_size[s];++ps,++p)
      {
        PosType u = cell.toUnit(P.R[p]-corner);
        //bool inside = true;
        //for(int d=0;d<DIM;++d)
        //  inside &= u[d]>0.0 && u[d]<1.0;
        //if(inside)
        //{
          int point=offset;
          for(int d=0;d<DIM;++d)
            point += gdims[d]*((int)(grid[d]*(u[d]-std::floor(u[d])))); //periodic only
          P.Collectables[point] += w;
        //}
      }
    return 0.0;
  }


  void SpinDensity::test(int moves,ParticleSet& P)
  {
    app_log()<<"  SpinDensity test"<< std::endl;
    RandomGenerator_t rng;
    int particles = P.getTotalNum();
    int pmin = std::numeric_limits<int>::max();
    int pmax = std::numeric_limits<int>::min();
    for(int m=0;m<moves;++m)
    {
      for(int p=0;p<particles;++p)
      {
        PosType u;
        for(int d=0;d<DIM;++d)
          u[d] = rng();
        P.R[p] = P.Lattice.toCart(u);
      }
      test_evaluate(P,pmin,pmax);
    }
    app_log()<<"  end SpinDensity test"<< std::endl;
    APP_ABORT("SpinDensity::test  test complete");
  }



  SpinDensity::Return_t SpinDensity::test_evaluate(ParticleSet& P,int& pmin,int&pmax)
  {
    RealType w=1.0;
    int p=0;
    int offset = 0;
    for(int s=0;s<nspecies;++s,offset+=npoints)
      for(int ps=0;ps<species_size[s];++ps,++p)
      {
        PosType u = cell.toUnit(P.R[p]-corner);
        bool inside = true;
        for(int d=0;d<DIM;++d)
          inside &= u[d]>0.0 && u[d]<1.0;
        if(inside)
        {
          int point=offset;
          for(int d=0;d<DIM;++d)
            point += gdims[d]*((int)(u[d]*grid[d]));
          pmin = std::min(pmin,point-offset);
          pmax = std::max(pmax,point-offset);
        }
      }
    app_log()<<"    pmin = "<<pmin<<" pmax = "<<pmax<<" npoints = "<<npoints<< std::endl;
    return 0.0;
  }


  void SpinDensity::
  postprocess_density(const std::string& infile,const std::string& species,
                      pts_t& points,dens_t& density,dens_t& density_err)
  {
    std::ifstream datafile;
    datafile.open(infile.c_str());
    if(!datafile.is_open())
      APP_ABORT("SpinDensity::postprocess_density\n  could not open file: "+infile);

    int n=0;
    dens_t dens;
    dens_t dens_err;
    RealType value;    
    while(datafile>>value)
    {
      if(n<npoints)
        dens.push_back(value);
      else
        dens_err.push_back(value);
      n++;
    }
    if(dens_err.size()!=npoints)
    {
      app_log()<<"SpinDensity::postprocess_density\n  file "<<infile<<"\n  contains "<<dens.size()+dens_err.size()<<" values\n  expected "<<2*npoints<<" values"<< std::endl;
      APP_ABORT("SpinDensity::postprocess_density");
    }

    //evaluate the density
    //  this error analysis neglects spatial (r,r') correlations

    //vector<int> dc;
    //dc.resize(npoints);
    //for(int p=0;p<npoints;++p)
    //  dc[p]=0;

    for(int p=0;p<points.size();++p)
    {
      PosType u = cell.toUnit(points[p]-corner);
      bool inside = true;
      for(int d=0;d<DIM;++d)
        inside &= u[d]>0.0 && u[d]<1.0;
      if(inside)
      {
        int point=0;
        for(int d=0;d<DIM;++d)
          point += gdims[d]*((int)(u[d]*grid[d]));
        density[p]     = dens[point];
        density_err[p] = dens_err[point];
        
        //dc[point]+=1;
      }
      else
      {
        density[p] = 0.0;
        density_err[p] = 0.0;
      }
    }
    
    //vector<int> count(20,0);
    //for(int p=0;p<npoints;++p)
    //  count[dc[p]]++;
    //for(int c=0;c<count.size();++c)
    //  app_log()<<"count "<<c<<" "<<count[c]<< std::endl;
    //
    //PosType pmin;
    //PosType pmax;
    //pmin =  1e99;
    //pmax = -1e99;
    //for(int p=0;p<points.size();++p)
    //  for(int d=0;d<DIM;++d)
    //  {
    //    pmin[d] = std::min(pmin[d],points[p][d]);
    //    pmax[d] = std::max(pmax[d],points[p][d]);
    //  }
    //app_log()<<" pmin "<<pmin<< std::endl;
    //app_log()<<" pmax "<<pmax<< std::endl;
    //app_log()<<" gridpoints = "<<points.size()<< std::endl;
    //
    //APP_ABORT("sdens check");

  }
  
  void SpinDensity::addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy)
  {
    int nw = W.WalkerList.size();
    for (int iw=0; iw<nw; iw++)
    {
      Walker_t &w = *W.WalkerList[iw];
      RealType weight=w.Weight/nw;
      int p=0;
      int offset = myIndex;
      for(int s=0; s<nspecies; ++s,offset+=npoints)
        for(int ps=0; ps<species_size[s]; ++ps,++p)
        {
          PosType u = cell.toUnit(w.R[p]-corner);
          //bool inside = true;
          //for(int d=0;d<DIM;++d)
          //  inside &= u[d]>0.0 && u[d]<1.0;
          //if(inside)
          //{
          int point=offset;
          for(int d=0;d<DIM;++d)
            point += gdims[d]*((int)(grid[d]*(u[d]-std::floor(u[d])))); //periodic only
          W.Collectables[point] += weight;
          //}
        }
    }
  }
}
