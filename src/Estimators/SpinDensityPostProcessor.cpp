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
    
    



#include <Utilities/string_utils.h>
#include <Utilities/unit_conversion.h>
#include <QMCHamiltonians/QMCHamiltonian.h>
#include <QMCHamiltonians/SpinDensity.h>
#include <QMCHamiltonians/StaticStructureFactor.h>
#include <QMCHamiltonians/DensityMatrices1B.h>
#include <Estimators/SpinDensityPostProcessor.h>

namespace qmcplusplus
{
  SpinDensityPostProcessor::
  SpinDensityPostProcessor(ParticleSet& pq,ParticleSet& pc,QMCHamiltonian& h)
    : Pq(pq),Pc(pc),H(h)
  {
  }


  void SpinDensityPostProcessor::put(xmlNodePtr cur)
  {
    using std::sqrt;
    using std::ceil;

    bool have_dr = false;
    bool have_grid = false;
    bool have_center = false;
    bool have_corner = false;
    bool have_cell = false;
    bool have_norm = false;

    TinyVector<int,DIM> gdims;
    PosType du;
    PosType dr;
    PosType center;
    Tensor<RealType,DIM> axes;
    format = "";
    normalization = "sum";

    xmlNodePtr element = cur->xmlChildrenNode;
    while(element!=NULL)
    {
      std::string ename((const char*)element->name);
      if(ename=="parameter")
      {
        std::string name((const char*)(xmlGetProp(element,(const xmlChar*)"name")));
        if(name=="sources")
        {
          putContent(sources,element);
        }
        else if(name=="format")
        {
          putContent(format,element);
        }
        else if(name=="dr")      
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
        else if(name=="norm") 
        {
          have_norm = true;
          putContent(normalization,element);        
        }
        else
          APP_ABORT("SpinDensityPostProcessor::put  "+name+" is not a valid parameter name\n  valid options are: sources, format, dr, grid, corner, center, cell");
      }
      element = element->next;
    }

    if(format=="")
      APP_ABORT("SpinDensityPostProcessor::put  format must be provided");

    if(normalization!="sum" && normalization!="max")
      APP_ABORT("SpinDensityPostProcessor::put  norm ("+normalization+") is invalid\n  valid options are: sum,max");
    
    if(have_dr && have_grid)
    {
      APP_ABORT("SpinDensityPostProcessor::put  dr and grid are provided, this is ambiguous");
    }
    else if(!have_dr && !have_grid)
      APP_ABORT("SpinDensityPostProcessor::put  must provide dr or grid");

    if(have_corner && have_center)
      APP_ABORT("SpinDensityPostProcessor::put  corner and center are provided, this is ambiguous");
    if(have_cell)
    {
      cell.set(axes);
      if(!have_corner && !have_center)
        APP_ABORT("SpinDensityPostProcessor::put  must provide corner or center");
    }
    else
      cell = Pq.Lattice;

    if(have_center)
      corner = center-cell.Center;

    if(have_dr)
      for(int d=0;d<DIM;++d)
        grid[d] = (int)ceil(sqrt(dot(cell.Rv[d],cell.Rv[d]))/dr[d]);

    int ncells = 1;
    for(int d=0;d<DIM;++d)
    {
      ncells *= grid[d];
      du[d] = 1.0/grid[d];
    }
    dV = cell.Volume/ncells;

    //make it a 'general' xcrysden grid
    // must have points covering all facets of the cell
    npoints = 1;
    for(int d=0;d<DIM;++d)
    {
      grid[d] += 1;
      npoints *= grid[d];
    }

    gdims[0] = 1;
    for(int d=1;d<DIM;++d)
      gdims[d] = gdims[d-1]*grid[d-1];

    gridpoints.resize(npoints);
    PosType u;
    for(int p=0;p<npoints;++p)
    {
      int nrem = p;
      for(int d=DIM-1;d>0;--d)
      {
        int ind = nrem/gdims[d];
        u[d] = ind*du[d];
        nrem-= ind*gdims[d];
      }
      u[0] = nrem*du[0];
      gridpoints[p] = cell.toCart(u) + corner;
    }

    SpeciesSet& species = Pq.getSpeciesSet();
    nspecies = species.size();
    for(int s=0;s<nspecies;++s)
      species_name.push_back(species.speciesName[s]);
    int isize  = species.addAttribute("membersize");
    if(isize==species.numAttributes())
      APP_ABORT("DensityMatrices1B::set_state  Species set does not have the required attribute 'membersize'");
    for(int s=0;s<nspecies;++s)
      species_size.push_back(species(isize,s));

    app_log()<<"      cell origin = "<<corner<< std::endl;
    for(int d=0;d<DIM;++d)
      app_log()<<"      cell axes "<<d<<" = "<<cell.Rv[d]<< std::endl;
    app_log()<<"      grid shape  = "<<grid<< std::endl;
    app_log()<<"      grid size   = "<<npoints<< std::endl;
    app_log()<<"      sources     = ";
    for(int i=0;i<sources.size();++i)
      app_log()<<sources[i]<<" ";
    app_log()<< std::endl;
    app_log()<<"      species     = ";
    for(int i=0;i<species_name.size();++i)
      app_log()<<species_name[i]<<" ";
    app_log()<< std::endl;
  }


  void SpinDensityPostProcessor::postprocess()
  {
    std::map<std::string,std::string> typemap;
    typemap["spindensity"]     = "sd";
    typemap["structurefactor"] = "sf";
    typemap["dm1b"]            = "dm";

    dens_t density;
    dens_t density_err;
    density.resize(npoints);
    density_err.resize(npoints);

    for(int series=series_start;series<series_end;++series)
    {
      char cseries[256];
      sprintf(cseries,"s%03d.",series);
      std::string sprefix(cseries);
      sprefix = project_id+"."+sprefix;

      for(int isrc=0;isrc<sources.size();++isrc)
      {
        std::string name           = sources[isrc];
        QMCHamiltonianBase* h = H.getHamiltonian(name);
        std::string type           = H.getOperatorType(name);
        if(h==0)
          APP_ABORT("SpinDensityPostProcessor::postprocess\n  could not find operator with name "+name);

        for(int ispec=0;ispec<nspecies;++ispec)
        {
          const std::string& species = species_name[ispec];
          std::string infile  = sprefix+type+"_"+species+".dat";
          std::string outfile = sprefix+"density_"+typemap[type]+"_"+species;

          app_log()<<"      evaluating & writing density for "<<outfile<< std::endl;

          if(type=="spindensity")
            get_density<SpinDensity>(infile,species,h,density,density_err);
          else if(type=="structurefactor")
            get_density<StaticStructureFactor>(infile,species,h,density,density_err);
          else if(type=="dm1b")
            get_density<DensityMatrices1B>(infile,species,h,density,density_err);            
          else
            APP_ABORT("SpinDensityPostProcessor::postprocess\n  "+name+" of type "+type+" is not a valid density source");

          normalize(species_size[ispec],density,density_err);

          // mean of density field
          write_density(outfile,density);
          
          // mean + error of density field
          for(int p=0;p<npoints;++p)
            density[p]+=density_err[p];
          write_density(outfile+"+err",density);
          
          // mean - error of density field
          for(int p=0;p<npoints;++p)
            density[p]-=2.0*density_err[p];
          write_density(outfile+"-err",density);
        }
      }
    }
  }

  
  template<typename SDO>
  void SpinDensityPostProcessor::get_density(
    const std::string& infile,const std::string& species,
    QMCHamiltonianBase* h,dens_t& density,dens_t& density_err)
  {
    SDO* sdo = dynamic_cast<SDO*>(h);
    if(sdo==0)
      APP_ABORT("SpinDensityPostProcessor::get_density  dynamic cast failed");
    sdo->postprocess_density(infile,species,gridpoints,density,density_err);
  }


  void SpinDensityPostProcessor::normalize(int nparticles,dens_t& density,dens_t& density_err)
  {
    RealType norm = 1.0;
    if(normalization=="sum")
    {
      RealType dsum = 0.0;
      for(int p=0;p<npoints;++p)
        dsum += density[p];
      if(dsum>1e-12)
        norm = nparticles/(dsum*cell.Volume/npoints);
    }
    else if(normalization=="max")
    {
      RealType dmax = -std::numeric_limits<RealType>::max();
      for(int p=0;p<npoints;++p)
        dmax = std::max(dmax,std::abs(density[p]));
      norm = 1.0/dmax;
    }
    for(int p=0;p<npoints;++p)
      density[p] *= norm;
    for(int p=0;p<npoints;++p)
      density_err[p] *= norm;
  }


  void SpinDensityPostProcessor::write_density(const std::string& outfile,dens_t& density)
  {
    if(format=="xsf")
      write_density_xsf(outfile,density);
    else
      APP_ABORT("SpinDensityPostProcessor::write_density\n  file format "+format+" is unknown");
  }


  void SpinDensityPostProcessor::write_density_xsf(const std::string& outfile,dens_t& density)
  {
    using Units::convert;
    using Units::B;
    using Units::A;

    std::ofstream file;
    std::string filename = outfile+".xsf";
    file.open(filename.c_str(),std::ios::out | std::ios::trunc);
    if(!file.is_open())
      APP_ABORT("SpinDensityPostProcessor::write_density_xsf\n  failed to open file for output: "+outfile+".xsf");

    file.precision(6);
    file<<std::scientific;
    int columns = 5;

    int natoms = Pc.getTotalNum();

    file<<" CRYSTAL"<< std::endl;
    file<<" PRIMVEC"<< std::endl;
    for(int i=0;i<DIM;++i)
    {
      file<<" ";
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(Pq.Lattice.Rv[i][d],B,A);
      file<< std::endl;
    }
    file<<" PRIMCOORD"<< std::endl;
    file<<"   "<<natoms<<"   1"<< std::endl;
    for(int i=0;i<natoms;++i)
    {
      file<<"   "<<pname(Pc,i);
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(Pc.R[i][d],B,A);
      file<< std::endl;
    }
    file<<" BEGIN_BLOCK_DATAGRID_3D"<< std::endl;
    file<<"   "<<outfile<< std::endl;
    file<<"   DATAGRID_3D_SPIN_DENSITY"<< std::endl;
    file<<"   ";
    for(int d=0;d<DIM;++d)
      file<<"  "<<grid[d];
    file<< std::endl;
    file<<"   ";
    for(int d=0;d<DIM;++d)
      file<<"  "<<convert(corner[d],B,A);
    file<< std::endl;
    for(int i=0;i<DIM;++i)
    {
      file<<"   ";
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(cell.Rv[i][d],B,A);
      file<< std::endl;
    }
    file<<"   ";
    for(int p=0;p<npoints;++p)
    {
      file<<"  "<<density[p];
      if((p+1)%columns==0)
        file<< std::endl<<"   ";
    }
    file<< std::endl;
    file<<"   END_DATAGRID_3D"<< std::endl;
    file<<" END_BLOCK_DATAGRID_3D"<< std::endl;
  }
}
