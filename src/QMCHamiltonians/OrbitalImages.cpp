//////////////////////////////////////////////////////////////////
// (c) Copyright 2014-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCHamiltonians/OrbitalImages.h>
#include <OhmmsData/AttributeSet.h>
#include <QMCWaveFunctions/BasisSetFactory.h>
#include <Utilities/unit_conversion.h>


namespace qmcplusplus
{

  OrbitalImages::OrbitalImages(ParticleSet& P,PSPool& PSP,Communicate* mpicomm)
    : psetpool(PSP)
  {
    //keep the electron particle to get the cell later, if necessary
    Peln = &P;

    //keep the communicator to select master task for file write
    comm = mpicomm;
  }

  
  QMCHamiltonianBase* OrbitalImages::makeClone(ParticleSet& P, TrialWaveFunction& Psi)
  {
    //cloning shouldn't strictly be necessary, but do it right just in case
    OrbitalImages* clone = new OrbitalImages(*this);
    clone->Peln = &P;
    for(int i=0;i<sposets.size();++i)
    {
      clone->sposet_indices[i] = new vector<int>(*sposet_indices[i]);
      clone->sposets[i]        = sposets[i]->makeClone();
    }
    return clone;
  }
 

  bool OrbitalImages::put(xmlNodePtr cur)
  {
    app_log()<<"OrbitalImages::put"<<endl;

    //set defaults
    myName = "OrbitalImages";
    corner = 0.0;
    batch_size = -1;

    //read simple attributes
    string write_report = "yes";
    string ion_psname = "ion0";
    OhmmsAttributeSet attrib;
    attrib.add(myName,"name");
    attrib.add(ion_psname,"ions");
    attrib.add(write_report,"report");
    attrib.put(cur);

    //read parameters
    bool have_grid = false;
    bool have_center = false;
    bool have_corner = false;
    bool have_cell = false;

    PosType center;
    Tensor<RealType,DIM> axes;
    vector<string> valtypes;
    string file_format = "xsf";

    xmlNodePtr element = cur->xmlChildrenNode;
    vector<xmlNodePtr> other_elements;
    while(element!=NULL)
    {
      string ename((const char*)element->name);
      if(ename=="parameter")
      {
        string name((const char*)(xmlGetProp(element,(const xmlChar*)"name")));
        if(name=="sposets")
          putContent(sposet_names,element);   
        else if(name=="batch_size") 
          putContent(batch_size,element);        
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
        else if(name=="value") 
          putContent(valtypes,element);        
        else if(name=="format") 
          putContent(file_format,element);        
        else
          other_elements.push_back(element);
      }
      element = element->next;
    }

    //second pass parameter read to get orbital indices
    //  each parameter is named after the corresponding sposet
    for(int i=0;i<sposet_names.size();++i)
      sposet_indices.push_back(new vector<int>);
    for(int n=0;n<other_elements.size();++n)
    {
      xmlNodePtr element = other_elements[n];
      string ename((const char*)element->name);
      if(ename=="parameter")
      {
        string name((const char*)(xmlGetProp(element,(const xmlChar*)"name")));
        for(int i=0;i<sposet_names.size();++i)
          if(name==sposet_names[i])
            putContent(*sposet_indices[i],element);
      }
    }
    
    //check that the value types requested are reasonable
    for(int i=0;i<valtypes.size();++i)
    {
      const string& valtype = valtypes[i];
      value_types_enum value_type;
      if(valtype=="real")
        value_type = real_val;
      else if(valtype=="imag")
        value_type = imag_val;
      else if(valtype=="abs")
      value_type = abs_val;
      else if(valtype=="abs2")
        value_type = abs2_val;
      else
      {
        APP_ABORT("OrbitalImages::put  value type "+valtype+" is unsupported\n  valid options are: value, abs, abs2");
      }
      value_types.push_back(value_type);
    }
    if(value_types.size()==0)
      value_types.push_back(real_val);

    //check the format
    if(file_format=="xsf")
      format = xsf;
    else
    {
      APP_ABORT("OrbitalImages::put  file format "+file_format+" is invalid\n  valid options are: xsf");
    }

    //get the ion particleset
    if(psetpool.find(ion_psname)==psetpool.end())
    {
      APP_ABORT("OrbitalImages::put  ParticleSet "+ion_psname+" does not exist");
    }
    Pion = psetpool[ion_psname];

    app_log()<<"  getting sposets"<<endl;

    // get the sposets for image output
    if(sposet_names.size()==0)
      APP_ABORT("OrbitalImages::put  must have at least one sposet");
    for(int i=0;i<sposet_names.size();++i)
    {
      SPOSetBase* sposet = get_sposet(sposet_names[i]);
      if(sposet==0)
        APP_ABORT("OrbitalImages::put  sposet "+sposet_names[i]+" does not exist");
      sposets.push_back(sposet);
      vector<int>& sposet_inds = *sposet_indices[i];
      if(sposet_inds.size()==0)
        for(int n=0;n<sposet->size();++n)
          sposet_inds.push_back(n);
      for(int n=0;n<sposet_inds.size();++n)
      {
        int index = sposet_inds[n];
        if(index<0 || index>sposet->size())
        {
          app_log()<<"\nindex for sposet "<<sposet_names[i]<<" is out of range\nindex must be between 0 and "<<sposet->size()-1<<"\nyou provided: "<<index<<endl;
          APP_ABORT("OrbitalImages::put  sposet index out of range, see message above");
        }
      }
    }

    //check that the grid and cell are properly defined
    if(!have_grid)
      APP_ABORT("OrbitalImages::put  must provide grid");
    if(have_corner && have_center)
      APP_ABORT("OrbitalImages::put  corner and center are provided, this is ambiguous");
    if(have_cell)
    {
      cell.set(axes);
      if(!have_corner && !have_center)
        APP_ABORT("OrbitalImages::put  must provide corner or center");
    }
    else
      cell = Peln->Lattice;

    //calculate the cell corner in the case that the cell center is provided
    if(have_center)
      corner = center-cell.Center;

    //get the number of grid points
    npoints = 1;
    for(int d=0;d<DIM;++d)
      npoints *= grid[d];

    //get grid strides in D dimensions
    gdims[0] = 1;
    for(int d=1;d<DIM;++d)
      gdims[d] = gdims[d-1]*grid[d-1];

    //write a brief report to the log file
    if(write_report=="yes")
      report("  ");

    app_log()<<"end OrbitalImages::put"<<endl;
    return true;
  }


  void OrbitalImages::report(const string& pad)
  {
    app_log()<<pad<<"OrbitalImages report"<<endl;
    app_log()<<pad<<"  nsposets = "<<sposets.size()<<" "<<sposet_names.size()<<" "<<sposet_indices.size()<<endl;
    for(int i=0;i<sposet_names.size();++i)
    {
      vector<int>& sposet_inds = *sposet_indices[i];
      SPOSetBase& sposet = *sposets[i];
      if(sposet_inds.size()==sposet.size())
        app_log()<<pad<<"  "<<sposet_names[i]<<" = all "<<sposet.size()<<" orbitals"<<endl;
      else
      {
        app_log()<<pad<<"  "<<sposet_names[i]<<" =";
        for(int n=0;n<sposet_inds.size();++n)
          app_log()<<" "<<sposet_inds[n];
        app_log()<<endl;
      }
    }
    app_log()<<pad<<"  npoints = "<< npoints <<endl;
    app_log()<<pad<<"  grid    = "<< grid <<endl;
    app_log()<<pad<<"  corner  = "<< corner <<endl;
    app_log()<<pad<<"  center  = "<< corner+cell.Center <<endl;
    app_log()<<pad<<"  cell " <<endl;
    for(int d=0;d<DIM;++d)
      app_log()<<pad<<"    "<< d <<" "<< cell.Rv[d] <<endl;
    app_log()<<pad<<"  end cell " <<endl;
    app_log()<<pad<<"end OrbitalImages report"<<endl;
  }




  OrbitalImages::Return_t OrbitalImages::evaluate(ParticleSet& P)
  {
    //only the first thread of the master task writes the orbitals
    if(comm->rank()==0 && omp_get_thread_num()==0)
    {
      app_log()<<endl;
      app_log()<<"OrbitalImages::evaluate  writing orbital images"<<endl;

      //create the grid points by mapping point index into D-dim points
      app_log()<<"  generating grid "<<grid<<endl;
      vector<PosType> rpoints;
      rpoints.resize(npoints);
      PosType u,du;
      for(int d=0;d<DIM;++d)
        du[d] += 1.0/grid[d];
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
        rpoints[p] = cell.toCart(u) + corner;
      }

      //evaluate and write the orbitals for each sposet
      app_log()<<"  evaluating all orbitals"<<endl;
      for(int i=0;i<sposets.size();++i)
      {
        //get sposet information
        const string& sposet_name = sposet_names[i];
        app_log()<<"  evaluating orbitals in "+sposet_name+" on the grid"<<endl;
        vector<int>& sposet_inds = *sposet_indices[i];
        SPOSetBase& sposet = *sposets[i];
        int nspo = sposet_inds.size();
        
        //set the batch size
        int bsize = batch_size;
        if(bsize==-1)
          bsize = nspo;

        //resize temporary evaluation arrays
        spo_vtmp.resize(sposet.size());
        batch_values.resize(npoints,bsize);
        orbital.resize(npoints);

        //loop over orbitals one batch at a time (batch_size orbitals)
        int bstart=0;
        int bend=min(bsize,nspo);
        while(bstart<nspo)
        {
          //fill up the temporary storage for a batch of orbitals
          app_log()<<"    evaluating orbital batch "<<bstart<<" to "<<bend<<" out of "<<nspo<<endl;
          for(int p=0;p<npoints;++p)
          {
            P.makeMove(0,rpoints[p]-P.R[0]);
            sposet.evaluate(P,0,spo_vtmp); //note that ALL orbitals are evaluated each time
            P.rejectMove(0);
            for(int b=bstart,ib=0;b<bend;++b,++ib)
              batch_values(p,ib) = spo_vtmp[sposet_inds[b]];
          }
          //write out the batch one orbital at a time for each value type requested
          app_log()<<"    writing all orbitals in the batch"<<endl;
          for(int b=bstart,ib=0;b<bend;++b,++ib)
          {
            //transfer strided orbital info into contiguous array
            for(int p=0;p<npoints;++p)
              orbital[p] = batch_values(p,ib);
            //write one file for each value type (real, imag, etc) selected
            for(int iv=0;iv<value_types.size();++iv)
              write_orbital(sposet_name,sposet_inds[b],orbital,value_types[iv]);
          }
          bstart = bend;
          bend += bsize;
          bend = min(bend,nspo);
        }
      }

      app_log()<<"OrbitalImages::evaluate  orbital images written successfully, exiting.\n"<<endl;
      APP_ABORT("  Not a fatal error, just exiting.");
    }

    //make sure no other process runs off
    comm->barrier();

    return 0.0;
  }


  void OrbitalImages::write_orbital(const string& sponame,int index,vector<ValueType>& orb,value_types_enum value_type)
  {
    //only xsf format is supported for now
    if(format==xsf)
      write_orbital_xsf(sponame,index,orb,value_type);
  }


  void OrbitalImages::write_orbital_xsf(const string& sponame,int index,vector<ValueType>& orb,value_types_enum value_type)
  {
    using Units::convert;
    using Units::B;
    using Units::A;

    //generate file name
    char filename[100];
    if(value_type==real_val)
      sprintf(filename,"%s_orbital_%04d.xsf",sponame.c_str(),index);
    else if(value_type==imag_val)
      sprintf(filename,"%s_orbital_%04d_imag.xsf",sponame.c_str(),index);
    else if(value_type==abs_val)
      sprintf(filename,"%s_orbital_%04d_abs.xsf",sponame.c_str(),index);
    else if(value_type==abs2_val)
      sprintf(filename,"%s_orbital_%04d_abs2.xsf",sponame.c_str(),index);

    app_log()<<"      writing file: "<<string(filename)<<endl;

    //get the cell containing the ion positions
    //  assume the evaluation cell if any boundaries are open
    ParticleSet& Pc = *Pion;
    Lattice_t* Lbox;
    if(Peln->Lattice.SuperCellEnum==SUPERCELL_BULK)
      Lbox = &Peln->Lattice; //periodic
    else
      Lbox = &cell; //at least partially open

    //open the file
    ofstream file;
    file.open(filename,ios::out | ios::trunc);
    if(!file.is_open())
      APP_ABORT("OrbitalImages::write_orbital\n  failed to open file for output: "+string(filename));

    //set the precision & number of columns
    file.precision(6);
    file<<std::scientific;
    int columns = 5;

    //get the number of atoms
    int natoms = Pc.getTotalNum();

    //write the cell, atomic positions, and orbital to the file in xsf format
    file<<" CRYSTAL"<<endl;
    file<<" PRIMVEC"<<endl;
    for(int i=0;i<DIM;++i)
    {
      file<<" ";
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(Lbox->Rv[i][d],B,A);
      file<<endl;
    }
    file<<" PRIMCOORD"<<endl;
    file<<"   "<<natoms<<"   1"<<endl;
    for(int i=0;i<natoms;++i)
    {
      file<<"   "<<Pc.species_from_index(i);
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(Pc.R[i][d],B,A);
      file<<endl;
    }
    file<<" BEGIN_BLOCK_DATAGRID_3D"<<endl;
    file<<"   orbital "<<index<<endl;
    file<<"   DATAGRID_3D_ORBITAL"<<endl;
    file<<"   ";
    for(int d=0;d<DIM;++d)
      file<<"  "<<grid[d];
    file<<endl;
    file<<"   ";
    for(int d=0;d<DIM;++d)
      file<<"  "<<convert(corner[d],B,A);
    file<<endl;
    for(int i=0;i<DIM;++i)
    {
      file<<"   ";
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(cell.Rv[i][d],B,A);
      file<<endl;
    }
    file<<"   ";
    if(value_type==real_val)      // real part
      for(int p=0;p<npoints;++p)
      {
        file<<"  "<<real(orb[p]);
        if((p+1)%columns==0)
          file<<endl<<"   ";
      }
    else if(value_type==imag_val) // imag part
      for(int p=0;p<npoints;++p)
      {
        file<<"  "<<imag(orb[p]);
        if((p+1)%columns==0)
          file<<endl<<"   ";
      }
    else if(value_type==abs_val)  // modulus
      for(int p=0;p<npoints;++p)
      {
        file<<"  "<<abs(orb[p]);
        if((p+1)%columns==0)
          file<<endl<<"   ";
      }
    else if(value_type==abs2_val) // modulus^2
      for(int p=0;p<npoints;++p)
      {
        file<<"  "<<abs(orb[p])*abs(orb[p]);
        if((p+1)%columns==0)
          file<<endl<<"   ";
      }
    file<<endl;
    file<<"   END_DATAGRID_3D"<<endl;
    file<<" END_BLOCK_DATAGRID_3D"<<endl;
  }

}
