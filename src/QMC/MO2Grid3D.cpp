//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Numerics/Spline3D/Config.h"
#include "QMC/MO2Grid3D.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleIO/XMLParticleIO.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/LCOrbitals.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFTriCubicSpline.h"
#include "Utilities/Clock.h"
#include "QMC/QMCUtilities.h"

namespace ohmmsqmc {

  MO2Grid3D::MO2Grid3D(int argc, char** argv): QMCApps(argc,argv) { 
    el.setName("e");
    int iu = el.Species.addSpecies("u");
    int id = el.Species.addSpecies("d");
    int icharge = el.Species.addAttribute("charge");
    el.Species(icharge,iu) = -1;
    el.Species(icharge,id) = -1;
  }

  ///destructor
  MO2Grid3D::~MO2Grid3D() {
    DEBUGMSG("MO2Grid3D::~MO2Grid3D")
  }

  bool MO2Grid3D::init() {

    if(!setParticleSets(m_root)) {
       ERRORMSG("Failed to initialize the ions and electrons. Exit now.")
       return false;
    }

    setWavefunctions(m_root);

    string fname(myProject.CurrentRoot());
    fname.append(".debug");
    ofstream fout(fname.c_str());
    fout << "Ionic configuration : " << ion.getName() << endl;
    ion.get(fout);
    fout << "Electronic configuration : " << el.getName() << endl;
    el.get(fout);
    Psi.VarList.print(fout);

    //return false to stop it
    return false;    
  }   

  bool  MO2Grid3D::setParticleSets(xmlNodePtr aroot) {

    bool init_els = determineNumOfElectrons(el,m_context);

    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"//particleset",m_context);

    xmlNodePtr el_ptr=NULL, ion_ptr=NULL;
    for(int i=0; i<result->nodesetval->nodeNr; i++) {
      xmlNodePtr cur=result->nodesetval->nodeTab[i];
      xmlChar* aname= xmlGetProp(cur,(const xmlChar*)"name");
      if(aname) {
	char fc = aname[0];
	if(fc == 'e') { el_ptr=cur;}
	else if(fc == 'i') {ion_ptr=cur;}
      }
    }

    bool donotresize = false;
    if(init_els) {
      el.setName("e");
      XMLReport("The configuration for electrons is already determined by the wave function")
      donotresize = true;
    } 

    if(ion_ptr) {
      XMLParticleParser pread(ion);
      pread.put(ion_ptr);
    }

    if(el_ptr) {
      XMLParticleParser pread(el,donotresize);
      pread.put(el_ptr);
    }

    xmlXPathFreeObject(result);

    if(!ion.getTotalNum()) {
      ion.setName("i");
      ion.create(1);
      ion.R[0] = 0.0;
    }

    return true;
  }

  bool MO2Grid3D::setWavefunctions(xmlNodePtr aroot) {

    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"//wavefunction/determinantset",m_context);

    ///make a temporary array to pass over JastrowBuilder
    vector<ParticleSet*> PtclSets;
    PtclSets.push_back(&ion);
    PtclSets.push_back(&el);

    bool foundwfs=true;
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      ERRORMSG("wavefunction/determinantset is missing. Exit." << endl)
      foundwfs=false;
    } else {
      generateNumericalOrbitals(result->nodesetval->nodeTab[0]);
    }
    xmlXPathFreeObject(result);
    return foundwfs;
  }

  bool MO2Grid3D::generateNumericalOrbitals(xmlNodePtr cur) {

    LOGMSG("MO2Grid3D::generateNumericalOrbitals from Molecular Orbitals")

    string troot(myProject.m_title);
    const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"target");
    if(a) {
      troot = (const char*)a;
    } 

    GridMolecularOrbitals Original(Psi,ion,el);

    typedef OrbitalBase::RealType RealType;
    typedef OrbitalBase::ValueType ValueType;
    typedef OrbitalBase::PosType PosType;

    typedef GridMolecularOrbitals::BasisSetType BasisSetType;
    typedef LCOrbitals<BasisSetType>            SPOSetType;

    int nels=0;

    std::vector<RealType> ri(3,-5.0);
    std::vector<RealType> rf(3,5.0);
    std::vector<int> npts(3,101);

    vector<SPOSetType*> InOrbs;
    map<string,int> DetCounter;

    BasisSetType *basisSet=0;
    cur = cur->xmlChildrenNode;
    int idir=0;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "spline") {
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
	  string tname((const char*)(tcur->name));
          if(tname == "grid") {
            a=xmlGetProp(tcur,(const xmlChar*)"dir");
            if(a) { idir=atoi((const char*)a);}
            a=xmlGetProp(tcur,(const xmlChar*)"ri");
            if(a) { ri[idir]=atof((const char*)a);}
            a=xmlGetProp(tcur,(const xmlChar*)"rf");
            if(a) { rf[idir]=atof((const char*)a);}
            a=xmlGetProp(tcur,(const xmlChar*)"npts");
            if(a) { npts[idir]=atoi((const char*)a);}
          }
          tcur = tcur->next;
        }
      }
      else if(cname == OrbitalBuilderBase::basisset_tag) {
        basisSet = Original.addBasisSet(cur);
      } else if(cname == OrbitalBuilderBase::sd_tag) {
        nels=0;
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
	  string tname((const char*)(tcur->name));
	  if(tname == OrbitalBuilderBase::det_tag) {
            string detname("det");
            a=xmlGetProp(tcur,(const xmlChar*)"id");
            if(a) {
              detname = (const char*)a;
            } 
            //else {
            //  ostringstream idassigned(detname);
            //  idassigned << "dummy";
            //}
            bool newset=true;
            a=xmlGetProp(tcur,(const xmlChar*)"ref");
            if(a) {
              string detref((const char*)a);
              if(DetCounter.find(detref) != DetCounter.end()) {
                newset=false;
              }
            }
            a=xmlGetProp(tcur,(const xmlChar*)"orbitals");
            int norbs = atoi((const char*)a);
            if(newset) {
              SPOSetType* psi=new SPOSetType(basisSet,nels);
              psi->put(tcur);
              psi->setName(detname);
              psi->resize(norbs);

              DetCounter[detname]=InOrbs.size();
              InOrbs.push_back(psi);

              //InOrbs[detname]=psi;
            }
            nels+=norbs;
	  }
	  tcur = tcur->next;
	}
      }
      cur = cur->next;
    }

    //resize with respect to the number of electrons
    basisSet->resize(nels);

    DistanceTable::create(1);
    el.setUpdateMode(MCWalkerConfiguration::Update_Particle);

    //Need only one electron to calculate this
    //map<string,SPOSetType*>::iterator oit(InOrbs.begin());
    SPOSetType* inorb=InOrbs[0];
    //Create one-dimensional grids for three orthogonal directions
    typedef LinearGrid<double> GridType;
    GridType gridX, gridY, gridZ;
    gridX.set(ri[0],rf[0],npts[0]);
    gridY.set(ri[1],rf[1],npts[1]);
    gridZ.set(ri[2],rf[2],npts[2]);

    cout << "XYZCubicGrid " << endl;
    cout << " x " << ri[0] << " " << rf[0] << " " << npts[0] << endl;
    cout << " y " << ri[1] << " " << rf[1] << " " << npts[1] << endl;
    cout << " z " << ri[2] << " " << rf[2] << " " << npts[2] << endl;

    typedef TriCubicSplineT<ValueType> NOType;

    XYZCubicGrid<ValueType> grid3(&gridX,&gridY,&gridZ);
    NOType torb(&grid3);

    int ntot = npts[0]*npts[1]*npts[2];
    vector<ValueType> phi(inorb->numOrbitals(),0.0);
    vector<ValueType> dat(ntot);

    Pooma::Clock timer;

    std::vector<RealType> lapsed_time(5,0.0);
    PosType pos(el.R[0]);
    for(int iorb=0; iorb<inorb->numOrbitals(); iorb++) {

      //evaluate the values on the grid points
      timer.start();
      for(int ix=1; ix<npts[0]-1; ix++) {
        RealType x(gridX(ix));
        for(int iy=1; iy<npts[1]-1; iy++) {
          RealType y(gridY(iy));
          for(int iz=1; iz<npts[2]-1; iz++) {
            PosType dr(x,y,gridZ(iz));
            PosType newpos(el.makeMove(0,dr-pos));
            torb(ix,iy,iz) = inorb->evaluate(el,0,iorb);
          }
        }
      }
      timer.stop();
      lapsed_time[0]+= timer.cpu_time();

      //spline
      timer.start();
      torb.reset(false);
      timer.stop();
      lapsed_time[1]+= timer.cpu_time();

      //test
      timer.start();
      for(int ix=0; ix<npts[0]-1; ix++) {
        double x(gridX(ix));
        for(int iy=0; iy<npts[1]-1; iy++) {
          double y(gridY(iy));
          int offset=npts[2]*(iy+npts[1]*ix);
          for(int iz=0; iz<npts[2]-1; iz++,offset++) {
             TinyVector<double,3> p(x,y,gridZ(iz));
             dat[offset]=torb.evaluate(p);
          }
        }
      }
      timer.stop();
      lapsed_time[2]+= timer.cpu_time();

      timer.start();
      char oname[128];
      sprintf(oname,"%s.wf%04d.h5",troot.c_str(),iorb);
      hid_t h_file = H5Fcreate(oname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      HDFAttribIO<std::vector<double> > dump(dat,npts);
      dump.write(h_file,"Orbital");
      HDFAttribIO<TriCubicSplineT<double> > dump1(torb);
      dump1.write(h_file,"CubicSpline");
      H5Fclose(h_file);
      timer.stop();
      lapsed_time[3]+= timer.cpu_time();
    }

    cout << "Timing results in sec" << endl;
    cout << "Function evaluation " << inorb->numOrbitals() << " orbitals = " << lapsed_time[0] << endl;
    cout << "Spline coefficients = " << lapsed_time[1] << endl;
    cout << "Testing spline      = " << lapsed_time[2] << endl;
    cout << "Writing hdf5 files  = " << lapsed_time[3] << endl;
    //HDFAttribIO<TriCubicSplineT<double> > dump1(torb);
    //dump1.write(h_file,"spline0000");
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
