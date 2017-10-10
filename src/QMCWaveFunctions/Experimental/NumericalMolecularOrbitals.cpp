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
    
    
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/MolecularOrbitals/NumericalMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/LCOrbitals.h"
#include "OhmmsData/AttributeSet.h"
#if defined(HAVE_LIBHDF5)
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFTriCubicSpline.h"
#else
#include "Numerics/TriCubicSplineT.h"
#endif
#include "Utilities/Clock.h"

namespace qmcplusplus
{

NumericalMolecularOrbitals::NumericalMolecularOrbitals(ParticleSet& els,
    TrialWaveFunction& wfs, ParticleSet& ions):
  OrbitalBuilderBase(els,wfs), Original(0)
{
  Original = new GridMolecularOrbitals(els,wfs,ions);
}

bool NumericalMolecularOrbitals::put(xmlNodePtr cur)
{
  LOGMSG("NumericMolecularOrbitals::put")
  typedef GridMolecularOrbitals::BasisSetType BasisSetType;
  typedef LCOrbitals<BasisSetType>            SPOSetType;
  int nels=0;
  //map<std::string,SPOSetType*> InOrbs;
  std::vector<SPOSetType*> InOrbs;
  std::map<std::string,int> DetCounter;
  BasisSetType *basisSet=0;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == basisset_tag)
    {
      basisSet = Original->addBasisSet(cur);
    }
    else
      if(cname == sd_tag)
      {
        nels=0;
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL)
        {
          std::string tname((const char*)(tcur->name));
          if(tname == det_tag)
          {
            std::string detname("det");
            const xmlChar* a=xmlGetProp(tcur,(const xmlChar*)"id");
            if(a)
            {
              detname = (const char*)a;
            }
            //else {
            //  std::ostringstream idassigned(detname);
            //  idassigned << "dummy";
            //}
            bool newset=true;
            a=xmlGetProp(tcur,(const xmlChar*)"ref");
            if(a)
            {
              std::string detref((const char*)a);
              if(DetCounter.find(detref) != DetCounter.end())
              {
                newset=false;
              }
            }
            a=xmlGetProp(tcur,(const xmlChar*)"orbitals");
            int norbs = atoi((const char*)a);
            if(newset)
            {
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
  targetPtcl.update(1);
  //Need only one electron to calculate this
  //map<std::string,SPOSetType*>::iterator oit(InOrbs.begin());
  SPOSetType* inorb=InOrbs[0];
  double ri = -5.0;
  double rf = 5.0;
  std::vector<int> npts(3);
  npts[0]=101;
  npts[1]=101;
  npts[2]=101;
  //Create one-dimensional grids for three orthogonal directions
  typedef LinearGrid<double> GridType;
  GridType gridX, gridY, gridZ;
  gridX.set(ri,rf,npts[0]);
  gridY.set(ri,rf,npts[1]);
  gridZ.set(ri,rf,npts[2]);
  XYZCubicGrid<ValueType> grid3(&gridX,&gridY,&gridZ);
  typedef TriCubicSplineT<ValueType> NOType;
  std::vector<NOType*> torb;
  for(int i=0; i<inorb->numOrbitals(); i++)
  {
    torb.push_back(new NOType(&grid3));
    if(i)
      torb[i]->setGridManager(false);
  }
  //TriCubicSplineT<ValueType> torb(&grid3);
  std::vector<ValueType> phi(inorb->numOrbitals(),0.0);
  std::vector<std::vector<ValueType>* > dat;
  int ntot = npts[0]*npts[1]*npts[2];
  for(int is=0; is<inorb->numOrbitals(); is++)
    dat.push_back(new std::vector<ValueType>(ntot,0));
  //int lastpsi=inorb->numOrbitals()-2;
  //int lastpsi=0;
  PosType pos(targetPtcl.R[0]);
  int ng=0;
  Pooma::Clock timer;
  timer.start();
  for(int ix=0; ix<npts[0]; ix++)
  {
    RealType x(gridX(ix));
    for(int iy=0; iy<npts[1]; iy++)
    {
      RealType y(gridY(iy));
      for(int iz=0; iz<npts[2]; iz++)
      {
        PosType dr(x,y,gridZ(iz));
        PosType newpos(targetPtcl.makeMove(0,dr-pos));
        int item=torb[0]->index(ix,iy,iz);
        inorb->evaluate(targetPtcl,0,phi);
        for(int is=0; is<phi.size(); is++)
          (*torb[is])(item)=phi[is];
        //for(int is=0; is<phi.size(); is++) dat[is]->operator[](offset)=phi[is];
        //torb(ix,iy,iz)=phi[lastpsi];
        //dat[offset]=phi[lastpsi];
      }
    }
  }
  timer.stop();
  timer.start();
  for(int is=0; is<torb.size(); is++)
    torb[is]->reset(true);
  timer.stop();
  for(int ix=0; ix<npts[0]-1; ix++)
  {
    double x(gridX(ix));
    for(int iy=0; iy<npts[1]-1; iy++)
    {
      double y(gridY(iy));
      int offset=npts[2]*(iy+npts[1]*ix);
      for(int iz=0; iz<npts[2]-1; iz++,offset++)
      {
        TinyVector<double,3> p(x,y,gridZ(iz));
        for(int is=0; is<dat.size(); is++)
          dat[is]->operator[](offset) = torb[is]->evaluate(p);
      }
    }
  }
#if defined(HAVE_LIBHDF5)
  for(int is=0; is<dat.size(); is++)
  {
    char oname[128];
    sprintf(oname,"out%04d.h5",is);
    hid_t h_file = H5Fcreate(oname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    HDFAttribIO<std::vector<double> > dump(*dat[is],npts);
    dump.write(h_file,"Eigvector");
    HDFAttribIO<TriCubicSplineT<double> > dump1(*torb[is]);
    dump1.write(h_file,"CubicSpline");
    H5Fclose(h_file);
  }
  //HDFAttribIO<TriCubicSplineT<double> > dump1(torb);
  //dump1.write(h_file,"spline0000");
  //stop for now
  //OHMMS::Controller->finalize();
#endif
  return true;
}
}
