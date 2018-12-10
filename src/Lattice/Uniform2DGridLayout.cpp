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
    
    



#include "Lattice/Uniform2DGridLayout.h"
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Utilities/OutputManager.h"
#include <limits>

namespace qmcplusplus
{

///Set LR_rc = radius of smallest sphere inside box and kc=dim/rc
void Uniform2DGridLayout::SetLRCutoffs()
{
  //Compute rc as the real-space cutoff of 1/2 the unit-cell.
  //Radius of maximum shere that fits in a...
  TinyVector<value_type,2> b,c,d,x; //Unit vector of each surface will be in here
  //Compute the coordinate of the box center
  c = 0.5*(a(0)+a(1));
  Scalar_t dotP= dot(a(0), a(1));
  SingleParticlePos_t L0 = a(0) - dotP*a(1);
  SingleParticlePos_t L1 = a(1) - dotP*a(0);
  LR_rc = 0.5*std::min(std::sqrt(dot(L0,L0)), std::sqrt(dot(L1,L1)));
  LR_kc = LR_dim_cutoff/LR_rc;
  LOGMSG("\tLong-range breakup parameters:");
  LOGMSG("\trc*kc = " << LR_dim_cutoff << "; rc = " << LR_rc << "; kc = " << LR_kc << "\n");
}

void Uniform2DGridLayout::makeShell(std::vector<SingleParticleIndex_t>& RS,
                                    std::vector<bool>& FS,
                                    const SingleParticleIndex_t& nc,
                                    value_type rr)
{
  std::map<int,std::vector<SingleParticleIndex_t>*> rs;
  std::map<int,bool> cb;
  Scalar_t scaleL = static_cast<value_type>(8*NP[0]*NP[1]);
  int ic=0,ih;
  Scalar_t dx=Delta[0]*0.5;
  Scalar_t dy=Delta[1]*0.5;
  for(int ix1=-nc[0]; ix1<=nc[0]; ix1++)
  {
    Scalar_t x= std::abs(ix1)*Delta[0];
    for(int ix2=-nc[1]; ix2<=nc[1]; ix2++)
    {
      if(ix1==0 && ix2 == 0)
        continue;
      Scalar_t y=std::abs(ix2)*Delta[1];
      ic++;
      SingleParticlePos_t tc(x-dx,y-dy);
      int ih = static_cast<int>(Dot(tc,tc)*scaleL);
      std::map<int,std::vector<SingleParticleIndex_t>*>::iterator it = rs.find(ih);
      if(it == rs.end())
      {
        std::vector<SingleParticleIndex_t>* ns = new std::vector<SingleParticleIndex_t>;
        ns->push_back(SingleParticleIndex_t(ix1,ix2));
        rs[ih] = ns;
      }
      else
      {
        (*it).second->push_back(SingleParticleIndex_t(ix1,ix2));
      }
    }
  }
  RS.resize(ic);
  FS.resize(ic);
  int irc = static_cast<int>(4.0*rr*scaleL);
  int ir=0;
  std::map<int,std::vector<SingleParticleIndex_t>*>::const_iterator cit = rs.begin();
  while(cit != rs.end())
  {
    bool inside=((*cit).first <= irc);
    const std::vector<SingleParticleIndex_t>& iv = *((*cit).second);
    for(int i=0; i<iv.size(); i++, ir++)
    {
      RS[ir] = iv[i];
      FS[ir] = inside;
    }
    cit++;
  }
  std::map<int,std::vector<SingleParticleIndex_t>*>::iterator it = rs.begin();
  while(it != rs.end())
  {
    delete (*it).second;
    it++;
  }
}

void Uniform2DGridLayout::checkGrid(value_type int_rad)
{
  if(int_rad<std::numeric_limits<Scalar_t>::epsilon())
  {
    WARNMSG("Invalid interaction radius.")
    return;
  }
  value_type RmaxSq = int_rad*int_rad;
  //calculate the grid based on the interaction radius
  int ncx=static_cast<int>(sqrt(dot(Rv[0],Rv[0])/RmaxSq));
  int ncy=static_cast<int>(sqrt(dot(Rv[1],Rv[1])/RmaxSq));
  ncx = (ncx)? ncx:ncx+1;
  ncy = (ncy)? ncy:ncy+1;
  if(Grid[SPATIAL_GRID][0]*Grid[SPATIAL_GRID][1]== 1)
  {
    Grid[SPATIAL_GRID][0] = ncx;
    Grid[SPATIAL_GRID][1] = ncy;
  }
  bool initOmpGrid=false;
  int np = Grid[OMP_GRID][0]*Grid[OMP_GRID][1];
  if(np != omp_get_max_threads())
  {
    WARNMSG("\tInput OpenMP grid " << Grid[OMP_GRID] << " does not fit to " << np << " threads")
    int idir = (ncx>ncy)? 0:1;
    int jdir = (ncx>ncy)? 1:0;
    Grid[OMP_GRID][idir] =1;
    Grid[OMP_GRID][jdir] = omp_get_max_threads();
    np = omp_get_max_threads();
  }
  if(np>1)
  {
    initOmpGrid=true;
    WARNMSG("\tInput OpenMP grid = " << Grid[OMP_GRID])
    Grid[SPATIAL_GRID][0]=(Grid[SPATIAL_GRID][0]/Grid[OMP_GRID][0])*Grid[OMP_GRID][0];
    Grid[SPATIAL_GRID][1]=(Grid[SPATIAL_GRID][1]/Grid[OMP_GRID][1])*Grid[OMP_GRID][1];
    WARNMSG("\tDetecting omp grid setting. Creating s supergrid")
    Grid[SPATIAL_GRID][0] =(Grid[SPATIAL_GRID][0])?Grid[SPATIAL_GRID][0]:1;
    Grid[SPATIAL_GRID][1] =(Grid[SPATIAL_GRID][1])?Grid[SPATIAL_GRID][1]:1;
  }
  setGrid(Grid[SPATIAL_GRID]);
  if(initOmpGrid && SuperGrid[OMP_GRID]==0)
  {
    Grid_t* agrid = addGrid(OMP_GRID,Grid[OMP_GRID]);
    initGrid(*agrid,*this);
  }
}

/** initialize a grid
 * @param agrid grid whose data will be created
 * @param subgrid subgrid to be partitioned to agrid
 *
 * Default method to create a supergrid. Subgrid is divided among the grid points of the
 * supergrid based on spatial locality.
 */
void Uniform2DGridLayout::initGrid(Grid_t& agrid, const Grid_t& subgrid)
{
  Scalar_t del0=subgrid.Delta[0];
  Scalar_t del1=subgrid.Delta[1];
  Scalar_t eps0=0.5*subgrid.Delta[0];
  Scalar_t eps1=0.5*subgrid.Delta[1];
  agrid.createData();
  for(int ig=0; ig<subgrid.NP[0]; ig++)
  {
    for(int jg=0; jg<subgrid.NP[1]; jg++)
    {
      int iloc = agrid.loc(ig*del0+eps0,jg*del1+eps1);
      agrid.addData(iloc,subgrid.loc(ig,jg));
    }
  }
}

int Uniform2DGridLayout::connectGrid(value_type int_rad, value_type con_rad)
{
  if(int_rad<std::numeric_limits<Scalar_t>::epsilon())
    return 1;
  //{  ///create the spatial grid
  //  setGrid(Grid[SPATIAL_GRID]);
  //  SingleParticlePos_t v0(Delta[0],0.0,0.0);
  //  SingleParticlePos_t v1(0.0,Delta[1],0.0);
  //  SingleParticlePos_t v2(0.0,0.0,Delta[2]);
  //  int_rad = Dot(v0,v0);
  //  int_rad = std::min(int_rad, Dot(v1,v1));
  //  int_rad = std::min(int_rad, Dot(v2,v2));
  //  int_rad = sqrt(int_rad);
  //  WARNMSG("The negative interaction radius is overwritten by the grid setting")
  //  WARNMSG("The new interaction radius = " << int_rad)
  //}
  value_type RmaxSq = int_rad*int_rad;
  if(con_rad<0)
    con_rad = int_rad;
  InteractionRadius = int_rad;
  ConnectionRadius = con_rad;
  SingleParticlePos_t u0(Delta[0],0.0);
  SingleParticlePos_t u1(0.0,Delta[1]);
  Scalar_t rc2 = con_rad*con_rad;
  int nx = static_cast<int>(sqrt(rc2/Dot(u0,u0)))+1;
  int ny = static_cast<int>(sqrt(rc2/Dot(u1,u1)))+1;
  if(nx == NCMax[0] && ny == NCMax[1])
  {
    for(int i=0; i<u_bc.size(); i++)
      c_bc[i]=toCart(u_bc[i]);
    return MaxConnections;
  }
  NCMax[0] = nx;
  NCMax[1] = ny;
  std::vector<SingleParticleIndex_t> RShell;
  std::vector<bool> FirstShell;
  makeShell(RShell,FirstShell,NCMax,RmaxSq);
  c_offset.resize(NumGrids+1);
  c_max.resize(NumGrids);
  int ntot = NumGrids*(2*nx+1)*(2*ny+1);
  if(c_id.capacity() < ntot)
    c_id.reserve(ntot);
  if(c_bc.capacity() < ntot)
    c_bc.reserve(ntot);
  if(u_bc.capacity() < ntot)
    u_bc.reserve(ntot);
  int maxnc = 0, gtot = 0;
  SingleParticlePos_t dx(Delta[0],Delta[1]),org,d,dp;
  c_offset[0] = 0;
  for(int ig=0; ig<NP[0]; ig++)
  {
    for(int jg=0; jg<NP[1]; jg++)
    {
      int nconnected = 0, ncinside=0;
      SingleParticleIndex_t org(ig,jg);
      for(int ixyz=0; ixyz<RShell.size(); ixyz++)
      {
        SingleParticleIndex_t d0 = org+RShell[ixyz],d;
        if(!BoxBConds[0] && d[0]<0 || d[0]>= NP[0])
          continue;
        if(!BoxBConds[1] && d[1]<0 || d[1]>= NP[1])
          continue;
        d[0]=d0[0]%NP[0];
        d[0]= (d[0]<0)?d[0]+NP[0] : d[0];
        d[1]=d0[1]%NP[1];
        d[1]= (d[1]<0)?d[1]+NP[1] : d[1];
        int iloc = d[1]+NP[1]*d[0];
        d0 -= d;
        d0[0]= (d0[0]-d[0])/NP[0];
        d0[1]= (d0[1]-d[1])/NP[1];
        c_id.push_back(iloc);
        u_bc.push_back(d0);
        c_bc.push_back(toCart(d0));
        if(FirstShell[ixyz])
          ncinside++;
        nconnected++;
      }
      c_offset[gtot+1] = c_offset[gtot]+nconnected;
      c_max[gtot] = c_offset[gtot]+ncinside;
      gtot++;
      //maxnc = std::max(maxnc,nconnected);
      maxnc = std::max(maxnc,ncinside);
    }
  }
  MaxConnections = maxnc;
  //print(std::cout);
  return maxnc; // return the maxmimum number of connected cells
}

void Uniform2DGridLayout::print(std::ostream& os) const
{
  os << "<unitcell>" << std::endl;
  Base_t::print(os);
  os << "<note>" << std::endl;
  os << "\tLong-range breakup parameters:" << std::endl;
  os << "\trc*kc = " << LR_dim_cutoff << "; rc = " << LR_rc << "; kc = " << LR_kc << "\n" << std::endl;
  os << "</note>" << std::endl;
  os << "</unitcell>" << std::endl;
  ////printGrid(os);
  //for(int ig=0; ig<c_offset.size()-1; ig++) {
  //  os << ig << " has neighboring cell "
  //     << c_max[ig]-c_offset[ig] << " "
  //     << c_offset[ig+1]-c_offset[ig]<< std::endl;
  //  //for(int ii=c_offset[ig]; ii<c_max[ig]; ii++) {
  //  //  os << c_id[ii] << " " << c_bc[ii] << std::endl;
  //  //}
  //}
}
}

