//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_RADIAL_NUMERICALGRIDORBITALBUILDER_H
#define QMCPLUSPLUS_RADIAL_NUMERICALGRIDORBITALBUILDER_H

#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"
#include <HDFVersion.h>
#include "Utilities/OhmmsInfo.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/GaussianBasisSet.h"
#include "Numerics/SlaterBasisSet.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCFactory/OneDimGridFactory.h"

namespace qmcplusplus
{

/** Build a set of radial orbitals at the origin
 *
 * For a center,
 *   - only one grid is used
 *   - any number of radial orbitals
 */
template<typename COT>
class RadialOrbitalSetBuilder: public QMCTraits
{

public:
  typedef typename COT::RadialOrbital_t RadialOrbitalType;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///the radial orbitals
  COT* m_orbitals;
  ///input grid in case transform is needed
  GridType *input_grid;
  ///maximum cutoff
  RealType m_rcut;
  ///the quantum number of this node
  QuantumNumberType m_nlms;
  ///the species
  std::string m_species;
  ///type of input function
  std::string m_infunctype;

  ///constructor
  RadialOrbitalSetBuilder(xmlNodePtr cur=NULL);
  ///destructor
  ~RadialOrbitalSetBuilder();

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(COT* oset, const std::string& acenter);

  ///add a grid
  bool addGrid(xmlNodePtr cur);

  /** add a radial functor
   * @param cur xml element
   * @param nlms quantum number
   */
  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

  /** put common element
   * @param cur xml element
   */
  bool putCommon(xmlNodePtr cur);

private:
  void addGaussian(xmlNodePtr cur);
  void addSlater(xmlNodePtr cur);
  void addNumerical(xmlNodePtr cur, const std::string& dsname);
  hid_t m_fileid;
};

  template<typename COT>
    RadialOrbitalSetBuilder<COT>::RadialOrbitalSetBuilder(xmlNodePtr cur)
    : Normalized(true),m_orbitals(0),input_grid(nullptr),m_rcut(-1.0), 
    m_infunctype("Gaussian"), m_fileid(-1)
  {
    if(cur != NULL)
      putCommon(cur);
  }

  template<typename COT>
    RadialOrbitalSetBuilder<COT>::~RadialOrbitalSetBuilder()
    {
      if(m_fileid>-1)
        H5Fclose(m_fileid);
      if(input_grid != nullptr)
        delete input_grid;
    }

  template<typename COT>
    bool RadialOrbitalSetBuilder<COT>::putCommon(xmlNodePtr cur)
    {
      std::string normin("yes");
      std::string afilename("0");
      OhmmsAttributeSet aAttrib;
      aAttrib.add(normin,"normalized");
      aAttrib.add(m_infunctype,"type");
      aAttrib.add(afilename,"href");
      aAttrib.add(afilename,"src");
      bool success=aAttrib.put(cur);
      //set the noramlization
      Normalized=(normin=="yes");
      if(afilename.find("h5")<afilename.size())
      {
        m_fileid = H5Fopen(afilename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        //current version
        HDFVersion res_version(0,1);
        HDFVersion in_version(0,0); //start using major=0 and minor=4
        herr_t status = H5Eset_auto(NULL, NULL);
        in_version.read(m_fileid,hdf::version);
        if(in_version<res_version)
        {
          APP_ABORT("NGOBuilder::putCommon Old format is not supported. Please rerun SQD if the output is produced by it");
        }
        app_log() << "  " << afilename <<  " version " << in_version << std::endl;
      }
      else
        m_fileid=-1;
      return success;
    }

  template<typename COT>
    void
    RadialOrbitalSetBuilder<COT>::setOrbitalSet(COT* oset, const std::string& acenter)
    {
      m_orbitals = oset;
      m_species = acenter;
    }

  template<typename COT>
    bool
    RadialOrbitalSetBuilder<COT>::addGrid(xmlNodePtr cur)
    {
      if(!m_orbitals)
      {
        APP_ABORT("NGOBuilder::addGrid SphericalOrbitals<ROT,GT>*, is not initialized");
      }
      if(m_fileid<0)
      {
        GridType *agrid = OneDimGridFactory::createGrid(cur);
        m_orbitals->Grids.push_back(agrid);
      }
      else
      {
        app_log() << "   Grid is created by the input paremters in h5" << std::endl;
        hid_t gid = H5Gopen(m_fileid,"radial_basis_states/grid");
        std::string gridtype;
        HDFAttribIO<std::string> gtypestr(gridtype);
        gtypestr.read(gid,"type");
        int npts=0;
        RealType ri=0.0,rf=10.0,rmax_safe=10;
        double tt=0;
        HDFAttribIO<double> ri_in(tt);
        ri_in.read(gid,"ri");
        ri=tt;
        ri_in.read(gid,"rf");
        rf=tt;
        ri_in.read(gid,"rmax_safe");
        rmax_safe=tt;
        HDFAttribIO<int> n_in(npts);
        n_in.read(gid,"npts");
        if(gridtype.empty())
        {
          APP_ABORT("Grid type is not specified.");
        }
        if(gridtype == "log")
        {
          app_log() << "    Using log grid ri = " << ri << " rf = " << rf << " npts = " << npts << std::endl;
          input_grid = new LogGrid<RealType>;
          input_grid->set(ri,rf,npts);
          //GridType *agrid = new LinearGrid<RealType>;
          //agrid->set(0.0,rmax_safe,rmax_safe/0.01);
          //m_orbitals->Grids.push_back(agrid);
          m_orbitals->Grids.push_back(input_grid);
          input_grid=0;
        }
        else
          if(gridtype == "linear")
          {
            app_log() << "    Using linear grid ri = " << ri << " rf = " << rf << " npts = " << npts << std::endl;
            input_grid = new LinearGrid<RealType>;
            input_grid->set(ri,rf,npts);
            m_orbitals->Grids.push_back(input_grid);
            input_grid=0;
          }
        //if(!input_grid)
        //{
        //  APP_ABORT("Grid is not defined.");
        //}
        //using 0.01 for the time being
      }
      return true;
    }

  /** Add a new Slater Type Orbital with quantum numbers \f$(n,l,m,s)\f$
   * \param cur  the current xmlNode to be processed
   * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
   * \return true is succeeds
   *
   This function puts the STO on a logarithmic grid and calculates the boundary
   conditions for the 1D Cubic Spline.  The derivates at the endpoint
   are assumed to be all zero.  Note: for the radial orbital we use
   \f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
   radial orbital and \f$ l \f$ is the angular momentum.
   */
  template<typename COT>
    bool
    RadialOrbitalSetBuilder<COT>::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms)
    {
      if(!m_orbitals)
      {
        ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
          return false;
      }
      std::string radtype(m_infunctype);
      std::string dsname("0");
      OhmmsAttributeSet aAttrib;
      aAttrib.add(radtype,"type");
      aAttrib.add(m_rcut,"rmax");
      aAttrib.add(dsname,"ds");
      aAttrib.put(cur);
      //const xmlChar *tptr = xmlGetProp(cur,(const xmlChar*)"type");
      //if(tptr) radtype = (const char*)tptr;
      //tptr = xmlGetProp(cur,(const xmlChar*)"rmax");
      //if(tptr) m_rcut = atof((const char*)tptr);
      int lastRnl = m_orbitals->Rnl.size();
      m_nlms = nlms;
      if(radtype == "Gaussian" || radtype == "GTO")
      {
        addGaussian(cur);
      }
      else if(radtype == "Slater" || radtype == "STO")
      {
        addSlater(cur);
      }
      else
      {
        addNumerical(cur,dsname);
      }
      if(lastRnl && m_orbitals->Rnl.size()> lastRnl)
      {
        app_log() << "\tSetting GridManager of " << lastRnl << " radial orbital to false" << std::endl;
        m_orbitals->Rnl[lastRnl]->setGridManager(false);
      }
      return true;
    }

  template<typename COT>
  void RadialOrbitalSetBuilder::addGaussian(xmlNodePtr cur)
  {
    int L= m_nlms[1];
    GaussianCombo<RealType> gset(L,Normalized);
    gset.putBasisGroup(cur);
    GridType* agrid = m_orbitals->Grids[0];
    RadialOrbitalType *radorb = new RadialOrbitalType(agrid);
    if(m_rcut<0)
      m_rcut = agrid->rmax();
    Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gset, *radorb);
    transform.generate(agrid->rmin(),m_rcut,agrid->size());
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(m_nlms);
  }

  template<typename COT>
  void RadialOrbitalSetBuilder::addSlater(xmlNodePtr cur)
  {
    ////pointer to the grid
    GridType* agrid = m_orbitals->Grids[0];
    RadialOrbitalType *radorb = new RadialOrbitalType(agrid);
    SlaterCombo<RealType> sto(m_nlms[1],Normalized);
    sto.putBasisGroup(cur);
    //spline the slater type orbital
    Transform2GridFunctor<SlaterCombo<RealType>,RadialOrbitalType> transform(sto, *radorb);
    transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());
    //transform.generate(agrid->rmax());
    //add the radial orbital to the list
    m_orbitals->Rnl.push_back(radorb);
    m_orbitals->RnlID.push_back(m_nlms);
  }

  template<typename COT>
  void RadialOrbitalSetBuilder::addNumerical(xmlNodePtr cur, const std::string& dsname)
  {
    int imin = 0;
    OhmmsAttributeSet aAttrib;
    aAttrib.add(imin,"imin");
    aAttrib.put(cur);
    char grpname[128];
    sprintf(grpname,"radial_basis_states/%s",dsname.c_str());
    hid_t group_id_orb=H5Gopen(m_fileid,grpname);
    double cusp_cond=0.0;
    HDFAttribIO<double> r_in(cusp_cond);
    r_in.read(group_id_orb,"cusp");
    //int rinv_p=0;
    //HDFAttribIO<int> ir(rinv_p);
    //ir.read(group_id_orb,"power");
    Vector<double> rad_orb;
    HDFAttribIO<Vector<double> > rin(rad_orb);
    rin.read(group_id_orb,"radial_orbital");
    H5Gclose(group_id_orb);
    //if(input_grid)
    //{
    //  int imax = rad_orb.size()-1;
    //  RadialOrbitalType torb(input_grid,rad_orb);
    //  RealType yprime_i = (rad_orb[imin+1]-rad_orb[imin])/((input_grid->r(imin+1)-input_grid->r(imin)));
    //  torb.spline(imin,yprime_i,imax,0.0);
    //  GridType* agrid = m_orbitals->Grids[0];
    //  std::vector<double> orb_linear(agrid->size(),0.0);
    //  for(int ig=0; ig<agrid->size()-2; ++ig) orb_linear[ig]=torb.f(agrid->r(ig));
    //  RadialOrbitalType *radorb = new RadialOrbitalType(agrid,orb_linear);
    //  radorb->spline(imin,yprime_i,agrid->size()-2,0.0);
    //  m_orbitals->Rnl.push_back(radorb);
    //}
    //else
    {
      GridType* agrid = m_orbitals->Grids[0];
      int imax = rad_orb.size()-1;
      RadialOrbitalType *radorb = new RadialOrbitalType(agrid,rad_orb);
      //calculate boundary condition, assume derivates at endpoint are 0.0
      RealType yprime_i = (rad_orb[imin+1]-rad_orb[imin])/((agrid->r(imin+1)-agrid->r(imin)));
      radorb->spline(imin,yprime_i,imax,0.0);
      m_orbitals->Rnl.push_back(radorb);
    }
    m_orbitals->RnlID.push_back(m_nlms);
    //ofstream dfile("spline.dat");
    //dfile.setf(std::ios::scientific, std::ios::floatfield);
    //for(int ig=imin; ig<radorb->size(); ig++) {
    //  RealType dr = (radorb->r(ig+1)- radorb->r(ig))/5.0;
    //  RealType _r = radorb->r(ig),y,dy,d2y;
    //  while(_r<radorb->r(ig+1)) {
    //    //Do not need this anymore
    //    //radorb->setgrid(_r);
    //    y = radorb->evaluate(_r,1.0/_r,dy,d2y);
    //    dfile << std::setw(15) << _r << std::setw(20) << setprecision(12) << y
    //          << std::setw(20) << dy << std::setw(20) << d2y
    //          << std::endl;
    //    _r+=dr;
    //  }
    //}
    //cout << " Power " << rinv_p << " size=" << rad_orb.size() << std::endl;
    //APP_ABORT("NGOBuilder::addNumerical");
  }
}
#endif
