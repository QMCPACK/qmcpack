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


#include <HDFVersion.h>
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/GaussianBasisSet.h"
#include "Numerics/SlaterBasisSet.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCFactory/OneDimGridFactory.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "io/hdf_archive.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "Message/CommOperators.h"


namespace qmcplusplus
{

OptimizableFunctorBase* NGOrbital::makeClone() const
{
  NGOrbital *myclone=new NGOrbital(*this);
  myclone->myFunc.m_grid=myFunc.m_grid->makeClone();
  myclone->setGridManager(true);
  return myclone;
}

NGOBuilder::NGOBuilder(xmlNodePtr cur)
  : Normalized(true),m_orbitals(0),input_grid(0)
  ,m_rcut(-1.0), m_infunctype("Gaussian")
  , m_fileid(-1)
{
  if(cur != NULL)
    putCommon(cur);
}

NGOBuilder::~NGOBuilder()
{
  if(m_fileid>-1)
    H5Fclose(m_fileid);
  if(input_grid)
    delete input_grid;
}

bool NGOBuilder::putCommon(xmlNodePtr cur)
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

void
NGOBuilder::setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter)
{
  m_orbitals = oset;
  m_species = acenter;
}

bool NGOBuilder::addGridH5(hdf_archive &hin)
{
  if(!m_orbitals)
  {
    APP_ABORT("NGOBuilder::addGrid SphericalOrbitals<ROT,GT>*, is not initialized");
  }

  app_log() << "   Grid is created by the input paremters in h5" << std::endl;



  std::string gridtype;

  if(hin.myComm->rank()==0){
     if(!hin.read(gridtype, "grid_type")){
         std::cerr<<"Could not read grid_type in H5; Probably Corrupt H5 file"<<std::endl;
         exit(0);
     }
  }
  hin.myComm->bcast(gridtype);

  int npts=0;
  RealType ri=0.0,rf=10.0,rmax_safe=10;

  if(hin.myComm->rank()==0){
      double tt=0;
      hin.read(tt,"grid_ri");
      ri=tt;
      hin.read(tt,"grid_rf");
      rf=tt;
      hin.read(tt,"rmax_safe");
      rmax_safe=tt;
      hin.read(npts,"grid_npts");
  }
  hin.myComm->bcast(ri);
  hin.myComm->bcast(rf);
  hin.myComm->bcast(rmax_safe);
  hin.myComm->bcast(npts);


  if(gridtype.empty())
  {
    APP_ABORT("Grid type is not specified.");
  }
  if(gridtype == "log")
  {
    app_log() << "    Using log grid ri = " << ri << " rf = " << rf << " npts = " << npts << std::endl;
    input_grid = new LogGrid<RealType>;
    input_grid->set(ri,rf,npts);
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
  return true;
}

bool NGOBuilder::addGrid(xmlNodePtr cur)
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
bool
NGOBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms)
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
  else
    if(radtype == "Slater" || radtype == "STO")
    {
      addSlater(cur);
    }
    else
      if(radtype == "Pade")
      {
        app_error() << "  Any2GridBuilder::addPade is disabled." << std::endl;
        APP_ABORT("NGOBuilder::addRadialOrbital");
        //addPade(cur);
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

bool
NGOBuilder::addRadialOrbitalH5(hdf_archive & hin, const QuantumNumberType& nlms)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  std::string radtype(m_infunctype);
  std::string dsname("0");
  int lastRnl = m_orbitals->Rnl.size();
  m_nlms = nlms;
  if(radtype == "Gaussian" || radtype == "GTO")
  {
    addGaussianH5(hin);
  }
  else
  {
    app_log()<<" RadType other than Gaussian not implemented to be read from H5 format"<<std::endl;
    exit(0);
  }
  if(lastRnl && m_orbitals->Rnl.size()> lastRnl)
  {
    app_log() << "\tSetting GridManager of " << lastRnl << " radial orbital to false" << std::endl;
    m_orbitals->Rnl[lastRnl]->setGridManager(false);
  }

  return true;
}

void NGOBuilder::addGaussian(xmlNodePtr cur)
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

void NGOBuilder::addGaussianH5(hdf_archive &hin)
{
  int L= m_nlms[1];

  GaussianCombo<RealType> gset(L,Normalized);
  gset.putBasisGroupH5(hin);
  GridType* agrid = m_orbitals->Grids[0];
  RadialOrbitalType *radorb = new RadialOrbitalType(agrid);
  if(m_rcut<0)
    m_rcut = agrid->rmax();
  Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gset, *radorb);
  transform.generate(agrid->rmin(),m_rcut,agrid->size());
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(m_nlms);
}


void NGOBuilder::addSlater(xmlNodePtr cur)
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

void NGOBuilder::addNumerical(xmlNodePtr cur, const std::string& dsname)
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

template<class T>
struct PadeOrbital: public OptimizableFunctorBase
{

  typedef T value_type;
  real_type a0,a1,a2,a3,rcut;
  std::string  nodeName;

  explicit
  PadeOrbital(const char* node_name="radfunc"):
    a0(1.0),a1(-0.5),a2(0.0),a3(-0.2),rcut(4.0),nodeName(node_name) {}

  ~PadeOrbital() { }

  OptimizableFunctorBase* makeClone() const
  {
    return new PadeOrbital<T>(*this);
  }

  void reset() {}

  inline real_type f(real_type r)
  {
    return a0*std::exp((a1+a2*r)*r/(1.0e0+a3*r));
  }

  inline real_type df(real_type r)
  {
    value_type t = 1.0/(1.0e0+a3*r);
    value_type z=(a1+a2*r)*r*t;
    value_type res = a0*std::exp(z);
    return res*(a1+2.0*a2*r-z*a3)*t;
  }

  bool putBasisGroup(xmlNodePtr cur)
  {
    cur = cur->children;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == nodeName)
      {
        OhmmsAttributeSet radAttrib;
        radAttrib.add(a0,"a0");
        radAttrib.add(a1,"a1");
        radAttrib.add(a2,"a2");
        radAttrib.add(a3,"a3");
        radAttrib.put(cur);
        rcut = -1.0/a3-std::numeric_limits<T>::epsilon();
        LOGMSG("Pade compoenent (a0,a1,a2,a3,rcut) = " << a0 << " " << a1 << " " << a2 << " " << a3 << " " << rcut)
      }
      cur=cur->next;
    }
    return true;
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  void addOptimizables( VarRegistry<real_type>& vlist) {}
};

void NGOBuilder::addPade(xmlNodePtr cur)
{
  //GridType* agrid = m_orbitals->Grids[0];
  //RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);
  //PadeOrbital<RealType> pade;
  //pade.putBasisGroup(cur);
  ////spline the slater type orbital
  //Transform2GridFunctor<PadeOrbital<RealType>,RadialOrbitalType> transform(pade, *radorb);
  //if(pade.rcut>0)
  //  transform.generate(pade.rcut);
  //else
  //  transform.generate(agrid->rmax());
  ////add the radial orbital to the list
  //m_orbitals->Rnl.push_back(radorb);
  //m_orbitals->RnlID.push_back(m_nlms);
}
}
