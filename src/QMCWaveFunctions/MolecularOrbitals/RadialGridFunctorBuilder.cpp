//////////////////////////////////////////////////////////////////
// (c) Copyright 2003 by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
#include "Utilities/OhmmsInfo.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/SlaterTypeOrbital.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/HDFNumericAttrib.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/RadialGridFunctorBuilder.h"

namespace ohmmsqmc {

 /*!
  * \param cur  the current xmlNode to be processed
  * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
  * \return true is succeeds
  * \brief Add a new Slater Type Orbital with quantum numbers 
  \f$(n,l,m,s)\f$ to the list of radial orbitals.  This function 
  puts the STO on a logarithmic grid and calculates the boundary 
  conditions for the 1D Cubic Spline.  The derivates at the endpoint 
  are assumed to be all zero.  Note: for the radial orbital we use
  /f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
  radial orbital and \f$ l \f$ is the angular momentum.
 */
  bool
  STO2GridBuilder::addRadialOrbital(xmlNodePtr cur, 
				    const QuantumNumberType& nmls) {
    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    RadialOrbitalType *radorb =  NULL;
    int n=nmls[0];
    int l=nmls[1];
    ValueType zeta = 1.0;
    xmlNodePtr s = cur->xmlChildrenNode;
    while(s != NULL) {
      string cname((const char*)(s->name));
      if(cname == "parameter" || cname =="Var") {
	putContent(zeta,s);
      }
      s=s->next;
    }
    XMLReport("Zeta = " << zeta)
    STONorm<RealType> anorm(n);
    
    GenericSTO<RealType> sto(n-l-1,zeta,anorm(n-1,zeta));
    XMLReport("Calculating 1D-Cubic spline.")

    GridType* agrid = m_orbitals->Grids[0];

    radorb = new OneDimCubicSpline<ValueType>(agrid);
    Transform2GridFunctor<GenericSTO<RealType>,RadialOrbitalType> transform(sto, *radorb);
    transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());
    m_orbitals->Rnl.push_back(radorb);

#ifdef PRINT_DEBUG
    ofstream dfile("check.dat");
    for(int ig=0; ig<agrid->size()-1; ig++) {
      RealType dr = ((*agrid)(ig+1)-(*agrid)(ig))/5.0;
      RealType _r =(*agrid)(ig)+dr*0.00131,y,dy,d2y;
      while(_r<(*agrid)(ig+1)) {
	radorb->setgrid(_r);
	y = radorb->evaluate(_r,1.0/_r,dy,d2y);
	dfile << setw(15) << _r << setw(20) << setprecision(12) << y 
	      << setw(20) << dy << setw(20) << d2y
	      << endl;
	_r+=dr;
      }
    }
#endif

    return true;
  }

 /*! 
    * \param cur the current xmlNode to be processed
    * \return true if succeeds
    * \brief Add a radial grid to the list of radial grids.
    The Slater Type Orbitals are automatically placed on a 
    logarithmic grid and interpolated.
  */

  bool 
  STO2GridBuilder::addGrid(xmlNodePtr cur) {
    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    XMLReport("Converting Slater-type orbitals to radial grid functions")
    RealType ri = 1e-5;
    RealType rf = 100.0;
    IndexType npts = 1001;
    GridType *agrid = new LogGrid<RealType>;
    agrid->set(ri,rf,npts);
    LOGMSG("Using log grid with default values: ri = " << ri << " rf = " 
	   << rf << " npts = " << npts)

    m_orbitals->Grids.push_back(agrid);
    return true;
  }

  /*! \fn RGFBuilder::~RGFBuilder()
  * \brief the destructor
  */

  RGFBuilder::~RGFBuilder() {
    //clean up
    H5Gclose(m_group_id);   
    H5Fclose(m_file_id);   
  }

  /*!
  * \param cur the current xmlNode to be processed
  * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
  * \return true is succeeds
  * \brief Add a new radial orbital with quantum numbers 
  \f$(n,l,m,s)\f$to the list of radial orbitals.
  This function finds the location of the corresponding data for 
  the radial orbital in the HDF5 file, reads in the data into a new
  RadialOrbitalType and calculates the boundary conditions for the 
  1D Cubic Spline.  The derivates at the endpoint are assumed to be 
  all zero.  Note: for the radial orbital we use
  \f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
  radial orbital and \f$ l \f$ is the angular momentum.
  */

  bool
  RGFBuilder::addRadialOrbital(xmlNodePtr cur,
				const QuantumNumberType& nlms) {
    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    int n=nlms[0];
    int l=nlms[1];
    int m=nlms[2];
    int s=nlms[3];

    //first valid index for the radial grid
    int imin = 0;
    xmlAttrPtr att = cur->properties;
    while(att != NULL) {
      string aname((const char*)(att->name));
      if(aname == "imin"){
	imin = atoi((const char*)(att->children->content));
	XMLReport("First valid index for radial grid = " << imin)
	  }
      att = att->next;
    }

    //open the group containing the proper orbital
    char grpname[128];
    sprintf(grpname,"orbital%04d",Counter++);
    hid_t group_id_orb = H5Gopen(m_group_id,grpname);

    //find the correct power such that:
    //R(r)/r^l = rad_orb(r)/r^{power} 
    Vector<int> power;
    HDFAttribIO<Vector<int> > PowerHDFIn(power);
    PowerHDFIn.read(group_id_orb,"power");

    //check that the quantum numbers are consistent
    Vector<int> quantum_numbers;
    HDFAttribIO<Vector<int> > QuantumNoHDFIn(quantum_numbers);
    QuantumNoHDFIn.read(group_id_orb,"quantum_numbers");

    if((n != quantum_numbers[0]) || (l != quantum_numbers[1])){
      ERRORMSG("Inconsistency between input file and HDF5 file!")
	return false;
    }
    XMLReport("Reading from HDF5 file for Zeta = " << quantum_numbers[2])

      //read the radial orbital data
    Vector<double> rad_orb;
    HDFAttribIO<Vector<double> > RadOrbIn(rad_orb);
    RadOrbIn.read(group_id_orb,"radial_orbital");
    H5Gclose(group_id_orb);   

    //pointer to the grid
    GridType* agrid = m_orbitals->Grids[0];

    int rinv_p = power[0];    
    if(rinv_p != 0)
      for(int ig=0; ig<rad_orb.size(); ig++) 
        rad_orb[ig] *= pow(agrid->r(ig),-rinv_p);
    XMLReport("Multiplying orbital by r^" << -rinv_p)


    //last valid index for radial grid
    int imax = rad_orb.size()-1;
    RadialOrbitalType *radorb = new OneDimCubicSpline<ValueType>(agrid,rad_orb);
    //calculate boundary condition, assume derivates at endpoint are 0.0
    RealType yprime_i = rad_orb[imin+1]-rad_orb[imin];
    if(fabs(yprime_i)<1e-10)  yprime_i = 0.0;
    yprime_i /= (agrid->r(imin+1)-agrid->r(imin)); 
    //set up 1D-Cubic Spline
    radorb->spline(imin,yprime_i,imax,0.0);
    //add radial orbital to list
    m_orbitals->Rnl.push_back(radorb);

    XMLReport("Calculating 1D-Cubic spline, cusp condtion = " << yprime_i)

      //#ifdef PRINT_DEBUG   
    string fname(grpname);
    fname.append(".dat");
    ofstream dfile(fname.c_str());
    for(int ig=imin; ig<radorb->size(); ig++) {
      RealType dr = (radorb->r(ig+1)- radorb->r(ig))/5.0;
      RealType _r = radorb->r(ig),y,dy,d2y;
      while(_r<radorb->r(ig+1)) {
        //gorb->setgrid(_r);
	//y = gorb->evaluate(_r,1.0/_r,dy,d2y);
	radorb->setgrid(_r);
	//y = aterm->evaluate(_r,1.0/_r,dy,d2y);
	y = radorb->evaluate(_r,1.0/_r,dy,d2y);
	dfile << setw(15) << _r << setw(20) << setprecision(12) << y 
	      << setw(20) << dy << setw(20) << d2y
	      << endl;
	_r+=dr;
      }
    }
    //#endif

    return true;
  }

  /*!
    * \param cur the current xmlNode to be processed
    * \return true if succeeds
    * \brief Add a radial grid to the list of radial grids.
    This function opens the HDF5 file named in the input file, 
    finds the location of the radial grid data and creates the 
    grid.  The possible grid types are:
    Linear \f[ r_i = r_0 + 
    i\left( \frac{r_f - r_0}{N-1} \right) \f]
    Logarithmic \f[ r_i = r_0 
    \left( \frac{r_f}{r_0} \right) ^{\frac{i}{N-1}} \f]
    or LogZero \f[ r_i = B 
    \left[ \exp(Ai)-1 \right] \f]
  */

  bool 
  RGFBuilder::addGrid(xmlNodePtr cur) {

    if(!m_orbitals) {
      ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
      return false;
    }

    XMLReport("Reading a table for radial grid functions")
      if(!xmlHasProp(cur, (const xmlChar *)"file")) {
      ERRORMSG("Missing a file name for basis functions")
      return false;
    }

    //open the HDF5 file
    string afilename((const char*)(xmlGetProp(cur, (const xmlChar *)"file")));

    m_file_id = H5Fopen(afilename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    XMLReport("Opening file: " << afilename)
      m_group_id = H5Gopen(m_file_id,"radial_basis_states");
    hid_t group_id_grid = H5Gopen(m_group_id,"grid");
    
    //create a grid and initialize it
    GridType* agrid = NULL;
    Vector<int> grid_type;
    Vector<double> grid_params;

    HDFAttribIO<Vector<int> > GridTypeIn(grid_type);
    GridTypeIn.read(group_id_grid,"grid_type");

    HDFAttribIO<Vector<double> > GridParamsIn(grid_params);
    GridParamsIn.read(group_id_grid,"params");

    /* documentation */
    enum {LINEAR_GRID=1, LOG_GRID, LOG_GRID_ZERO};

    H5Gclose(group_id_grid);
    if(grid_type[0] == LINEAR_GRID){
      agrid = new LinearGrid<RealType>;
      LOGMSG("Using a linear grid.")
	RealType ri = grid_params(0);
      RealType rf = grid_params(1);
      IndexType npts = static_cast<int>(grid_params(2));
      agrid->set(ri,rf,npts);
    } else if(grid_type[0] == LOG_GRID) {
      agrid = new LogGrid<RealType>;
      LOGMSG("Using a log grid.")
	RealType ri = grid_params(0);
      RealType rf = grid_params(1);
      IndexType npts = static_cast<int>(grid_params(2));
      agrid->set(ri,rf,npts);
    } else if(grid_type[0] == LOG_GRID_ZERO){
      agrid = new LogGridZero<RealType>;
      LOGMSG("Using a log zero grid.")
	RealType a = grid_params(0);
      RealType b = grid_params(1);
      IndexType npts = static_cast<int>(grid_params(2));
      agrid->set(a,b,npts);
      
    }

    LOGMSG("Grid Values: ri = " << agrid->rmin() << " rf = "
	   << agrid->rmax() << " npts = " << agrid->size())
      //add the radial grid to the list   
      m_orbitals->Grids.push_back(agrid);
    
    return true;
  }
  
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
