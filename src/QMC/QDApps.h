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
#ifndef OHMMS_QMC_QDAPP_H
#define OHMMS_QMC_QDAPP_H
#include "QMC/QMCApps.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/Config.h"
#include <string>
#include <vector>

namespace ohmmsqmc {


  class QDApps: public QMCApps {

    int idir;
    string v_file;

    double ulength;      /// length factor for conversion to a_B*
    double epsbym;
    double mbyepsq;

    double minv_ref;     /// inverse of effective mass in quantum region
    double eps_ref;      /// dielectric constant in reference region

    std::vector<int> m_M,mat_m,prior_m;
    std::vector<double> inv_meff,e_gap,eps_m,phi_s;

    void setRefParams();

  public:

    ///constructor
    QDApps(int argc, char** argv);
    ~QDApps();
    ///initialization with a file
    bool init();

  protected:
    bool setParticleSets(xmlpp::Node*);
    bool setElectrons(xmlpp::Node*);
    bool setIons(xmlpp::Node*);
    bool setWavefunctions(xmlpp::Node*);
    bool setHamiltonian(xmlpp::Node*);  
    bool setSimulation(xmlpp::Node*); 
    bool setQDot(xmlpp::Node*);
    bool rearrange(std::vector<int>&);
    bool rearrange(std::vector<double>&);

    void elocal(const gridvec_t&);
    bool readGrid(xmlpp::Node*);

    Grid3D* DeviceGrid;

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
