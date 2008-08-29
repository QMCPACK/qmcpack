//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
// -*- C++ -*-
/** @file ECPComponentBuilderBuilder.h
 * @brief Declaration of a builder class for an ECP component for an ionic type
 */
#ifndef QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#define QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPotential.h"

namespace qmcplusplus {

  struct ECPComponentBuilder: public MPIObjectBase, public QMCTraits {

    typedef LocalECPotential::GridType GridType;
    typedef LocalECPotential::RadialPotentialType RadialPotentialType;

    int NumNonLocal;
    int Lmax, Llocal, Nrule;
    RealType Zeff;
    RealType RcutMax;
    string Species;
    GridType *grid_global;
    map<string,GridType*> grid_inp;
    RadialPotentialType* pp_loc;
    NonLocalECPComponent* pp_nonloc;
    map<string,int> angMon;

    ECPComponentBuilder(const string& aname, Communicate* c);

    bool parse(const string& fname, xmlNodePtr cur);
    bool put(xmlNodePtr cur);
    void addSemiLocal(xmlNodePtr cur);
    void buildLocal(xmlNodePtr cur);
    void buildSemiLocalAndLocal(vector<xmlNodePtr>& semiPtr);

    bool parseCasino(const string& fname, xmlNodePtr cur); //std::string& fname, RealType rc);
    //bool parseCasino(std::string& fname, RealType rc);
    // This sets the spherical quadrature rule used to apply the
    // projection operators.  rule can be 1 to 7.  See
    // J. Chem. Phys. 95 (3467) (1991)
    // Rule     # points     lexact
    //  1           1          0
    //  2           4          2
    //  3           6          3
    //  4          12          5
    //  5          18          5
    //  6          26          7
    //  7          50         11
    void SetQuadratureRule(int rule);
    void CheckQuadratureRule(int lexact);

    GridType* createGrid(xmlNodePtr cur, bool useLinear=false);
    RadialPotentialType* createVrWithBasisGroup(xmlNodePtr cur, GridType* agrid);
    RadialPotentialType* createVrWithData(xmlNodePtr cur, GridType* agrid, int rCorrection=0);

    void doBreakUp(const vector<int>& angList, const Matrix<RealType>& vnn, 
        RealType rmax, RealType Vprefactor=1.0);

    void printECPTable();
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
