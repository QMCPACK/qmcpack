//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef FITTING_H
#define FITTING_H

#include "Blitz.h"
#include "MatrixOps.h"

/// LitFit performs a least-squares fit to data given in y with the
/// errors given by sigma.  It performs a fit to a function of the
/// form \f[ y_{\mathrm{fit}}(x) \approx \sum_{j=0}^M a_j F_j(x) \f]. 
/// \f$ F_{ij} = F_j(x_i) \f$.  
void LinFitLU (Array<double,1> &y, Array<double,1> &sigma,   // inputs
	       Array<double,2> &F,                           // input
	       Array<double,1> &a, Array<double,1> &errors); // outputs

void LinFitSVD (Array<double,1> &y, Array<double,1> &sigma,   // inputs
		Array<double,2> &F,                           // input
		Array<double,1> &a, Array<double,1> &error,   // outputs
		double tolerance);

double LinFitSVD (Array<double,1> &y, Array<double,1> &sigma,   // inputs
		  Array<double,2> &F, Array<bool,1> &adjust,    // inputs
		  Array<double,1> &a, Array<double,1> &error,   // outputs
		  double tolerance);

#endif
