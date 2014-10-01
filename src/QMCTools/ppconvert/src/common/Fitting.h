/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef FITTING_H
#define FITTING_H

#include "Blitz.h"
#include "MatrixOps.h"

/// LitFit performs a least-squares fit to data given in y with the
/// errors given by sigma.  It performs a fit to a function of the
/// form \f[ y_{\text{fit}}(x) \approx \sum_{j=0}^M a_j F_j(x) \f]. 
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
