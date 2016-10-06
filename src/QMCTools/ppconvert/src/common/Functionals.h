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

#ifndef FUNCTIONALS_H
#define FUNCTIONALS_H

void ExchangePotential (double nup, double ndown,
			double &Vup, double &Vdown);
void CorrelationPotential(double  nup, double ndown,
			  double &Vup, double &Vdown);

void CPPExCorr(double nup, double ndown,
	       double &Vup, double &Vdown);

void FortranExCorr(double  nup, double  ndown,
		   double &Vup, double &Vdown);

void FortranExCorr(double n, double &Exc, double &Vxc);


double FortranXCE (double nup, double ndown);

#endif
