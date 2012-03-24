/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//                                                                         //
//  This program is free software; you can redistribute it and/or modify   //
//  it under the terms of the GNU General Public License as published by   //
//  the Free Software Foundation; either version 2 of the License, or      //
//  (at your option) any later version.                                    //
//                                                                         //
//  This program is distributed in the hope that it will be useful,        //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//  GNU General Public License for more details.                           //
//                                                                         //
//  You should have received a copy of the GNU General Public License      //
//  along with this program; if not, write to the Free Software            //
//  Foundation, Inc., 51 Franklin Street, Fifth Floor,                     //
//  Boston, MA  02110-1301  USA                                            //
/////////////////////////////////////////////////////////////////////////////

#ifndef BLIP_CREATE_H
#define BLIP_CREATE_H

#include "bspline_base.h"
#include "bspline_structs.h"
#include <stdbool.h>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Blip creation functions               ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

UBspline_3d_s*
create_blip_3d_s (double *lattice, double *Gvecs, 
		  complex_float *coefs, int numG,
		  double factor, bool useReal);

UBspline_3d_d*
create_blip_3d_d (double *lattice, double *Gvecs, 
		  complex_double *coefs, int numG,
		  double factor, bool useReal);

UBspline_3d_c*
create_blip_3d_c (double *lattice, double *Gvecs, 
		  complex_float *coefs, int numG,
		  double factor);

UBspline_3d_z*
create_blip_3d_z (double *lattice, double *Gvecs, 
		  complex_double *coefs, int numG,
		  double factor);



#endif
