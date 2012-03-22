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

#ifndef MULTI_NUBSPLINE_H
#define MULTI_NUBSPLINE_H

#include "bspline_base.h"
#include "multi_nubspline_structs.h"

// Currently, some of the single-precision routines use SSE2 instructions
#include "multi_nubspline_eval_s.h"
#include "multi_nubspline_eval_c.h"
#include "multi_nubspline_eval_d.h"
#include "multi_nubspline_eval_z.h"

#include "nubspline_create.h"
#include "multi_nubspline_create.h"
#endif
