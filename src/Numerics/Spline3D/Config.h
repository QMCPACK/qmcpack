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
    
    



#ifndef GUARD_SPLINE3DCONFIGURATION_H
#define GUARD_SPLINE3DCONFIGURATION_H

#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"

typedef TinyVector<int,3> gridvec_t;
typedef TinyVector<double,3> posvec_t;
typedef Vector<double> scalar_array_t;
typedef Vector<posvec_t> posarray_t;
typedef std::vector<posvec_t> rarray_t;

#endif
