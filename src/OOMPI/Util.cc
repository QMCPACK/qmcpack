// -*- c++ -*-
//
// Copyright (c) 2002-2003 Indiana University.  All rights reserved.
// Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
//                         All rights reserved.
// 
// This file is part of the OOMPI software package.  For license
// information, see the LICENSE file in the top level directory of the
// OOMPI source distribution.
//
// $Id$
//
// OOMPI Class library
// Utility functions
//

#include "oompi-config.h"
#include "Util.h"

int OOMPI_Util::util_sem = 1;
int OOMPI_Util::sem_array[SEM_MAX];

