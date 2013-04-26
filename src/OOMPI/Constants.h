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
// OOMPI Constants
// Constants used by OOMPI
//


#ifndef _OOMPI_CONSTANTS_H_
#define _OOMPI_CONSTANTS_H_

//
// Other definitions
//

#include <mpi.h>
#include <limits.h>
#include "oompi-config.h"


//
// OOMPI error actions
//

typedef enum
{
  OOMPI_ERRORS_ARE_FATAL,
  OOMPI_ERRORS_EXCEPTION,
  OOMPI_ERRORS_RETURN
} OOMPI_Error_action;


//
// Environment object
//

class OOMPI_Environment;
extern OOMPI_Environment OOMPI_ENV;


//
// Return types
//

typedef int OOMPI_Error_type;

extern const OOMPI_Error_type OOMPI_SUCCESS;
extern const OOMPI_Error_type OOMPI_ERR_BUFFER;
extern const OOMPI_Error_type OOMPI_ERR_COUNT;
extern const OOMPI_Error_type OOMPI_ERR_TYPE;
extern const OOMPI_Error_type OOMPI_ERR_TAG;
extern const OOMPI_Error_type OOMPI_ERR_COMM;
extern const OOMPI_Error_type OOMPI_ERR_RANK;
extern const OOMPI_Error_type OOMPI_ERR_REQUEST;
extern const OOMPI_Error_type OOMPI_ERR_ROOT;
extern const OOMPI_Error_type OOMPI_ERR_GROUP;
extern const OOMPI_Error_type OOMPI_ERR_OP;
extern const OOMPI_Error_type OOMPI_ERR_TOPOLOGY;
extern const OOMPI_Error_type OOMPI_ERR_DIMS;
extern const OOMPI_Error_type OOMPI_ERR_ARG;
extern const OOMPI_Error_type OOMPI_ERR_UNKNOWN;
extern const OOMPI_Error_type OOMPI_ERR_TRUNCATE;
extern const OOMPI_Error_type OOMPI_ERR_OTHER;
extern const OOMPI_Error_type OOMPI_ERR_INTERN;
extern const OOMPI_Error_type OOMPI_ERR_PENDING;
extern const OOMPI_Error_type OOMPI_ERR_IN_STATUS;
extern const OOMPI_Error_type OOMPI_ERR_LASTCODE;

extern OOMPI_Error_type OOMPI_errno;
extern OOMPI_Error_type &OOMPI_Errno;


//
// Error Handling Specifiers
//

typedef MPI_Errhandler OOMPI_Errhandler;

class OOMPI_Error_table;
extern OOMPI_Error_table OOMPI_ERROR;


//
// Default type tags
//

extern const int OOMPI_RESERVED_TAGS;
extern const int OOMPI_TAG_UB;

extern const int OOMPI_CHAR_TAG;
extern const int OOMPI_SHORT_TAG;
extern const int OOMPI_INT_TAG;
extern const int OOMPI_LONG_TAG;
extern const int OOMPI_UNSIGNED_CHAR_TAG;
extern const int OOMPI_UNSIGNED_SHORT_TAG;
extern const int OOMPI_UNSIGNED_TAG;
extern const int OOMPI_UNSIGNED_LONG_TAG;
extern const int OOMPI_FLOAT_TAG;
extern const int OOMPI_DOUBLE_TAG;
#if OOMPI_HAVE_LONG_DOUBLE
extern const int OOMPI_LONG_DOUBLE_TAG;
#endif
#if OOMPI_HAVE_LONG_LONG_INT
extern const int OOMPI_LONG_LONG_INT_TAG;
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
extern const int OOMPI_COMPLEX_FLOAT_TAG;
extern const int OOMPI_COMPLEX_DOUBLE_TAG;
extern const int OOMPI_COMPLEX_LONG_DOUBLE_TAG;
#endif
extern const int OOMPI_BYTE_TAG;
extern const int OOMPI_MESSAGE_TAG;
extern const int OOMPI_PACKED_TAG;

extern const int OOMPI_MPI_DATATYPE_TAG;
extern const int OOMPI_INTERCOMM_CREATE_TAG;

extern const int OOMPI_NO_TAG;
extern const int OOMPI_NO_COUNT;


//
// Some OOMPI equivalents of MPI constants
//

extern const int OOMPI_ANY_TAG;
extern const int OOMPI_ANY_SOURCE;
extern const int OOMPI_PROC_NULL;
extern const int OOMPI_UNDEFINED;


//
// result of communicator and group comparisons
//

typedef int OOMPI_Compare;

extern const OOMPI_Compare OOMPI_IDENT;
extern const OOMPI_Compare OOMPI_CONGRUENT;
extern const OOMPI_Compare OOMPI_SIMILAR;
extern const OOMPI_Compare OOMPI_UNEQUAL;


//
// Topologies
//

typedef int OOMPI_Topology;

extern const OOMPI_Topology OOMPI_GRAPH;
extern const OOMPI_Topology OOMPI_CART;


//
// Reduction operations
//

class OOMPI_Op;

extern const OOMPI_Op OOMPI_MAX;
extern const OOMPI_Op OOMPI_MIN;
extern const OOMPI_Op OOMPI_SUM;
extern const OOMPI_Op OOMPI_PROD;
extern const OOMPI_Op OOMPI_MAXLOC;
extern const OOMPI_Op OOMPI_MINLOC;
extern const OOMPI_Op OOMPI_BAND;
extern const OOMPI_Op OOMPI_BOR;
extern const OOMPI_Op OOMPI_BXOR;
extern const OOMPI_Op OOMPI_LAND;
extern const OOMPI_Op OOMPI_LOR;
extern const OOMPI_Op OOMPI_LXOR;


//
// OOMPI_Datatypes
//

class OOMPI_Datatype;

extern const OOMPI_Datatype OOMPI_CHAR;
extern const OOMPI_Datatype OOMPI_SHORT;
extern const OOMPI_Datatype OOMPI_INT;
extern const OOMPI_Datatype OOMPI_LONG;
extern const OOMPI_Datatype OOMPI_UNSIGNED_CHAR;
extern const OOMPI_Datatype OOMPI_UNSIGNED_SHORT;
extern const OOMPI_Datatype OOMPI_UNSIGNED;
extern const OOMPI_Datatype OOMPI_UNSIGNED_LONG;
extern const OOMPI_Datatype OOMPI_FLOAT;
extern const OOMPI_Datatype OOMPI_DOUBLE;
extern const OOMPI_Datatype OOMPI_LONG_DOUBLE;
extern const OOMPI_Datatype OOMPI_BYTE;
#if OOMPI_HAVE_LONG_DOUBLE
extern const OOMPI_Datatype OOMPI_LONG_DOUBLE;
#endif
#if OOMPI_HAVE_LONG_LONG_INT
extern const OOMPI_Datatype OOMPI_LONG_DOUBLE;
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
extern const OOMPI_Datatype OOMPI_COMPLEX_FLOAT;
extern const OOMPI_Datatype OOMPI_COMPLEX_DOUBLE;
extern const OOMPI_Datatype OOMPI_COMPLEX_LONG_DOUBLE;
#endif
extern const OOMPI_Datatype OOMPI_PACKED;
extern const OOMPI_Datatype OOMPI_MESSAGE;


//
// Miscellaneous
//

typedef MPI_Aint OOMPI_Aint;

extern int OOMPI_HOST;
extern int OOMPI_IO;
extern bool OOMPI_WTIME_IS_GLOBAL;
extern const int OOMPI_MAX_ERROR_STRING;
extern const int OOMPI_MAX_PROCESSOR_NAME;


#endif
