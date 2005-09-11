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

//
// Other definitions
//

#include <mpi.h>
#include "oompi-config.h"
#include "Wrapper_ptr.cct"
#include "Datatype.h"
#include "Op.h"
#include "Hidden.h"
#include "Constants.h"


//
// Return types
//

const OOMPI_Error_type OOMPI_SUCCESS = MPI_SUCCESS;
const OOMPI_Error_type OOMPI_ERR_BUFFER = MPI_ERR_BUFFER;
const OOMPI_Error_type OOMPI_ERR_COUNT = MPI_ERR_COUNT;
const OOMPI_Error_type OOMPI_ERR_TYPE = MPI_ERR_TYPE;
const OOMPI_Error_type OOMPI_ERR_TAG = MPI_ERR_TAG;
const OOMPI_Error_type OOMPI_ERR_COMM= MPI_ERR_COMM;
const OOMPI_Error_type OOMPI_ERR_RANK = MPI_ERR_RANK;
const OOMPI_Error_type OOMPI_ERR_REQUEST = MPI_ERR_REQUEST;
const OOMPI_Error_type OOMPI_ERR_ROOT = MPI_ERR_ROOT;
const OOMPI_Error_type OOMPI_ERR_GROUP = MPI_ERR_GROUP;
const OOMPI_Error_type OOMPI_ERR_OP = MPI_ERR_OP;
const OOMPI_Error_type OOMPI_ERR_TOPOLOGY = MPI_ERR_TOPOLOGY;
const OOMPI_Error_type OOMPI_ERR_DIMS = MPI_ERR_DIMS;
const OOMPI_Error_type OOMPI_ERR_ARG = MPI_ERR_ARG;
const OOMPI_Error_type OOMPI_ERR_UNKNOWN = MPI_ERR_UNKNOWN;
const OOMPI_Error_type OOMPI_ERR_TRUNCATE = MPI_ERR_TRUNCATE;
const OOMPI_Error_type OOMPI_ERR_OTHER = MPI_ERR_OTHER;
const OOMPI_Error_type OOMPI_ERR_INTERN = MPI_ERR_INTERN;
#if LSC_MPI_ERR_PENDING
const OOMPI_Error_type OOMPI_ERR_PENDING = MPI_PENDING;
#else
const OOMPI_Error_type OOMPI_ERR_PENDING = MPI_ERR_PENDING;
#endif
const OOMPI_Error_type OOMPI_ERR_IN_STATUS = MPI_ERR_IN_STATUS;
const OOMPI_Error_type OOMPI_ERR_LASTCODE = MPI_ERR_LASTCODE;


OOMPI_Error_type OOMPI_errno = OOMPI_SUCCESS;
OOMPI_Error_type &OOMPI_Errno = OOMPI_errno;


//
// Default type tags
//

const int OOMPI_RESERVED_TAGS = 256;
const int OOMPI_TAG_UB = 32767 - OOMPI_RESERVED_TAGS;

const int OOMPI_CHAR_TAG = OOMPI_TAG_UB + 1;
const int OOMPI_SHORT_TAG = OOMPI_TAG_UB + 2;
const int OOMPI_INT_TAG = OOMPI_TAG_UB + 3;
const int OOMPI_LONG_TAG = OOMPI_TAG_UB + 4;
const int OOMPI_UNSIGNED_CHAR_TAG = OOMPI_TAG_UB + 5;
const int OOMPI_UNSIGNED_SHORT_TAG = OOMPI_TAG_UB + 6;
const int OOMPI_UNSIGNED_TAG = OOMPI_TAG_UB + 7;
const int OOMPI_UNSIGNED_LONG_TAG = OOMPI_TAG_UB + 8;
const int OOMPI_FLOAT_TAG = OOMPI_TAG_UB + 9;
const int OOMPI_DOUBLE_TAG = OOMPI_TAG_UB + 10;
#if OOMPI_HAVE_LONG_DOUBLE
const int OOMPI_LONG_DOUBLE_TAG = OOMPI_TAG_UB + 11;
#endif
#if OOMPI_HAVE_LONG_LONG_INT
const int OOMPI_LONG_LONG_INT_TAG = OOMPI_TAG_UB + 12;
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
const int OOMPI_COMPLEX_FLOAT_TAG = OOMPI_TAG_UB + 13;
const int OOMPI_COMPLEX_DOUBLE_TAG = OOMPI_TAG_UB + 14;
const int OOMPI_COMPLEX_LONG_DOUBLE_TAG = OOMPI_TAG_UB + 15;
#endif
const int OOMPI_BYTE_TAG = OOMPI_TAG_UB + 16;
const int OOMPI_MESSAGE_TAG = OOMPI_TAG_UB + 17;
const int OOMPI_PACKED_TAG = OOMPI_TAG_UB + 18;

const int OOMPI_MPI_DATATYPE_TAG = OOMPI_TAG_UB + 30;
const int OOMPI_INTERCOMM_CREATE_TAG = OOMPI_TAG_UB + 31;

const int OOMPI_NO_TAG = OOMPI_TAG_UB + 40;
const int OOMPI_NO_COUNT = OOMPI_TAG_UB + 41;


//
// Some OOMPI equivalents of MPI constants
//

const int OOMPI_ANY_TAG = MPI_ANY_TAG;
const int OOMPI_ANY_SOURCE = MPI_ANY_SOURCE;
const int OOMPI_PROC_NULL = MPI_PROC_NULL;
const int OOMPI_UNDEFINED = MPI_UNDEFINED;


//
// result of communicator and group comparisons
//

const OOMPI_Compare OOMPI_IDENT = MPI_IDENT;
const OOMPI_Compare OOMPI_CONGRUENT = MPI_CONGRUENT;
const OOMPI_Compare OOMPI_SIMILAR = MPI_SIMILAR;
const OOMPI_Compare OOMPI_UNEQUAL = MPI_UNEQUAL;


//
// Topologies
//

const OOMPI_Topology OOMPI_GRAPH = MPI_GRAPH;
const OOMPI_Topology OOMPI_CART = MPI_CART;


//
// Reduction operations
//

const OOMPI_Op OOMPI_SUM = OOMPI_Hidden::Create_op(MPI_SUM);
const OOMPI_Op OOMPI_MAX = OOMPI_Hidden::Create_op(MPI_MAX);
const OOMPI_Op OOMPI_MIN = OOMPI_Hidden::Create_op(MPI_MIN);
const OOMPI_Op OOMPI_PROD = OOMPI_Hidden::Create_op(MPI_PROD);
const OOMPI_Op OOMPI_MAXLOC = OOMPI_Hidden::Create_op(MPI_MAXLOC);
const OOMPI_Op OOMPI_MINLOC = OOMPI_Hidden::Create_op(MPI_MINLOC);
const OOMPI_Op OOMPI_BAND = OOMPI_Hidden::Create_op(MPI_BAND);
const OOMPI_Op OOMPI_BOR = OOMPI_Hidden::Create_op(MPI_BOR);
const OOMPI_Op OOMPI_BXOR = OOMPI_Hidden::Create_op(MPI_BXOR);
const OOMPI_Op OOMPI_LAND = OOMPI_Hidden::Create_op(MPI_LAND);
const OOMPI_Op OOMPI_LOR = OOMPI_Hidden::Create_op(MPI_LOR);
const OOMPI_Op OOMPI_LXOR = OOMPI_Hidden::Create_op(MPI_LXOR);


//
// OOMPI_Datatypes
//

const OOMPI_Datatype OOMPI_CHAR = 
  OOMPI_Hidden::Create_datatype(MPI_CHAR, OOMPI_CHAR_TAG);
const OOMPI_Datatype OOMPI_SHORT = 
  OOMPI_Hidden::Create_datatype(MPI_SHORT, OOMPI_SHORT_TAG);
const OOMPI_Datatype OOMPI_INT = 
  OOMPI_Hidden::Create_datatype(MPI_INT, OOMPI_INT_TAG);
const OOMPI_Datatype OOMPI_LONG = 
  OOMPI_Hidden::Create_datatype(MPI_LONG, OOMPI_LONG_TAG);
const OOMPI_Datatype OOMPI_UNSIGNED_CHAR = 
  OOMPI_Hidden::Create_datatype(MPI_UNSIGNED_CHAR, OOMPI_UNSIGNED_CHAR_TAG);
const OOMPI_Datatype OOMPI_UNSIGNED_SHORT = 
  OOMPI_Hidden::Create_datatype(MPI_UNSIGNED_SHORT, OOMPI_UNSIGNED_SHORT_TAG);
const OOMPI_Datatype OOMPI_UNSIGNED = 
  OOMPI_Hidden::Create_datatype(MPI_UNSIGNED, OOMPI_UNSIGNED_TAG);
const OOMPI_Datatype OOMPI_UNSIGNED_LONG = 
  OOMPI_Hidden::Create_datatype(MPI_UNSIGNED_LONG, OOMPI_UNSIGNED_LONG_TAG);
const OOMPI_Datatype OOMPI_FLOAT = 
  OOMPI_Hidden::Create_datatype(MPI_FLOAT, OOMPI_FLOAT_TAG);
const OOMPI_Datatype OOMPI_DOUBLE = 
  OOMPI_Hidden::Create_datatype(MPI_DOUBLE, OOMPI_DOUBLE_TAG);
#if OOMPI_HAVE_LONG_DOUBLE
const OOMPI_Datatype OOMPI_LONG_DOUBLE = 
  OOMPI_Hidden::Create_datatype(MPI_LONG_DOUBLE, OOMPI_LONG_DOUBLE_TAG);
#endif
#if OOMPI_HAVE_LONG_LONG_INT
const OOMPI_Datatype OOMPI_LONG_LONG_INT = 
  OOMPI_Hidden::Create_datatype(MPI_LONG_LONG_INT, OOMPI_LONG_LONG_INT_TAG);
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
const OOMPI_Datatype OOMPI_COMPLEX_FLOAT = 
  OOMPI_Hidden::Create_datatype(MPI_COMPLEX_FLOAT, OOMPI_COMPLEX_FLOAT_TAG);
const OOMPI_Datatype OOMPI_COMPLEX_DOUBLE = 
  OOMPI_Hidden::Create_datatype(MPI_COMPLEX_DOUBLE, OOMPI_COMPLEX_DOUBLE_TAG);
const OOMPI_Datatype OOMPI_COMPLEX_LONG_DOUBLE = 
  OOMPI_Hidden::Create_datatype(MPI_COMPLEX_LONG_DOUBLE, 
				OOMPI_COMPLEX_LONG_DOUBLE_TAG);
#endif
const OOMPI_Datatype OOMPI_PACKED = 
  OOMPI_Hidden::Create_datatype(MPI_PACKED, OOMPI_PACKED_TAG);
const OOMPI_Datatype OOMPI_MESSAGE;


//
// Miscellaneous
//

int OOMPI_HOST;
int OOMPI_IO;
bool OOMPI_WTIME_IS_GLOBAL;
const int OOMPI_MAX_ERROR_STRING = MPI_MAX_ERROR_STRING;
const int OOMPI_MAX_PROCESSOR_NAME = MPI_MAX_PROCESSOR_NAME;




