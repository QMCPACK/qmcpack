//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_COMMUNICATION_OPERATORS_H
#define OHMMS_COMMUNICATION_OPERATORS_H
#include <Message/Communicate.h>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/Tensor.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#if defined(HAVE_MPI)
#include <Message/CommOperatorsMPI.h>
#else
#include <Message/CommOperatorsSingle.h>
#endif
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
