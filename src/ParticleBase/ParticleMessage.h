//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_PARTICLEMESSAGE_H
#define OHMMS_PARTICLEMESSAGE_H

template<class T> void gsum(T&, int) { }

#ifdef USE_MPI
template<>
void gsum(ParticleAttrib<TinyVector<double,3> >& f0, int gid) { 
  typedef ParticleAttrib<TinyVector<double,3> > ParticlePos_t;
  int nat =f0.size();
  ParticlePos_t fg(nat);
  MPI_Allreduce(&(f0[0][0]), &(fg[0][0]), nat*3, 
  	        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  f0 = fg;
}

template<>
void gsum(ParticleAttrib<double>& q, int gid) { 
  typedef ParticleAttrib<double> ParticleScalar_t;
  int nat =q.size();
  ParticleScalar_t qg(nat);
  MPI_Allreduce(&(q[0]), &(qg[0]), nat, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  q = qg;
}
#endif

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
