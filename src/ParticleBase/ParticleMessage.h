//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_PARTICLEMESSAGE_H
#define OHMMS_PARTICLEMESSAGE_H

template<class T> void gsum(T&, int) { }

#ifdef USE_MPI
template<>
void gsum(ParticleAttrib<TinyVector<double,3> >& f0, int gid)
{
  typedef ParticleAttrib<TinyVector<double,3> > ParticlePos_t;
  int nat =f0.size();
  ParticlePos_t fg(nat);
  MPI_Allreduce(&(f0[0][0]), &(fg[0][0]), nat*3,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  f0 = fg;
}

template<>
void gsum(ParticleAttrib<double>& q, int gid)
{
  typedef ParticleAttrib<double> ParticleScalar_t;
  int nat =q.size();
  ParticleScalar_t qg(nat);
  MPI_Allreduce(&(q[0]), &(qg[0]), nat, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  q = qg;
}
#endif

#endif

