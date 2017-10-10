//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <cstdio>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <iomanip>
#include "Optimize/GSLMinimize.h"

double MCost (const gsl_vector *v, void *params)
{
  typedef CostFunctionBase<double> MinimizeFunction;
  GSLConjugateGradient &CG = *(GSLConjugateGradient *)params;
  MinimizeFunction &MF = *CG.MinFunc;
  for (int i=0; i<MF.NumParams(); i++)
    MF.Params(i) = gsl_vector_get(v, i);
  return (MF.Cost());
}


void MGradCost(const gsl_vector *v, void *params, gsl_vector *df)
{
  typedef double scalar;
  typedef CostFunctionBase<double> MinimizeFunction;
  GSLConjugateGradient &CG = *(GSLConjugateGradient *)params;
  MinimizeFunction &MF = *CG.MinFunc;
  double epsilon = CG.epsilon;
  std::vector<scalar> gradient(MF.NumParams());
  //blitz::Array<scalar,1> gradient(MF.NumParams());
  // std::vector<double> gradient;
  //gradient.resize(MF.NumParams());
  for (int i=0; i<MF.NumParams(); i++)
  {
    for (int j=0; j<MF.NumParams(); j++)
      MF.Params(j) = gsl_vector_get(v,j);
    scalar CostPlus, CostMinus;
    MF.Params(i) = gsl_vector_get(v,i) + epsilon;
    CostPlus = MF.Cost();
    MF.Params(i) = gsl_vector_get(v,i) - epsilon;
    CostMinus = MF.Cost();
    gradient[i] = (CostPlus-CostMinus)/(2.0*epsilon);
    gsl_vector_set(df, i, (CostPlus-CostMinus)/(2.0*epsilon));
  }
}

void MBoth (const gsl_vector *v, void *params, double *f,
            gsl_vector *df)
{
  *f = MCost(v, params);
  MGradCost(v, params, df);
}


bool
GSLConjugateGradient::optimize(ObjectFuncType *atarget)
{
  typedef double scalar;
  size_t iter = 0;
  int status;
  MinFunc = atarget;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *Minimizer;
  gsl_multimin_function_fdf gslMinFunc;
  gslMinFunc.f = &MCost;
  gslMinFunc.df = &MGradCost;
  gslMinFunc.fdf = &MBoth;
  gslMinFunc.n = atarget->NumParams();
  gslMinFunc.params = this;
  gsl_vector *StartPoint;
  StartPoint = gsl_vector_alloc(atarget->NumParams());
  for (int i=0; i<atarget->NumParams(); i++)
    gsl_vector_set(StartPoint, i, atarget->Params(i));
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  Minimizer = gsl_multimin_fdfminimizer_alloc(T,atarget->NumParams());
  gsl_multimin_fdfminimizer_set(Minimizer, &gslMinFunc, StartPoint, StepSize, Tolerance);
  do
  {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(Minimizer);
    if (status)
    {
      if(msg_stream)
        *msg_stream << "CG: Problem with iteration. " << std::endl;
      break;
    }
    status = gsl_multimin_test_gradient (Minimizer->gradient, Tolerance);
    if (status == GSL_SUCCESS && msg_stream)
      *msg_stream << "Minimum found at:" << std::endl;
    for (int i=0; i<atarget->NumParams(); i++)
      atarget->Params(i) = gsl_vector_get(Minimizer->x,i);
    atarget->Report();
  }
  while (status == GSL_CONTINUE && iter < MaxCGStep && atarget->IsValid);
  for (int i=0; i<atarget->NumParams(); i++)
    atarget->Params(i) = gsl_vector_get(Minimizer->x,i);
  gsl_multimin_fdfminimizer_free (Minimizer);
  gsl_vector_free (StartPoint);
  return true;
}

GSLConjugateGradient::GSLConjugateGradient():
  MaxCGStep(100), epsilon(1e-6), StepSize(1e-3),
  Tolerance(1.e-6)
{
}

bool GSLConjugateGradient::put(xmlNodePtr cur)
{
  ParameterSet p;
  p.add(epsilon,"epsilon","scalar");
  p.add(StepSize,"stepsize","scalar");
  p.add(Tolerance,"tolerance","scalar");
  p.add(MaxCGStep,"max_steps","int");
  p.put(cur);
  return true;
}

//void
//MinimizeFunction::ReadParameters(char *FileName)
//{
//  typedef double scalar;
//  FILE *fin;
//  if ((fin = fopen (FileName, "r")) == NULL)
//    {
//      std::cerr << "Can't open parmeters file.  Exitting.\n";
//      exit(1);
//    }
//
//  for (int i=0; i<NumParams(); i++)
//    {
//      scalar temp;
//      fscanf (fin, " %lf ", &temp);
//      Params(i) = temp;
//    }
//}

//class TestMinim : public CostFunctionBase<double>
//{
//public:
//  typedef double scalar;
//  std::vector<scalar> x;
//  int NumParams()
//  {
//    return (2);
//  }
//  scalar &Params(int i)
//  {
//    return(x[i]);
//  }
//  scalar Params(int i) const
//  {
//    return (x[i]);
//  }
//  scalar Cost()
//  {
//    return ((x[0]-1.0)*(x[0]-1.0) + (x[1]-2.0)*(x[1]-2.0));
//  }
//};
//  main()
//  {
//    TestMinim Minim;

//    Minim.Tolerance = 1.0e-10;
//    Minim.StepSize = 0.01;
//    Minim.epsilon = 1.0e-6;
//    Minim.x.resize(2);
//    Minim.x(0) = 3.19;
//    Minim.x(1) = -2.56;

//    Minim.Minimize();
//  }
