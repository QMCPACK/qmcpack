#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <iostream>
#include "Optimize/Minimize.h"
//using namespace std;

double MCost (const gsl_vector *v, void *params)
{
  ConjugateGradient &CG = *(ConjugateGradient *)params;
  MinimizeFunction &MF = *CG.MinFunc;
  for (int i=0; i<MF.NumParams(); i++)
    MF.Params(i) = gsl_vector_get(v, i);
  return (MF.Cost());
}


void MGradCost(const gsl_vector *v, void *params, gsl_vector *df)
{
  typedef double scalar;
  ConjugateGradient &CG = *(ConjugateGradient *)params;
  MinimizeFunction &MF = *CG.MinFunc;
  double epsilon = CG.epsilon;  
  blitz::Array<scalar,1> gradient(MF.NumParams());
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
      gradient(i) = (CostPlus-CostMinus)/(2.0*epsilon);
      gsl_vector_set(df, i, (CostPlus-CostMinus)/(2.0*epsilon));
    }
  cout << "Gradient = " << endl;
  for (int i=0;i<MF.NumParams(); i++) cout << setw(14) << gradient(i) << endl;
}

void MBoth (const gsl_vector *v, void *params, double *f,
	    gsl_vector *df)
{
  *f = MCost(v, params);
  MGradCost(v, params, df);
}


void
ConjugateGradient::Minimize(MinimizeFunction &MinimFunc)
{
  typedef double scalar;
  size_t iter = 0;
  int status;

  MinFunc = &MinimFunc;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *Minimizer;

  gsl_multimin_function_fdf gslMinFunc;
 
  gslMinFunc.f = &MCost;
  gslMinFunc.df = &MGradCost;
  gslMinFunc.fdf = &MBoth;
  gslMinFunc.n = MinimFunc.NumParams();
  gslMinFunc.params = this;

  gsl_vector *StartPoint;
  StartPoint = gsl_vector_alloc(MinimFunc.NumParams());
  for (int i=0; i<MinimFunc.NumParams(); i++)
    gsl_vector_set(StartPoint, i, MinimFunc.Params(i));

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  Minimizer = gsl_multimin_fdfminimizer_alloc(T,MinimFunc.NumParams());

  gsl_multimin_fdfminimizer_set(Minimizer, &gslMinFunc, StartPoint,
				StepSize, Tolerance);
  
  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(Minimizer);
      
      if (status) {
	cout << "CG: Problem with iteration. " << endl;
	break;
      }
      
      status = gsl_multimin_test_gradient (Minimizer->gradient,
					   Tolerance);

      if (status == GSL_SUCCESS)
	cout << "Minimum found at:" << endl;

      //    fprintf (stderr, "%5d ", iter);
      cout << "CG: " << iter << endl;
      for (int i=0; i<MinimFunc.NumParams(); i++)
	cout << "CG: " << setw(15) <<  gsl_vector_get(Minimizer->x, i);
      cout << endl;
	//	fprintf (stderr, "%15.12e ", gsl_vector_get(Minimizer->x, i));
      //  fprintf (stderr, "\n");
      for (int i=0; i<MinimFunc.NumParams(); i++)
	MinimFunc.Params(i) = gsl_vector_get(Minimizer->x,i);
      MinimFunc.WriteStuff();
    }
  while (status == GSL_CONTINUE && iter < 100); // JNKIM reduces this number from 100

  for (int i=0; i<MinimFunc.NumParams(); i++)
    MinimFunc.Params(i) = gsl_vector_get(Minimizer->x,i);

  gsl_multimin_fdfminimizer_free (Minimizer);
  gsl_vector_free (StartPoint);

}

void
MinimizeFunction::ReadParameters(char *FileName)
{
  typedef double scalar;
  FILE *fin;
  if ((fin = fopen (FileName, "r")) == NULL)
    {
      std::cerr << "Can't open parmeters file.  Exitting.\n";
      exit(1);
    }

  for (int i=0; i<NumParams(); i++)
    {
      scalar temp;
      fscanf (fin, " %lf ", &temp);
      Params(i) = temp;
    }
}


class TestMinim : public MinimizeFunction
{
public:
  typedef double scalar;
  blitz::Array<scalar, 1> x;
  int NumParams()
  {
    return (2);
  }
  scalar &Params(int i)
  {
    return(x(i));
  }
  scalar Params(int i) const
  {
    return (x(i));
  }
  scalar Cost()
  {
    return ((x(0)-1.0)*(x(0)-1.0) + (x(1)-2.0)*(x(1)-2.0));
  }
};



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
