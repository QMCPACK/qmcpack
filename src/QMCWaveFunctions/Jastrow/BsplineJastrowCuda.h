#ifndef BSPLINE_JASTROW_CUDA_H
#define BSPLINE_JASTROW_CUDA_H

template <typename S>
struct NLjobGPU
{
  int Elec, NumQuadPoints;
  S *R, *QuadPoints, *Ratios;
};

///////////////////////
// Two-Body routines //
///////////////////////

void
two_body_sum (float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
	      float spline_coefs[], int numCoefs, float rMax,  
	      float lattice[], float latticeInv[], float sum[], int numWalkers);

void
two_body_sum (double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
	      double spline_coefs[], int numCoefs, double rMax,  
	      double lattice[], double latticeInv[], double sum[], int numWalkers);

void
two_body_ratio (float *R[], int first, int last,
		float Rnew[], int inew,
		float spline_coefs[], int numCoefs, float rMax,  
		float lattice[], float latticeInv[], float sum[], int numWalkers);

void
two_body_ratio (double *R[], int first, int last,
		double Rnew[], int inew,
		double spline_coefs[], int numCoefs, double rMax,  
		double lattice[], double latticeInv[], double sum[], int numWalkers);

void
two_body_ratio_grad(float *R[], int first, int last,
		    float  Rnew[], int inew,
		    float spline_coefs[], int numCoefs, float rMax,  
		    float lattice[], float latticeInv[], bool zero,
		    float ratio_grad[], int numWalkers, bool use_fast_image);

void
two_body_ratio_grad(double *R[], int first, int last,
		    double  Rnew[], int inew,
		    double spline_coefs[], int numCoefs, double rMax,  
		    double lattice[], double latticeInv[], bool zero,
		    double ratio_grad[], int numWalkers);

void
two_body_NLratios(NLjobGPU<float> jobs[], int first, int last,
		  float* spline_coefs[], int numCoefs[], float rMax[], 
		  float lattice[], float latticeInv[], float sim_cell_radius,
		  int numjobs);

void
two_body_NLratios(NLjobGPU<double> jobs[], int first, int last,
		  double* spline_coefs[], int numCoefs[], double rMax[], 
		  double lattice[], double latticeInv[], 
		  double sim_cell_radius, int numjobs);


void
two_body_update(float *R[], int N, int iat, int numWalkers);

void
two_body_update(double *R[], int N, int iat, int numWalkers);


void
two_body_grad_lapl(float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
		   float spline_coefs[], int numCoefs, float rMax,  
		   float lattice[], float latticeInv[], float sim_cell_radius,
		   float gradLapl[], int row_stride, int numWalkers);

void
two_body_grad_lapl(double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
		   double spline_coefs[], int numCoefs, double rMax,  
		   double lattice[], double latticeInv[], 
		   double gradLapl[], int row_stride, int numWalkers);


void
two_body_gradient (float *R[], int first, int last, int iat, 
		   float spline_coefs[], int numCoefs, float rMax,
		   float lattice[], float latticeInv[], float sim_cell_radius,
		   bool zeroOut, float grad[], int numWalkers);

void
two_body_gradient (double *R[], int first, int last, int iat, 
		   double spline_coefs[], int numCoefs, double rMax,
		   double lattice[], double latticeInv[], bool zeroOut,
		   double grad[], int numWalkers);

void
two_body_derivs(float *R[], float *gradLogPsi[], int e1_first, int e1_last, 
		int e2_first, int e2_last,
		int numCoefs, float rMax,  
		float lattice[], float latticeInv[], float sim_cell_radius,
		float *derivs[], int numWalkers);
void
two_body_derivs(double *R[], double *gradLogPsi[], int e1_first, int e1_last, 
		int e2_first, int e2_last,
		int numCoefs, double rMax,  
		double lattice[], double latticeInv[], double sim_cell_radius,
		double *derivs[], int numWalkers);


///////////////////////
// One-Body routines //
///////////////////////

void
one_body_sum (float C[], float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
	      float spline_coefs[], int numCoefs, float rMax,  
	      float lattice[], float latticeInv[], float sum[], int numWalkers);

void
one_body_sum (double C[], double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
	      double spline_coefs[], int numCoefs, double rMax,  
	      double lattice[], double latticeInv[], double sum[], int numWalkers);

void
one_body_ratio (float C[], float *R[], int first, int last,
		float Rnew[], int inew,
		float spline_coefs[], int numCoefs, float rMax,  
		float lattice[], float latticeInv[], float sum[], int numWalkers);

void
one_body_ratio (double C[], double *R[], int first, int last,
		double Rnew[], int inew,
		double spline_coefs[], int numCoefs, double rMax,  
		double lattice[], double latticeInv[], double sum[], int numWalkers);

void
one_body_ratio_grad (float C[], float *R[], int first, int last,
		     float Rnew[], int inew,
		     float spline_coefs[], int numCoefs, float rMax,  
		     float lattice[], float latticeInv[], bool zero,
		     float ratio_grad[], int numWalkers, bool use_fast_image);
void
one_body_ratio_grad (double C[], double *R[], int first, int last,
		     double Rnew[], int inew,
		     double spline_coefs[], int numCoefs, double rMax,  
		     double lattice[], double latticeInv[], bool zero,
		     double ratio_grad[], int numWalkers);

void
one_body_NLratios(NLjobGPU<float> jobs[], float C[], int first, int last,
		  float spline_coefs[], int numCoefs, float rMax, 
		  float lattice[], float latticeInv[], int numjobs);

void
one_body_NLratios(NLjobGPU<float> jobs[], float C[], int first, int last,
		  float spline_coefs[], int numCoefs, float rMax, 
		  float lattice[], float latticeInv[], float sim_cell_radius,
		  int numjobs);


void
one_body_NLratios(NLjobGPU<double> jobs[], double C[], int first, int last,
		  double spline_coefs[], int numCoefs, double rMax, 
		  double lattice[], double latticeInv[], int numjobs);

void
one_body_update(float *R[], int N, int iat, int numWalkers);

void
one_body_update(double *R[], int N, int iat, int numWalkers);


void
one_body_grad_lapl(float C[], float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
		   float spline_coefs[], int numCoefs, float rMax,  
		   float lattice[], float latticeInv[], 
		   float gradLapl[], int row_stride, int numWalkers);

void
one_body_grad_lapl(double C[], double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
		   double spline_coefs[], int numCoefs, double rMax,  
		   double lattice[], double latticeInv[], 
		   double gradLapl[], int row_stride, int numWalkers);

void
one_body_gradient (float *Rlist[], int iat, float C[], int first, int last,
		   float spline_coefs[], int num_coefs, float rMax,
		   float L[], float Linv[], float sim_cell_radius,
		   bool zeroSum, float grad[], int numWalkers);

void
one_body_gradient (double *Rlist[], int iat, double C[], int first, int last,
		   double spline_coefs[], int num_coefs, double rMax,
		   double L[], double Linv[], bool zeroSum,
		   double grad[], int numWalkers);


void
one_body_derivs(float C[], float *R[], float *gradLogPsi[], 
		int cfirst, int clast, 
		int efirst, int elast,
		int numCoefs, float rMax,  
		float lattice[], float latticeInv[], float sim_cell_radius,
		float *derivs[], int numWalkers);

void
one_body_derivs(double C[], double *R[], double *gradLogPsi[], 
		int cfirst, int clast, 
		int efirst, int elast,
		int numCoefs, double rMax,  
		double lattice[], double latticeInv[], double sim_cell_radius,
		double *derivs[], int numWalkers);

#endif
