#ifndef BSPLINE_JASTROW_CUDA_PBC_H
#define BSPLINE_JASTROW_CUDA_PBC_H

#include "NLjobGPU.h"

///////////////////////
// Two-Body routines //
///////////////////////

void
two_body_sum_PBC (float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                  float spline_coefs[], int numCoefs, float rMax,
                  float lattice[], float latticeInv[], float sum[], int numWalkers);

void
two_body_sum_PBC (double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                  double spline_coefs[], int numCoefs, double rMax,
                  double lattice[], double latticeInv[], double sum[], int numWalkers);

void
two_body_ratio_PBC (float *R[], int first, int last,
                    float Rnew[], int inew,
                    float spline_coefs[], int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sum[], int numWalkers);

void
two_body_ratio_PBC (double *R[], int first, int last,
                    double Rnew[], int inew,
                    double spline_coefs[], int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sum[], int numWalkers);

void
two_body_ratio_grad_PBC(float *R[], int first, int last,
                        float  Rnew[], int inew,
                        float spline_coefs[], int numCoefs, float rMax,
                        float lattice[], float latticeInv[], bool zero,
                        float ratio_grad[], int numWalkers, bool use_fast_image);

void
two_body_ratio_grad_PBC(double *R[], int first, int last,
                        double  Rnew[], int inew,
                        double spline_coefs[], int numCoefs, double rMax,
                        double lattice[], double latticeInv[], bool zero,
                        double ratio_grad[], int numWalkers);

void
two_body_NLratios_PBC(NLjobGPU<float> jobs[], int first, int last,
                      float* spline_coefs[], int numCoefs[], float rMax[],
                      float lattice[], float latticeInv[], float sim_cell_radius,
                      int numjobs);

void
two_body_NLratios_PBC(NLjobGPU<double> jobs[], int first, int last,
                      double* spline_coefs[], int numCoefs[], double rMax[],
                      double lattice[], double latticeInv[],
                      double sim_cell_radius, int numjobs);


void
two_body_update(float *R[], int N, int iat, int numWalkers);

void
two_body_update(double *R[], int N, int iat, int numWalkers);


void
two_body_grad_lapl_PBC(float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                       float spline_coefs[], int numCoefs, float rMax,
                       float lattice[], float latticeInv[], float sim_cell_radius,
                       float gradLapl[], int row_stride, int numWalkers);

void
two_body_grad_lapl_PBC(double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                       double spline_coefs[], int numCoefs, double rMax,
                       double lattice[], double latticeInv[],
                       double gradLapl[], int row_stride, int numWalkers);


void
two_body_gradient_PBC (float *R[], int first, int last, int iat,
                       float spline_coefs[], int numCoefs, float rMax,
                       float lattice[], float latticeInv[], float sim_cell_radius,
                       bool zeroOut, float grad[], int numWalkers);

void
two_body_gradient_PBC (double *R[], int first, int last, int iat,
                       double spline_coefs[], int numCoefs, double rMax,
                       double lattice[], double latticeInv[], bool zeroOut,
                       double grad[], int numWalkers);

void
two_body_derivs_PBC(float *R[], float *gradLogPsi[], int e1_first, int e1_last,
                    int e2_first, int e2_last,
                    int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sim_cell_radius,
                    float *derivs[], int numWalkers);
void
two_body_derivs_PBC(double *R[], double *gradLogPsi[], int e1_first, int e1_last,
                    int e2_first, int e2_last,
                    int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sim_cell_radius,
                    double *derivs[], int numWalkers);


///////////////////////
// One-Body routines //
///////////////////////

void
one_body_sum_PBC (float C[], float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                  float spline_coefs[], int numCoefs, float rMax,
                  float lattice[], float latticeInv[], float sum[], int numWalkers);

void
one_body_sum_PBC (double C[], double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                  double spline_coefs[], int numCoefs, double rMax,
                  double lattice[], double latticeInv[], double sum[], int numWalkers);

void
one_body_ratio_PBC (float C[], float *R[], int first, int last,
                    float Rnew[], int inew,
                    float spline_coefs[], int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sum[], int numWalkers);

void
one_body_ratio_PBC (double C[], double *R[], int first, int last,
                    double Rnew[], int inew,
                    double spline_coefs[], int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sum[], int numWalkers);

void
one_body_ratio_grad_PBC (float C[], float *R[], int first, int last,
                         float Rnew[], int inew,
                         float spline_coefs[], int numCoefs, float rMax,
                         float lattice[], float latticeInv[], bool zero,
                         float ratio_grad[], int numWalkers, bool use_fast_image);
void
one_body_ratio_grad_PBC (double C[], double *R[], int first, int last,
                         double Rnew[], int inew,
                         double spline_coefs[], int numCoefs, double rMax,
                         double lattice[], double latticeInv[], bool zero,
                         double ratio_grad[], int numWalkers);

void
one_body_NLratios_PBC(NLjobGPU<float> jobs[], float C[], int first, int last,
                      float spline_coefs[], int numCoefs, float rMax,
                      float lattice[], float latticeInv[], int numjobs);

void
one_body_NLratios_PBC(NLjobGPU<float> jobs[], float C[], int first, int last,
                      float spline_coefs[], int numCoefs, float rMax,
                      float lattice[], float latticeInv[], float sim_cell_radius,
                      int numjobs);


void
one_body_NLratios_PBC(NLjobGPU<double> jobs[], double C[], int first, int last,
                      double spline_coefs[], int numCoefs, double rMax,
                      double lattice[], double latticeInv[], int numjobs);

void
one_body_update(float *R[], int N, int iat, int numWalkers);

void
one_body_update(double *R[], int N, int iat, int numWalkers);


void
one_body_grad_lapl_PBC(float C[], float *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                       float spline_coefs[], int numCoefs, float rMax,
                       float lattice[], float latticeInv[],
                       float gradLapl[], int row_stride, int numWalkers);

void
one_body_grad_lapl_PBC(double C[], double *R[], int e1_first, int e1_last, int e2_first, int e2_last,
                       double spline_coefs[], int numCoefs, double rMax,
                       double lattice[], double latticeInv[],
                       double gradLapl[], int row_stride, int numWalkers);

void
one_body_gradient_PBC (float *Rlist[], int iat, float C[], int first, int last,
                       float spline_coefs[], int num_coefs, float rMax,
                       float L[], float Linv[], float sim_cell_radius,
                       bool zeroSum, float grad[], int numWalkers);

void
one_body_gradient_PBC (double *Rlist[], int iat, double C[], int first, int last,
                       double spline_coefs[], int num_coefs, double rMax,
                       double L[], double Linv[], bool zeroSum,
                       double grad[], int numWalkers);


void
one_body_derivs_PBC(float C[], float *R[], float *gradLogPsi[],
                    int cfirst, int clast,
                    int efirst, int elast,
                    int numCoefs, float rMax,
                    float lattice[], float latticeInv[], float sim_cell_radius,
                    float *derivs[], int numWalkers);

void
one_body_derivs_PBC(double C[], double *R[], double *gradLogPsi[],
                    int cfirst, int clast,
                    int efirst, int elast,
                    int numCoefs, double rMax,
                    double lattice[], double latticeInv[], double sim_cell_radius,
                    double *derivs[], int numWalkers);

#endif
