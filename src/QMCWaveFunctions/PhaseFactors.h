#ifndef PHASE_FACTORS_H
#define PHASE_FACTORS_H

#include<complex>

void apply_phase_factors(float kPoints[], int makeTwoCopies[],
                         float pos[], float *phi_in[], float *phi_out[],
                         int num_splines, int num_walkers);

void apply_phase_factors(float kPoints[], int makeTwoCopies[],
                         float pos[], float *phi_in[], float *phi_out[],
                         float *GL_in[], float *GL_out[],
                         int num_splines, int num_walkers, int row_stride);

void apply_phase_factors(float kPoints[], int makeTwoCopies[], int TwoCopiesIndex[],
                         float pos[], float *phi_in[], float *phi_out[],
                         float *GL_in[], float *GL_out[],
                         int num_splines, int num_walkers, int row_stride);

void apply_phase_factors(double kPoints[], int makeTwoCopies[],
                         double pos[], double *phi_in[], double *phi_out[],
                         int num_splines, int num_walkers);

void apply_phase_factors(double kPoints[], int makeTwoCopies[],
                         double pos[], double *phi_in[], double *phi_out[],
                         double *GL_in[], double *GL_out[],
                         int num_splines, int num_walkers, int row_stride);

void apply_phase_factors(double kPoints[], int makeTwoCopies[], int TwoCopiesIndex[],
                         double pos[], double *phi_in[], double *phi_out[],
                         double *GL_in[], double *GL_out[],
                         int num_splines, int num_walkers, int row_stride);


#ifdef QMC_COMPLEX
// YingWai's complex implementation
void apply_phase_factors(float kPoints[], float pos[],
                         std::complex<float>* phi_in[], std::complex<float>* phi_out[],
                         int num_splines, int num_walkers);

void apply_phase_factors(double kPoints[], double pos[],
                         std::complex<double>* phi_in[], std::complex<double>* phi_out[],
                         int num_splines, int num_walkers);

void apply_phase_factors(float kPoints[], float pos[],
                         std::complex<float>* phi_in[], std::complex<float>* phi_out[],
                         std::complex<float>* GL_in[], std::complex<float>* GL_out[],
                         int num_splines, int num_walkers, int row_stride);

void apply_phase_factors(double kPoints[], double pos[],
                         std::complex<double>* phi_in[], std::complex<double>* phi_out[],
                         std::complex<double>* GL_in[], std::complex<double>* GL_out[],
                         int num_splines, int num_walkers, int row_stride);
#endif


#endif
