#ifndef CUDA_NLPP_H
#define CUDA_NLPP_H

void
find_core_electrons_PBC (float *R[], int numElec,
                         float I[], int firstIon, int lastIon,
                         float rcut, float L[], float Linv[],
                         int2 *pairs[], float *dist[],
                         int numPairs[], int numWalkers);
void
find_core_electrons_PBC (double *R[], int numElec,
                         double I[], int firstIon, int lastIon,
                         double rcut, double L[], double Linv[],
                         int2 *pairs[], double *dist[],
                         int numPairs[], int numWalkers);

void
find_core_electrons_PBC (float *R[], int numElec,
                         float I[], int firstIon, int lastIon,
                         float rcut, float L[], float Linv[],
                         float quadPoints[], int numQuadPoints,
                         int *elecs[], float *ratioPos[],
                         float *dist_list[], float *cosTheta_list[],
                         int numPairs[], int numWalkers);

void
find_core_electrons_PBC (double *R[], int numElec,
                         double I[], int firstIon, int lastIon,
                         double rcut, double L[], double Linv[],
                         double quadPoints[], int numQuadPoints,
                         int *elecs[], double *ratioPos[],
                         double *dist_list[], double *cosTheta_list[],
                         int numPairs[], int numWalkers);



void
find_core_electrons (float *R[], int numElec,
                     float I[], int firstIon, int lastIon,
                     float rcut,
                     int2 *pairs[], float *dist[],
                     int numPairs[], int numWalkers);
void
find_core_electrons (double *R[], int numElec,
                     double I[], int firstIon, int lastIon,
                     double rcut, double L[], double Linv[],
                     int2 *pairs[], double *dist[],
                     int numPairs[], int numWalkers);

void
find_core_electrons (float *R[], int numElec,
                     float I[], int firstIon, int lastIon,
                     float rcut,
                     float quadPoints[], int numQuadPoints,
                     int *elecs[], float *ratioPos[],
                     float *dist_list[], float *cosTheta_list[],
                     int numPairs[], int numWalkers);

void
find_core_electrons (double *R[], int numElec,
                     double I[], int firstIon, int lastIon,
                     double rcut,
                     double quadPoints[], int numQuadPoints,
                     int *elecs[], double *ratioPos[],
                     double *dist_list[], double *cosTheta_list[],
                     int numPairs[], int numWalkers);

#endif
