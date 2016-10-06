//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20
#define UMFPACK_A       (0)             /* Ax=b    */
#define UMFPACK_DROPTOL 18                /* drop tolerance for entries in L,U */
#define UMFPACK_SCALE 16                /* what row scaling to do */
#define UMFPACK_SCALE_NONE 0
//extern "C" {

int umfpack_di_symbolic
(
  int n_row,
  int n_col,
  const int Ap [ ],
  const int Ai [ ],
  const double Ax [ ],
  void **Symbolic,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
) ;

int umfpack_di_numeric
(
  const int Ap [ ],
  const int Ai [ ],
  const double Ax [ ],
  void *Symbolic,
  void **Numeric,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
) ;

void umfpack_di_free_symbolic
(
  void **Symbolic
) ;

int umfpack_di_solve
(
  int sys,
  const int Ap [ ],
  const int Ai [ ],
  const double Ax [ ],
  double X [ ],
  const double B [ ],
  void *Numeric,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
) ;

void umfpack_di_free_numeric
(
  void **Numeric
) ;

int umfpack_di_get_lunz
(
  int *lnz,
  int *unz,
  int *n_row,
  int *n_col,
  int *nz_udiag,
  void *Numeric
) ;

int umfpack_di_get_numeric
(
  int Lp [ ],
  int Lj [ ],
  double Lx [ ],
  int Up [ ],
  int Ui [ ],
  double Ux [ ],
  int P [ ],
  int Q [ ],
  double Dx [ ],
  int *do_recip,
  double Rs [ ],
  void *Numeric
) ;

int umfpack_di_transpose
(
  int n_row,
  int n_col,
  const int Ap [ ],
  const int Ai [ ],
  const double Ax [ ],
  const int P [ ],
  const int Q [ ],
  int Rp [ ],
  int Ri [ ],
  double Rx [ ]
) ;
//extern "C" {
int gmresm(int, int, int, double, double *,
           double *, int *, int *,
           double *, int *, int *, double *, double *,
           double *, int *, int *, FILE *);

void ilutp_(int *, double *, int *, int *, int *, double *,
            double *, int *, double *, int *, int *, int *,
            double *, double *, int *, int *, int *, int *, int *);

void lusol_(int *, double *, double *, double *, int *, int *);
//}
double TestMe();

void calcDeterminantILUGMRES(int *movedParticle, int *size, int *passedArrLength, double u[], int Ap[], int Ai[], double Ax[], int Arp[], int Ari[], double Arx[], double *detRatio_ILU);

void calcDeterminantILUGMRESNew();
//}
