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
    
    
#include <stdio.h>
#include "ILUGMRESInterface.h"
#include <stdlib.h>

// As ILUTP messes the original array so two sets passed.
void calcDeterminantILUGMRES(int *movedParticle, int *size, int *passedArrLength, double u[], int Arp[], int Ari[], double Arx[], int ArpNew[], int AriNew[], double ArxNew[], double *detRatio_ILU)
{
	fprintf(stderr,"check\n");
        // ------------------------------ COMMON CODE FOR GMRES CALL ---------------------
        ///////////////////////////////////////////////////////////////////////////////////
        int n = *size;
        int i;
        int status;

        double b[n];
        double x[n];
        for (i = 0; i < n; i++)
        {
            b[i] = 0;
            x[i] = 0;
        }
	fprintf(stderr,"Moved particel si %i\n",*movedParticle);
	if (*movedParticle>=n){
	  fprintf(stderr,"Above n error");
	  exit(1);
	}
	    
        b[*movedParticle] = 1;
        int arrLength = *passedArrLength;
	fprintf(stderr,"GMRES done\n");
        //------------------------------------ILUTP CODE FOR ILU ------------------------
        ///////////////////////////////////////////////////////////////////////////////////

        // We also need to convert to the fact that here is C and ilutp is FORTRAN
	for (i = 0; i < n+1; i++)
		Arp[i] = Arp[i] + 1;
	fprintf(stderr,"A done %i\n",arrLength);
	for (i = 0; i < arrLength; i++)
	       	Ari[i] = Ari[i] + 1;
	fprintf(stderr,"B done\n");
        // Calling ILUTP from ilut.f
        int lfil = n;
        double droptol = 0.01;
        double permtol = 0.0;
        int mbloc = n;
        int iwk = arrLength + 2*lfil*n + 2;
	fprintf(stderr,"C done %i %i %i %i\n",iwk,arrLength,lfil,n);
	double testMe[1411790];
	int testMe2[1411790];
	fprintf(stderr,"Cp\n");
        double alu[iwk];
        int jlu[iwk];
	fprintf(stderr,"CC done\n");
        int ju[n];
        double wu[n+1];
        double wl[n];
	fprintf(stderr,"CCC done\n");
        int jr[n];
        int jwl[n];
        int jwu[n];
	fprintf(stderr,"D done\n");
        int iperm[2*n];
        int ierr;
	fprintf(stderr,"hi\n");
        ilutp_(&n, Arx, Ari, Arp, &lfil, &droptol, &permtol, &mbloc, alu, jlu, ju, &iwk, wu, wl, jr, jwl, jwu, iperm, &ierr);

        //For data collection purposes, we need the nnz of L+U. See work out in NoteBook III.
        //Coming from FORTRAN but goes back in FORTRAN (in gmres) so NO need to remove 1 back.
        int nnzLU = jlu[n] - 2 + n;
	printf("nnzLU = %d \n", nnzLU);
	fprintf(stderr,"ilutp done\n");
        //------------------------------------GMRES CALL ---------------------------------
        ///////////////////////////////////////////////////////////////////////////////////
        int restart = n;
        int maxItn = 1;
        double tol = 1e-12;
        double res[n];
        double resnorm[maxItn+1];
        int nit;
        int nmv;
        FILE *outFile = NULL;
        outFile = fopen("gmres-results.txt", "w");
        int flag = gmresm(n, restart, maxItn, tol, x, ArxNew, AriNew, ArpNew, alu, jlu, ju, b, res, resnorm, &nit, &nmv, outFile);
        fclose(outFile);
	printf("GMRES itn = %d \n", nmv);

        // Now need to go to lusol_ from ilut.f as in fmgr_drv.c and gmr_drv.c
        lusol_(&n, x, x, alu, jlu, ju);

        //Computing the dot product (easier to do myself then look for methods).
        *detRatio_ILU = 0;
        for (i = 0; i < n; i++)
                *detRatio_ILU = *detRatio_ILU + x[i]*u[i];
}

double TestMe()
{
  return 5;
}
