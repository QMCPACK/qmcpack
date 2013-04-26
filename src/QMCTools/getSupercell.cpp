#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
using namespace std;


template<typename T>
T getDet(T * mat)
{
  return mat[0]*(mat[4]*mat[8] - mat[7]*mat[5]) - mat[1]*(mat[3]*mat[8] - mat[5]*mat[6]) + mat[2]*(mat[3]*mat[7]-mat[4]*mat[6]);
}

void getSupercell(double* prim, int* tile, double* super)
{
  super[0] = tile[0]*prim[0]+tile[1]*prim[3]+tile[2]*prim[6];
  super[1] = tile[0]*prim[1]+tile[1]*prim[4]+tile[2]*prim[7];
  super[2] = tile[0]*prim[2]+tile[1]*prim[5]+tile[2]*prim[8];
  super[3] = tile[3]*prim[0]+tile[4]*prim[3]+tile[5]*prim[6];
  super[4] = tile[3]*prim[1]+tile[4]*prim[4]+tile[5]*prim[7];
  super[5] = tile[3]*prim[2]+tile[4]*prim[5]+tile[5]*prim[8];
  super[6] = tile[6]*prim[0]+tile[7]*prim[3]+tile[8]*prim[6];
  super[7] = tile[6]*prim[1]+tile[7]*prim[4]+tile[8]*prim[7];
  super[8] = tile[6]*prim[2]+tile[7]*prim[5]+tile[8]*prim[8];
}

template<typename T>
T dot(T* a, T* b)
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

template<typename T>
void cross(T* a, T* b, T* c)
{
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = -a[0]*b[2]+a[2]*b[0];
  c[2] = a[0]*b[1]-a[1]*b[0];
}

double SimCellRad(double* mat)
{
  double A[3];
  double B[3];
  double C[3];
  double BxC[3];
  double radius = 5000000000000000.0;
  for (int i = 0; i < 3; i++)
  {
    const int astart = i*3;
    const int bstart = ((i+1)%3)*3;
    const int cstart = ((i+2)%3)*3;
    for (int j = 0; j < 3; j++)
    {
      A[j] = mat[astart+j];
      B[j] = mat[bstart+j];
      C[j] = mat[cstart+j];
    }
    cross(B,C,BxC);
    double val = abs(0.5*dot(A,BxC) / sqrt(dot(BxC,BxC)));
    if (val < radius)
      radius = val;
  }
  return radius;
}

double getScore(double *mat)
{
  double score = 0;
  static const double tol = 0.001;
  // Highest preference for diagonal matrices
  const double abssumoffidag = abs(mat[1]) + abs(mat[2]) + abs(mat[3]) + abs(mat[5]) + abs(mat[6]) + abs(mat[7]);
  if (abssumoffidag < tol)
  {
    //    cout << "Found a diagonal supercell!" << endl;
    score += 50000;
  }
  // Now prefer positive elements that come as early as possible
  for (int i = 0; i < 9; i++)
  {
    if (mat[i] > 0.0)
    {
      score += 10;
      double v = (9.0-static_cast<double>(i))*0.1;
      score += v*v;
    }
  }
  return score;
}

double WigSeitzRad(double *mat)
{
  double rmin = 1000000000000000;
  for (int i = -1; i <= 1; i++)
  {
    for (int j = -1; j <= 1; j++)
    {
      for (int k = -1; k <= 1; k++)
      {
        if ((i != 0) || (j != 0) || (k != 0))
        {
          double d[3];
          d[0] = i*mat[0]+j*mat[3]+k*mat[6];
          d[1] = i*mat[1]+j*mat[4]+k*mat[7];
          d[2] = i*mat[2]+j*mat[5]+k*mat[8];
          double dist = 0.5 * sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
          if (dist < rmin)
            rmin = dist;
        }
      }
    }
  }
  return rmin;
}


void getBestTile(double* primcell, int target, int* tilemat, double& radius, int range = 7)
{
  double largest = 0.0;
  double bestScore = 0.0;
  static const double tol = 0.0000001;
  double detPrim = getDet(primcell);
  //  cout << "detPrim = " << detPrim << endl;
  if (detPrim < 0)
  {
    target *= -1;
  }
  //  target *= -1;
  #pragma omp parallel
  {
    double my_largest = 0.0;
    int my_besttile[9];
    double localBestScore = 0.0;
    int flag = 0;
    #pragma omp for
    for (int i = -range; i <= range; i++)
    {
      int d[9];
      double super[9];
      d[0] = i;
      for (int j = -range; j <= range; j++)
      {
        d[1] = j;
        for (int k = -range; k <= range; k++)
        {
          d[2] = k;
          for (int l = -range; l <= range; l++)
          {
            d[3] = l;
            for (int m = -range; m <= range; m++)
            {
              d[4] = m;
              int denominator = j*l-i*m;
              for (int n = -range; n <= range; n++)
              {
                d[5] = n;
                int fpp = k*l-i*n;
                for (int o = -range; o <= range; o++)
                {
                  d[6] = o;
                  int sp = o*(n*j-k*m);
                  for (int p = -range; p <= range; p++)
                  {
                    d[7] = p;
                    int numpart = p*fpp+sp;
                    if (denominator != 0)
                    {
                      int rem = 5;
                      rem = (numpart - target) % denominator;
                      if (rem == 0)
                      {
                        d[8] = (numpart - target) / denominator;
                        getSupercell(primcell, d, super);
                        double score = getScore(super);
                        double rad = SimCellRad(super);
                        //double rad = WigSeitzRad(super);
                        if (rad > my_largest+tol || (rad > my_largest-tol && score > localBestScore))
                        {
                          my_largest = rad;
                          localBestScore = score;
                          memcpy(my_besttile, d, 9*sizeof(int));
                        }
                      }
                    }
                    else
                      // Handle case where denominator is 0
                    {
                      if (numpart == target)
                      {
                        //			cout << "Got here" << endl;
                        for (int q = -range; q <= range; q++)
                        {
                          d[8] = q;
                          getSupercell(primcell, d, super);
                          double score = getScore(super);
                          double rad = SimCellRad(super);
                          //double rad = WigSeitzRad(super);
                          if (rad > my_largest+tol || (rad > my_largest-tol && score > localBestScore))
                          {
                            my_largest = rad;
                            localBestScore = score;
                            memcpy(my_besttile, d, 9*sizeof(int));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (my_largest > largest+tol || (my_largest > largest-tol && localBestScore > bestScore))
    {
      #pragma omp critical
      {
        if (my_largest > largest+tol || (my_largest > largest-tol && localBestScore > bestScore) )
        {
          largest = my_largest;
          bestScore = localBestScore;
          memcpy(tilemat, my_besttile, 9*sizeof(int));
        }
      }
    }
  }
  radius = largest;
}



int main(int argc, char* argv[])
{
  double prim[9];
  int target;
  char* pend;
  int verbose = 0;
  int maxentry = 4;
  for (int i =1; i < argc; i++)
  {
    if (i <= argc)
    {
      if (!strcmp(argv[i],"--ptvs"))
      {
        for (int j = 0; j < 9; j++)
        {
          prim[j] = strtod(argv[i+j+1], &pend);
        }
        i += 9;
      }
      if (!strcmp(argv[i],"--target"))
      {
        target = strtol(argv[i+1], &pend, 10);
        i++;
      }
      if (!strcmp(argv[i],"--verbose"))
      {
        verbose=1;
      }
      if (!strcmp(argv[i],"--maxentry"))
      {
        maxentry = strtol(argv[i+1], &pend, 10);
        i++;
      }
    }
  }
  /*
  cout << "The input primitive translation vectors are: \n";
  cout << prim[0] << "   " << prim[1] << "   " << prim[2] << endl;
  cout << prim[3] << "   " << prim[4] << "   " << prim[5] << endl;
  cout << prim[6] << "   " << prim[7] << "   " << prim[8] << endl;
  cout << "The target is " << target << endl;
  */
  int besttile[9];
  double super[9];
  double radius;
  getBestTile(prim, target, besttile, radius, maxentry);
  getSupercell(prim, besttile, super);
  if (verbose)
  {
    cout << "Best Tiling Matrix = " << endl;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        cout << besttile[i*3+j] << "   ";
      }
      cout << endl;
    }
    cout << "Best Supercell = " << endl;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        cout << super[i*3+j] << "   ";
      }
      cout << endl;
    }
    cout << "radius = " << radius << endl;
    double score = getScore(super);
    cout << "score = " << score << endl;
    int trialtile[9];
    trialtile[0] = -1;
    trialtile[1] =  3;
    trialtile[2] = -1;
    trialtile[3] = -3;
    trialtile[4] =  1;
    trialtile[5] =  1;
    trialtile[6] =  1;
    trialtile[7] =  1;
    trialtile[8] =  1;
    double ts[9];
    getSupercell(prim, trialtile, ts);
    cout << "Trial Supercell = " << endl;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        cout << ts[i*3+j] << "   ";
      }
      cout << endl;
    }
    cout << "radius = " << SimCellRad(ts) << endl;
    score = getScore(ts);
    cout << "score = " << score << endl;
  }
  else
  {
    cout << radius << "   ";
    for (int i = 0; i < 9; i++)
    {
      cout << besttile[i] << "   ";
    }
    for (int i = 0; i < 9; i++)
    {
      cout << super[i] << "   ";
    }
    cout << endl;
  }
}


/**************************************************************************
 * $RCSfile$   $Author: lshulenburger $
 * $Revision: 5115 $   $Date: 2011-9-15 10:19:08 -0700 (Thu, 15 Sep 2011) $
 * $Id: getSupercell.cpp 5115 2011-9-15 10:19:08 lshulenburger $
*************************************************************************/
