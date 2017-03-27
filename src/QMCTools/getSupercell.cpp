//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <cstring>


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
    double val = std::abs(0.5*dot(A,BxC) / sqrt(dot(BxC,BxC)));
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
  const double abssumoffidag = std::abs(mat[1]) + std::abs(mat[2]) + std::abs(mat[3]) + std::abs(mat[5]) + std::abs(mat[6]) + std::abs(mat[7]);
  if (abssumoffidag < tol)
  {
    //    std::cout << "Found a diagonal supercell!" << std::endl;
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
  //  std::cout << "detPrim = " << detPrim << std::endl;
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
                          std::memcpy(my_besttile, d, 9*sizeof(int));
                        }
                      }
                    }
                    else
                      // Handle case where denominator is 0
                    {
                      if (numpart == target)
                      {
                        //			cout << "Got here" << std::endl;
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
                            std::memcpy(my_besttile, d, 9*sizeof(int));
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
          std::memcpy(tilemat, my_besttile, 9*sizeof(int));
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
      if (!std::strcmp(argv[i],"--ptvs"))
      {
        for (int j = 0; j < 9; j++)
        {
          prim[j] = strtod(argv[i+j+1], &pend);
        }
        i += 9;
      }
      if (!std::strcmp(argv[i],"--target"))
      {
        target = strtol(argv[i+1], &pend, 10);
        i++;
      }
      if (!std::strcmp(argv[i],"--verbose"))
      {
        verbose=1;
      }
      if (!std::strcmp(argv[i],"--maxentry"))
      {
        maxentry = strtol(argv[i+1], &pend, 10);
        i++;
      }
    }
  }
  /*
  std::cout << "The input primitive translation vectors are: \n";
  std::cout << prim[0] << "   " << prim[1] << "   " << prim[2] << std::endl;
  std::cout << prim[3] << "   " << prim[4] << "   " << prim[5] << std::endl;
  std::cout << prim[6] << "   " << prim[7] << "   " << prim[8] << std::endl;
  std::cout << "The target is " << target << std::endl;
  */
  int besttile[9];
  double super[9];
  double radius;
  getBestTile(prim, target, besttile, radius, maxentry);
  getSupercell(prim, besttile, super);
  if (verbose)
  {
    std::cout << "Best Tiling Matrix = " << std::endl;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        std::cout << besttile[i*3+j] << "   ";
      }
      std::cout << std::endl;
    }
    std::cout << "Best Supercell = " << std::endl;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        std::cout << super[i*3+j] << "   ";
      }
      std::cout << std::endl;
    }
    std::cout << "radius = " << radius << std::endl;
    double score = getScore(super);
    std::cout << "score = " << score << std::endl;
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
    std::cout << "Trial Supercell = " << std::endl;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        std::cout << ts[i*3+j] << "   ";
      }
      std::cout << std::endl;
    }
    std::cout << "radius = " << SimCellRad(ts) << std::endl;
    score = getScore(ts);
    std::cout << "score = " << score << std::endl;
  }
  else
  {
    std::cout << radius << "   ";
    for (int i = 0; i < 9; i++)
    {
      std::cout << besttile[i] << "   ";
    }
    for (int i = 0; i < 9; i++)
    {
      std::cout << super[i] << "   ";
    }
    std::cout << std::endl;
  }
}


