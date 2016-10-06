//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef GRID_H
#define GRID_H

#include "IO.h"
using namespace IO;
#include "Blitz.h"

//Ken's Grid Class

/// The different types of grids that we currently allow
typedef enum {NONE, LINEAR, OPTIMAL, OPTIMAL2, LOG, 
	      CLUSTER, GENERAL, CENTER} GridType;


/// Parent class for all grids
class Grid
{
 protected:
  /// Contains the grid points 
  Array<double,1> grid;
 public:
  /// First and last grid points
  double Start, End;

  /// Number of points in the grid
  int NumPoints;

  /// The i'th point in the grid
  inline double operator()(int i) const
  {
    return (grid(i));
  }
  inline double* data()
  {
    return grid.data();
  }
  inline Array<double,1>& Points()
  {
    return grid;
  }

  /// Returns the type of the grid (i.e. linear, optimal, etc)
  virtual GridType Type() = 0;
   
  ///Returns the index of the nearest point below r. 
  virtual int ReverseMap (double r) = 0;
  virtual void Write (IOSectionClass &out) = 0;
  virtual void Read  (IOSectionClass &inSection) = 0;
};


/// Linear Grid inherets from Grid.  
class LinearGrid : public Grid
{
 private:
  /// The value between successive grid points.
  double delta, deltainv;
  inline void CheckRoundingMode();
 public:
  /// Returns the type of the grid (in this case LINEAR)
  GridType Type()
  { return (LINEAR); }  

  /// Returns the index of the nearest point below r. 
  int ReverseMap(double r)
  {
    return ((int)nearbyint((r-Start)*deltainv-0.5));
  }

  /// Initializes the linear grid.
  inline void Init(double start, double end, int numpoints)
  {
    Start=start; End=end; NumPoints=numpoints;
    grid.resize(NumPoints);
    delta = (End-Start)/(double)(NumPoints-1);
    deltainv = 1.0/delta;
    for (int i=0; i<NumPoints; i++)
      grid(i) = Start + (double)i*delta;
    CheckRoundingMode();
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", std::string("Linear"));
    outSection.WriteVar ("Start", Start);
    outSection.WriteVar ("End", End);
    outSection.WriteVar ("NumPoints", NumPoints);
  }

  void Read (IOSectionClass &inSection)
  {
    assert(inSection.ReadVar("Start", Start));
    assert(inSection.ReadVar("End", End));
    assert(inSection.ReadVar("NumPoints", NumPoints));
    Init (Start, End, NumPoints);
  }

  /// Useless constructor
  LinearGrid ()
  { /*  Do nothing */ }

  LinearGrid& operator= (const LinearGrid& lin)
  {
    grid.resize(lin.grid.shape());
    Start    = lin.Start;
    End      = lin.End;
    grid     = lin.grid;
    delta    = lin.delta;
    deltainv = lin.deltainv;
    return *this;
  }

  /// Constructor that sets the number of points, start and end point
  /// of the original grid 
  LinearGrid (double start, double end, int numpoints)
  {
    Init (start, end, numpoints);
  }
};


/// General Grid inherets from Grid.  
class GeneralGrid : public Grid
{
 public:
  /// Returns the type of the grid (in this case GENERAL)
  GridType Type()
  { return (GENERAL); }  

  /// Returns the index of the nearest point below r. 
  int ReverseMap(double r)
  {
    if (r <= grid(0))
      return (0);
    else if (r >= grid(NumPoints-1))
      return (NumPoints-1);
    else {
      int hi = NumPoints-1;
      int lo = 0;
      bool done = false;
      while (!done) {
	int i = (hi+lo)>>1;
	if (grid(i) > r)
	  hi = i;
	else
	  lo = i;
	done = (hi-lo)<2;
      }
      return min (lo, NumPoints-2);
    }
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", std::string("General"));
  }

  void Read (IOSectionClass &inSection)
  {
    assert (inSection.ReadVar("Points", grid));
    Start = grid(0);
    End = grid(grid.size()-1);
    NumPoints = grid.size();
  }

  void Init (Array<double,1> &points)
  {
    NumPoints = points.size();
    grid.resize(NumPoints);
    grid = points;
    Start = grid(0);
    End = grid(NumPoints-1);
  }

  /// Useless constructor
  GeneralGrid ()
  { /*  Do nothing */ }

};



/// The OptimalGrid class stores a grid which has linear spacing at
/// the origin and exponential spacing further out.  It has the
/// analytic form \f[r_k = a\left(e^{kb}-1\right)\f].
class OptimalGrid : public Grid
{
 private:
  double a, b;
  
 public:

  GridType Type()
  { return (OPTIMAL); }

  int ReverseMap(double r)
  {
    if ((r/a) < 1e-6)
      return ((int)floor(r/(a*b)+0.5)-1);
    else
      return((int)floor(log(r/a + 1.0)/b + 0.5) -1);
  }
  
  /// Returns a parameter
  double Geta() const
  { return (a); }
  
  /// Returns b parameter
  double Getb() const
  { return (b); }
  
  OptimalGrid ()
  {
    // Do nothing
  }
  
  /// This form of the constructor takes the number of points, the
  /// maximum radius and the value of b.
  void Init (int numpoints, double rmax, double bval)
  {
    NumPoints = numpoints;
    b = bval;
    End = rmax;
    a = End / (exp(b*(double)NumPoints) - 1.0);  
    Start = a * (exp(b) - 1.0);
    grid.resize(NumPoints);
    
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*(i+1))-1.0);
  }
  
  OptimalGrid (int numPoints, double rmax, double bval)
  { Init (numPoints, rmax, bval); }
  
  void Init (double aval, double bval, int numPoints)
  {
    a = aval;
    b = bval;
    NumPoints = numPoints;
    Start = a * (exp(b) - 1.0);
    End   = a * (exp(b*NumPoints) - 1.0);
    
    grid.resize(NumPoints);
    
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*(i+1))-1.0);
  }
  
  
  /// This form of the constructor takes a, b, and the number of points.
  OptimalGrid (double aval, double bval, int numPoints)
  { 
    Init (aval, bval, numPoints);
  }
  
  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", std::string("Optimal"));
    outSection.WriteVar ("a", a);
    outSection.WriteVar ("b", b);
    outSection.WriteVar ("NumPoints", NumPoints);
  }
  
  void Read (IOSectionClass &inSection)
  {
    double aval, bval;
    int numPoints;
    if (inSection.ReadVar("a", aval)) {
      assert(inSection.ReadVar("b", bval));
      assert(inSection.ReadVar("NumPoints", numPoints));
      Init (aval,bval,numPoints);
    }
    else {
      double Z, rmax;
      assert(inSection.ReadVar("Z", Z));
      assert(inSection.ReadVar("rmax", rmax));
      Init (Z, rmax);
    }
  }
  
  void InitRatio (double end, double ratio, int numpoints)
  {
    End = end; 
    NumPoints = numpoints;
    
    b = log(ratio)/(double)(numpoints-2);
    a = end/(exp(b*(double)(numpoints-1)) - 1);
      
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*i)-1.0);
  }



  /// This form of the constructor takes a nuclear charge and a
  /// maxmimum radius and chooses an appropriate number of points for
  /// that atom.
  void Init (double Z, double rmax)
  {
    a = 4.34e-6/Z;
    //a = 4.0e-2;
    b = 0.002304;
    //b = 0.004;
    
    NumPoints = (int)ceil(log(rmax/a+1.0)/b);
    b = log(rmax/a+1.0)/(double)NumPoints;
    Start = a * (exp(b) - 1.0);
    End = rmax;
    //End   = a * (exp(b*NumPoints) - 1.0);
    
    grid.resize(NumPoints);
    
    for (int i=0; i<NumPoints; i++) {
      grid(i) = a*(exp(b*(i+1))-1.0);
      //fprintf (stdout, "%1.12e\n", grid(i));
    }
  }

  inline OptimalGrid& operator= (const OptimalGrid& opt)
  {
    grid.resize(opt.grid.shape());
    a         = opt.a;
    b         = opt.b;
    NumPoints = opt.NumPoints;
    Start     = opt.Start;
    End       = opt.End;
    grid      = opt.grid;
    return *this;
  }

  OptimalGrid (double Z, double rmax)
  { Init (Z, rmax); }

};


/// The OptimalGrid class stores a grid which has linear spacing at
/// the origin and exponential spacing further out.  It has the
/// analytic form \f[r_k = a\left(e^{kb}-1\right)\f].
class OptimalGrid2 : public Grid
{
 private:
  double a, b, c;
  double Ratio;
 public:

  GridType Type()
  { return (OPTIMAL2); }

  int ReverseMap(double r)
  {
//     if ((r/a) < 1e-6)
//       return ((int)floor(r/(a*b)));
//     else
    return((int)floor(log1p((r-c)/a)/b));
  }
  
  /// Returns a parameter
  double Geta() const
  { return (a); }
  
  /// Returns b parameter
  double Getb() const
  { return (b); }
  
  OptimalGrid2 ()
  {
    // Do nothing
  }

  /// This form of the constructor takes the number of points, the
  /// maximum radius and the value of b.
  OptimalGrid2 (int numpoints, double rmax, double bval)
  {
    NumPoints = numpoints;
    b = bval;
    End = rmax;
    a = End / (exp(b*(double)NumPoints) - 1.0);  
    Start = a * (exp(b) - 1.0);
    grid.resize(NumPoints);
    c = 0.0;
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = c + a*expm1(b*i);
  }

  void Init (double start, double end, double ratio, int numpoints)
  {
    Start = start;
    End = end; 
    Ratio = ratio;
    NumPoints = numpoints;
    
    b = log(ratio)/(double)(numpoints-2);
    c = Start;
    a = (end - c)/expm1(b*(double)(numpoints-1));
      
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = c + a*expm1(b*i);
  }


  /// This form of the constructor takes a, b, and the number of points.
  OptimalGrid2 (double start, double end, double ratio, int numpoints)
  { 
    Init (start, end, ratio, numpoints);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", std::string("Optimal2"));
    outSection.WriteVar ("Start", Start);
    outSection.WriteVar ("End", End);
    outSection.WriteVar ("Ratio", Ratio);
    outSection.WriteVar ("NumPoints", NumPoints);
  }

  void Read (IOSectionClass &inSection)
  {
    double start, end, ratio;
    int numPoints;
    assert(inSection.ReadVar("Start", start));
    assert(inSection.ReadVar("End", end));
    assert(inSection.ReadVar("Ratio", ratio));
    assert(inSection.ReadVar("NumPoints", numPoints));
    Init (start,end,ratio,numPoints);
  }
};


class CenterGrid : public Grid
{
private:
  double a, aInv, b, bInv, center;
  int HalfPoints;
  bool Odd;
  double EvenHalf;
  int OddOne;
public:
  // ratio gives approximately the largest grid spacing divided by the
  // smallest. 
  GridType Type()
  { return CENTER; }

  int ReverseMap (double x)
  {
    x -= center;
    double index = copysign (log1p(std::abs(x)*aInv)*bInv, x);
    return (int)floor(HalfPoints + index - EvenHalf);
  }
  void Write (IOSectionClass &out) {
  }
  void Read  (IOSectionClass &in) {
  }
  void Init (double start, double end, double ratio, int numPoints) {
    assert (ratio > 1.0);
    Start      = start;
    End        = end;
    center     = 0.5*(start + end);
    NumPoints  = numPoints;
    HalfPoints = numPoints/2;
    Odd = ((numPoints % 2) == 1);
    b = log(ratio) / (double)(HalfPoints-1);
    bInv = 1.0/b;
    grid.resize(numPoints);
    if (Odd) {
      EvenHalf = 0.0;  OddOne = 1;
      a = 0.5*(end-start)/expm1(b*HalfPoints);
      aInv = 1.0/a;
      for (int i=-HalfPoints; i<=HalfPoints; i++) {
	double sign = (i<0) ? -1.0 : 1.0;
	grid(i+HalfPoints) = sign * a*expm1(b*std::abs(i))+center;
      }
    }
    else {
      EvenHalf = 0.5;  OddOne = 0;
      a = 0.5*(end-start)/expm1(b*(-0.5+HalfPoints));
      aInv = 1.0/a;
      for (int i=-HalfPoints; i<HalfPoints; i++) {
	double sign = (i<0) ? -1.0 : 1.0;
	grid(i+HalfPoints) = sign * a*expm1(b*std::abs(0.5+i)) + center;
      }
    }
  }

};



/// LogGrid is a function whose gridpoints increase exponentially with
/// the index.  That is, it has the analytic form
/// \f[ r_k = \frac{r_0}{Z} \Delta^k.\f]  It is appropriate for
/// functions which change rapidly near the origin but vary smoothly
/// further out.
class LogGrid : public Grid
{
 public:
  double Z, r0, Spacing;

  GridType Type()
    { return (LOG); }

  int ReverseMap(double r)
  {
    return ((int)(floor(log(Z*r/r0)/log(Spacing))));
  }
  
  LogGrid ()
  {
    // Do nothing
  }
  
  void Init (double R0, double spacing, int numpoints)
  {
    NumPoints = numpoints;
    Z = 1.0; r0 = R0; Spacing = spacing;
    Start = r0;
    End = r0 * pow(Spacing, (double) NumPoints-1);
    grid.resize (NumPoints);
    
    for (int i=0; i<NumPoints; i++)
      grid(i) = r0 * pow(Spacing, (double) i);
  }



  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", std::string("Log"));
    outSection.WriteVar ("r0", r0);
    outSection.WriteVar ("Spacing", Spacing);
  }

  void Read (IOSectionClass &inSection)
  {
    double tempr0, tempSpacing;
    int  tempNumPoints;
    assert (inSection.ReadVar("r0", tempr0));
    assert (inSection.ReadVar("Spacing", tempSpacing));
    assert (inSection.ReadVar("NumPoints", tempNumPoints));
    Init (tempr0, tempSpacing, tempNumPoints);
  }

  LogGrid (double R0, double spacing, int numpoints)
  {
    Init  (R0, spacing, numpoints);
  }

  LogGrid (int numpoints, double z, double R0, double spacing)
  {
    NumPoints = numpoints;
    Z = z; r0 = R0; Spacing = spacing;
    
    Start = r0 / Z;
    End = r0/Z * pow(Spacing, (double) (NumPoints-1));
    
    grid.resize(NumPoints);
    
    for (int i=0; i<NumPoints; i++)
      grid(i) = r0/Z * pow (Spacing, (double) i);
  }
};




/// ClusterGrid is a function whose gridpoints are clustered tightly
/// around the origin.
class ClusterGrid : public Grid
{
private:
  double x0, dri, rr;

 public:
  double Start, End, Cluster;

  GridType Type()
    { return (CLUSTER); }

  int ReverseMap(double r)
    {
      return ((int)floor (dri/(r-rr) -1.0 + x0));
    }

  void Init (double start, double end, double cluster, int numpoints)
  {
    Start = start; End = end; Cluster = cluster;
    NumPoints = numpoints;
    
    x0 = (NumPoints - Cluster)/(1.0-Cluster);
    dri = -(End-Start)*((double) NumPoints-x0)*(1.0-x0) / 
      ((double)NumPoints-1.0);
    rr = Start - dri/(1.0-x0);
    grid.resize(NumPoints);
    for (int i=0; i<NumPoints; i++)
      grid(i) = rr+dri/(i+1.0-x0);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", std::string("Cluster"));
    outSection.WriteVar ("Start", Start);
    outSection.WriteVar ("End", End);
    outSection.WriteVar ("Cluster", Cluster);
    outSection.WriteVar ("NumPoints", NumPoints);
  }

  void Read (IOSectionClass &inSection)
  {
    double start, end, cluster;
    int numpoints;
    assert (inSection.ReadVar("Start", start));
    assert (inSection.ReadVar("End", end));
    assert (inSection.ReadVar("Cluster", cluster));
    assert (inSection.ReadVar("NumPoints", numpoints));
    Init (start, end, cluster, numpoints);
  }

  ClusterGrid (double start, double end, double cluster, int numpoints)
  {
    Init  (start, end, cluster, numpoints);
  }
  ClusterGrid () 
  { /* Do nothing */ }
};


inline Grid* ReadGrid (IOSectionClass &inSection)
{
  std::string Type;
  assert (inSection.ReadVar ("Type", Type));

  Grid *newGrid;
  if (Type == "Linear")
    newGrid = new LinearGrid;
  else if (Type == "General")
    newGrid = new GeneralGrid;
  else if (Type == "Optimal")
    newGrid = new OptimalGrid;
  else if (Type == "Optimal2")
    newGrid = new OptimalGrid2;
  else if (Type == "Log")
    newGrid = new LogGrid;
  else if (Type == "Cluster")
    newGrid = new ClusterGrid;
  else
    {
      std::cerr << "Unrecognized Grid type " << Type << "\n";
      exit(1);
    }
  newGrid->Read(inSection);
  return (newGrid);
}


inline void
LinearGrid::CheckRoundingMode()
{
  for (int i=0; i<100; i++) {
    double x = 100.0*drand48()-50.0;
    if (nearbyint(x) != round(x)) {
      std::cerr << "Error in rounding mode detected in LinearGrid.  Abort.\n";
      abort();
    }
  }
}
      




#endif
