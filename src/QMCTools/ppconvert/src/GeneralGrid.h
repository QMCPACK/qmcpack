#ifndef SIMPLE_GRID_H
#define SIMPLE_GRID_H

#include <vector>

using namespace std;

class SimpleGrid
{
private:
  vector<double> grid;
  
public:
  inline int NumPoints()
  { return grid.size();  }

  inline vector<double>& Points()
  { return grid; }

  inline double operator[](int i)
  { return grid[i]; }

  /// Returns the index of the nearest point below r. 
  int ReverseMap(double r)
  {
    int n = grid.size();
    if (r <= grid[0])
      return (0);
    else if (r >= grid[n-1])
      return n-1;
    else {
      int hi = n-1;
      int lo = 0;
      bool done = false;
      while (!done) {
	int i = (hi+lo)>>1;
	if (grid[i] > r)
	  hi = i;
	else
	  lo = i;
	done = (hi-lo)<2;
      }
      if (grid[lo] >= r)
	lo--;
      return (lo);
    }
  }
  
  inline double Start()
  { return grid[0]; }

  inline double End()
  { return grid[grid.size()-1]; }

  void Init (vector<double> &points)
  {
    grid.resize(points.size());
    grid = points;
  }

  /// Useless constructor
  SimpleGrid ()
  { /*  Do nothing */ }
};

#endif
