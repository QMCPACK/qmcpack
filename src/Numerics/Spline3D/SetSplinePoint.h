#ifndef GUARD_SETSPLINEPOINT_H
#define GUARD_SETSPLINEPOINT_H

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/Grid3D.h"

class SetSplinePoint
{

public:

  bool ifout;
  int ix,iy,iz;
  double h,k,l;
  double hinv,kinv,linv;
  double u,v,w;

  void set_point(const posvec_t& r,
                 Grid3D* agrid)
  {
    gridvec_t ir = agrid->ptl(r);
    ix = ir[0];
    iy = ir[1];
    iz = ir[2];
    ifout = false;
    if( r[2] < agrid->m_axis[2].m_start || r[0] > agrid->m_axis[0].m_end ||
        r[2] > agrid->m_axis[2].m_end   || r[0] < agrid->m_axis[0].m_start ||
        r[1] < agrid->m_axis[1].m_start || r[1] > agrid->m_axis[1].m_end )
    {
      ifout = true;
      /*
      cout << "OUTSIDE" << endl;
      cout << r[0] << '\t' << agrid->m_axis[0].m_start << '\t'
      << agrid->m_axis[0].m_end << endl;
      cout << r[1] << '\t' << agrid->m_axis[1].m_start << '\t'
      << agrid->m_axis[1].m_end << endl;
      cout << r[2] << '\t' << agrid->m_axis[2].m_start << '\t'
      << agrid->m_axis[2].m_end << endl;
      */
    }
    h = agrid->m_axis[0].h(ix);
    k = agrid->m_axis[1].h(iy);
    l = agrid->m_axis[2].h(iz);
    hinv = 1.0/h;
    kinv = 1.0/k;
    linv = 1.0/l;
    u = (r[0] - agrid->m_axis[0].m_coord[ix])*hinv;
    v = (r[1] - agrid->m_axis[1].m_coord[iy])*kinv;
    w = (r[2] - agrid->m_axis[2].m_coord[iz])*linv;
  }










};
#endif
