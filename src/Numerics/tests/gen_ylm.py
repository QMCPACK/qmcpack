
# Generate spherical harmonic values to test against

# Redirect output to ylm.inc ("python gen_ylm.py > ylm.inc") and adjust
# the #ifdef in test_ylm.cpp

from sympy import mpmath
import math

def gen_spherical_harmonics():

  # Set up a grid of points in (theta, phi) to compare against
  n_theta = 5
  a_theta = 0.0
  b_theta = math.pi
  h_theta = (b_theta - a_theta)/n_theta

  n_phi = 5
  a_phi = 0.0
  b_phi = 2*math.pi
  h_phi = (b_phi - a_phi)/n_phi

  vals = []
  for l in range(0,5):
    for m in range(-l,l+1):
      for i in range (n_theta+1):
        theta = a_theta + h_theta*i
        for j in range (n_phi+1):
          phi = a_phi + h_phi*j
          x = math.sin(theta)*math.cos(phi)
          y = math.sin(theta)*math.sin(phi)
          z = math.cos(theta)
          vec = (x,y,z)
          y = mpmath.spherharm(l,m,theta,phi)
          #print theta,vec,y
          vals.append( (vec,l,m,y) )

 # Assumes
 # struct Point
 #  {
 #    double x;
 #    double y;
 #    double z;
 #  };
 # struct YlmValue
 # {
 #   Point p;
 #   int l;
 #   int m;
 #   double y_re;
 #   double y_im;
 # };
 #   
  tmpl = """
  const int N = {N};
  YlmValue Vals[N] = {{
    {values}
  }};
"""
  fmt_values = ',\n  '.join(["{ {%g, %g, %g}, %d, %d, %g, %g}"%(p[0],p[1],p[2],l,m,y.real,y.imag) for (p,l,m,y) in vals])
  s = tmpl.format(N=len(vals), values=fmt_values)
  print s

if __name__ == '__main__':
  gen_spherical_harmonics()
