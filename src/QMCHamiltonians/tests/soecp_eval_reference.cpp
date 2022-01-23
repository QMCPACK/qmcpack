#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

double PI   = 4. * std::atan(1);
using dcomp = std::complex<double>;
using vec   = std::vector<double>;

//Same quadrature grid as in QMCPACK
std::vector<vec> quad = {{1, 0, 0},
                         {-1, 1.224646853e-16, 0},
                         {0.4472135901, 0.8944271803, 0},
                         {-0.4472135901, 0.7236068249, 0.5257310867},
                         {0.4472135901, 0.2763932049, 0.8506507874},
                         {-0.4472135901, -0.2763932049, 0.8506507874},
                         {0.4472135901, -0.7236068249, 0.5257310867},
                         {-0.4472135901, -0.8944271803, 1.095357398e-16},
                         {0.4472135901, -0.7236068249, -0.5257310867},
                         {-0.4472135901, -0.2763932049, -0.8506507874},
                         {0.4472135901, 0.2763932049, -0.8506507874},
                         {-0.4472135901, 0.7236068249, -0.5257310867}};

vec wt = {0.08333333582, 0.08333333582, 0.08333333582, 0.08333333582, 0.08333333582, 0.08333333582,
          0.08333333582, 0.08333333582, 0.08333333582, 0.08333333582, 0.08333333582, 0.08333333582};

dcomp I = dcomp(0, 1);


dcomp Ylm(int l, int m, const vec& sph)
{
  //From wiki
  double pref;
  switch (l)
  {
  case 0:
    return dcomp(0.5 * std::sqrt(1. / PI), 0.0);
    break;
  case 1:
    switch (m)
    {
    case -1:
      pref = 0.5 * std::sqrt(3. / (2 * PI));
      return pref * std::exp(-I * sph[2]) * std::sin(sph[1]);
      break;
    case 0:
      pref = 0.5 * std::sqrt(3. / (PI));
      return pref * std::cos(sph[1]);
      break;
    case 1:
      pref = -0.5 * std::sqrt(3. / (2 * PI));
      return pref * std::exp(I * sph[2]) * std::sin(sph[1]);
      break;
    default:
      exit(1);
      break;
    }
    break;
  case 2:
    switch (m)
    {
    case -2:
      pref = 0.25 * std::sqrt(15. / (2 * PI));
      return pref * std::exp(-2.0 * I * sph[2]) * std::sin(sph[1]) * std::sin(sph[1]);
      break;
    case -1:
      pref = 0.5 * std::sqrt(15. / (2 * PI));
      return pref * std::exp(-I * sph[2]) * std::sin(sph[1]) * std::cos(sph[1]);
      break;
    case 0:
      pref = 0.25 * std::sqrt(5. / (PI));
      return pref * (3. * std::cos(sph[1]) * std::cos(sph[1]) - 1.);
      break;

    case 1:
      pref = -0.5 * std::sqrt(15. / (2 * PI));
      return pref * std::exp(I * sph[2]) * std::sin(sph[1]) * std::cos(sph[1]);
      break;
    case 2:
      pref = 0.25 * std::sqrt(15. / (2 * PI));
      return pref * std::exp(2.0 * I * sph[2]) * std::sin(sph[1]) * std::sin(sph[1]);
      break;
    }
  case 3:
    switch (m)
    {
    case -3:
      pref = 0.125 * std::sqrt(35. / PI);
      return pref * std::exp(-3. * I * sph[2]) * std::pow(std::sin(sph[1]), 3);
      break;
    case -2:
      pref = 0.25 * std::sqrt(105. / (2 * PI));
      return pref * std::exp(-2. * I * sph[2]) * std::pow(std::sin(sph[1]), 2) * std::cos(sph[1]);
      break;
    case -1:
      pref = 0.125 * std::sqrt(21. / PI);
      return pref * std::exp(-I * sph[2]) * std::sin(sph[1]) * (5 * std::pow(std::cos(sph[1]), 2) - 1);
      break;
    case 0:
      pref = 0.25 * std::sqrt(7. / PI);
      return pref * (5 * std::pow(std::cos(sph[1]), 3) - 3. * std::cos(sph[1]));
      break;
    case 1:
      pref = -0.125 * std::sqrt(21. / PI);
      return pref * std::exp(I * sph[2]) * std::sin(sph[1]) * (5 * std::pow(std::cos(sph[1]), 2) - 1);
      break;
    case 2:
      pref = 0.25 * std::sqrt(105. / (2 * PI));
      return pref * std::exp(2. * I * sph[2]) * std::pow(std::sin(sph[1]), 2) * std::cos(sph[1]);
      break;
    case 3:
      pref = -0.125 * std::sqrt(35. / PI);
      return pref * std::exp(3. * I * sph[2]) * std::pow(std::sin(sph[1]), 3);
      break;
    }
    break;
  default:
    exit(1);
    break;
  }
}


vec cart2sph(const vec& cart)
{
  //Physics convention
  vec sph(3);
  sph[0] = std::sqrt(cart[0] * cart[0] + cart[1] * cart[1] + cart[2] * cart[2]);
  sph[1] = std::acos(cart[2] / sph[0]);
  sph[2] = std::atan2(cart[1], cart[0]);
  return sph;
}

vec sph2cart(const vec& sph)
{
  //Physics convention
  vec cart(3);
  cart[0] = sph[0] * std::sin(sph[1]) * std::cos(sph[2]);
  cart[1] = sph[0] * std::sin(sph[1]) * std::sin(sph[2]);
  cart[2] = sph[0] * std::cos(sph[1]);
  return cart;
}

//spin representation
dcomp chiu(double s) { return std::exp(I * s); }
dcomp chid(double s) { return std::exp(-I * s); }

//simple plane-wave
dcomp orbu(const vec& sph)
{
  vec cart = sph2cart(sph);
  return std::exp(I * (cart[0] + cart[1] + cart[2]));
}
dcomp orbd(const vec& sph)
{
  vec cart = sph2cart(sph);
  return std::exp(2. * I * (cart[0] + cart[1] + cart[2]));
}

//single particle spinor
dcomp spinor(const vec& sph, double s) { return orbu(sph) * chiu(s) + orbd(sph) * chid(s); }

// <s1 | s_d | s2>
dcomp sMatrixElement(double s1, double s2, int d)
{
  switch (d)
  {
  case 0:
    return dcomp(std::cos(s1 + s2), 0.0);
    break;
  case 1:
    return dcomp(std::sin(s1 + s2), 0.0);
    break;
  case 2:
    return dcomp(0.0, std::sin(s1 - s2));
    break;
  default:
    exit(1);
    break;
  }
}

int kroneckerDelta(int x, int y) { return (x == y) ? 1 : 0; }

//<l m1 | l_d | l m2>
dcomp lMatrixElement(int l, int m1, int m2, int d)
{
  double pref1, pref2, val;
  switch (d)
  {
  case 0:
    pref1 = std::sqrt(l * (l + 1) - m2 * (m2 + 1));
    pref2 = std::sqrt(l * (l + 1) - m2 * (m2 - 1));
    val   = 0.5 * (pref1 * kroneckerDelta(m1, m2 + 1) + pref2 * kroneckerDelta(m1, m2 - 1));
    return dcomp(val, 0.0);
    break;
  case 1:
    pref1 = std::sqrt(l * (l + 1) - m2 * (m2 - 1));
    pref2 = std::sqrt(l * (l + 1) - m2 * (m2 + 1));
    val   = 0.5 * (pref1 * kroneckerDelta(m1, m2 - 1) - pref2 * kroneckerDelta(m1, m2 + 1));
    return dcomp(0.0, val);
    break;
  case 2:
    return dcomp(m2 * kroneckerDelta(m1, m2), 0.0);
    break;
  default:
    exit(1);
    break;
  }
}

//simple radial dependence. coded in so_ecp_test.xml
//used in so_ecp_test.xml for each spin channel
double Wso(int l, double r) { return exp(-l * r * r); }

dcomp calcAngInt(const vec& sph1, double s1, double s2)
{
  dcomp angint(0.0, 0.0);
  for (int i = 0; i < quad.size(); i++)
  {
    vec sph2 = cart2sph(quad[i]);
    sph2[0] *= sph1[0]; //now scaled to appropriate distance

    dcomp integrand(0.0, 0.0);
    for (int l = 1; l <= 3; l++)
    {
      dcomp msum(0.0, 0.0);
      for (int m1 = -l; m1 <= l; m1++)
      {
        for (int m2 = -l; m2 <= l; m2++)
        {
          dcomp ldots(0.0, 0.0);
          for (int d = 0; d < 3; d++)
            ldots += lMatrixElement(l, m1, m2, d) * sMatrixElement(s1, s2, d);
          msum += Ylm(l, m1, sph1) * std::conj(Ylm(l, m2, sph2)) * ldots;
        }
      }
      integrand += Wso(l, sph1[0]) * msum;
    }
    integrand *= spinor(sph2, s2) / spinor(sph1, s1);
    angint += integrand * wt[i] * 4.0 * PI;
  }
  return angint;
}

void calcVal(int npts)
{
  vec cart1 = {0.138, -0.24, 0.216};
  vec sph1  = cart2sph(cart1);
  double s1 = 0.0;

  dcomp sint(0.0, 0.0);
  double smin = 0.0;
  double smax = 2 * PI;
  double h    = (smax - smin) / npts;
  for (int k = 1; k <= npts - 1; k += 2)
  {
    double s2    = smin + k * h;
    dcomp angint = calcAngInt(sph1, s1, s2);
    sint += 4 * h / 3. * angint;
  }
  for (int k = 2; k <= npts - 2; k += 2)
  {
    double s2    = smin + k * h;
    dcomp angint = calcAngInt(sph1, s1, s2);
    sint += 2 * h / 3. * angint;
  }
  sint += h / 3. * calcAngInt(sph1, s1, smin);
  sint += h / 3. * calcAngInt(sph1, s1, smax);
  sint /= (2.0 * PI);
  std::cout << npts << " " << std::setprecision(10) << std::real(sint) << " " << std::imag(sint) << std::endl;
}

int main()
{
  for (int n = 2; n <= 100; n += 2)
    calcVal(n);
}
