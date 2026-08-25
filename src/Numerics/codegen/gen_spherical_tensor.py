#!/usr/bin/env python3
# Generate SoaSphericalTensor using the solid-harmonic derivative recurrence.

import subprocess
from functools import lru_cache
from pathlib import Path

from sympy import (Integer, Poly, Symbol, cos, diff, expand, expand_complex,
                   expand_func, expand_trig, im, powdenest, re, simplify, sin,
                   sqrt, sympify)
from sympy.functions.special.spherical_harmonics import Ynm
from sympy.printing.cxx import cxxcode


X, Y, Z = Symbol("x", real=True), Symbol("y", real=True), Symbol("z", real=True)
R, RHO = Symbol("r", positive=True), Symbol("rho", positive=True)
THETA, PHI = Symbol("theta", real=True), Symbol("phi", real=True)


def index(l, m):
    """Return the QMCPACK storage index for a real solid harmonic."""
    return l * (l + 1) + m


@lru_cache(maxsize=None)
def create_solid_harmonic_symbolic(l, m):
    """Construct r**l S_l^m in the Gaussian real-harmonic convention.

    This follows the addsign=true convention documented by
    SoaSphericalTensor. It is constructed directly from SymPy's complex
    spherical harmonic rather than from QMCPACK's runtime recurrence.
    """
    if m == 0:
        angular = expand_func(Ynm(l, 0, THETA, PHI))
    else:
        complex_ylm = expand_complex(expand_func(Ynm(l, abs(m), THETA, PHI)))
        real_part = re(complex_ylm) if m > 0 else im(complex_ylm)
        angular = Integer(-1) ** m * sqrt(2) * real_part

    angular = expand_trig(angular)
    regular = expand(angular * R**l)
    regular = regular.subs(
        {
            cos(PHI): X / RHO,
            sin(PHI): Y / RHO,
            cos(THETA): Z / R,
            sin(THETA): RHO / R,
        },
        simultaneous=True,
    )
    regular = expand(regular)
    regular = regular.subs(RHO**2, X**2 + Y**2).subs(R**2, X**2 + Y**2 + Z**2)
    regular = powdenest(regular, force=True)
    regular = simplify(
        regular.subs(RHO, sqrt(X**2 + Y**2)).subs(R, sqrt(X**2 + Y**2 + Z**2))
    )
    regular = expand(regular)

    # A regular solid harmonic must be a Cartesian polynomial.
    Poly(regular, X, Y, Z)
    return regular


@lru_cache(maxsize=None)
def create_raw_solid_harmonic_symbolic(l, m):
    """Construct the solid harmonic before norm_factor is applied."""
    solid_harmonic = create_solid_harmonic_symbolic(l, m)
    if m:
        solid_harmonic /= Integer(-1) ** m * sqrt(2)
    return expand(solid_harmonic)


def apply_gradient_recurrence(l, m, lower_value):
    """Apply the runtime gradient recurrence to lower-l symbolic values."""
    ma = abs(m)
    fac = Integer(2 * l + 1) / Integer(2 * l - 1)
    cp = sqrt(fac * (l - ma - 1) * (l - ma))
    cm = sqrt(fac * (l + ma - 1) * (l + ma))
    c0 = sqrt(fac * (l - ma) * (l + ma))

    dz = c0 * lower_value(m) if l > ma else Integer(0)

    if l > ma + 1:
        dpr = cp * lower_value(ma + 1)
        dpi = cp * lower_value(-ma - 1)
    else:
        dpr = Integer(0)
        dpi = Integer(0)

    if l > 1:
        if ma == 0:
            dmr = -cm * lower_value(1)
            dmi = cm * lower_value(-1)
        elif ma == 1:
            dmr = cm * lower_value(0)
            dmi = Integer(0)
        else:
            dmr = cm * lower_value(ma - 1)
            dmi = cm * lower_value(-ma + 1)
    else:
        dmr = cm * lower_value(0)
        dmi = Integer(0)

    half = Integer(1) / 2
    if m < 0:
        dx = half * (dpi - dmi)
        dy = -half * (dpr + dmr)
    else:
        dx = half * (dpr - dmr)
        dy = half * (dpi + dmi)

    return tuple(simplify(value) for value in (dx, dy, dz))


def verify_hessian_recurrence(lmax=6):
    """Verify the compact Hessian recurrence against symbolic harmonics."""
    coordinates = (X, Y, Z)
    for l in range(1, lmax + 1):
        for m in range(-l, l + 1):
            target = create_raw_solid_harmonic_symbolic(l, m)
            for source_coordinate in coordinates:
                lower_value = lambda lower_m: diff(
                    create_raw_solid_harmonic_symbolic(l - 1, lower_m), source_coordinate
                )
                recurrence_derivatives = apply_gradient_recurrence(l, m, lower_value)
                for derivative_coordinate, recurrence_derivative in zip(
                    coordinates, recurrence_derivatives
                ):
                    expected = diff(target, source_coordinate, derivative_coordinate)
                    difference = simplify(recurrence_derivative - expected)
                    if difference != 0:
                        raise RuntimeError(
                            f"Hessian recurrence verification failed for l={l}, m={m}: "
                            f"{difference}"
                        )


def gen_recurrence_coefficients():
    """Emit the three recurrence coefficients of the gradient recurrence.

    The factor association is pinned to the order the original evaluateVGL used, so
    that factoring the recurrence out of evaluateVGL_impl is bitwise identical: a
    three-way product reassociates under IEEE754, and reordering it perturbs the
    stored gradients by a few ulp. SymPy's printer canonicalizes factor order and
    cannot preserve the association, so the C++ is written out literally and SymPy
    is used to check it against the symbolic form instead of to format it.
    """
    l, ma, fac = Symbol("l"), Symbol("ma"), Symbol("fac")
    coefficients = (
        ("cp", "fac * (l - ma - 1) * (l - ma)", sqrt(fac * (l - ma - 1) * (l - ma))),
        ("cm", "fac * (l + ma - 1) * (l + ma)", sqrt(fac * (l + ma - 1) * (l + ma))),
        ("c0", "fac * (l - ma) * (l + ma)", sqrt(fac * (l - ma) * (l + ma))),
    )
    lines = []
    for name, emitted, expression in coefficients:
        if simplify(sqrt(sympify(emitted)) - expression) != 0:
            raise RuntimeError(f"emitted coefficient {name} does not match its symbolic form")
        lines.append(f"  const T {name:<2} = std::sqrt({emitted});")
    return "\n".join(lines)


def gen_soa_gradient_recurrence():
    """Generate the shared first-derivative recurrence used by VGL and VGH."""
    return """template<typename T>
template<typename Accessor>
inline void SoaSphericalTensor<T>::gradient_recurrence(const int l,
                                                       const int m,
                                                       const T fac,
                                                       Accessor lower,
                                                       T& gx,
                                                       T& gy,
                                                       T& gz)
{
  constexpr T czero(0);
  constexpr T ahalf(0.5);
  const int ma = std::abs(m);
  const int lm = index(l - 1, 0);
%s

  T dpr, dpi, dmr, dmi;
  gz = (l > ma) ? c0 * lower(lm + m) : czero;
  if (l > ma + 1)
  {
    dpr = cp * lower(lm + ma + 1);
    dpi = cp * lower(lm - ma - 1);
  }
  else
  {
    dpr = czero;
    dpi = czero;
  }
  if (l > 1)
  {
    switch (ma)
    {
    case 0:
      dmr = -cm * lower(lm + 1);
      dmi = cm * lower(lm - 1);
      break;
    case 1:
      dmr = cm * lower(lm);
      dmi = czero;
      break;
    default:
      dmr = cm * lower(lm + ma - 1);
      dmi = cm * lower(lm - ma + 1);
    }
  }
  else
  {
    dmr = cm * lower(lm);
    dmi = czero;
  }
  if (m < 0)
  {
    gx = ahalf * (dpi - dmi);
    gy = -ahalf * (dpr + dmr);
  }
  else
  {
    gx = ahalf * (dpr - dmr);
    gy = ahalf * (dpi + dmi);
  }
}
""" % gen_recurrence_coefficients()


def gen_soa_evaluate_vgh():
    """Generate evaluateVGH by reapplying the shared gradient recurrence."""
    return """template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGH(T x, T y, T z)
{
  evaluateVGL(x, y, z);

  constexpr T czero(0);
  constexpr T ahalf(0.5);
  const T* restrict gYlmX = cYlm.data(1);
  const T* restrict gYlmY = cYlm.data(2);
  const T* restrict gYlmZ = cYlm.data(3);
  T* restrict hYlmXX      = cYlm.data(4);
  T* restrict hYlmXY      = cYlm.data(5);
  T* restrict hYlmXZ      = cYlm.data(6);
  T* restrict hYlmYY      = cYlm.data(7);
  T* restrict hYlmYZ      = cYlm.data(8);
  T* restrict hYlmZZ      = cYlm.data(9);

  hYlmXX[0] = czero;
  hYlmXY[0] = czero;
  hYlmXZ[0] = czero;
  hYlmYY[0] = czero;
  hYlmYZ[0] = czero;
  hYlmZZ[0] = czero;

  for (int l = 1; l <= Lmax; ++l)
  {
    const T fac = factor2L_[l];
    for (int m = -l; m <= l; ++m)
    {
      const int lm = index(l, m);

      // The recurrence coefficients do not depend on position, so differentiating
      // a second time is the same recurrence applied to the l-1 gradients that
      // evaluateVGL just stored. Those carry norm_factor_, so strip it on the way
      // in and apply the target harmonic's norm_factor_ on the way out.
      const auto differentiate_gradient = [&](const T* restrict lower_derivative) -> TinyVector<T, 3> {
        T dx, dy, dz;
        gradient_recurrence(
            l, m, fac, [&](const int lower_lm) { return lower_derivative[lower_lm] / norm_factor_[lower_lm]; }, dx, dy,
            dz);
        return {dx, dy, dz};
      };

      // d_gx[i] is the i-th derivative of the stored x gradient, so the mixed
      // second derivatives each come out of two separate recurrence passes.
      const TinyVector<T, 3> d_gx = differentiate_gradient(gYlmX);
      const TinyVector<T, 3> d_gy = differentiate_gradient(gYlmY);
      const TinyVector<T, 3> d_gz = differentiate_gradient(gYlmZ);

      // The Hessian is symmetric analytically, so d_gx[1] and d_gy[0] differ
      // only by roundoff. Both are already in hand, so averaging them is free
      // and keeps the stored Hessian exactly symmetric.
      const T norm = norm_factor_[lm];
      hYlmXX[lm]   = norm * d_gx[0];
      hYlmXY[lm]   = norm * ahalf * (d_gx[1] + d_gy[0]);
      hYlmXZ[lm]   = norm * ahalf * (d_gx[2] + d_gz[0]);
      hYlmYY[lm]   = norm * d_gy[1];
      hYlmYZ[lm]   = norm * ahalf * (d_gy[2] + d_gz[1]);
      hYlmZZ[lm]   = norm * d_gz[2];
    }
  }
}
"""


def gen_soa_evaluate_vghgh():
    return """template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGHGH(T x, T y, T z)
{
  throw std::runtime_error("SoaSphericalTensor<T>::evaluateVGHGH(x,y,z):  Not implemented\\n");
}
"""


def run_template(template_path, output_path, bodies):
    """Apply the same line-oriented template convention as the Cartesian generator."""
    output = ""
    with template_path.open("r", encoding="utf-8") as template_file:
        for line in template_file:
            if line.startswith("%"):
                key = line.strip()[1:]
                if key not in bodies:
                    raise KeyError(f"Template item not found: {key}")
                line = bodies[key]
            output += line

    output_path.write_text(output, encoding="utf-8")


def create_soa_spherical_tensor_h():
    script_dir = Path(__file__).resolve().parent
    template_path = script_dir / "SoaSphericalTensor.h.in"
    output_path = script_dir.parent / "SoaSphericalTensor.h"

    verify_hessian_recurrence()

    warning = """/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/Numerics/codegen/gen_spherical_tensor.py and
 SoaSphericalTensor.h.in.

 Edit the template or generator and rerun gen_spherical_tensor.py.
*/
"""
    bodies = {
        "dire_codegen_warning": warning,
        "gradient_recurrence": gen_soa_gradient_recurrence(),
        "evaluate_vgh": gen_soa_evaluate_vgh(),
        "evaluate_vghgh": gen_soa_evaluate_vghgh(),
    }
    run_template(template_path, output_path, bodies)
    subprocess.run(["clang-format", "-i", output_path], check=True)


if __name__ == "__main__":
    create_soa_spherical_tensor_h()
