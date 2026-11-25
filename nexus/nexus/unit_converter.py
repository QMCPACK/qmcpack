##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

"""Simple unit conversion utilities."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import Self

import numpy as np
import numpy.typing as npt


@dataclass
class Unit:
    val: float
    shift: float = 0.0


class UnitEnum(Unit, Enum):

    def __new__(cls, val: float, shift: float = 0.0) -> Self:
        unit = Unit.__new__(cls)
        unit._value_ = Unit(val, shift)
        return unit

    @classmethod
    def _missing_(cls, value) -> Self:
        value = value.replace(" ", "_")
        if (val := cls.__members__.get(value)) is not None:
            return val
        else:
            return None


class Constant(UnitEnum):
    """Constant values that are not necessarily a base unit. All values are in SI units.

    Attributes
    ----------
    kB : float
        Boltzmann constant in Joules per Kelvin [1]_.
    mol : float
        Avogadro's number [2]_.

    References
    ----------
    .. [1] https://physics.nist.gov/cgi-bin/cuu/Value?k
    .. [2] https://physics.nist.gov/cgi-bin/cuu/Value?na
    """

    kB = 1.380_649e-23
    kb = kB
    mole = 6.022_140_76e23
    mol = mole


class Energy(UnitEnum):
    """Energy conversions. All values are in units of joules.

    Attributes
    ----------
    Joule, J : float
        One joule.
    electron_volt, eV : float
        One electron volt [1]_.
    Hartree, Ha, au : float
        One Hartree [2]_.
    Rydberg, Ry : float
        One Rydberg [3]_.
    kJ_mol : float
        One kilojoule per mole.
    kcal_mol : float
        One kilocalorie per mole [4]_. See also :attr:`Constant.mole`.
    Kelvin, K : float
        Thermal energy per Kelvin. See :attr:`Constant.kB`.
    Celsius, degC : float
        Thermal energy per degree Celsius. Equal to :attr:`K`.
    Fahrenheit, degF : float
        Thermal energy per degree Fahrenheit [5]_.

    References
    ----------
    .. [1] https://physics.nist.gov/cgi-bin/cuu/Value?evj
    .. [2] https://physics.nist.gov/cgi-bin/cuu/Value?hrev
    .. [3] https://physics.nist.gov/cgi-bin/cuu/Value?rydhcj
    .. [4] https://www.nist.gov/glossary-term/26261
    .. [5] https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9#ENERGY
    """

    Joule = 1.0
    J = Joule
    electron_volt = 1.602_176_634e-19
    eV = electron_volt
    Hartree = 4.359_744_722_2060e-18
    Ha = Hartree
    au = Hartree
    Rydberg = 2.179_872_361_1030e-18
    Ry = Rydberg
    kJ_mol = 1_000.0 / Constant.mol.val
    kcal_mol = kJ_mol * 4.184
    Kelvin = Constant.kb.val
    K = Constant.kb.val
    Celsius = Kelvin, -273.15
    degC = Celsius
    Fahrenheit = (Kelvin * 5.0 / 9.0), (9.0 / 5.0 * -273.15) + 32
    degF = Fahrenheit


class Length(UnitEnum):
    """Length conversions. All values are in units of meters.

    Attributes
    ----------
    meter, m : float
        One meter.
    Angstrom, A, ang : float
        One Angstrom [1]_.
    Bohr, B, a0, au : float
        One Bohr radius [2]_.
    millimeter, mm : float
        One millimeter [3]_.
    micrometer, um : float
        One micrometer [3]_.
    nanometer, nm : float
        One nanometer [3]_.
    picometer, pm : float
        One picometer [3]_.
    femtometer, fm : float
        One femtometer [3]_.

    References
    ----------
    .. [1] https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9#LENGTH
    .. [2] https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
    .. [3] https://www.nist.gov/pml/owm/metric-si-prefixes#prefixes
    """

    meter = 1.0
    m = meter
    Angstrom = 1.0e-10
    A = Angstrom
    ang = Angstrom
    Bohr = 0.529_177_208_59e-10
    B = Bohr
    a0 = Bohr
    au = Bohr
    millimeter = 1.0e-3
    mm = millimeter
    micrometer = 1.0e-6
    um = micrometer
    nanometer = 1.0e-9
    nm = nanometer
    picometer = 1.0e-12
    pm = picometer
    femtometer = 1.0e-15
    fm = femtometer


class Time(UnitEnum):
    """Time conversions. All values are in units of seconds.

    Attributes
    ----------
    second, s : float
        One second.
    millisecond, ms : float
        One millisecond [1]_.
    microsecond, us : float
        One microsecond [1]_.
    nanosecond, ns : float
        One nanosecond [1]_.
    picosecond, ps : float
        One picosecond [1]_.
    femtosecond, fs : float
        One femtosecond [1]_.
    au : float
        Atomic unit of time [2]_.

    References
    ----------
    .. [1] https://www.nist.gov/pml/owm/metric-si-prefixes#prefixes
    .. [2] https://physics.nist.gov/cgi-bin/cuu/Value?aut
    """

    second = 1.0
    s = second
    millisecond = 1.0e-3
    ms = millisecond
    microsecond = 1.0e-6
    us = microsecond
    nanosecond = 1.0e-9
    ns = nanosecond
    picosecond = 1.0e-12
    ps = picosecond
    femtosecond = 1.0e-15
    fs = femtosecond
    au = 2.418_884_326_5864e-17


class Mass(UnitEnum):
    """Mass conversions. All values are in units of kilograms.

    Attributes
    ----------
    kilogram, kg : float
        One kilogram.
    electron_mass, me, au : float
        Mass of one electron [2]_.
    proton_mass, mp : float
        Mass of one proton [3]_.
    atomic_mass_unit, amu, Dalton, Da : float
        Atomic mass constant [4]_.

    References
    ----------
    .. [1] https://www.nist.gov/pml/owm/metric-si-prefixes#prefixes
    .. [2] https://physics.nist.gov/cgi-bin/cuu/Value?me
    .. [3] https://physics.nist.gov/cgi-bin/cuu/Value?mp
    .. [4] https://physics.nist.gov/cgi-bin/cuu/Value?u
    """
    kilogram = 1.0
    kg = kilogram
    gram = 1.0e-3
    g = gram
    electron_mass = 9.109_383_7139e-31
    me = electron_mass
    au = me
    proton_mass = 1.672_621_925_95e-27
    mp = proton_mass
    atomic_mass_unit = 1.660_539_068_92e-27
    amu = atomic_mass_unit
    Dalton = atomic_mass_unit
    Da = atomic_mass_unit


class Charge(UnitEnum):
    """Charge conversions. All values are in units of Coulombs.

    Attributes
    ----------
    Coulomb, C : float
        One Coulomb.
    electron_charge, e, au : float
        Elementary charge [1]_.

    References
    ----------
    .. [1] https://physics.nist.gov/cgi-bin/cuu/Value?e
    """

    Coulomb = 1.0
    C = Coulomb
    electron_charge = 1.602_176_634e-19
    e = electron_charge
    au = electron_charge


class Pressure(UnitEnum):
    """Pressure conversions. All values are in units of pascals.

    Attributes
    ----------
    Pascal, Pa : float
        One pascal.
    kiloPascal, kPa : float
        One kilopascal [1]_.
    megaPascal, MPa : float
        One megapascal [1]_.
    gigaPascal, GPa : float
        One gigapascal [1]_.
    bar : float
        One bar [2]_.
    kilobar, kbar : float
        One kilobar [2]_.
    megabar, Mbar : float
        One megabar [2]_.
    gigabar, Gbar : float
        One gigabar [2]_.
    atmosphere, atm : float
        One atmosphere [2]_.
    mmHg : float
        One millimeter of Mercury [2]_.

    References
    ----------
    .. [1] https://www.nist.gov/pml/owm/metric-si-prefixes#prefixes
    .. [2] https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9#PRESSURE
    """

    Pascal = 1.0
    Pa = Pascal
    kiloPascal = 1.0e3
    kPa = kiloPascal
    megaPascal = 1.0e6
    MPa = megaPascal
    gigaPascal = 1.0e9
    GPa = gigaPascal
    bar = 1.0e5
    kilobar = 1.0e8
    kbar = kilobar
    megabar = 1.0e11
    Mbar = megabar
    gigabar = 1.0e14
    Gbar = gigabar
    atmosphere = 1.013_25e5
    atm = atmosphere
    mmHg = 1.333_224e2


class Force(UnitEnum):
    """Force conversions. All values are in units of Newtons.

    Attributes
    ----------
    Newton, N : float
        One newton.
    kiloNewton, kN : float
        One kilonewton [1]_.
    megaNewton, MN : float
        One meganewton [1]_.
    gigaNewton, GN : float
        One giganewton [1]_.
    milliNewton, mN : float
        One millinewton [1]_.
    microNewton, uN : float
        One micronewton [1]_.
    nanoNewton, nN : float
        One nanonewton [1]_.
    picoNewton, pN : float
        One piconewton [1]_.
    Hartree_Bohr, au : float
        One Hartree per bohr radius [2]_.
    Rydberg_Bohr, Ry_a0 : float
        One Rydberg per bohr radius [2]_.

    References
    ----------
    .. [1] https://www.nist.gov/pml/owm/metric-si-prefixes#prefixes
    .. [2] https://physics.nist.gov/cgi-bin/cuu/Value?auforce
    """

    Newton = 1.0
    N = Newton
    kiloNewton = 1e3
    kN = kiloNewton
    megaNewton = 1e6
    MN = megaNewton
    gigaNewton = 1e9
    GN = gigaNewton
    milliNewton = 1e-3
    mN = milliNewton
    microNewton = 1e-6
    uN = microNewton
    nanoNewton = 1e-9
    nN = nanoNewton
    picoNewton = 1e-12
    pN = picoNewton
    Hartree_Bohr = 8.238_723_5038e-8
    au = Hartree_Bohr
    Rydberg_Bohr = 4.119_361_7519e-8
    Ry_a0 = Rydberg_Bohr


class ThermalConductivity(UnitEnum):
    """Thermal conductivity conversions. All values are in units of W/mK.

    Attributes
    ----------
    W_mK : float
        Watts per meter Kelvin.
    """
    W_mK = 1.0


ALL_UNITS: MappingProxyType[str, UnitEnum] = MappingProxyType(
    {
        name: val for unit_enum in UnitEnum.__subclasses__()
        for name, val in unit_enum.__members__.items() if name != "au"
    }
)


type Convertable = float | int | npt.NDArray[np.number]
type UnitStr = str

def convert[V: Convertable, Unit: UnitStr | UnitEnum](
    value: V, source_unit: Unit, target_unit: Unit,
) -> V:
    """Convert a value from a source unit to a target unit.

    Examples
    --------
    Simple float conversion.

    >>> convert(1.0, "A", "B")
    1.8897261328856432

    Array conversion.

    >>> convert(np.eye(3), "Ha", "eV")
    array([[27.21138625,  0.        ,  0.        ],
           [ 0.        , 27.21138625,  0.        ],
           [ 0.        ,  0.        , 27.21138625]])
    """
    if isinstance(source_unit, str):
        source_unit = ALL_UNITS[source_unit.replace(" ", "_")]

    if isinstance(target_unit, str):
        target_unit = ALL_UNITS[target_unit.replace(" ", "_")]

    if source_unit is target_unit:
        return value

    if (su_typ := type(source_unit)) is not (tu_typ := type(target_unit)):
        msg = f"Can not convert between {su_typ.__name__} units and {tu_typ.__name__} units!"
        raise ValueError(msg)

    value_out = (
        (value - source_unit.shift) * source_unit.val
        / target_unit.val
        + target_unit.shift
    )

    return value_out


def convert_scalar_to_all[T: Convertable](
    units: str | UnitEnum, value_orig: T
) -> dict[str, T]:
    if isinstance(units, str):
        unit_type = type(ALL_UNITS[units.replace(" ", "_")])
    else:
        unit_type = type(units)

    value = {"orig": value_orig}
    for unit in unit_type:
        value[unit.name] = convert(value_orig, units, unit)

    return value
