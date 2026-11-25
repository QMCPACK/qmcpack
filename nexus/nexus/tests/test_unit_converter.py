import numpy as np
import pytest

from ..testing import value_eq
from ..unit_converter import Energy, Length, convert, convert_scalar_to_all
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.UNIT_CONVERTER)

def test_convert():
    # distance
    meters_per_angstrom = convert(1.0,'A','m')
    meters_per_bohr     = convert(1.0,'B','m')
    angstrom_per_bohr   = convert(1.0,'B','A')

    B = .52917720859e-10
    B_per_A = B*1e10
    assert(value_eq(meters_per_angstrom,1e-10))
    assert(value_eq(meters_per_bohr,B,atol=1e-10))
    assert(value_eq(angstrom_per_bohr,B_per_A,atol=1e-10))

    v  = 2.34
    vc = convert(v,'A','B')
    assert(value_eq(vc,v/B_per_A))

    vc = convert(v,'B','A')
    assert(value_eq(vc,v*B_per_A))

    v  = 2.34*np.arange(5)
    vc = convert(v,'A','B')
    assert(value_eq(vc,v/B_per_A))

    vc = convert(v,'B','A')
    assert(value_eq(vc,v*B_per_A))


    # mass
    kg_per_electron = convert(1.0,'me' ,'kg')
    kg_per_proton   = convert(1.0,'mp' ,'kg')
    kg_per_amu      = convert(1.0,'amu','kg')

    assert(value_eq(kg_per_electron,9.10938371e-31 ,rtol=1e-8))
    assert(value_eq(kg_per_proton  ,1.672621926e-27,rtol=1e-8))
    assert(value_eq(kg_per_amu     ,1.660539069e-27,rtol=1e-8))

    electrons_per_proton = convert(1.0,'mp','me')
    assert(value_eq(electrons_per_proton,1836.15267195))


    # energy
    joules_per_eV   = convert(1.0,'eV','J')
    joules_per_K    = convert(1.0,'K' ,'J')
    eV_per_rydberg  = convert(1.0,'Ry','eV')
    eV_per_hartree  = convert(1.0,'Ha','eV')
    eV_per_kJ_mol   = convert(1.0,'kJ_mol','eV')
    eV_per_kcal_mol = convert(1.0,'kcal_mol','eV')
    eV_per_kelvin   = convert(1.0,'K','eV')
    kelvin_per_eV   = convert(1.0,'eV','K')

    assert(value_eq(joules_per_eV  ,1.60217663e-19   ,rtol=1e-8))
    assert(value_eq(joules_per_K   ,1.380649e-23     ,rtol=1e-8))
    assert(value_eq(eV_per_rydberg ,13.605693123     ,rtol=1e-8))
    assert(value_eq(eV_per_hartree ,27.2113862       ,rtol=1e-8))
    assert(value_eq(eV_per_kJ_mol  ,0.0103642696563  ,rtol=1e-8))
    assert(value_eq(eV_per_kcal_mol,0.04336410424    ,rtol=1e-8))
    assert(value_eq(eV_per_kelvin  ,8.61733326215e-05,rtol=1e-8))
    assert(value_eq(kelvin_per_eV  ,11604.5181216    ,rtol=1e-8))


    # temperature
    assert(value_eq(convert(100,'degC','K'   ),373.15))
    assert(value_eq(convert(  0,'degC','K'   ),273.15))
    assert(value_eq(convert(212,'degF','K'   ),373.15))
    assert(value_eq(convert( 32,'degF','K'   ),273.15))
    assert(value_eq(convert(100,'degC','degF'),212.0 ))
    assert(value_eq(convert(  0,'degC','degF'), 32.0 ))
    assert(value_eq(convert(212,'degF','degC'),100.0 ))
    assert(value_eq(convert( 32,'degF','degC'),  0.0,atol=1e-8))

    Ha_to_eV = convert(1.0, "Ha", "eV")
    assert(convert(1.0, "Ha", "electron volt") == Ha_to_eV)
    assert(convert(1.0, Energy.Ha, Energy.eV) == Ha_to_eV)
    assert(convert(1.0, Energy("Ha"), Energy("electron volt")) == Ha_to_eV)

    with pytest.raises(
        ValueError,
        match="Can not convert between Energy units and Length units!",
    ):
        convert(1.0, "Ha", "meter")

    with pytest.raises(
        ValueError,
        match="'bean' is not a valid Energy",
    ):
        Energy("bean")
#end def test_convert



def test_convert_scalar_to_all():
    eV_to = {
        'Joule':         1.602176634e-19,
        'Hartree':       0.03674932217566422,
        'Rydberg':       0.07349864435132844,
        'electron_volt': 1.0,
        'kcal_mol':      23.06054783061903,
        'kJ_mol':        96.48533212331002,
        'Kelvin':        11604.518121550082,
        'Celsius':       11331.368121550082,
        'Fahrenheit':    20428.46261879015,
    }

    v = convert_scalar_to_all('eV', 1.0)

    for unit, val in eV_to.items():
        assert(value_eq(v[unit], val))

    meter_to = {
        'meter': 1.0,
        'Angstrom': 10000000000.0,
        'Bohr': 18897261328.85643,
        'millimeter': 1000.0,
        'micrometer': 1000000.0,
        'nanometer':1000000000.0,
        'picometer': 1000000000000.0,
        'femtometer': 1000000000000000.0,
    }

    v = convert_scalar_to_all(Length.meter, 1.0)

    for unit, val in meter_to.items():
        assert(value_eq(v[unit], val))
#end def test_convert_scalar_to_all
