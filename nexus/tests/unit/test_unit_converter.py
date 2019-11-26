

def test_import():
    from unit_converter import convert
#end def test_import



def test_convert():
    import numpy as np
    from testing import value_eq
    from unit_converter import convert


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

    assert(value_eq(kg_per_electron,9.10938291e-31 ,rtol=1e-8))
    assert(value_eq(kg_per_proton  ,1.672621777e-27,rtol=1e-8))
    assert(value_eq(kg_per_amu     ,1.660538921e-27,rtol=1e-8))

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

    assert(value_eq(joules_per_eV  ,1.60217646e-19   ,rtol=1e-8))
    assert(value_eq(eV_per_rydberg ,13.6056923       ,rtol=1e-8))
    assert(value_eq(eV_per_hartree ,27.2113846       ,rtol=1e-8))
    assert(value_eq(eV_per_kJ_mol  ,0.0103642695083  ,rtol=1e-8))
    assert(value_eq(eV_per_kcal_mol,0.04336411531    ,rtol=1e-8))
    assert(value_eq(eV_per_kelvin  ,8.61734231197e-05,rtol=1e-8))
    assert(value_eq(kelvin_per_eV  ,11604.5059346    ,rtol=1e-8))


    # temperature
    assert(value_eq(convert(100,'degC','K'   ),373.15))
    assert(value_eq(convert(  0,'degC','K'   ),273.15))
    assert(value_eq(convert(212,'degF','K'   ),373.15))
    assert(value_eq(convert( 32,'degF','K'   ),273.15))
    assert(value_eq(convert(100,'degC','degF'),212.0 ))
    assert(value_eq(convert(  0,'degC','degF'), 32.0 ))
    assert(value_eq(convert(212,'degF','degC'),100.0 ))
    assert(value_eq(convert( 32,'degF','degC'),  0.0,atol=1e-8))

#end def test_convert



def test_convert_scalar_to_all():
    from testing import value_eq
    from unit_converter import UnitConverter

    eV_to = {
        'J'       : 1.60217646e-19, 
        'Ha'      : 0.03674932439858279, 
        'Ry'      : 0.07349864879716558, 
        'eV'      : 1.0,
        'kcal_mol': 23.0605419446755, 
        'kJ_mol'  : 96.48533350089092, 
        'K'       : 11604.505934630948, 
        'degC'    : 11331.355934630948, 
        'degF'    : 20428.440682335706, 
        }

    v = UnitConverter.convert_scalar_to_all('eV',1.0)

    for unit in eV_to.keys():
        assert(value_eq(v[unit],eV_to[unit]))
    #end for

#end def test_convert_scalar_to_all
