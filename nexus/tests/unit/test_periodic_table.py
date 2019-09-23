

def test_import():
    import periodic_table
    from periodic_table import pt,is_element
#end def test_import



def test_periodic_table():
    from testing import value_eq
    from periodic_table import pt

    elements = set('''
        H He 
        Li Be B C N O F Ne 
        Na Mg Al Si P S Cl Ar 
        K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr 
        Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe 
        Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn 
        Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr 
        '''.split())

    atomic_numbers = set(range(1,len(elements)+1))
    
    missing = elements - set(pt.keys())
    assert(len(missing)==0)
    
    missing = elements - set(pt.elements.keys())
    assert(len(missing)==0)

    missing = atomic_numbers - set(pt.simple_elements.keys())
    assert(len(missing)==0)

    fields = '''
        atomic_weight
        atomic_radius
        nuclear_charge
        electron_affinity
        ionization_energy
        ionic_radius
        thermal_cond
        melting_point
        boiling_point
        '''.split()

    for e in elements:
        elem = pt.elements[e]
        assert(id(elem)==id(pt[e]))
        selem = pt.simple_elements[elem.atomic_number]
        for f in fields:
            assert(value_eq(selem[f],elem[f].orig))
        #end for
    #end for

    C = pt.C
    assert(C.atomic_number==6)
    assert(C.symbol=='C')
    assert(value_eq(C.atomic_radius.pm        ,77.2))
    assert(value_eq(C.atomic_weight.amu       ,12.011))
    assert(value_eq(C.boiling_point.degC      ,4827.0))
    assert(value_eq(C.electron_affinity.kJ_mol,121.9))
    assert(value_eq(C.ionic_radius.pm,260.0))
    assert(value_eq(C.ionization_energy.eV,11.26))
    assert(value_eq(C.melting_point.degC,3550.0))
    assert(value_eq(C.nuclear_charge.e,3.25))
    assert(value_eq(C.thermal_cond.W_mK,1960.0))

#end def test_periodic_table



def test_is_element():
    from periodic_table import is_element

    elements = '''
        H He 
        Li Be B C N O F Ne 
        Na Mg Al Si P S Cl Ar 
        K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr 
        Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe 
        Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn 
        Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr 
        '''.split()


    for e in elements:
        assert(is_element(e))
        is_elem,symbol = is_element(e,symbol=True)
        assert(is_elem)
        assert(symbol==e)
    #end for
    
    valid = '''
        C
        C1
        C2
        C12
        C123
        C_1
        C_2
        C_12
        C_123
        C_a
        C_abc
        '''.split()
    for e in valid:
        assert(is_element(e))
        is_elem,symbol = is_element(e,symbol=True)
        assert(is_elem)
        assert(symbol=='C')
    #end for
    
    valid = '''
        Co
        Co1
        Co2
        Co12
        Co123
        Co_1
        Co_2
        Co_12
        Co_123
        Co_a
        Co_abc
        '''.split()
    for e in valid:
        assert(is_element(e))
        is_elem,symbol = is_element(e,symbol=True)
        assert(is_elem)
        assert(symbol=='Co')
    #end for
        
#end def test_is_element
