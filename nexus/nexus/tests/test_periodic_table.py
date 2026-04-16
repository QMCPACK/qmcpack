try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.PERIODIC_TABLE)
except ImportError:
    pass


def test_import():
    from .. import periodic_table
    from ..periodic_table import Elements
#end def test_import



def test_periodic_table():
    from ..testing import value_eq
    from ..periodic_table import Elements

    ref_element_symbols = (
        "Xx", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
        "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
        "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
        "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
        "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
        "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
        "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
        "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    )

    ref_atomic_numbers = tuple(range(0,len(ref_element_symbols)))

    element_symbols = [e.symbol for e in Elements]
    atomic_numbers = [e.atomic_number for e in Elements]

    assert(len(element_symbols) == len(ref_element_symbols))
    assert(len(atomic_numbers) == len(ref_atomic_numbers))

    for elem in element_symbols:
        assert(elem in ref_element_symbols)

    for number in atomic_numbers:
        assert(number in ref_atomic_numbers)
    
    ref_carbon_name     = "Carbon"
    ref_carbon_symbol   = "C"
    ref_carbon_number   = 6
    ref_carbon_weight   = 12.011
    ref_carbon_group    = 14
    ref_carbon_isotopes = {12:12.0000000, 13:13.00335483507, 14:14.0032419884}

    assert(Elements.Carbon.name          == ref_carbon_name)
    assert(Elements.Carbon.symbol        == ref_carbon_symbol)
    assert(str(Elements.Carbon)          == ref_carbon_symbol)
    assert(Elements.Carbon.atomic_number == ref_carbon_number)
    assert(Elements.Carbon.atomic_weight == ref_carbon_weight)
    assert(Elements.Carbon.group         == ref_carbon_group)
    assert(Elements.Carbon.isotopes      == ref_carbon_isotopes)

    assert(Elements.C.name          == ref_carbon_name)
    assert(Elements.C.symbol        == ref_carbon_symbol)
    assert(str(Elements.C)          == ref_carbon_symbol)
    assert(Elements.C.atomic_number == ref_carbon_number)
    assert(Elements.C.atomic_weight == ref_carbon_weight)
    assert(Elements.C.group         == ref_carbon_group)
    assert(Elements.C.isotopes      == ref_carbon_isotopes)

    assert(Elements.Carbon is Elements.C)

    assert(Elements.num_elements() == 118)
#end def test_periodic_table


def test_call_elements():
    from ..periodic_table import Elements

    # Good calls
    assert(Elements("Hydrogen") is Elements.Hydrogen)
    assert(Elements("H") is Elements.Hydrogen)
    assert(Elements(1) is Elements.Hydrogen)

    # Calls that need to go through `_missing_`
    # Improper case
    assert(Elements("hydrogen") is Elements.Hydrogen)
    assert(Elements("h") is Elements.Hydrogen)
    # Trailing and leading whitespace
    assert(Elements("Hydrogen ") is Elements.Hydrogen)
    assert(Elements(" H") is Elements.Hydrogen)
    # One step from good
    assert(Elements("1") is Elements.Hydrogen)
    assert(Elements(1.0) is Elements.Hydrogen)

    # Another spot check, same situations
    assert(Elements("Ruthenium") is Elements.Ruthenium)
    assert(Elements("Ru") is Elements.Ruthenium)
    assert(Elements(44) is Elements.Ruthenium)

    assert(Elements("ruthenium") is Elements.Ruthenium)
    assert(Elements("ru") is Elements.Ruthenium)
    assert(Elements("Ruthenium ") is Elements.Ruthenium)
    assert(Elements(" Ru") is Elements.Ruthenium)
    assert(Elements("44") is Elements.Ruthenium)
    assert(Elements(44.0) is Elements.Ruthenium)
#end def test_call_elements


def test_is_element():
    from ..periodic_table import Elements

    ref_symbols = (
        "Xx", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
        "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
        "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
        "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
        "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
        "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
        "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
        "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    )

    ref_elements = (
        Elements.Xx,
        Elements.H,  Elements.He, Elements.Li, Elements.Be, Elements.B,
        Elements.C,  Elements.N,  Elements.O,  Elements.F,  Elements.Ne,
        Elements.Na, Elements.Mg, Elements.Al, Elements.Si, Elements.P,
        Elements.S,  Elements.Cl, Elements.Ar, Elements.K,  Elements.Ca,
        Elements.Sc, Elements.Ti, Elements.V,  Elements.Cr, Elements.Mn,
        Elements.Fe, Elements.Co, Elements.Ni, Elements.Cu, Elements.Zn,
        Elements.Ga, Elements.Ge, Elements.As, Elements.Se, Elements.Br,
        Elements.Kr, Elements.Rb, Elements.Sr, Elements.Y,  Elements.Zr,
        Elements.Nb, Elements.Mo, Elements.Tc, Elements.Ru, Elements.Rh,
        Elements.Pd, Elements.Ag, Elements.Cd, Elements.In, Elements.Sn,
        Elements.Sb, Elements.Te, Elements.I,  Elements.Xe, Elements.Cs,
        Elements.Ba, Elements.La, Elements.Ce, Elements.Pr, Elements.Nd,
        Elements.Pm, Elements.Sm, Elements.Eu, Elements.Gd, Elements.Tb,
        Elements.Dy, Elements.Ho, Elements.Er, Elements.Tm, Elements.Yb,
        Elements.Lu, Elements.Hf, Elements.Ta, Elements.W,  Elements.Re,
        Elements.Os, Elements.Ir, Elements.Pt, Elements.Au, Elements.Hg,
        Elements.Tl, Elements.Pb, Elements.Bi, Elements.Po, Elements.At,
        Elements.Rn, Elements.Fr, Elements.Ra, Elements.Ac, Elements.Th,
        Elements.Pa, Elements.U,  Elements.Np, Elements.Pu, Elements.Am,
        Elements.Cm, Elements.Bk, Elements.Cf, Elements.Es, Elements.Fm,
        Elements.Md, Elements.No, Elements.Lr, Elements.Rf, Elements.Db,
        Elements.Sg, Elements.Bh, Elements.Hs, Elements.Mt, Elements.Ds,
        Elements.Rg, Elements.Cn, Elements.Nh, Elements.Fl, Elements.Mc,
        Elements.Lv, Elements.Ts, Elements.Og,
    )

    for symbol, element in zip(ref_symbols, ref_elements):
        assert(Elements.is_element(symbol)) # True for symbols
        assert(Elements.is_element(element)) # True for members

        is_elem, elem = Elements.is_element(element, return_element=True)
        assert(is_elem)
        assert(elem.symbol is symbol)
        assert(elem is element)
        is_elem, elem = Elements.is_element(symbol, return_element=True)
        assert(is_elem)
        assert(elem.symbol is symbol)
        assert(elem is element)
    #end for
    
    carbon_strs = (
        "C",
        "C1",
        "C2",
        "C12",
        "C123",
        "C_1",
        "C_2",
        "C_12",
        "C_123",
        "C_a",
        "C_abc",
        "C-1",
        "C-2",
        "C-12",
        "C-123",
        "C-a",
        "C-abc",
        "c",
        "c1",
        "c2",
        "c12",
        "c123",
        "c_1",
        "c_2",
        "c_12",
        "c_123",
        "c_a",
        "c_abc",
        "c-1",
        "c-2",
        "c-12",
        "c-123",
        "c-a",
        "c-abc",
    )

    for string in carbon_strs:
        assert(Elements.is_element(string))

        is_elem, elem = Elements.is_element(string, return_element=True)
        assert(is_elem)
        assert(elem is Elements.Carbon)
        assert(elem.symbol == "C")
        assert(elem.name == "Carbon")
    #end for

    cobalt_strs = [
        "Co",
        "Co1",
        "Co2",
        "Co12",
        "Co123",
        "Co_1",
        "Co_2",
        "Co_12",
        "Co_123",
        "Co_a",
        "Co_abc",
        "Co-1",
        "Co-2",
        "Co-12",
        "Co-123",
        "Co-a",
        "Co-abc",
        "co",
        "co1",
        "co2",
        "co12",
        "co123",
        "co_1",
        "co_2",
        "co_12",
        "co_123",
        "co_a",
        "co_abc",
        "co-1",
        "co-2",
        "co-12",
        "co-123",
        "co-a",
        "co-abc",
    ]
    for string in cobalt_strs:
        assert(Elements.is_element(string))

        is_elem, element = Elements.is_element(string, return_element=True)
        assert(is_elem)
        assert(element is Elements.Cobalt)
        assert(element.symbol == "Co")
        assert(element.name == "Cobalt")
    #end for
#end def test_is_element


def test_element_set():
    from ..periodic_table import Elements
    ref_set = set([
        Elements.Xx,
        Elements.H,
        Elements.Dy,
        Elements.U,
        Elements.Nh,
    ])

    element_set = set([
        Elements.Xx,
        Elements.H,  Elements.H,  Elements.H,
        Elements.Dy,
        Elements.U,  Elements.U,  Elements.U,  Elements.U,  Elements.U,
        Elements.Nh, Elements.Nh, Elements.Nh, Elements.Nh,
    ])

    assert(ref_set == element_set)


def test_representation():
    from ..periodic_table import Elements

    ref_repr = "<Elements.Carbon: symbol='C', atomic_number=6, atomic_weight=12.011, group=14>"
    assert(repr(Elements.Carbon) == ref_repr)
