import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PHYSICAL_SYSTEM)

from ..generic import generic_settings
generic_settings.raise_error = True

import numpy as np
import numpy.typing as npt
from ..physical_system import IonSpecies, PhysicalSystem, Electrons, Positrons
from ..physical_system import generate_physical_system
from ..periodic_table import Elements
from ..structure import Structure
from ..unit_converter import convert


def get_LaAlO3_references() -> dict[str, list | npt.NDArray[np.floating] | dict[str, IonSpecies] | Electrons | Positrons]:
    """Get reference elem, pos, axes, ions, electrons, and positrons for LaAlO3.

    This includes both folded and 2x2x2-tiled versions.
    """
    ref_folded_elem = [
        "La", "Al", "O1", "O2", "O3",
        ]
    ref_tiled_elem = [
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        "La", "Al", "O1", "O2", "O3",
        ]
    ref_folded_axes = np.array([
        [3.780, 0.000, 0.000],
        [0.000, 3.780, 0.000],
        [0.000, 0.000, 3.780],
        ], dtype=float)
    ref_tiled_axes = np.array([
        [7.560,  0.000,  0.000],
        [0.000,  7.560,  0.000],
        [0.000,  0.000,  7.560],
        ], dtype=float)
    ref_folded_pos = np.array([
        [0.000,  0.000,  0.000],
        [1.890,  1.890,  1.890],
        [1.890,  1.890,  0.000],
        [0.000,  1.890,  1.890],
        [1.890,  0.000,  1.890],
        ], dtype=float)
    ref_tiled_pos = np.array([
        [0.000,  0.000,  0.000],
        [1.890,  1.890,  1.890],
        [1.890,  1.890,  0.000],
        [0.000,  1.890,  1.890],
        [1.890,  0.000,  1.890],
        [3.780,  0.000,  0.000],
        [5.670,  1.890,  1.890],
        [5.670,  1.890,  0.000],
        [3.780,  1.890,  1.890],
        [5.670,  0.000,  1.890],
        [0.000,  3.780,  0.000],
        [1.890,  5.670,  1.890],
        [1.890,  5.670,  0.000],
        [0.000,  5.670,  1.890],
        [1.890,  3.780,  1.890],
        [3.780,  3.780,  0.000],
        [5.670,  5.670,  1.890],
        [5.670,  5.670,  0.000],
        [3.780,  5.670,  1.890],
        [5.670,  3.780,  1.890],
        [0.000,  0.000,  3.780],
        [1.890,  1.890,  5.670],
        [1.890,  1.890,  3.780],
        [0.000,  1.890,  5.670],
        [1.890,  0.000,  5.670],
        [3.780,  0.000,  3.780],
        [5.670,  1.890,  5.670],
        [5.670,  1.890,  3.780],
        [3.780,  1.890,  5.670],
        [5.670,  0.000,  5.670],
        [0.000,  3.780,  3.780],
        [1.890,  5.670,  5.670],
        [1.890,  5.670,  3.780],
        [0.000,  5.670,  5.670],
        [1.890,  3.780,  5.670],
        [3.780,  3.780,  3.780],
        [5.670,  5.670,  5.670],
        [5.670,  5.670,  3.780],
        [3.780,  5.670,  5.670],
        [5.670,  3.780,  5.670],
        ], dtype=float)

    ref_folded_ions = {
        'Al': IonSpecies(element=Elements.Aluminum,  count=1, label="Al", formal_charge=0, magnetization=2.5, Zeff=3),
        'La': IonSpecies(element=Elements.Lanthanum, count=1, label="La", formal_charge=0, magnetization=3.5, Zeff=11),
        'O1': IonSpecies(element=Elements.Oxygen,    count=1, label="O1", formal_charge=0, magnetization=0,   Zeff=6),
        'O2': IonSpecies(element=Elements.Oxygen,    count=1, label="O2", formal_charge=0, magnetization=0,   Zeff=6),
        'O3': IonSpecies(element=Elements.Oxygen,    count=1, label="O3", formal_charge=0, magnetization=0,   Zeff=6),
        }
    ref_tiled_ions = {
        'Al': IonSpecies(element=Elements.Aluminum,  count=8, label="Al", formal_charge=0, magnetization=2.5, Zeff=3),
        'La': IonSpecies(element=Elements.Lanthanum, count=8, label="La", formal_charge=0, magnetization=3.5, Zeff=11),
        'O1': IonSpecies(element=Elements.Oxygen,    count=8, label="O1", formal_charge=0, magnetization=0,   Zeff=6),
        'O2': IonSpecies(element=Elements.Oxygen,    count=8, label="O2", formal_charge=0, magnetization=0,   Zeff=6),
        'O3': IonSpecies(element=Elements.Oxygen,    count=8, label="O3", formal_charge=0, magnetization=0,   Zeff=6),
        }
    ref_folded_electrons = Electrons(count=32,  n_unpaired=0, spin_orbit=False)
    ref_tiled_electrons  = Electrons(count=256, n_unpaired=0, spin_orbit=False)

    ref_folded_positrons = Positrons(count=8,  n_unpaired=0, spin_orbit=False)
    ref_tiled_positrons  = Positrons(count=64, n_unpaired=0, spin_orbit=False)

    return {
        "ref_folded_elem"      : ref_folded_elem,
        "ref_tiled_elem"       : ref_tiled_elem,
        "ref_folded_axes"      : ref_folded_axes,
        "ref_tiled_axes"       : ref_tiled_axes,
        "ref_folded_pos"       : ref_folded_pos,
        "ref_tiled_pos"        : ref_tiled_pos,
        "ref_folded_ions"      : ref_folded_ions,
        "ref_tiled_ions"       : ref_tiled_ions,
        "ref_folded_electrons" : ref_folded_electrons,
        "ref_tiled_electrons"  : ref_tiled_electrons,
        "ref_folded_positrons" : ref_folded_positrons,
        "ref_tiled_positrons"  : ref_tiled_positrons,
        }
#end def get_LaAlO3_references


def get_d2_d8_references()-> dict[str, list[str] | npt.NDArray[np.floating] | tuple[str]]:
    # Primitive Diamond Cell
    d2_ref_axes = np.array([
        [1.785, 1.785, 0.000],
        [0.000, 1.785, 1.785],
        [1.785, 0.000, 1.785],
        ], dtype=float)
    d2_ref_elem = ["C", "C"]
    d2_ref_pos  = np.array([
        [0.0000, 0.0000, 0.0000],
        [0.8925, 0.8925, 0.8925],
        ], dtype=float)
    d2_ref_units  = "A"
    d2_ref_bconds = ("p", "p", "p")

    d2_structure = Structure(
        axes   = d2_ref_axes,
        elem   = d2_ref_elem,
        pos    = d2_ref_pos,
        units  = d2_ref_units,
        bconds = d2_ref_bconds,
    )

    d2_ref_ions = {
        "C": IonSpecies(
            element       = Elements.Carbon,
            count         = 2,
            label         = "C",
            formal_charge = 0,
            magnetization     = 0,
            Zeff          = 4,
            ),
    }
    d2_ref_electrons = Electrons(count=8, n_unpaired=0, spin_orbit=False)

    # Conventional Diamond Cell
    d8_ref_axes = np.array([
        [3.57, 0.00, 0.00],
        [0.00, 3.57, 0.00],
        [0.00, 0.00, 3.57],
        ], dtype=float)
    d8_ref_elem = ["C", "C", "C", "C", "C", "C", "C", "C"]
    d8_ref_pos  = np.array([
        [0.0000, 0.0000, 0.0000],
        [0.8925, 0.8925, 0.8925],
        [1.7850, 1.7850, 0.0000],
        [2.6775, 2.6775, 0.8925],
        [0.0000, 1.7850, 1.7850],
        [0.8925, 2.6775, 2.6775],
        [1.7850, 0.0000, 1.7850],
        [2.6775, 0.8925, 2.6775],
        ], dtype=float)
    # The conventional cell in `known_crystals` and the cell produced
    # by `Structure.tile()` can have slightly different orderings
    d8_ref_pos_known = np.array([
        [0.0000, 0.0000, 0.0000],
        [0.8925, 0.8925, 0.8925],
        [0.0000, 1.7850, 1.7850],
        [0.8925, 2.6775, 2.6775],
        [1.7850, 0.0000, 1.7850],
        [2.6775, 0.8925, 2.6775],
        [1.7850, 1.7850, 0.0000],
        [2.6775, 2.6775, 0.8925],
        ], dtype=float)
    d8_ref_units = "A"
    d8_ref_bconds = ("p", "p", "p")

    d8_structure = Structure(
        axes   = d8_ref_axes,
        elem   = d8_ref_elem,
        pos    = d8_ref_pos,
        units  = d8_ref_units,
        bconds = d8_ref_bconds,
    )
    d8_ref_ions = {
        "C": IonSpecies(
            element       = Elements.Carbon,
            count         = 8,
            label         = "C",
            formal_charge = 0,
            magnetization     = 0,
            Zeff          = 4,
            ),
    }
    d8_ref_electrons = Electrons(count=32, n_unpaired=0, spin_orbit=False)

    return {
        "d2_ref_axes":      d2_ref_axes,
        "d2_ref_elem":      d2_ref_elem,
        "d2_ref_pos":       d2_ref_pos,
        "d2_ref_units":     d2_ref_units,
        "d2_ref_bconds":    d2_ref_bconds,
        "d2_structure":     d2_structure,
        "d2_ref_ions":      d2_ref_ions,
        "d2_ref_electrons": d2_ref_electrons,
        "d8_ref_axes":      d8_ref_axes,
        "d8_ref_elem":      d8_ref_elem,
        "d8_ref_pos":       d8_ref_pos,
        "d8_ref_pos_known": d8_ref_pos_known,
        "d8_ref_units":     d8_ref_units,
        "d8_ref_bconds":    d8_ref_bconds,
        "d8_structure":     d8_structure,
        "d8_ref_ions":      d8_ref_ions,
        "d8_ref_electrons": d8_ref_electrons,
        }
#end def get_d2_d8_references


def test_d2_d8_references():
    # We don't test the ions and electrons because they aren't put into
    # another object. The structures, however, should be verified.
    refs = get_d2_d8_references()

    d2_ref_axes   = refs["d2_ref_axes"]
    d2_ref_elem   = refs["d2_ref_elem"]
    d2_ref_pos    = refs["d2_ref_pos"]
    d2_ref_units  = refs["d2_ref_units"]
    d2_ref_bconds = refs["d2_ref_bconds"]
    d2_structure  = refs["d2_structure"]
    np.testing.assert_allclose(d2_structure.axes, d2_ref_axes)
    np.testing.assert_allclose(d2_structure.pos, d2_ref_pos)
    assert(d2_structure.elem.tolist()  == d2_ref_elem)
    assert(d2_structure.units          == d2_ref_units)
    assert(tuple(d2_structure.bconds)  == d2_ref_bconds)

    d8_ref_axes   = refs["d8_ref_axes"]
    d8_ref_elem   = refs["d8_ref_elem"]
    d8_ref_pos    = refs["d8_ref_pos"]
    d8_ref_units  = refs["d8_ref_units"]
    d8_ref_bconds = refs["d8_ref_bconds"]
    d8_structure  = refs["d8_structure"]
    np.testing.assert_allclose(d8_structure.axes, d8_ref_axes)
    np.testing.assert_allclose(d8_structure.pos, d8_ref_pos)
    assert(d8_structure.elem.tolist()  == d8_ref_elem)
    assert(d8_structure.units          == d8_ref_units)
    assert(tuple(d8_structure.bconds)  == d8_ref_bconds)
#end def test_d2_d8_references


def test_electrons():
    ref_charge       = -16
    ref_multiplicity = 3
    ref_n_up         = 9
    ref_n_down       = 7

    electrons = Electrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(electrons.total_charge    == ref_charge)
    assert(electrons.multiplicity    == ref_multiplicity)
    assert(electrons.n_up            == ref_n_up)
    assert(electrons.n_down          == ref_n_down)
    assert(not electrons.is_fractional())

    ref_charge       = -15
    ref_multiplicity = 3
    ref_n_up         = 8.5
    ref_n_down       = 6.5

    electrons = Electrons(count=15, n_unpaired=2, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)

    ref_charge       = -16
    ref_multiplicity = 2
    ref_n_up         = 8.5
    ref_n_down       = 7.5

    electrons = Electrons(count=16, n_unpaired=1, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)

    ref_charge       = -15
    ref_multiplicity = 2
    ref_n_up         = 8
    ref_n_down       = 7

    electrons = Electrons(count=15, n_unpaired=1, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)

    ref_charge       = -15
    ref_multiplicity = 2
    ref_n_up         = 7
    ref_n_down       = 8

    electrons = Electrons(count=15, n_unpaired=-1, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)
#end def test_electrons


def test_electrons_eq():
    ref_charge       = -16
    ref_multiplicity = 3
    ref_n_up         = 9
    ref_n_down       = 7

    electrons1 = Electrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(electrons1.total_charge == ref_charge)
    assert(electrons1.multiplicity == ref_multiplicity)
    assert(electrons1.n_up         == ref_n_up)
    assert(electrons1.n_down       == ref_n_down)
    assert(not electrons1.is_fractional())

    electrons2 = Electrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(electrons2.total_charge == ref_charge)
    assert(electrons2.multiplicity == ref_multiplicity)
    assert(electrons2.n_up         == ref_n_up)
    assert(electrons2.n_down       == ref_n_down)
    assert(not electrons2.is_fractional())

    assert(electrons1 == electrons2)
#end def test_electrons_eq


def test_electrons_repr():
    ref_repr = "Electrons(count=16, n_unpaired=2, spin_orbit=False)"
    electrons = Electrons(count=16, n_unpaired=2, spin_orbit=False)
    assert(repr(electrons) == ref_repr)
#end def test_electrons_repr


def test_positrons():
    """More minimal test since the code is the same as for ``Electrons``."""

    ref_charge       = 16
    ref_multiplicity = 3
    ref_n_up         = 9
    ref_n_down       = 7

    positrons1 = Positrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(positrons1.total_charge == ref_charge)
    assert(positrons1.multiplicity == ref_multiplicity)
    assert(positrons1.n_up         == ref_n_up)
    assert(positrons1.n_down       == ref_n_down)
    assert(not positrons1.is_fractional())

    positrons2 = Positrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(positrons2.total_charge == ref_charge)
    assert(positrons2.multiplicity == ref_multiplicity)
    assert(positrons2.n_up         == ref_n_up)
    assert(positrons2.n_down       == ref_n_down)
    assert(not positrons2.is_fractional())

    assert(positrons1 == positrons2)
#end def test_positrons


def test_positrons_repr():
    ref_repr = "Positrons(count=16, n_unpaired=2, spin_orbit=False)"
    positrons = Positrons(count=16, n_unpaired=2, spin_orbit=False)
    assert(repr(positrons) == ref_repr)
#end def test_positrons_repr


def test_electron_positron_neq():

    ref_positron_charge       = 16
    ref_positron_multiplicity = 3
    ref_positron_n_up         = 9
    ref_positron_n_down       = 7

    positrons = Positrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(positrons.total_charge == ref_positron_charge)
    assert(positrons.multiplicity == ref_positron_multiplicity)
    assert(positrons.n_up         == ref_positron_n_up)
    assert(positrons.n_down       == ref_positron_n_down)
    assert(not positrons.is_fractional())

    ref_electron_charge       = -16
    ref_electron_multiplicity = 3
    ref_electron_n_up         = 9
    ref_electron_n_down       = 7

    electrons = Electrons(count=16, n_unpaired=2, spin_orbit=False)

    assert(electrons.total_charge == ref_electron_charge)
    assert(electrons.multiplicity == ref_electron_multiplicity)
    assert(electrons.n_up         == ref_electron_n_up)
    assert(electrons.n_down       == ref_electron_n_down)
    assert(not electrons.is_fractional())

    assert(electrons != positrons)
#end def test_electron_positron_neq


def test_custom_ion_species():
    ion = IonSpecies(
        element       = Elements.Iron,
        count         = 12,
        label         = "Fe1",
        formal_charge = 2,
        magnetization = 1,
        Zeff          = 16,
        )

    assert(ion.element              is Elements.Iron)
    assert(ion.count                == 12)
    assert(ion.label                == "Fe1")
    assert(ion.formal_charge        == 2)
    assert(ion.magnetization        == 1)
    assert(ion.Zeff                 == 16)
    assert(ion.pseudized()          is True)
    assert(ion.is_ghost()           is False)
    assert(ion.symbol               == "Fe")
    assert(ion.total_mass           == 670.14)
    assert(ion.total_magnetization  == 12)
    assert(ion.total_charge_deficit == 168)
#end def test_custom_ion_species


def test_minimal_ion_species():
    """Test to make sure the defaults are populated correctly."""
    ion = IonSpecies(element="Fe", count=12)

    assert(ion.element              is Elements.Iron)
    assert(ion.label                == "Fe")
    assert(ion.count                == 12)
    assert(ion.formal_charge        == 0)
    assert(ion.magnetization        == 0)
    assert(ion.Zeff                 == Elements.Iron.atomic_number)
    assert(ion.pseudized()          is False)
    assert(ion.is_ghost()           is False)
    assert(ion.symbol               == Elements.Iron.symbol)
    assert(ion.total_mass           == 670.14)
    assert(ion.total_magnetization  == 0)
    assert(ion.total_charge_deficit == 312)
#end def test_minimal_ion_species


def test_ion_species_eq():
    ref_element               = Elements.Iron
    ref_label                 = "Fe1"
    ref_count                 = 12
    ref_formal_charge         = 0
    ref_magnetization         = 0
    ref_Zeff                  = Elements.Iron.atomic_number
    ref_is_pseudo             = False
    ref_is_ghost              = False
    ref_symbol                = Elements.Iron.symbol
    ref_total_mass            = 670.14
    ref_total_charge_deficit  = 312
    ref_total_spin            = 0

    ion1 = IonSpecies(element="Fe", label="Fe1", count=12)

    assert(ion1.element               is ref_element)
    assert(ion1.label                 == ref_label)
    assert(ion1.count                 == ref_count)
    assert(ion1.formal_charge         == ref_formal_charge)
    assert(ion1.magnetization         == ref_magnetization)
    assert(ion1.Zeff                  == ref_Zeff)
    assert(ion1.pseudized()           is ref_is_pseudo)
    assert(ion1.is_ghost()            is ref_is_ghost)
    assert(ion1.symbol                == ref_symbol)
    assert(ion1.total_mass            == ref_total_mass)
    assert(ion1.total_charge_deficit  == ref_total_charge_deficit)
    assert(ion1.total_magnetization   == ref_total_spin)

    ion2 = IonSpecies(element=Elements.Iron, label="Fe1", count=12)

    assert(ion2.element               is ref_element)
    assert(ion2.label                 == ref_label)
    assert(ion2.count                 == ref_count)
    assert(ion2.formal_charge         == ref_formal_charge)
    assert(ion2.magnetization         == ref_magnetization)
    assert(ion2.Zeff                  == ref_Zeff)
    assert(ion2.pseudized()           is ref_is_pseudo)
    assert(ion2.is_ghost()            is ref_is_ghost)
    assert(ion2.symbol                == ref_symbol)
    assert(ion2.total_mass            == ref_total_mass)
    assert(ion2.total_charge_deficit  == ref_total_charge_deficit)
    assert(ion2.total_magnetization   == ref_total_spin)

    assert(ion1 == ion2)
#end def test_ion_species_eq


def test_ion_species_repr():
    ref_repr = "IonSpecies(element=Fe, count=12, label=Fe, formal_charge=0, magnetization=0, Zeff=26)"
    ion = IonSpecies(element="Fe", count=12)
    assert(repr(ion) == ref_repr)
#end def test_ion_species_repr


def test_ion_species_hash():
    ion1 = IonSpecies(
        element       = Elements.Carbon,
        count         = 1,
        label         = "C2",
        formal_charge = 0,
        magnetization     = 0,
        Zeff          = 6,
        )
    ion2 = IonSpecies(
        element       = "C", # Use a slightly different element specifier to make sure it resolves.
        count         = 1,
        label         = "C2",
        formal_charge = 0,
        magnetization     = 0,
        Zeff          = 6,
        )
    assert(ion1 == ion2) # Make sure they're actually equal.
    assert(hash(ion1) == hash(ion2))
#end def test_ion_species_hash


def test_ion_species_from_structure():
    ref_ions = {
        "C1": IonSpecies(element=Elements.Carbon,   count=1, label="C1", formal_charge=0, magnetization=0, Zeff=6),
        "C2": IonSpecies(element=Elements.Carbon,   count=1, label="C2", formal_charge=0, magnetization=0, Zeff=6),
        "H" : IonSpecies(element=Elements.Hydrogen, count=5, label="H",  formal_charge=0, magnetization=0, Zeff=1),
        "N" : IonSpecies(element=Elements.Nitrogen, count=1, label="N",  formal_charge=0, magnetization=0, Zeff=7),
        "O" : IonSpecies(element=Elements.Oxygen,   count=2, label="O",  formal_charge=0, magnetization=0, Zeff=8),
        }

    structure = Structure(
        elem = ["N", "C1", "C2", "O", "O", "H", "H", "H", "H", "H"],
        pos = np.zeros((10,3)),
        )
    ions = IonSpecies.from_structure(
        structure=structure,
        )

    assert(ions == ref_ions)

    ref_ions = {
        "C1": IonSpecies(element=Elements.Carbon,   count=1, label="C1", formal_charge=0, magnetization=0,   Zeff=4),
        "C2": IonSpecies(element=Elements.Carbon,   count=1, label="C2", formal_charge=0, magnetization=0.5, Zeff=6),
        "H" : IonSpecies(element=Elements.Hydrogen, count=5, label="H",  formal_charge=0, magnetization=0.5, Zeff=1),
        "N" : IonSpecies(element=Elements.Nitrogen, count=1, label="N",  formal_charge=0, magnetization=1,   Zeff=5),
        "O" : IonSpecies(element=Elements.Oxygen,   count=2, label="O",  formal_charge=0, magnetization=0,   Zeff=6),
        }

    structure = Structure(
        elem = ["N", "C1", "C2", "O", "O", "H", "H", "H", "H", "H"],
        pos = np.zeros((10,3)),
        )
    ions = IonSpecies.from_structure(
        structure   = structure,
        elem_charge = dict(N=0, C1=0, C2=0,   O=0, H=0),
        elem_mag    = dict(N=1, C1=0, C2=0.5, O=0, H=0.5),
        elem_Zeff   = dict(N=5, C1=4, C2=6,   O=6, H=1),
        )

    assert(ions == ref_ions)
#end def test_ion_species_from_structure


def test_electrons_neutralize_to():
    ions = {
        "C1": IonSpecies(element=Elements.Carbon,   count=1, label="C1", formal_charge=0, magnetization=0,   Zeff=4),
        "C2": IonSpecies(element=Elements.Carbon,   count=1, label="C2", formal_charge=0, magnetization=0.5, Zeff=6),
        "H" : IonSpecies(element=Elements.Hydrogen, count=5, label="H",  formal_charge=0, magnetization=0.5, Zeff=1),
        "N" : IonSpecies(element=Elements.Nitrogen, count=1, label="N",  formal_charge=0, magnetization=1,   Zeff=5),
        "O" : IonSpecies(element=Elements.Oxygen,   count=2, label="O",  formal_charge=0, magnetization=0,   Zeff=6),
        }

    ref_electrons = Electrons(count=32, n_unpaired=0, spin_orbit=False)

    electrons = Electrons.neutralize_to(
        ions         = ions,
        total_charge = 0,
        n_unpaired   = 0,
        spin_orbit   = False,
    )
    assert(electrons == ref_electrons)
    assert(electrons.n_up   == 16)
    assert(electrons.n_down == 16)
    assert(isinstance(electrons.count, int))
    assert(isinstance(electrons.spin,  int))

    ref_electrons = Electrons(count=31, n_unpaired=1, spin_orbit=False)

    electrons = Electrons.neutralize_to(
        ions          = ions,
        total_charge  = 1,
        n_unpaired    = 1,
        spin_orbit    = False,
    )
    assert(electrons == ref_electrons)
    assert(electrons.n_up   == 16)
    assert(electrons.n_down == 15)
    assert(isinstance(electrons.count, int))
    assert(isinstance(electrons.spin,  int))

    # Test with integer-value floats for count and spin.
    # Should be turned into ints.
    ref_electrons = Electrons(count=32, n_unpaired=2, spin_orbit=False)

    electrons = Electrons.neutralize_to(
        ions          = ions,
        total_charge  = 0.0,
        n_unpaired    = 2.0,
        spin_orbit    = False,
    )
    assert(electrons == ref_electrons)
    assert(electrons.n_up   == 17)
    assert(electrons.n_down == 15)
    assert(isinstance(electrons.count, int))
    assert(isinstance(electrons.spin,  int))

    # See if it can make the right spin with floats
    ref_electrons = Electrons(count=31.5, n_unpaired=1.5, spin_orbit=False)

    electrons = Electrons.neutralize_to(
        ions          = ions,
        total_charge  = 0.5,
        n_unpaired    = None,
        spin_orbit    = False,
    )
    assert(electrons == ref_electrons)
    assert(electrons.n_up   == 16.5)
    assert(electrons.n_down == 15)
    assert(isinstance(electrons.count, float))
    assert(isinstance(electrons.spin,  float))
#end def test_electrons_neutralize_to


def test_molecular_system():
    ref_ions = {
        'C':  IonSpecies(element=Elements.Carbon,   count=2, label="C",  formal_charge= 0, magnetization=0,   Zeff=6),
        'H':  IonSpecies(element=Elements.Hydrogen, count=4, label="H",  formal_charge= 0, magnetization=0.5, Zeff=1),
        'N':  IonSpecies(element=Elements.Nitrogen, count=1, label="N",  formal_charge= 0, magnetization=1,   Zeff=7),
        'O1': IonSpecies(element=Elements.Oxygen,   count=1, label="O1", formal_charge= 0, magnetization=0,   Zeff=8),
        'O2': IonSpecies(element=Elements.Oxygen,   count=1, label="O2", formal_charge=-1, magnetization=0,   Zeff=8)
        }
    ref_electrons = Electrons(count=40, n_unpaired=0, spin_orbit=False)
    ref_structure_elem = ["N", "C", "C", "O1", "O2", "H", "H", "H", "H"]
    ref_structure_pos = np.array([
        [ 0.711045, 1.361274, 3.966292],
        [ 1.848745, 2.865874, 3.041292],
        [ 0.679145, 1.977474, 3.067692],
        [ 1.827545, 3.514674, 3.813792],
        [-0.580355, 2.805074, 3.070592],
        [-0.510755, 4.011174, 3.052592],
        [ 2.706445, 2.334474, 3.038892],
        [ 0.690245, 1.335874, 2.186592],
        [-1.779455, 2.202274, 3.093292],
        ], dtype=float)
    ref_structure_units = "A"

    # Glycinate
    structure = Structure(
        elem  = ref_structure_elem,
        pos   = ref_structure_pos,
        units = ref_structure_units,
        )
    system = PhysicalSystem(
        structure     = structure,
        total_charge  = None,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_charge   = dict(O2=-1),
        elem_mag      = dict(N=1, C=0, O1=0, O2= 0, H=0.5),
        )

    # Structure comparison
    assert(system.structure.elem.tolist() == ref_structure_elem)
    np.testing.assert_allclose(system.structure.pos, ref_structure_pos)
    assert(system.structure.units == ref_structure_units)

    # Attribute comparison
    assert(system.ions            == ref_ions)
    assert(system.electrons       == ref_electrons)
    assert(system.positrons       is None)
    assert(system.folded_system   is None)

    # Property comparison
    assert(system.ion_charge      == 39)
    assert(system.ion_mag         == 3.0)
    assert(system.electron_charge == -40)
    assert(system.electron_charge == system.electrons.total_charge)
    assert(system.net_charge      == -1)
    assert(system.net_charge      == system.ion_charge + system.electron_charge)


    # Test again, but use `total_charge` instead of `elem_charge`

    ref_ions = {
        'C':  IonSpecies(element=Elements.Carbon,   count=2, label="C",  formal_charge=0, magnetization=0,   Zeff=6),
        'H':  IonSpecies(element=Elements.Hydrogen, count=4, label="H",  formal_charge=0, magnetization=0.5, Zeff=1),
        'N':  IonSpecies(element=Elements.Nitrogen, count=1, label="N",  formal_charge=0, magnetization=1,   Zeff=7),
        'O1': IonSpecies(element=Elements.Oxygen,   count=1, label="O1", formal_charge=0, magnetization=0,   Zeff=8),
        'O2': IonSpecies(element=Elements.Oxygen,   count=1, label="O2", formal_charge=0, magnetization=0,   Zeff=8)
        }

    system = PhysicalSystem(
        structure     = structure,
        total_charge  = -1,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_mag      = dict(N=1, C=0, O1=0, O2= 0, H=0.5),
        )

    # Structure comparison
    assert(system.structure.elem.tolist() == ref_structure_elem)
    assert(system.structure.units == ref_structure_units)
    np.testing.assert_allclose(system.structure.pos, ref_structure_pos)

    # Attribute comparison
    assert(system.ions            == ref_ions)
    assert(system.electrons       == ref_electrons)
    assert(system.positrons       is None)
    assert(system.folded_system   is None)

    # Property comparison
    assert(system.ion_charge      == 39)
    assert(system.electron_charge == -40)
    assert(system.electron_charge == system.electrons.total_charge)
    assert(system.net_charge      == -1)
    assert(system.net_charge      == system.ion_charge + system.electron_charge)

    # Check some simple functions
    assert(system.num_ions()      == 9)
    assert(system.num_ion_types() == 5)
#end def test_molecular_system


def test_pretiled_system():
    refs = get_LaAlO3_references()
    ref_folded_elem      = refs["ref_folded_elem"]
    ref_tiled_elem       = refs["ref_tiled_elem"]
    ref_folded_axes      = refs["ref_folded_axes"]
    ref_tiled_axes       = refs["ref_tiled_axes"]
    ref_folded_pos       = refs["ref_folded_pos"]
    ref_tiled_pos        = refs["ref_tiled_pos"]
    ref_folded_ions      = refs["ref_folded_ions"]
    ref_tiled_ions       = refs["ref_tiled_ions"]
    ref_folded_electrons = refs["ref_folded_electrons"]
    ref_tiled_electrons  = refs["ref_tiled_electrons"]
    ref_folded_positrons = refs["ref_folded_positrons"]
    ref_tiled_positrons  = refs["ref_tiled_positrons"]

    folded_structure = Structure(
        axes=ref_folded_axes,
        elem=ref_folded_elem,
        pos=ref_folded_pos,
        units="A",
        )

    # 8x everything
    tiled_structure = folded_structure.tile(2,2,2)

    # Sanity check on tiling.
    assert(tiled_structure.elem.tolist() == ref_tiled_elem)
    np.testing.assert_allclose(tiled_structure.axes, ref_tiled_axes)
    np.testing.assert_allclose(tiled_structure.pos,  ref_tiled_pos)

    assert(tiled_structure.folded_structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(tiled_structure.folded_structure.axes, ref_folded_axes)
    np.testing.assert_allclose(tiled_structure.folded_structure.pos,  ref_folded_pos)

    tiled_ps = PhysicalSystem(
        structure     = tiled_structure,
        total_charge  = 0,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_mag      = dict(La=3.5, Al=2.5, O1=0, O2=0, O3=0),
        elem_Zeff     = dict(La=11,  Al=3,   O1=6, O2=6, O3=6),
        positrons     = ref_tiled_positrons,
        )

    # Make sure the folded system exists, and its structure is the folded structure
    assert(tiled_ps.structure.elem.tolist() == ref_tiled_elem)
    np.testing.assert_allclose(tiled_ps.structure.axes, ref_tiled_axes)
    np.testing.assert_allclose(tiled_ps.structure.pos,  ref_tiled_pos)

    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.axes, ref_folded_axes)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  ref_folded_pos)

    assert(tiled_ps.ions      == ref_tiled_ions)
    assert(tiled_ps.electrons == ref_tiled_electrons)
    assert(tiled_ps.positrons == ref_tiled_positrons)

    assert(tiled_ps.folded_system.ions      == ref_folded_ions)
    assert(tiled_ps.folded_system.electrons == ref_folded_electrons)
    assert(tiled_ps.folded_system.positrons == ref_folded_positrons)

    assert(tiled_ps.net_charge      ==   64) # 64 positrons
    assert(tiled_ps.ion_charge      ==  256)
    assert(tiled_ps.ion_mag         ==   48)
    assert(tiled_ps.electron_charge == -256)
    assert(tiled_ps.electron_spin   ==    0)
    assert(tiled_ps.positron_charge ==   64)
    assert(tiled_ps.positron_spin   ==    0)

    assert(tiled_ps.folded_system.net_charge      ==   8) # 8 positrons
    assert(tiled_ps.folded_system.ion_charge      ==  32)
    assert(tiled_ps.folded_system.ion_mag         ==   6)
    assert(tiled_ps.folded_system.electron_charge == -32)
    assert(tiled_ps.folded_system.electron_spin   ==   0)
    assert(tiled_ps.folded_system.positron_charge ==   8)
    assert(tiled_ps.folded_system.positron_spin   ==   0)

    assert(tiled_ps.is_magnetic())
    assert(tiled_ps.folded_system.is_magnetic())

    # Now that we've checked everything manually, the last
    # check is to make sure the built-in checks work
    assert(tiled_ps.check_folded_system())
    assert(tiled_ps.check_consistent())
    assert(tiled_ps.is_valid())
    assert(tiled_ps.has_folded())
    assert(tiled_ps.get_smallest() is tiled_ps.folded_system)
    tiled_ps.remove_folded_system()
    assert(tiled_ps.folded_system is None)
#end def test_pretiled_system


def test_tile():
    refs = get_LaAlO3_references()
    ref_folded_elem      = refs["ref_folded_elem"]
    ref_tiled_elem       = refs["ref_tiled_elem"]
    ref_folded_axes      = refs["ref_folded_axes"]
    ref_tiled_axes       = refs["ref_tiled_axes"]
    ref_folded_pos       = refs["ref_folded_pos"]
    ref_tiled_pos        = refs["ref_tiled_pos"]
    ref_folded_ions      = refs["ref_folded_ions"]
    ref_tiled_ions       = refs["ref_tiled_ions"]
    ref_folded_electrons = refs["ref_folded_electrons"]
    ref_tiled_electrons  = refs["ref_tiled_electrons"]
    ref_folded_positrons = refs["ref_folded_positrons"]
    ref_tiled_positrons  = refs["ref_tiled_positrons"]

    folded_structure = Structure(
        axes=ref_folded_axes,
        elem=ref_folded_elem,
        pos=ref_folded_pos,
        units="A",
        )

    assert(folded_structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(folded_structure.axes, ref_folded_axes)
    np.testing.assert_allclose(folded_structure.pos,  ref_folded_pos)

    folded_ps = PhysicalSystem(
        structure     = folded_structure,
        total_charge  = 0,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_mag      = dict(La=3.5, Al=2.5, O1=0, O2=0, O3=0),
        elem_Zeff     = dict(La=11,  Al=3,   O1=6, O2=6, O3=6),
        positrons     = ref_folded_positrons,
        )

    assert(folded_ps.folded_system is None)
    assert(folded_ps.ions          == ref_folded_ions)
    assert(folded_ps.electrons     == ref_folded_electrons)
    assert(folded_ps.positrons     == ref_folded_positrons)

    tiled_ps = folded_ps.tile(2,2,2)

    assert(tiled_ps.structure.elem.tolist() == ref_tiled_elem)
    np.testing.assert_allclose(tiled_ps.structure.axes, ref_tiled_axes)
    np.testing.assert_allclose(tiled_ps.structure.pos,  ref_tiled_pos)

    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.axes, ref_folded_axes)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  ref_folded_pos)

    assert(tiled_ps.ions      == ref_tiled_ions)
    assert(tiled_ps.electrons == ref_tiled_electrons)
    assert(tiled_ps.positrons == ref_tiled_positrons)

    assert(tiled_ps.folded_system.ions      == ref_folded_ions)
    assert(tiled_ps.folded_system.electrons == ref_folded_electrons)
    assert(tiled_ps.folded_system.positrons == ref_folded_positrons)
#end def test_tile


def test_init_with_path(tmp_path):
    refs = get_d2_d8_references()
    d8_ref_axes   = refs["d8_ref_axes"]
    d8_ref_pos    = refs["d8_ref_pos"]
    d8_ref_elem   = refs["d8_ref_elem"]
    d8_ref_units  = refs["d8_ref_units"]
    d8_ref_bconds = refs["d8_ref_bconds"]
    d8_structure  = refs["d8_structure"]

    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    d8_path = tmp_path / 'diamond8.xsf'
    d8_structure.write(d8_path)

    system = PhysicalSystem(
        structure   = d8_path,
        elem_Zeff   = dict(C=4),
        )

    assert(system.net_charge    == 0)
    assert(system.electron_spin == 0)
    assert(system.pseudized)
    np.testing.assert_allclose(system.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(system.structure.pos,  d8_ref_pos)
    assert(system.structure.elem.tolist() == d8_ref_elem)
    assert(system.structure.units         == d8_ref_units)
    assert(tuple(system.structure.bconds) == d8_ref_bconds)
    assert(system.ions      == d8_ref_ions)
    assert(system.electrons == d8_ref_electrons)
#end def test_init_with_path


def test_change_units():
    refs = get_LaAlO3_references()
    ref_folded_elem      = refs["ref_folded_elem"]
    ref_tiled_elem       = refs["ref_tiled_elem"]
    ref_folded_axes      = refs["ref_folded_axes"]
    ref_tiled_axes       = refs["ref_tiled_axes"]
    ref_folded_pos       = refs["ref_folded_pos"]
    ref_tiled_pos        = refs["ref_tiled_pos"]
    ref_folded_ions      = refs["ref_folded_ions"]
    ref_tiled_ions       = refs["ref_tiled_ions"]
    ref_folded_electrons = refs["ref_folded_electrons"]
    ref_tiled_electrons  = refs["ref_tiled_electrons"]
    ref_folded_positrons = refs["ref_folded_positrons"]
    ref_tiled_positrons  = refs["ref_tiled_positrons"]

    # This is already tested in `test_unit_converter.py`
    ref_folded_axes_bohr = convert(ref_folded_axes, "A", "B")
    ref_tiled_axes_bohr  = convert(ref_tiled_axes, "A", "B")
    ref_folded_pos_bohr  = convert(ref_folded_pos, "A", "B")
    ref_tiled_pos_bohr   = convert(ref_tiled_pos, "A", "B")

    folded_structure = Structure(
        axes=ref_folded_axes,
        elem=ref_folded_elem,
        pos=ref_folded_pos,
        units="A",
        )

    # 8x everything
    tiled_structure = folded_structure.tile(2,2,2)
    tiled_ps = PhysicalSystem(
        structure     = tiled_structure,
        total_charge  = 0,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_mag      = dict(La=3.5, Al=2.5, O1=0, O2=0, O3=0),
        elem_Zeff     = dict(La=11,  Al=3,   O1=6, O2=6, O3=6),
        positrons     = ref_tiled_positrons,
        )

    # Test everything beforehand, make sure it's in working order.
    assert(tiled_ps.structure.elem.tolist() == ref_tiled_elem)
    np.testing.assert_allclose(tiled_ps.structure.axes, ref_tiled_axes)
    np.testing.assert_allclose(tiled_ps.structure.pos,  ref_tiled_pos)

    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.axes, ref_folded_axes)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  ref_folded_pos)

    assert(tiled_ps.ions      == ref_tiled_ions)
    assert(tiled_ps.electrons == ref_tiled_electrons)
    assert(tiled_ps.positrons == ref_tiled_positrons)

    assert(tiled_ps.folded_system.ions      == ref_folded_ions)
    assert(tiled_ps.folded_system.electrons == ref_folded_electrons)
    assert(tiled_ps.folded_system.positrons == ref_folded_positrons)

    # Change units
    tiled_ps.change_units("B")

    # Make sure nothing else was touched.
    assert(tiled_ps.structure.elem.tolist() == ref_tiled_elem)
    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem)

    assert(tiled_ps.ions      == ref_tiled_ions)
    assert(tiled_ps.electrons == ref_tiled_electrons)
    assert(tiled_ps.positrons == ref_tiled_positrons)

    assert(tiled_ps.folded_system.ions      == ref_folded_ions)
    assert(tiled_ps.folded_system.electrons == ref_folded_electrons)
    assert(tiled_ps.folded_system.positrons == ref_folded_positrons)

    # These are the only things that should have changed
    np.testing.assert_allclose(tiled_ps.structure.axes, ref_tiled_axes_bohr)
    np.testing.assert_allclose(tiled_ps.structure.pos,  ref_tiled_pos_bohr)

    np.testing.assert_allclose(tiled_ps.folded_system.structure.axes, ref_folded_axes_bohr)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  ref_folded_pos_bohr)
#end def test_change_units


def test_rename():
    refs = get_LaAlO3_references()
    ref_folded_elem      = refs["ref_folded_elem"]
    ref_tiled_elem       = refs["ref_tiled_elem"]
    ref_folded_ions      = refs["ref_folded_ions"]
    ref_tiled_ions       = refs["ref_tiled_ions"]

    # "La" -> "La2"
    # "Al" -> "Al2"
    ref_folded_elem_new = [
        "La2", "Al2", "O1", "O2", "O3",
        ]
    ref_tiled_elem_new  = [
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        "La2", "Al2", "O1", "O2", "O3",
        ]
    ref_folded_ions_new = {
        'Al2': IonSpecies(element=Elements.Aluminum,  count=1, label="Al2", formal_charge=0, magnetization=2.5, Zeff=3),
        'La2': IonSpecies(element=Elements.Lanthanum, count=1, label="La2", formal_charge=0, magnetization=3.5, Zeff=11),
        'O1': IonSpecies(element=Elements.Oxygen,    count=1, label="O1", formal_charge=0, magnetization=0,   Zeff=6),
        'O2': IonSpecies(element=Elements.Oxygen,    count=1, label="O2", formal_charge=0, magnetization=0,   Zeff=6),
        'O3': IonSpecies(element=Elements.Oxygen,    count=1, label="O3", formal_charge=0, magnetization=0,   Zeff=6),
        }
    ref_tiled_ions_new = {
        'Al2': IonSpecies(element=Elements.Aluminum,  count=8, label="Al2", formal_charge=0, magnetization=2.5, Zeff=3),
        'La2': IonSpecies(element=Elements.Lanthanum, count=8, label="La2", formal_charge=0, magnetization=3.5, Zeff=11),
        'O1': IonSpecies(element=Elements.Oxygen,    count=8, label="O1", formal_charge=0, magnetization=0,   Zeff=6),
        'O2': IonSpecies(element=Elements.Oxygen,    count=8, label="O2", formal_charge=0, magnetization=0,   Zeff=6),
        'O3': IonSpecies(element=Elements.Oxygen,    count=8, label="O3", formal_charge=0, magnetization=0,   Zeff=6),
        }

    folded_structure = Structure(
        axes=refs["ref_folded_axes"],
        elem=ref_folded_elem,
        pos=refs["ref_folded_pos"],
        units="A",
        )

    folded_ps = PhysicalSystem(
        structure     = folded_structure,
        total_charge  = 0,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_mag      = dict(La=3.5, Al=2.5, O1=0, O2=0, O3=0),
        elem_Zeff     = dict(La=11,  Al=3,   O1=6, O2=6, O3=6),
        positrons     = refs["ref_folded_positrons"],
        )

    tiled_ps = folded_ps.tile(2,2,2)

    assert(tiled_ps.structure.elem.tolist()               == ref_tiled_elem)
    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem)

    assert(tiled_ps.ions               == ref_tiled_ions)
    assert(tiled_ps.folded_system.ions == ref_folded_ions)

    tiled_ps.rename(
        La="La2",
        Al="Al2",
    )

    assert(tiled_ps.structure.elem.tolist()               == ref_tiled_elem_new)
    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem_new)

    assert(tiled_ps.ions               == ref_tiled_ions_new)
    assert(tiled_ps.folded_system.ions == ref_folded_ions_new)
#end def test_rename


def test_group_atoms():
    folded_ungrouped_elem = ["B1", "N", "O", "He", "Al"]
    tiled_ungrouped_elem  = ["B1", "N", "O", "He", "Al", "B1", "N", "O", "He", "Al"]
    folded_grouped_elem   = ["Al", "B1", "He", "N", "O"]
    tiled_grouped_elem    = ["Al", "Al", "B1", "B1", "He", "He", "N", "N", "O", "O"]

    folded_ungrouped_pos = np.array([
        [1, 1, 1], # B1 = 1
        [3, 3, 3], # N  = 3
        [4, 4, 4], # O  = 4
        [2, 2, 2], # He = 2
        [0, 0, 0], # Al = 0
    ], dtype=float)
    tiled_ungrouped_pos = np.array([
        [ 1,  1,  1], # B1 = 1
        [ 3,  3,  3], # N  = 3
        [ 4,  4,  4], # O  = 4
        [ 2,  2,  2], # He = 2
        [ 0,  0,  0], # Al = 0
        [11,  1,  1], # B1 = 1
        [13,  3,  3], # N  = 3
        [14,  4,  4], # O  = 4
        [12,  2,  2], # He = 2
        [10,  0,  0], # Al = 0
    ], dtype=float)
    folded_grouped_pos = np.array([
        [0, 0, 0], # Al = 0
        [1, 1, 1], # B1 = 1
        [2, 2, 2], # He = 2
        [3, 3, 3], # N  = 3
        [4, 4, 4], # O  = 4
    ], dtype=float)
    tiled_grouped_pos = np.array([
        [ 0, 0, 0], # Al = 0
        [10, 0, 0], # Al = 0
        [ 1, 1, 1], # B1 = 1
        [11, 1, 1], # B1 = 1
        [ 2, 2, 2], # He = 2
        [12, 2, 2], # He = 2
        [ 3, 3, 3], # N  = 3
        [13, 3, 3], # N  = 3
        [ 4, 4, 4], # O  = 4
        [14, 4, 4], # O  = 4
    ], dtype=float)

    folded_structure = Structure(
        axes=np.array([
            [10,  0,  0],
            [ 0, 10,  0],
            [ 0,  0, 10],
        ], dtype=float),
        elem=folded_ungrouped_elem,
        pos=folded_ungrouped_pos,
        units="A",
        )

    assert(folded_structure.elem.tolist() == folded_ungrouped_elem)
    np.testing.assert_allclose(folded_structure.pos, folded_ungrouped_pos)

    folded_ps = PhysicalSystem(structure=folded_structure)

    tiled_ps = folded_ps.tile(2,1,1)

    assert(tiled_ps.structure.elem.tolist() == tiled_ungrouped_elem)
    np.testing.assert_allclose(tiled_ps.structure.pos,  tiled_ungrouped_pos)

    assert(tiled_ps.folded_system.structure.elem.tolist() == folded_ungrouped_elem)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  folded_ungrouped_pos)

    tiled_ps.group_atoms()

    assert(tiled_ps.structure.elem.tolist() == tiled_grouped_elem)
    np.testing.assert_allclose(tiled_ps.structure.pos,  tiled_grouped_pos)

    assert(tiled_ps.folded_system.structure.elem.tolist() == folded_grouped_elem)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  folded_grouped_pos)
#end def test_group_atoms


def test_generate_struct():
    refs = get_d2_d8_references()
    d2_ref_axes   = refs["d2_ref_axes"]
    d2_ref_pos    = refs["d2_ref_pos"]
    d2_ref_elem   = refs["d2_ref_elem"]
    d2_ref_units  = refs["d2_ref_units"]
    d2_ref_bconds = refs["d2_ref_bconds"]
    d2_structure  = refs["d2_structure"]
    d8_ref_axes   = refs["d8_ref_axes"]
    d8_ref_pos    = refs["d8_ref_pos"]
    d8_ref_elem   = refs["d8_ref_elem"]
    d8_ref_units  = refs["d8_ref_units"]
    d8_ref_bconds = refs["d8_ref_bconds"]
    d8_structure  = refs["d8_structure"]

    d2_ref_ions      = refs["d2_ref_ions"]
    d2_ref_electrons = refs["d2_ref_electrons"]
    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    struct_notile = generate_physical_system(
        structure = d8_structure,
        C         = 4,
        )

    struct_tile = generate_physical_system(
        structure = d2_structure,
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    actual = struct_notile
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    actual        = struct_tile
    actual_folded = struct_tile.folded_system
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    assert(actual_folded.net_charge    == 0)
    assert(actual_folded.electron_spin == 0)
    assert(actual_folded.pseudized)
    np.testing.assert_allclose(actual_folded.structure.axes, d2_ref_axes)
    np.testing.assert_allclose(actual_folded.structure.pos,  d2_ref_pos)
    assert(actual_folded.structure.elem.tolist()  == d2_ref_elem)
    assert(actual_folded.structure.units          == d2_ref_units)
    assert(tuple(actual_folded.structure.bconds)  == d2_ref_bconds)
    assert(actual_folded.ions      == d2_ref_ions)
    assert(actual_folded.electrons == d2_ref_electrons)

    assert(struct_notile.is_valid())
    assert(struct_tile.is_valid())
    assert(not struct_notile.has_folded())
    assert(struct_tile.has_folded())
#end def test_generate_struct


def test_generate_direct():
    refs = get_d2_d8_references()
    d2_ref_axes   = refs["d2_ref_axes"]
    d2_ref_pos    = refs["d2_ref_pos"]
    d2_ref_elem   = refs["d2_ref_elem"]
    d2_ref_units  = refs["d2_ref_units"]
    d2_ref_bconds = refs["d2_ref_bconds"]
    d8_ref_axes   = refs["d8_ref_axes"]
    d8_ref_pos    = refs["d8_ref_pos"]
    d8_ref_elem   = refs["d8_ref_elem"]
    d8_ref_units  = refs["d8_ref_units"]
    d8_ref_bconds = refs["d8_ref_bconds"]

    d2_ref_ions      = refs["d2_ref_ions"]
    d2_ref_electrons = refs["d2_ref_electrons"]
    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    direct_notile = generate_physical_system(
        units = 'A',
        axes  = [[3.57, 0.00, 0.00],
                 [0.00, 3.57, 0.00],
                 [0.00, 0.00, 3.57]],
        elem  = 8*['C'],
        posu  = [[0.00, 0.00, 0.00],
                 [0.25, 0.25, 0.25],
                 [0.50, 0.50, 0.00],
                 [0.75, 0.75, 0.25],
                 [0.00, 0.50, 0.50],
                 [0.25, 0.75, 0.75],
                 [0.50, 0.00, 0.50],
                 [0.75, 0.25, 0.75]],
        C     = 4,
        )

    direct_tile = generate_physical_system(
        units  = 'A',
        axes   = [[1.785, 1.785, 0.   ],
                  [0.   , 1.785, 1.785],
                  [1.785, 0.   , 1.785]],
        elem   = 2*['C'],
        posu   = [[0.00, 0.00, 0.00],
                  [0.25, 0.25, 0.25]],
        tiling = [[ 1, -1,  1],
                  [ 1,  1, -1],
                  [-1,  1,  1]],
        C      = 4,
        )

    actual = direct_notile
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    actual        = direct_tile
    actual_folded = direct_tile.folded_system
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    assert(actual_folded.net_charge    == 0)
    assert(actual_folded.electron_spin == 0)
    assert(actual_folded.pseudized)
    np.testing.assert_allclose(actual_folded.structure.axes, d2_ref_axes)
    np.testing.assert_allclose(actual_folded.structure.pos,  d2_ref_pos)
    assert(actual_folded.structure.elem.tolist()  == d2_ref_elem)
    assert(actual_folded.structure.units          == d2_ref_units)
    assert(tuple(actual_folded.structure.bconds)  == d2_ref_bconds)
    assert(actual_folded.ions      == d2_ref_ions)
    assert(actual_folded.electrons == d2_ref_electrons)

    assert(direct_notile.is_valid())
    assert(direct_tile.is_valid())
    assert(not direct_notile.has_folded())
    assert(direct_tile.has_folded())
#end def test_generate_direct


def test_generate_with_path(tmp_path):
    refs = get_d2_d8_references()
    d2_ref_axes   = refs["d2_ref_axes"]
    d2_ref_pos    = refs["d2_ref_pos"]
    d2_ref_elem   = refs["d2_ref_elem"]
    d2_ref_units  = refs["d2_ref_units"]
    d2_ref_bconds = refs["d2_ref_bconds"]
    d2_structure  = refs["d2_structure"]
    d8_ref_axes   = refs["d8_ref_axes"]
    d8_ref_pos    = refs["d8_ref_pos"]
    d8_ref_elem   = refs["d8_ref_elem"]
    d8_ref_units  = refs["d8_ref_units"]
    d8_ref_bconds = refs["d8_ref_bconds"]
    d8_structure  = refs["d8_structure"]

    d2_ref_ions      = refs["d2_ref_ions"]
    d2_ref_electrons = refs["d2_ref_electrons"]
    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    d2_path = tmp_path / 'diamond2.xsf'
    d2_structure.write(d2_path)

    d8_path = tmp_path / 'diamond8.xsf'
    d8_structure.write(d8_path)

    read_notile = generate_physical_system(
        structure = d8_path,
        C         = 4,
        )

    read_tile = generate_physical_system(
        structure = d2_path,
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    actual = read_notile
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    actual        = read_tile
    actual_folded = read_tile.folded_system
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    assert(actual_folded.net_charge    == 0)
    assert(actual_folded.electron_spin == 0)
    assert(actual_folded.pseudized)
    np.testing.assert_allclose(actual_folded.structure.axes, d2_ref_axes)
    np.testing.assert_allclose(actual_folded.structure.pos,  d2_ref_pos)
    assert(actual_folded.structure.elem.tolist()  == d2_ref_elem)
    assert(actual_folded.structure.units          == d2_ref_units)
    assert(tuple(actual_folded.structure.bconds)  == d2_ref_bconds)
    assert(actual_folded.ions      == d2_ref_ions)
    assert(actual_folded.electrons == d2_ref_electrons)

    assert(read_notile.is_valid())
    assert(read_tile.is_valid())
    assert(not read_notile.has_folded())
    assert(read_tile.has_folded())
#end def test_generate_with_path


def test_generate_gen():
    refs = get_d2_d8_references()
    d2_ref_axes      = refs["d2_ref_axes"]
    d2_ref_pos       = refs["d2_ref_pos"]
    d2_ref_elem      = refs["d2_ref_elem"]
    d2_ref_units     = refs["d2_ref_units"]
    d2_ref_bconds    = refs["d2_ref_bconds"]
    d8_ref_axes      = refs["d8_ref_axes"]
    d8_ref_pos       = refs["d8_ref_pos"]
    d8_ref_pos_known = refs["d8_ref_pos_known"] # Use known_crystals version
    d8_ref_elem      = refs["d8_ref_elem"]
    d8_ref_units     = refs["d8_ref_units"]
    d8_ref_bconds    = refs["d8_ref_bconds"]

    d2_ref_ions      = refs["d2_ref_ions"]
    d2_ref_electrons = refs["d2_ref_electrons"]
    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    gen_notile = generate_physical_system(
        lattice   = 'cubic',        # cubic tetragonal orthorhombic rhombohedral
                                    # hexagonal triclinic monoclinic
        cell      = 'conventional', # primitive, conventional
        centering = 'F',            # P A B C I R F
        constants = 3.57,           # a,b,c,alpha,beta,gamma
        units     = 'A',            # A or B
        atoms     = 'C',            # species in primitive cell
        basis     = [[0,0,0],       # basis vectors (optional)
                     [.25,.25,.25]],
        C         = 4,
        )

    gen_tile = generate_physical_system(
        lattice   = 'cubic',        # cubic tetragonal orthorhombic rhombohedral
                                    # hexagonal triclinic monoclinic
        cell      = 'primitive',    # primitive, conventional
        centering = 'F',            # P A B C I R F
        constants = 3.57,           # a,b,c,alpha,beta,gamma
        units     = 'A',            # A or B
        atoms     = 'C',            # species in primitive cell
        basis     = [[0,0,0],       # basis vectors (optional)
                     [.25,.25,.25]],
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    actual = gen_notile
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes,      atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos_known, atol=1e-12)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    actual        = gen_tile
    actual_folded = gen_tile.folded_system
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    assert(actual_folded.net_charge    == 0)
    assert(actual_folded.electron_spin == 0)
    assert(actual_folded.pseudized)
    np.testing.assert_allclose(actual_folded.structure.axes, d2_ref_axes)
    np.testing.assert_allclose(actual_folded.structure.pos,  d2_ref_pos)
    assert(actual_folded.structure.elem.tolist()  == d2_ref_elem)
    assert(actual_folded.structure.units          == d2_ref_units)
    assert(tuple(actual_folded.structure.bconds)  == d2_ref_bconds)
    assert(actual_folded.ions      == d2_ref_ions)
    assert(actual_folded.electrons == d2_ref_electrons)

    assert(gen_notile.is_valid())
    assert(gen_tile.is_valid())
    assert(not gen_notile.has_folded())
    assert(gen_tile.has_folded())
#end def test_generate_gen


def test_generate_lookup():
    refs = get_d2_d8_references()
    d2_ref_axes      = refs["d2_ref_axes"]
    d2_ref_pos       = refs["d2_ref_pos"]
    d2_ref_elem      = refs["d2_ref_elem"]
    d2_ref_units     = refs["d2_ref_units"]
    d2_ref_bconds    = refs["d2_ref_bconds"]
    d8_ref_axes      = refs["d8_ref_axes"]
    d8_ref_pos       = refs["d8_ref_pos"]
    d8_ref_pos_known = refs["d8_ref_pos_known"] # Use known_crystals version
    d8_ref_elem      = refs["d8_ref_elem"]
    d8_ref_units     = refs["d8_ref_units"]
    d8_ref_bconds    = refs["d8_ref_bconds"]

    d2_ref_ions      = refs["d2_ref_ions"]
    d2_ref_electrons = refs["d2_ref_electrons"]
    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    lookup_notile = generate_physical_system(
        structure = 'diamond',
        cell      = 'conv',
        C         = 4,
        )

    lookup_tile = generate_physical_system(
        structure = 'diamond',
        cell      = 'prim',
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    actual = lookup_notile
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes,      atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos_known, atol=1e-12)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    actual        = lookup_tile
    actual_folded = lookup_tile.folded_system
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    assert(actual_folded.net_charge    == 0)
    assert(actual_folded.electron_spin == 0)
    assert(actual_folded.pseudized)
    np.testing.assert_allclose(actual_folded.structure.axes, d2_ref_axes)
    np.testing.assert_allclose(actual_folded.structure.pos,  d2_ref_pos)
    assert(actual_folded.structure.elem.tolist()  == d2_ref_elem)
    assert(actual_folded.structure.units          == d2_ref_units)
    assert(tuple(actual_folded.structure.bconds)  == d2_ref_bconds)
    assert(actual_folded.ions      == d2_ref_ions)
    assert(actual_folded.electrons == d2_ref_electrons)

    assert(lookup_notile.is_valid())
    assert(lookup_tile.is_valid())
    assert(not lookup_notile.has_folded())
    assert(lookup_tile.has_folded())
#end def test_generate_lookup


@pytest.mark.skip
def test_save_load(tmp_path):
    refs = get_d2_d8_references()
    d2_ref_axes   = refs["d2_ref_axes"]
    d2_ref_pos    = refs["d2_ref_pos"]
    d2_ref_elem   = refs["d2_ref_elem"]
    d2_ref_units  = refs["d2_ref_units"]
    d2_ref_bconds = refs["d2_ref_bconds"]
    d8_ref_axes   = refs["d8_ref_axes"]
    d8_ref_pos    = refs["d8_ref_pos"]
    d8_ref_elem   = refs["d8_ref_elem"]
    d8_ref_units  = refs["d8_ref_units"]
    d8_ref_bconds = refs["d8_ref_bconds"]

    d2_ref_ions      = refs["d2_ref_ions"]
    d2_ref_electrons = refs["d2_ref_electrons"]
    d8_ref_ions      = refs["d8_ref_ions"]
    d8_ref_electrons = refs["d8_ref_electrons"]

    lookup_notile = generate_physical_system(
        structure = 'diamond',
        cell      = 'conv',
        C         = 4,
        )

    lookup_tile = generate_physical_system(
        structure = 'diamond',
        cell      = 'prim',
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    actual = lookup_notile
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    actual        = lookup_tile
    actual_folded = lookup_tile.folded_system
    assert(actual.net_charge    == 0)
    assert(actual.electron_spin == 0)
    assert(actual.pseudized)
    np.testing.assert_allclose(actual.structure.axes, d8_ref_axes, atol=1e-12)
    np.testing.assert_allclose(actual.structure.pos,  d8_ref_pos)
    assert(actual.structure.elem.tolist()  == d8_ref_elem)
    assert(actual.structure.units          == d8_ref_units)
    assert(tuple(actual.structure.bconds)  == d8_ref_bconds)
    assert(actual.ions      == d8_ref_ions)
    assert(actual.electrons == d8_ref_electrons)

    assert(actual_folded.net_charge    == 0)
    assert(actual_folded.electron_spin == 0)
    assert(actual_folded.pseudized)
    np.testing.assert_allclose(actual_folded.structure.axes, d2_ref_axes)
    np.testing.assert_allclose(actual_folded.structure.pos,  d2_ref_pos)
    assert(actual_folded.structure.elem.tolist()  == d2_ref_elem)
    assert(actual_folded.structure.units          == d2_ref_units)
    assert(tuple(actual_folded.structure.bconds)  == d2_ref_bconds)
    assert(actual.ions      == d2_ref_ions)
    assert(actual.electrons == d2_ref_electrons)

    assert(lookup_notile.is_valid())
    assert(lookup_tile.is_valid())
    assert(not lookup_notile.has_folded())
    assert(lookup_tile.has_folded())
    # test load
    for i,sys in enumerate(systems):
        path = tmp_path / 'system_{}'.format(i)
        sys.save(path)
        sys2 = PhysicalSystem()
        sys2.load(path)
        assert(sys2.is_valid())
        assert(system_same(sys2,sys,tiled=sys.has_folded()))
#end def test_save_load


def test_kf_rpa():
    alat = 3.3521298178767225

    structure = Structure(
        axes = np.eye(3) * alat,
        elem = ["H", "H", "H", "H"],
        pos  = np.array([
            [   0.0,    0.0,    0.0],
            [alat/2,    0.0,    0.0],
            [   0.0, alat/2,    0.0],
            [   0.0,    0.0, alat/2],
            ], dtype=float),
        units = "B",
        )
    ps = PhysicalSystem(
        structure     = structure,
        total_charge  = 1,
        unpaired_elns = 1,
        elem_Zeff     = {"H": 1},
    )
    kfs = ps.kf_rpa()
    np.testing.assert_allclose(kfs[0], 1.465,           atol=1e-3)
    np.testing.assert_allclose(kfs[1], 1.465/2**(1./3), atol=1e-3)
#end def test_kf_rpa


def test_ae_pp_species():
    # Glycinate
    structure = Structure(
        elem  = ["N", "C", "C", "O1", "O2", "H", "H", "H", "H"],
        pos   = np.empty((9, 3), dtype=float),
        units = "A",
        )
    system = PhysicalSystem(
        structure     = structure,
        total_charge  = None,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_Zeff     = dict(N=5, O1=6),
        )
    
    ae_species, pp_species = system.ae_pp_species()
    assert(ae_species == {"C", "O2", "H"})
    assert(pp_species == {"N", "O1"})
#end def test_ae_pp_species


def test_large_Zeff_elem():
        # Glycinate
    structure = Structure(
        elem  = ["N", "C", "C", "O1", "O2", "H", "H", "H", "H"],
        pos   = np.empty((9, 3), dtype=float),
        units = "A",
        )

    system = PhysicalSystem(
        structure     = structure,
        total_charge  = None,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_Zeff     = dict(N=5, O1=6),
        )

    large_Zeff_elem = system.large_Zeff_elem(Zmin=5)
    assert(set(large_Zeff_elem) == {"C", "O1", "O2"})
#end def test_large_Zeff_elem


def test_copy():
    refs = get_LaAlO3_references()
    ref_folded_elem      = refs["ref_folded_elem"]
    ref_tiled_elem       = refs["ref_tiled_elem"]
    ref_folded_axes      = refs["ref_folded_axes"]
    ref_tiled_axes       = refs["ref_tiled_axes"]
    ref_folded_pos       = refs["ref_folded_pos"]
    ref_tiled_pos        = refs["ref_tiled_pos"]
    ref_folded_ions      = refs["ref_folded_ions"]
    ref_tiled_ions       = refs["ref_tiled_ions"]
    ref_folded_electrons = refs["ref_folded_electrons"]
    ref_tiled_electrons  = refs["ref_tiled_electrons"]
    ref_folded_positrons = refs["ref_folded_positrons"]
    ref_tiled_positrons  = refs["ref_tiled_positrons"]

    tiled_structure = Structure(
        axes=ref_folded_axes,
        elem=ref_folded_elem,
        pos=ref_folded_pos,
        units="A",
        tiling=(2,2,2),
        )

    tiled_ps = PhysicalSystem(
        structure     = tiled_structure,
        total_charge  = 0,
        unpaired_elns = 0,
        spin_orbit    = False,
        elem_mag      = dict(La=3.5, Al=2.5, O1=0, O2=0, O3=0),
        elem_Zeff     = dict(La=11,  Al=3,   O1=6, O2=6, O3=6),
        positrons     = ref_tiled_positrons,
        )

    # Make sure the folded system exists, and its structure
    # is the folded structure. We use `is` because it compares
    # memory addresses, which is what we want to be the same.
    assert(tiled_ps.folded_system.structure is tiled_ps.structure.folded_structure)

    assert(tiled_ps.structure.elem.tolist() == ref_tiled_elem)
    np.testing.assert_allclose(tiled_ps.structure.axes, ref_tiled_axes)
    np.testing.assert_allclose(tiled_ps.structure.pos,  ref_tiled_pos)

    assert(tiled_ps.folded_system.structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.axes, ref_folded_axes)
    np.testing.assert_allclose(tiled_ps.folded_system.structure.pos,  ref_folded_pos)

    assert(tiled_ps.ions      == ref_tiled_ions)
    assert(tiled_ps.electrons == ref_tiled_electrons)
    assert(tiled_ps.positrons == ref_tiled_positrons)

    assert(tiled_ps.folded_system.ions      == ref_folded_ions)
    assert(tiled_ps.folded_system.electrons == ref_folded_electrons)
    assert(tiled_ps.folded_system.positrons == ref_folded_positrons)

    # Make a copy of the tiled physical system, re-check everything.
    tiled_ps_copy = tiled_ps.copy()
    # Copied system is in sync with itself
    assert(tiled_ps_copy.folded_system.structure is tiled_ps_copy.structure.folded_structure)
    # Copied system's folded system's structure is not the original folded structure
    assert(tiled_ps_copy.folded_system.structure is not tiled_ps.structure.folded_structure)
    # Original system's folded system's structure is not the copy's folded structure
    assert(tiled_ps.folded_system.structure is not tiled_ps_copy.structure.folded_structure)

    assert(tiled_ps_copy.structure.elem.tolist() == ref_tiled_elem)
    np.testing.assert_allclose(tiled_ps_copy.structure.axes, ref_tiled_axes)
    np.testing.assert_allclose(tiled_ps_copy.structure.pos,  ref_tiled_pos)

    assert(tiled_ps_copy.folded_system.structure.elem.tolist() == ref_folded_elem)
    np.testing.assert_allclose(tiled_ps_copy.folded_system.structure.axes, ref_folded_axes)
    np.testing.assert_allclose(tiled_ps_copy.folded_system.structure.pos,  ref_folded_pos)

    assert(tiled_ps_copy.ions      == ref_tiled_ions)
    assert(tiled_ps_copy.electrons == ref_tiled_electrons)
    assert(tiled_ps_copy.positrons == ref_tiled_positrons)

    assert(tiled_ps_copy.folded_system.ions      == ref_folded_ions)
    assert(tiled_ps_copy.folded_system.electrons == ref_folded_electrons)
    assert(tiled_ps_copy.folded_system.positrons == ref_folded_positrons)
#end def test_copy


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Can not split into up- and down-spin with spin_orbit=True",
)
def test_spin_orbit_fail_updown():
    electrons = Electrons(count=16, n_unpaired=1, spin_orbit=True)
    electrons.n_up_down()
#end def test_spin_orbit_fail_updown


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Can not split into up-spin with spin_orbit=True",
)
def test_spin_orbit_fail_up():
    electrons = Electrons(count=16, n_unpaired=1, spin_orbit=True)
    electrons.n_up
#end def test_spin_orbit_fail_up


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Can not split into down-spin with spin_orbit=True",
)
def test_spin_orbit_fail_down():
    electrons = Electrons(count=16, n_unpaired=1, spin_orbit=True)
    electrons.n_down
#end def test_spin_orbit_fail_down


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Multiplicity undefined with spin_orbit=True",
)
def test_spin_orbit_fail_multiplicity():
    electrons = Electrons(count=16, n_unpaired=1, spin_orbit=True)
    electrons.multiplicity
#end def test_spin_orbit_fail_multiplicity


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Multiplicity undefined with fractional count",
)
def test_fractional_count_fail_multiplicity():
    electrons = Electrons(count=16.5, n_unpaired=1, spin_orbit=False)
    electrons.multiplicity
#end def test_fractional_count_multiplicity
