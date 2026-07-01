##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

"""Representations of electrons, ions, and complete physical systems."""

from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Mapping, Iterable
from copy import deepcopy
import os
from os import PathLike
from pathlib import Path
from typing import Self, TypeAlias, Literal

import numpy as np
import numpy.typing as npt

from .developer import obj, warn, error
from .periodic_table import Elements, ElementLike
from .structure import Structure, generate_structure, read_structure

LabelNumMap: TypeAlias = Mapping[str, int | float]
"""Mapping (e.g. ``dict`` or ``obj``) from an ion label to a number."""

ElemNumMap: TypeAlias = Mapping[ElementLike, int | float]
"""Mapping (e.g. ``dict`` or ``obj``) from an element to a number."""


class ElectronsPositronsBase(ABC):
    """Base class for ``Electrons`` and ``Positrons``.

    Note that this class does not make guarantees about having an
    integer amount of particles, but provides ``is_fractional`` to check
    for non-integer numbers of particles. If the class is created with
    integer-value floats for ``count`` and ``spin``, they will be
    converted into ints.

    Attributes
    ----------
    count : int or float, property
        The total number of particles
    spin : int or float, property
        The total spin of the system.

        An up-spin particle has a spin of +1/2, a down-spin particle has
        a spin of -1/2.
    n_up : int or float, read-only property
        The number of up-spin particles.
        Not defined for spin-orbit systems.
    n_down : int or float, read-only property
        The number of down-spin particles.
        Not defined for spin-orbit systems.
    total_charge : int or float, read-only property
        The total charge of the particles, equal to
        ``self.unit_charge * self.count``.
    multiplicity : int or float, read-only property
        The spin multiplicity of the particles, equal to :math:`2S+1`
        where :math:`S` is the spin.
    """

    @property
    @abstractmethod
    def unit_charge(self) -> int:
        ...

    @property
    def mass(self) -> float:
        return 1.0

    def __init__(
        self,
        count     : int | float,
        spin      : int | float,
        spin_orbit: bool = False,
        ):
        self.count      = count
        self.spin       = spin
        self.spin_orbit = spin_orbit

    @property
    def count(self) -> int | float:
        return self._count

    @count.setter
    def count(self, new_count: int | float) -> None:
        if abs((int_count := int(new_count)) - new_count) < 1e-8:
            self._count = int_count
        else:
            self._count = new_count

    @property
    def spin(self) -> int | float:
        return self._spin

    @spin.setter
    def spin(self, new_spin: int | float) -> None:
        if abs((int_spin := int(new_spin)) - new_spin) < 1e-8:
            self._spin = int_spin
        else:
            self._spin = new_spin

    def n_up_down(self) -> tuple[int, int] | tuple[float, float]:
        """Return a tuple representing the number of up- and down-spin electrons.

        Examples
        --------
        >>> Electrons(count=16, spin=1).n_up_down()
        (9, 7) # (up, down)
        >>> Electrons(count=15, spin=1).n_up_down()
        (8.5, 6.5)
        >>> Electrons(count=16, spin=1/2).n_up_down()
        (8.5, 7.5)
        >>> Electrons(count=15, spin=1/2).n_up_down()
        (8, 7)
        >>> Electrons(count=15, spin=-1/2).n_up_down()
        (7, 8)
        """
        if self.spin_orbit:
            # Use self.__class__.__name__ to get name of subclass, not base class
            raise RuntimeError(
                f"{self.__class__.__name__} can not be split into up- and down-spin with a spin-orbit system!"
                )
        if isinstance(self.count, int):
            if isinstance(self.spin, int):
                if self.count % 2 == 0:
                    n_up   = (self.count // 2) + self.spin
                    n_down = (self.count // 2) - self.spin
                else:
                    n_up   = (self.count / 2) + self.spin
                    n_down = (self.count / 2) - self.spin
            else:
                if self.count % 2 == 0:
                    n_up   = (self.count / 2) + self.spin
                    n_down = (self.count / 2) - self.spin
                else:
                    if self.spin > 0:
                        n_down = self.count // 2
                        n_up   = self.count - n_down
                    else:
                        n_up   = self.count // 2
                        n_down = self.count - n_up
        else:
            n_up   = (self.count / 2) + self.spin
            n_down = (self.count / 2) - self.spin

        return n_up, n_down

    @property
    def n_up(self) -> int | float:
        if self.spin_orbit:
            raise RuntimeError(
                f"Up-spin {self.__class__.__name__} are not defined with a spin-orbit system!"
                )
        return self.n_up_down()[0]

    @property
    def n_down(self) -> int | float:
        if self.spin_orbit:
            raise RuntimeError(
                f"Down-spin {self.__class__.__name__} are not defined with a spin-orbit system!"
                )
        return self.n_up_down()[1]

    def is_fractional(self) -> bool:
        """Returns ``True`` if the count of particles is not an int."""
        return isinstance(self.count, float)

    @property
    def total_charge(self) -> int | float:
        return self.unit_charge * self.count

    @property
    def multiplicity(self) -> int | float:
        """Defined as :math:`2S+1` where :math:`S` is ``self.spin``.

        Undefined for spin-orbit systems and fractional counts.
        """
        if self.spin_orbit:
            raise RuntimeError("Multiplicity is undefined for spin-orbit systems!")
        elif self.is_fractional():
            raise RuntimeError("Multiplicity is undefined for fractional counts!")
        else:
            return (2 * abs(self.spin)) + 1

    def __eq__(self, other: Self) -> bool:
        return (
            self.unit_charge    == other.unit_charge
            and self.count      == other.count
            and self.spin       == other.spin
            and self.spin_orbit is other.spin_orbit
            )

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}("
            f"count={self.count}, "
            f"spin={self.spin}, "
            f"spin_orbit={self.spin_orbit})"
            )
#end class ElectronsPositronsBase


class Electrons(ElectronsPositronsBase):
    """Class representing a collection of electrons."""
    @property
    def unit_charge(self) -> int:
        return -1

    @classmethod
    def neutralize_to(
        cls,
        ions         : Iterable[IonSpecies],
        total_charge : int | float,
        electron_spin: int | float | None = None,
        spin_orbit   : bool = False,
    ) -> Self:
        """Neutralize the charge of ``ions`` to ``total_charge``.

        This will prioritize creating an integer number of electrons
        with the specified spin, however it will fall back to a
        fractional number of electrons if necessary.

        Parameters
        ----------
        ions : Iterable of IonSpecies
            A ``dict`` or ``obj`` with ``IonSpecies`` as values or a
            list of ``IonSpecies``.
        total_charge : int or float
            The desired total charge of the combined ion-electron system.
        electron_spin : int or float, optional
            The desired total spin of the electrons. If this is not
            specified, then this will be set to the smallest, positive,
            half-integer or integer spin state that is compatible with
            the number of electrons.
        spin_orbit : bool, default=False
            Specify whether or not the system is a spin-orbit system.
            Passed to the class constructor.
        """
        if isinstance(ions, Mapping):
            ions = ions.values()

        ions_charge = 0
        for ion in ions:
            ions_charge += ion.total_charge_deficit

        n_electrons = ions_charge - total_charge
        if electron_spin is None:
            if isinstance(ions_charge, int) and isinstance(n_electrons, int):
                if n_electrons % 2 == 0:
                    electron_spin = 0
                else:
                    electron_spin = 0.5
            else:
                electron_spin = (n_electrons % 2) / 2

        return Electrons(
            count      = n_electrons,
            spin       = electron_spin,
            spin_orbit = spin_orbit,
            )
#end class Electrons


class Positrons(ElectronsPositronsBase):
    """Class representing a collection of positrons."""
    @property
    def unit_charge(self) -> int:
        return 1
#end class Positrons


class IonSpecies:
    """Class representing a collection of ions of the same type.

    Attributes
    ----------
    element : Elements
        The element for this ion collection.
    count : int or float
        The number of ions in this collection.
    label : str
        The label for the ion collection.
    formal_charge : int or float
        The formal charge associated with a single one of the ions.
    unit_spin : int or float
        The spin of a single one of the ions.
    Zeff : int or float
        The effective nuclear charge of one of the ions.
    symbol : str, read-only property
        The atomic symbol of the element.
    total_formal_charge : int or float, read-only property
        Formal charge multiplied by count.
    total_nuclear_charge : int or float, read-only property
        Zeff multiplied by count.
    electron_deficit : int or float, read-only property
        Number of electrons required to reach the formal charge.
    total_electron_deficit : int or float, read-only property
        Number of electrons required to reach the formal charge
        multiplied by count.
    total_spin : int or float, read-only property
        Unit spin multiplied by count.
    total_mass : float, read-only property
        Atomic weight multiplied by count.

    Parameters
    ----------
    element : ElementLike
        A member of the ``Elements`` enum, atomic symbol, or atomic
        number.
    count : int or float
        The number of ions in this collection.
    label : str, optional
        The label for the ion. If not given, defaults to
        ``element.symbol``.
    formal_charge : int or float, default=0
        The formal charge associated with a single one of the ions.
    unit_spin : int or float, default=0
        The spin of a single one of the ions.
    Zeff : int, optional
        The effective nuclear charge of the ion. Defaults to the atomic
        number (a.k.a. all-electron).
    """

    def __init__(
        self,
        element      : ElementLike,
        count        : int | float,
        label        : str | None         = None,
        formal_charge: int | float        = 0,
        unit_spin    : int | float        = 0,
        Zeff         : int | float | None = None,
        ):
        self.element       = Elements(element)
        self.count         = count
        self.label         = label if label is not None else self.element.symbol
        self.formal_charge = formal_charge
        self.unit_spin     = unit_spin
        self.Zeff          = Zeff if Zeff is not None else self.element.atomic_number

    def pseudized(self) -> bool:
        """Check if this ion is pseudized."""
        return self.Zeff != self.element.atomic_number

    def pseudize(self, Zeff: int):
        """Equivalent to setting ``self.Zeff = Zeff``."""
        self.Zeff = Zeff

    def is_ghost(self) -> bool:
        """Check if this collection of ions represents ghost atoms."""
        return self.element is Elements.Unknown

    @property
    def symbol(self) -> str:
        """Atomic symbol of the element."""
        return self.element.symbol

    @property
    def total_mass(self) -> float:
        """Total mass of all ions in the collection."""
        return self.element.atomic_weight * self.count

    @property
    def charge_deficit(self) -> int | float:
        """The required number of electrons to get to the formal charge.

        Equal to the effective nuclear charge minus the formal charge.

        For example, a single all-electron oxygen (``Zeff=8``) that has 
        ``formal_charge=-1`` needs :math:`8 - (-1) = 8 + 1 = 9`
        electrons to achieve the desired format charge.
        """
        return self.Zeff - self.formal_charge

    @property
    def total_nuclear_charge(self) -> int | float:
        """Nuclear charge times count."""
        return self.Zeff * self.count

    @property
    def total_formal_charge(self) -> int | float:
        """Unit charge times count."""
        return self.formal_charge * self.count

    @property
    def total_charge_deficit(self) -> int | float:
        """Number of electrons to achieve the desired formal charge times count."""
        return self.charge_deficit * self.count

    @property
    def total_spin(self) -> int | float:
        """Total spin of all ions in the collection."""
        return self.unit_spin * self.count

    def __eq__(self, other: Self) -> bool:
        return (
            self.element is other.element
            and self.count         == other.count
            and self.label         == other.label
            and self.formal_charge == other.formal_charge
            and self.unit_spin     == other.unit_spin
            and self.Zeff          == other.Zeff
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"element={self.element}, "
            f"count={self.count}, "
            f"label={self.label}, "
            f"formal_charge={self.formal_charge}, "
            f"unit_spin={self.unit_spin}, "
            f"Zeff={self.Zeff})"
            )

    def __hash__(self) -> int: # Enables making unordered sets
        return hash((
            self.element,
            self.count,
            self.label,
            self.formal_charge,
            self.unit_spin,
            self.Zeff,
            ))

    @classmethod
    def from_structure(
        cls,
        structure  : Structure,
        elem_charge: LabelNumMap = dict(),
        elem_spin  : LabelNumMap = dict(),
        elem_Zeff  : LabelNumMap = dict(),
        ) -> dict[str, Self]:
        """Create a dict with ``IonSpecies`` from a ``Structure`` object.

        It is important to note that this class only represents the ions
        in a structure, not any other particles (e.g. electrons). Thus,
        this will not capture the background charge of the structure if
        it is defined.

        Parameters
        ----------
        structure : Structure
            The structure from which to pull ions.
        elem_charge : Mapping[str, int or float], optional
            A dict or ``obj`` mapping the elements to formal charges.
            Defaults to 0 if not given.
        elem_spin : Mapping[str, int or float], optional
            A dict or ``obj`` mapping the elements to spin states.
            Defaults to 0 if not given.
        elem_Zeff : Mapping[str, int or float], optional
            A dict or ``obj`` mapping the elements to effective nuclear
            charges.
            Defaults to the atomic number if not supplied.

        Returns
        -------
        ions : dict[str, IonSpecies]
            A dictionary of the ion species found in the structure.
            The ion labels are used for the keys, sorted in alphabetical
            order.

        Examples
        --------
        Minimal call signature, populated with defaults.

        >>> structure = Structure(
        ...     elem = ["N", "C1", "C2", "O", "O", "H", "H", "H", "H", "H"],
        ...     pos = np.zeros((10,3)),
        ...     )
        >>> ions = IonSpecies.from_structure(
        ...     structure=structure,
        ...     )
        >>> for label, ion in ions.items():
        ...     print(f"{label:2}: {repr(ion)}")
        C1: IonSpecies(element=C, count=1, label=C1, formal_charge=0, unit_spin=0, Zeff=6)
        C2: IonSpecies(element=C, count=1, label=C2, formal_charge=0, unit_spin=0, Zeff=6)
         H: IonSpecies(element=H, count=5, label=H, formal_charge=0, unit_spin=0, Zeff=1)
         N: IonSpecies(element=N, count=1, label=N, formal_charge=0, unit_spin=0, Zeff=7)
         O: IonSpecies(element=O, count=2, label=O, formal_charge=0, unit_spin=0, Zeff=8)

        Full call signature, all values specified.

        >>> structure = Structure(
        ...     elem = ["N", "C1", "C2", "O", "O", "H", "H", "H", "H", "H"],
        ...     pos = np.zeros((10,3)),
        ...     )
        >>> ions = IonSpecies.from_structure(
        ...     structure   = structure,
        ...     elem_charge = dict(N=3, C1=2, C2=4,   O=2, H=1  ),
        ...     elem_spin   = dict(N=1, C1=0, C2=0.5, O=0, H=0.5),
        ...     elem_Zeff   = dict(N=5, C1=4, C2=6,   O=6, H=1  ),
        ...     )
        >>> for label, ion in ions.items():
        ...     print(f"{label:2}: {repr(ion)}")
        C1: IonSpecies(element=C, count=1, label=C1, formal_charge=2, unit_spin=0, Zeff=4)
        C2: IonSpecies(element=C, count=1, label=C2, formal_charge=4, unit_spin=0.5, Zeff=6)
         H: IonSpecies(element=H, count=5, label=H, formal_charge=1, unit_spin=0.5, Zeff=1)
         N: IonSpecies(element=N, count=1, label=N, formal_charge=3, unit_spin=1, Zeff=5)
         O: IonSpecies(element=O, count=2, label=O, formal_charge=2, unit_spin=0, Zeff=6)
        """
        ions = {}
        elem_list = list(structure.elem)
        elem_set = set(elem_list)
        for label in elem_set:
            is_elem, element = Elements.is_element(label, return_element=True)
            if not is_elem:
                raise ValueError(
                    f"Can not determine element from label {label}!"
                    )

            ion = cls(
                element       = element,
                count         = elem_list.count(label),
                label         = label,
                formal_charge = elem_charge.get(label, 0),
                unit_spin     = elem_spin.get(label, 0),
                Zeff          = elem_Zeff.get(label, element.atomic_number),
                )
            ions[label] = ion
        # Sort so keys (ion labels) are in alphabetical order.
        # This is ever so slightly better than the random order
        # that comes from using set(elem_list) above.
        ions = {lbl: ions[lbl] for lbl in sorted(ions.keys())}
        return ions
#end class IonSpecies


class PhysicalSystem:
    """A system of ions and electrons with a structure.

    The ``PhysicalSystem`` is used to create inputs for all simulations.
    The difference between a ``PhysicalSystem`` and a ``Structure`` is
    that a ``Structure`` contains no information about charge, spin,
    pseudopotentials, elements, or electrons. The ``PhysicalSystem``
    class contains a ``Structure``, but it also includes the additional
    data required to define a real system that would be used in a
    calculation.

    Attributes
    ----------
    structure : Structure
        The structure of the ions in the system.
    ions : Mapping[str, IonSpecies]
        The unique ions of the system.
    electrons : Electrons
        The electrons in the system.
    positrons : Positrons or None
        The positrons in the system.
    folded : PhysicalSystem or None
        If ``structure`` has a folded structure, this will be the
        ``PhysicalSystem`` that corresponds to the folded structure.

    Parameters
    ----------
    structure : Structure or PathLike
        The structure to use for creating the ``PhysicalSystem``, or a
        path to a file with structural information.
    total_charge : int or float, optional
        The total charge of the system. Will use
        ``Structure.background_charge`` if it is not ``None``,
        otherwise it will default to zero.
    electron_spin : int or float, optional
        Set the total spin of the electrons.

        Note that up-spin electrons are +1/2 and down-spin are -1/2.
    spin_orbit : bool, default=False
        Signal that the system has spin-orbit coupling.
    elem_charge : Mapping[str, int or float], optional
        Mapping from element labels to their formal charges.
    elem_spin : Mapping[str, int or float], optional
        Mapping from element labels to their spin states.
    elem_Zeff : Mapping[str, int or float], optional
        Mapping from element labels to effective nuclear charges.
    positrons : Positrons, optional
        Additional positrons to add to the system.

    See Also
    --------
    IonSpecies.from_structure()
    Electrons.neutralize_to()
    """

    def __init__(
        self,
        structure    : Structure | PathLike,
        total_charge : int | float | None = None,
        electron_spin: int | float | None = None,
        spin_orbit   : bool               = False,
        elem_charge  : LabelNumMap        = dict(),
        elem_spin    : LabelNumMap        = dict(),
        elem_Zeff    : LabelNumMap        = dict(),
        positrons    : Positrons | None   = None,
        ):
        if isinstance(structure, str | bytes | Path):
            s = Structure()
            s.read(structure)
            self.structure = s
        elif isinstance(structure, Structure):
            self.structure = structure
        else:
            error("The `structure` parameter for PhysicalSystem must be a file path or a Structure!")

        self.ions = IonSpecies.from_structure(
            structure   = self.structure,
            elem_charge = elem_charge,
            elem_spin   = elem_spin,
            elem_Zeff   = elem_Zeff,
            )

        if total_charge is None:
            if self.structure.background_charge is not None:
                total_charge = self.structure.background_charge
            else:
                total_charge = 0
                for ion in self.ions.values():
                    total_charge += ion.total_formal_charge

        self.electrons = Electrons.neutralize_to(
            ions          = self.ions,
            total_charge  = total_charge,
            electron_spin = electron_spin,
            spin_orbit    = spin_orbit,
            )
        self.positrons = positrons
        self.folded_system = self._process_folded_structure()

    def _process_folded_structure(self) -> Self | None:
        """Get the appropriate folded system from ``self.structure.folded_structure``.

        If ``self.structure.folded_structure`` is ``None``, this will
        just return ``None``.
        """
        if not self.structure.has_folded():
            return None

        if (msg := self.structure.check_tiling(exit=False)) != "":
            raise RuntimeError(
                "Folded structure is not consistent, can not initialize PhysicalSystem!\n"
                f"Reason:\n{msg}"
                )

        if not hasattr(self.structure, "tmatrix") or self.structure.tmatrix is None:
            # Folded molecule
            n_cells_tiled = 1
        else:
            # Folded crystal
            n_cells_tiled = int(np.round(np.abs(np.linalg.det(self.structure.tmatrix))))

        # Ignore positrons for this, since we use this to get the number of electrons
        folded_total_charge = (self.ion_charge + self.electron_charge) / n_cells_tiled
        if abs((int_fold_chg := int(folded_total_charge)) - folded_total_charge) < 1e-8:
            folded_total_charge = int_fold_chg

        folded_elec_spin = self.electron_spin / n_cells_tiled
        if abs((int_fold_spin := int(folded_elec_spin)) - folded_elec_spin) < 1e-8:
            folded_elec_spin = int_fold_spin

        elem_chg_map  = {}
        elem_spin_map = {}
        elem_Zeff_map = {}
        for label, species in self.ions.items():
            elem_chg_map[label]  = species.formal_charge
            elem_spin_map[label] = species.unit_spin
            elem_Zeff_map[label] = species.Zeff

        if self.positrons is not None:
            if self.positrons.count % n_cells_tiled == 0:
                folded_positron_count = self.positrons.count // n_cells_tiled
            else:
                folded_positron_count = self.positrons.count / n_cells_tiled

            if self.positrons.spin % n_cells_tiled == 0:
                folded_positron_spin = self.positrons.spin // n_cells_tiled
            else:
                folded_positron_spin = self.positrons.spin / n_cells_tiled

            folded_positrons = Positrons(
                count      = folded_positron_count,
                spin       = folded_positron_spin,
                spin_orbit = self.positrons.spin_orbit,
                )
        else:
            folded_positrons = None

        return PhysicalSystem(
            structure     = self.structure.folded_structure,
            total_charge  = folded_total_charge,
            electron_spin = folded_elec_spin,
            spin_orbit    = self.electrons.spin_orbit,
            elem_charge   = elem_chg_map,
            elem_spin     = elem_spin_map,
            elem_Zeff     = elem_Zeff_map,
            positrons     = folded_positrons,
            )

    @property
    def ion_charge(self) -> int | float:
        """The total charge of all ions in the system."""
        return sum([ion.total_nuclear_charge for ion in self.ions.values()])

    @property
    def ion_spin(self) -> int | float:
        """The total spin of all ions in the system."""
        return sum([ion.total_spin for ion in self.ions.values()])

    @property
    def electron_charge(self) -> int | float:
        """Alias for ``PhysicalSystem.electrons.total_charge``"""
        return self.electrons.total_charge

    @property
    def electron_spin(self) -> int | float:
        """Alias for ``PhysicalSystem.electrons.spin``"""
        return self.electrons.spin

    @property
    def positron_charge(self) -> int | float | None:
        """Returns ``None`` if ``PhysicalSystem.positrons`` is ``None``."""
        if self.positrons is not None:
            return self.positrons.total_charge
        else:
            return None

    @property
    def positron_spin(self) -> int | float | None:
        """Returns ``None`` if ``PhysicalSystem.positrons`` is ``None``."""
        if self.positrons is not None:
            return self.positrons.spin
        else:
            return None

    @property
    def net_charge(self) -> int | float:
        """Ion charge plus electron charge plus positron charge."""
        if self.positrons is not None:
            return self.ion_charge + self.electron_charge + self.positron_charge
        else:
            return self.ion_charge + self.electron_charge

    @property
    def pseudized(self) -> bool:
        """Return ``True`` if **any** of the ions are pseudized."""
        pseudized = False
        for ion in self.ions.values():
            pseudized |= ion.pseudized()

        return pseudized

    def num_ions(self) -> int | float:
        return sum([ion.count for ion in self.ions.values()])

    def num_ion_types(self) -> int:
        return len(self.ions)

    def spin_polarized(self) -> bool:
        spn = self.electron_spin + sum([ion.total_spin for ion in self.ions.values()])
        return spn != 0 or self.structure.is_magnetic()

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}("
            f"structure=..., " # Not sure how to do this yet.
            f"ions={self.ions!s}, "
            f"electrons={self.electrons!s}, "
            f"positrons={self.positrons!s})"
            )

    def check_folded_system(self, exit=True, message=False) -> bool | tuple[bool, str]:
        """Check if the current folded system matches the folded structure.

        Parameters
        ----------
        exit : bool, default=True
            Whether or not to raise an error on mismatch.
        message : bool, default=False
            Whether or not to include the message with the output.
            Even if the systems match this will return the message along
            with the success state, in which case the message will be an
            empty string.

        Returns
        -------
        success : bool
            ``True`` if the systems match, ``False`` otherwise.
        msg : str, optional
            Only returned if ``message=True``.

        Raises
        ------
        ValueError
            If the systems do not match and ``exit=True``.
        """
        msg = ""
        # Either they are both None or they both point to the same object.
        if self.folded_system is None:
            if message:
                return True, msg
            else:
                return True

        if self.folded_system.structure is not self.structure.folded_structure:
            msg += "The folded structure is not the folded system's structure!\n"

        tmp_folded_system = self._process_folded_structure()
        if self.folded_system.ions != tmp_folded_system.ions:
            msg += "The folded system's ions do not match the expected ions!\n"
            msg += "Current Folded Ions:\n\n"
            for label, ion in self.folded_system.ions.items():
                msg += f"    {label}: {ion!r}\n"

            msg += "\nExpected Folded Ions:\n\n"
            for label, ion in tmp_folded_system.ions.items():
                msg += f"    {label}: {ion!r}\n"

        if self.folded_system.electrons != tmp_folded_system.electrons:
            msg += "The folded system's electrons do not match the expected electrons!\n"
            msg += f"Current Folded Electrons:  {self.folded_system.electrons!r}\n"
            msg += f"Expected Folded Electrons: {tmp_folded_system.electrons!r}\n"
        
        if self.folded_system.positrons != tmp_folded_system.positrons:
            msg += "The folded system's positrons do not match the expected positrons!\n"
            msg += f"Current Folded Positrons:  {self.folded_system.positrons!r}\n"
            msg += f"Expected Folded Positrons: {tmp_folded_system.positrons!r}\n"

        success = len(msg) == 0
        if exit and not success:
            raise ValueError(msg)
        elif message:
            return success, message
        else:
            return success

    def check_consistent(
        self,
        tol    : float = 1e-8,
        exit   : bool  = True,
        message: bool  = False,
        ) -> bool | tuple[bool, str]:
        """Check that the folded system and structure match their tiled versions.

        See ``PhysicalSystem.check_folded_system`` and
        ``Structure.check_consistent`` for a description of the
        parameters.
        """
        fs, fm = self.check_folded_system(exit=False, message=True)
        cs, cm = self.structure.check_consistent(tol, exit=False, message=True)
        msg = ""
        if not fs:
            msg += fm + "\n"

        if not cs:
            msg += cm + "\n"

        consistent = len(msg) == 0
        if not consistent and exit:
            raise ValueError(msg)

        if not message:
            return consistent
        else:
            return consistent, msg

    def is_valid(self) -> bool:
        """Check if the system's folded system and structure are valid."""
        return self.check_consistent(exit=False)

    def change_units(self, units: str) -> None:
        """Alias for ``self.structure.change_units()``."""
        self.structure.change_units(units, folded=False)
        if self.folded_system is not None:
            self.folded_system.change_units(units)

    def group_atoms(self) -> None:
        """Alias for ``self.structure.group_atoms()``."""
        self.structure.group_atoms(folded=False)
        if self.folded_system is not None:
            self.folded_system.group_atoms()

    def rename(self, **name_pairs):
        """Replace the labels for the ions.

        This will also rename ``self.structure.elem``.
        """
        if name_pairs.pop("folded", None) is not None:
            warn(
                "Passing `folded` to PhysicalSystem.rename() is no longer supported.\n"
                "This function will now always rename the folded system to maintain parity."
            )
        old_keys = list(self.ions.keys())
        for old, new in name_pairs.items():
            if old not in self.ions:
                warn(
                    f"The label '{old}' does not exist in the current system.\n"
                    f"Known labels are {old_keys}"
                    )
            else:
                self.ions[new] = self.ions[old]
                self.ions[new].label = new
                del self.ions[old]

        # folded=False because we will do it through the folded system.
        self.structure.rename(folded=False, **name_pairs)
        if self.folded_system is not None:
            self.folded_system.rename(**name_pairs)

    def tile(self, *td) -> Self:
        """Tile an existing system.

        Note that this function returns a brand new ``PhysicalSystem``
        instance, with a folded system that matches the old system. The
        old system is not changed.
        """
        supercell = self.structure.tile(*td)

        n_cells_tiled = int(np.round(np.abs(np.linalg.det(supercell.tmatrix))))
        # Copy over data from the old system
        elem_chg_map  = {}
        elem_spin_map = {}
        elem_Zeff_map = {}
        for label, species in self.ions.items():
            elem_chg_map[label]  = species.formal_charge
            elem_spin_map[label] = species.unit_spin
            elem_Zeff_map[label] = species.Zeff

        if self.positrons is not None:
            tiled_positrons = Positrons(
                count      = self.positrons.count * n_cells_tiled,
                spin       = self.positrons.spin * n_cells_tiled,
                spin_orbit = self.positrons.spin_orbit,
                )
        else:
            tiled_positrons = None

        tiled_system = PhysicalSystem(
            structure     = supercell,
            total_charge  = (self.ion_charge + self.electron_charge) * n_cells_tiled,
            electron_spin = self.electron_spin * n_cells_tiled,
            elem_charge   = elem_chg_map,
            elem_spin     = elem_spin_map,
            elem_Zeff     = elem_Zeff_map,
            positrons     = tiled_positrons,
            )
        # TODO: Check that the old system is equal to the new folded system.
        return tiled_system

    def has_folded(self) -> bool:
        """Check if the system has a folded system."""
        return self.folded_system is not None

    def remove_folded_system(self) -> None:
        """Remove the existing folded system, if it exists."""
        self.folded_system = None
        self.structure.remove_folded_structure()

    def get_smallest(self) -> Self:
        """Return the folded system if it exists, otherwise return itself."""
        if self.has_folded():
            return self.folded_system
        else:
            return self

    def is_magnetic(self) -> bool:
        return self.electron_spin != 0 or self.ion_spin != 0

    def large_Zeff_elem(self, Zmin: int | float) -> list[str]:
        """Get the element labels with a ``Zeff`` greater than ``Zmin``."""
        elem = []
        for label, ion in self.ions.items():
            if ion.Zeff > Zmin:
                elem.append(label)

        return elem

    def ae_pp_species(self) -> tuple[set[str], set[str] | set]:
        """Get the all-electron and pseudized species, respectively."""
        if self.pseudized:
            pp_species = set([lbl for lbl, ion in self.ions.items() if ion.pseudized()])
            ae_species = set(self.ions.keys()) - pp_species
        else: # All species are all-electron
            pp_species = set()
            ae_species = set(self.ions.keys())

        return ae_species, pp_species

    def kf_rpa(self) -> npt.NDArray[np.floating]:
        n_elecs = self.electrons.n_up_down()
        n_elecs = [float(e) for e in n_elecs]

        # k-space volume per particle
        kvol1 = (2 * np.pi)**3 / self.structure.volume()
        kfs = []
        for n_elec in n_elecs:
            kf = (3 * n_elec * kvol1 / (4 * np.pi))**(1./3)
            kfs.append(kf)

        return np.array(kfs, dtype=float)

    def copy(self) -> Self:
        return deepcopy(self)
#end class PhysicalSystem


ps_defaults = dict(
    type='crystal',
    kshift=(0,0,0),
    net_charge=0,
    net_spin=0,
    pretile=None,
    tiling=None,
    tiled_spin=None,
    extensive=True,
    )
def generate_physical_system(**kwargs) -> PhysicalSystem:
    for var,val in ps_defaults.items():
        if var not in kwargs:
            kwargs[var] = val

    type = kwargs['type']
    if type=='atom' or type=='dimer' or type=='trimer':
        del kwargs['kshift']
        del kwargs['tiling']
        tiling = None
    else:
        tiling = kwargs['tiling']

    if 'structure' in kwargs:
        s = kwargs['structure']
        if isinstance(s, str | Path):
            if os.path.exists(s):
                if 'elem' in kwargs:
                    s = read_structure(s,elem=kwargs['elem'])
                else:
                    s = read_structure(s)

                if 'axes' in kwargs:
                    s.reset_axes(kwargs['axes'])

                kwargs['structure'] = s
            else:
                s_low = s.lower()
                format = None
                if '.' in s_low:
                    format = s_low.rsplit('.')[1]
                elif 'poscar' in s_low:
                    format = 'poscar'

                if '/' in s or format in set(["xyz", "xsf", "poscar", "cif", "fhi-aims"]):
                    error(
                        'User provided structure file does not exist\n'
                        'Structure file path: '+s
                        )

    generation_info = obj()
    generation_info.transfer_from(deepcopy(kwargs))

    net_charge = kwargs.pop("net_charge")
    net_spin   = kwargs.pop("net_spin")
    tiled_spin = kwargs.pop("tiled_spin")
    extensive  = kwargs.pop("extensive")

    if 'particles' in kwargs:
        warn("generate_physical_system no longer supports the `particles` parameter.")
        del kwargs['particles']

    elem_Zeff_map = dict()
    elems_to_delete = []
    for var in kwargs:
        if Elements.is_element(var):
            elem_Zeff_map[var] = kwargs[var]
            elems_to_delete.append(var)

    for var in elems_to_delete:
        del kwargs[var]

    pretile = kwargs.pop("pretile", None)
    if pretile is None:
        structure = generate_structure(**kwargs)
    else:
        for d in range(len(pretile)):
            if tiling[d] % pretile[d] != 0:
                error(
                    'pretile does not divide evenly into tiling\n'
                    '  tiling provided: {0}\n'
                    '  pretile provided: {1}'
                    .format(tiling, pretile)
                    )

        tiling = tuple(np.array(tiling) // np.array(pretile))
        kwargs['tiling'] = pretile
        pre = generate_structure(**kwargs)
        pre.remove_folded_structure()
        structure = pre.tile(tiling)

    if isinstance(tiling, tuple):
        tiling_mat = np.diag(tiling)
    elif tiling is None:
        tiling_mat = np.eye(3)
    else:
        tiling_mat = tiling

    if not np.array_equal(tiling_mat, np.eye(3)) and structure.has_folded():
        # Has some supercell tiling
        fps = PhysicalSystem(
            structure     = structure.folded_structure,
            total_charge  = net_charge,
            electron_spin = net_spin,
            elem_Zeff     = elem_Zeff_map,
            )
        structure.remove_folded()
        folded_structure = fps.structure
        if extensive:
            ncells = int(round(structure.volume()/folded_structure.volume()))
            net_charge = ncells * net_charge
            if not isinstance(net_spin, str):
                net_spin = ncells * net_spin

        if tiled_spin is not None:
            net_spin = tiled_spin

        ps = PhysicalSystem(
            structure     = structure,
            total_charge  = net_charge,
            electron_spin = net_spin,
            elem_Zeff     = elem_Zeff_map,
            )
        structure.set_folded(folded_structure)
        ps.folded_system = fps
    else:
        # No supercell tiling
        ps = PhysicalSystem(
            structure     = structure,
            total_charge  = net_charge,
            electron_spin = net_spin,
            elem_Zeff     = elem_Zeff_map,
            )

    # For now we will store it.
    ps.generation_info = generation_info

    return ps
#end def generate_physical_system



# test needed
def ghost_atoms(*particles):
    for particle in particles:
        ...
    #end for
#end def ghost_atoms
