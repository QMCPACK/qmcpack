from __future__ import annotations
import os
from os import PathLike
from collections.abc import Sequence
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Self, TypeVar
from enum import Enum
import builtins

import numpy as np
import numpy.typing as npt
from numpy.linalg import inv

from .physical_system import PhysicalSystem
from .structure import kmesh


PwscfInputType = str | bool | int | float | Sequence
# For structured array inputs (k-points, weights, reciprocal lattice, etc.)
PwscfArrayInput = npt.NDArray[np.floating] | Sequence[Sequence[float]] | Sequence[float]


@dataclass(frozen=True) # We generally don't want any of these to be modified
class NamelistParamDefinition:
    """Base class for all namelist variables."""

    datatype:        tuple[PwscfInputType]
    required:        bool
    shape:           tuple
    allowed_values:  tuple[PwscfInputType] | None = None
    dependencies:    tuple[NamelistParamDefinition] | None = None
    version_added:   tuple[float] | None = None
    version_removed: tuple[float] | None = None
#end class NamelistDefinition

class NamelistEnumBase(Enum, ABC):
    """Abstract base class for the namelist enumerations.
    Provides the ``__new__`` method for the enums.
    """
    def __new__(
        cls,
        datatype:        PwscfInputType,
        required:        bool,
        allowed_values:  tuple[PwscfInputType] | None = None,
        dependencies:    tuple[NamelistParamDefinition] | None = None,
        version_added:   tuple[float] | None = None,
        version_removed: tuple[float] | None = None,
    ):
        definition = NamelistParamDefinition.__new__(cls)
        definition._value_ = NamelistParamDefinition(
            datatype,
            required,
            allowed_values,
            dependencies,
            version_added,
            version_removed,
        )
        return definition
#end class NamelistEnumBase

class ControlDefinitions(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &CONTROL input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    calculation   = str,   False, (1,), ('scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md')
    title         = str,   False, (1,), None
    verbosity     = str,   False, (1,), ('high', 'low')
    restart_mode  = str,   False, (1,), ('from_scratch', 'restart')
    wf_collect    = bool,  False, (1,), None
    nstep         = int,   False, (1,), None
    iprint        = int,   False, (1,), None
    tstress       = bool,  False, (1,), None
    tprnfor       = bool,  False, (1,), None
    dt            = float, False, (1,), None
    outdir        = str,   False, (1,), None
    wfcdir        = str,   False, (1,), None
    prefix        = str,   False, (1,), None
    lkpoint_dir   = bool,  False, (1,), None
    max_seconds   = float, False, (1,), None
    etot_conv_thr = float, False, (1,), None
    forc_conv_thr = float, False, (1,), None
    disk_io       = str,   False, (1,), ('high', 'medium', 'low', 'nowf', 'minimal', 'none')
    pseudo_dir    = str,   False, (1,), None
    tefield       = bool,  False, (1,), None
    dipfield      = bool,  False, (1,), None
    lelfield      = bool,  False, (1,), None
    nberrycyc     = int,   False, (1,), None
    lorbm         = bool,  False, (1,), None
    lberry        = bool,  False, (1,), None
    gdir          = int,   False, (1,), None
    nppstr        = int,   False, (1,), None
    gate          = bool,  False, (1,), None
    twochem       = bool,  False, (1,), None
    lfcp          = bool,  False, (1,), None
    trism         = bool,  False, (1,), None
#end class ControlDefinitions

class SystemVariables(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &SYSTEM input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    ibrav                     = int,   True,  (1,), None
    nat                       = int,   True,  (1,), None
    ntyp                      = int,   True,  (1,), None
    nbnd                      = int,   False, (1,), None
    nbnd_cond                 = int,   False, (1,), None
    tot_charge                = float, False, (1,), None
    tot_magnetization         = float, False, (1,), None
    ecutwfc                   = float, True,  (1,), None
    ecutrho                   = float, False, (1,), None
    ecutfock                  = float, False, (1,), None
    nosym                     = bool,  False, (1,), None
    nosym_evc                 = bool,  False, (1,), None
    noinv                     = bool,  False, (1,), None
    no_t_rev                  = bool,  False, (1,), None
    force_symmorphic          = bool,  False, (1,), None
    use_all_frac              = bool,  False, (1,), None
    occupations               = str,   False, (1,), ('smearing', 'tetrahedra', 'tetrahedra_lin', 'tetrahedra_opt', 'fixed', 'from_input')
    one_atom_occupations      = bool,  False, (1,), None
    starting_spin_angle       = bool,  False, (1,), None
    degauss_cond              = float, False, (1,), None
    nelec_cond                = float, False, (1,), None
    degauss                   = float, False, (1,), None
    smearing                  = str,   False, (1,), ('gaussian', 'gauss', 'methfessel-paxton', 'm-p', 'mp', 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'fermi-dirac', 'f-d', 'fd')
    nspin                     = int,   False, (1,), None
    sic_gamma                 = float, False, (1,), None
    pol_type                  = str,   False, (1,), ('e', 'h')
    sic_energy                = bool,  False, (1,), None
    sci_vb                    = float, False, (1,), None
    sci_cb                    = float, False, (1,), None
    noncolin                  = bool,  False, (1,), None
    ecfixed                   = float, False, (1,), None
    qcutz                     = float, False, (1,), None
    q2sigma                   = float, False, (1,), None
    input_dft                 = str,   False, (1,), None
    ace                       = bool,  False, (1,), None
    exx_fraction              = float, False, (1,), None
    screening_parameter       = float, False, (1,), None
    exxdiv_treatment          = str,   False, (1,), ('gygi-baldereschi', 'vcut_spherical', 'vcut_ws', 'none')
    x_gamma_extrapolation     = bool,  False, (1,), None
    ecutvcut                  = float, False, (1,), None
    localization_thr          = float, False, (1,), None
    dmft                      = bool,  False, (1,), None
    dmft_prefix               = str,   False, (1,), None
    ensemble_energies         = bool,  False, (1,), None
    edir                      = int,   False, (1,), None
    emaxpos                   = float, False, (1,), None
    eopreg                    = float, False, (1,), None
    eamp                      = float, False, (1,), None
    lforcet                   = bool,  False, (1,), None
    constrained_magnetization = str,   False, (1,), ('none', 'total', 'atomic', 'total direction', 'atomic direction')
    # Lambda is reserved by Python, so it gets an underscore.
    # Will need to handle the underscore at some point.
    lambda_                   = float, False, (1,), None
    report                    = int,   False, (1,), None
    lspinorb                  = bool,  False, (1,), None
    assume_isolated           = str,   False, (1,), ('none', 'makov-payne', 'm-p', 'mp', 'martyna-tuckerman', 'm-t', 'mt', 'esm', '2D')
    esm_bc                    = str,   False, (1,), ('pbc', 'bc1', 'bc2', 'bc3')
    esm_w                     = float, False, (1,), None
    esm_efield                = float, False, (1,), None
    esm_nfit                  = int,   False, (1,), None
    lgcscf                    = bool,  False, (1,), None
    gcscf_mu                  = float, True,  (1,), None
    gcscf_conv_thr            = float, False, (1,), None
    gcscf_beta                = float, False, (1,), None
    vdw_corr                  = str,   False, (1,), ('grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d', 'grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3', 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler', 'MBD', 'mbd', 'many-body-dispersion', 'mbd_vdw', 'XDM', 'xdm')
    london                    = bool,  False, (1,), None
    london_s6                 = float, False, (1,), None
    london_rcut               = float, False, (1,), None
    dftd3_version             = int,   False, (1,), None
    dftd3_threebody           = bool,  False, (1,), None
    ts_vdw_econv_thr          = float, False, (1,), None
    ts_vdw_isolated           = bool,  False, (1,), None
    xdm                       = bool,  False, (1,), None
    xdm_a1                    = float, False, (1,), None
    xdm_a2                    = float, False, (1,), None
    space_group               = int,   False, (1,), None
    uniqueb                   = bool,  False, (1,), None
    origin_choice             = int,   False, (1,), None
    rhombohedral              = bool,  False, (1,), None
    nextffield                = int,   False, (1,), None
    celldm                    = float, False, ((1,), (6,)), None
    A                         = float, False, (1,), None
    B                         = float, False, (1,), None
    C                         = float, False, (1,), None
    cosAB                     = float, False, (1,), None
    cosAC                     = float, False, (1,), None
    cosBC                     = float, False, (1,), None
    zgate                     = float, False, (1,), None
    relaxz                    = bool,  False, (1,), None
    block                     = bool,  False, (1,), None
    block_1                   = float, False, (1,), None
    block_2                   = float, False, (1,), None
    block_height              = float, False, (1,), None
    starting_charge           = float, False, ((1,), ('ntyp',)), None
    starting_magnetization    = float, False, ((1,), ('ntyp',)), None
    Hubbard_beta              = float, False, ((1,), ('ntyp',)), None
    angle1                    = float, False, ((1,), ('ntyp',)), None
    angle2                    = float, False, ((1,), ('ntyp',)), None
    fixed_magnetization       = float, False, ((1,), (3,)), None
    london_c6                 = float, False, ((1,), ('ntyp',)), None
    london_rvdw               = float, False, ((1,), ('ntyp',)), None
    nr1                       = int,   False, (1,), None
    nr2                       = int,   False, (1,), None
    nr3                       = int,   False, (1,), None
    nr1s                      = int,   False, (1,), None
    nr2s                      = int,   False, (1,), None
    nr3s                      = int,   False, (1,), None
    nqx1                      = int,   False, (1,), None
    nqx2                      = int,   False, (1,), None
    nqx3                      = int,   False, (1,), None
    Hubbard_occ               = float, False, ((1,), ('ntyp,3',)), None
    starting_ns_eigenvalue    = float, False, ((1,), ('2*lmax+1,nspin or npol,ntyp',)), None
#end class SystemVariables

class ElectronsVariables(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &ELECTRONS input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    electron_maxstep  = int,   False, (1,), None
    exx_maxstep       = int,   False, (1,), None
    scf_must_converge = bool,  False, (1,), None
    conv_thr          = float, False, (1,), None
    adaptive_thr      = bool,  False, (1,), None
    conv_thr_init     = float, False, (1,), None
    conv_thr_multi    = float, False, (1,), None
    mixing_mode       = str,   False, (1,), ('plain', 'TF', 'local-TF')
    mixing_beta       = float, False, (1,), None
    mixing_ndim       = int,   False, (1,), None
    mixing_fixed_ns   = int,   False, (1,), None
    diagonalization   = str,   False, (1,), ('david', 'cg', 'ppcg', 'paro', 'ParO', 'rmm-davidson', 'rmm-paro')
    diago_thr_init    = float, False, (1,), None
    diago_cg_maxiter  = int,   False, (1,), None
    diago_david_ndim  = int,   False, (1,), None
    diago_rmm_ndim    = int,   False, (1,), None
    diago_rmm_conv    = bool,  False, (1,), None
    diago_gs_nblock   = int,   False, (1,), None
    diago_full_acc    = bool,  False, (1,), None
    efield            = float, False, (1,), None
    efield_phase      = str,   False, (1,), ('read', 'write', 'none')
    startingpot       = str,   False, (1,), ('atomic', 'file')
    startingwfc       = str,   False, (1,), ('atomic', 'atomic+random', 'random', 'file')
    tqr               = bool,  False, (1,), None
    real_space        = bool,  False, (1,), None
    efield_cart       = float, False, ((1,), (3,)), None
#end class ElectronsVariables

class IonsVariables(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &IONS input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    ion_positions     = str,   False, (1,), ('default', 'from_input')
    ion_velocities    = str,   False, (1,), ('default', 'from_input')
    ion_dynamics      = str,   False, (1,), ('bfgs', 'damp', 'fire', 'verlet', 'velocity-verlet', 'langevin', 'langevin-smc', 'bfgs', 'damp', 'beeman')
    pot_extrapolation = str,   False, (1,), ('none', 'atomic', 'first_order', 'second_order')
    wfc_extrapolation = str,   False, (1,), ('none', 'first_order', 'second_order')
    remove_rigid_rot  = bool,  False, (1,), None
    ion_temperature   = str,   False, (1,), ('rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'nose', 'berendsen', 'andersen', 'svr', 'initial', 'not_controlled')
    tempw             = float, False, (1,), None
    fnosep            = float, False, (1,), None
    nhpcl             = int,   False, (1,), None
    nhptyp            = int,   False, (1,), None
    ndega             = int,   False, (1,), None
    tolp              = float, False, (1,), None
    delta_t           = float, False, (1,), None
    nraise            = int,   False, (1,), None
    refold_pos        = bool,  False, (1,), None
    nhgrp             = int,   False, ((1,), ('ntyp',)), None
    fnhscl            = float, False, ((1,), ('ntyp',)), None
    upscale           = float, False, (1,), None
    bfgs_ndim         = int,   False, (1,), None
    tgdiis_step       = bool,  False, (1,), None
    trust_radius_max  = float, False, (1,), None
    trust_radius_min  = float, False, (1,), None
    trust_radius_ini  = float, False, (1,), None
    w_1               = float, False, (1,), None
    w_2               = float, False, (1,), None
    fire_alpha_init   = float, False, (1,), None
    fire_falpha       = float, False, (1,), None
    fire_nmin         = int,   False, (1,), None
    fire_f_inc        = float, False, (1,), None
    fire_f_dec        = float, False, (1,), None
    fire_dtmax        = float, False, (1,), None
#end class IonsVariables

class CellVariables(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &CELL input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    cell_dynamics  = str,   False, (1,), ('none', 'sd', 'damp-pr', 'damp-w', 'bfgs', 'none', 'pr', 'w')
    press          = float, False, (1,), None
    wmass          = float, False, (1,), None
    cell_factor    = float, False, (1,), None
    press_conv_thr = float, False, (1,), None
    cell_dofree    = str,   False, (1,), ('all', 'ibrav', 'a', 'b', 'c', 'fixa', 'fixb', 'fixc', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape', 'epitaxial_ab', 'epitaxial_ac', 'epitaxial_bc')
#end class CellVariables

class FCPVariables(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &FCP input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    fcp_mu           = float, True,  (1,), None
    fcp_dynamics     = str,   False, (1,), ('bfgs', 'newton', 'damp', 'lm', 'velocity-verlet', 'verlet')
    fcp_conv_thr     = float, False, (1,), None
    fcp_ndiis        = int,   False, (1,), None
    freeze_all_atoms = bool,  False, (1,), None
    fcp_mass         = float, False, (1,), None
    fcp_velocity     = float, False, (1,), None
    fcp_temperature  = str,   False, (1,), ('rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'berendsen', 'andersen', 'initial', 'not_controlled')
    fcp_tempw        = float, False, (1,), None
    fcp_tolp         = float, False, (1,), None
    fcp_delta_t      = float, False, (1,), None
    fcp_nraise       = int,   False, (1,), None
#end class FCPVariables

class RismVariables(NamelistParamDefinition, NamelistEnumBase):
    """Class representing all variables that belongs to the
    &RISM input namelist for ``pw.x``.

    Notes
    -----
    If you wish to add a new variable to the &CONTROL namelist, you
    can simply append the variable to the end of the definitions.
    The only required attributes that you must define are ``datatype``
    and ``default``, however it is ***strongly encouraged*** to specify
    the remaining information so that the various checks present
    throughout the rest of the code base can be utilized properly.
    """

    nsolv                 = int,   True,  (1,), None
    closure               = str,   False, (1,), ('kh', 'hnc')
    tempv                 = float, False, (1,), None
    ecutsolv              = float, False, (1,), None
    starting1d            = str,   False, (1,), ('zero', 'file', 'fix')
    starting3d            = str,   False, (1,), ('zero', 'file')
    smear1d               = float, False, (1,), None
    smear3d               = float, False, (1,), None
    rism1d_maxstep        = int,   False, (1,), None
    rism3d_maxstep        = int,   False, (1,), None
    rism1d_conv_thr       = float, False, (1,), None
    rism3d_conv_thr       = float, False, (1,), None
    mdiis1d_size          = int,   False, (1,), None
    mdiis3d_size          = int,   False, (1,), None
    mdiis1d_step          = float, False, (1,), None
    mdiis3d_step          = float, False, (1,), None
    rism1d_bond_width     = float, False, (1,), None
    rism1d_dielectric     = float, False, (1,), None
    rism1d_molesize       = float, False, (1,), None
    rism1d_nproc          = int,   False, (1,), None
    rism3d_conv_level     = float, False, (1,), None
    rism3d_planar_average = bool,  False, (1,), None
    laue_nfit             = int,   False, (1,), None
    laue_expand_right     = float, False, (1,), None
    laue_expand_left      = float, False, (1,), None
    laue_starting_right   = float, False, (1,), None
    laue_starting_left    = float, False, (1,), None
    laue_buffer_right     = float, False, (1,), None
    laue_buffer_left      = float, False, (1,), None
    laue_both_hands       = bool,  False, (1,), None
    laue_wall             = str,   False, (1,), ('none', 'auto', 'manual')
    laue_wall_z           = float, False, (1,), None
    laue_wall_rho         = float, False, (1,), None
    laue_wall_epsilon     = float, False, (1,), None
    laue_wall_sigma       = float, False, (1,), None
    laue_wall_lj6         = bool,  False, (1,), None
    solute_lj             = str,   False, ((1,), ('ntyp',)), None
    solute_epsilon        = float, False, ((1,), ('ntyp',)), None
    solute_sigma          = float, False, ((1,), ('ntyp',)), None
#end class RismVariables


class PWscfNamelistBase(ABC):
    """Abstract base class for real instances of PWscf namelists."""

    @abstractmethod
    def write(self) -> str: ...

    @abstractmethod
    def meets_requirements(self) -> bool: ...

    @abstractmethod
    def __setattr__(self, name, value): ...

    @property
    def value(self) -> PwscfInputType:
        return self._value

    @value.setter
    def value(self, val: PwscfInputType):
        if val is None and self.required:
            raise ValueError(
                f"Can not set {self.name} to None, it is required!"
            )
        elif not isinstance(val, self.datatype):
            return TypeError(
                f"Expected a value of type {self.datatype} for {self.name} but got {type(val)} instead!"
            )
        else:
            self._value = val


    def write_value(self) -> str:
        match self.datatype:
            case builtins.bool:
                if self._value:
                    return f"   {self.name:<15} = .true."
                else:
                    return f"   {self.name:<15} = .false."

            case builtins.str:
                return f"   {self.name:<15} = '{self.value}'"

            case builtins.int:
                return f"   {self.name:<15} = {self.value}"

            case builtins.float:
                return f"   {self.name:<15} = {self.value}"

            case builtins.list:
                string = ""
                if isinstance(self.value[0], str):
                    for index, val in enumerate(self.value):
                        string += f"   {self.name+f'({index+1})':<15} = '{val}'"
                else:
                    for index, val in enumerate(self.value):
                        string += f"   {self.name+f'({index+1})':<15} = {val}"
                return string
            case _:
                raise TypeError("Value is not a valid Pwscf Input Type!")
    #end def write_value
#end class PWscfNamelistBase



class PWscfCardBaseIO(ABC):
    """Abstract base for PWscf card read/write IO."""

    @abstractmethod
    def write(self, card: PWscfCardBase) -> str:
        """Write the card body (excluding header line)."""
        ...

    @abstractmethod
    def read(self, body: str, **kwargs) -> PWscfCardBase:
        """Read the card body (excluding header line). Returns the card instance."""
        ...


class PWscfCardBase(ABC):
    """Abstract base class for PWscf cards."""
    allowed_specifiers: frozenset[str] = frozenset()
    _io: dict[str, 'PWscfCardBaseIO'] = {}

    @abstractmethod
    def write(self) -> str: ...

    @classmethod
    @abstractmethod
    def read(cls, card: str): ...

#end class PWscfCardBase


class AtomicSpecies(PWscfCardBase):
    """Class for the ATOMIC_SPECIES card."""

    def __init__(
        self,
        elements: Sequence[str],
        pseudos: Sequence[str | PathLike],
    ):
        self.elements = list(elements)
        self.pseudos = [Path(ps) for ps in pseudos]
        #+ Will need to add masses, contingent on PR #5832
        self.masses = None


    @classmethod
    def from_physical_system(cls, system: PhysicalSystem, pseudos: str) -> Self:
        """Generate the information for the ATOMIC_SPECIES card from a
        Nexus ``PhysicalSystem``. Must have elements and be pseudized.
        """
        ...

_PWSCF_ARRAY_FORMAT = '{0:16.8f}' # could be imported from somewhere else to be more general


def _array_to_string(
    a: PwscfArrayInput,
    pad: str = '   ',
    fmt: str = _PWSCF_ARRAY_FORMAT,
    rowsep: str = '\n',
) -> str:
    """Format a numpy array for QE card output."""
    a = np.asarray(a)
    if a.ndim == 1:
        return pad + ' '.join(fmt.format(v) for v in a) + rowsep
    return rowsep.join(
        pad + ' '.join(fmt.format(v) for v in row) for row in a
    )

# Kpoints Card
class KPointsGammaIO(PWscfCardBaseIO):
    """IO for K_POINTS {gamma} - no additional lines."""

    def write(self, card: KPoints) -> str:
        return ''
    @staticmethod
    def read(lines: list[str]) -> KPoints:
        return KPoints(specifier='gamma')

class KPointsAutomaticIO(PWscfCardBaseIO):
    """IO for K_POINTS {automatic} - Monkhorst-Pack grid."""

    def write(self, card: KPoints) -> str:
        grid = np.asarray(card.grid, dtype=int)
        shift = np.asarray(card.shift, dtype=int)
        s = '   '
        s += ' '.join(str(g) for g in grid)
        s += ' '
        s += ' '.join(str(sv) for sv in shift)
        return s

    def read(self, lines: list[str]) -> KPoints:
        a = np.fromstring(lines[0], sep=' ')
        grid = a[0:3]
        shift = a[3:]
        return KPoints(specifier='automatic', grid=grid, shift=shift)


class KPointsExplicitIO(PWscfCardBaseIO):
    """IO for K_POINTS {tpiba,crystal,tpiba_b,crystal_b,tpiba_c,crystal_c}.

    The card must store kpoints in the coordinate system matching the specifier:
      - tpiba*: Cartesian coordinates in 2π/alat units
      - crystal*: Crystal (reciprocal lattice) coordinates
    The writer outputs the stored values as-is; no unit conversion is performed.
    """

    def write(self, card: KPoints) -> str:
        kpoints = np.asarray(card.kpoints)
        weights = np.asarray(card.weights)
        nkpoints = len(kpoints)
        a = np.empty((nkpoints, 4))
        a[:, 0:3] = kpoints
        a[:, 3] = weights
        return f'   {nkpoints}\n' + _array_to_string(a)

    def read(self, lines: list[str], specifier: str = 'crystal', **kwargs) -> KPoints:
        nkpoints = int(lines[0].strip())
        if len(lines) < 1 + nkpoints:
            raise ValueError(
                f"K_POINTS explicit: expected {nkpoints} data rows, got {len(lines) - 1}"
            )
        arr = np.loadtxt(lines[1:1 + nkpoints])
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        kpoints = arr[:, :3]
        weights = arr[:, 3] if arr.shape[1] > 3 else np.ones(len(kpoints))
        return KPoints(specifier=specifier, kpoints=kpoints, weights=weights)

class KPoints(PWscfCardBase):
    """K_POINTS card with strategy-pattern writers for specifier-dependent formats.

    Specifiers:
      - gamma: Gamma-point only
      - automatic: Monkhorst-Pack grid (nk1 nk2 nk3 sk1 sk2 sk3)
      - tpiba, crystal, tpiba_b, crystal_b, tpiba_c, crystal_c: require explicit list

    For explicit specifiers, kpoints must be in the matching coordinate system, no unit conversion is performed.
      - tpiba*: Cartesian, 2pi/alat units
      - crystal*: Crystal (reciprocal lattice) coordinates
    
    Use change_specifier() to convert between explicit specifiers
    Use from_physical_system() to create to inherit kpoints from a Nexus PhysicalSystem
    """

    _allowed_specifiers = frozenset({
        'gamma', 'automatic',
        'tpiba', 'crystal', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c',
    })
    _explicit_specifiers = frozenset({
        'tpiba', 'crystal', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c',
    })

    _io: dict[str, PWscfCardBaseIO] = {
        'gamma': KPointsGammaIO(),
        'automatic': KPointsAutomaticIO(),
        'tpiba': KPointsExplicitIO(),
        'crystal': KPointsExplicitIO(),
        'tpiba_b': KPointsExplicitIO(),
        'crystal_b': KPointsExplicitIO(),
        'tpiba_c': KPointsExplicitIO(),
        'crystal_c': KPointsExplicitIO(),
    }

    def __init__(
        self,
        specifier: str = 'gamma',
        grid: tuple[int, int, int] | None = None,
        shift: tuple[int, int, int] | None = None,
        kpoints: PwscfArrayInput | None = None,
        weights: PwscfArrayInput | None = None,
    ):
        spec = specifier.lower() if specifier else 'gamma'
        if spec not in self._allowed_specifiers:
            raise ValueError(
                f"K_POINTS specifier '{specifier}' not recognized. "
                f"Allowed: {sorted(self._allowed_specifiers)}"
            )
        self.specifier = spec

        if spec == 'gamma':
            self.grid = None
            self.shift = None
            self.kpoints = np.empty((0, 3))
            self.weights = np.empty((0,))
        elif spec == 'automatic':
            if grid is None or shift is None:
                raise ValueError("grid and shift required for specifier 'automatic'")
            self.grid = tuple(int(g) for g in grid)
            self.shift = tuple(int(sv) for sv in shift)
            self.kpoints = np.empty((0, 3))
            self.weights = np.empty((0,))
        else:
            if kpoints is None or weights is None:
                raise ValueError(
                    f"kpoints and weights required for specifier '{spec}'"
                )
            kpoints = np.asarray(kpoints)
            weights = np.asarray(weights)
            if len(kpoints) != len(weights):
                raise ValueError(
                    f"kpoints and weights length mismatch: {len(kpoints)} vs {len(weights)}"
                )
            self.grid = None
            self.shift = None
            self.kpoints = kpoints
            self.weights = weights

    def write(self) -> str:
        header = f"K_POINTS {self.specifier}\n"
        body = self._io[self.specifier].write(self)
        return header + body + ('\n' if body else '')

    @classmethod
    def read(cls, lines: list[str], specifier: str = 'crystal') -> Self:
        if specifier in cls._explicit_specifiers:
            return cls._io[specifier].read(lines, specifier=specifier.lower())
        else:
            return cls._io[specifier].read(lines)
    
    def change_specifier(
        self,
        new_specifier: str,
        scale: float,
        kaxes: PwscfArrayInput,
    ) -> None:
        """Convert k-points to a different coordinate specifier.

        Reuses the legacy logic. Converts via cartesian reciprocal space as intermediate.

        Parameters
        ----------
        new_specifier : str
            One of: gamma, automatic, tpiba, crystal (tpiba_b, crystal_b not yet implemented).
        scale : float
            Lattice parameter alat (e.g. celldm(1) in Bohr).
        kaxes : ndarray
            Reciprocal lattice vectors, shape (3,3). kaxes = 2π * inv(axes).T
            for axes in Bohr.
        """
        pi = np.pi
        spec = self.specifier
        new_spec = new_specifier.lower()

        # Convert from current specifier to Cartesian reciprocal
        if spec in ('tpiba', ''):
            kpoints = self.kpoints * (2 * pi) / scale
        elif spec == 'gamma':
            kpoints = np.array([[0.0, 0.0, 0.0]])
        elif spec == 'crystal':
            kpoints = np.dot(self.kpoints, kaxes)
        elif spec == 'automatic':
            grid = np.array(self.grid, dtype=int)
            shift = 0.5 * np.array(self.shift)
            kpoints = kmesh(kaxes, grid, shift)
            # automatic -> explicit: need weights (uniform for MP mesh)
            nk = len(kpoints)
            self.weights = np.full(nk, 1.0 / nk)
        elif spec in ('tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c'):
            raise NotImplementedError(
                f"change_specifier from '{spec}' not yet implemented"
            )
        else:
            raise ValueError(
                f"Invalid current specifier '{spec}'. "
                f"Valid: tpiba, gamma, crystal, automatic, tpiba_b, crystal_b"
            )

        # Convert from Cartesian to new specifier
        if new_spec in ('tpiba', 'tpiba_b', 'tpiba_c'):
            kpoints = kpoints / ((2 * pi) / scale)
        elif new_spec == 'gamma':
            kpoints = np.array([[0.0, 0.0, 0.0]])
        elif new_spec in ('crystal', 'crystal_b', 'crystal_c'):
            kpoints = np.dot(kpoints, inv(kaxes))
        elif new_spec == 'automatic':
            if spec != 'automatic':
                raise ValueError(
                    "cannot map arbitrary k-points into a Monkhorst-Pack mesh"
                )
        else:
            raise ValueError(
                f"Invalid new specifier '{new_specifier}'. "
                f"Valid: {sorted(self._explicit_specifiers)}"
            )

        self.kpoints = np.asarray(kpoints)
        self.specifier = new_spec

    @classmethod
    def from_gamma(cls) -> Self:
        """Create K_POINTS for Gamma-point-only calculation."""
        return cls(specifier='gamma')

    @classmethod
    def from_automatic(
        cls,
        kgrid: Sequence[int],
        kshift: Sequence[int] = (0, 0, 0),
    ) -> Self:
        """Create K_POINTS using uniform grid and shift."""
        kgrid = tuple(kgrid)
        kshift = tuple(kshift)
        if len(kgrid) != 3 or len(kshift) != 3:
            raise ValueError("kgrid and kshift must have 3 elements each")
        return cls(specifier='automatic', grid=kgrid, shift=kshift)

    @classmethod
    def from_explicit(
        cls,
        kpoints: PwscfArrayInput,
        weights: PwscfArrayInput,
        specifier: str = 'crystal',
    ) -> Self:
        """Create K_POINTS from an explicit k-point list."""
        spec = specifier.lower()
        if spec not in cls._explicit_specifiers:
            raise ValueError(f"Specifier '{spec}' must be one of {sorted(cls.explicit_specifiers)}")
        kpoints = np.asarray(kpoints)
        weights = np.asarray(weights)
        return cls(
            specifier=spec,
            kpoints=kpoints,
            weights=weights,
        )

    @classmethod
    def from_physical_system(
        cls,
        system: PhysicalSystem,
        specifier: str = 'crystal',
    ) -> Self:
        """Create K_POINTS from a Nexus PhysicalSystem.

        Uses structure.kpoints (Cartesian reciprocal) and structure.kweights.
        Builds in crystal coordinates first, then uses change_specifier to convert
        if the requested specifier is tpiba or another format.
        """
        structure = system.structure
        kpoints = structure.kpoints
        kweights = structure.kweights

        if len(kpoints) == 0 or len(kweights) == 0:
            raise ValueError("PhysicalSystem has no k-points or k-weights")

        spec = specifier.lower()
        if spec not in cls._explicit_specifiers:
            raise ValueError(f"Specifier '{spec}' must be one of {sorted(cls._explicit_specifiers)}")

        # Build in crystal coordinates (from structure.kpoints_unit)
        card = cls.from_explicit(
            kpoints=structure.kpoints_unit(),
            weights=kweights,
            specifier='crystal',
        )
        if spec in ('tpiba', 'tpiba_b', 'tpiba_c'):
            scale = structure.scale
            kaxes = structure.kaxes
            card.change_specifier('tpiba', scale, kaxes)
            card.specifier = spec  # preserve tpiba_b, tpiba_c if requested
        elif spec in ('crystal_b', 'crystal_c'):
            card.specifier = spec  # same coords as crystal, different specifier
        return card

# Hubbard Card
@dataclass(frozen=True)
class HubbardOnSite:
    """On site Hubbard parameters.
    Format: param label value
    """
    _allowed_params: ClassVar[frozenset[str]] = frozenset({'U', 'ALPHA', 'J0', 'J', 'B', 'E2', 'E3'})
    param: str
    label: str  # e.g. "Ni-3d"
    value: float

    def __post_init__(self):
        if self.param not in self._allowed_params:
            raise ValueError(f"Invalid on-site parameter: {self.param}")

@dataclass(frozen=True)
class HubbardOrbitalResolved(HubbardOnSite):
    """On site Hubbard parameters for orbital-resolved calculations.
    Format: param label value orb1 orb2 ... orbn, n should be in the range [1, 2*l+1]
    """
    param: str
    orbitals: tuple[int, ...]  # orbital indices in manifold

@dataclass(frozen=True)
class HubbardInterSite:
    """Inter-site Hubbard parameters.
    Format: V label1 label2 ind1 ind2 value
    """
    _allowed_params: ClassVar[frozenset[str]] = frozenset({'V'})
    param: str
    label1: str  # e.g. "Co-3d"
    label2: str  # e.g. "O-2p"
    ind1: int  # atom index (1-based, ATOMIC_POSITIONS order)
    ind2: int
    value: float

    def __post_init__(self):
        if self.param not in self._allowed_params:
            raise ValueError(f"Invalid inter-site parameter: {self.param}")


HubbardEntry = HubbardOnSite | HubbardInterSite | HubbardOrbitalResolved


class HubbardIO(PWscfCardBaseIO):
    """IO for HUBBARD card. Body format is the same for all specifiers."""

    _param_width = 5
    _label_width = 5
    _value_fmt = '.3f'

    def write(self, card: 'Hubbard') -> str:
        entries = card.entries
        on_site = [e for e in entries if isinstance(e, HubbardOnSite) and not isinstance(e, HubbardOrbitalResolved)]
        orbital = [e for e in entries if isinstance(e, HubbardOrbitalResolved)]
        inter_site = [e for e in entries if isinstance(e, HubbardInterSite)]
        ordered = on_site + orbital + inter_site

        parts = []
        for e in ordered:
            p = f"{e.param:<{self._param_width}}"
            v = format(e.value, self._value_fmt)
            if isinstance(e, HubbardOrbitalResolved):
                lb = f"{e.label:<{self._label_width}}"
                orb_str = ' '.join(str(o) for o in e.orbitals)
                parts.append(f"{p} {lb} {v}  {orb_str}")
            elif isinstance(e, HubbardOnSite):
                lb = f"{e.label:<{self._label_width}}"
                parts.append(f"{p} {lb} {v}")
            elif isinstance(e, HubbardInterSite):
                l1 = f"{e.label1:<{self._label_width}}"
                l2 = f"{e.label2:<{self._label_width}}"
                parts.append(f"{p} {l1} {l2} {e.ind1} {e.ind2} {v}")
        return '\n'.join(parts)
    
    def read(self, lines: list[str], specifier: str = 'atomic', **kwargs) -> Hubbard:
        entries: list[HubbardEntry] = []
        on_site_params = HubbardOnSite._allowed_params
        for line in lines:
            tokens = line.split()
            if len(tokens) == 3:
                param, label, val = tokens[0], tokens[1], float(tokens[2])
                if param in on_site_params:
                    entries.append(HubbardOnSite(param, label, val))
            elif len(tokens) == 6 and tokens[0].upper() == 'V':
                entries.append(HubbardInterSite(
                    tokens[0], tokens[1], tokens[2],
                    int(tokens[3]), int(tokens[4]), float(tokens[5]),
                ))
            elif len(tokens) >= 4 and tokens[0] in on_site_params:
                param, label, val = tokens[0], tokens[1], float(tokens[2])
                orbitals = tuple(int(t) for t in tokens[3:])
                entries.append(HubbardOrbitalResolved(param, label, val, orbitals))
        return Hubbard(entries=entries, specifier=specifier)


class Hubbard(PWscfCardBase):
    """HUBBARD card with strategy-pattern writers for line types 1, 2, 3.

    Specifiers (projection type): atomic, ortho-atomic, norm-atomic, wf, pseudo.
    Line types:
      1: param label value (on-site shell-averaged)
      2: V label1 label2 ind1 ind2 value (inter-site)
      3: param label value orb1 orb2 ... (on-site orbital-resolved)
    """

    allowed_specifiers = frozenset({
        'atomic', 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo',
    })
    default_specifier = 'atomic'
    _io = HubbardIO()

    def __init__(
        self,
        entries: Sequence[HubbardEntry],
        specifier: str = 'atomic',
    ):
        spec = specifier.lower() if specifier else self.default_specifier
        if spec not in self.allowed_specifiers:
            raise ValueError(
                f"HUBBARD specifier '{specifier}' not recognized. "
                f"Allowed: {sorted(self.allowed_specifiers)}"
            )
        self.specifier = spec
        self.entries = list(entries)

    def check_compatible(self, physical_system: PhysicalSystem) -> bool:
        """Check if the HUBBARD card is compatible with the PhysicalSystem.

        Validates that:
        - Species in Hubbard labels (e.g. 'Ni' from 'Ni-3d') exist in structure.elem
        - Inter-site indices are not checked, as these are typically provided by hp.x 
        """
        if len(self.entries) == 0:
            return True

        valid_species = set(physical_system.structure.elem)
        def _extract_species(label: str) -> str:
            """Extract species from Hubbard label (e.g. 'Ni-3d' -> 'Ni')."""
            return label.split('-')[0] if '-' in label else label

        for e in self.entries:
            if isinstance(e, (HubbardOnSite, HubbardOrbitalResolved)):
                sp = _extract_species(e.label)
                if sp not in valid_species:
                    raise ValueError(
                        f"Hubbard label '{e.label}' references unknown species '{sp}'. "
                        f"Valid species: {sorted(valid_species)}"
                    )
            elif isinstance(e, HubbardInterSite):
                for lbl in (e.label1, e.label2):
                    sp = _extract_species(lbl)
                    if sp not in valid_species:
                        raise ValueError(
                            f"Hubbard inter-site label '{lbl}' references unknown species '{sp}'. "
                            f"Valid species: {sorted(valid_species)}"
                        )
        return True

    def write(self) -> str:
        header = f"HUBBARD {self.specifier}\n"
        body = self._io.write(self)
        return header + body + '\n'
    
    @classmethod
    def read(cls, lines: list[str], specifier: str = 'atomic') -> Self:
        return cls._io.read(lines, specifier=specifier)

    @classmethod
    def from_entries(
        cls,
        entries: Sequence[HubbardEntry],
        specifier: str = 'atomic',
    ) -> Self:
        """Create HUBBARD card from a sequence of HubbardEntry objects."""
        return cls(entries=entries, specifier=specifier)

    @classmethod
    def from_records(
        cls,
        records: Sequence[dict],
        specifier: str = 'atomic',
    ) -> Self:
        """Create HUBBARD card from a list of dicts.

        Options:
          1. on_site:  {type: 'on_site', param, label, value}
          2. orbital:  {type: 'orbital', param, label, value, orbitals: tuple|list}
          3. inter_site: {type: 'inter_site', param, label1, label2, ind1, ind2, value}
        """
        entries: list[HubbardEntry] = []
        for r in records:
            kind = r.get('type')
            if kind == 'on_site':
                entries.append(HubbardOnSite(
                    r['param'], r['label'], float(r['value'])
                ))
            elif kind == 'orbital':
                orb = r['orbitals']
                entries.append(HubbardOrbitalResolved(
                    r['param'], r['label'], float(r['value']),
                    tuple(orb) if not isinstance(orb, tuple) else orb,
                ))
            elif kind == 'inter_site':
                entries.append(HubbardInterSite(
                    r['param'],
                    r['label1'], r['label2'],
                    int(r['ind1']), int(r['ind2']),
                    float(r['value']),
                ))
            else:
                raise ValueError(
                    f"Hubbard record type '{kind}' not recognized. "
                    "Use 'on_site', 'orbital', or 'inter_site'."
                )
        return cls(entries=entries, specifier=specifier)

    @classmethod
    def from_hp_output(cls, text: str) -> Self:
        """Create HUBBARD card from hp.x HUBBARD.dat output.

        Parses the block after '# Copy this data in the pw.x input file...'
        """
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        specifier = cls.default_specifier
        entries: list[HubbardEntry] = []
        on_site_params = HubbardOnSite._allowed_params

        for line in lines:
            if line.startswith('#'):
                continue
            tokens = line.split()
            if not tokens:
                continue
            if tokens[0].upper() == 'HUBBARD':
                for t in tokens[1:]:
                    s = t.strip('{}()')
                    if s and s.lower() in cls.allowed_specifiers:
                        specifier = s.lower()
                        break
                continue
            if len(tokens) == 3:
                param, label, val = tokens[0], tokens[1], float(tokens[2])
                if param in on_site_params:
                    entries.append(HubbardOnSite(param, label, val))
            elif len(tokens) == 6 and tokens[0].upper() == 'V':
                entries.append(HubbardInterSite(
                    tokens[0],
                    tokens[1], tokens[2],
                    int(tokens[3]), int(tokens[4]),
                    float(tokens[5]),
                ))
            elif len(tokens) >= 4 and tokens[0] in on_site_params:
                param, label, val = tokens[0], tokens[1], float(tokens[2])
                orbitals = tuple(int(t) for t in tokens[3:])
                entries.append(HubbardOrbitalResolved(param, label, val, orbitals))
        return cls(entries=entries, specifier=specifier)

    def from_hubbard_dat_file(self, file: PathLike) -> Self:
        """Create HUBBARD card from a hubbard.dat file."""
        with Path(file).open() as f:
            return self.from_hp_output(f.read())