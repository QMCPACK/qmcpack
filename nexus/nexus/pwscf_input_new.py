from __future__ import annotations
import os
import numpy as np
import numpy.typing as npt
import builtins
from dataclasses import dataclass
from abc import ABC, abstractmethod
from enum import Enum
from collections.abc import Sequence


type PwscfInputType = (
    str
    | bool
    | int
    | float
    | Sequence
)


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

    calculation               = str,   False, (1,), ('scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md')
    title                     = str,   False, (1,), None
    verbosity                 = str,   False, (1,), ('high', 'low')
    restart_mode              = str,   False, (1,), ('from_scratch', 'restart')
    wf_collect                = bool,  False, (1,), None
    nstep                     = int,   False, (1,), None
    iprint                    = int,   False, (1,), None
    tstress                   = bool,  False, (1,), None
    tprnfor                   = bool,  False, (1,), None
    dt                        = float, False, (1,), None
    outdir                    = str,   False, (1,), None
    wfcdir                    = str,   False, (1,), None
    prefix                    = str,   False, (1,), None
    lkpoint_dir               = bool,  False, (1,), None
    max_seconds               = float, False, (1,), None
    etot_conv_thr             = float, False, (1,), None
    forc_conv_thr             = float, False, (1,), None
    disk_io                   = str,   False, (1,), ('high', 'medium', 'low', 'nowf', 'minimal', 'none')
    pseudo_dir                = str,   False, (1,), None
    tefield                   = bool,  False, (1,), None
    dipfield                  = bool,  False, (1,), None
    lelfield                  = bool,  False, (1,), None
    nberrycyc                 = int,   False, (1,), None
    lorbm                     = bool,  False, (1,), None
    lberry                    = bool,  False, (1,), None
    gdir                      = int,   False, (1,), None
    nppstr                    = int,   False, (1,), None
    gate                      = bool,  False, (1,), None
    twochem                   = bool,  False, (1,), None
    lfcp                      = bool,  False, (1,), None
    trism                     = bool,  False, (1,), None
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

    electron_maxstep          = int,   False, (1,), None
    exx_maxstep               = int,   False, (1,), None
    scf_must_converge         = bool,  False, (1,), None
    conv_thr                  = float, False, (1,), None
    adaptive_thr              = bool,  False, (1,), None
    conv_thr_init             = float, False, (1,), None
    conv_thr_multi            = float, False, (1,), None
    mixing_mode               = str,   False, (1,), ('plain', 'TF', 'local-TF')
    mixing_beta               = float, False, (1,), None
    mixing_ndim               = int,   False, (1,), None
    mixing_fixed_ns           = int,   False, (1,), None
    diagonalization           = str,   False, (1,), ('david', 'cg', 'ppcg', 'paro', 'ParO', 'rmm-davidson', 'rmm-paro')
    diago_thr_init            = float, False, (1,), None
    diago_cg_maxiter          = int,   False, (1,), None
    diago_david_ndim          = int,   False, (1,), None
    diago_rmm_ndim            = int,   False, (1,), None
    diago_rmm_conv            = bool,  False, (1,), None
    diago_gs_nblock           = int,   False, (1,), None
    diago_full_acc            = bool,  False, (1,), None
    efield                    = float, False, (1,), None
    efield_phase              = str,   False, (1,), ('read', 'write', 'none')
    startingpot               = str,   False, (1,), ('atomic', 'file')
    startingwfc               = str,   False, (1,), ('atomic', 'atomic+random', 'random', 'file')
    tqr                       = bool,  False, (1,), None
    real_space                = bool,  False, (1,), None
    efield_cart               = float, False, ((1,), (3,)), None
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

    ion_positions             = str,   False, (1,), ('default', 'from_input')
    ion_velocities            = str,   False, (1,), ('default', 'from_input')
    ion_dynamics              = str,   False, (1,), ('bfgs', 'damp', 'fire', 'verlet', 'velocity-verlet', 'langevin', 'langevin-smc', 'bfgs', 'damp', 'beeman')
    pot_extrapolation         = str,   False, (1,), ('none', 'atomic', 'first_order', 'second_order')
    wfc_extrapolation         = str,   False, (1,), ('none', 'first_order', 'second_order')
    remove_rigid_rot          = bool,  False, (1,), None
    ion_temperature           = str,   False, (1,), ('rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'nose', 'berendsen', 'andersen', 'svr', 'initial', 'not_controlled')
    tempw                     = float, False, (1,), None
    fnosep                    = float, False, (1,), None
    nhpcl                     = int,   False, (1,), None
    nhptyp                    = int,   False, (1,), None
    ndega                     = int,   False, (1,), None
    tolp                      = float, False, (1,), None
    delta_t                   = float, False, (1,), None
    nraise                    = int,   False, (1,), None
    refold_pos                = bool,  False, (1,), None
    nhgrp                     = int,   False, ((1,), ('ntyp',)), None
    fnhscl                    = float, False, ((1,), ('ntyp',)), None
    upscale                   = float, False, (1,), None
    bfgs_ndim                 = int,   False, (1,), None
    tgdiis_step               = bool,  False, (1,), None
    trust_radius_max          = float, False, (1,), None
    trust_radius_min          = float, False, (1,), None
    trust_radius_ini          = float, False, (1,), None
    w_1                       = float, False, (1,), None
    w_2                       = float, False, (1,), None
    fire_alpha_init           = float, False, (1,), None
    fire_falpha               = float, False, (1,), None
    fire_nmin                 = int,   False, (1,), None
    fire_f_inc                = float, False, (1,), None
    fire_f_dec                = float, False, (1,), None
    fire_dtmax                = float, False, (1,), None
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

    cell_dynamics             = str,   False, (1,), ('none', 'sd', 'damp-pr', 'damp-w', 'bfgs', 'none', 'pr', 'w')
    press                     = float, False, (1,), None
    wmass                     = float, False, (1,), None
    cell_factor               = float, False, (1,), None
    press_conv_thr            = float, False, (1,), None
    cell_dofree               = str,   False, (1,), ('all', 'ibrav', 'a', 'b', 'c', 'fixa', 'fixb', 'fixc', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape', 'epitaxial_ab', 'epitaxial_ac', 'epitaxial_bc')
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

    fcp_mu                    = float, True,  (1,), None
    fcp_dynamics              = str,   False, (1,), ('bfgs', 'newton', 'damp', 'lm', 'velocity-verlet', 'verlet')
    fcp_conv_thr              = float, False, (1,), None
    fcp_ndiis                 = int,   False, (1,), None
    freeze_all_atoms          = bool,  False, (1,), None
    fcp_mass                  = float, False, (1,), None
    fcp_velocity              = float, False, (1,), None
    fcp_temperature           = str,   False, (1,), ('rescaling', 'rescale-v', 'rescale-T', 'reduce-T', 'berendsen', 'andersen', 'initial', 'not_controlled')
    fcp_tempw                 = float, False, (1,), None
    fcp_tolp                  = float, False, (1,), None
    fcp_delta_t               = float, False, (1,), None
    fcp_nraise                = int,   False, (1,), None
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

    nsolv                     = int,   True,  (1,), None
    closure                   = str,   False, (1,), ('kh', 'hnc')
    tempv                     = float, False, (1,), None
    ecutsolv                  = float, False, (1,), None
    starting1d                = str,   False, (1,), ('zero', 'file', 'fix')
    starting3d                = str,   False, (1,), ('zero', 'file')
    smear1d                   = float, False, (1,), None
    smear3d                   = float, False, (1,), None
    rism1d_maxstep            = int,   False, (1,), None
    rism3d_maxstep            = int,   False, (1,), None
    rism1d_conv_thr           = float, False, (1,), None
    rism3d_conv_thr           = float, False, (1,), None
    mdiis1d_size              = int,   False, (1,), None
    mdiis3d_size              = int,   False, (1,), None
    mdiis1d_step              = float, False, (1,), None
    mdiis3d_step              = float, False, (1,), None
    rism1d_bond_width         = float, False, (1,), None
    rism1d_dielectric         = float, False, (1,), None
    rism1d_molesize           = float, False, (1,), None
    rism1d_nproc              = int,   False, (1,), None
    rism3d_conv_level         = float, False, (1,), None
    rism3d_planar_average     = bool,  False, (1,), None
    laue_nfit                 = int,   False, (1,), None
    laue_expand_right         = float, False, (1,), None
    laue_expand_left          = float, False, (1,), None
    laue_starting_right       = float, False, (1,), None
    laue_starting_left        = float, False, (1,), None
    laue_buffer_right         = float, False, (1,), None
    laue_buffer_left          = float, False, (1,), None
    laue_both_hands           = bool,  False, (1,), None
    laue_wall                 = str,   False, (1,), ('none', 'auto', 'manual')
    laue_wall_z               = float, False, (1,), None
    laue_wall_rho             = float, False, (1,), None
    laue_wall_epsilon         = float, False, (1,), None
    laue_wall_sigma           = float, False, (1,), None
    laue_wall_lj6             = bool,  False, (1,), None
    solute_lj                 = str,   False, ((1,), ('ntyp',)), None
    solute_epsilon            = float, False, ((1,), ('ntyp',)), None
    solute_sigma              = float, False, ((1,), ('ntyp',)), None
#end class RismVariables


class PWscfNamelist(ABC):
    """Abstract base class for real instances of PWscf namelists."""

    @abstractmethod
    def write_card(self) -> str:
        pass

    @abstractmethod
    def meets_requirements(self) -> bool:
        pass

    @abstractmethod
    def __setattr__(self, name, value):
        pass


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
#end class PWscfNamelist
