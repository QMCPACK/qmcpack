##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf_analyzer.py                                                 #
#    Supports data analysis for PWSCF output.  Can handle log file   #
#    and XML output.                                                 #
#                                                                    #
#  Content summary:                                                  #
#    PwscfOutData                                                    #
#      Reads and stores physical data from PWSCF log output.         #
#                                                                    #
#    PwscfXmlData                                                    #
#      Reads and stores schema-based PWSCF XML output.               #
#                                                                    #
#    Pw2CasinoAnalyzer                                               #
#      Reads and stores physical data from PW2CASINO output.         #
#                                                                    #
#    PwscfAnalyzer                                                   #
#      SimulationAnalyzer class for PWSCF.                           #
#      Coordinates text and XML output analysis.                     #
#      Can also read data-file.xml.  See pwscf_data_reader.py.       #
#                                                                    #
#====================================================================#


import os
import re
import xml.etree.ElementTree as ET
from copy import deepcopy
from glob import glob
from types import MappingProxyType
import numpy as np
from .developer import DevBase,obj,dotdict
from .unit_converter import convert
from .numerics import simstats, simplestats
from .simulation import SimulationAnalyzer, Simulation
from .structure import Structure, get_kpath
from .pwscf_input import PwscfInput
from .pwscf_data_reader import read_qexml
from .utilities import path_string
from . import numpy_extensions as npe


# Match one complete decimal or scientific-notation number as PWSCF writes it.
# Examples include ``-168.12345678``, ``.5000000``, ``6.3E-09``, and
# ``-1.250D+02``.  The pattern intentionally excludes nonnumeric XML values
# such as ``true``, non-finite spellings such as ``NaN``, and incomplete
# exponents such as ``1.0E``.  It has no capturing groups so it can be safely
# embedded in each of the field-specific expressions below.
number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


# Match the numeric prefix of a PWSCF band or occupation line.  PWSCF usually
# separates values with spaces (``-5.3120  1.2470  4.9810``), but fixed-width
# output can join a negative value to its predecessor (``0.0488-0.0345``).
# A trailing diagnostic is permitted (``1.0 2.0  convergence achieved``),
# while a line beginning with a label (``bands: 1.0 2.0``) is rejected.
leading_number_list_pattern = (
    rf'^\s*(?P<values>{number_pattern}'
    rf'(?:(?:\s+|(?=[+-])){number_pattern})*)'
    )


# Match an entire whitespace-separated numeric text field.  This is tailored
# to schema XML arrays such as ``0.0 0.5 -0.5`` and ``1.0D+00\n0.0D+00``.
# Unlike the band-line pattern, it rejects adjacent values (``0.5-0.5``),
# units (``1.0 Ry``), labels (``weight=1.0``), and trailing prose.
numeric_text_pattern = (
    rf'^\s*(?P<values>{number_pattern}(?:\s+{number_pattern})*)\s*$'
    )


# Match the first three entries of a whitespace-separated vector row.  This
# accepts realistic CELL_PARAMETERS rows such as ``1.0 0.0 0.0`` and
# ``-2.5D-01  .500000  0.``; extra printed columns or comments are harmless.
# It rejects a short row (``1.0 0.0``) and concatenated entries
# (``1.0-0.5 0.0``), which are not valid PWSCF cell-vector formatting.
vector3_pattern = (
    rf'^\s*(?P<x>{number_pattern})\s+(?P<y>{number_pattern})\s+'
    rf'(?P<z>{number_pattern})(?=\s|$)'
    )


# Match one entry in each Cartesian or crystal k-point table.  Examples are
# ``k( 1) = ( 0.0 0.0 0.0), wk = 0.25`` and
# ``k(12)=(.5 -.5 5.0D-1), wk=1.0D+00``.  The required parentheses, comma,
# ``wk`` label, three coordinates, and weight prevent unrelated vectors,
# incomplete k-points, and malformed ``weight =`` variants from matching.
kpoint_table_pattern = (
    rf'\bk\(\s*\d+\s*\)\s*=\s*\(\s*'
    rf'(?P<kx>{number_pattern})\s+(?P<ky>{number_pattern})\s+'
    rf'(?P<kz>{number_pattern})\s*\)\s*,\s*wk\s*=\s*'
    rf'(?P<weight>{number_pattern})(?=\s|$)'
    )


# Match the k-vector in a ``bands (ev)`` header.  Normal output such as
# ``k = 0.0000 0.0000 -0.7071 (2138 PWs)`` and fixed-width output such as
# ``k = 0.0488-0.0345 0.0345 (42052 PWs)`` are both accepted.  Requiring
# exactly three coordinates followed by ``(`` rejects truncated vectors and
# labels such as ``k-point =`` that belong to other output sections.
band_kpoint_pattern = (
    rf'\bk\s*=\s*(?P<kx>{number_pattern})(?:\s+|(?=[+-]))'
    rf'(?P<ky>{number_pattern})(?:\s+|(?=[+-]))'
    rf'(?P<kz>{number_pattern})\s*\('
    )


# Match a PWSCF total-energy assignment including its Rydberg unit.  Both SCF
# iteration lines (``total energy = -168.12345678 Ry``) and completed lines
# (``! total energy = -1.0D+02 Ry``) match.  Lines in eV, missing ``=``, or
# reporting ``one-electron contribution`` do not match this field.
total_energy_pattern = (
    rf'\btotal\s+energy\s*=\s*(?P<energy>{number_pattern})\s+Ry\b'
    )


# Match the accuracy estimate immediately associated with an SCF iteration.
# Real forms include ``estimated scf accuracy < 6.3E-09 Ry`` and
# ``estimated scf accuracy = 1.2D-10 Ry``.  A missing comparison operator,
# another unit, or a shortened ``estimated accuracy`` label is rejected.
scf_accuracy_pattern = (
    rf'\bestimated\s+scf\s+accuracy\s*[<=>]\s*'
    rf'(?P<accuracy>{number_pattern})\s+Ry\b'
    )


# Match the ionic kinetic energy printed during molecular dynamics.  It
# accepts ``kinetic energy = 0.00285013 Ry`` and the parenthesized-label form
# ``kinetic energy (Ekin) = 2.85D-03 Ry``.  Temperature lines, eV values, and
# a bare ``Ekin =`` summary are intentionally handled by other patterns.
kinetic_energy_pattern = (
    rf'\bkinetic\s+energy(?:\s*\(Ekin\))?\s*=\s*'
    rf'(?P<kinetic_energy>{number_pattern})\s+Ry\b'
    )


# Match the instantaneous ionic temperature, for example
# ``temperature = 300.00000000 K`` or ``temperature=2.5D+02 K``.
# It rejects target-temperature prose, a missing Kelvin unit, and the compact
# ``Ekin ... T =`` summary that is parsed as a coupled record below.
temperature_pattern = (
    rf'^\s*temperature\s*=\s*(?P<temperature>{number_pattern})\s+K\b'
    )


# Match the coupled kinetic-energy/temperature record used by variable-cell
# dynamics.  Examples include ``Ekin = .00285 Ry T = 300.0 K Etot = ...`` and
# ``Ekin=2.8D-03 Ry  T=0.0 K``.  A missing Ry or K unit, reversed field order,
# or an isolated Ekin or T value is rejected so partial records are not used.
ekin_temperature_pattern = (
    rf'^\s*Ekin\s*=\s*(?P<kinetic_energy>{number_pattern})\s+Ry\s+'
    rf'T\s*=\s*(?P<temperature>{number_pattern})\s+K\b'
    )


# Match the dynamics time assignment in either standard MD output
# (``time = 0.50000 pico-seconds``) or the ``Entering Dynamics`` line
# (``Entering Dynamics; it = 1 time = 0.00000 pico-seconds``).  Requiring a
# time unit avoids taking the integer from ``it =`` or an unrelated timer.
md_time_pattern = (
    rf'\btime\s*=\s*(?P<time>{number_pattern})\s*(?:pico-seconds|ps)\b'
    )


# Match pressure from a PWSCF total-stress heading.  It accepts forms such as
# ``total stress (Ry/bohr**3) (kbar) P= -170.96`` and ``P = 1.2D+03``.
# Lowercase ``p=``, the word ``pressure``, and a missing assignment operator
# are rejected because they do not have the exact PWSCF field label.
pressure_pattern = rf'\bP\s*=\s*(?P<pressure>{number_pattern})(?=\s|$)'


# Match the unit-cell volume assignment, including realistic output such as
# ``unit-cell volume = 380.6210 (a.u.)^3`` and a scientific-notation value.
# A ``new unit-cell volume:`` message, a missing value, or ``cell volume`` is
# rejected so only the canonical PWSCF result line updates the volume.
volume_pattern = (
    rf'\bunit-cell\s+volume\s*=\s*(?P<volume>{number_pattern})(?=\s|$)'
    )


# Match the optional lattice scale in CELL_PARAMETERS headers.  Examples are
# ``CELL_PARAMETERS (alat= 10.20)`` and ``CELL_PARAMETERS (alat = 1.0D+01)``.
# The expression rejects ``celldm(1)=`` and ``alat:`` variants, which do not
# carry the same local scaling semantics for this structure block.
alat_pattern = rf'\balat\s*=\s*(?P<alat>{number_pattern})(?=\s|\))'


# Match one atomic-force row.  Representative lines are
# ``atom 1 type 2 force = -0.001 0.002 0.000`` and
# ``atom 12 type 1 force=1.0D-03 -2.0D-03 .0``.  The atom/type indices and all
# three components are required; total-force summaries and short rows fail.
atomic_force_pattern = (
    rf'\batom\s+(?P<atom>\d+)\s+type\s+(?P<type>\d+)\s+force\s*=\s*'
    rf'(?P<fx>{number_pattern})\s+(?P<fy>{number_pattern})\s+'
    rf'(?P<fz>{number_pattern})(?=\s|$)'
    )


# Match the aggregate force magnitude, for example
# ``Total force = 0.173046 Total SCF correction = 0.001`` or
# ``Total force=1.25D-04``.  Atomic ``force =`` rows, lowercase labels, and a
# value with no assignment operator are rejected.
total_force_pattern = (
    rf'\bTotal\s+force\s*=\s*(?P<total_force>{number_pattern})(?=\s|$)'
    )


# Match a complete six-column stress row: three Ry/bohr**3 entries followed by
# three kbar entries.  Examples include ``-.001 0.0 .001 -147.1 0.0 147.1``
# and exponent-form rows.  Five-column rows, labeled rows, and concatenated
# columns are rejected because they cannot form the two complete tensors.
stress_row_pattern = (
    rf'^\s*(?P<sxx>{number_pattern})\s+(?P<sxy>{number_pattern})\s+'
    rf'(?P<sxz>{number_pattern})\s+(?P<kxx>{number_pattern})\s+'
    rf'(?P<kxy>{number_pattern})\s+(?P<kxz>{number_pattern})(?=\s|$)'
    )


# Match one component of a PWSCF timing value.  This covers attached forms
# such as ``4m33.69s`` and spaced forms such as ``1h 2m 3.5s``.  Units other
# than h/m/s, unitless values, and malformed decimals do not match.  The final
# lookahead permits the next attached component while rejecting ``ms`` units.
timing_value_pattern = (
    rf'(?P<value>{number_pattern})\s*(?P<unit>[hms])(?=\s|[-+.\d]|$)'
    )


# Match the Fermi-energy result without collecting unrelated numbers earlier
# on the line.  Real singular forms include ``the Fermi energy is 10.1198 ev``
# and ``the Fermi energy = -3.22772442 eV``; spin-polarized output may report
# ``the spin up/dw Fermi energies are 5.1 5.2 ev``.  Missing eV units, three
# energies, ``highest occupied level``, and prose merely mentioning Fermi
# energy fail.
fermi_energies_pattern = (
    rf'(?i:\b(?:the\s+)?(?:spin\s+up/dw\s+)?Fermi\s+energ(?:y|ies)\s*'
    rf'(?:is|are|=)\s*)'
    rf'(?P<values>{number_pattern}(?:\s+{number_pattern})?)\s+(?i:eV)\b'
    )


def match_float(pattern,line):
    """Return the first captured floating-point value matching a pattern."""
    match = re.search(pattern,line)
    if match is None:
        return None
    return float(match.group(1).replace('D','E').replace('d','e'))
#end def match_float





class PwscfOutData(DevBase):
    """Read and store physical data from PWSCF text output.

    The calculation type and available result fields are determined solely
    from the output log.  Each reader tolerates absent or incomplete records
    and retains any complete physical data that precedes them.

    Parameters
    ----------
    filepath : str or os.PathLike
        Path to the PWSCF text-output file.
    md_only : bool, default=False
        For molecular-dynamics calculations, stop after reading the dynamics
        history and leave the remaining applicable attributes as ``None``.

    Attributes
    ----------
    calculation : {'scf', 'nscf', 'bands', 'relax', 'vc-relax', 'md', 'vc-md'}
        Calculation type inferred from the output text.  This attribute is
        present for every calculation type.
    Ef : float or None
        Final Fermi energy in eV.  Present for every calculation type.
    fermi_energies : numpy.ndarray or None
        One-dimensional history of Fermi energies in eV.  Present for every
        calculation type.
    bands : obj or None
        Spin-resolved band records.  The ``up`` and ``down`` members map
        k-point indices to objects containing ``eigs`` and ``occs`` arrays,
        k-point coordinates, polarization, and optional band-edge data.
        Present for every calculation type.
    volume : float or None
        Final unit-cell volume in bohr cubed.  Present for every calculation
        type.
    cputime, walltime : float or None
        Total CPU and wall-clock time in hours.  Present for every calculation
        type.
    kpoints_cart, kpoints_unit : numpy.ndarray or None
        Cartesian and crystal k-point arrays with shape ``(nkpoints, 3)``.
        Present for every calculation type.
    kweights : numpy.ndarray or None
        One-dimensional k-point weight array.  Present for every calculation
        type.
    E : float or None
        Final total energy in Ry.  Present for ``scf``, ``relax``,
        ``vc-relax``, ``md``, and ``vc-md`` calculations.
    relax_energies : numpy.ndarray or None
        One-dimensional sequence of completed SCF total energies in Ry.
        Present for ``scf``, ``relax``, ``vc-relax``, ``md``, and ``vc-md``.
    scf_conv_energy, scf_conv_accuracy : numpy.ndarray or None
        SCF-iteration energies and estimated accuracies in Ry.  Present for
        ``scf``, ``relax``, ``vc-relax``, ``md``, and ``vc-md``.
    pressure : float or None
        Final pressure in kbar.  Present for ``scf``, ``relax``,
        ``vc-relax``, ``md``, and ``vc-md`` calculations.
    stress : list of list of float or None
        Reported stress rows, with three rows per stress tensor.  Present for
        ``scf``, ``relax``, ``vc-relax``, ``md``, and ``vc-md``.
    forces : numpy.ndarray or None
        Atomic-force history with shape ``(nsteps, natoms, 3)`` in Ry/bohr.
        Present for ``scf``, ``relax``, ``vc-relax``, ``md``, and ``vc-md``.
    tot_forces, max_forces : numpy.ndarray or None
        One-dimensional histories of reported total forces and maximum atomic
        force magnitudes.  Present for ``scf``, ``relax``, ``vc-relax``,
        ``md``, and ``vc-md``.
    relax_structures : obj or None
        Integer-indexed structure records containing atom labels, Cartesian
        positions, and, when reported, cell axes.  Present for ``relax``,
        ``vc-relax``, ``md``, and ``vc-md`` calculations.
    md_data : obj or None
        Molecular-dynamics histories.  Its ``total_energy``, ``pressure``,
        ``time``, ``kinetic_energy``, ``temperature``, and
        ``potential_energy`` members are one-dimensional NumPy arrays.
        Present for ``md`` and ``vc-md`` calculations.
    md_stats : obj or None
        Molecular-dynamics summary statistics.  Each member is a
        ``(mean, error)`` tuple of floats.  Present for ``md`` and ``vc-md``.

    Notes
    -----
    The expected attribute groups by calculation type are:

    ``scf``
        Common attributes plus energy, convergence, pressure, stress, and
        force attributes.
    ``nscf`` and ``bands``
        Common attributes only.
    ``relax`` and ``vc-relax``
        The ``scf`` attributes plus ``relax_structures``.
    ``md`` and ``vc-md``
        The relaxation attributes plus ``md_data`` and ``md_stats``.

    An applicable attribute remains ``None`` when its record is absent or
    cannot be parsed.  Attributes that do not apply to the inferred
    calculation type are removed from the object.
    """

    def __init__(self,filepath,*,md_only=False):
        """Read a PWSCF log and initialize its accessible physical data."""
        self.calculation = None

        # all calculation types
        self.Ef                = None
        self.fermi_energies    = None
        self.bands             = None
        self.volume            = None
        self.cputime           = None
        self.walltime          = None
        self.kpoints_cart      = None
        self.kpoints_unit      = None
        self.kweights          = None
        # scf/relax/vc-relax/md/vc-md
        self.E                 = None
        self.relax_energies    = None
        self.scf_conv_energy   = None
        self.scf_conv_accuracy = None
        self.pressure          = None
        self.stress            = None
        self.forces            = None
        self.tot_forces        = None
        self.max_forces        = None
        # md/vc-md
        self.md_data           = None
        self.md_stats          = None
        # md/vc-md/relax/vc-relax
        self.relax_structures  = None

        with open(filepath,'r') as fobj:
            lines = fobj.read().splitlines()
        # read the calculation type
        self.read_calculation(lines)
        # remove unused attributes, depending on the calculation type
        if self.calculation in ('nscf','bands'):
            for name in (
                'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
                'pressure','stress','forces','tot_forces','max_forces',
                ):
                del self[name]
            #end for
        if self.calculation in ('scf','nscf','bands','relax','vc-relax'):
            del self.md_data
            del self.md_stats
        if self.calculation in ('scf','nscf','bands'):
            del self.relax_structures
        # read md
        if self.calculation in ('md','vc-md'):
            self.read_md(lines)
            if md_only:
                return
        # all calculations
        self.read_fermi_energies(lines)
        self.read_kpoints(lines)
        self.read_bands(lines)
        self.read_pressure_volume(lines)
        # all but nscf/bands
        if self.calculation in ('scf','md','vc-md','relax','vc-relax'):
            self.read_scf_convergence(lines)
            self.read_energies(lines)
            self.read_stress(lines)
            self.read_forces(lines)
        # all but scf/nscf/bands
        if self.calculation in ('md','vc-md','relax','vc-relax'):
            self.read_structures(lines)

        self.read_timing(lines)
    #end def __init__


    def read_calculation(self,lines):
        """Infer and bind the PWSCF calculation type from log records."""
        has_cell        = any(line.strip().startswith('CELL_PARAMETERS') for line in lines)
        has_bfgs        = any('BFGS Geometry Optimization' in line for line in lines)
        has_band_run    = any('Band Structure Calculation' in line for line in lines)
        has_vc_dynamics = any('Entering Dynamics;' in line for line in lines)
        has_dynamics    = has_vc_dynamics or any(
            'Molecular Dynamics Calculation' in line or 'Entering Dynamics:' in line
            for line in lines
            )
        if has_dynamics:
            calculation = 'vc-md' if has_vc_dynamics or has_cell else 'md'
        elif has_bfgs:
            calculation = 'vc-relax' if has_cell else 'relax'
        elif has_band_run:
            # QE uses the same heading for nscf and bands, but suppresses
            # electronic-reference and occupation records for bands runs.
            has_reference = any(
                'Fermi energ' in line
                or 'highest occupied' in line
                or 'occupation numbers' in line
                for line in lines
                )
            calculation = 'nscf' if has_reference else 'bands'
        else:
            calculation = 'scf'
        self.calculation = calculation
    #end def read_calculation


    def read_md(self,lines):
        """Read and bind molecular-dynamics histories.

        Complete ionic steps are stored in ``md_data`` as an ``obj`` of
        one-dimensional NumPy arrays named ``total_energy``, ``pressure``,
        ``time``, ``kinetic_energy``, ``temperature``, and
        ``potential_energy``.  ``md_stats`` contains a ``(mean, error)``
        tuple for each array.  Fixed-cell ``md`` and variable-cell ``vc-md``
        outputs provide time, kinetic energy, and temperature in different
        line formats, but produce the same stored data structure.

        Only steps containing all required quantities are retained.
        """
        def record_value(name,pattern,line):
            """Add a matched value to the current dynamics record."""
            value = match_float(pattern,line)
            if value is not None and 'total_energy' in record:
                record[name] = value
        #end def record_value

        if self.calculation not in ('md','vc-md'):
            return

        records  = []
        record   = dotdict()
        required = ('total_energy','pressure','time','kinetic_energy','temperature')
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                energy = match_float(total_energy_pattern,line)
                record = dotdict()
                if energy is not None:
                    record.total_energy = energy
            elif 'total   stress' in line and 'P=' in line:
                record_value('pressure',pressure_pattern,line)
            if self.calculation=='md':
                if 'time      =' in line:
                    record_value('time',md_time_pattern,line)
                elif 'kinetic energy' in line and '=' in line:
                    record_value('kinetic_energy',kinetic_energy_pattern,line)
                elif line.strip().startswith('temperature') and '=' in line:
                    record_value('temperature',temperature_pattern,line)
            else:
                if 'Entering Dynamics;' in line and 'time' in line:
                    record_value('time',md_time_pattern,line)
                elif line.strip().startswith('Ekin'):
                    match = re.search(ekin_temperature_pattern,line)
                    if match is not None and 'total_energy' in record:
                        record.kinetic_energy = float(
                            match.group('kinetic_energy').replace('D','E').replace('d','e')
                            )
                        record.temperature    = float(
                            match.group('temperature').replace('D','E').replace('d','e')
                            )
            if all(name in record for name in required):
                records.append(record)
                record = dotdict()
        if len(records)==0:
            return
        md = obj({
            name:np.array([record[name] for record in records],dtype=float)
            for name in required
            })
        md.potential_energy = md.total_energy - md.kinetic_energy
        self.md_data        = md
        self.md_stats       = self.md_statistics()
    #end def read_md


    def read_fermi_energies(self,lines):
        """Read and bind the sequence of reported Fermi energies.

        ``fermi_energies`` is a one-dimensional NumPy array containing every
        successfully parsed value in eV, and ``Ef`` is its final value.  Both
        remain ``None`` when no Fermi-energy record is available.
        """
        fermi_energies = []
        for line in lines:
            if 'Fermi energ' in line:
                match = re.search(fermi_energies_pattern,line)
                if match is not None:
                    values = re.findall(number_pattern,match.group('values'))
                    fermi_energies.extend(
                        float(value.replace('D','E').replace('d','e'))
                        for value in values
                        )
        if len(fermi_energies)>0:
            self.Ef             = fermi_energies[-1]
            self.fermi_energies = np.array(fermi_energies,dtype=float)
    #end def read_fermi_energies


    def read_energies(self,lines):
        """Read and bind completed SCF total energies.

        ``relax_energies`` is a one-dimensional NumPy array of the total
        energies marked with ``!`` in the output, in Ry, and ``E`` is its
        final value.  For relaxation and dynamics calculations the array is
        the energy history over ionic steps; for ``scf`` it normally contains
        one entry.
        """
        relax_energies = []
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                energy = match_float(total_energy_pattern,line)
                if energy is not None:
                    relax_energies.append(energy)
        if len(relax_energies)>0:
            self.E              = relax_energies[-1]
            self.relax_energies = np.array(relax_energies,dtype=float)
    #end def read_energies


    def read_scf_convergence(self,lines):
        """Read and bind electronic SCF convergence histories.

        ``scf_conv_energy`` and ``scf_conv_accuracy`` are one-dimensional
        NumPy arrays in Ry containing iterative, non-final total energies and
        their following estimated accuracies.  Relaxation and dynamics runs
        can contribute multiple SCF cycles to the same flattened histories.
        Missing accuracy lines are skipped, so the two arrays can have
        different lengths.
        """
        scf_conv_energy   = []
        scf_conv_accuracy = []
        capture_accuracy  = False
        for line in lines:
            if 'total energy' in line and '=' in line:
                capture_accuracy = False
                if not line.lstrip().startswith('!'):
                    energy = match_float(total_energy_pattern,line)
                    if energy is not None:
                        scf_conv_energy.append(energy)
                        capture_accuracy = True
            elif capture_accuracy and 'estimated scf accuracy' in line:
                accuracy = match_float(scf_accuracy_pattern,line)
                if accuracy is not None:
                    scf_conv_accuracy.append(accuracy)
                capture_accuracy = False
        if len(scf_conv_energy)>0:
            self.scf_conv_energy = np.array(scf_conv_energy,dtype=float)
        if len(scf_conv_accuracy)>0:
            self.scf_conv_accuracy = np.array(scf_conv_accuracy,dtype=float)
    #end def read_scf_convergence


    def read_bands(self,lines):
        """Read and bind band data for each reported k-point.

        ``bands`` is an ``obj`` with ``up`` and ``down`` members.  Each member
        maps a zero-based k-point index to an ``obj`` containing ``index``,
        ``kpoint_2pi_alat``, ``kpoint_rel``, ``eigs``, ``occs``, and ``pol``.
        Eigenvalues and occupations are one-dimensional NumPy arrays in eV
        and electrons, respectively; either k-point representation can be
        ``None`` when its source table was not printed.

        Non-spin-polarized output places all records in ``bands.up`` and sets
        ``pol`` to ``'none'``.  Spin-polarized output separates records into
        ``up`` and ``down`` and labels them accordingly.  For ``bands`` and
        some ``nscf`` outputs, occupation arrays can be empty.  When complete
        occupations are present, :meth:`read_band_edges` adds band-edge and
        gap information to the same object.
        """
        def leading_numbers(line):
            """Return the leading sequence of numeric values from a line."""
            match = re.match(leading_number_list_pattern,line)
            if match is None:
                return np.array([],dtype=float)
            return np.array(
                re.findall(number_pattern,match.group('values').replace('D','E').replace('d','e')),
                dtype=float,
                )
        #end def leading_numbers

        def read_values(start,markers=()):
            """Read a contiguous block of numeric values."""
            values = []
            i      = start
            while i<len(lines):
                text = lines[i].strip()
                if len(text)==0:
                    if len(values)>0:
                        break
                elif any(marker in text for marker in markers):
                    break
                else:
                    numbers = leading_numbers(text)
                    if len(numbers)==0:
                        break
                    values.extend(numbers)
                i+=1
            return values,i
        #end def read_values

        table_cart  = self.kpoints_cart
        table_unit  = self.kpoints_unit
        polarized    = any('- SPIN UP -' in line or '- SPIN DOWN -' in line for line in lines)
        bands        = obj(up=obj(),down=obj())
        band_channel = bands.up
        up_spin      = True
        for i,line in enumerate(lines):
            if ('End of self-consistent calculation' in line
                and len(bands.up)+len(bands.down)>0
                ):
                bands        = obj(up=obj(),down=obj())
                band_channel = bands.up
                up_spin      = True
                continue
            elif '- SPIN UP -' in line:
                up_spin      = True
                band_channel = bands.up
                continue
            elif '- SPIN DOWN -' in line:
                up_spin      = False
                band_channel = bands.down
                continue
            elif 'bands (ev)' not in line:
                continue

            eigs,j = read_values(i+1,('occupation numbers','bands (ev)'))
            if len(eigs)==0:
                continue
            while j<len(lines) and len(lines[j].strip())==0:
                j+=1
            occs = []
            if j<len(lines) and 'occupation numbers' in lines[j]:
                occs,_ = read_values(j+1)

            match       = re.search(band_kpoint_pattern,line)
            kpoint_cart = None
            if match is not None:
                values = [
                    float(match.group(name).replace('D','E').replace('d','e'))
                    for name in ('kx','ky','kz')
                    ]
                kpoint_cart = np.array(values,dtype=float)
            index      = len(band_channel)
            kpoint_rel = kpoint_cart
            if table_unit is not None and index<len(table_unit):
                kpoint_rel = table_unit[index]
            if table_cart is not None and index<len(table_cart):
                kpoint_cart = table_cart[index]
            pol = ('up' if up_spin else 'down') if polarized else 'none'
            band_channel[index] = obj(
                index           = index,
                kpoint_2pi_alat = kpoint_cart,
                kpoint_rel      = kpoint_rel,
                eigs            = np.array(eigs,dtype=float),
                occs            = np.array(occs,dtype=float),
                pol             = pol,
                )
        if len(bands.up)+len(bands.down)==0:
            return
        self.bands = bands
        self.read_band_edges()
    #end def read_bands


    def read_band_edges(self):
        """Add band-edge and gap records to a parsed bands object.

        Occupied and unoccupied eigenvalues are identified independently at
        every k-point and spin.  When usable occupations exist, ``bands`` is
        augmented with ``vbm``, ``cbm``, ``direct_gap``, and
        ``electronic_structure``.  Edge objects contain the energy, k-point
        representations, k-point index, polarization, and, for VBM/CBM,
        band number.  Insulating systems whose extrema occur at different
        k-points also receive ``indirect_gap``, whose ``kpoints`` member holds
        the VBM and CBM records.  No members are added when occupations are
        absent or cannot distinguish occupied and unoccupied states.
        """
        def edge_data(band,energy,band_number):
            """Build a band-edge record for one k-point."""
            return obj(
                energy          = energy,
                kpoint_rel      = band.kpoint_rel,
                kpoint_2pi_alat = band.kpoint_2pi_alat,
                index           = band.index,
                pol             = band.pol,
                band_number     = band_number,
                )
        #end def edge_data

        bands      = self.bands
        vbm        = None
        cbm        = None
        direct_gap = None
        for band_channel in (bands.up,bands.down):
            for band in band_channel.values():
                if len(band.occs)!=len(band.eigs) or len(band.occs)==0:
                    continue
                occ   = band.occs > 0.5
                unocc = band.occs < 0.5
                if not occ.any() or not unocc.any():
                    continue
                e_val  = np.max(band.eigs[occ])
                e_cond = np.min(band.eigs[unocc])
                if vbm is None or e_val>vbm.energy:
                    vbm = edge_data(band,e_val,np.max(np.where(occ)))
                if cbm is None or e_cond<cbm.energy:
                    cbm = edge_data(band,e_cond,np.min(np.where(unocc)))
                if direct_gap is None or e_cond-e_val<direct_gap.energy:
                    direct_gap = obj(
                        energy          = e_cond-e_val,
                        kpoint_rel      = band.kpoint_rel,
                        kpoint_2pi_alat = band.kpoint_2pi_alat,
                        index           = band.index,
                        pol             = [vbm.pol,cbm.pol],
                        )
        if vbm is None:
            return
        if vbm.energy+0.025>=cbm.energy:
            electronic_structure = (
                'metallic' if vbm.band_number==cbm.band_number else 'semi-metal'
                )
        else:
            electronic_structure = 'insulating'
            if (vbm.kpoint_rel is not None
                and cbm.kpoint_rel is not None
                and not np.equal(vbm.kpoint_rel,cbm.kpoint_rel).all()
                ):
                bands.indirect_gap = obj(
                    energy  = round(cbm.energy-vbm.energy,3),
                    kpoints = obj(vbm=vbm,cbm=cbm),
                    )
        bands.update(
            electronic_structure = electronic_structure,
            vbm                  = vbm,
            cbm                  = cbm,
            direct_gap           = direct_gap,
            )
    #end def read_band_edges


    def read_structures(self,lines):
        """Read and bind structures from ionic-step output blocks.

        ``relax_structures`` is an ``obj`` mapping zero-based step indices to
        configuration objects.  Each configuration contains ``atoms`` as a
        list of element labels and ``positions`` as an ``(natoms, 3)`` NumPy
        array.  An ``axes`` ``(3, 3)`` array is included when a preceding
        ``CELL_PARAMETERS`` block is available.  Crystal positions are
        converted to Cartesian coordinates when those axes are known.

        The layout is the same for ``relax``, ``vc-relax``, ``md``, and
        ``vc-md``.  Fixed-cell output can omit cell blocks, while variable-cell
        output normally supplies new axes with each structure.
        """
        structures = obj()
        conf       = None
        i          = 0
        while i<len(lines):
            line = lines[i]
            if line.strip().startswith('CELL_PARAMETERS'):
                axes = []
                if i+3<len(lines):
                    for axis_line in lines[i+1:i+4]:
                        match = re.match(vector3_pattern,axis_line)
                        if match is None:
                            axes = []
                            break
                        values = [
                            float(match.group(name).replace('D','E').replace('d','e'))
                            for name in ('x','y','z')
                            ]
                        axes.append(values)
                if len(axes)==3:
                    conf = obj()
                    axes = np.array(axes,dtype=float)
                    alat = match_float(alat_pattern,line)
                    if alat is not None:
                        axes *= alat
                    conf.axes = axes
                    i+=3
                else:
                    conf = None
            elif 'ATOMIC_POSITIONS' in line:
                if conf is None:
                    conf = obj()
                atoms     = []
                positions = []
                i+=1
                while i<len(lines):
                    tokens = lines[i].split()
                    if len(tokens)<4 or tokens[0].lower()=='end':
                        break
                    coordinates = tokens[1:4]
                    # Check whether every coordinate is a complete numeric value.
                    if not all(
                        re.fullmatch(number_pattern,value.strip()) is not None
                        for value in coordinates
                        ):
                        break
                    atoms.append(tokens[0])
                    positions.append([
                        float(value.replace('D','E').replace('d','e'))
                        for value in coordinates
                        ])
                    i+=1
                if len(positions)==0:
                    conf = None
                    continue
                conf.atoms     = atoms
                conf.positions = np.array(positions,dtype=float)
                if 'crystal' in line.lower() and 'axes' in conf:
                    conf.positions = np.dot(conf.positions,conf.axes)
                structures[len(structures)] = conf
                conf = None
                continue
            i+=1
        if len(structures)>0:
            self.relax_structures = structures
    #end def read_structures


    def read_pressure_volume(self,lines):
        """Read the final reported pressure and unit-cell volume."""
        pressure = None
        volume   = None
        for line in lines:
            if 'unit-cell volume' in line:
                value = match_float(volume_pattern,line)
                if value is not None:
                    volume = value
            if 'total' in line and 'stress' in line and 'P=' in line:
                value = match_float(pressure_pattern,line)
                if value is not None:
                    pressure = value
        if pressure is not None and 'pressure' in self:
            self.pressure = pressure
        if volume is not None:
            self.volume = volume
    #end def read_pressure_volume


    def read_stress(self,lines):
        """Read and bind the sequence of reported stress tensors.

        ``stress`` is a list of numeric rows, with three consecutive rows per
        complete tensor.  Each row contains the three stress components in
        Ry/bohr cubed followed by the three values printed in kbar.  Tensors
        from successive ionic steps are appended to the same flat list.
        """
        stress = []
        for i,line in enumerate(lines):
            if 'total   stress' in line:
                rows = []
                if i+3<len(lines):
                    for stress_line in lines[i+1:i+4]:
                        match = re.match(stress_row_pattern,stress_line)
                        if match is None:
                            rows = []
                            break
                        names  = ('sxx','sxy','sxz','kxx','kxy','kxz')
                        values = [
                            float(match.group(name).replace('D','E').replace('d','e'))
                            for name in names
                            ]
                        rows.append(values)
                if len(rows)==3:
                    stress.extend(rows)
        if len(stress)>0:
            self.stress = stress
    #end def read_stress


    def read_forces(self,lines):
        """Read and bind atomic and total force histories.

        ``forces`` is a NumPy array with shape ``(nsteps, natoms, 3)`` in
        Ry/bohr.  ``max_forces`` is the maximum atomic-force norm for each
        retained step, and ``tot_forces`` is a one-dimensional array of every
        separately reported total-force value.  Atomic-force blocks with a
        known atom count are retained only when all atoms are present; total
        forces remain available even when their atomic block is incomplete.
        """
        forces     = []
        tot_forces = []
        nat        = None
        for line in lines:
            if 'number of atoms/cell' in line:
                match = re.search(r'number of atoms/cell\s*=\s*(\d+)',line)
                if match is not None:
                    nat = int(match.group(1))
                    break
        for i,line in enumerate(lines):
            if 'Forces acting on atoms' in line:
                aforces = []
                j       = i+1
                while j<len(lines):
                    match = re.search(atomic_force_pattern,lines[j])
                    if match is not None:
                        values = [
                            float(match.group(name).replace('D','E').replace('d','e'))
                            for name in ('fx','fy','fz')
                            ]
                        aforces.append(values)
                    elif len(aforces)>0:
                        break
                    j+=1
                if len(aforces)>0 and (nat is None or len(aforces)==nat):
                    forces.append(aforces)
            if 'Total force' in line:
                match = re.search(total_force_pattern,line)
                if match is not None:
                    value = match.group('total_force').replace('D','E').replace('d','e')
                    tot_forces.append(float(value))
        if len(forces)>0:
            forces          = np.array(forces,dtype=float)
            self.forces     = forces
            self.max_forces = np.linalg.norm(forces,axis=2).max(axis=1)
        if len(tot_forces)>0:
            self.tot_forces = np.array(tot_forces,dtype=float)
    #end def read_forces


    def read_timing(self,lines):
        """Read total CPU and wall-clock times."""
        def pwscf_time(text):
            """Convert a PWSCF timing string to hours."""
            scales = dict(h=1.,m=60.,s=3600.)
            return sum(
                float(match.group('value').replace('D','E').replace('d','e'))
                / scales[match.group('unit')]
                for match in re.finditer(timing_value_pattern,text)
                )
        #end def pwscf_time

        for line in lines:
            if 'PWSCF        :' in line:
                match = re.search(r'PWSCF\s*:\s*(.*?)\s+CPU\s+(.*?)\s+WALL',line)
                if match is not None:
                    self.cputime  = pwscf_time(match.group(1))
                    self.walltime = pwscf_time(match.group(2))
                    return
    #end def read_timing


    def read_kpoints(self,lines):
        """Read and bind paired k-point tables and their weights.

        A complete Cartesian table and its following crystal-coordinate table
        are required.  ``kpoints_cart`` and ``kpoints_unit`` are NumPy arrays
        with shape ``(nkpoints, 3)`` in units of ``2 pi/alat`` and crystal
        reciprocal coordinates, respectively.  ``kweights`` is the matching
        one-dimensional weight array.  No member is updated when either table
        is incomplete.
        """
        for i,line in enumerate(lines):
            if 'number of k points=' not in line:
                continue
            match = re.search(r'number of k points=\s*(\d+)',line)
            if match is None:
                continue
            nkpoints = int(match.group(1))
            if i+1>=len(lines) or 'cart. coord.' not in lines[i+1]:
                continue
            cart    = []
            weights = []
            valid   = len(lines[i+2:i+2+nkpoints])==nkpoints
            for kline in lines[i+2:i+2+nkpoints]:
                match = re.search(kpoint_table_pattern,kline)
                if match is None:
                    valid = False
                    break
                coordinates = [
                    float(match.group(name).replace('D','E').replace('d','e'))
                    for name in ('kx','ky','kz')
                    ]
                cart.append(coordinates)
                weights.append(
                    float(match.group('weight').replace('D','E').replace('d','e'))
                    )
            if not valid:
                continue
            j = i+2+nkpoints
            while j<len(lines) and 'cryst. coord.' not in lines[j]:
                j+=1
            if j>=len(lines):
                continue
            unit  = []
            valid = len(lines[j+1:j+1+nkpoints])==nkpoints
            for kline in lines[j+1:j+1+nkpoints]:
                match = re.search(kpoint_table_pattern,kline)
                if match is None:
                    valid = False
                    break
                coordinates = [
                    float(match.group(name).replace('D','E').replace('d','e'))
                    for name in ('kx','ky','kz')
                    ]
                unit.append(coordinates)
            if not valid:
                continue
            self.kpoints_cart = np.array(cart,dtype=float)
            self.kpoints_unit = np.array(unit,dtype=float)
            self.kweights     = np.array(weights,dtype=float)
            return
    #end def read_kpoints


    def md_statistics(self,equil=None,autocorr=None):
        """Build summary statistics for each molecular-dynamics history.

        The returned ``obj`` has the same quantity names as ``md_data``.  Each
        value is a ``(mean, error)`` tuple.  ``equil`` discards an initial
        number of samples; ``autocorr`` optionally groups the remaining data
        into blocks of that length before estimating the error.  This method
        is used to bind ``md_stats`` after a complete MD history is read.
        """
        mds = obj()
        for quantity,values in self.md_data.items():
            if equil is not None:
                values = values[equil:]
            if autocorr is None:
                mean,_,error,_ = simstats(values)
            else:
                nvalues = len(values)
                nblocks = int(np.floor(float(nvalues)/autocorr))
                values  = values[nvalues-nblocks*autocorr:]
                npe.reshape_inplace(values,(nblocks,autocorr))
                mean,error = simplestats(values.mean(axis=1))
            mds[quantity] = mean,error
        return mds
    #end def md_statistics

#end class PwscfOutData



class PwscfXmlData(DevBase):
    """Read physical results from schema-based Quantum ESPRESSO XML output.

    The complete converted XML hierarchy is retained while commonly used
    materials data are exposed directly as scalars and NumPy arrays.  Results
    that are absent or incomplete remain ``None``.

    Parameters
    ----------
    filepath : str or os.PathLike
        Path to a ``data-file-schema.xml`` file.

    Attributes
    ----------
    data : obj or None
        Complete converted XML hierarchy rooted at ``data.root``.  It is
        ``None`` only when the XML file cannot be read.
    parse_failed : bool
        Whether direct reading or XML syntax parsing failed.
    version : str or None
        Quantum ESPRESSO version reported by the XML creator record.
    calculation : {'scf', 'nscf', 'bands', 'relax', 'vc-relax', 'md', 'vc-md'} or None
        Calculation type reported in the input section of the XML file.
    exit_status : int or None
        Program exit status.
    total_energy, band_energy, hartree_energy : float or None
        Final total, band, and Hartree energies.
    xc_energy, xc_potential_energy, ewald_energy : float or None
        Final exchange-correlation, exchange-correlation potential, and Ewald
        energy contributions.
    scf_converged : bool or None
        Whether the final self-consistent calculation converged.
    scf_steps : int or None
        Number of electronic iterations in the final SCF cycle.
    scf_error : float or None
        Estimated error from the final SCF cycle.
    initial_atoms, atoms : numpy.ndarray or None
        One-dimensional arrays of initial and final atomic species labels.
    initial_positions, positions : numpy.ndarray or None
        Initial and final Cartesian atomic positions with shape ``(natoms, 3)``.
    initial_axes, axes : numpy.ndarray or None
        Initial and final cell axes with shape ``(3, 3)``.
    initial_alat, alat : float or None
        Initial and final lattice parameters reported by PWSCF.
    initial_volume, volume : float or None
        Initial and final cell volumes computed from the cell axes.
    forces : numpy.ndarray or None
        Final atomic forces with shape ``(natoms, 3)``.  This attribute is not
        present for ``nscf`` and ``bands`` calculations.
    stress : numpy.ndarray or None
        Final stress tensor with shape ``(3, 3)``.  This attribute is not
        present for ``nscf`` and ``bands`` calculations.
    functional : str or None
        Exchange-correlation functional reported by PWSCF.
    spin_polarized, noncollinear, spin_orbit : bool or None
        Collinear spin-polarization, noncollinearity, and spin-orbit flags.
    absolute_magnetization : float or None
        Final absolute magnetization.
    nelectrons : float or None
        Number of electrons.
    nbands, nkpoints : int or None
        Numbers of electronic bands and k-points.
    kgrid, kshift : numpy.ndarray or None
        Three-element integer arrays containing the Monkhorst-Pack grid and
        grid shift.
    fermi_energy : float or None
        Final Fermi energy.
    highest_occupied_level, lowest_unoccupied_level : float or None
        Highest occupied and lowest unoccupied electronic energy levels.
    occupations_kind : str or None
        Occupation method reported in the band-structure record.
    kpoints_rel : numpy.ndarray or None
        K-point coordinates in units of ``2 pi/alat`` with shape
        ``(nkpoints, 3)``.
    kweights : numpy.ndarray or None
        K-point weights with shape ``(nkpoints,)``.
    eigenvalues, occupations : numpy.ndarray or None
        Electronic eigenvalues and occupations.  The shape is
        ``(nkpoints, nbands)`` for non-spin-polarized calculations and
        ``(2, nkpoints, nbands)`` for collinear spin-polarized calculations.
    plane_waves : numpy.ndarray or None
        Number of plane waves at each k-point.  The shape is ``(nkpoints,)``
        for non-spin-polarized calculations and ``(2, nkpoints)`` for
        collinear spin-polarized calculations.
    species : numpy.ndarray or None
        One-dimensional array of unique atomic species labels.
    atomic_masses : numpy.ndarray or None
        Atomic mass for each species.
    pseudopotentials : numpy.ndarray or None
        Pseudopotential filename for each species.
    valence_charges : numpy.ndarray or None
        Pseudopotential valence charge for each species.
    ecutwfc, ecutrho : float or None
        Wavefunction and charge-density plane-wave cutoffs.
    gamma_only : bool or None
        Whether the plane-wave basis uses gamma-only storage.
    ngm, ngms, npwx : int or None
        Plane-wave and smooth-grid vector counts and the maximum number of
        wavefunction plane waves.
    fft_grid, fft_smooth, fft_box : numpy.ndarray or None
        Three-element integer arrays containing FFT-grid dimensions.
    reciprocal_lattice : numpy.ndarray or None
        Reciprocal-lattice vectors with shape ``(3, 3)``.
    nsym, nrot, space_group : int or None
        Numbers of active symmetries and rotations and the space-group index.
    symmetry_names, symmetry_classes : numpy.ndarray or None
        One-dimensional string arrays describing each parsed symmetry.
    symmetry_rotations : numpy.ndarray or None
        Symmetry rotation matrices with shape ``(nsym, 3, 3)``.
    symmetry_translations : numpy.ndarray or None
        Fractional translations with shape ``(nsym, 3)``.
    equivalent_atoms : numpy.ndarray or None
        Equivalent-atom indices with shape ``(nsym, natoms)``.
    cpu_time, wall_time : float or None
        Total CPU and wall-clock times reported in the XML timing record.
    optimization_converged : bool or None
        Whether ionic optimization converged.  Present only for ``relax`` and
        ``vc-relax`` calculations.
    optimization_steps : int or None
        Number of ionic optimization steps.  Present only for ``relax`` and
        ``vc-relax`` calculations.
    gradient_norm : float or None
        Final ionic gradient norm.  Present only for ``relax`` and
        ``vc-relax`` calculations.
    trajectory_energies, trajectory_band_energies : numpy.ndarray or None
        Total- and band-energy histories over ionic steps.
    trajectory_hartree_energies, trajectory_xc_energies : numpy.ndarray or None
        Hartree- and exchange-correlation-energy histories over ionic steps.
    trajectory_xc_potential_energies, trajectory_ewald_energies : numpy.ndarray or None
        Exchange-correlation-potential and Ewald-energy histories over ionic
        steps.
    trajectory_positions : numpy.ndarray or None
        Atomic-position history with shape ``(nsteps, natoms, 3)``.
    trajectory_axes : numpy.ndarray or None
        Cell-axis history with shape ``(nsteps, 3, 3)``.
    trajectory_volumes : numpy.ndarray or None
        Cell-volume history with shape ``(nsteps,)``.
    trajectory_forces : numpy.ndarray or None
        Atomic-force history with shape ``(nsteps, natoms, 3)``.
    trajectory_stress : numpy.ndarray or None
        Stress-tensor history with shape ``(nsteps, 3, 3)``.
    trajectory_scf_converged : numpy.ndarray or None
        Boolean SCF-convergence history over ionic steps.
    trajectory_scf_steps : numpy.ndarray or None
        Number of SCF iterations at each ionic step.
    trajectory_scf_error : numpy.ndarray or None
        Estimated SCF error at each ionic step.

    Notes
    -----
    Trajectory attributes are present only for ``relax``, ``vc-relax``,
    ``md``, and ``vc-md`` calculations.  Applicable attributes remain
    ``None`` when their XML records are absent or incomplete; attributes that
    do not apply to the reported calculation type are removed.  Physical
    quantities retain the native units of the Quantum ESPRESSO schema,
    generally atomic units.
    """

    def __init__(self,filepath):
        """Read schema XML and initialize extracted PWSCF results."""
        self.data                         = None
        self.parse_failed                 = False
        self.version                      = None
        self.calculation                  = None
        self.exit_status                  = None
        self.total_energy                 = None
        self.band_energy                  = None
        self.hartree_energy               = None
        self.xc_energy                    = None
        self.xc_potential_energy          = None
        self.ewald_energy                 = None
        self.scf_converged                = None
        self.scf_steps                    = None
        self.scf_error                    = None
        self.initial_atoms                = None
        self.initial_positions            = None
        self.initial_axes                 = None
        self.initial_alat                 = None
        self.initial_volume               = None
        self.atoms                        = None
        self.positions                    = None
        self.axes                         = None
        self.alat                         = None
        self.volume                       = None
        self.forces                       = None
        self.stress                       = None
        self.functional                   = None
        self.spin_polarized               = None
        self.noncollinear                 = None
        self.spin_orbit                   = None
        self.absolute_magnetization       = None
        self.nelectrons                   = None
        self.nbands                       = None
        self.nkpoints                     = None
        self.kgrid                        = None
        self.kshift                       = None
        self.fermi_energy                 = None
        self.highest_occupied_level       = None
        self.lowest_unoccupied_level      = None
        self.occupations_kind             = None
        self.kpoints_rel                  = None
        self.kweights                     = None
        self.eigenvalues                  = None
        self.occupations                  = None
        self.plane_waves                  = None
        self.species                      = None
        self.atomic_masses                = None
        self.pseudopotentials             = None
        self.valence_charges              = None
        self.ecutwfc                      = None
        self.ecutrho                      = None
        self.gamma_only                   = None
        self.ngm                          = None
        self.ngms                         = None
        self.npwx                         = None
        self.fft_grid                     = None
        self.fft_smooth                   = None
        self.fft_box                      = None
        self.reciprocal_lattice           = None
        self.nsym                         = None
        self.nrot                         = None
        self.space_group                  = None
        self.symmetry_names               = None
        self.symmetry_classes             = None
        self.symmetry_rotations           = None
        self.symmetry_translations        = None
        self.equivalent_atoms             = None
        self.cpu_time                     = None
        self.wall_time                    = None
        self.optimization_converged       = None
        self.optimization_steps           = None
        self.gradient_norm                = None
        self.trajectory_energies              = None
        self.trajectory_band_energies         = None
        self.trajectory_hartree_energies      = None
        self.trajectory_xc_energies           = None
        self.trajectory_xc_potential_energies = None
        self.trajectory_ewald_energies        = None
        self.trajectory_positions             = None
        self.trajectory_axes                  = None
        self.trajectory_volumes               = None
        self.trajectory_forces                = None
        self.trajectory_stress                = None
        self.trajectory_scf_converged         = None
        self.trajectory_scf_steps             = None
        self.trajectory_scf_error             = None

        def xml_value(text):
            """Convert XML text to its natural scalar or array value."""
            text = text.strip()
            if len(text)==0:
                return ''
            elif text.lower() in ('true','false'):
                return text.lower()=='true'
            match = re.fullmatch(numeric_text_pattern,text)
            if match is not None:
                tokens = re.findall(number_pattern,match.group('values'))
                values = np.array([
                    float(token.replace('D','E').replace('d','e'))
                    for token in tokens
                    ],dtype=float)
                if len(values)==1:
                    value = values[0]
                    if value.is_integer() and all(c not in text.lower() for c in ('.','e')):
                        return int(value)
                    return value
                return values
            return text
        #end def xml_value

        def xml_element(element):
            """Convert an XML element tree to an object hierarchy."""
            node = obj()
            for name,value in element.attrib.items():
                node[name.lower()] = xml_value(value)
            groups = {}
            for child in element:
                name = child.tag.rsplit('}',1)[-1].lower()
                groups.setdefault(name,[]).append(child)
            for name,group in groups.items():
                node[name] = (
                    xml_element(group[0])
                    if len(group)==1
                    else obj({i+1:xml_element(child) for i,child in enumerate(group)})
                    )
            text = (element.text or '').strip()
            if len(text)>0:
                value = xml_value(text)
                if len(node)==0:
                    return value
                node.value = value
            return node
        #end def xml_element

        try:
            root = ET.parse(filepath).getroot()
        except (OSError,UnicodeError,ET.ParseError):
            self.parse_failed = True
            return

        self.data = obj(root=xml_element(root))
        self.extract_results(root)
    #end def __init__


    def extract_results(self,root):
        """Extract user-facing physical quantities from the schema XML tree."""
        def child(element,name):
            """Return a named direct child element."""
            if element is None:
                return None
            name = name.lower()
            return next(
                (item for item in element
                 if item.tag.rsplit('}',1)[-1].lower()==name),
                None,
                )
        #end def child

        def children(element,name):
            """Return all named direct child elements."""
            if element is None:
                return []
            name = name.lower()
            return [
                item for item in element
                if item.tag.rsplit('}',1)[-1].lower()==name
                ]
        #end def children

        def path(element,*names):
            """Follow a direct-child path through the XML tree."""
            for name in names:
                element = child(element,name)
                if element is None:
                    break
            return element
        #end def path

        def value(element):
            """Convert a scalar or numeric-array element value."""
            if element is None:
                return None
            text = (element.text or '').strip()
            if len(text)==0:
                return None
            if text.lower() in ('true','false'):
                return text.lower()=='true'
            match = re.fullmatch(numeric_text_pattern,text)
            if match is None:
                return text
            tokens = re.findall(number_pattern,match.group('values'))
            values = np.array([
                float(token.replace('D','E').replace('d','e'))
                for token in tokens
                ],dtype=float)
            if len(values)>1:
                return values
            number = values[0]
            if number.is_integer() and all(c not in text.lower() for c in ('.','e')):
                return int(number)
            return float(number)
        #end def value

        def array(element,dtype=float):
            """Return element text as a one-dimensional array."""
            result = value(element)
            if result is None or isinstance(result,(str,bool)):
                return None
            return np.asarray(result,dtype=dtype).reshape(-1)
        #end def array

        def matrix(element,ncols):
            """Return element text as a matrix with a fixed column count."""
            result = array(element)
            if result is None or len(result)%ncols!=0:
                return None
            return result.reshape(-1,ncols)
        #end def matrix

        def structure(element):
            """Return atom labels, positions, and cell axes from a structure."""
            atoms     = None
            positions = None
            axes      = None
            position_element = child(element,'atomic_positions')
            atom_elements    = children(position_element,'atom')
            if len(atom_elements)>0:
                labels = []
                rows   = []
                for atom in atom_elements:
                    coordinates = array(atom)
                    if coordinates is None or len(coordinates)<3:
                        rows = []
                        break
                    labels.append(atom.attrib.get('name',''))
                    rows.append(coordinates[:3])
                if len(rows)>0:
                    atoms     = np.array(labels,dtype=str)
                    positions = np.array(rows,dtype=float)
            cell = child(element,'cell')
            rows = [array(child(cell,name)) for name in ('a1','a2','a3')]
            if all(row is not None and len(row)>=3 for row in rows):
                axes = np.array([row[:3] for row in rows],dtype=float)
            return atoms,positions,axes
        #end def structure

        def complete_array(values,dtype=float):
            """Stack a sequence only when every entry is available."""
            if len(values)==0 or any(item is None for item in values):
                return None
            try:
                return np.array(values,dtype=dtype)
            except ValueError:
                return None
        #end def complete_array

        # version and calculation type
        creator = path(root,'general_info','creator')
        if creator is not None:
            self.version = creator.attrib.get('VERSION',creator.attrib.get('version'))
        self.calculation = value(path(root,'input','control_variables','calculation'))

        # principal XML sections
        input_element  = child(root,'input')
        output         = child(root,'output')
        band_structure = child(output,'band_structure')

        # exit status and total-energy components
        self.exit_status = value(child(root,'exit_status'))
        energy = child(output,'total_energy')
        energy_members = (
            ('total_energy','etot'),
            ('band_energy','eband'),
            ('hartree_energy','ehart'),
            ('xc_energy','etxc'),
            ('xc_potential_energy','vtxc'),
            ('ewald_energy','ewald'),
            )
        for member,name in energy_members:
            self[member] = value(child(energy,name))

        # electronic and structural convergence
        convergence = child(output,'convergence_info')
        scf_conv     = child(convergence,'scf_conv')
        self.scf_converged = value(child(scf_conv,'convergence_achieved'))
        self.scf_steps     = value(child(scf_conv,'n_scf_steps'))
        self.scf_error     = value(child(scf_conv,'scf_error'))
        opt_conv = child(convergence,'opt_conv')
        self.optimization_converged = value(child(opt_conv,'convergence_achieved'))
        self.optimization_steps     = value(child(opt_conv,'n_opt_steps'))
        self.gradient_norm          = value(child(opt_conv,'grad_norm'))

        # initial and final structures, forces, and stress
        initial_structure = child(input_element,'atomic_structure')
        final_structure   = child(output,'atomic_structure')
        self.initial_atoms,self.initial_positions,self.initial_axes = structure(initial_structure)
        self.atoms,self.positions,self.axes = structure(final_structure)
        for member,element in (
            ('initial_alat',initial_structure),
            ('alat',final_structure),
            ):
            text = None if element is None else element.attrib.get('alat')
            if text is not None and re.fullmatch(number_pattern,text.strip()) is not None:
                self[member] = float(text.replace('D','E').replace('d','e'))
        if self.initial_axes is not None:
            self.initial_volume = abs(np.linalg.det(self.initial_axes))
        if self.axes is not None:
            self.volume = abs(np.linalg.det(self.axes))
        self.forces = matrix(child(output,'forces'),3)
        self.stress = matrix(child(output,'stress'),3)

        # exchange-correlation functional and magnetization
        dft = child(output,'dft')
        self.functional = value(child(dft,'functional'))
        magnetization = child(output,'magnetization')
        self.spin_polarized         = value(child(magnetization,'lsda'))
        self.noncollinear           = value(child(magnetization,'noncolin'))
        self.spin_orbit             = value(child(magnetization,'spinorbit'))
        self.absolute_magnetization = value(child(magnetization,'absolute'))

        # electron, band, occupation, and Fermi-level counts and values
        self.nelectrons              = value(child(band_structure,'nelec'))
        self.nbands                  = value(child(band_structure,'nbnd'))
        self.nkpoints                = value(child(band_structure,'nks'))
        self.fermi_energy            = value(child(band_structure,'fermi_energy'))
        self.highest_occupied_level  = value(child(band_structure,'highestOccupiedLevel'))
        self.lowest_unoccupied_level = value(child(band_structure,'lowestUnoccupiedLevel'))
        self.occupations_kind        = value(child(band_structure,'occupations_kind'))
        lsda = value(child(band_structure,'lsda'))
        if lsda is not None:
            self.spin_polarized = lsda

        # k-points, weights, eigenvalues, occupations, and plane-wave counts
        records = children(band_structure,'ks_energies')
        coordinate_map = {}
        electronic     = []
        for record in records:
            k_element  = child(record,'k_point')
            coordinates = array(k_element)
            eigenvalues = array(child(record,'eigenvalues'))
            occupations = array(child(record,'occupations'))
            weight_text = '' if k_element is None else k_element.attrib.get('weight','')
            if (coordinates is None
                or len(coordinates)<3
                or eigenvalues is None
                or occupations is None
                or len(eigenvalues)!=len(occupations)
                or re.fullmatch(number_pattern,weight_text.strip()) is None
                ):
                continue
            coordinates = coordinates[:3]
            key = (
                tuple(np.round(coordinates,12))
                if self.spin_polarized else len(electronic)
                )
            if key not in coordinate_map:
                coordinate_map[key] = len(electronic)
                electronic.append(dotdict(
                    kpoint      = coordinates,
                    weight      = float(weight_text.replace('D','E').replace('d','e')),
                    eigenvalues = [],
                    occupations = [],
                    plane_waves = [],
                    ))
            entry = electronic[coordinate_map[key]]
            entry.eigenvalues.append(eigenvalues)
            entry.occupations.append(occupations)
            entry.plane_waves.append(value(child(record,'npw')))
        if len(electronic)>0:
            self.kpoints_rel = np.array([entry.kpoint for entry in electronic],dtype=float)
            self.kweights    = np.array([entry.weight for entry in electronic],dtype=float)
            nspin = 2 if self.spin_polarized else 1
            if all(len(entry.eigenvalues)==nspin for entry in electronic):
                eigenvalues = complete_array([
                    [entry.eigenvalues[spin] for entry in electronic]
                    for spin in range(nspin)
                    ])
                occupations = complete_array([
                    [entry.occupations[spin] for entry in electronic]
                    for spin in range(nspin)
                    ])
                plane_waves = complete_array([
                    [entry.plane_waves[spin] for entry in electronic]
                    for spin in range(nspin)
                    ])
                if nspin==1:
                    if eigenvalues is not None:
                        eigenvalues = eigenvalues[0]
                    if occupations is not None:
                        occupations = occupations[0]
                    if plane_waves is not None:
                        plane_waves = plane_waves[0]
                self.eigenvalues = eigenvalues
                self.occupations = occupations
                self.plane_waves = plane_waves

        # atomic species, masses, pseudopotentials, and valence charges
        species_elements = children(child(output,'atomic_species'),'species')
        if len(species_elements)>0:
            names       = []
            masses      = []
            pseudos     = []
            charges     = []
            complete_z  = True
            for species in species_elements:
                names.append(species.attrib.get('name',''))
                masses.append(value(child(species,'mass')))
                pseudo = child(species,'pseudo_file')
                pseudos.append(value(pseudo))
                zval = None if pseudo is None else pseudo.attrib.get('Zval',pseudo.attrib.get('zval'))
                if zval is None or re.fullmatch(number_pattern,zval.strip()) is None:
                    complete_z = False
                else:
                    charges.append(float(zval.replace('D','E').replace('d','e')))
            #end for
            self.species          = np.array(names,dtype=str)
            self.atomic_masses    = complete_array(masses)
            self.pseudopotentials = np.array(pseudos,dtype=str)
            if complete_z:
                self.valence_charges = np.array(charges,dtype=float)

        # plane-wave cutoffs, FFT grids, and reciprocal lattice
        basis = child(output,'basis_set')
        self.ecutwfc    = value(child(basis,'ecutwfc'))
        self.ecutrho    = value(child(basis,'ecutrho'))
        self.gamma_only = value(child(basis,'gamma_only'))
        self.ngm        = value(child(basis,'ngm'))
        self.ngms       = value(child(basis,'ngms'))
        self.npwx       = value(child(basis,'npwx'))
        for member,name in (
            ('fft_grid','fft_grid'),
            ('fft_smooth','fft_smooth'),
            ('fft_box','fft_box'),
            ):
            element = child(basis,name)
            if element is not None:
                dimensions = [
                    element.attrib.get(label)
                    for label in ('nr1','nr2','nr3')
                    ]
                if all(item is not None for item in dimensions):
                    self[member] = np.array(dimensions,dtype=int)
        reciprocal = child(basis,'reciprocal_lattice')
        rows = [array(child(reciprocal,name)) for name in ('b1','b2','b3')]
        if all(row is not None and len(row)>=3 for row in rows):
            self.reciprocal_lattice = np.array([row[:3] for row in rows],dtype=float)

        # Monkhorst-Pack grid and shift
        starting_kpoints = child(band_structure,'starting_k_points')
        if starting_kpoints is None:
            starting_kpoints = child(basis,'starting_k_points')
        monkhorst_pack = child(starting_kpoints,'monkhorst_pack')
        if monkhorst_pack is not None:
            grid  = [monkhorst_pack.attrib.get(f'nk{i}') for i in range(1,4)]
            shift = [monkhorst_pack.attrib.get(f'k{i}') for i in range(1,4)]
            if all(item is not None for item in grid):
                self.kgrid = np.array(grid,dtype=int)
            if all(item is not None for item in shift):
                self.kshift = np.array(shift,dtype=int)

        # crystal symmetry operations and equivalent atoms
        symmetries         = child(output,'symmetries')
        self.nsym          = value(child(symmetries,'nsym'))
        self.nrot          = value(child(symmetries,'nrot'))
        self.space_group   = value(child(symmetries,'space_group'))
        symmetry_elements = children(symmetries,'symmetry')
        names        = []
        classes      = []
        rotations    = []
        translations = []
        equivalents  = []
        for symmetry in symmetry_elements:
            info        = child(symmetry,'info')
            rotation    = matrix(child(symmetry,'rotation'),3)
            translation = array(child(symmetry,'fractional_translation'))
            equivalent  = array(child(symmetry,'equivalent_atoms'),dtype=int)
            if (info is None
                or rotation is None
                or rotation.shape!=(3,3)
                or translation is None
                or len(translation)<3
                or equivalent is None
                ):
                continue
            names.append(info.attrib.get('name',''))
            classes.append(info.attrib.get('class',''))
            rotations.append(rotation)
            translations.append(translation[:3])
            equivalents.append(equivalent)
        if len(rotations)>0:
            self.symmetry_names        = np.array(names,dtype=str)
            self.symmetry_classes      = np.array(classes,dtype=str)
            self.symmetry_rotations    = np.array(rotations,dtype=float)
            self.symmetry_translations = np.array(translations,dtype=float)
            self.equivalent_atoms      = complete_array(equivalents,dtype=int)

        # total execution timing
        timing = path(root,'timing_info','total')
        self.cpu_time  = value(child(timing,'cpu'))
        self.wall_time = value(child(timing,'wall'))

        # ionic-step energy, structure, force, stress, and convergence histories
        steps = children(root,'step')
        if len(steps)>0:
            trajectory = dotdict(
                energies              = [],
                band_energies         = [],
                hartree_energies      = [],
                xc_energies           = [],
                xc_potential_energies = [],
                ewald_energies        = [],
                positions             = [],
                axes                  = [],
                volumes               = [],
                forces                = [],
                stress                = [],
                scf_converged         = [],
                scf_steps             = [],
                scf_error             = [],
                )
            for step in steps:
                step_energy = child(step,'total_energy')
                trajectory.energies.append(value(child(step_energy,'etot')))
                trajectory.band_energies.append(value(child(step_energy,'eband')))
                trajectory.hartree_energies.append(value(child(step_energy,'ehart')))
                trajectory.xc_energies.append(value(child(step_energy,'etxc')))
                trajectory.xc_potential_energies.append(value(child(step_energy,'vtxc')))
                trajectory.ewald_energies.append(value(child(step_energy,'ewald')))
                _,positions,axes = structure(child(step,'atomic_structure'))
                trajectory.positions.append(positions)
                trajectory.axes.append(axes)
                trajectory.volumes.append(
                    None if axes is None else abs(np.linalg.det(axes))
                    )
                trajectory.forces.append(matrix(child(step,'forces'),3))
                trajectory.stress.append(matrix(child(step,'stress'),3))
                step_scf = child(step,'scf_conv')
                trajectory.scf_converged.append(value(child(step_scf,'convergence_achieved')))
                trajectory.scf_steps.append(value(child(step_scf,'n_scf_steps')))
                trajectory.scf_error.append(value(child(step_scf,'scf_error')))
            for member in trajectory.keys():
                self['trajectory_'+member] = complete_array(trajectory[member])

        # calculation-specific result membership updates
        if self.calculation in ('nscf','bands'):
            del self.forces
            del self.stress
        if self.calculation not in ('relax','vc-relax'):
            del self.optimization_converged
            del self.optimization_steps
            del self.gradient_norm
        if self.calculation not in ('relax','vc-relax','md','vc-md'):
            for name in tuple(self.keys()):
                if name.startswith('trajectory_'):
                    del self[name]
    #end def extract_results

#end class PwscfXmlData



class Pw2CasinoAnalyzer(DevBase):
    """Read and store physical data from PW2CASINO output.

    The auxiliary analyzer is kept separate from PWSCF text and XML results.
    Its data members remain ``None`` when the requested information is absent
    or the output file cannot be read.
    """

    def __init__(self,filepath):
        """Initialize empty PW2CASINO results."""
        self.K = None

        with open(filepath,'r') as fobj:
            lines = fobj.readlines()
        for line in lines:
            if 'Kinetic' in line:
                tokens = line.split()
                # Check whether the kinetic energy token is a complete numeric value.
                if len(tokens)>5 and re.fullmatch(number_pattern,tokens[5].strip()) is not None:
                    self.K = float(tokens[5].replace('D','E').replace('d','e'))
    #end def __init__

#end class Pw2CasinoAnalyzer



class PwscfAnalyzer(SimulationAnalyzer):
    """Analyze output produced by Quantum ESPRESSO PWscf calculations.

    The analyzer coordinates PWSCF text, auxiliary, and XML readers across
    common run modes. It also provides summaries and visualizations useful
    for inspecting molecular dynamics and electronic structure.

    Parameters
    ----------
    arg0 : Simulation or str or os.PathLike or None, optional
        PWSCF simulation to analyze, or path to a calculation directory,
        input file, or output file. If ``None``, an unconfigured analyzer is
        created.
    infile_name : str, optional
        Name of the PWSCF input file within the calculation directory.
    outfile_name : str, optional
        Name of the PWSCF output file. It is inferred from ``infile_name``
        when possible.
    pw2c_outfile_name : str, optional
        Name of an accompanying PW2CASINO output file.
    analyze : bool, optional
        If ``True``, parse the available log, XML, and auxiliary output during
        initialization.
    md_only : bool, optional
        Limit text-output analysis to molecular-dynamics data. Available XML
        output is still parsed.

    Attributes
    ----------
    path : str
        Directory containing the PWSCF input and output files.
    abspath : str
        Absolute path to the calculation directory.
    infile_name, outfile_name : str or None
        Names of the PWSCF input and text-output files.
    pw2c_outfile_name : str or None
        Name of the optional PW2CASINO output file.
    info : obj
        General analyzer metadata, including the ``md_only`` setting.
    input : PwscfInput or None
        Parsed PWSCF input when an input file is available.
    simulation_structure : Structure
        Input structure supplied by a ``Simulation``. This member is bound
        only when the analyzer is constructed from a simulation.
    results_out : PwscfOutData or None
        Parsed PWSCF text-output data. It is ``None`` until analysis is
        performed.
    results_xml : PwscfXmlData or obj or None
        Parsed schema or legacy XML data. It remains ``None`` when XML output
        is absent or cannot be read.
    pw2casino : Pw2CasinoAnalyzer or None
        Parsed PW2CASINO data when an auxiliary output file is requested.

    Methods
    -------
    initial_structure(units='A') : Structure or None
        Initial structure in Angstrom (``'A'``) or bohr (``'B'``).
    energy(units='Ha') : float or numpy.floating or None
        Final total energy in ``'eV'``, ``'Ha'``, or ``'Ry'``.
    kpoints(units='B') : numpy.ndarray or None
        Cartesian k-points in inverse Angstrom or inverse bohr with shape
        ``(nkpoints, 3)``. The unit argument is ``'A'`` or ``'B'``.
    kweights() : numpy.ndarray or None
        Dimensionless k-point weights with shape ``(nkpoints,)``.
    eigenvalues(units='eV') : numpy.ndarray or None
        Kohn--Sham eigenvalues in ``'eV'``, ``'Ha'``, or ``'Ry'``. The leading
        dimension has length ``nkpoints``; remaining dimensions represent
        spin, when present, and bands.
    occupations() : numpy.ndarray or None
        Dimensionless Kohn--Sham occupations with the same layout as the
        eigenvalue array.
    Ef(units='eV') : float or numpy.floating or None
        Final Fermi energy in ``'eV'``, ``'Ha'``, or ``'Ry'``.
    Evbm(units='eV') : float or numpy.floating or None
        Final valence-band maximum in selected energy units.
    Ecbm(units='eV') : float or numpy.floating or None
        Final conduction-band minimum in selected energy units.
    band_gap(units='eV') : float or numpy.floating or None
        Fundamental electronic band gap in selected energy units.
    fractional_occs() : bool or None
        Whether any occupation differs from both empty and full by more than
        ``1e-3``.
    relaxed_structure(units='A') : Structure or None
        Final relaxed structure in Angstrom (``'A'``) or bohr (``'B'``).
    forces(units='eV/A') : numpy.ndarray or None
        Ionic-force history with shape ``(nsteps, natoms, 3)``. Available
        units are ``'eV/A'``, ``'Ry/B'``, and ``'Ha/B'``.
    stress(units='GPa') : numpy.ndarray or None
        Stress-tensor history with shape ``(nsteps, 3, 3)``. Available units
        are ``'Pa'``, ``'bar'``, ``'kbar'``, ``'Mbar'``, ``'GPa'``, and
        ``'atm'``.
    pressure(units='GPa') : float or numpy.floating or None
        Final hydrostatic pressure in the units accepted by ``stress``.

    Raises
    ------
    FileNotFoundError
        If a supplied path, input file, output file, or requested PW2CASINO
        file does not exist.
    RuntimeError
        If a supplied file cannot be identified as input or output.

    Notes
    -----
    Log and XML output are both parsed automatically when present. Query
    methods prefer schema XML data and fall back to text-output data. A query
    returns ``None`` when its quantity applies to the detected calculation but
    was not parsed. Calling a query before analysis or for an unsupported
    calculation raises ``RuntimeError``. Supplying unsupported units raises
    ``ValueError``.

    Initial structure, k-point, and electronic quantities apply to all
    supported calculation modes. Energy applies to all modes except
    ``bands``; relaxed structures apply only to ``relax`` and ``vc-relax``;
    forces, stress, and pressure apply to ``scf``, ``relax``, ``vc-relax``,
    ``md``, and ``vc-md``.
    """

    all_modes        = frozenset({'scf','nscf','bands','relax','vc-relax','md','vc-md'})
    energy_modes     = frozenset({'scf','nscf','relax','vc-relax','md','vc-md'})
    electronic_modes = all_modes
    relaxation_modes = frozenset({'relax','vc-relax'})
    force_modes      = frozenset({'scf','relax','vc-relax','md','vc-md'})
    pressure_units   = MappingProxyType({
        'Pa'   : 1.0,
        'bar'  : 1e5,
        'kbar' : 1e8,
        'Mbar' : 1e11,
        'GPa'  : 1e9,
        'atm'  : 1.01325e5,
        })


    def _require_supported(self,quantity,modes):
        """Require analyzed output and a calculation supporting the quantity."""
        if 'results_out' not in self or self.results_out is None:
            msg = f'PWSCF quantity "{quantity}" is unavailable because output has not been analyzed'
            raise RuntimeError(msg)
        calculation = self.results_out.calculation
        if isinstance(self.results_xml,PwscfXmlData):
            if self.results_xml.calculation is not None:
                calculation = self.results_xml.calculation
        if calculation not in modes:
            msg = f'PWSCF quantity "{quantity}" is not supported for calculation "{calculation}"'
            raise RuntimeError(msg)
    #end def _require_supported


    def initial_structure(self,units='A'):
        """Return the initial ``Structure`` in Angstrom or bohr."""
        self._require_supported('initial_structure',self.all_modes)
        if units not in {'A','B'}:
            raise ValueError('initial_structure units must be one of: A, B')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if (xml is not None
            and xml.initial_atoms is not None
            and xml.initial_positions is not None
            and xml.initial_axes is not None
            ):
            structure = Structure(
                axes    = xml.initial_axes.copy(),
                elem    = xml.initial_atoms.copy(),
                pos     = xml.initial_positions.copy(),
                units   = 'B',
                rescale = False,
                )
        elif 'simulation_structure' in self and self.simulation_structure is not None:
            structure = deepcopy(self.simulation_structure)
        elif (self.input is not None
            and 'system' in self.input
            and 'ibrav' in self.input.system
            and self.input.system.ibrav==0
            and 'atomic_positions' in self.input
            and 'cell_parameters' in self.input
            and 'k_points' in self.input
            ):
            input_data = deepcopy(self.input)
            system     = input_data.system
            if ('celldm(1)' not in system
                and 'celldm' in system
                and 1 in system.celldm
                ):
                system['celldm(1)'] = system.celldm[1]
            structure = input_data.return_system(structure_only=True)
        else:
            return None
        structure.change_units(units)
        return structure
    #end def initial_structure


    def energy(self,units='Ha'):
        """Return the final total energy in eV, Hartree, or Rydberg."""
        self._require_supported('energy',self.energy_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('energy units must be one of: eV, Ha, Ry')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.total_energy is not None:
            return convert(xml.total_energy,'Ha',units)
        if 'E' in self.results_out and self.results_out.E is not None:
            return convert(self.results_out.E,'Ry',units)
        return None
    #end def energy


    def kpoints(self,units='B'):
        """Return Cartesian k-points in inverse Angstrom or inverse bohr."""
        self._require_supported('kpoints',self.all_modes)
        if units not in {'A','B'}:
            raise ValueError('kpoints units must be one of: A, B')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        kpoints = None
        if xml is not None and xml.kpoints_rel is not None:
            scale = xml.alat if xml.alat is not None else xml.initial_alat
            if scale is not None and scale>0:
                kpoints = xml.kpoints_rel*2*np.pi/scale
        if (kpoints is None
            and self.results_out.kpoints_cart is not None
            and self.input is not None
            and 'system' in self.input
            ):
            system = self.input.system
            scale  = None
            if 'celldm(1)' in system:
                scale = system['celldm(1)']
            elif 'celldm' in system and 1 in system.celldm:
                scale = system.celldm[1]
            if scale is not None:
                kpoints = self.results_out.kpoints_cart*2*np.pi/scale
        if kpoints is None and self.results_out.kpoints_unit is not None:
            structure = self.initial_structure('B')
            if structure is not None:
                kpoints = np.dot(self.results_out.kpoints_unit,structure.kaxes)
        if kpoints is None:
            return None
        return kpoints*convert(1.0,units,'B')
    #end def kpoints


    def kweights(self):
        """Return dimensionless k-point weights, or ``None`` if unavailable."""
        self._require_supported('kweights',self.all_modes)
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.kweights is not None:
            return xml.kweights
        return self.results_out.kweights
    #end def kweights


    def _log_band_values(self,name):
        """Return complete k-point-major band values from text output."""
        bands = self.results_out.bands
        if bands is None:
            return None
        channels = []
        for channel in (bands.up,bands.down):
            if len(channel)==0:
                continue
            values = [
                np.asarray(channel[index][name],dtype=float)
                for index in sorted(channel.keys())
                ]
            if len(values)==0 or len(values[0])==0:
                return None
            shape = values[0].shape
            if any(value.shape!=shape for value in values):
                return None
            channels.append(np.array(values,dtype=float))
        if len(channels)==0:
            return None
        if len(channels)==1:
            return channels[0]
        if channels[0].shape!=channels[1].shape:
            return None
        return np.stack(channels,axis=1)
    #end def _log_band_values


    def eigenvalues(self,units='eV'):
        """Return k-point-major eigenvalues in eV, Hartree, or Rydberg."""
        self._require_supported('eigenvalues',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('eigenvalues units must be one of: eV, Ha, Ry')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.eigenvalues is not None:
            values = xml.eigenvalues
            if values.ndim==3:
                values = np.transpose(values,(1,0,2))
            return convert(values,'Ha',units)
        values = self._log_band_values('eigs')
        if values is None:
            return None
        return convert(values,'eV',units)
    #end def eigenvalues


    def occupations(self):
        """Return the dimensionless k-point-major occupation array."""
        self._require_supported('occupations',self.electronic_modes)
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.occupations is not None:
            values = xml.occupations
            if values.ndim==3:
                values = np.transpose(values,(1,0,2))
            return values
        return self._log_band_values('occs')
    #end def occupations


    def Ef(self,units='eV'):
        """Return the final Fermi energy in eV, Hartree, or Rydberg."""
        self._require_supported('Ef',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Ef units must be one of: eV, Ha, Ry')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.fermi_energy is not None:
            return convert(xml.fermi_energy,'Ha',units)
        if self.results_out.Ef is not None:
            return convert(self.results_out.Ef,'eV',units)
        return None
    #end def Ef


    def Evbm(self,units='eV'):
        """Return the final valence-band maximum in selected energy units."""
        self._require_supported('Evbm',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Evbm units must be one of: eV, Ha, Ry')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.highest_occupied_level is not None:
            return convert(xml.highest_occupied_level,'Ha',units)
        bands = self.results_out.bands
        if bands is not None and 'vbm' in bands:
            return convert(bands.vbm.energy,'eV',units)
        return None
    #end def Evbm


    def Ecbm(self,units='eV'):
        """Return the final conduction-band minimum in selected energy units."""
        self._require_supported('Ecbm',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Ecbm units must be one of: eV, Ha, Ry')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if xml is not None and xml.lowest_unoccupied_level is not None:
            return convert(xml.lowest_unoccupied_level,'Ha',units)
        bands = self.results_out.bands
        if bands is not None and 'cbm' in bands:
            return convert(bands.cbm.energy,'eV',units)
        return None
    #end def Ecbm


    def band_gap(self,units='eV'):
        """Return the fundamental band gap in eV, Hartree, or Rydberg."""
        self._require_supported('band_gap',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('band_gap units must be one of: eV, Ha, Ry')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if (xml is not None
            and xml.highest_occupied_level is not None
            and xml.lowest_unoccupied_level is not None
            ):
            gap = xml.lowest_unoccupied_level-xml.highest_occupied_level
            return convert(gap,'Ha',units)
        bands = self.results_out.bands
        if bands is not None and 'vbm' in bands and 'cbm' in bands:
            gap = bands.cbm.energy-bands.vbm.energy
            return convert(gap,'eV',units)
        return None
    #end def band_gap


    def fractional_occs(self):
        """Whether any occupation differs from empty or full by over ``1e-3``."""
        self._require_supported('fractional_occs',self.electronic_modes)
        occupations = self.occupations()
        if occupations is None:
            return None
        tolerance = 1e-3
        empty     = np.isclose(occupations,0.0,rtol=0.0,atol=tolerance)
        full      = np.isclose(occupations,1.0,rtol=0.0,atol=tolerance)
        return bool(np.any(~(empty|full)))
    #end def fractional_occs


    def relaxed_structure(self,units='A'):
        """Return the final relaxed ``Structure`` in Angstrom or bohr."""
        self._require_supported('relaxed_structure',self.relaxation_modes)
        if units not in {'A','B'}:
            raise ValueError('relaxed_structure units must be one of: A, B')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        if (xml is not None
            and xml.atoms is not None
            and xml.positions is not None
            and xml.axes is not None
            ):
            structure = Structure(
                axes    = xml.axes.copy(),
                elem    = xml.atoms.copy(),
                pos     = xml.positions.copy(),
                units   = 'B',
                rescale = False,
                )
        elif ('relax_structures' in self.results_out
            and self.results_out.relax_structures is not None
            and len(self.results_out.relax_structures)>0
            ):
            structures = self.results_out.relax_structures
            result     = structures[max(structures.keys())]
            initial    = self.initial_structure('B')
            axes       = result.axes if 'axes' in result else None
            if axes is None and initial is not None:
                axes = initial.axes
            if axes is None:
                return None
            structure = Structure(
                axes    = np.asarray(axes,dtype=float),
                elem    = np.asarray(result.atoms,dtype=str),
                pos     = np.asarray(result.positions,dtype=float),
                units   = 'B',
                rescale = False,
                )
        else:
            return None
        structure.change_units(units)
        return structure
    #end def relaxed_structure


    def forces(self,units='eV/A'):
        """Return ionic forces in ``eV/A``, ``Ry/B``, or ``Ha/B``."""
        self._require_supported('forces',self.force_modes)
        if units not in {'eV/A','Ry/B','Ha/B'}:
            raise ValueError('forces units must be one of: eV/A, Ry/B, Ha/B')
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        values        = None
        source_energy = None
        if xml is not None:
            if 'trajectory_forces' in xml and xml.trajectory_forces is not None:
                values = xml.trajectory_forces
            elif 'forces' in xml and xml.forces is not None:
                values = xml.forces[np.newaxis,:,:]
            if values is not None:
                source_energy = 'Ha'
        if values is None and self.results_out.forces is not None:
            values        = self.results_out.forces
            source_energy = 'Ry'
        if values is None:
            return None
        energy_units,length_units = units.split('/')
        factor = convert(1.0,source_energy,energy_units)/convert(1.0,'B',length_units)
        return values*factor
    #end def forces


    def stress(self,units='GPa'):
        """Return the stress-tensor history in selected pressure units."""
        self._require_supported('stress',self.force_modes)
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            raise ValueError('stress units must be one of: '+supported)
        xml = self.results_xml if isinstance(self.results_xml,PwscfXmlData) else None
        values          = None
        source_pressure = None
        if xml is not None:
            if 'trajectory_stress' in xml and xml.trajectory_stress is not None:
                values = xml.trajectory_stress
            elif 'stress' in xml and xml.stress is not None:
                values = xml.stress[np.newaxis,:,:]
            if values is not None:
                hartree = convert(1.0,'Ha','J')
                bohr    = convert(1.0,'B','m')
                source_pressure = hartree/bohr**3
        if values is None and self.results_out.stress is not None:
            rows = np.asarray(self.results_out.stress,dtype=float)
            if rows.ndim!=2 or rows.shape[1]<6 or len(rows)%3!=0:
                return None
            values          = rows[:,3:6].reshape(-1,3,3)
            source_pressure = 1e8
        if values is None:
            return None
        return values*source_pressure/self.pressure_units[units]
    #end def stress


    def pressure(self,units='GPa'):
        """Return the final hydrostatic pressure in selected pressure units."""
        self._require_supported('pressure',self.force_modes)
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            raise ValueError('pressure units must be one of: '+supported)
        stress = self.stress(units)
        if stress is not None and len(stress)>0:
            return np.trace(stress[-1])/3.0
        if self.results_out.pressure is not None:
            return self.results_out.pressure*1e8/self.pressure_units[units]
        return None
    #end def pressure


    def __init__(
        self,
        arg0              = None,
        infile_name       = None,
        outfile_name      = None,
        pw2c_outfile_name = None,
        *,
        analyze           = False,
        md_only           = False,
        ):
        """Initialize an analyzer for a PWSCF simulation or output path."""
        if isinstance(arg0,Simulation):
            sim                       = arg0
            path                      = sim.locdir
            infile_name               = sim.infile
            outfile_name              = sim.outfile
            self.simulation_structure = sim.system.structure
        elif arg0 is not None:
            path = path_string(arg0)
            if not os.path.exists(path):
                msg = (
                    'path to QE data does not exist\n'
                    f'path provided: {path}'
                    )
                raise FileNotFoundError(msg)
            if os.path.isfile(path):
                filepath = path
                path,filename = os.path.split(filepath)
                if filename.endswith('.in'):
                    infile_name = filename
                elif filename.endswith('.out'):
                    outfile_name = filename
                else:
                    msg = (
                        'could not determine whether file is QE input or output\n'
                        f'file provided: {filepath}'
                        )
                    raise RuntimeError(msg)
            if outfile_name is None and infile_name is not None:
                outfile_name = f"{infile_name.rsplit('.',1)[0]}.out"
        else:
            return

        inp = None
        if infile_name is not None:
            infile = os.path.join(path,infile_name)
            if os.path.isfile(infile):
                inp = PwscfInput(infile)
            else:
                msg = f'PWSCF input file is not available\nfile not found: {infile}'
                raise FileNotFoundError(msg)

        self.info              = obj(md_only=md_only)
        self.infile_name       = infile_name
        self.outfile_name      = outfile_name
        self.path              = path
        self.abspath           = os.path.abspath(path)
        self.pw2c_outfile_name = pw2c_outfile_name
        self.input             = inp
        self.results_out       = None
        self.results_xml       = None
        self.pw2casino         = None
        if analyze:
            self.analyze()
    #end def __init__


    def analyze(self):
        """Analyze available PWSCF text, XML, and auxiliary output."""
        self.results_out = None
        self.results_xml = None
        self.pw2casino   = None
        if ('path' not in self
            or 'outfile_name' not in self
            or self.outfile_name is None
            ):
            raise RuntimeError('PWSCF output file name is not available')
        outfile = os.path.join(self.path,self.outfile_name)
        if not os.path.isfile(outfile):
            msg = f'PWSCF output file is not available\nfile not found: {outfile}'
            raise FileNotFoundError(msg)
        self.results_out = PwscfOutData(outfile,md_only=self.info.md_only)
        self.analyze_xml()
        if self.info.md_only:
            return
        if self.pw2c_outfile_name is not None:
            filepath = os.path.join(self.path,self.pw2c_outfile_name)
            if os.path.isfile(filepath):
                self.pw2casino = Pw2CasinoAnalyzer(filepath)
            else:
                msg = f'PW2CASINO output file is not available\nfile not found: {filepath}'
                raise FileNotFoundError(msg)
    #end def analyze


    def analyze_xml(self):
        """Locate and parse schema or legacy PWscf XML output."""
        self.results_xml = None
        schema_file = None
        legacy_file = None
        legacy_dir  = None
        if ('input' in self
            and self.input is not None
            and 'control' in self.input
            and 'outdir' in self.input.control
            and 'prefix' in self.input.control
            ):
            cont         = self.input.control
            savedir      = os.path.join(self.path,cont.outdir,f'{cont.prefix}.save')
            schema_path  = os.path.join(savedir,'data-file-schema.xml')
            legacy_path  = os.path.join(savedir,'data-file.xml')
            fallback_dir = os.path.join(self.path,cont.outdir)
            fallback     = os.path.join(fallback_dir,f'{cont.prefix}.xml')
            if os.path.isfile(schema_path):
                schema_file = schema_path
            elif os.path.isfile(legacy_path):
                legacy_file = legacy_path
                legacy_dir  = savedir
            elif os.path.isfile(fallback):
                legacy_file = fallback
                legacy_dir  = fallback_dir

        if schema_file is None and legacy_file is None:
            schema_candidates = sorted(set(
                glob(os.path.join(self.path,'*.save','data-file-schema.xml'))
                + glob(os.path.join(self.path,'*','*.save','data-file-schema.xml'))
                ))
            if len(schema_candidates)==1:
                schema_file = schema_candidates[0]
            else:
                legacy_candidates = sorted(set(
                    glob(os.path.join(self.path,'*.save','data-file.xml'))
                    + glob(os.path.join(self.path,'*','*.save','data-file.xml'))
                    + glob(os.path.join(self.path,'*.xml'))
                    + glob(os.path.join(self.path,'*','*.xml'))
                    ))
                if len(legacy_candidates)==1:
                    legacy_file = legacy_candidates[0]
                    legacy_dir  = os.path.dirname(legacy_file)

        if schema_file is not None:
            results_xml = PwscfXmlData(schema_file)
            if results_xml.parse_failed:
                return
            self.results_xml = results_xml
        else:
            if legacy_file is None:
                return
            results_xml = PwscfXmlData(legacy_file)
            if results_xml.parse_failed:
                return
            if results_xml.calculation is not None:
                self.results_xml = results_xml
                return
            data = read_qexml(legacy_file)
            self.results_xml = obj(data=None,kpoints=None,failed=False)
            self.analyze_legacy_xml(data,legacy_dir)
            if self.results_xml.failed:
                self.results_xml = None
    #end def analyze_xml


    def analyze_legacy_xml(self,data,datadir):
        """Extract k-point and orbital data from legacy PWscf XML."""
        def object_path(value,*names):
            """Return a nested value, or None when its path is incomplete."""
            for name in names:
                if value is None or name not in value:
                    return None
                value = value[name]
            return value
        #end def object_path

        kpdata = object_path(data,'root','eigenvalues','k_point')
        if kpdata is None:
            self.results_xml.update(data=data,kpoints=obj())
            self.results_xml.failed = True
            return
        kpoints = obj()
        for ki,kpd in kpdata.items():
            if 'k_point_coords' not in kpd or 'weight' not in kpd or 'datafile' not in kpd:
                self.results_xml.failed = True
                continue
            kp = obj(kpoint=kpd.k_point_coords,weight=kpd.weight)
            kpoints[ki] = kp
            for si,dfile in kpd.datafile.items():
                efilepath = os.path.join(datadir,dfile.iotk_link)
                if not os.path.isfile(efilepath):
                    self.results_xml.failed = True
                    continue
                try:
                    edata = read_qexml(efilepath)
                except Exception:
                    self.results_xml.failed = True
                    continue
                eunits      = object_path(edata,'root','units_for_energies','units')
                eigenvalues = object_path(edata,'root','eigenvalues')
                occupations = object_path(edata,'root','occupations')
                if eunits is None or eigenvalues is None or occupations is None:
                    self.results_xml.failed = True
                    continue
                units = dict(ha='Ha',ry='Ry',ev='eV').get(eunits.lower()[:2],'Ha')
                spin  = obj(units=units,eigenvalues=eigenvalues,occupations=occupations)
                if si==1:
                    kp.up = spin
                elif si==2:
                    kp.down = spin
        self.results_xml.update(data=data,kpoints=kpoints)
    #end def analyze_legacy_xml






    def md_statistics(self,equil=None,autocorr=None):
        """Calculate summary statistics for molecular-dynamics histories."""
        if ('results_out' not in self
            or self.results_out is None
            or 'md_data' not in self.results_out
            or self.results_out.md_data is None
            ):
            return
        return self.results_out.md_statistics(equil,autocorr)
    #end def md_statistics


    def md_plots(self,*,show=True):
        """Plot energy, temperature, and pressure histories from dynamics."""
        if ('results_out' not in self
            or self.results_out is None
            or 'md_data' not in self.results_out
            or self.results_out.md_data is None
            ):
            return
        md = self.results_out.md_data

        import matplotlib.pyplot as plt
        fig = plt.figure()

        plt.subplot(3,1,1)
        plt.plot(md.time,md.total_energy-md.total_energy[0],label='Etot')
        plt.plot(md.time,md.kinetic_energy-md.kinetic_energy[0],label='Ekin')
        plt.plot(md.time,md.potential_energy-md.potential_energy[0],label='Epot')
        plt.ylabel('E (Ryd)')
        plt.legend()

        plt.subplot(3,1,2)
        plt.plot(md.time,md.temperature)
        plt.ylabel('T (K)')

        plt.subplot(3,1,3)
        plt.plot(md.time,md.pressure)
        plt.ylabel('P (kbar)')
        plt.xlabel('time (ps)')

        if show:
            plt.show()

        return fig
    #end def md_plots


    def make_movie(self,filename,filepath=None):
        """Write relaxed or dynamic structures as a tiled XYZ movie."""
        if ('results_out' not in self
            or self.results_out is None
            or 'relax_structures' not in self.results_out
            or self.results_out.relax_structures is None
            or 'input' not in self
            or self.input is None
            ):
            return
        structures = self.results_out.relax_structures
        filepath = self.abspath if filepath is None else filepath
        filepath = os.path.join(filepath,filename)

        movie         = ''
        alat_angstrom = convert(self.input.system['celldm(1)'],'B','A')
        cell          = self.input.cell_parameters.vectors
        for structure in structures.values():
            tiled = Structure(
                elem  = structure.atoms,
                pos   = structure.positions,
                axes  = cell,
                scale = alat_angstrom,
                units = 'A',
                ).tile(2,2,2)
            movie += tiled.write_xyz()
        with open(filepath,'w') as fobj:
            fobj.write(movie)
    #end def make_movie


    def plot_bandstructure(
        self,
        filename      = None,
        filepath      = None,
        max_min_e     = None,
        *,
        show          = False,
        save          = True,
        show_vbm_cbm  = True,
        k_labels      = None,
        ):
        """Plot the analyzed band structure along its reciprocal-space path."""
        import matplotlib.pyplot as plt
        if ('results_out' not in self
            or self.results_out is None
            or 'input' not in self
            or self.input is None
            or 'system' not in self.input
            or 'nbnd' not in self.input.system
            ):
            return
        bands = self.results_out.bands
        if bands is None:
            return
        params = {
            'legend.fontsize'      : 14,
            'figure.facecolor'     : 'white',
            'figure.subplot.hspace': 0.,
            'axes.labelsize'       : 16,
            'xtick.labelsize'      : 14,
            'ytick.labelsize'      : 14,
            }
        plt.rcParams.update(params)
        filename = 'band_structure.pdf' if filename is None else filename
        filepath = os.path.join(self.abspath if filepath is None else filepath,filename)

        plt.figure()
        ax     = plt.gca()
        nbands = self.input.system.nbnd

        if k_labels is None:
            structure = self.initial_structure()
            if structure is None:
                return
            kpath  = get_kpath(structure=structure,check_standard=False)
            x      = kpath['explicit_path_linearcoords']
            labels = list(kpath['explicit_kpoints_labels'])
        else:
            if self.results_out.kpoints_cart is None:
                return
            labels = list(k_labels)
            # Calculate linear coordinates from self.results_out.kpoints_cart
            x          = []
            prev_label = ''
            ref_kpt    = self.results_out.kpoints_cart[0]
            lincoord   = 0.0
            for kpt_idx,kpt in enumerate(self.results_out.kpoints_cart):
                curr_label = labels[kpt_idx]
                if (curr_label != '' and prev_label == '') or curr_label == '':
                    lincoord += np.linalg.norm(kpt-ref_kpt)
                ref_kpt = kpt
                x.append(lincoord)
                prev_label = curr_label
        for band_channel,color in ((bands.up,'k'),(bands.down,'r')):
            if len(band_channel)==0:
                continue
            for nb in range(nbands):
                eigs = np.array([band.eigs[nb] for band in band_channel],dtype=float)
                y = eigs - bands.vbm.energy
                plt.plot(x,y,color)
        for ln,li in enumerate(labels):
            if li != '':
                plt.axvline(x[ln],ymin=-100,ymax=100,linewidth=3,color='k')
                labels[ln] = r'$\Gamma$' if li=='GAMMA' else f'${li}$'
                if ln>0 and labels[ln-1]!='':
                    labels[ln] = f'{labels[ln-1]}|{labels[ln]}'
                    labels[ln-1] = ''

        plt.xlim([np.min(x),np.max(x)])
        plt.ylim((-5,+5) if max_min_e is None else max_min_e)
        plt.ylabel('Energy (eV)')
        plt.xticks(x,labels)
        ax.tick_params(axis='x',which='both',length=0)
        ax.tick_params(axis='x',which='both',pad=10)
        if show_vbm_cbm:
            vbm = bands.vbm
            cbm = bands.cbm
            for kn,ki in enumerate(bands.up):
                if (vbm.kpoint_rel==ki.kpoint_rel).all():
                    plt.scatter(x[kn],0,c='green',s=100)
                if (cbm.kpoint_rel==ki.kpoint_rel).all():
                    plt.scatter(x[kn],cbm.energy-vbm.energy,c='r',s=100)
        if save:
            plt.savefig(filepath,format='pdf',bbox_inches='tight')
        if show:
            plt.show()
    #end def plot_bandstructure

#end class PwscfAnalyzer
        
