##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf_analyzer.py                                                 #
#    Supports data analysis for PWSCF output.  Can handle log file   #
#    and legacy XML output.                                          #
#                                                                    #
#  Content summary:                                                  #
#    PwscfOutData                                                    #
#      Reads and stores physical data from PWSCF log output.         #
#                                                                    #
#    PwscfAnalyzer                                                   #
#      SimulationAnalyzer class for PWSCF.                           #
#      Coordinates text and legacy XML output analysis.              #
#      Can also read data-file.xml.  See pwscf_data_reader.py.       #
#                                                                    #
#====================================================================#


import os
import re
from copy import deepcopy
from glob import glob
from types import MappingProxyType

import numpy as np

from .developer import DevBase, obj
from .pwscf_data_reader import read_qexml
from .pwscf_input import PwscfInput
from .simulation import Simulation, SimulationAnalyzer
from .structure import Structure
from .unit_converter import convert
from .utilities import path_string

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


# Match a PWSCF total-energy assignment including its Rydberg unit.  Both SCF
# iteration lines (``total energy = -168.12345678 Ry``) and completed lines
# (``! total energy = -1.0D+02 Ry``) match.  Lines in eV, missing ``=``, or
# reporting ``one-electron contribution`` do not match this field.
total_energy_pattern = (
    rf'\btotal\s+energy\s*=\s*(?P<energy>{number_pattern})\s+Ry\b'
    )


# Match pressure from a PWSCF total-stress heading.  It accepts forms such as
# ``total stress (Ry/bohr**3) (kbar) P= -170.96`` and ``P = 1.2D+03``.
# Lowercase ``p=``, the word ``pressure``, and a missing assignment operator
# are rejected because they do not have the exact PWSCF field label.
pressure_pattern = rf'\bP\s*=\s*(?P<pressure>{number_pattern})(?=\s|$)'


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


# Match a complete six-column stress row: three Ry/bohr**3 entries followed by
# three kbar entries.  Examples include ``-.001 0.0 .001 -147.1 0.0 147.1``
# and exponent-form rows.  Five-column rows, labeled rows, and concatenated
# columns are rejected because they cannot form the two complete tensors.
stress_row_pattern = (
    rf'^\s*(?P<sxx>{number_pattern})\s+(?P<sxy>{number_pattern})\s+'
    rf'(?P<sxz>{number_pattern})\s+(?P<kxx>{number_pattern})\s+'
    rf'(?P<kxy>{number_pattern})\s+(?P<kxz>{number_pattern})(?=\s|$)'
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

    Parameters
    ----------
    filepath : str or os.PathLike
        Path to the PWSCF text-output file.

    Attributes
    ----------
    calculation : {'scf', 'nscf', 'relax', 'vc-relax'}
        Calculation type inferred from the output text.
    Ef : float or None
        Final Fermi energy in eV.
    fermi_energies : numpy.ndarray or None
        One-dimensional history of Fermi energies in eV.
    bands : obj or None
        Spin-resolved band records.  The ``up`` and ``down`` members map
        k-point indices to objects containing ``eigs`` and ``occs`` arrays.
        Optional ``vbm`` and ``cbm`` members contain their respective
        energies.
    kpoints_cart, kpoints_unit : numpy.ndarray or None
        Cartesian and crystal k-point arrays with shape ``(nkpoints, 3)``.
    kweights : numpy.ndarray or None
        One-dimensional k-point weight array.  Present for every calculation
        type.
    E : float or None
        Final total energy in Ry for ``scf``, ``relax``, and ``vc-relax``.
    pressure : float or None
        Final pressure in kbar.
    stress : numpy.ndarray or None
        Stress-tensor history in kbar with shape ``(nsteps, 3, 3)``.
    forces : numpy.ndarray or None
        Atomic-force history with shape ``(nsteps, natoms, 3)`` in Ry/bohr.
    relax_structures : obj or None
        Integer-indexed structure records containing atom labels, Cartesian
        positions, and, when reported, cell axes.  Present for relaxation
        calculations.

    Notes
    -----
    An applicable attribute remains ``None`` when its record is absent or
    cannot be parsed.  Attributes that do not apply to the inferred
    calculation type are removed from the object.
    """

    def __init__(self,filepath):
        """Read a PWSCF log and initialize its accessible physical data."""
        self.calculation = None

        # all calculation types
        self.Ef                = None
        self.fermi_energies    = None
        self.bands             = None
        self.kpoints_cart      = None
        self.kpoints_unit      = None
        self.kweights          = None
        # scf/relax/vc-relax
        self.E                 = None
        self.pressure          = None
        self.stress            = None
        self.forces            = None
        # relax/vc-relax
        self.relax_structures  = None

        with open(filepath,'r') as fobj:
            lines = fobj.read().splitlines()
        # read the calculation type
        self.read_calculation(lines)
        # remove unused attributes, depending on the calculation type
        if self.calculation=='nscf':
            for name in {  # noqa: PLC0208
                'E','pressure','stress','forces',
                }:
                del self[name]
            #end for
        if self.calculation in {'scf','nscf'}:
            del self.relax_structures
        # all calculations
        self.read_fermi_energies(lines)
        self.read_kpoints(lines)
        self.read_bands(lines)
        # all but nscf
        if self.calculation in {'scf','relax','vc-relax'}:
            self.read_energies(lines)
            self.read_pressure(lines)
            self.read_stress(lines)
            self.read_forces(lines)
        # relaxation calculations
        if self.calculation in {'relax','vc-relax'}:
            self.read_structures(lines)

    #end def __init__


    def read_calculation(self,lines):
        """Infer and bind the PWSCF calculation type from log records."""
        has_cell     = any(line.strip().startswith('CELL_PARAMETERS') for line in lines)
        has_bfgs     = any('BFGS Geometry Optimization' in line for line in lines)
        has_band_run = any('Band Structure Calculation' in line for line in lines)
        has_dynamics = any(
            'Entering Dynamics' in line or 'Molecular Dynamics Calculation' in line
            for line in lines
            )
        if has_dynamics:
            raise RuntimeError('PWSCF molecular-dynamics calculations are not supported')
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
            if not has_reference:
                raise RuntimeError('PWSCF bands calculations are not supported')
            calculation = 'nscf'
        else:
            calculation = 'scf'
        self.calculation = calculation
    #end def read_calculation



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

        ``E`` is the final completed total energy marked with ``!`` in the
        output, in Ry.
        """
        energy = None
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                value = match_float(total_energy_pattern,line)
                if value is not None:
                    energy = value
        if energy is not None:
            self.E = energy
    #end def read_energies



    def read_bands(self,lines):
        """Read and bind band data for each reported k-point.

        ``bands`` is an ``obj`` with ``up`` and ``down`` members.  Each member
        maps a zero-based k-point index to an ``obj`` containing ``eigs`` and
        ``occs``.  Eigenvalues and occupations are one-dimensional NumPy
        arrays in eV and electrons, respectively.

        Non-spin-polarized output places all records in ``bands.up``.
        Spin-polarized output separates records into ``up`` and ``down``.
        Occupation arrays can be empty.  When complete occupations are
        present, :meth:`read_band_edges` adds VBM and CBM energies.
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

        bands        = obj(up=obj(),down=obj())
        band_channel = bands.up
        for i,line in enumerate(lines):
            if ('End of self-consistent calculation' in line
                and len(bands.up)+len(bands.down)>0
                ):
                bands        = obj(up=obj(),down=obj())
                band_channel = bands.up
                continue
            if '- SPIN UP -' in line:
                band_channel = bands.up
                continue
            if '- SPIN DOWN -' in line:
                band_channel = bands.down
                continue
            if 'bands (ev)' not in line:
                continue

            eigs,j = read_values(i+1,('occupation numbers','bands (ev)'))
            if len(eigs)==0:
                continue
            while j<len(lines) and len(lines[j].strip())==0:
                j+=1
            occs = []
            if j<len(lines) and 'occupation numbers' in lines[j]:
                occs,_ = read_values(j+1)

            index = len(band_channel)
            band_channel[index] = obj(
                eigs = np.array(eigs,dtype=float),
                occs = np.array(occs,dtype=float),
                )
        if len(bands.up)+len(bands.down)==0:
            return
        self.bands = bands
        self.read_band_edges()
    #end def read_bands


    def read_band_edges(self):
        """Add VBM and CBM energies to a parsed bands object."""
        bands = self.bands
        vbm   = None
        cbm   = None
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
                if vbm is None or e_val>vbm:
                    vbm = e_val
                if cbm is None or e_cond<cbm:
                    cbm = e_cond
        if vbm is None:
            return
        bands.vbm = obj(energy=vbm)
        bands.cbm = obj(energy=cbm)
    #end def read_band_edges


    def read_structures(self,lines):
        """Read and bind structures from ionic-step output blocks.

        ``relax_structures`` is an ``obj`` mapping zero-based step indices to
        configuration objects.  Each configuration contains ``atoms`` as a
        list of element labels and ``positions`` as an ``(natoms, 3)`` NumPy
        array.  An ``axes`` ``(3, 3)`` array is included when a preceding
        ``CELL_PARAMETERS`` block is available.  Crystal positions are
        converted to Cartesian coordinates when those axes are known.

        Fixed-cell output can omit cell blocks, while variable-cell output
        normally supplies new axes with each structure.
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


    def read_pressure(self,lines):
        """Read the final reported pressure."""
        pressure = None
        for line in lines:
            if 'total' in line and 'stress' in line and 'P=' in line:
                value = match_float(pressure_pattern,line)
                if value is not None:
                    pressure = value
        if pressure is not None and 'pressure' in self:
            self.pressure = pressure
    #end def read_pressure


    def read_stress(self,lines):
        """Read and bind the sequence of reported stress tensors.

        ``stress`` is a NumPy array containing complete stress tensors in
        kbar, with shape ``(nsteps, 3, 3)``.
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
                        names  = ('kxx','kxy','kxz')
                        values = [
                            float(match.group(name).replace('D','E').replace('d','e'))
                            for name in names
                            ]
                        rows.append(values)
                if len(rows)==3:
                    stress.append(rows)
        if len(stress)>0:
            self.stress = np.array(stress,dtype=float)
    #end def read_stress


    def read_forces(self,lines):
        """Read and bind atomic-force histories.

        ``forces`` is a NumPy array with shape ``(nsteps, natoms, 3)`` in
        Ry/bohr. Atomic-force blocks with a known atom count are retained only
        when all atoms are present.
        """
        forces = []
        nat    = None
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
        if len(forces)>0:
            self.forces = np.array(forces,dtype=float)
    #end def read_forces



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


#end class PwscfOutData







class PwscfAnalyzer(SimulationAnalyzer):
    """Analyze output produced by Quantum ESPRESSO PWscf calculations.

    The analyzer coordinates PWSCF text and legacy XML readers
    for SCF, NSCF, relaxation, and variable-cell relaxation calculations.

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
    analyze : bool, optional
        If ``True``, parse the available log and legacy XML during
        initialization.

    Attributes
    ----------
    path : str
        Directory containing the PWSCF input and output files.
    abspath : str
        Absolute path to the calculation directory.
    infile_name, outfile_name : str or None
        Names of the PWSCF input and text-output files.
    input : PwscfInput or None
        Parsed PWSCF input when an input file is available.
    simulation_structure : Structure
        Input structure supplied by a ``Simulation``. This member is bound
        only when the analyzer is constructed from a simulation.
    results_out : PwscfOutData or None
        Parsed PWSCF text-output data. It is ``None`` until analysis is
        performed.
    results_xml : obj or None
        Parsed legacy XML data. It remains ``None`` when legacy XML output is
        absent or cannot be read.

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
        If a supplied path, input file, or output file does not exist.
    RuntimeError
        If a supplied file cannot be identified as input or output.

    Notes
    -----
    Log output is parsed automatically, and legacy XML is retained when
    present. Physical quantity methods use only the text-output results. A
    query returns ``None`` when its quantity applies to the detected
    calculation but was not parsed. Calling a query before analysis or for an
    unsupported calculation raises ``RuntimeError``. Supplying unsupported
    units raises ``ValueError``.

    Initial structure, k-point, and electronic quantities apply to all
    supported calculation modes. Relaxed structures apply only to ``relax``
    and ``vc-relax``; forces, stress, and pressure apply to ``scf``,
    ``relax``, and ``vc-relax``.
    """

    all_modes        = frozenset({'scf','nscf','relax','vc-relax'})
    energy_modes     = all_modes
    electronic_modes = all_modes
    relaxation_modes = frozenset({'relax','vc-relax'})
    force_modes      = frozenset({'scf','relax','vc-relax'})
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
        if calculation not in modes:
            msg = f'PWSCF quantity "{quantity}" is not supported for calculation "{calculation}"'
            raise RuntimeError(msg)
    #end def _require_supported


    def initial_structure(self,units='A'):
        """Return the initial ``Structure`` in Angstrom or bohr."""
        self._require_supported('initial_structure',self.all_modes)
        if units not in {'A','B'}:
            raise ValueError('initial_structure units must be one of: A, B')
        if 'simulation_structure' in self and self.simulation_structure is not None:
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
        if 'E' in self.results_out and self.results_out.E is not None:
            return convert(self.results_out.E,'Ry',units)
        return None
    #end def energy


    def kpoints(self,units='B'):
        """Return Cartesian k-points in inverse Angstrom or inverse bohr."""
        self._require_supported('kpoints',self.all_modes)
        if units not in {'A','B'}:
            raise ValueError('kpoints units must be one of: A, B')
        kpoints = None
        if (self.results_out.kpoints_cart is not None
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
        values = self._log_band_values('eigs')
        if values is None:
            return None
        return convert(values,'eV',units)
    #end def eigenvalues


    def occupations(self):
        """Return the dimensionless k-point-major occupation array."""
        self._require_supported('occupations',self.electronic_modes)
        return self._log_band_values('occs')
    #end def occupations


    def Ef(self,units='eV'):
        """Return the final Fermi energy in eV, Hartree, or Rydberg."""
        self._require_supported('Ef',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Ef units must be one of: eV, Ha, Ry')
        if self.results_out.Ef is not None:
            return convert(self.results_out.Ef,'eV',units)
        return None
    #end def Ef


    def Evbm(self,units='eV'):
        """Return the final valence-band maximum in selected energy units."""
        self._require_supported('Evbm',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Evbm units must be one of: eV, Ha, Ry')
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
        if ('relax_structures' in self.results_out
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
        values = self.results_out.forces
        if values is None:
            return None
        energy_units,length_units = units.split('/')
        factor = convert(1.0,'Ry',energy_units)/convert(1.0,'B',length_units)
        return values*factor
    #end def forces


    def stress(self,units='GPa'):
        """Return the stress-tensor history in selected pressure units."""
        self._require_supported('stress',self.force_modes)
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            raise ValueError('stress units must be one of: '+supported)
        values = None
        if self.results_out.stress is not None:
            values = np.asarray(self.results_out.stress,dtype=float)
            if values.ndim!=3 or values.shape[1:]!=(3,3):
                return None
        if values is None:
            return None
        return values*1e8/self.pressure_units[units]
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
        arg0         = None,
        infile_name  = None,
        outfile_name = None,
        *,
        analyze           = False,
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

        self.infile_name = infile_name
        self.outfile_name = outfile_name
        self.path         = path
        self.abspath      = os.path.abspath(path)
        self.input        = inp
        self.results_out  = None
        self.results_xml  = None
        if analyze:
            self.analyze()
    #end def __init__


    def analyze(self):
        """Analyze available PWSCF text and legacy XML output."""
        self.results_out = None
        self.results_xml = None
        if ('path' not in self
            or 'outfile_name' not in self
            or self.outfile_name is None
            ):
            raise RuntimeError('PWSCF output file name is not available')
        outfile = os.path.join(self.path,self.outfile_name)
        if not os.path.isfile(outfile):
            msg = f'PWSCF output file is not available\nfile not found: {outfile}'
            raise FileNotFoundError(msg)
        self.results_out = PwscfOutData(outfile)
        self.analyze_xml()
    #end def analyze


    def analyze_xml(self):
        """Locate and parse legacy PWscf XML output."""
        self.results_xml = None
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
            legacy_path  = os.path.join(savedir,'data-file.xml')
            if os.path.isfile(legacy_path):
                legacy_file = legacy_path
                legacy_dir  = savedir

        if legacy_file is None:
            legacy_candidates = sorted(set(
                glob(os.path.join(self.path,'*.save','data-file.xml'))
                + glob(os.path.join(self.path,'*','*.save','data-file.xml'))
                ))
            if len(legacy_candidates)==1:
                legacy_file = legacy_candidates[0]
                legacy_dir  = os.path.dirname(legacy_file)
        if legacy_file is None:
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
                except Exception:  # noqa: BLE001
                    self.results_xml.failed = True
                    continue
                eunits      = object_path(edata,'root','units_for_energies','units')
                eigenvalues = object_path(edata,'root','eigenvalues')
                occupations = object_path(edata,'root','occupations')
                if eunits is None or eigenvalues is None or occupations is None:
                    self.results_xml.failed = True
                    continue
                units = {'ha':'Ha','ry':'Ry','ev':'eV'}.get(eunits.lower()[:2],'Ha')
                spin  = obj(units=units,eigenvalues=eigenvalues,occupations=occupations)
                if si==1:
                    kp.up = spin
                elif si==2:
                    kp.down = spin
        self.results_xml.update(data=data,kpoints=kpoints)
    #end def analyze_legacy_xml







#end class PwscfAnalyzer
        
