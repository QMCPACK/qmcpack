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
import xml.etree.ElementTree as ET
from copy import deepcopy
from glob import glob
from types import MappingProxyType

import numpy as np

from .developer import DevBase, obj
from .pwscf_data_reader import read_qexml
from .pwscf_input import PwscfInput
from .simulation import Simulation, SimulationAnalyzer
from .structure import Structure, get_kpath
from .unit_converter import UnitConverter, convert
from .utilities import path_string


def parse_float(text):
    """Return a finite floating-point value from a complete numeric token."""
    if '_' in text:
        return None
    try:
        value = float(text.lower().replace('d','e'))
    except ValueError:
        return None
    if not np.isfinite(value):
        return None
    return value
#end def parse_float



class PwscfOutData(DevBase):
    """Read and store physical data from PWSCF text output.

    Parameters
    ----------
    filepath : str or os.PathLike
        Path to the PWSCF text-output file.
    calculation : str, optional
        Explicit calculation type. When omitted, the type is inferred only
        from output-text records. ``PwscfAnalyzer`` always uses inference so
        that input and output calculation types remain independent.

    Attributes
    ----------
    calculation : {'scf', 'nscf', 'bands', 'relax', 'vc-relax', 'md', 'vc-md'} or None
        Calculation type inferred from the output text, or ``None`` when no
        recognized run marker is present.
    run_type_detected : bool
        Whether ``calculation`` was explicitly supplied or inferred from a
        recognized output-log marker.
    Ef : float or None
        Final Fermi energy in eV.
    fermi_energies : numpy.ndarray or None
        One-dimensional history of Fermi energies in eV.
    bands : obj or None
        Spin-resolved band records.  The ``up`` and ``down`` members map
        k-point indices to objects containing eigenvalues, occupations,
        k-point coordinates, index, and polarization.  When complete
        occupations are present, ``vbm``, ``cbm``, ``direct_gap``,
        ``indirect_gap``, and ``electronic_structure`` may also be present.
    kpoints_cart, kpoints_unit : numpy.ndarray or None
        Cartesian and crystal k-point arrays with shape ``(nkpoints, 3)``.
    initial_structure_data : obj or None
        Initial QE-generated cell axes and Cartesian ion positions in bohr.
    kweights : numpy.ndarray or None
        One-dimensional k-point weight array.
    volume : float or None
        Final unit-cell volume in bohr cubed.
    cputime, walltime : float or None
        Total CPU and wall-clock time in hours.
    E : float or None
        Final total energy in Ry.
    relax_energies : numpy.ndarray or None
        Completed SCF total-energy history in Ry.
    scf_conv_energy, scf_conv_accuracy : numpy.ndarray or None
        Electronic SCF iteration energies and accuracy estimates in Ry.
    pressure : float or None
        Final pressure in kbar.
    stress : numpy.ndarray or None
        Stress-tensor history in kbar with shape ``(nsteps, 3, 3)``.
    forces : numpy.ndarray or None
        Atomic-force history with shape ``(nsteps, natoms, 3)`` in Ry/bohr.
    tot_forces, max_forces : numpy.ndarray or None
        Histories of reported total-force and maximum atomic-force magnitudes.
    relax_structures : list or None
        Structure records containing atom labels, Cartesian positions, and,
        when reported, cell axes.

    Notes
    -----
    Every result attribute is retained for every calculation type. Readers
    appropriate to a detected type are called selectively; when no type can
    be detected, every reader is attempted permissively. An attribute remains
    ``None`` when its reader is not selected or its record cannot be parsed.
    """

    def __init__(self,filepath,calculation=None,*,md_only=False):
        """Read a PWSCF log and initialize its accessible physical data."""
        self.calculation = None
        self.run_type_detected = False

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
        self.initial_structure_data = None
        # scf/relax/vc-relax
        self.E                 = None
        self.relax_energies    = None
        self.scf_conv_energy   = None
        self.scf_conv_accuracy = None
        self.pressure          = None
        self.stress            = None
        self.forces            = None
        self.tot_forces        = None
        self.max_forces        = None
        self.md_data           = None
        self.md_stats          = None
        # relax/vc-relax
        self.relax_structures  = None

        with open(filepath,'r') as fobj:
            lines = fobj.read().splitlines()
        # read the calculation type
        self.read_calculation(lines,calculation)
        if self.run_type_detected:
            if self.calculation in {'md','vc-md'}:
                self.read_md(lines)
                if md_only:
                    return
            # all calculations
            self.read_initial_structure(lines)
            self.read_fermi_energies(lines)
            self.read_kpoints(lines)
            self.read_bands(lines)
            self.read_volume(lines)
            # all but nscf/bands
            if self.calculation in {'scf','relax','vc-relax','md','vc-md'}:
                self.read_energies(lines)
                self.read_scf_convergence(lines)
                self.read_pressure(lines)
                self.read_stress(lines)
                self.read_forces(lines)
            # relaxation and molecular-dynamics calculations
            if self.calculation in {'relax','vc-relax','md','vc-md'}:
                self.read_structures(lines)
            self.read_timing(lines)
        else:
            # With no trustworthy run type, retain as much recognizable data
            # as possible instead of assuming an SCF-specific layout.
            self.read_md(lines)
            if md_only:
                return
            self.read_initial_structure(lines)
            self.read_fermi_energies(lines)
            self.read_kpoints(lines)
            self.read_bands(lines)
            self.read_volume(lines)
            self.read_energies(lines)
            self.read_scf_convergence(lines)
            self.read_pressure(lines)
            self.read_stress(lines)
            self.read_forces(lines)
            self.read_structures(lines)
            self.read_timing(lines)
    #end def __init__


    def read_calculation(self,lines,calculation=None):
        """Infer and bind the PWSCF calculation type from log records."""
        if calculation is not None:
            calculation = calculation.lower()
            if calculation not in {'scf','nscf','bands','relax','vc-relax','md','vc-md'}:
                msg = f'PWSCF calculation "{calculation}" is not supported'
                raise RuntimeError(msg)
            self.calculation = calculation
            self.run_type_detected = True
            return
        has_cell      = False
        has_bfgs      = False
        has_band_run  = False
        has_dynamics  = False
        has_scf_run   = False
        has_reference = False
        for line in lines:
            if not has_cell and line.strip().startswith('CELL_PARAMETERS'):
                has_cell = True
            if not has_bfgs and 'BFGS Geometry Optimization' in line:
                has_bfgs = True
            if not has_band_run and 'Band Structure Calculation' in line:
                has_band_run = True
            if not has_scf_run and 'Self-consistent Calculation' in line:
                has_scf_run = True
            if (
                not has_dynamics
                and (
                    'Entering Dynamics' in line
                    or 'Molecular Dynamics Calculation' in line
                    )
                ):
                has_dynamics = True
            if (
                not has_reference
                and (
                    'Fermi energ' in line
                    or 'highest occupied' in line
                    or 'occupation numbers' in line
                    )
                ):
                has_reference = True
        if has_dynamics:
            calculation = 'vc-md' if has_cell or any('Entering Dynamics;' in line for line in lines) else 'md'
        elif has_bfgs:
            calculation = 'vc-relax' if has_cell else 'relax'
        elif has_band_run:
            # QE uses the same heading for nscf and bands, but suppresses
            # electronic-reference and occupation records for bands runs.
            calculation = 'nscf' if has_reference else 'bands'
        elif has_scf_run:
            calculation = 'scf'
        else:
            calculation = None
        self.calculation = calculation
        self.run_type_detected = calculation is not None
    #end def read_calculation


    def read_md(self,lines):
        """Read complete molecular-dynamics records from text output."""
        records = []
        record = None
        for line in lines:
            text = ' '.join(line.split())
            if line.lstrip().startswith('!') and 'total energy' in line:
                tokens = line.replace('=',' = ').split()
                index = tokens.index('=') if '=' in tokens else len(tokens)
                value = parse_float(tokens[index+1]) if index+1<len(tokens) else None
                record = {} if value is None else {'total_energy':value}
            elif record is not None and 'total stress' in text and 'P=' in text:
                tokens = line.replace('=',' = ').split()
                if 'P' in tokens:
                    index = tokens.index('P')
                    value = parse_float(tokens[index+2]) if index+2<len(tokens) and tokens[index+1]=='=' else None
                    if value is not None:
                        record['pressure'] = value
            elif record is not None and 'time' in line and ('Entering Dynamics' in line or line.strip().startswith('time')):
                tokens = line.replace('=',' = ').split()
                if 'time' in tokens:
                    index = tokens.index('time')
                    value = parse_float(tokens[index+2]) if index+2<len(tokens) and tokens[index+1]=='=' else None
                    if value is not None:
                        record['time'] = value
            elif record is not None and ('kinetic energy' in line or line.strip().startswith('Ekin')):
                tokens = line.replace('=',' = ').split()
                if 'kinetic energy' in line and '=' in tokens:
                    index = tokens.index('=')
                    value = parse_float(tokens[index+1]) if index+1<len(tokens) else None
                    if value is not None:
                        record['kinetic_energy'] = value
                elif 'Ekin' in tokens and 'T' in tokens:
                    eindex = tokens.index('Ekin')
                    tindex = tokens.index('T')
                    evalue = parse_float(tokens[eindex+2]) if eindex+2<len(tokens) and tokens[eindex+1]=='=' else None
                    tvalue = parse_float(tokens[tindex+2]) if tindex+2<len(tokens) and tokens[tindex+1]=='=' else None
                    if evalue is not None and tvalue is not None:
                        record['kinetic_energy'],record['temperature'] = evalue,tvalue
            elif record is not None and line.strip().startswith('temperature'):
                tokens = line.replace('=',' = ').split()
                if '=' in tokens:
                    index = tokens.index('=')
                    value = parse_float(tokens[index+1]) if index+1<len(tokens) else None
                    if value is not None:
                        record['temperature'] = value
            if record is not None and all(name in record for name in ('total_energy','pressure','time','kinetic_energy','temperature')):
                records.append(record)
                record = None
        if records:
            self.md_data = obj({name:np.array([r[name] for r in records],dtype=float) for name in records[0]})
            self.md_data.potential_energy = self.md_data.total_energy-self.md_data.kinetic_energy
            self.md_stats = self.md_statistics()
    #end def read_md


    def md_statistics(self,equil=None):
        """Return mean and standard error for each MD history."""
        if self.md_data is None:
            return None
        stats = obj()
        for name,values in self.md_data.items():
            values = values[equil:] if equil is not None else values
            if len(values):
                stats[name] = (float(np.mean(values)),float(np.std(values,ddof=1)/np.sqrt(len(values))) if len(values)>1 else 0.0)
        return stats
    #end def md_statistics


    def read_fermi_energies(self,lines):
        """Read and bind the sequence of reported Fermi energies.

        ``fermi_energies`` is a one-dimensional NumPy array containing every
        successfully parsed value in eV, and ``Ef`` is its final value.  Both
        remain ``None`` when no Fermi-energy record is available.
        """
        fermi_energies = []
        for line in lines:
            tokens = line.replace('=',' = ').split()
            lower  = [token.lower() for token in tokens]
            if 'fermi' not in lower:
                continue
            index = lower.index('fermi')
            if index+2>=len(tokens) or lower[index+1] not in {'energy','energies'}:
                continue
            if lower[index+2] not in {'is','are','='}:
                continue
            values = []
            for token in tokens[index+3:index+5]:
                value = parse_float(token)
                if value is None:
                    break
                values.append(value)
            unit_index = index+3+len(values)
            if len(values)>0 and unit_index<len(tokens) and lower[unit_index]=='ev':
                fermi_energies.extend(values)
        if len(fermi_energies)>0:
            self.Ef             = fermi_energies[-1]
            self.fermi_energies = np.array(fermi_energies,dtype=float)
    #end def read_fermi_energies


    def read_energies(self,lines):
        """Read and bind completed SCF total energies.

        ``E`` is the final completed total energy marked with ``!`` in the
        output, in Ry.
        """
        energies = []
        for line in lines:
            tokens = line.replace('=',' = ').split()
            if (
                len(tokens)>=6
                and tokens[:4]==['!','total','energy','=']
                and tokens[5]=='Ry'
                ):
                value = parse_float(tokens[4])
                if value is not None:
                    energies.append(value)
        if len(energies)>0:
            self.E              = energies[-1]
            self.relax_energies = np.array(energies,dtype=float)
    #end def read_energies


    def read_scf_convergence(self,lines):
        """Read electronic SCF iteration energies and accuracy estimates."""
        energies   = []
        accuracies = []
        capture    = False
        for line in lines:
            tokens = line.replace('=',' = ').replace('<',' < ').split()
            if 'total' in tokens and 'energy' in tokens and '=' in tokens:
                capture = False
                if tokens[0]!='!':
                    index = tokens.index('=')
                    if index+2<len(tokens) and tokens[index+2]=='Ry':
                        value = parse_float(tokens[index+1])
                        if value is not None:
                            energies.append(value)
                            capture = True
            elif capture and 'estimated scf accuracy' in ' '.join(tokens):
                for operator in ('<','=','>'):
                    if operator in tokens:
                        index = tokens.index(operator)
                        if index+2<len(tokens) and tokens[index+2]=='Ry':
                            value = parse_float(tokens[index+1])
                            if value is not None:
                                accuracies.append(value)
                        break
                capture = False
        if len(energies)>0:
            self.scf_conv_energy = np.array(energies,dtype=float)
        if len(accuracies)>0:
            self.scf_conv_accuracy = np.array(accuracies,dtype=float)
    #end def read_scf_convergence


    def read_bands(self,lines):
        """Read and bind band data for each reported k-point.

        ``bands`` is an ``obj`` with ``up`` and ``down`` members.  Each member
        maps a zero-based k-point index to an ``obj`` containing ``eigs`` and
        ``occs``.  Eigenvalues and occupations are one-dimensional NumPy
        arrays in eV and electrons, respectively.

        Non-spin-polarized output places all records in ``bands.up``.
        Spin-polarized output separates records into ``up`` and ``down``.
        Occupation arrays can be empty.  When complete occupations are
        present, band-edge and gap metadata is added to ``bands``.
        """
        # Match a numeric prefix, including joined fixed-width negatives.
        # This cannot be parsed by whitespace tokenization alone.
        number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
        leading_number_list_pattern = (
            rf'^\s*(?P<values>{number_pattern}'
            rf'(?:(?:\s+|(?=[+-])){number_pattern})*)'
            )
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
                    match = re.match(leading_number_list_pattern,text)
                    if match is None:
                        break
                    numbers = [
                        parse_float(number)
                        for number in re.findall(
                            number_pattern,
                            match.group('values'),
                            )
                        ]
                    if any(number is None for number in numbers):
                        return [],i
                    values.extend(numbers)
                i+=1
            return values,i
        #end def read_values

        band_kpoint_pattern = (
            rf'\bk\s*=\s*(?P<kx>{number_pattern})(?:\s+|(?=[+-]))'
            rf'(?P<ky>{number_pattern})(?:\s+|(?=[+-]))'
            rf'(?P<kz>{number_pattern})\s*\('
            )
        bands        = obj(up=obj(),down=obj())
        band_channel = bands.up
        polarized    = any('- SPIN ' in line for line in lines)
        up_spin      = True
        for i,line in enumerate(lines):
            if (
                'End of self-consistent calculation' in line
                and len(bands.up)+len(bands.down)>0
                ):
                bands        = obj(up=obj(),down=obj())
                band_channel = bands.up
                up_spin      = True
                continue
            if '- SPIN UP -' in line:
                band_channel = bands.up
                up_spin      = True
                continue
            if '- SPIN DOWN -' in line:
                band_channel = bands.down
                up_spin      = False
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

            index       = len(band_channel)
            kpoint_cart = None
            match       = re.search(band_kpoint_pattern,line)
            if match is not None:
                coordinates = [
                    parse_float(match.group(name))
                    for name in ('kx','ky','kz')
                    ]
                if all(value is not None for value in coordinates):
                    kpoint_cart = np.array(coordinates,dtype=float)
            kpoint_rel = kpoint_cart
            if self.kpoints_cart is not None and index<len(self.kpoints_cart):
                kpoint_cart = self.kpoints_cart[index]
            if self.kpoints_unit is not None and index<len(self.kpoints_unit):
                kpoint_rel = self.kpoints_unit[index]
            band_channel[index] = obj(
                index           = index,
                kpoint_2pi_alat = kpoint_cart,
                kpoint_rel      = kpoint_rel,
                eigs            = np.array(eigs,dtype=float),
                occs            = np.array(occs,dtype=float),
                pol             = ('up' if up_spin else 'down') if polarized else 'none',
                )
        if len(bands.up)+len(bands.down)==0:
            return
        self.bands = bands

        def read_band_edges():
            """Add band edges, gaps, and electronic classification to bands."""
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
                        vbm = edge_data(band,e_val,int(np.max(np.where(occ)[0])))
                    if cbm is None or e_cond<cbm.energy:
                        cbm = edge_data(band,e_cond,int(np.min(np.where(unocc)[0])))
                    if direct_gap is None or e_cond-e_val<direct_gap.energy:
                        direct_gap = obj(
                            energy          = e_cond-e_val,
                            kpoint_rel      = band.kpoint_rel,
                            kpoint_2pi_alat = band.kpoint_2pi_alat,
                            index           = band.index,
                            pol             = band.pol,
                            )
            if vbm is None:
                return
            electronic_structure = 'insulating'
            if vbm.energy+.025>=cbm.energy:
                electronic_structure = 'metallic' if vbm.band_number==cbm.band_number else 'semi-metal'
            elif (
                vbm.kpoint_rel is not None
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

        def edge_data(band,energy,band_number):
            """Return identifying data for a valence or conduction band edge."""
            return obj(
                energy          = energy,
                kpoint_rel      = band.kpoint_rel,
                kpoint_2pi_alat = band.kpoint_2pi_alat,
                index           = band.index,
                pol             = band.pol,
                band_number     = band_number,
                )
        #end def edge_data

        read_band_edges()
    #end def read_bands


    def read_initial_structure(self,lines):
        """Read QE-generated initial axes and ion positions in bohr.

        This is the authoritative geometry fallback for Bravais lattices
        generated by QE from a nonzero ``ibrav`` value.
        """
        number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'

        def numbers(text):
            return [
                parse_float(value)
                for value in re.findall(number_pattern,text.lower().replace('d','e'))
                ]
        #end def numbers

        alat = None
        for line in lines:
            if 'lattice parameter (alat)' not in line:
                continue
            _,separator,text = line.partition('=')
            values = numbers(text) if separator else []
            if len(values)>0 and values[0] is not None:
                alat = values[0]
                break
        if alat is None:
            return

        axes = None
        for i,line in enumerate(lines):
            if 'crystal axes:' not in line.lower() or i+3>=len(lines):
                continue
            candidate = []
            for axis_line in lines[i+1:i+4]:
                _,separator,text = axis_line.partition('=')
                values = numbers(text) if separator else []
                if len(values)!=3 or any(value is None for value in values):
                    candidate = []
                    break
                candidate.append(values)
            if len(candidate)==3:
                axes = alat*np.array(candidate,dtype=float)
                break
        if axes is None:
            return

        positions = None
        atoms     = None
        for coordinate_type in ('cryst. coord.','alat units'):
            for i,line in enumerate(lines):
                if 'site n.' not in line.lower() or coordinate_type not in line.lower():
                    continue
                candidate_atoms = []
                candidate_pos   = []
                for position_line in lines[i+1:]:
                    label,separator,text = position_line.partition('=')
                    label_tokens = label.split()
                    if (
                        not separator
                        or len(label_tokens)<2
                        or not label_tokens[0].isdecimal()
                        ):
                        break
                    values = numbers(text)
                    if len(values)!=3 or any(value is None for value in values):
                        candidate_atoms = []
                        break
                    candidate_atoms.append(label_tokens[1])
                    candidate_pos.append(values)
                if len(candidate_pos)>0:
                    atoms     = candidate_atoms
                    positions = np.array(candidate_pos,dtype=float)
                    if coordinate_type=='cryst. coord.':
                        positions = np.dot(positions,axes)
                    else:
                        positions *= alat
                    break
            if positions is not None:
                break
        if positions is not None:
            self.initial_structure_data = obj(
                axes      = axes,
                atoms     = atoms,
                positions = positions,
                )
    #end def read_initial_structure


    def read_structures(self,lines):
        """Read and bind structures from ionic-step output blocks.

        ``relax_structures`` is a list of configuration objects.  Each
        configuration contains ``atoms`` as a
        list of element labels and ``positions`` as an ``(natoms, 3)`` NumPy
        array.  An ``axes`` ``(3, 3)`` array is included when a preceding
        ``CELL_PARAMETERS`` block is available.  Crystal positions are
        converted to Cartesian coordinates when those axes are known.

        Fixed-cell output can omit cell blocks, while variable-cell output
        normally supplies new axes with each structure.
        """
        def card_option(line,name):
            """Return the lower-case unit option from a QE card header."""
            text = line.strip()
            if not text.startswith(name):
                return None
            text = text[len(name):].strip().lower()
            if len(text)==0:
                return None
            if text[0] in '({':
                end = ')' if text[0]=='(' else '}'
                text = text[1:text.find(end)] if end in text else text[1:]
            return text.split()[0] if len(text)>0 else None
        #end def card_option

        def alat_from_header(line):
            """Return the bohr lattice parameter given in a card header."""
            tokens = line.replace('(',' ').replace(')',' ').replace('=',' = ').split()
            lower  = [token.lower() for token in tokens]
            if 'alat' not in lower:
                return None
            index = lower.index('alat')
            if index+2<len(tokens) and tokens[index+1]=='=':
                return parse_float(tokens[index+2])
            return None
        #end def alat_from_header

        nat = None
        for line in lines:
            if 'number of atoms/cell' not in line:
                continue
            label,separator,text = line.partition('=')
            tokens = text.split()
            if (
                separator
                and label.rstrip().endswith('number of atoms/cell')
                and len(tokens)>0
                ):
                try:
                    candidate = int(tokens[0])
                except ValueError:
                    candidate = None
                if candidate is not None and candidate>0:
                    nat = candidate
                    break

        structures = []
        conf       = None
        i          = 0
        while i<len(lines):
            line = lines[i]
            if line.strip().startswith('CELL_PARAMETERS'):
                axes = []
                if i+3<len(lines):
                    for axis_line in lines[i+1:i+4]:
                        tokens = axis_line.split()
                        values = [parse_float(token) for token in tokens[:3]]
                        if len(tokens)<3 or any(value is None for value in values):
                            axes = []
                            break
                        axes.append(values)
                if len(axes)==3:
                    conf = obj()
                    axes = np.array(axes,dtype=float)
                    option = card_option(line,'CELL_PARAMETERS')
                    alat   = alat_from_header(line)
                    if option is not None and option.startswith('ang'):
                        axes *= convert(1.0,'A','B')
                    elif option is not None and option.startswith('alat'):
                        if alat is None:
                            conf = None
                            i += 3
                            continue
                        axes *= alat
                    conf.axes = axes
                    if alat is not None:
                        conf.alat = alat
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
                    values = [parse_float(value) for value in coordinates]
                    if any(value is None for value in values):
                        break
                    atoms.append(tokens[0])
                    positions.append(values)
                    i+=1
                if len(positions)==0 or (nat is not None and len(positions)!=nat):
                    conf = None
                    if (
                        i<len(lines)
                        and (
                            lines[i].strip().startswith('CELL_PARAMETERS')
                            or 'ATOMIC_POSITIONS' in lines[i]
                            )
                        ):
                        continue
                else:
                    conf.atoms     = atoms
                    conf.positions = np.array(positions,dtype=float)
                    option = card_option(line,'ATOMIC_POSITIONS')
                    if option is not None and option.startswith('crystal') and 'axes' in conf:
                        conf.positions = np.dot(conf.positions,conf.axes)
                    elif option is not None and option.startswith('ang'):
                        conf.positions *= convert(1.0,'A','B')
                    elif option is not None and option.startswith('alat'):
                        alat = conf.alat if 'alat' in conf else None
                        if alat is not None:
                            conf.positions *= alat
                        else:
                            conf.position_units = 'alat'
                    elif option is not None and option.startswith('crystal'):
                        conf.position_units = 'crystal'
                    structures.append(conf)
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
            tokens = line.replace('=',' = ').split()
            if 'total' in tokens and 'stress' in tokens and 'P' in tokens:
                index = tokens.index('P')
                if index+2<len(tokens) and tokens[index+1]=='=':
                    value = parse_float(tokens[index+2])
                    if value is not None:
                        pressure = value
        if pressure is not None and 'pressure' in self:
            self.pressure = pressure
    #end def read_pressure


    def read_volume(self,lines):
        """Read the final reported unit-cell volume in bohr cubed."""
        volume = None
        for line in lines:
            tokens = line.replace('=',' = ').split()
            if tokens[:3]==['unit-cell','volume','='] and len(tokens)>3:
                value = parse_float(tokens[3])
                if value is not None:
                    volume = value
        if volume is not None:
            self.volume = volume
    #end def read_volume


    def read_stress(self,lines):
        """Read and bind the sequence of reported stress tensors.

        ``stress`` is a NumPy array containing complete stress tensors in
        kbar, with shape ``(nsteps, 3, 3)``.
        """
        stress = []
        for i,line in enumerate(lines):
            if 'total stress' in ' '.join(line.split()):
                rows = []
                if i+3<len(lines):
                    for stress_line in lines[i+1:i+4]:
                        tokens = stress_line.split()
                        values = [parse_float(token) for token in tokens[:6]]
                        if len(tokens)<6 or any(value is None for value in values):
                            rows = []
                            break
                        rows.append(values[3:6])
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
        nat = None
        for line in lines:
            if 'number of atoms/cell' not in line:
                continue
            label,separator,text = line.partition('=')
            tokens = text.split()
            if (
                not separator
                or not label.rstrip().endswith('number of atoms/cell')
                or len(tokens)==0
                ):
                continue
            try:
                candidate = int(tokens[0])
            except ValueError:
                continue
            if candidate>0:
                nat = candidate
                break
        force_blocks = []
        tot_forces = []
        for i,line in enumerate(lines):
            if 'Forces acting on atoms' not in line:
                continue
            aforces = []
            j       = i+1
            while j<len(lines):
                tokens = lines[j].replace('=',' = ').split()
                values = []
                if (
                    len(tokens)>=9
                    and tokens[0]=='atom'
                    and tokens[2]=='type'
                    and tokens[4:6]==['force','=']
                    ):
                    for token in tokens[6:9]:
                        value = parse_float(token)
                        if value is None:
                            break
                        values.append(value)
                if len(values)==3:
                    aforces.append(values)
                elif len(aforces)>0:
                    break
                j+=1
            if len(aforces)>0:
                force_blocks.append(aforces)
        force_count = nat
        if force_count is None and len(force_blocks)>0:
            force_counts = [len(aforces) for aforces in force_blocks]
            force_count = max(
                set(force_counts),
                key=lambda count:(force_counts.count(count),count),
                )
        forces = [
            aforces for aforces in force_blocks
            if len(aforces)==force_count
            ]
        for line in lines:
            tokens = line.replace('=',' = ').split()
            if tokens[:3]==['Total','force','='] and len(tokens)>3:
                value = parse_float(tokens[3])
                if value is not None:
                    tot_forces.append(value)
        if len(forces)>0:
            self.forces     = np.array(forces,dtype=float)
            self.max_forces = np.linalg.norm(self.forces,axis=2).max(axis=1)
        if len(tot_forces)>0:
            self.tot_forces = np.array(tot_forces,dtype=float)
    #end def read_forces


    def read_timing(self,lines):
        """Read total PWSCF CPU and wall-clock time in hours."""
        number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
        timing_value_pattern = (
            rf'(?P<value>{number_pattern})\s*(?P<unit>[hms])(?=\s|[-+.\d]|$)'
            )
        def pwscf_time(text):
            scales = {'h':1.0,'m':60.0,'s':3600.0}
            values = []
            for match in re.finditer(timing_value_pattern,text):
                value = parse_float(match.group('value'))
                if value is not None:
                    values.append(value/scales[match.group('unit')])
            return sum(values) if len(values)>0 else None
        #end def pwscf_time

        for line in lines:
            if 'PWSCF' not in line or 'CPU' not in line or 'WALL' not in line:
                continue
            label,separator,text = line.partition(':')
            if not separator or label.strip()!='PWSCF':
                continue
            cpu,separator,wall = text.partition('CPU')
            if not separator:
                continue
            wall,separator,_ = wall.partition('WALL')
            if not separator:
                continue
            cputime  = pwscf_time(cpu)
            walltime = pwscf_time(wall)
            if cputime is not None:
                self.cputime = cputime
            if walltime is not None:
                self.walltime = walltime
            if cputime is not None or walltime is not None:
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
        def read_kpoint(line):
            """Parse a complete QE k-point table row."""
            tokens = line.translate(str.maketrans('(),=','    ')).split()
            if (
                len(tokens)<7
                or tokens[0]!='k'
                or not tokens[1].isdecimal()
                or tokens[5]!='wk'
                ):
                return None
            values = [parse_float(token) for token in tokens[2:5]]
            weight = parse_float(tokens[6])
            if any(value is None for value in values) or weight is None:
                return None
            return values,weight
        #end def read_kpoint
        for i,line in enumerate(lines):
            if 'number of k points' not in line:
                continue
            label,separator,text = line.partition('=')
            tokens = text.split()
            if (
                not separator
                or not label.rstrip().endswith('number of k points')
                or len(tokens)==0
                or not tokens[0].isdecimal()
                ):
                continue
            nkpoints = int(tokens[0])
            cart_header = None
            for j in range(i+1,min(i+9,len(lines))):
                if 'cart. coord.' in lines[j]:
                    cart_header = j
                    break
            if cart_header is None:
                continue
            cart    = []
            weights = []
            cart_start = cart_header+1
            valid   = len(lines[cart_start:cart_start+nkpoints])==nkpoints
            for kline in lines[cart_start:cart_start+nkpoints]:
                kpoint = read_kpoint(kline)
                if kpoint is None:
                    valid = False
                    break
                coordinates,weight = kpoint
                cart.append(coordinates)
                weights.append(weight)
            if not valid:
                continue
            j = cart_start+nkpoints
            cryst_end = min(j+9,len(lines))
            while j<cryst_end and 'cryst. coord.' not in lines[j]:
                j+=1
            if j>=cryst_end:
                continue
            unit  = []
            valid = len(lines[j+1:j+1+nkpoints])==nkpoints
            for kline in lines[j+1:j+1+nkpoints]:
                kpoint = read_kpoint(kline)
                if kpoint is None:
                    valid = False
                    break
                coordinates,_ = kpoint
                unit.append(coordinates)
            if not valid:
                continue
            self.kpoints_cart = np.array(cart,dtype=float)
            self.kpoints_unit = np.array(unit,dtype=float)
            self.kweights     = np.array(weights,dtype=float)
            return
    #end def read_kpoints


#end class PwscfOutData



class PwscfXmlData(DevBase):
    """Read primary physical results from QE schema XML output."""

    def __init__(self,filepath):
        self.data = None
        self.parse_failed = False
        for name in ('version','calculation','total_energy','initial_atoms',
                     'initial_positions','initial_axes','initial_alat','atoms',
                     'positions','axes','alat','volume','forces','stress',
                     'spin_polarized','kpoints_rel','kweights','eigenvalues',
                     'occupations','fermi_energy'):
            self[name] = None
        try:
            root = ET.parse(filepath).getroot()
        except (OSError,LookupError,ET.ParseError):
            self.parse_failed = True
            return
        self.data = obj(root=root)
        self.extract_results(root)
    #end def __init__


    def extract_results(self,root):
        """Extract structures, energies, forces, and electronic arrays."""
        # Navigate XML elements and normalize their textual data.
        def tag(element):
            """Return an XML tag without its optional namespace."""
            return element.tag.rsplit('}',1)[-1]
        #end def tag

        def child(element,name):
            """Return a named child element, or ``None`` when absent."""
            if element is None:
                return None
            return next(
                (item for item in element if tag(item)==name),
                None,
                )
        #end def child

        def children(element,name):
            """Return all child elements with a given name."""
            if element is None:
                return []
            return [item for item in element if tag(item)==name]
        #end def children

        def element_path(element,*names):
            """Follow a sequence of named children from an XML element."""
            for name in names:
                element = child(element,name)
            return element
        #end def element_path

        def scalar(element,*,allow_text=False):
            """Return scalar text as a bool, float, array, or optional string."""
            text = '' if element is None else (element.text or '').strip()
            if len(text)==0:
                return None
            if text.lower() in {'true','false'}:
                return text.lower()=='true'
            values = [parse_float(value) for value in text.split()]
            if any(value is None for value in values):
                return text if allow_text else None
            if len(values)==1:
                return values[0]
            return np.array(values)
        #end def scalar

        def number(element):
            """Return a numeric XML scalar, excluding booleans and arrays."""
            value = scalar(element)
            return value if isinstance(value,(float,np.floating)) else None
        #end def number

        def vector(element):
            """Return numeric element text as a one-dimensional float array."""
            value = scalar(element)
            if value is None or isinstance(value,(str,bool)):
                return None
            return np.asarray(value,dtype=float).reshape(-1)
        #end def vector

        def structure(element):
            """Return atoms, positions, and axes from an atomic-structure node."""
            atom_nodes = children(child(element,'atomic_positions'),'atom')
            positions  = [vector(atom) for atom in atom_nodes]
            cell       = child(element,'cell')
            axes       = [vector(child(cell,name)) for name in ('a1','a2','a3')]
            valid_positions = (
                len(positions)>0
                and all(position is not None and len(position)>=3 for position in positions)
                )
            if valid_positions:
                atoms = np.array(
                    [atom.attrib.get('name','') for atom in atom_nodes],
                    dtype = str,
                    )
                positions = np.array([position[:3] for position in positions])
            else:
                atoms     = None
                positions = None
            valid_axes = all(axis is not None and len(axis)>=3 for axis in axes)
            axes = np.array([axis[:3] for axis in axes]) if valid_axes else None
            return atoms,positions,axes
        #end def structure

        def stack(values):
            """Stack finite arrays with a common shape, or return ``None``."""
            if len(values)==0 or any(value is None for value in values):
                return None
            arrays = [np.asarray(value,dtype=float) for value in values]
            shape  = arrays[0].shape
            if any(array.shape!=shape or not np.isfinite(array).all() for array in arrays):
                return None
            return np.stack(arrays)
        #end def stack

        # Read general metadata and the initial/final atomic structures.
        creator = element_path(root,'general_info','creator')
        if creator is not None:
            self.version = creator.attrib.get(
                'VERSION',
                creator.attrib.get('version'),
                )
        self.calculation = scalar(
            element_path(root,'input','control_variables','calculation'),
            allow_text = True,
            )
        output  = child(root,'output')
        initial = element_path(root,'input','atomic_structure')
        final   = child(output,'atomic_structure')
        iatoms,ipos,iaxes      = structure(initial)
        self.initial_atoms     = iatoms
        self.initial_positions = ipos
        self.initial_axes      = iaxes
        fatoms,fpos,faxes = structure(final)
        self.atoms        = fatoms
        self.positions    = fpos
        self.axes         = faxes
        for name,structure_data in (('initial_alat',initial),('alat',final)):
            if structure_data is not None:
                self[name] = parse_float(structure_data.attrib.get('alat',''))
        if self.axes is not None:
            self.volume = abs(np.linalg.det(self.axes))

        # Extract scalar output quantities and final force/stress tensors.
        self.total_energy = number(element_path(output,'total_energy','etot'))
        forces            = vector(child(output,'forces'))
        stress            = vector(child(output,'stress'))
        if forces is not None and len(forces)>0 and len(forces)%3==0:
            self.forces = forces.reshape(-1,3)
        if stress is not None and len(stress)==9:
            self.stress = stress.reshape(3,3)

        # Collect complete k-point records with matching band occupations.
        band           = child(output,'band_structure')
        spin_polarized = scalar(child(band,'lsda'))
        self.spin_polarized = (
            spin_polarized if isinstance(spin_polarized,bool) else None)
        self.fermi_energy = number(child(band,'fermi_energy'))
        kpoints = []
        weights = []
        eigs    = []
        occs    = []
        for record in children(band,'ks_energies'):
            kpoint      = child(record,'k_point')
            point       = vector(kpoint)
            eigenvalues = vector(child(record,'eigenvalues'))
            occupations = vector(child(record,'occupations'))
            weight = None
            if kpoint is not None:
                weight = parse_float(kpoint.attrib.get('weight',''))
            valid_record = (
                point is not None
                and len(point)>=3
                and eigenvalues is not None
                and occupations is not None
                and len(eigenvalues)==len(occupations)
                and weight is not None
                )
            if valid_record:
                kpoints.append(point[:3])
                weights.append(weight)
                eigs.append(eigenvalues)
                occs.append(occupations)
        if len(kpoints)>0:
            if self.spin_polarized:
                # Group separate spin records by their common k-point.
                groups = obj()
                for point,weight,eigenvalues,occupation in zip(kpoints,weights,eigs,occs):
                    key = tuple(np.round(point,12))
                    if key not in groups:
                        groups[key] = obj(
                            point  = point,
                            weight = weight,
                            eigs   = [],
                            occs   = [],
                            )
                    groups[key].eigs.append(eigenvalues)
                    groups[key].occs.append(occupation)
                records = list(groups.values())
                if all(len(record.eigs)==2 for record in records):
                    # Store one explicitly paired up/down record per k-point.
                    eigenvalues = stack([stack(record.eigs) for record in records])
                    occupations = stack([stack(record.occs) for record in records])
                    if eigenvalues is not None and occupations is not None:
                        self.kpoints_rel = np.array([record.point for record in records])
                        self.kweights    = np.array([record.weight for record in records])
                        self.eigenvalues = eigenvalues
                        self.occupations = occupations
                        return
                # Accommodate schema variants with adjacent spin records.
                eigenvalues = stack(eigs)
                occupations = stack(occs)
                if (
                    eigenvalues is not None
                    and occupations is not None
                    and eigenvalues.shape==occupations.shape
                    and eigenvalues.shape[1]%2==0
                    ):
                    self.kpoints_rel = np.array(kpoints)
                    self.kweights    = np.array(weights)
                    self.eigenvalues = eigenvalues.reshape(len(eigs),2,-1)
                    self.occupations = occupations.reshape(len(occs),2,-1)
                    return
            # Store a consistent non-spin band table.
            eigenvalues = stack(eigs)
            occupations = stack(occs)
            if eigenvalues is not None and occupations is not None and eigenvalues.shape==occupations.shape:
                self.kpoints_rel = np.array(kpoints)
                self.kweights    = np.array(weights)
                self.eigenvalues = eigenvalues
                self.occupations = occupations
    #end def extract_results

#end class PwscfXmlData


class Pw2CasinoAnalyzer(DevBase):
    """Read kinetic energy reported by a PW2CASINO output file.

    Parameters
    ----------
    filepath : str or os.PathLike
        Path to the PW2CASINO text-output file.

    Attributes
    ----------
    K : float or None
        Kinetic energy parsed from the last applicable ``Kinetic ... =``
        record.  It remains ``None`` when no valid record is present.
    """

    def __init__(self,filepath):
        self.K = None
        with open(filepath,'r') as fobj:
            for line in fobj:
                if 'Kinetic' not in line:
                    continue
                label,separator,text = line.partition('=')
                if not separator or not label.strip().startswith('Kinetic'):
                    continue
                tokens = text.split()
                if len(tokens)==0:
                    continue
                value = parse_float(tokens[0])
                if value is not None:
                    self.K = value
    #end def __init__

#end class Pw2CasinoAnalyzer




class PwscfAnalyzer(SimulationAnalyzer):
    """Analyze output produced by Quantum ESPRESSO PWscf calculations.

    The analyzer coordinates PWSCF text, legacy XML, and optional PW2CASINO
    readers for SCF, NSCF, relaxation, and variable-cell relaxation
    calculations.

    Parameters
    ----------
    input : Simulation, PwscfInput, str, os.PathLike, or None, optional
        PWSCF simulation, parsed input, input-file path, or directory in which
        a single ``*.in`` input can be discovered.
    outfile : str, os.PathLike, or None, optional
        Path to PWSCF text output. Relative paths are resolved below ``path``.
    analyze : bool, default=False
        If ``True``, parse the available log, legacy XML, and requested
        PW2CASINO output during initialization.
    path : str, os.PathLike, or None, optional
        Base directory for relative paths and file discovery.
    xmlfile : str, os.PathLike, or None, optional
        Explicit modern ``data-file-schema.xml`` path.
    pw2c_outfile_name : str, os.PathLike, or None, optional
        Path to optional PW2CASINO text output.
    read_all : bool, default=True
        Parse all available modern XML and text output. If ``False``, parse
        XML first and parse text output only when a constructor-required
        quantity remains unavailable.
    strict : bool, default=True
        Require explicitly supplied files, reject ambiguous discovery, and
        require at least one modern XML or text-output source. Parsed
        top-level calculation types must also agree. If ``False``, missing or
        ambiguous files are skipped and disagreements are retained as a tuple.
    required : str, iterable of str, or None, optional
        Query quantities whose absence should raise ``RuntimeError`` instead
        of returning ``None``. These quantities also control text-output
        fallback when ``read_all=False``.
    md_only : bool, default=False
        For molecular-dynamics text output, stop after parsing MD histories.

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
    input : PwscfInput or None
        Parsed PWSCF input when an input file is available.
    simulation_structure : Structure
        Input structure supplied by a ``Simulation``. This member is bound
        only when the analyzer is constructed from a simulation.
    results_out : PwscfOutData or None
        Parsed PWSCF text-output data. It is ``None`` until analysis is
        performed.
    results_xml : PwscfXmlData or obj or None
        Parsed modern schema XML, or browse-only legacy XML when schema XML is
        absent. It remains ``None`` when XML is unavailable or cannot be
        parsed.
    pw2casino : Pw2CasinoAnalyzer or None
        Parsed PW2CASINO data when an auxiliary output file is requested.
    calculation : str, tuple, or None
        Reconciled calculation type. A string is stored when all known parsed
        sources agree. Disagreement is represented by
        ``(input_type, xml_type, output_type)``, with ``None`` in unavailable
        or unparsed positions.
    required : set of str
        Validated query quantities subject to required-data behavior.
    source_status : obj
        Diagnostic input, modern-XML, and text-output states, such as
        ``'parsed'``, ``'missing'``, ``'ambiguous'``, ``'parse_failed'``, or
        ``'skipped'``.

    Methods
    -------
    initial_structure(units='A') : Structure or None
        Initial structure in Angstrom (``'A'``) or bohr (``'B'``).
    energy(units='Ha') : float or None
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
    Ef(units='eV') : float or None
        Final Fermi energy in ``'eV'``, ``'Ha'``, or ``'Ry'``.
    Evbm(units='eV') : float or None
        Final valence-band maximum in selected energy units.
    Ecbm(units='eV') : float or None
        Final conduction-band minimum in selected energy units.
    band_gap(units='eV') : float or None
        Fundamental electronic band gap in selected energy units.
    fractional_occs(tol=1e-3) : bool or None
        Whether any occupation differs from both empty and full by more than
        ``tol``.
    relaxed_structure(units='A') : Structure or None
        Final relaxed structure in Angstrom (``'A'``) or bohr (``'B'``).
    forces(units='eV/A') : numpy.ndarray or None
        Final ionic forces with shape ``(natoms, 3)``. Available units are
        ``'eV/A'``, ``'Ry/B'``, and ``'Ha/B'``.
    stress(units='GPa') : numpy.ndarray or None
        Final stress tensor with shape ``(3, 3)``. Available units are
        ``'Pa'``, ``'bar'``, ``'kbar'``, ``'Mbar'``, ``'GPa'``, and
        ``'atm'``, ``'eV/A^3'``, ``'Ha/Bohr^3'``, and ``'Ry/Bohr^3'``.
    pressure(units='GPa') : float or None
        Final hydrostatic pressure in the units accepted by ``stress``.
    require(*quantities)
        Add query quantities to the required-data policy without parsing.
    available(*quantities) : bool
        Report whether every named quantity is currently queryable.
    make_movie(filename, filepath=None)
        Write the parsed relaxation trajectory as a tiled XYZ movie.
    plot_bandstructure(...)
        Plot analyzed band energies relative to the valence-band maximum.

    Raises
    ------
    FileNotFoundError
        During strict analysis, if a supplied input, output, modern XML, or
        requested PW2CASINO file does not exist.
    RuntimeError
        During strict analysis, if discovery is ambiguous or parsed top-level
        calculation types disagree.

    Notes
    -----
    Modern schema XML is the primary source for query data, with text output
    used as a field-by-field fallback. Legacy XML is browse-only. A query
    returns ``None`` when applicable data was not parsed unless the quantity
    is required. Calling a query before completed analysis or for a
    definitively unsupported calculation raises ``RuntimeError``. Data that
    is present is returned even when unusual for the reconciled calculation
    type. Supplying unsupported units raises ``ValueError``.

    Initial structure, k-point, and electronic quantities apply to all
    supported calculation modes. Relaxed structures apply to relaxation and
    molecular-dynamics modes; forces, stress, and pressure apply to ``scf``,
    relaxation, and molecular-dynamics modes.
    """

    all_modes        = frozenset({'scf','nscf','bands','relax','vc-relax','md','vc-md'})
    energy_modes     = frozenset({'scf','nscf','relax','vc-relax','md','vc-md'})
    electronic_modes = all_modes
    relaxation_modes = frozenset({'relax','vc-relax','md','vc-md'})
    force_modes      = frozenset({'scf','relax','vc-relax','md','vc-md'})
    quantity_names   = frozenset({
        'initial_structure','energy','kpoints','kweights','eigenvalues',
        'occupations','Ef','Evbm','Ecbm','band_gap','fractional_occs',
        'relaxed_structure','forces','stress','pressure',
        })
    pressure_units   = MappingProxyType({
        'Pa'   : 1.0,
        'bar'  : 1e5,
        'kbar' : 1e8,
        'Mbar' : 1e11,
        'GPa'  : 1e9,
        'atm'  : 1.01325e5,
        'eV/A^3'    : UnitConverter.eV/UnitConverter.A**3,
        'Ha/Bohr^3' : UnitConverter.Ha/UnitConverter.B**3,
        'Ry/Bohr^3' : UnitConverter.Ry/UnitConverter.B**3,
        })


    def initial_structure(self,units='A'):
        """Return the initial structure.

        Parameters
        ----------
        units : {'A', 'B'}, default='A'
            Requested length unit: Angstrom or bohr.

        Returns
        -------
        Structure or None
            Initial structure with ``axes`` and ``pos`` arrays of shape
            ``(3, 3)`` and ``(natoms, 3)``, respectively, or ``None``.
        """
        if units not in {'A','B'}:
            msg = 'initial_structure units must be one of: A, B'
            raise ValueError(msg)
        self._require_analyzed('initial_structure')
        xml = self._schema_results()
        if xml is not None and xml.initial_atoms is not None and xml.initial_positions is not None and xml.initial_axes is not None:
            structure = Structure(
                axes    = np.asarray(xml.initial_axes,dtype=float),
                elem    = np.asarray(xml.initial_atoms,dtype=str),
                pos     = np.asarray(xml.initial_positions,dtype=float),
                units   = 'B',
                rescale = False,
                )
        elif 'simulation_structure' in self and self.simulation_structure is not None:
            structure = deepcopy(self.simulation_structure)
        elif self.results_out is not None and self.results_out.initial_structure_data is not None:
            data = self.results_out.initial_structure_data
            structure = Structure(
                axes    = np.asarray(data.axes,dtype=float),
                elem    = np.asarray(data.atoms,dtype=str),
                pos     = np.asarray(data.positions,dtype=float),
                units   = 'B',
                rescale = False,
                )
        elif (
            self.input is not None
            and 'system' in self.input
            and 'ibrav' in self.input.system
            and self.input.system.ibrav==0
            and 'atomic_positions' in self.input
            and 'cell_parameters' in self.input
            ):
            input_data = deepcopy(self.input)
            cell       = deepcopy(input_data.cell_parameters)
            positions  = np.asarray(input_data.atomic_positions.positions,dtype=float)
            cell.change_specifier('bohr',input_data)
            axes = np.asarray(cell.vectors,dtype=float)
            specifier = input_data.atomic_positions.specifier
            if specifier in {'alat',''}:
                scale = input_data.get_common_vars('scale')
                positions = positions*scale
            elif specifier=='angstrom':
                positions = positions*convert(1.0,'A','B')
            elif specifier=='crystal':
                positions = np.dot(positions,axes)
            elif specifier!='bohr':
                return self._unavailable('initial_structure',self.all_modes)
            structure = Structure(
                axes    = axes,
                elem    = np.asarray(input_data.atomic_positions.atoms,dtype=str),
                pos     = positions,
                units   = 'B',
                rescale = False,
                )
        else:
            return self._unavailable('initial_structure',self.all_modes)
        structure.change_units(units)
        return structure
    #end def initial_structure


    def energy(self,units='Ha'):
        """Return the final total energy.

        Parameters
        ----------
        units : {'eV', 'Ha', 'Ry'}, default='Ha'
            Requested energy unit.

        Returns
        -------
        float or None
            Final total energy, or ``None`` if it was not reported.
        """
        if units not in {'eV','Ha','Ry'}:
            msg = 'energy units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        self._require_analyzed('energy')
        xml = self._schema_results()
        if xml is not None and xml.total_energy is not None:
            return float(convert(xml.total_energy,'Ha',units))
        if self.results_out is not None and 'E' in self.results_out and self.results_out.E is not None:
            return float(convert(self.results_out.E,'Ry',units))
        return self._unavailable('energy',self.energy_modes)
    #end def energy


    def kpoints(self,units='B'):
        """Return Cartesian k-points.

        Parameters
        ----------
        units : {'A', 'B'}, default='B'
            Reciprocal length unit, inverse Angstrom or inverse bohr.

        Returns
        -------
        numpy.ndarray or None
            Float array with shape ``(nkpoints, 3)``, or ``None``.
        """
        if units not in {'A','B'}:
            msg = 'kpoints units must be one of: A, B'
            raise ValueError(msg)
        self._require_analyzed('kpoints')
        kpoints = None
        xml = self._schema_results()
        if xml is not None and xml.kpoints_rel is not None:
            alat = xml.alat if xml.alat is not None else xml.initial_alat
            if alat is not None and alat>0:
                kpoints = xml.kpoints_rel*2*np.pi/alat
        if kpoints is None and (
            self.results_out is not None
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
        if kpoints is None and self.results_out is not None and self.results_out.kpoints_unit is not None:
            structure = self._query_value('initial_structure','B')
            if structure is not None:
                kpoints = np.dot(self.results_out.kpoints_unit,structure.kaxes)
        if kpoints is None:
            return self._unavailable('kpoints',self.all_modes)
        return kpoints*convert(1.0,units,'B')
    #end def kpoints


    def kweights(self):
        """Return k-point integration weights.

        Returns
        -------
        numpy.ndarray or None
            One-dimensional float array of shape ``(nkpoints,)``, or ``None``.
        """
        self._require_analyzed('kweights')
        xml = self._schema_results()
        if xml is not None and xml.kweights is not None:
            return xml.kweights
        values = None if self.results_out is None else self.results_out.kweights
        if values is None:
            return self._unavailable('kweights',self.all_modes)
        return values
    #end def kweights


    def eigenvalues(self,units='eV'):
        """Return Kohn--Sham eigenvalues.

        Parameters
        ----------
        units : {'eV', 'Ha', 'Ry'}, default='eV'
            Requested energy unit.

        Returns
        -------
        numpy.ndarray or None
            Float array shaped ``(nkpoints, nbands)`` for non-spin runs or
            ``(nkpoints, 2, nbands)`` for collinear spin runs, or ``None``.
        """
        if units not in {'eV','Ha','Ry'}:
            msg = 'eigenvalues units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        self._require_analyzed('eigenvalues')
        xml = self._schema_results()
        values = None if xml is None else xml.eigenvalues
        if values is not None:
            return convert(values,'Ha',units)
        values = self._log_band_values('eigs')
        if values is None:
            return self._unavailable('eigenvalues',self.electronic_modes)
        return convert(values,'eV',units)
    #end def eigenvalues


    def occupations(self):
        """Return Kohn--Sham occupations.

        Returns
        -------
        numpy.ndarray or None
            Dimensionless float array with the shape returned by
            :meth:`eigenvalues`, or ``None``.
        """
        self._require_analyzed('occupations')
        xml = self._schema_results()
        if xml is not None and xml.occupations is not None:
            return xml.occupations
        values = self._log_band_values('occs')
        if values is None:
            return self._unavailable('occupations',self.electronic_modes)
        return values
    #end def occupations


    def Ef(self,units='eV'):
        """Return the final Fermi energy.

        Parameters
        ----------
        units : {'eV', 'Ha', 'Ry'}, default='eV'
            Requested energy unit.

        Returns
        -------
        float or None
            Fermi energy, or ``None`` when unavailable.
        """
        if units not in {'eV','Ha','Ry'}:
            msg = 'Ef units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        self._require_analyzed('Ef')
        xml = self._schema_results()
        if xml is not None and xml.fermi_energy is not None:
            return float(convert(xml.fermi_energy,'Ha',units))
        if self.results_out is not None and self.results_out.Ef is not None:
            return float(convert(self.results_out.Ef,'eV',units))
        return self._unavailable('Ef',self.electronic_modes)
    #end def Ef


    def Evbm(self,units='eV'):
        """Return the valence-band maximum.

        Parameters
        ----------
        units : {'eV', 'Ha', 'Ry'}, default='eV'
            Requested energy unit.

        Returns
        -------
        float or None
            Valence-band maximum, or ``None`` when band edges are unavailable.
        """
        if units not in {'eV','Ha','Ry'}:
            msg = 'Evbm units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        self._require_analyzed('Evbm')
        vbm,_ = self._schema_band_edges()
        if vbm is not None:
            return float(convert(vbm,'Ha',units))
        bands = None if self.results_out is None else self.results_out.bands
        if bands is not None and 'vbm' in bands:
            return float(convert(bands.vbm.energy,'eV',units))
        return self._unavailable('Evbm',self.electronic_modes)
    #end def Evbm


    def Ecbm(self,units='eV'):
        """Return the conduction-band minimum.

        Parameters
        ----------
        units : {'eV', 'Ha', 'Ry'}, default='eV'
            Requested energy unit.

        Returns
        -------
        float or None
            Conduction-band minimum, or ``None`` when unavailable.
        """
        if units not in {'eV','Ha','Ry'}:
            msg = 'Ecbm units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        self._require_analyzed('Ecbm')
        _,cbm = self._schema_band_edges()
        if cbm is not None:
            return float(convert(cbm,'Ha',units))
        bands = None if self.results_out is None else self.results_out.bands
        if bands is not None and 'cbm' in bands:
            return float(convert(bands.cbm.energy,'eV',units))
        return self._unavailable('Ecbm',self.electronic_modes)
    #end def Ecbm


    def band_gap(self,units='eV'):
        """Return the fundamental band gap.

        Parameters
        ----------
        units : {'eV', 'Ha', 'Ry'}, default='eV'
            Requested energy unit.

        Returns
        -------
        float or None
            Conduction-minus-valence energy, or ``None`` when unavailable.
        """
        if units not in {'eV','Ha','Ry'}:
            msg = 'band_gap units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        self._require_analyzed('band_gap')
        vbm,cbm = self._schema_band_edges()
        if vbm is not None and cbm is not None:
            return float(convert(cbm-vbm,'Ha',units))
        bands = None if self.results_out is None else self.results_out.bands
        if bands is not None and 'vbm' in bands and 'cbm' in bands:
            gap = bands.cbm.energy-bands.vbm.energy
            return float(convert(gap,'eV',units))
        return self._unavailable('band_gap',self.electronic_modes)
    #end def band_gap


    def fractional_occs(self,tol=1e-3):
        """Determine whether any occupation is fractional.

        Parameters
        ----------
        tol : float, default=1e-3
            Absolute tolerance for identifying empty and full occupations.

        Returns
        -------
        bool or None
            ``True`` when any occupation is fractional, or ``None`` when the
            occupation array is unavailable.
        """
        self._require_analyzed('fractional_occs')
        occupations = self._query_value('occupations')
        if occupations is None:
            return self._unavailable('fractional_occs',self.electronic_modes)
        empty = np.isclose(occupations,0.0,rtol=0.0,atol=tol)
        full  = np.isclose(occupations,1.0,rtol=0.0,atol=tol)
        return bool(np.any(~(empty|full)))
    #end def fractional_occs


    def relaxed_structure(self,units='A'):
        """Return the final ionic structure.

        Parameters
        ----------
        units : {'A', 'B'}, default='A'
            Requested length unit: Angstrom or bohr.

        Returns
        -------
        Structure or None
            Final structure with axes ``(3, 3)`` and positions
            ``(natoms, 3)``, or ``None``.
        """
        if units not in {'A','B'}:
            msg = 'relaxed_structure units must be one of: A, B'
            raise ValueError(msg)
        self._require_analyzed('relaxed_structure')
        xml = self._schema_results()
        if xml is not None and xml.atoms is not None and xml.positions is not None and xml.axes is not None:
            structure = Structure(
                axes    = np.asarray(xml.axes,dtype=float),
                elem    = np.asarray(xml.atoms,dtype=str),
                pos     = np.asarray(xml.positions,dtype=float),
                units   = 'B',
                rescale = False,
                )
        elif (
            self.results_out is not None
            and 'relax_structures' in self.results_out
            and self.results_out.relax_structures is not None
            and len(self.results_out.relax_structures)>0
            ):
            structures = self.results_out.relax_structures
            result     = structures[-1]
            initial    = self._query_value('initial_structure','B')
            axes       = result.axes if 'axes' in result else None
            if axes is None and initial is not None:
                axes = initial.axes
            if axes is None:
                return self._unavailable('relaxed_structure',self.relaxation_modes)
            positions = np.asarray(result.positions,dtype=float)
            position_units = result.position_units if 'position_units' in result else 'B'
            if position_units=='crystal':
                positions = np.dot(positions,axes)
            elif position_units=='alat':
                alat = result.alat if 'alat' in result else None
                if alat is None:
                    return self._unavailable('relaxed_structure',self.relaxation_modes)
                positions *= alat
            structure = Structure(
                axes    = np.asarray(axes,dtype=float),
                elem    = np.asarray(result.atoms,dtype=str),
                pos     = positions,
                units   = 'B',
                rescale = False,
                )
        else:
            return self._unavailable('relaxed_structure',self.relaxation_modes)
        structure.change_units(units)
        return structure
    #end def relaxed_structure


    def forces(self,units='eV/A'):
        """Return final ionic forces.

        Parameters
        ----------
        units : {'eV/A', 'Ry/B', 'Ha/B'}, default='eV/A'
            Requested energy-per-length unit.

        Returns
        -------
        numpy.ndarray or None
            Float array of shape ``(natoms, 3)`` containing final Cartesian
            forces, or ``None``.
        """
        if units not in {'eV/A','Ry/B','Ha/B'}:
            msg = 'forces units must be one of: eV/A, Ry/B, Ha/B'
            raise ValueError(msg)
        self._require_analyzed('forces')
        xml = self._schema_results()
        values = None
        source_energy = None
        if xml is not None and xml.forces is not None:
            values = xml.forces
            source_energy = 'Ha'
        if values is None:
            values = None if self.results_out is None else self.results_out.forces
            source_energy = 'Ry'
        if values is None:
            return self._unavailable('forces',self.force_modes)
        energy_units,length_units = units.split('/')
        factor                    = convert(1.0,source_energy,energy_units)/convert(1.0,'B',length_units)
        values = np.asarray(values,dtype=float)
        if values.ndim==3:
            values = values[-1]
        if values.ndim!=2 or values.shape[-1]!=3:
            return self._unavailable('forces',self.force_modes)
        return values*factor
    #end def forces


    def stress(self,units='GPa'):
        """Return the final stress tensor.

        Parameters
        ----------
        units : {'Pa', 'bar', 'kbar', 'Mbar', 'GPa', 'atm', 'eV/A^3', 'Ha/Bohr^3', 'Ry/Bohr^3'}, default='GPa'
            Requested pressure or energy-density unit.

        Returns
        -------
        numpy.ndarray or None
            Float array of shape ``(3, 3)``, or ``None``.
        """
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            msg = f'stress units must be one of: {supported}'
            raise ValueError(msg)
        self._require_analyzed('stress')
        xml = self._schema_results()
        values = None
        source_pressure = None
        if xml is not None and xml.stress is not None:
            values = np.asarray(xml.stress,dtype=float)
            source_pressure = convert(1.0,'Ha','J')/convert(1.0,'B','m')**3
        elif self.results_out is not None and self.results_out.stress is not None:
            values = np.asarray(self.results_out.stress,dtype=float)
            source_pressure = 1e8
        if values is None:
            return self._unavailable('stress',self.force_modes)
        if values.ndim==3:
            values = values[-1]
        if values.shape!=(3,3):
            return self._unavailable('stress',self.force_modes)
        return values*source_pressure/self.pressure_units[units]
    #end def stress


    def pressure(self,units='GPa'):
        """Return final hydrostatic pressure.

        Parameters
        ----------
        units : {'Pa', 'bar', 'kbar', 'Mbar', 'GPa', 'atm', 'eV/A^3', 'Ha/Bohr^3', 'Ry/Bohr^3'}, default='GPa'
            Requested pressure or energy-density unit.

        Returns
        -------
        float or None
            Trace of the final stress tensor divided by three, or ``None``.
        """
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            msg = f'pressure units must be one of: {supported}'
            raise ValueError(msg)
        self._require_analyzed('pressure')
        stress = self._query_value('stress',units)
        if stress is not None:
            return float(np.trace(stress)/3.0)
        if self.results_out is not None and self.results_out.pressure is not None:
            return float(self.results_out.pressure*1e8/self.pressure_units[units])
        return self._unavailable('pressure',self.force_modes)
    #end def pressure


    def _require_analyzed(self,quantity):
        """Require completed analysis unless a private query is in progress."""
        if self._query_depth>0:
            return
        if self.analysis_state!='analyzed':
            msg = f'PWSCF quantity "{quantity}" is unavailable because output has not been analyzed'
            raise RuntimeError(msg)
    #end def _require_analyzed


    def _unavailable(self,quantity,modes):
        """Apply unsupported and required-quantity policy to absent data."""
        if self._query_depth>0:
            return None
        calculation = self.calculation
        if isinstance(calculation,str) and calculation not in modes:
            msg = f'PWSCF quantity "{quantity}" is not supported for calculation "{calculation}"'
            raise RuntimeError(msg)
        if quantity in self.required:
            msg = f'required PWSCF quantity "{quantity}" is not available'
            raise RuntimeError(msg)
        return None
    #end def _unavailable


    def _query_value(self,quantity,*args):
        """Return a query value without applying public missing-data policy."""
        self._query_depth += 1
        try:
            return getattr(self,quantity)(*args)
        finally:
            self._query_depth -= 1
    #end def _query_value


    def _validate_quantities(self,quantities):
        """Validate quantity names atomically and return them as a tuple."""
        names = tuple(quantities)
        unknown = [
            name for name in names
            if not isinstance(name,str) or name not in self.quantity_names
            ]
        if len(unknown)>0:
            names_text = ', '.join(repr(name) for name in unknown)
            msg = f'unknown PWSCF quantity name(s): {names_text}'
            raise ValueError(msg)
        return names
    #end def _validate_quantities


    def require(self,*quantities):
        """Add query quantities to the required-data policy.

        Parameters
        ----------
        *quantities : str
            Case-sensitive names from :attr:`quantity_names`.

        Returns
        -------
        None

        Notes
        -----
        This method only updates :attr:`required`. It does not parse files,
        reanalyze existing results, or check current availability.
        """
        names = self._validate_quantities(quantities)
        self.required.update(names)
    #end def require


    def available(self,*quantities):
        """Return whether all named query quantities are available.

        Parameters
        ----------
        *quantities : str
            Case-sensitive names from :attr:`quantity_names`.

        Returns
        -------
        bool
            ``True`` only when every named query returns a value other than
            ``None``. With no names, returns ``True``.

        Raises
        ------
        RuntimeError
            If analysis has not completed.
        ValueError
            If any quantity name is invalid.

        Notes
        -----
        Required-data and calculation-applicability errors are suppressed for
        this check. The method never parses or reanalyzes data.
        """
        names = self._validate_quantities(quantities)
        if self.analysis_state!='analyzed':
            raise RuntimeError('PWSCF output has not been analyzed')
        return all(self._query_value(name) is not None for name in names)
    #end def available


    def _schema_results(self):
        """Return modern XML results, excluding browse-only legacy XML."""
        if 'results_xml' in self and isinstance(self.results_xml,PwscfXmlData):
            return self.results_xml
        return None
    #end def _schema_results


    def _log_band_values(self,name):
        """Return complete k-point-major band values from text output."""
        if self.results_out is None:
            return None
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


    def _schema_band_edges(self):
        """Return schema-XML valence and conduction edges in Hartree."""
        xml = self._schema_results()
        if xml is None or xml.eigenvalues is None or xml.occupations is None:
            return None,None
        eigenvalues = np.asarray(xml.eigenvalues,dtype=float)
        occupations = np.asarray(xml.occupations,dtype=float)
        if eigenvalues.shape!=occupations.shape:
            return None,None
        occupied = occupations>0.5
        unoccupied = occupations<0.5
        if not occupied.any() or not unoccupied.any():
            return None,None
        return np.max(eigenvalues[occupied]),np.min(eigenvalues[unoccupied])
    #end def _schema_band_edges


    def _schema_file(self):
        """Resolve the modern schema file and report discovery status."""
        if self.xmlfile is not None:
            filepath = path_string(self.xmlfile)
            if not os.path.isabs(filepath):
                filepath = os.path.join(self.path,filepath)
            status = 'found' if os.path.isfile(filepath) else 'missing'
            return filepath,status
        if self.input is not None and 'control' in self.input:
            control = self.input.control
            if 'outdir' in control and 'prefix' in control:
                savedir = f'{control.prefix}.save'
                filepath = os.path.join(
                    self.path,
                    control.outdir,
                    savedir,
                    'data-file-schema.xml',
                    )
                status = 'found' if os.path.isfile(filepath) else 'missing'
                return filepath,status
        candidates = sorted(set(
            glob(os.path.join(self.path,'*.save','data-file-schema.xml'))
            + glob(os.path.join(self.path,'*','*.save','data-file-schema.xml'))
            ))
        if len(candidates)==1:
            return candidates[0],'found'
        if len(candidates)>1:
            return candidates,'ambiguous'
        return None,'missing'
    #end def _schema_file


    def _input_file(self):
        """Resolve the PWSCF input file and report discovery status."""
        if isinstance(self.input,PwscfInput):
            return None,'parsed'
        if self.infile_name is not None:
            filepath = os.path.join(self.path,self.infile_name)
            status = 'found' if os.path.isfile(filepath) else 'missing'
            return filepath,status
        if not self.input_requested:
            return None,'omitted'
        candidates = sorted(glob(os.path.join(self.path,'*.in')))
        if len(candidates)==1:
            return candidates[0],'found'
        if len(candidates)>1:
            return candidates,'ambiguous'
        return None,'missing'
    #end def _input_file


    def _output_file(self):
        """Resolve the text-output file and report discovery status."""
        if self.outfile_name is not None:
            filepath = os.path.join(self.path,self.outfile_name)
            status = 'found' if os.path.isfile(filepath) else 'missing'
            return filepath,status
        candidates = sorted(glob(os.path.join(self.path,'*.out')))
        if self.pw2c_outfile_name is not None:
            auxiliary = os.path.abspath(
                os.path.join(self.path,self.pw2c_outfile_name),
                )
            candidates = [
                path for path in candidates
                if os.path.abspath(path)!=auxiliary
                ]
        if len(candidates)==1:
            return candidates[0],'found'
        if len(candidates)>1:
            return candidates,'ambiguous'
        return None,'missing'
    #end def _output_file


    def _set_calculation(self):
        """Reconcile calculation types from the input, XML, and output."""
        def normalized(value):
            if value is None:
                return None
            return str(value).strip().lower().replace('_','-')
        #end def normalized

        input_calculation = None
        if self.input is not None and 'control' in self.input:
            control = self.input.control
            if 'calculation' in control:
                input_calculation = control.calculation
        xml = self._schema_results()
        xml_calculation = None if xml is None else xml.calculation
        out_calculation = (
            None if self.results_out is None else self.results_out.calculation
            )
        calculations = tuple(normalized(value) for value in (
            input_calculation,
            xml_calculation,
            out_calculation,
            ))
        known = [value for value in calculations if value is not None]
        if len(known)==0:
            self.calculation = None
        elif len(set(known))==1:
            self.calculation = known[0]
        else:
            self.calculation = calculations
    #end def _set_calculation


    def __init__(
        self,
        input             = None,
        outfile           = None,
        *,
        analyze           = False,
        path              = None,
        xmlfile           = None,
        pw2c_outfile_name = None,
        read_all          = True,
        strict            = True,
        required          = None,
        md_only           = False,
        ):
        """Initialize a PWSCF output analyzer.

        Parameters
        ----------
        input : Simulation, PwscfInput, str, or os.PathLike, optional
            PWSCF simulation, parsed input, input-file path, or directory in
            which a single ``*.in`` input can be discovered.
        outfile : str or os.PathLike, optional
            Text-output path.  Relative paths are resolved below ``path``.
        analyze : bool, default=False
            Perform analysis during construction.
        path : str or os.PathLike, optional
            Base directory for relative paths and file discovery.
        xmlfile : str or os.PathLike, optional
            Explicit modern Quantum ESPRESSO schema-XML path.
        pw2c_outfile_name : str or os.PathLike, optional
            Optional PW2CASINO text-output path.
        read_all : bool, default=True
            Parse all available modern XML and text output.  If ``False``,
            parse text output only when a constructor-required quantity is
            unavailable after XML parsing.
        strict : bool, default=True
            Require explicitly supplied files, reject ambiguous discovery,
            require at least one modern XML or text-output source, and require
            parsed top-level calculation types to agree.
        required : str or iterable of str, optional
            Quantity names whose query functions must raise when unavailable.
        md_only : bool, default=False
            Restrict text parsing to molecular-dynamics data.
        """
        if not isinstance(read_all,bool):
            raise TypeError('read_all must be a bool')
        if not isinstance(strict,bool):
            raise TypeError('strict must be a bool')
        if input is not None and not isinstance(
            input,(Simulation,PwscfInput,str,os.PathLike)
            ):
            raise TypeError(
                'input must be a Simulation, PwscfInput, path, or None'
                )
        for name,value in (
            ('outfile',outfile),
            ('path',path),
            ('xmlfile',xmlfile),
            ('pw2c_outfile_name',pw2c_outfile_name),
            ):
            if value is not None and not isinstance(value,(str,os.PathLike)):
                raise TypeError(f'{name} must be a path or None')
        if required is None:
            required = ()
        elif isinstance(required,str):
            required = (required,)
        else:
            try:
                required = tuple(required)
            except TypeError as error:
                raise TypeError('required must be a quantity name or iterable') from error

        self.path              = None
        self.abspath           = None
        self.infile_name       = None
        self.outfile_name      = None
        self.pw2c_outfile_name = pw2c_outfile_name
        self.xmlfile           = xmlfile
        self.input_requested   = input is not None
        self.outfile_requested = outfile is not None
        self.xmlfile_requested = xmlfile is not None
        self.input             = None
        self.results_out       = None
        self.results_xml       = None
        self.pw2casino         = None
        self.calculation       = None
        self.read_all          = read_all
        self.strict            = strict
        self.required          = set()
        self.md_only           = md_only
        self.analysis_state    = 'not_analyzed'
        self.source_status     = obj(
            input = 'not_analyzed',
            xml   = 'not_analyzed',
            out   = 'not_analyzed',
            )
        self._query_depth      = 0
        self.require(*required)

        if isinstance(input,Simulation):
            sim                       = input
            path                      = sim.locdir if path is None else path
            input                     = sim.input
            outfile                   = sim.outfile if outfile is None else outfile
            self.simulation_structure = sim.system.structure
            self.input_requested      = True
            self.outfile_requested    = True
        if input is not None and isinstance(input,(str,os.PathLike)):
            input_path = path_string(input)
            if path is not None and not os.path.isabs(input_path):
                input_path = os.path.join(path_string(path),input_path)
            input_path = os.path.abspath(input_path)
            if os.path.isdir(input_path):
                path  = input_path
                input = None
            else:
                if path is None:
                    path = os.path.dirname(input_path)
                input = input_path
        if path is None and outfile is not None:
            path = os.path.dirname(path_string(outfile))
        if path is None and xmlfile is not None:
            path = os.path.dirname(path_string(xmlfile))
        if path is None:
            path = '.'
        self.path    = os.path.abspath(path_string(path) or '.')
        self.abspath = os.path.abspath(self.path)
        if isinstance(input, PwscfInput):
            self.input = input
        elif input is not None:
            infile = path_string(input)
            if not os.path.isabs(infile):
                infile = os.path.join(self.path,infile)
            self.infile_name = os.path.relpath(infile,self.path)
        if outfile is not None:
            outpath = path_string(outfile)
            if not os.path.isabs(outpath):
                outpath = os.path.join(self.path,outpath)
            self.outfile_name = os.path.relpath(outpath,self.path)
        if analyze:
            self.analyze()
    #end def __init__


    def analyze(self):
        """Analyze available PWSCF text, legacy XML, and PW2CASINO output."""
        self.results_out = None
        self.results_xml = None
        self.pw2casino   = None
        self.calculation = None
        self.analysis_state = 'analyzing'
        self.source_status = obj(
            input = 'parsed' if isinstance(self.input,PwscfInput) else 'omitted',
            xml   = 'missing',
            out   = 'missing',
            )
        if self.path is None:
            self.analysis_state = 'not_analyzed'
            msg = 'PWSCF output file name is not available'
            raise RuntimeError(msg)
        try:
            errors = []
            input_file,input_status = self._input_file()
            self.source_status.input = input_status
            if input_status=='found':
                self.input = PwscfInput(input_file)
                self.source_status.input = 'parsed'
            elif self.strict and self.input_requested:
                if input_status=='ambiguous':
                    paths = '\n'.join(input_file)
                    errors.append(f'multiple PWSCF input files were found\n{paths}')
                elif not isinstance(self.input,PwscfInput):
                    errors.append('PWSCF input file is not available')

            schema_file,schema_status = self._schema_file()
            output_file,output_status = self._output_file()
            self.source_status.xml = schema_status
            self.source_status.out = output_status
            if self.strict:
                if self.xmlfile_requested and schema_status=='missing':
                    errors.append('PWSCF schema XML file is not available')
                if self.outfile_requested and output_status=='missing':
                    errors.append('PWSCF output file is not available')
                if (
                    not self.xmlfile_requested
                    and not self.outfile_requested
                    and schema_status=='missing'
                    and output_status=='missing'
                    ):
                    errors.append(
                        'PWSCF schema XML file and output file are not available'
                        )
                if schema_status=='ambiguous':
                    paths = '\n'.join(schema_file)
                    errors.append(f'multiple PWSCF schema XML files were found\n{paths}')
                if output_status=='ambiguous':
                    paths = '\n'.join(output_file)
                    errors.append(f'multiple PWSCF output files were found\n{paths}')
            auxiliary = None
            if self.pw2c_outfile_name is not None:
                auxiliary = os.path.join(self.path,self.pw2c_outfile_name)
                if not os.path.isfile(auxiliary):
                    if self.strict:
                        errors.append(
                            'PW2CASINO output file is not available\n'
                            f'file not found: {auxiliary}'
                            )
                    auxiliary = None
            if len(errors)>0:
                message = '\n\n'.join(errors)
                ambiguous = any('multiple ' in error for error in errors)
                error_type = RuntimeError if ambiguous else FileNotFoundError
                raise error_type(message)

            if schema_status=='found':
                self.analyze_xml(schema_file)
                self.source_status.xml = (
                    'parsed' if self._schema_results() is not None
                    else 'parse_failed'
                    )
            else:
                # Legacy XML remains browse-only and does not change the
                # modern XML resolution status.
                self.analyze_xml(discover=False)
            self._set_calculation()

            parse_output = self.read_all
            if not self.read_all:
                parse_output = False
                if not parse_output:
                    parse_output = len(self.required)>0 and not all(
                        self._query_value(name) is not None
                        for name in self.required
                        )
            if parse_output and output_status=='found':
                self.results_out = PwscfOutData(
                    output_file,
                    md_only = self.md_only,
                    )
                self.source_status.out = 'parsed'
            elif output_status=='found':
                self.source_status.out = 'skipped'
            self._set_calculation()
            calculation_types = self.calculation
            if isinstance(calculation_types,tuple):
                input_type,xml_type,out_type = calculation_types
                known_types = {
                    calculation_type
                    for calculation_type in calculation_types
                    if calculation_type is not None
                    }
                outer_types = {'relax','vc-relax','md','vc-md'}
                if known_types & outer_types:
                    # SCF cycles are encapsulated within outer ionic runs and
                    # do not represent a conflicting calculation type.
                    known_types.discard('scf')
            else:
                known_types = set()
            if self.strict and len(known_types)>1:
                msg = (
                    'PWSCF top-level calculation types disagree.\n'
                    f'input: {input_type}\n'
                    f'xml: {xml_type}\n'
                    f'output: {out_type}'
                    )
                raise RuntimeError(msg)

            if auxiliary is not None:
                self.pw2casino = Pw2CasinoAnalyzer(auxiliary)
            self.analysis_state = 'analyzed'
        except Exception:
            self.analysis_state = 'not_analyzed'
            raise
    #end def analyze


    def analyze_xml(self,schema_file=None,*,discover=True):
        """Locate schema XML first, falling back to legacy PWscf XML."""
        self.results_xml = None

        if (
            discover
            and schema_file is None
            and self.input is not None
            and 'control' in self.input
            ):
            control = self.input.control
            if 'outdir' in control and 'prefix' in control:
                savedir = f'{control.prefix}.save'
                candidate = os.path.join(
                    self.path,
                    control.outdir,
                    savedir,
                    'data-file-schema.xml',
                    )
                if os.path.isfile(candidate):
                    schema_file = candidate
        if discover and schema_file is None:
            candidates = sorted(set(glob(os.path.join(self.path,'*.save','data-file-schema.xml')) + glob(os.path.join(self.path,'*','*.save','data-file-schema.xml'))))
            if len(candidates)==1:
                schema_file = candidates[0]
        if schema_file is not None:
            results = PwscfXmlData(schema_file)
            if not results.parse_failed:
                self.results_xml = results
                return

        legacy_file = None
        legacy_dir  = None
        if (
            'input' in self
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
        try:
            data = read_qexml(legacy_file)
            self.results_xml = obj(data=None,kpoints=None,failed=False)
            self.analyze_legacy_xml(data,legacy_dir)
        except Exception:  # noqa: BLE001
            self.results_xml = None
            return
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


    def md_statistics(self,equil=None):
        """Return summary statistics for parsed molecular-dynamics histories."""
        if self.results_out is None or 'md_data' not in self.results_out:
            return None
        return self.results_out.md_statistics(equil)
    #end def md_statistics


    def md_plots(self,*,show=True):
        """Plot molecular-dynamics energy, temperature, and pressure histories."""
        if self.results_out is None or 'md_data' not in self.results_out or self.results_out.md_data is None:
            return None
        import matplotlib.pyplot as plt
        md = self.results_out.md_data
        fig,axes = plt.subplots(3,1,sharex=True)
        axes[0].plot(md.time,md.total_energy-md.total_energy[0],label='Etot')
        axes[0].plot(md.time,md.kinetic_energy-md.kinetic_energy[0],label='Ekin')
        axes[0].plot(md.time,md.potential_energy-md.potential_energy[0],label='Epot')
        axes[0].set_ylabel('E (Ry)'); axes[0].legend()
        axes[1].plot(md.time,md.temperature); axes[1].set_ylabel('T (K)')
        axes[2].plot(md.time,md.pressure); axes[2].set_ylabel('P (kbar)'); axes[2].set_xlabel('time (ps)')
        if show: plt.show()
        return fig
    #end def md_plots


    def make_movie(self,filename,filepath=None):
        """Write the parsed relaxation trajectory as a tiled XYZ movie."""
        if 'results_out' not in self or self.results_out is None:
            msg = 'PWSCF output has not been analyzed'
            raise RuntimeError(msg)
        if (
            'relax_structures' not in self.results_out
            or self.results_out.relax_structures is None
            ):
            return
        initial = self.initial_structure('B') if self.input is not None else None
        frames  = []
        for result in self.results_out.relax_structures:
            axes = result.axes if 'axes' in result else None
            if axes is None and initial is not None:
                axes = initial.axes
            if axes is None:
                return
            positions = np.asarray(result.positions,dtype=float)
            position_units = result.position_units if 'position_units' in result else 'B'
            if position_units=='crystal':
                positions = np.dot(positions,axes)
            elif position_units=='alat':
                alat = result.alat if 'alat' in result else None
                if alat is None:
                    return
                positions *= alat
            frame = Structure(
                axes    = np.asarray(axes,dtype=float),
                elem    = np.asarray(result.atoms,dtype=str),
                pos     = positions,
                units   = 'B',
                rescale = False,
                ).tile(2,2,2)
            frames.append(frame.write_xyz())
        target_dir = self.abspath if filepath is None else filepath
        with open(os.path.join(target_dir,filename),'w') as fobj:
            fobj.write(''.join(frames))
    #end def make_movie


    def plot_bandstructure(
        self,
        filename     = None,
        filepath     = None,
        max_min_e    = None,
        *,
        show         = False,
        save         = True,
        show_vbm_cbm = True,
        k_labels     = None,
        ):
        """Plot the analyzed band structure along a reciprocal-space path."""
        import matplotlib.pyplot as plt

        if self.results_out is None or self.results_out.bands is None:
            return
        bands = self.results_out.bands
        if 'vbm' not in bands:
            return
        channels = [channel for channel in (bands.up,bands.down) if len(channel)>0]
        if len(channels)==0:
            return
        nkpoints = len(channels[0])
        if any(len(channel)!=nkpoints for channel in channels):
            return
        if k_labels is None:
            structure = self.initial_structure()
            if structure is None:
                return
            kpath  = get_kpath(structure=structure,check_standard=False)
            x      = np.asarray(kpath['explicit_path_linearcoords'],dtype=float)
            labels = list(kpath['explicit_kpoints_labels'])
            if len(x)!=nkpoints:
                return
        else:
            if self.results_out.kpoints_cart is None or len(k_labels)!=nkpoints:
                return
            labels = list(k_labels)
            kpoints = self.results_out.kpoints_cart
            x       = np.zeros(nkpoints,dtype=float)
            for index in range(1,nkpoints):
                x[index] = x[index-1]+np.linalg.norm(kpoints[index]-kpoints[index-1])
        plt.figure()
        ax = plt.gca()
        for channel,color in zip(channels,('k','r')):
            records = list(channel.values())
            nbands  = min(len(record.eigs) for record in records)
            for band_index in range(nbands):
                values = [record.eigs[band_index]-bands.vbm.energy for record in records]
                plt.plot(x,values,color=color)
        for index,label in enumerate(labels):
            if label:
                plt.axvline(x[index],linewidth=1,color='k')
                labels[index] = r'$\Gamma$' if label=='GAMMA' else f'${label}$'
        plt.xlim([np.min(x),np.max(x)])
        plt.ylim((-5,5) if max_min_e is None else max_min_e)
        plt.ylabel('Energy (eV)')
        plt.xticks(x,labels)
        ax.tick_params(axis='x',which='both',length=0,pad=10)
        if show_vbm_cbm:
            for edge,color in ((bands.vbm,'green'),(bands.cbm,'red')):
                if edge.kpoint_rel is None:
                    continue
                channel = bands.up if edge.pol!='down' else bands.down
                for index,record in channel.items():
                    if record.kpoint_rel is not None and np.equal(
                        edge.kpoint_rel,record.kpoint_rel,
                        ).all():
                        plt.scatter(x[index],edge.energy-bands.vbm.energy,c=color,s=100)
        if save:
            name   = 'band_structure.pdf' if filename is None else filename
            target = self.abspath if filepath is None else filepath
            plt.savefig(os.path.join(target,name),format='pdf',bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()
    #end def plot_bandstructure

#end class PwscfAnalyzer
