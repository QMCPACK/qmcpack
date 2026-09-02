##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os
import re
from copy import deepcopy
from types import MappingProxyType

import numpy as np

from .developer import DevBase, obj
from .simulation import Simulation, SimulationAnalyzer
from .structure import generate_structure
from .unit_converter import UnitConverter, convert


class RmgOutData(DevBase):
    """Read an RMG output file and collect results appropriate to its run mode.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to the RMG log output file.

    Attributes
    ----------
    path : str
        Directory containing the output file.
    abspath : str
        Absolute path to the output directory.
    outfile_name : str
        Name of the RMG output file.
    setup_info : obj
        Run mode and the initial structure, when available.
    run_mode : str or None
        Short RMG calculation mode: ``"scf"``, ``"nscf"``, or ``"relax"``.
    geometry : obj or None
        Cartesian k-points and k-point weights.
    energy : float or numpy.floating or None
        Last total energy obtained from the eigenvalue sum.
    energy_units : str or None
        Units associated with ``energy``.
    electronic : obj or None
        Fermi energies, band edges, gaps, k-points, eigenvalues, and
        occupations when reported.
    forces : numpy.ndarray or None
        Ionic forces with shape ``(nsteps, natoms, 3)``.
    structures : obj or None
        Mapping from ionic-step index to a :class:`Structure` instance.
    stress : numpy.ndarray or None
        Stress tensors with shape ``(nsteps, 3, 3)``.
    pressure : float or numpy.floating or None
        Last hydrostatic pressure.

    Notes
    -----
    The supported modes are ``scf``, ``nscf``, and ``relax``. A
    mode-applicable member is initialized to ``None`` and remains ``None``
    when its data cannot be obtained.

    Raises
    ------
    TypeError
        If ``filepath`` is not a string or path-like object.
    FileNotFoundError
        If ``filepath`` does not exist.
    IsADirectoryError
        If ``filepath`` does not identify a regular file.
    """

    # Match a signed integer or decimal with an optional E- or D-exponent.
    # Example: -1.2345D+02
    number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


    def __init__(self,filepath):
        """Initialize the parsed data by reading an RMG output file."""
        if not isinstance(filepath,(str,os.PathLike)):
            msg = (
                'invalid type provided for filepath\n'
                'Type expected: str or os.PathLike\n'
                f'Type provided: {type(filepath).__name__}'
                )
            raise TypeError(msg)
        filepath = os.fspath(filepath)
        if not os.path.exists(filepath):
            msg = (
                'RMG log output file does not exist.\n'
                f'Path provided: {filepath}'
                )
            raise FileNotFoundError(msg)
        elif not os.path.isfile(filepath):
            msg = (
                'Path provided for RMG log output is not a file.\n'
                f'Path provided: {filepath}'
                )
            raise IsADirectoryError(msg)
        path,outfile_name  = os.path.split(filepath)
        self.path          = path
        self.abspath       = os.path.abspath(path)
        self.outfile_name  = outfile_name
        self.setup_info    = None

        with open(filepath,'r') as input_file:
            lines = input_file.read().splitlines()
        self.read_setup_info(lines)

        # modes: scf, nscf, relax
        if self.run_mode in {'scf','nscf','relax'}:
            self.geometry     = None
            self.energy       = None
            self.energy_units = None
            self.electronic   = None
            self.forces       = None
            self.structures   = None
            self.stress       = None
            self.pressure     = None

            self.read_geometry()
            self.read_energies(lines)
            self.read_ions(lines)
            self.read_stress(lines)
            self.read_electronic(lines)
    #end def __init__


    def read_setup_info(self,lines):
        """Read the run mode and initial structure from the setup report."""
        setup_info = obj()

        def split_assignment(line):
            """Split a colon- or equals-delimited output assignment."""
            indices = [i for i in (line.find(':'),line.find('=')) if i>=0]
            if len(indices)==0:
                return None,None
            index = min(indices)
            return line[:index].strip(),line[index+1:].strip()
        #end def split_assignment

        def as_float(token):
            """Convert an RMG-formatted numeric token when possible."""
            try:
                value = float(token.lower().replace('d','e'))
            except (AttributeError,ValueError):
                return None
            return value if np.isfinite(value) else None
        #end def as_float

        # Match an initial-position table heading and capture its units.
        # Example: Initial Ionic Positions And Displacements (Bohr)
        position_header = re.compile(
            r'^\s*initial\s+ionic\s+positions\s+and\s+displacements\s*'
            r'\(\s*(bohr|angstrom)\s*\)',re.IGNORECASE)

        run_mode = None
        for line in lines:
            label,value = split_assignment(' '.join(line.split()))
            if label is None or label.lower()!='calculation type':
                continue
            mode_words = tuple(
                word.strip('.,:;()[]{}')
                for word in value.lower().replace('-',' ').split()
                )
            mode_pairs = {
                mode_words[i:i+2] for i in range(len(mode_words)-1)
                }
            if ('quench','electrons') in mode_pairs:
                run_mode = 'scf'
            elif 'nscf' in mode_words:
                run_mode = 'nscf'
            elif (
                ('structure','optimization') in mode_pairs
                or ('relax','structure') in mode_pairs
                ):
                run_mode = 'relax'
            break
        self.run_mode       = run_mode
        setup_info.run_mode = run_mode

        axes      = {}
        axis_unit = None
        for line in lines:
            label,value = split_assignment(' '.join(line.split()))
            label_words = label.lower().split() if label is not None else []
            if len(label_words)!=3 or label_words[1:]!=['basis','vector']:
                continue
            axis = label_words[0]
            if axis not in {'x','y','z'}:
                continue
            tokens = value.split()
            if len(tokens)<3:
                continue
            values = [as_float(token) for token in tokens[:3]]
            if any(v is None for v in values):
                continue
            axes[axis] = values
            if len(tokens)>=4:
                axis_unit = tokens[3].strip(',;')

        position_tables = []
        i               = 0
        while i<len(lines):
            match = position_header.match(lines[i])
            if match is None:
                i += 1
                continue
            units     = 'B' if match.group(1).lower()=='bohr' else 'A'
            atoms     = []
            positions = []
            i += 1
            while i<len(lines):
                line = lines[i]
                if position_header.match(line) is not None:
                    break
                tokens = line.split()
                atom   = tokens[0] if len(tokens)>0 else ''
                valid_atom = (
                    len(atom)>0
                    and atom[0].isalpha()
                    and all(c.isalnum() or c=='_' for c in atom)
                    and atom.lower()!='species'
                    )
                values = (
                    [as_float(token) for token in tokens[1:4]]
                    if len(tokens)>=4 else []
                    )
                if valid_atom and len(values)==3 and None not in values:
                    atoms.append(atom)
                    positions.append(values)
                elif len(atoms)>0 and len(line.strip())==0:
                    break
                i += 1
            if len(atoms)>0:
                position_tables.append(obj(
                    units     = units,
                    atoms     = np.array(atoms,dtype=object),
                    positions = np.array(positions,dtype=float),
                    ))

        if set(axes)=={'x','y','z'} and len(position_tables)>0:
            ion_positions = next(
                (table for table in position_tables if table.units=='B'),
                position_tables[0],
                )
            aunits = 'B' if axis_unit in {None,'a0','B','bohr'} else 'A'
            axes_array = np.array(
                [axes[c] for c in ('x','y','z')],dtype=float)
            axes_array = convert(axes_array,aunits,'B')
            positions  = convert(ion_positions.positions,ion_positions.units,'B')
            valid      = (
                axes_array.shape==(3,3)
                and positions.ndim==2
                and positions.shape[1:]==(3,)
                and len(ion_positions.atoms)==len(positions)
                )
            if valid:
                setup_info.structure = generate_structure(
                    units = 'B',
                    axes  = axes_array,
                    elem  = ion_positions.atoms,
                    pos   = positions,
                    )

        kpoints  = []
        kweights = []
        for i,line in enumerate(lines):
            header_words = {
                word.strip('.,:;()[]{}').lower() for word in line.split()
                }
            if not {'kx','ky','kz','weight','crystal'}<=header_words:
                continue
            for row_line in lines[i+1:]:
                tokens = row_line.split()
                values = (
                    [as_float(token) for token in tokens[:4]]
                    if len(tokens)>=4 else []
                    )
                if len(values)!=4 or None in values:
                    if len(kpoints)>0:
                        break
                    continue
                kpoints.append(values[:3])
                kweights.append(values[3])
            break
        if len(kpoints)>0:
            setup_info.k_points = obj(
                kpoints_crystal = np.array(kpoints,dtype=float),
                kweights        = np.array(kweights,dtype=float),
                )
        self.setup_info = setup_info
    #end def read_setup_info

    def read_energies(self,lines):
        """Read the final total energy obtained from the eigenvalue sum."""
        label = 'final total energy from eig sum'
        for line in lines:
            text  = ' '.join(line.split())
            lower = text.lower()
            if label not in lower:
                continue
            remainder = text[lower.index(label)+len(label):].lstrip()
            if len(remainder)==0 or remainder[0] not in {':','='}:
                continue
            tokens = remainder[1:].split()
            if len(tokens)==0:
                continue
            try:
                value = float(tokens[0].lower().replace('d','e'))
            except ValueError:
                continue
            if not np.isfinite(value):
                continue
            self.energy       = value
            self.energy_units = tokens[1].strip(',;') if len(tokens)>=2 else 'Ha'
    #end def read_energies


    def read_electronic(self,lines):
        """Parse electronic quantities exposed by ``RmgAnalyzer``.

        Binds ``electronic`` to an ``obj`` containing Fermi energies, band
        edges, gaps, k-point-major eigenvalues, occupations, and k-points.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        def assigned_value(text,lower,*labels):
            """Return a numeric value following a labeled assignment."""
            for label in labels:
                index = lower.find(label)
                if index<0:
                    continue
                remainder = text[index+len(label):].lstrip()
                if len(remainder)==0 or remainder[0] not in {':','='}:
                    continue
                tokens = remainder[1:].split()
                if len(tokens)==0:
                    continue
                try:
                    value = float(
                        tokens[0].strip(',;').lower().replace('d','e'))
                except ValueError:
                    continue
                if np.isfinite(value):
                    return value
            return None
        #end def assigned_value

        def read_eigenvalue_data():
            """Return the final complete k-point eigenvalue and occupation table."""
            npat = self.number_pattern
            # Match an eigenvalue-table heading with its k-point index and coordinates.
            # Example: KOHN SHAM EIGENVALUES [eV] AT K-POINT [ 1]: 0.0 0.5 0.0
            header_pattern = re.compile(
                r'KOHN\s+SHAM\s+EIGENVALUES\s*\[\s*eV\s*\]\s+AT\s+'
                r'K-POINT\s*\[\s*(\d+)\s*\]\s*:\s*'
                r'('+npat+r')\s+('+npat+r')\s+('+npat+r')',re.IGNORECASE)
            # Match an eigenvalue-table row and capture its k-point index and values.
            # Example: [kpt 1 0 0] -6.4 [2.000] -1.9 [0.000]
            row_pattern = re.compile(
                r'^\s*\[\s*kpt\s+(\d+)\s+[-+]?\d+\s+\d+\s*\]\s*(.*)$',
                re.IGNORECASE)
            # Match one eigenvalue followed by its bracketed occupation.
            # Example: -6.4238 [2.000]
            pair_pattern = re.compile(
                r'('+npat+r')\s*\[\s*('+npat+r')\s*\]',re.IGNORECASE)

            datasets  = []
            dataset   = {}
            kpoint    = None
            spin      = 'none'
            for line in lines:
                match = header_pattern.search(line)
                if match is not None:
                    index = int(match.group(1))
                    if index in dataset:
                        datasets.append(dataset)
                        dataset = {}
                    coordinates = np.array([
                        float(match.group(i).lower().replace('d','e'))
                        for i in range(2,5)
                        ],dtype=float)
                    dataset[index] = {
                        'kpoint'   : coordinates,
                        'channels' : {},
                        }
                    kpoint = index
                    spin   = 'none'
                    continue
                if kpoint is None:
                    continue
                normalized = ' '.join(line.split()).lower()
                if 'spin up' in normalized:
                    spin = 'up'
                    continue
                elif 'spin down' in normalized:
                    spin = 'down'
                    continue
                match = row_pattern.match(line)
                if match is None or int(match.group(1))!=kpoint:
                    continue
                pairs = pair_pattern.findall(match.group(2))
                if len(pairs)==0:
                    continue
                channels = dataset[kpoint]['channels']
                eigs,occs = channels.setdefault(spin,[[],[]])
                for eigenvalue,occupation in pairs:
                    eigs.append(float(eigenvalue.lower().replace('d','e')))
                    occs.append(float(occupation.lower().replace('d','e')))
            if len(dataset)>0:
                datasets.append(dataset)

            expected_kpoints = None
            if 'k_points' in self.setup_info:
                expected_kpoints = len(self.setup_info.k_points.kpoints_crystal)
            for dataset in reversed(datasets):
                indices = sorted(dataset)
                if indices!=list(range(len(indices))):
                    continue
                if expected_kpoints is not None and len(indices)!=expected_kpoints:
                    continue
                spin_channels = set()
                for record in dataset.values():
                    spin_channels.update(record['channels'])
                if spin_channels=={'none'}:
                    spins = ['none']
                elif spin_channels=={'up','down'}:
                    spins = ['up','down']
                else:
                    continue
                channels = [
                    dataset[index]['channels'].get(spin)
                    for index in indices for spin in spins
                ]
                if any(
                    channel is None
                    or len(channel[0])==0
                    or len(channel[0])!=len(channel[1])
                    for channel in channels
                    ):
                    continue
                if len({len(channel[0]) for channel in channels})!=1:
                    continue
                kpoints     = np.array(
                    [dataset[i]['kpoint'] for i in indices],dtype=float)
                eigenvalues = np.array([
                    [dataset[i]['channels'][spin][0] for spin in spins]
                    for i in indices
                    ],dtype=float)
                occupations = np.array([
                    [dataset[i]['channels'][spin][1] for spin in spins]
                    for i in indices
                    ],dtype=float)
                if spins==['none']:
                    eigenvalues = eigenvalues[:,0,:]
                    occupations = occupations[:,0,:]
                return {
                    'kpoints_crystal' : kpoints,
                    'eigenvalues'     : eigenvalues,
                    'occupations'     : occupations,
                    }
            return None
        #end def read_eigenvalue_data

        data = obj(
            fermi_energies         = [],
            valence_band_maxima    = [],
            conduction_band_minima = [],
            band_gaps              = [],
            )

        for line in lines:
            text  = ' '.join(line.split())
            lower = text.lower()
            fermi = assigned_value(text,lower,'fermi energy')
            vbm   = assigned_value(text,lower,'valence band maximum')
            cbm   = assigned_value(
                text,
                lower,
                'conduction band minimum',
                'conduction band minumm',
                )
            gap = assigned_value(text,lower,'band gap')
            if fermi is not None:
                data.fermi_energies.append(fermi)
            elif vbm is not None and cbm is not None:
                data.valence_band_maxima.append(vbm)
                data.conduction_band_minima.append(cbm)
            elif gap is not None:
                data.band_gaps.append(gap)
        for name,values in data.items():
            data[name] = np.array(values,dtype=float)

        eigenvalue_data = read_eigenvalue_data()
        if eigenvalue_data is not None:
            data.kpoints_crystal = eigenvalue_data['kpoints_crystal']
            data.eigenvalues     = eigenvalue_data['eigenvalues']
            data.occupations     = eigenvalue_data['occupations']
            if 'structure' in self.setup_info:
                data.kpoints = np.dot(
                    data.kpoints_crystal,self.setup_info.structure.kaxes)

        nfound = sum(len(v) for v in data.values() if isinstance(v,np.ndarray))
        if nfound>0:
            self.electronic = data
    #end def read_electronic


    def read_ions(self,lines):
        """Read ionic forces and structures from ionic-step tables.

        Binds ``forces`` to a trajectory array and ``structures`` to a mapping
        from ionic-step indices to ``Structure`` instances.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        records = []
        i       = 0
        while i<len(lines):
            header_tokens = lines[i].split()
            is_header     = (
                len(header_tokens)>=3
                and header_tokens[0].upper()=='@ION'
                and header_tokens[1].lower()=='ion'
                and header_tokens[2].lower()=='species'
                )
            if not is_header:
                i += 1
                continue
            atoms     = []
            positions = []
            forces    = []
            i += 1
            while i<len(lines):
                tokens = lines[i].split()
                if len(tokens)==0 or tokens[0].upper()!='@ION':
                    break
                i += 1
                if len(tokens)<14:
                    continue
                numeric_tokens = tokens[3:14]
                try:
                    values = [
                        float(v.lower().replace('d','e'))
                        for v in numeric_tokens
                        ]
                except ValueError:
                    continue
                if not all(np.isfinite(values)):
                    continue
                atoms.append(tokens[2])
                positions.append(values[:3])
                forces.append(values[5:8])
            if len(atoms)>0:
                records.append(obj(
                    atoms     = np.array(atoms,dtype=object),
                    positions = np.array(positions,dtype=float),
                    forces    = np.array(forces,dtype=float),
                    ))
        if len(records)>0 and len({len(r.atoms) for r in records})==1:
            self.forces = np.array([r.forces for r in records],dtype=float)
            if 'structure' in self.setup_info:
                structures = obj()
                initial    = self.setup_info.structure
                for n,record in enumerate(records):
                    structures[n] = generate_structure(
                        units = 'B',
                        axes  = initial.axes,
                        elem  = record.atoms,
                        pos   = record.positions,
                        )
                self.structures = structures
    #end def read_ions


    def read_geometry(self):
        """Collect k-point data derived from the setup report.

        Binds ``geometry`` to an ``obj`` containing Cartesian k-points and
        k-point weights.
        """
        geometry = obj()
        if 'k_points' in self.setup_info:
            kpoints                  = self.setup_info.k_points
            geometry.kweights        = kpoints.kweights
            if 'structure' in self.setup_info and len(kpoints.kpoints_crystal)>0:
                geometry.kpoints_cart = np.dot(
                    kpoints.kpoints_crystal,self.setup_info.structure.kaxes)
        if len(geometry)>0:
            self.geometry = geometry
    #end def read_geometry


    def read_stress(self,lines):
        """Parse stress tensors and derive hydrostatic pressures.

        Binds ``stress`` to a NumPy array of shape ``(nsteps, 3, 3)`` and
        ``pressure`` to the final hydrostatic pressure. Values are in kbar.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        tensors = []
        for i,line in enumerate(lines):
            normalized = ' '.join(line.split()).lower()
            # Match the heading for a total stress tensor reported in kbar.
            # Example: stress total in unit of kbar
            if 'stress total' not in normalized or 'kbar' not in normalized:
                continue
            rows = []
            j    = i+1
            while j<len(lines) and len(rows)<3:
                row_line = lines[j]
                j += 1
                if len(row_line.strip())==0:
                    continue
                values = []
                for token in row_line.replace(',',' ').split():
                    try:
                        value = float(token.lower().replace('d','e'))
                    except ValueError:
                        break
                    if not np.isfinite(value):
                        break
                    values.append(value)
                if len(values)<3:
                    if len(rows)>0:
                        rows = []
                        break
                    continue
                if len(values)>=4 and int(values[0])==len(rows)+1:
                    values = values[1:]
                rows.append(values[:3])
            if len(rows)==3:
                tensors.append(rows)
        if len(tensors)>0:
            stress        = np.array(tensors,dtype=float)
            pressures     = -np.trace(stress,axis1=1,axis2=2)/3.0
            self.stress   = stress
            self.pressure = pressures[-1]
    #end def read_stress



class RmgAnalyzer(SimulationAnalyzer):
    """Analyze an RMG simulation or log output file.

    Parameters
    ----------
    arg0 : Simulation or str or pathlib.Path or None, optional
        RMG simulation to analyze or path to an RMG log output file. If
        ``None``, an unconfigured analyzer is created.
    analyze : bool, optional
        If ``True``, parse the RMG output during initialization.

    Attributes
    ----------
    path : str or None
        Directory containing the RMG log output file.
    abspath : str or None
        Absolute path to the output directory.
    outfile_name : str or None
        Name of the RMG log output file.
    info : obj
        General analyzer metadata.
    input : RmgInput or None
        Reserved analyzer input member; this reduced implementation leaves it
        as ``None``.
    run_mode : str or None
        Short RMG calculation mode determined during analysis: ``"scf"``,
        ``"nscf"``, or ``"relax"``.
    results : RmgOutData or None
        Parsed RMG output data. ``None`` until analysis is performed.

    Methods
    -------
    initial_structure(units='A') : Structure or None
        Input atomic structure in Angstrom (``'A'``) or bohr (``'B'``).
    energy(units='Ha') : float or numpy.floating or None
        Final total energy in ``'eV'``, ``'Ha'``, or ``'Ry'``.
    kpoints(units='B') : numpy.ndarray or None
        Cartesian k-points in inverse Angstrom or inverse bohr with shape
        ``(nkpoints, 3)``. The unit argument is ``'A'`` or ``'B'``.
    kweights() : numpy.ndarray or None
        Dimensionless k-point weights with shape ``(nkpoints,)``.
    eigenvalues(units='eV') : numpy.ndarray or None
        Kohn--Sham eigenvalues in ``'eV'``, ``'Ha'``, or ``'Ry'``. The
        leading dimension has length ``nkpoints``; remaining dimensions
        represent spin, when present, and bands.
    occupations() : numpy.ndarray or None
        Dimensionless Kohn--Sham occupations. The leading dimension has
        length ``nkpoints``; remaining dimensions represent spin, when
        present, and bands.
    Ef(units='eV') : float or numpy.floating or None
        Final Fermi energy in ``'eV'``, ``'Ha'``, or ``'Ry'``.
    Evbm(units='eV') : float or numpy.floating or None
        Final reported valence-band maximum in selected energy units.
    Ecbm(units='eV') : float or numpy.floating or None
        Final reported conduction-band minimum in selected energy units.
    band_gap(units='eV') : float or numpy.floating or None
        Final reported electronic band gap in selected energy units.
    fractional_occs() : bool or None
        Whether any occupation is farther than ``1e-3`` from both empty
        and full occupation.
    relaxed_structure(units='A') : Structure or None
        Final relaxed structure in Angstrom (``'A'``) or bohr (``'B'``).
    forces(units='eV/A') : numpy.ndarray or None
        Ionic-force history with shape ``(nsteps, natoms, 3)``. Available
        units are ``'eV/A'``, ``'Ry/B'``, and ``'Ha/B'``.
    stress(units='GPa') : numpy.ndarray or None
        Stress-tensor history with shape ``(nsteps, 3, 3)``. Available
        units are ``'Pa'``, ``'bar'``, ``'kbar'``, ``'Mbar'``, ``'GPa'``,
        ``'atm'``, ``'eV/A^3'``, ``'Ha/Bohr^3'``, and ``'Ry/Bohr^3'``.
    pressure(units='GPa') : float or numpy.floating or None
        Final hydrostatic pressure in the units accepted by ``stress``.

    Notes
    -----
    A physical query method returns ``None`` when its quantity is supported
    by the detected run mode but was not successfully parsed. Calling a
    query before analysis, or for a run mode that does not support the
    quantity, raises ``RuntimeError``. Supplying unsupported units raises
    ``ValueError``.

    Raises
    ------
    TypeError
        If ``arg0`` is neither a ``Simulation``, a string, a path-like object,
        nor ``None``.
    FileNotFoundError
        If a supplied output path does not exist.
    IsADirectoryError
        If a supplied output path does not identify a regular file.
    """

    all_modes        = frozenset({'scf','nscf','relax'})
    relaxation_modes = frozenset({'relax'})
    pressure_units   = MappingProxyType({
        'Pa'        : 1e8,
        'bar'       : 1e3,
        'kbar'      : 1.0,
        'Mbar'      : 1e-3,
        'GPa'       : 1e-1,
        'atm'       : 1e8/UnitConverter.atm,
        'eV/A^3'    : 1e8*UnitConverter.A**3/UnitConverter.eV,
        'Ha/Bohr^3' : 1e8*UnitConverter.B**3/UnitConverter.Ha,
        'Ry/Bohr^3' : 1e8*UnitConverter.B**3/UnitConverter.Ry,
        })


    def _require_supported(self,quantity,modes):
        """Require analyzed output and a run mode supporting the quantity."""
        if self.results is None:
            msg = (
                f'RMG quantity "{quantity}" is unavailable because output has not been analyzed'
                )
            raise RuntimeError(msg)
        if self.run_mode not in modes:
            msg = (
                f'RMG quantity "{quantity}" is not supported for run mode "{self.run_mode}"'
                )
            raise RuntimeError(msg)
    #end def _require_supported


    def initial_structure(self,units='A'):
        """Return the input ``Structure`` in Angstrom or bohr."""
        self._require_supported('initial_structure',self.all_modes)
        if units not in {'A','B'}:
            msg = 'initial_structure units must be one of: A, B'
            raise ValueError(msg)
        if 'structure' not in self.results.setup_info:
            return None
        structure = deepcopy(self.results.setup_info.structure)
        structure.change_units(units)
        return structure
    #end def initial_structure


    def energy(self,units='Ha'):
        """Return the final total energy in eV, Hartree, or Rydberg."""
        self._require_supported('energy',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            msg = 'energy units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        value        = self.results.energy
        source_units = self.results.energy_units
        if value is None or source_units is None:
            return None
        try:
            return convert(value,source_units,units)
        except (KeyError,TypeError,ValueError):
            return None
    #end def energy


    def kpoints(self,units='B'):
        """Return Cartesian k-points in inverse Angstrom or inverse bohr."""
        self._require_supported('kpoints',self.all_modes)
        if units not in {'A','B'}:
            msg = 'kpoints units must be one of: A, B'
            raise ValueError(msg)
        electronic = self.results.electronic if 'electronic' in self.results else None
        if electronic is not None and 'kpoints' in electronic:
            kpoints = electronic.kpoints
        else:
            geometry = self.results.geometry
            if geometry is None or 'kpoints_cart' not in geometry:
                return None
            kpoints = geometry.kpoints_cart
        return kpoints*convert(1.0,units,'B')
    #end def kpoints


    def kweights(self):
        """Return dimensionless k-point weights, or ``None`` if unavailable."""
        self._require_supported('kweights',self.all_modes)
        geometry = self.results.geometry
        if geometry is None or 'kweights' not in geometry:
            return None
        return geometry.kweights
    #end def kweights


    def eigenvalues(self,units='eV'):
        """Return K-point-major eigenvalues in eV, Hartree, or Rydberg."""
        self._require_supported('eigenvalues',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            msg = 'eigenvalues units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        electronic = self.results.electronic
        if electronic is None or 'eigenvalues' not in electronic:
            return None
        return convert(electronic.eigenvalues,'eV',units)
    #end def eigenvalues


    def occupations(self):
        """Return the dimensionless K-point-major occupation array."""
        self._require_supported('occupations',self.all_modes)
        electronic = self.results.electronic
        if electronic is None or 'occupations' not in electronic:
            return None
        return electronic.occupations
    #end def occupations


    def Ef(self,units='eV'):
        """Return the final Fermi energy in eV, Hartree, or Rydberg."""
        self._require_supported('Ef',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            msg = 'Ef units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        electronic = self.results.electronic
        if electronic is None or len(electronic.fermi_energies)==0:
            return None
        return convert(electronic.fermi_energies[-1],'eV',units)
    #end def Ef


    def Evbm(self,units='eV'):
        """Return the final valence-band maximum in selected energy units."""
        self._require_supported('Evbm',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            msg = 'Evbm units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        electronic = self.results.electronic
        if electronic is None or len(electronic.valence_band_maxima)==0:
            return None
        return convert(electronic.valence_band_maxima[-1],'eV',units)
    #end def Evbm


    def Ecbm(self,units='eV'):
        """Return the final conduction-band minimum in selected energy units."""
        self._require_supported('Ecbm',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            msg = 'Ecbm units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        electronic = self.results.electronic
        if electronic is None or len(electronic.conduction_band_minima)==0:
            return None
        return convert(electronic.conduction_band_minima[-1],'eV',units)
    #end def Ecbm


    def band_gap(self,units='eV'):
        """Return the final band gap in eV, Hartree, or Rydberg."""
        self._require_supported('band_gap',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            msg = 'band_gap units must be one of: eV, Ha, Ry'
            raise ValueError(msg)
        electronic = self.results.electronic
        if electronic is None or len(electronic.band_gaps)==0:
            return None
        return convert(electronic.band_gaps[-1],'eV',units)
    #end def band_gap


    def fractional_occs(self):
        """Whether any occupation differs from empty or full by over ``1e-3``."""
        self._require_supported('fractional_occs',self.all_modes)
        occupations = self.occupations()
        if occupations is None:
            return None
        tolerance       = 1e-3
        full_occupation = 2.0 if occupations.ndim==2 else 1.0
        empty           = np.isclose(occupations,0.0,rtol=0.0,atol=tolerance)
        full            = np.isclose(
            occupations,
            full_occupation,
            rtol = 0.0,
            atol = tolerance,
            )
        return bool(np.any(~(empty|full)))
    #end def fractional_occs


    def relaxed_structure(self,units='A'):
        """Return the final relaxed ``Structure`` in Angstrom or bohr."""
        self._require_supported('relaxed_structure',self.relaxation_modes)
        if units not in {'A','B'}:
            msg = 'relaxed_structure units must be one of: A, B'
            raise ValueError(msg)
        structures = self.results.structures
        if structures is None or len(structures)==0:
            return None
        structure = deepcopy(structures[max(structures.keys())])
        structure.change_units(units)
        return structure
    #end def relaxed_structure


    def forces(self,units='eV/A'):
        """Return ionic forces in ``eV/A``, ``Ry/B``, or ``Ha/B``."""
        self._require_supported('forces',self.all_modes)
        if units not in {'eV/A','Ry/B','Ha/B'}:
            msg = 'forces units must be one of: eV/A, Ry/B, Ha/B'
            raise ValueError(msg)
        forces = self.results.forces
        if forces is None:
            return None
        energy_units,length_units = units.split('/')
        factor = (
            convert(1.0,'Ha',energy_units)/convert(1.0,'B',length_units))
        return forces*factor
    #end def forces


    def stress(self,units='GPa'):
        """Return the stress-tensor history in selected pressure units."""
        self._require_supported('stress',self.all_modes)
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            msg = f'stress units must be one of: {supported}'
            raise ValueError(msg)
        stress = self.results.stress
        if stress is None:
            return None
        return stress*self.pressure_units[units]
    #end def stress


    def pressure(self,units='GPa'):
        """Return the final hydrostatic pressure in selected pressure units."""
        self._require_supported('pressure',self.all_modes)
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            msg = f'pressure units must be one of: {supported}'
            raise ValueError(msg)
        pressure = self.results.pressure
        if pressure is None:
            return None
        return pressure*self.pressure_units[units]
    #end def pressure


    def __init__(self,arg0=None,*,analyze=False):
        """Initialize analyzer state and optionally parse the RMG output."""
        self.path         = None
        self.abspath      = None
        self.outfile_name = None
        self.info         = obj()
        self.input        = None
        self.run_mode     = None
        self.results      = None

        if arg0 is None:
            return
        if isinstance(arg0,Simulation):
            path     = arg0.locdir
            filename = arg0.outfile
        else:
            if not isinstance(arg0,(str,os.PathLike)):
                msg = (
                    'invalid type provided for log_file\n'
                    'Type expected: str or os.PathLike\n'
                    f'Type provided: {type(arg0).__name__}'
                    )
                raise TypeError(msg)
            arg0 = os.fspath(arg0)
            if not os.path.exists(arg0):
                msg = (
                    'RMG log output file does not exist.\n'
                    f'Path provided: {arg0}'
                    )
                raise FileNotFoundError(msg)
            elif not os.path.isfile(arg0):
                msg = (
                    'Path provided for RMG log output is not a file.\n'
                    f'Path provided: {arg0}'
                    )
                raise IsADirectoryError(msg)
            path,filename = os.path.split(arg0)

        self.path         = path
        self.abspath      = os.path.abspath(path)
        self.outfile_name = filename

        if analyze:
            self.analyze()
    #end def __init__


    def analyze(self):
        """Parse the configured RMG output into an ``RmgOutData`` instance."""
        filepath      = os.path.join(self.path,self.outfile_name)
        results       = RmgOutData(filepath)
        self.results  = results
        self.run_mode = results.run_mode
    #end def analyze

#end class RmgAnalyzer
