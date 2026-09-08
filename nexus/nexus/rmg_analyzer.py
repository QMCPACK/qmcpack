##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os
import re
from copy import deepcopy
from types import MappingProxyType

import numpy as np

from .developer import DevBase, dotdict, obj
from .generic import NexusError
from .rmg_input import RmgInput
from .simulation import Simulation, SimulationAnalyzer
from .structure import generate_structure
from .unit_converter import UnitConverter, convert


def as_float(text):
    """Convert a finite RMG-formatted number to a float when possible."""
    try:
        value = float(text.lower().replace('d','e'))
    except (AttributeError,ValueError):
        return None
    return value if np.isfinite(value) else None
#end def as_float


def normalize_line(line):
    """Collapse whitespace in an output line."""
    return ' '.join(line.split())
#end def normalize_line


def line_numbers(line,number_pattern):
    """Return all finite RMG-formatted numbers found in a line."""
    # Match each signed integer or decimal with an optional E- or D-exponent.
    # Example: SUM FORCE = 0.1 0.2 0.3
    values = re.findall(number_pattern,line)
    values = [as_float(value) for value in values]
    return np.array([value for value in values if value is not None],dtype=float)
#end def line_numbers


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
    input : RmgInput or None
        Parsed control input when the referenced input file is available.
    setup_info : obj
        Parsed setup sections, run mode, structure, and k-point information.
    run_mode : str or None
        Short RMG calculation mode: ``"scf"``, ``"nscf"``, or ``"relax"``.
    geometry : obj or None
        Cell volume and crystal/Cartesian k-point information.
    convergence : obj or None
        Electronic and ionic convergence indicators and event counts.
    timing : obj or None
        Total, per-step, and section-resolved timing data in seconds.
    energy : float or numpy.floating or None
        Last total energy obtained from the eigenvalue sum.
    energy_units : str or None
        Units associated with ``energy``.
    energies : numpy.ndarray or None
        History of total energies obtained from eigenvalue sums.
    energy_units_history : numpy.ndarray or None
        Units corresponding to ``energies``.
    direct_energies : numpy.ndarray or None
        History of directly evaluated total energies.
    direct_energy_units : numpy.ndarray or None
        Units corresponding to ``direct_energies``.
    electronic : obj or None
        Fermi energies, band edges, gaps, k-points, eigenvalues, and
        occupations, charge, magnetization, force, volume, and per-atom energy
        data when reported.
    scf : obj or None
        SCF energy components, iteration indices, residuals, and timing data.
    ionic_steps : obj or None
        Detailed per-step ionic records.
    position_units : str or None
        Units associated with ionic positions.
    force_units : str or None
        Units associated with ionic forces.
    positions : numpy.ndarray or None
        Ionic positions with shape ``(nsteps, natoms, 3)``.
    forces : numpy.ndarray or None
        Ionic forces with shape ``(nsteps, natoms, 3)``.
    charges : numpy.ndarray or None
        Ionic charges with shape ``(nsteps, natoms)``.
    magnetizations : numpy.ndarray or None
        Ionic magnetizations with shape ``(nsteps, natoms)``.
    max_forces : numpy.ndarray or None
        Maximum ionic force magnitude at each ionic step.
    structures : obj or None
        Mapping from ionic-step index to a :class:`Structure` instance.
    stress : numpy.ndarray or None
        Stress tensors with shape ``(nsteps, 3, 3)``.
    stress_units : str or None
        Units associated with stress and pressure values.
    pressures : numpy.ndarray or None
        Hydrostatic pressure at each reported stress step.
    pressure : float or numpy.floating or None
        Last hydrostatic pressure.
    produced_files : obj or None
        Paths to recognized files produced by an SCF run.

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
        path,outfile_name = os.path.split(filepath)
        self.path         = path
        self.abspath      = os.path.abspath(path)
        self.outfile_name = outfile_name
        self.input        = None
        self.setup_info   = None

        with open(filepath,'r') as input_file:
            lines = input_file.read().splitlines()
        self.read_setup_info(lines)

        # modes: scf, nscf, relax
        if self.run_mode in {'scf','nscf','relax'}:
            self.geometry             = None
            self.convergence          = None
            self.timing               = None
            self.energy               = None
            self.energy_units         = None
            self.energies             = None
            self.energy_units_history = None
            self.direct_energies      = None
            self.direct_energy_units  = None
            self.electronic           = None
            self.scf                  = None
            self.ionic_steps          = None
            self.position_units       = None
            self.force_units          = None
            self.positions            = None
            self.forces               = None
            self.charges              = None
            self.magnetizations       = None
            self.max_forces           = None
            self.structures           = None
            self.stress               = None
            self.stress_units         = None
            self.pressures            = None
            self.pressure             = None

            self.read_geometry()
            self.read_convergence(lines)
            self.read_timing(lines)
            self.read_energies(lines)
            self.read_scf(lines)
            self.read_ions(lines)
            self.read_stress(lines)
            self.read_electronic(lines)

        # modes: scf
        if self.run_mode=='scf':
            self.produced_files = None
            self.read_produced_files()
    #end def __init__


    def read_setup_info(self,lines):
        """Read setup sections, the run mode, and the initial structure.

        Binds ``setup_info`` to an ``obj`` containing normalized setup blocks
        and derived grid, lattice, ion, k-point, and structure data. When the
        referenced control file is available, also binds ``input``.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        setup_info = obj(
            run_mode  = None,
            structure = None,
            k_points  = None,
            files     = None,
            )
        position_heading = 'initial ionic positions and displacements'

        def process_name(text):
            """Convert an RMG setup label to a normalized member name."""
            # Remove parenthetical annotations from setup labels.
            # Example: Grid Points (Linear Anisotropy: 1.000)
            text = re.sub(r'\([^)]*\)','',text)
            name = '_'.join(text.strip().lower().split())
            return name.replace('/','_').replace('-','_')
        #end def process_name

        # Parse the indented setup report into named persistent sections.
        sections      = obj()
        current       = None
        files         = None
        grid_points   = None
        lattice_setup = None
        in_setup      = False
        section_added = False
        for raw_line in lines:
            stripped = raw_line.strip()
            if not in_setup:
                if stripped.lower()!='files':
                    continue
                in_setup = True
            if normalize_line(stripped).lower().startswith(position_heading):
                break
            if len(stripped)==0:
                continue
            if not raw_line[0].isspace():
                section_name  = process_name(stripped.rstrip(':'))
                section_added = False
                if section_name=='files':
                    current = obj(
                        control_input_file = None,
                        data_output_file   = None,
                        )
                    files = current
                elif section_name=='grid_points':
                    current = obj(
                        equivalent_energy_cutoffs = None,
                        units                     = None,
                        )
                    grid_points = current
                elif section_name=='lattice_setup':
                    current       = obj()
                    lattice_setup = current
                else:
                    current = obj()
                continue
            if current is None or ':' not in stripped:
                continue
            if not section_added and section_name!='k_points':
                sections[section_name] = current
                section_added = True
            label,value = stripped.split(':',1)
            name  = process_name(label)
            value = value.strip()
            units = None

            # Convert simple Boolean, integer, and floating-point fields.
            upper_value = value.upper()
            number      = as_float(value)
            if upper_value in {'ON','OFF'}:
                value = upper_value=='ON'
            elif number is not None:
                is_integer = (
                    number.is_integer()
                    and not any(c in value.lower() for c in '.e')
                    )
                value = int(number) if is_integer else number
            else:
                # Convert numeric sequences, separating a trailing unit label.
                tokens = value.replace(',',' ').split()
                if len(tokens)>1 and as_float(tokens[-1]) is None:
                    numeric = [as_float(token) for token in tokens[:-1]]
                    if None not in numeric:
                        units  = tokens[-1]
                        tokens = tokens[:-1]
                values = [as_float(token) for token in tokens]
                if (
                    (len(values)>1 or units is not None and len(values)>0)
                    and None not in values
                    ):
                    value = np.array(values,dtype=float)
            current[name] = value
            if units is not None:
                current.units = units
        setup_info.update(sections)

        # Read the calculation mode from the run setup block.
        run_mode = None
        for line in lines:
            label,separator,value = normalize_line(line).partition(':')
            if len(separator)==0 or label.strip().lower()!='calculation type':
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

        # Read the lattice vectors used to construct the input structure.
        axes      = {}
        axis_unit = None
        for line in lines:
            label,separator,value = normalize_line(line).partition(':')
            if len(separator)==0:
                continue
            label_words = label.lower().split()
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
        i = 0
        # Read each initial-position table, retaining its reported units.
        while i<len(lines):
            line = normalize_line(lines[i]).lower()
            if not line.startswith(position_heading) or '(' not in line:
                i += 1
                continue
            units = line.split('(')[-1].partition(')')[0].strip()
            if units not in {'bohr','angstrom'}:
                i += 1
                continue
            units     = 'B' if units=='bohr' else 'A'
            atoms     = []
            positions = []
            i += 1
            while i<len(lines):
                line = lines[i]
                if len(line.strip())==0 and len(atoms)>0:
                    break
                tokens     = line.split()
                atom       = tokens[0] if len(tokens)>0 else ''
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
                i += 1
            if len(atoms)>0:
                position_tables.append(
                    obj(
                        units     = units,
                        atoms     = np.array(atoms,dtype=object),
                        positions = np.array(positions,dtype=float),
                        ),
                    )

        if grid_points is not None:
            grid         = []
            grid_pe      = []
            grid_spacing = []
            for direction in ('x','y','z'):
                if direction not in grid_points:
                    break
                text   = normalize_line(str(grid_points[direction]))
                lower  = text.lower()
                values = []
                for label in ('total','per pe','spacing'):
                    start = lower.find(label)
                    if start<0:
                        break
                    remainder = text[start+len(label):].lstrip(' :')
                    tokens    = remainder.split()
                    value     = as_float(tokens[0]) if len(tokens)>0 else None
                    if value is None:
                        break
                    values.append(value)
                if len(values)!=3:
                    break
                grid.append(values[0])
                grid_pe.append(values[1])
                grid_spacing.append(values[2])
            if len(grid)==3:
                grid_points.grid         = np.array(grid,dtype=int)
                grid_points.grid_pe      = np.array(grid_pe,dtype=int)
                grid_points.grid_spacing = np.array(grid_spacing,dtype=float)
                grid_points.grid_units   = 'a0'
            cutoffs = grid_points.equivalent_energy_cutoffs
            if cutoffs is not None:
                cutoff_values = line_numbers(str(cutoffs),self.number_pattern)
                if len(cutoff_values)>=2:
                    grid_points.ecut        = cutoff_values[0]
                    grid_points.ecut_charge = cutoff_values[1]
                    grid_points.ecut_units  = grid_points.units

        if set(axes)=={'x','y','z'} and len(position_tables)>0:
            ion_positions = next(
                (table for table in position_tables if table.units=='B'),
                position_tables[0],
                )
            setup_info.ion_positions = ion_positions
            aunits        = 'B' if axis_unit in {None,'a0','B','bohr'} else 'A'
            reported_axes = np.array(
                [axes[c] for c in ('x','y','z')],
                dtype=float,
                )
            axes_array = convert(reported_axes,aunits,'B')
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
                if lattice_setup is not None:
                    lattice_setup.axes = reported_axes

        kpoints  = []
        kweights = []
        # Read crystal k-points and attach them to the input structure.
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
            kpoints  = np.array(kpoints,dtype=float)
            kweights = np.array(kweights,dtype=float)
            setup_info.k_points = obj(
                kpoints_crystal = kpoints,
                kweights        = kweights,
                )
            if setup_info.structure is not None:
                setup_info.structure.add_kpoints(
                    kpoints,
                    kweights,
                    recenter  = False,
                    cell_unit = True,
                    )

        if files is not None and files.control_input_file is not None:
            control_file = str(files.control_input_file)
            filepaths    = (
                os.path.join(self.path,control_file),
                os.path.join(self.path,os.path.basename(control_file)),
                os.path.join(os.path.dirname(self.path),control_file),
                )
            filepath = next(
                (path for path in filepaths if os.path.isfile(path)),None)
            if filepath is not None:
                try:
                    self.input = RmgInput(filepath)
                except (NexusError,OSError,TypeError,ValueError):
                    pass
        self.setup_info = setup_info
    #end def read_setup_info

    def read_energies(self,lines):
        """Read eigenvalue-sum and direct total-energy histories.

        Binds the history arrays and their unit arrays, as well as the final
        eigenvalue-sum energy and its units.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        energies            = []
        energy_units        = []
        direct_energies     = []
        direct_energy_units = []
        for line in lines:
            text         = normalize_line(line)
            lower        = text.lower()
            label         = None
            target_values = None
            target_units  = None
            if 'final total energy from eig sum' in lower:
                label         = 'final total energy from eig sum'
                target_values = energies
                target_units  = energy_units
            elif 'final total energy from direct' in lower:
                label         = 'final total energy from direct'
                target_values = direct_energies
                target_units  = direct_energy_units
            if label is None:
                continue
            remainder = text[lower.index(label)+len(label):].lstrip()
            if len(remainder)==0 or remainder[0] not in {':','='}:
                continue
            tokens = remainder[1:].split()
            if len(tokens)==0:
                continue
            value = as_float(tokens[0])
            if value is None:
                continue
            target_values.append(value)
            target_units.append(
                tokens[1].strip(',;') if len(tokens)>=2 else None,
                )
        if len(energies)>0:
            self.energies             = np.array(energies,dtype=float)
            self.energy_units_history = np.array(energy_units,dtype=object)
            self.energy               = self.energies[-1]
            self.energy_units         = self.energy_units_history[-1] or 'Ha'
        if len(direct_energies)>0:
            self.direct_energies     = np.array(direct_energies,dtype=float)
            self.direct_energy_units = np.array(
                direct_energy_units,
                dtype=object,
                )
    #end def read_energies


    def read_electronic(self,lines):
        """Parse electronic quantities exposed by ``RmgAnalyzer``.

        Binds ``electronic`` to an ``obj`` containing Fermi energies, band
        edges, gaps, charge and magnetization values, summed forces, per-atom
        volume and energy, k-point-major eigenvalues and occupations, and
        k-points.

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
                value = as_float(tokens[0].strip(',;'))
                if value is not None:
                    return value
            return None
        #end def assigned_value

        data = obj(
            fermi_energies          = [],
            valence_band_maxima     = [],
            conduction_band_minima  = [],
            band_gaps               = [],
            total_charges           = [],
            total_magnetizations    = [],
            absolute_magnetizations = [],
            sum_forces              = [],
            volume_per_atom         = [],
            energy_per_atom         = [],
            kpoints_crystal         = None,
            kpoints                 = None,
            eigenvalues             = None,
            occupations             = None,
            )

        # Match one eigenvalue followed by its bracketed occupation.
        # Example: -6.4238 [2.000]
        npat         = self.number_pattern
        pair_pattern = re.compile(
            r'('+npat+r')\s*\[\s*('+npat+r')\s*\]',
            re.IGNORECASE,
            )
        datasets = []
        dataset  = dotdict()
        kpoint   = None
        spin     = 'none'
        # Collect scalar results and candidate eigenvalue tables in one pass.
        for line in lines:
            text  = normalize_line(line)
            lower = text.lower()
            fermi = assigned_value(text,lower,'fermi energy')
            vbm   = assigned_value(text,lower,'valence band maximum')
            cbm   = assigned_value(
                text,
                lower,
                'conduction band minimum',
                'conduction band minumm',
                )
            gap                    = assigned_value(text,lower,'band gap')
            total_charge           = assigned_value(
                text,
                lower,
                'total charge in supercell',
                )
            total_magnetization    = assigned_value(
                text,
                lower,
                'total magnetization',
                )
            absolute_magnetization = assigned_value(
                text,
                lower,
                'absolute magnetization',
                )
            if fermi is not None:
                data.fermi_energies.append(fermi)
            elif vbm is not None and cbm is not None:
                data.valence_band_maxima.append(vbm)
                data.conduction_band_minima.append(cbm)
            elif gap is not None:
                data.band_gaps.append(gap)
            elif total_charge is not None:
                data.total_charges.append(total_charge)
            elif total_magnetization is not None:
                data.total_magnetizations.append(total_magnetization)
            elif absolute_magnetization is not None:
                data.absolute_magnetizations.append(absolute_magnetization)
            elif lower.startswith('sum force'):
                values = line_numbers(
                    text.partition('=')[2],
                    self.number_pattern,
                    )
                if len(values)>=3:
                    data.sum_forces.append(values[:3])
            elif 'volume and energy per atom' in lower:
                values = line_numbers(
                    text.partition('=')[2],
                    self.number_pattern,
                    )
                if len(values)>=2:
                    data.volume_per_atom.append(values[0])
                    data.energy_per_atom.append(values[1])

            if 'kohn sham eigenvalues' in lower and 'k-point' in lower:
                kpoint_start = lower.rfind('k-point')+len('k-point')
                kpoint_text  = text[kpoint_start:]
                index_text,separator,coordinates_text = kpoint_text.partition(']')
                if len(separator)==0 or '[' not in index_text:
                    continue
                try:
                    index = int(index_text.rsplit('[',1)[1].strip())
                except ValueError:
                    continue
                coordinate_tokens = coordinates_text.lstrip(' :').split()
                coordinates       = [as_float(v) for v in coordinate_tokens[:3]]
                if len(coordinates)!=3 or None in coordinates:
                    continue
                if index in dataset:
                    datasets.append(dataset)
                    dataset = dotdict()
                dataset[index] = dotdict(
                    kpoint   = np.array(coordinates,dtype=float),
                    channels = dotdict(),
                    )
                kpoint = index
                spin   = 'none'
                continue
            if kpoint is None:
                continue
            if 'spin up' in lower:
                spin = 'up'
                continue
            elif 'spin down' in lower:
                spin = 'down'
                continue
            row_prefix,row_separator,row_text = line.lstrip().partition(']')
            if (
                len(row_separator)==0
                or not row_prefix.lower().startswith('[kpt')
                ):
                continue
            pairs = pair_pattern.findall(row_text)
            if len(pairs)==0:
                continue
            channels = dataset[kpoint].channels
            eigs,occs = channels.setdefault(spin,[[],[]])
            for eigenvalue,occupation in pairs:
                eigenvalue = as_float(eigenvalue)
                occupation = as_float(occupation)
                if eigenvalue is not None and occupation is not None:
                    eigs.append(eigenvalue)
                    occs.append(occupation)

        for name,values in data.items():
            if values is not None:
                data[name] = np.array(values,dtype=float)
        data.energy_units          = 'eV'
        data.magnetization_units   = 'Bohr mag/cell'
        data.sum_force_units       = 'Ha/a0'
        data.energy_per_atom_units = 'eV'

        if len(dataset)>0:
            datasets.append(dataset)
        # Retain the final complete table with consistent spin and band counts.
        expected_kpoints = None
        if self.setup_info.k_points is not None:
            expected_kpoints = len(self.setup_info.k_points.kpoints_crystal)
        for candidate in reversed(datasets):
            indices = sorted(candidate)
            if indices!=list(range(len(indices))):
                continue
            if expected_kpoints is not None and len(indices)!=expected_kpoints:
                continue
            spin_channels = set()
            for record in candidate.values():
                spin_channels.update(record.channels)
            if spin_channels=={'none'}:
                spins = ['none']
            elif spin_channels=={'up','down'}:
                spins = ['up','down']
            else:
                continue
            channels = [
                candidate[index].channels.get(spin)
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
            data.kpoints_crystal = np.array(
                [candidate[i].kpoint for i in indices],
                dtype=float,
                )
            data.eigenvalues = np.array(
                [[candidate[i].channels[spin][0] for spin in spins]
                 for i in indices],
                dtype=float,
                )
            data.occupations = np.array(
                [[candidate[i].channels[spin][1] for spin in spins]
                 for i in indices],
                dtype=float,
                )
            if spins==['none']:
                data.eigenvalues = data.eigenvalues[:,0,:]
                data.occupations = data.occupations[:,0,:]
            if self.setup_info.structure is not None:
                data.kpoints = np.dot(
                    data.kpoints_crystal,
                    self.setup_info.structure.kaxes,
                    )
            break

        nfound = sum(v.size for v in data.values() if isinstance(v,np.ndarray))
        if nfound>0:
            self.electronic = data
    #end def read_electronic


    def read_scf(self,lines):
        """Read SCF energy components, iteration indices, residuals, and times.

        Binds ``scf`` to an ``obj`` containing NumPy histories. Energies are
        in Hartree and times are in seconds.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        component_names = {
            'eigenvalue sum'  : 'eigenvalue_sum',
            'ion_ion'         : 'ion_ion',
            'electrostatic'   : 'electrostatic',
            'vxc'             : 'vxc',
            'exc'             : 'exc',
            'total energy'    : 'total_energy',
            'estimated error' : 'estimated_error',
            }
        values = dotdict(
            eigenvalue_sum  = [],
            ion_ion         = [],
            electrostatic   = [],
            vxc             = [],
            exc             = [],
            total_energy    = [],
            estimated_error = [],
            )

        # Match an SCF energy-component label and its value or overflow stars.
        # Example: @@ TOTAL ENERGY = -1.250000 Ha
        component_pattern = re.compile(
            r'^\s*@@\s*(?P<label>eigenvalue\s+sum|ion_ion|electrostatic|'
            r'vxc|exc|total\s+energy|estimated\s+error)\s*[:=]\s*'
            r'(?P<value>'+self.number_pattern+r'|\*+)',
            re.IGNORECASE,
            )
        # Match fields within a detailed SCF-iteration summary.
        # Example: quench: [md: 0/2 scf: 3/20 step time: 0.20 RMS[dV]: 2e-5]
        detail_pattern = re.compile(
            r'\bmd\s*:\s*(?P<md>\d+)\s*/|'
            r'\bscf\s*:\s*(?P<scf>\d+)\s*/|'
            r'\bstep\s+time\s*:\s*(?P<step>'+self.number_pattern+r')|'
            r'\bscf\s+time\s*:\s*(?P<time>'+self.number_pattern+r')|'
            r'\brms\s*\[\s*dv\s*\]\s*:\s*(?P<rms>[^\]\s]+)',
            re.IGNORECASE,
            )
        md_steps   = []
        scf_steps  = []
        step_times = []
        scf_times  = []
        rms_dv     = []
        for line in lines:
            match = component_pattern.search(line)
            if match is not None:
                label = normalize_line(match.group('label')).lower()
                name  = component_names.get(label)
                token = match.group('value')
                if name is not None:
                    value = np.nan if '*' in token else as_float(token)
                    if value is not None:
                        values[name].append(value)
            summary_label,separator,_ = normalize_line(line).partition(':')
            if len(separator)==0 or summary_label.lower()!='quench':
                continue
            details = dotdict(
                md   = None,
                scf  = None,
                step = None,
                time = None,
                rms  = None,
                )
            for match in detail_pattern.finditer(line):
                details[match.lastgroup] = match.group(match.lastgroup)
            md_steps.append(int(details.md) if details.md is not None else -1)
            scf_steps.append(int(details.scf) if details.scf is not None else -1)
            step_times.append(
                as_float(details.step) if details.step is not None else np.nan,
                )
            scf_times.append(
                as_float(details.time) if details.time is not None else np.nan,
                )
            rms = as_float(details.rms) if details.rms is not None else None
            rms_dv.append(rms if rms is not None else np.nan)

        if len(values.total_energy)>0:
            scf = obj()
            for name,array in values.items():
                scf[name] = np.array(array,dtype=float)
            scf.md_steps     = np.array(md_steps,dtype=int)
            scf.scf_steps    = np.array(scf_steps,dtype=int)
            scf.step_times   = np.array(step_times,dtype=float)
            scf.scf_times    = np.array(scf_times,dtype=float)
            scf.rms_dv       = np.array(rms_dv,dtype=float)
            scf.energy_units = 'Ha'
            scf.time_units   = 's'
            self.scf         = scf
    #end def read_scf


    def read_ions(self,lines):
        """Read detailed ionic records and construct trajectory-level data.

        Binds ``ionic_steps`` to detailed per-step records and, for consistent
        atom counts, binds trajectory arrays and corresponding structures.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        records    = []
        structures = obj()
        initial    = self.setup_info.structure
        i       = 0
        # Collect complete ionic rows and construct each reported structure.
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
            atoms          = []
            positions      = []
            charges        = []
            magnetizations = []
            forces         = []
            movable        = []
            i += 1
            while i<len(lines):
                tokens = lines[i].split()
                if len(tokens)==0 or tokens[0].upper()!='@ION':
                    break
                i += 1
                if len(tokens)<14:
                    continue
                values = [as_float(v) for v in tokens[3:14]]
                if None in values:
                    continue
                move_values = values[8:11]
                if not all(value.is_integer() for value in move_values):
                    continue
                atoms.append(tokens[2])
                positions.append(values[:3])
                charges.append(values[3])
                magnetizations.append(values[4])
                forces.append(values[5:8])
                movable.append([int(value) for value in move_values])
            if len(atoms)>0:
                record = obj(
                    atoms          = np.array(atoms,dtype=object),
                    positions      = np.array(positions,dtype=float),
                    position_units = 'a0',
                    charges        = np.array(charges,dtype=float),
                    magnetizations = np.array(magnetizations,dtype=float),
                    forces         = np.array(forces,dtype=float),
                    force_units    = 'Ha/a0',
                    movable        = np.array(movable,dtype=int),
                    )
                records.append(record)
                if initial is not None:
                    structure = generate_structure(
                        units = 'B',
                        axes  = initial.axes,
                        elem  = record.atoms,
                        pos   = record.positions,
                        )
                    structure.add_kpoints(
                        initial.kpoints,
                        initial.kweights,
                        recenter=False,
                        )
                    structures[len(records)-1] = structure
        # Bind trajectories only when every ionic step has a consistent size.
        if len(records)>0:
            self.ionic_steps    = obj(dict(enumerate(records)))
            self.position_units = 'a0'
            self.force_units    = 'Ha/a0'
        if len(records)>0 and len({len(record.atoms) for record in records})==1:
            self.positions      = np.array(
                [record.positions for record in records],
                dtype=float,
                )
            self.forces         = np.array(
                [record.forces for record in records],
                dtype=float,
                )
            self.charges        = np.array(
                [record.charges for record in records],
                dtype=float,
                )
            self.magnetizations = np.array(
                [record.magnetizations for record in records],
                dtype=float,
                )
            self.max_forces     = np.array(
                [np.linalg.norm(record.forces,axis=1).max()
                 for record in records],
                dtype=float,
                )
            if len(structures)==len(records):
                self.structures = structures
    #end def read_ions


    def read_geometry(self):
        """Collect cell volume and k-point data from the setup report.

        Binds ``geometry`` to an ``obj`` containing volume and crystal and
        Cartesian k-point data.
        """
        geometry = obj(
            volume          = None,
            volume_units    = None,
            kpoints_crystal = None,
            kpoints_cart    = None,
            kweights        = None,
            )
        structure = self.setup_info.structure
        if structure is not None:
            geometry.volume       = abs(np.linalg.det(structure.axes))
            geometry.volume_units = 'a0^3'
            if len(structure.kpoints)>0:
                geometry.kpoints_cart = structure.kpoints
                geometry.kweights     = structure.kweights
        if self.setup_info.k_points is not None:
            kpoints                  = self.setup_info.k_points
            geometry.kpoints_crystal = kpoints.kpoints_crystal
            geometry.kweights        = kpoints.kweights
            if structure is not None and len(kpoints.kpoints_crystal)>0:
                geometry.kpoints_cart = np.dot(
                    kpoints.kpoints_crystal,
                    structure.kaxes,
                    )
        if any(value is not None for value in geometry.values()):
            self.geometry = geometry
    #end def read_geometry


    def read_stress(self,lines):
        """Parse stress tensors and derive hydrostatic pressures.

        Binds ``stress`` and ``pressures`` histories, their shared unit label,
        and ``pressure`` as the final hydrostatic pressure. Values are in kbar.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        tensors = []
        for i,line in enumerate(lines):
            normalized = normalize_line(line).lower()
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
                    value = as_float(token)
                    if value is None:
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
            stress    = np.array(tensors,dtype=float)
            pressures = -np.trace(stress,axis1=1,axis2=2)/3.0
            self.stress       = stress
            self.stress_units = 'kbar'
            self.pressures    = pressures
            self.pressure     = pressures[-1]
    #end def read_stress


    def read_convergence(self,lines):
        """Read electronic and ionic convergence messages.

        Binds ``convergence`` to an ``obj`` containing nullable status values
        and electronic success and failure counts.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        electronic_successes = 0
        electronic_failures  = 0
        ionic_converged      = None
        for line in lines:
            text         = normalize_line(line).lower()
            not_achieved = (
                'not achieved' in text or 'not been achieved' in text)
            electronic_failure = (
                'potential convergence' in text
                and not_achieved
                or 'convergence criterion' in text
                and 'not met' in text
                )
            electronic_success = (
                'potential convergence' in text
                and 'achieved' in text
                and not electronic_failure
                )
            ionic_failure = (
                'force convergence' in text
                and not_achieved
                )
            ionic_success = (
                'force convergence' in text
                and 'achieved' in text
                and not ionic_failure
                )
            if electronic_failure:
                electronic_failures += 1
            elif electronic_success:
                electronic_successes += 1
            elif ionic_failure:
                ionic_converged = False
            elif ionic_success:
                ionic_converged = True
        electronic_converged = None
        if electronic_successes+electronic_failures>0:
            electronic_converged = (
                electronic_successes>0 and electronic_failures==0)
        if electronic_converged is not None or ionic_converged is not None:
            self.convergence = obj(
                electronic_converged = electronic_converged,
                electronic_successes = electronic_successes,
                electronic_failures  = electronic_failures,
                ionic_converged      = ionic_converged,
                )
    #end def read_convergence


    def read_timing(self,lines):
        """Read total, per-step, and section-resolved timing information.

        Binds ``timing`` to an ``obj`` measured in seconds, with individual
        timing rows stored under ``timing.sections``.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.
        """
        # Match a numbered timing section followed by total and per-step times.
        # Example: 1-TOTAL 3.00 0.50
        time_value = r'(?:'+self.number_pattern+r'|inf|nan)'
        pattern    = re.compile(
            r'^\s*(\d+\s*-\s*.*?)\s+('+time_value+r')\s+('+time_value+r')'
            r'(?:\s+.*)?$',
            re.IGNORECASE,
            )
        timing   = None
        sections = obj()
        for line in lines:
            match = pattern.match(line)
            if match is None:
                continue
            name                    = normalize_line(match.group(1))
            prefix,separator,suffix = name.partition('-')
            if len(separator)>0:
                name = prefix.strip()+'-'+suffix.strip()
            total    = float(match.group(2).lower().replace('d','e'))
            per_step = float(match.group(3).lower().replace('d','e'))
            # Replace non-alphanumeric runs to form a stable section key.
            # Example: 1-TOTAL
            key = re.sub(r'[^a-z0-9]+','_',name.lower()).strip('_')
            sections[key] = obj(total=total,per_step=per_step)
            if name.lower()=='1-total':
                timing = obj(
                    total    = total,
                    per_step = per_step,
                    units    = 's',
                    )
        if timing is not None:
            timing.sections = sections
            self.timing     = timing
    #end def read_timing


    def read_produced_files(self):
        """Locate recognized files produced by QMCPACK-interface SCF runs.

        Binds ``produced_files`` to an ``obj`` containing the QMCPACK restart
        path when the file exists.
        """
        produced_files = obj()
        files          = self.setup_info.files
        if files is not None and files.data_output_file is not None:
            qmcpack_file = os.path.join(
                self.path,
                str(files.data_output_file)+'.h5',
                )
            if os.path.isfile(qmcpack_file):
                produced_files.qmcpack_restart = qmcpack_file
        if len(produced_files)>0:
            self.produced_files = produced_files
    #end def read_produced_files
#end class RmgOutData



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
        RMG input reconstructed from the referenced control file when it is
        available.
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
        if self.results.setup_info.structure is None:
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
        electronic = self.results.electronic
        if electronic is not None and electronic.kpoints is not None:
            kpoints = electronic.kpoints
        else:
            geometry = self.results.geometry
            if geometry is None or geometry.kpoints_cart is None:
                return None
            kpoints = geometry.kpoints_cart
        return kpoints*convert(1.0,units,'B')
    #end def kpoints


    def kweights(self):
        """Return dimensionless k-point weights, or ``None`` if unavailable."""
        self._require_supported('kweights',self.all_modes)
        geometry = self.results.geometry
        if geometry is None or geometry.kweights is None:
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
        if electronic is None or electronic.eigenvalues is None:
            return None
        return convert(electronic.eigenvalues,'eV',units)
    #end def eigenvalues


    def occupations(self):
        """Return the dimensionless K-point-major occupation array."""
        self._require_supported('occupations',self.all_modes)
        electronic = self.results.electronic
        if electronic is None or electronic.occupations is None:
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
        empty           = np.isclose(
            occupations,
            0.0,
            rtol = 0.0,
            atol = tolerance,
            )
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
            msg       = f'stress units must be one of: {supported}'
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
            msg       = f'pressure units must be one of: {supported}'
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
        filepath = os.path.join(self.path,self.outfile_name)
        results  = RmgOutData(filepath)
        self.results  = results
        self.run_mode = results.run_mode
        self.input    = results.input
    #end def analyze

#end class RmgAnalyzer
