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
from .rmg_input import RmgInput, rmg_modes
from .simulation import Simulation, SimulationAnalyzer
from .structure import generate_structure
from .unit_converter import convert


class RmgOutData(DevBase):
    """Read an RMG output file and collect results appropriate to its run mode.

    Parameters
    ----------
    filepath : str
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
        Parsed control input when the referenced input file can be found and
        read.
    setup_info : obj
        Parsed setup sections, including the run mode and, when available,
        lattice, ion, k-point, grid, and file information.
    run_mode : str or None
        Short RMG calculation mode: ``"scf"``, ``"nscf"``, or ``"relax"``.
    geometry : obj or None
        Cell volume and k-point arrays.
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
        Unit label corresponding to each entry in ``energies``.
    direct_energies : numpy.ndarray or None
        History of directly evaluated total energies.
    direct_energy_units : numpy.ndarray or None
        Unit label corresponding to each direct energy.
    electronic : obj or None
        Electronic observables represented by NumPy arrays, including Fermi
        energies, band edges, gaps, charges, magnetizations, forces, volumes,
        per-atom energies, k-points, Kohn--Sham eigenvalues, and occupations
        when reported.
    scf : obj or None
        SCF energy components, iteration indices, residuals, and timing data;
        numerical histories are stored as NumPy arrays.
    ionic_steps : obj or None
        Mapping from ionic-step index to an ``obj`` containing atom labels and
        position, charge, magnetization, force, and movable-flag arrays.
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
        Mapping from ionic-step index to a ``Structure`` instance.
    stress : numpy.ndarray or None
        Stress tensors with shape ``(nsteps, 3, 3)``.
    stress_units : str or None
        Units associated with stress and pressure values.
    pressures : numpy.ndarray or None
        Hydrostatic pressure at each reported stress step.
    pressure : float or numpy.floating or None
        Last hydrostatic pressure.
    produced_files : obj or None
        Paths for SCF output files, such as a QMCPACK restart file.

    Notes
    -----
    The supported modes are ``scf``, ``nscf``, and ``relax``. A
    mode-applicable member is initialized to ``None`` and remains ``None``
    when its data cannot be obtained. ``produced_files`` applies only to
    ``scf``.

    Raises
    ------
    TypeError
        If ``filepath`` is not a string.
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
        if not isinstance(filepath,str):
            raise TypeError(
                'invalid type provided for filepath\n'
                'Type expected: str\n'
                f'Type provided: {filepath.__class__.__name__}'
                )
        elif not os.path.exists(filepath):
            raise FileNotFoundError(
                'RMG log output file does not exist.\n'
                f'Path provided: {filepath}'
                )
        elif not os.path.isfile(filepath):
            raise IsADirectoryError(
                'Path provided for RMG log output is not a file.\n'
                f'Path provided: {filepath}'
                )
        path,outfile_name  = os.path.split(filepath)
        self.path          = path
        self.abspath       = os.path.abspath(path)
        self.outfile_name  = outfile_name
        self.input         = None
        self.setup_info    = None

        with open(filepath,'r') as fobj:
            lines = fobj.read().splitlines()
        self.read_setup_info(lines)

        supported_modes = {'scf','nscf','relax'}

        # modes: scf, nscf, relax
        if self.run_mode in supported_modes:
            self.geometry    = None
            self.convergence = None
            self.timing      = None
            self.energy               = None
            self.energy_units         = None
            self.energies             = None
            self.energy_units_history = None
            self.direct_energies      = None
            self.direct_energy_units  = None
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
            self.electronic            = None

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
        """Parse the RMG setup report and referenced control input.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines without newline characters.

        Notes
        -----
        Binds ``setup_info`` as an ``obj`` containing named setup blocks and
        derived NumPy arrays for grids, lattice vectors, k-points, and ionic
        positions. It also binds ``run_mode`` and may bind ``input`` to an
        ``RmgInput`` instance.
        """
        setup_info = obj()
        log_text = '\n'.join(lines)

        def normalize_line(line):
            """Normalize case and whitespace for tolerant text matching."""
            return ' '.join(line.split()).lower()
        #end def normalize_line

        def line_numbers(line):
            """Return all RMG-formatted numbers in a line as a float array."""
            # Find every RMG-formatted number occurring anywhere in a line.
            # Example: SUM FORCE = 0.1 0.2 0.3
            values = re.findall(
                self.number_pattern,line.replace('D','E').replace('d','e'))
            return np.array(values,dtype=float)
        #end def line_numbers

        def leading_numbers(line):
            """Return the leading whitespace-separated numbers as a float array."""
            # Match a whitespace-separated numeric sequence at the start of a line.
            # Example: 1 1.0 0.1 0.2 trailing
            npat    = self.number_pattern
            pattern = r'^\s*((?:'+npat+r')(?:\s+(?:'+npat+r'))*)'
            match   = re.match(pattern,line)
            if match is None:
                return np.array([],dtype=float)
            return line_numbers(match.group(1))
        #end def leading_numbers

        def match_float(pattern,line):
            """Return the first floating-point value captured by a regex."""
            # Apply a caller-provided expression that captures one physical value.
            # Example: Grid spacing: 0.30 a0
            match = re.search(pattern,line,re.IGNORECASE)
            if match is None:
                return None
            return float(match.group(1).replace('D','E').replace('d','e'))
        #end def match_float

        def identify_mode(text):
            """Identify the short RMG run-mode name from descriptive text."""
            text          = normalize_line(text)
            mode_patterns = (
                # Match the self-consistent electronic-quench calculation label.
                # Example: Quench electrons
                (r'quench\s+electrons','scf'),
                # Match the standalone abbreviation for a non-self-consistent run.
                # Example: NSCF
                (r'\bnscf\b','nscf'),
                # Match either supported structural-relaxation label.
                # Example: Structure optimization
                (r'(?:structure\s+optimization|relax\s+structure)','relax'),
                )
            for pattern,run_mode in mode_patterns:
                if re.search(pattern,text,re.IGNORECASE):
                    return run_mode
            run_mode = rmg_modes.mode_match(text,short=True)
            return run_mode if run_mode in {'scf','nscf','relax'} else None
        #end def identify_mode

        mode = None
        for line in lines:
            normalized = normalize_line(line)
            # Match the calculation-type heading and capture its value.
            # Example: Calculation type: Quench electrons
            match = re.match(r'^calculation\s*type\s*[:=]\s*(.*)$',normalized)
            if match is not None:
                mode = identify_mode(match.group(1))
                break
        self.run_mode       = mode
        setup_info.run_mode = mode
        unit_set = {'a0'}
        on_off   = {
            'ON'  : True,
            'OFF' : False,
            }
        def process_name(s):
            """Convert an RMG setup label into a normalized object key."""
            # Match parenthetical annotations so they can be removed from names.
            # Example: Grid spacing (a0)
            s = re.sub(r'\([^)]*\)','',s)
            name = '_'.join(s.strip().lower().split())
            return name.replace('/','_').replace('-','_')
        #end def process_name
        def process_value(v,*,list=False):
            """Convert setup text to a scalar, array, list, or string and units."""
            v     = v.strip()
            units = None
            try:
                v = int(v)
            except ValueError:
                try:
                    v = float(v)
                except ValueError:
                    if ' ' in v or ',' in v:
                        tokens = v.replace(',',' ').split()
                        if len(tokens)>0:
                            if tokens[-1] in unit_set:
                                units  = tokens[-1]
                                tokens = tokens[:-1]
                            try:
                                if not list:
                                    v = np.array(tokens,dtype=float)
                                else:
                                    v = [process_value(t,list=True)[0] for t in tokens]
                            except ValueError:
                                units = None
                    elif v in on_off:
                        v = on_off[v]
            return v,units
        #end def process_value
        if mode is not None:
            # Match the standalone Files setup-section heading.
            # Example: Files:
            setup_flags = re.IGNORECASE|re.MULTILINE
            start_match = re.search(
                r'^\s*files\s*(?::\s*)?$',log_text,setup_flags)
            istart      = start_match.start() if start_match is not None else -1
            if istart!=-1:
                # Match the normal heading that terminates the setup section.
                # Example: Initial Ionic Positions And Displacements (Angstrom)
                end_pattern = (
                    r'^\s*initial\s+ionic\s+positions\s+and\s+'
                    r'displacements\s*\(\s*angstroms?\s*\)'
                    )
                end_match   = re.search(
                    end_pattern,log_text[istart:],setup_flags)
                if end_match is not None:
                    iend = istart+end_match.start()
                else:
                    # Match alternate post-setup headings when ionic positions are absent.
                    # Example: Davidson converged
                    fallback_pattern = (
                        r'^\s*(?:diagonalization\s+using|davidson\s+'
                        r'(?:converged|incomplete)|-+\s*timing\s+information)'
                        )
                    end_match        = re.search(
                        fallback_pattern,log_text[istart:],setup_flags)
                    iend             = istart+end_match.start() if end_match is not None else len(log_text)
                if iend!=-1:
                    blocks       = []
                    body_started = False
                    for line in log_text[istart:iend].expandtabs().splitlines():
                        if len(line)==0:
                            continue
                        if line[0]!=' ' and (len(blocks)==0 or body_started):
                            blocks.append([line])
                        else:
                            blocks[-1].append(line)
                        body_started = line[0]==' '
                    other_blocks = dotdict()
                    table_blocks = {
                        'k_points',
                        'initial_ionic_positions_and_displacements',
                        }
                    for block in blocks:
                        if len(block)<2:
                            continue
                        header = block[0]
                        bname  = process_name(header)
                        lines  = block[1:]
                        if (
                            bname not in table_blocks and
                            any(':' in line for line in lines)
                            ):
                            bvalues = obj()
                            for line in lines:
                                if ':' not in line:
                                    continue
                                name,value    = line.split(':',1)
                                name          = process_name(name)
                                value,units   = process_value(value)
                                bvalues[name] = value
                                if units is not None:
                                    bvalues.units = units
                            setup_info[bname] = bvalues
                        elif bname in table_blocks:
                            other_blocks[bname] = header,lines
                    # additional processing for specific blocks
                    if 'grid_points' in setup_info:
                        b          = setup_info.grid_points
                        grid       = []
                        grid_pe    = []
                        spacing    = []
                        grid_units = None
                        for c in sorted({'x','y','z'}):
                            if c not in b:
                                continue
                            s = str(b[c])
                            # Match and capture the total grid-point count.
                            # Example: Total: 48 Per PE: 12 Spacing: 0.30 a0
                            total  = match_float(
                                r'\btotal\s*:\s*('+self.number_pattern+r')',s)
                            # Match and capture the per-process grid-point count.
                            # Example: Total: 48 Per PE: 12 Spacing: 0.30 a0
                            per_pe = match_float(
                                r'\bper\s*pe\s*:\s*('+self.number_pattern+r')',s)
                            # Match and capture the real-space grid spacing.
                            # Example: Total: 48 Per PE: 12 Spacing: 0.30 a0
                            space  = match_float(
                                r'\bspacing\s*:\s*('+self.number_pattern+r')',s)
                            if total is None or per_pe is None or space is None:
                                continue
                            # Match the unit following a grid-spacing value.
                            # Example: Spacing: 0.30 a0
                            unit_match = re.search(
                                r'\bspacing\s*:\s*'+self.number_pattern+
                                r'\s*([^\s,;]+)',s,re.IGNORECASE)
                            if unit_match is not None:
                                grid_units = unit_match.group(1)
                            grid.append(total)
                            grid_pe.append(per_pe)
                            spacing.append(space)
                        if len(grid)==3:
                            grid    = np.array(grid,dtype=int)
                            grid_pe = np.array(grid_pe,dtype=int)
                            spacing = np.array(spacing,dtype=float)
                            b.update(
                                grid         = grid,
                                grid_pe      = grid_pe,
                                grid_spacing = spacing,
                                grid_units   = grid_units,
                                )
                        if 'equivalent_energy_cutoffs' in b:
                            cutoff_text   = str(b.equivalent_energy_cutoffs)
                            cutoff_values = line_numbers(cutoff_text)
                            if len(cutoff_values)>=2:
                                # Match the unit following the two equivalent cutoffs.
                                # Example: 60.0 240.0 Ry
                                unit_match = re.search(
                                    self.number_pattern+r'\s+'+self.number_pattern+
                                    r'\s*([^\s,;]+)',cutoff_text)
                                ecut_units = (
                                    unit_match.group(1) if unit_match is not None else None)
                                b.update(
                                    ecut         = cutoff_values[0],
                                    ecut_charge  = cutoff_values[1],
                                    ecut_units   = ecut_units,
                                    )
                    if 'lattice_setup' in setup_info:
                        b    = setup_info.lattice_setup
                        axes = []
                        for name in sorted({
                            'x_basis_vector','y_basis_vector','z_basis_vector',
                            }):
                            if name not in b:
                                axes = []
                                break
                            values = line_numbers(str(b[name]))
                            if len(values)<3:
                                axes = []
                                break
                            axes.append(values[:3])
                        if len(axes)==3:
                            b.axes = np.array(axes,dtype=float)
                    if 'k_points' in other_blocks:
                        _,lines = other_blocks.k_points
                        first_row = None
                        for i,line in enumerate(lines):
                            normalized = normalize_line(line)
                            if 'weight' in normalized and 'crystal' in normalized:
                                first_row = i+1
                                break
                        if first_row is not None:
                            kp = []
                            kw = []
                            for line in lines[first_row:]:
                                normalized = normalize_line(line)
                                if 'weight' in normalized and len(kp)>0:
                                    break
                                values = leading_numbers(line)
                                if len(values)<4:
                                    continue
                                kp.append(values[:3])
                                kw.append(values[3])
                            if len(kp)>0:
                                setup_info.k_points = obj(
                                    kpoints_crystal = np.array(kp,dtype=float),
                                    kweights        = np.array(kw,dtype=float),
                                    )
                    k = 'initial_ionic_positions_and_displacements'
                    if k in other_blocks:
                        header,lines = other_blocks[k]
                        h      = header.lower()
                        punits = None
                        if 'bohr' in h:
                            punits = 'B'
                        elif 'angstrom' in h:
                            punits = 'A'
                        pos       = []
                        spec      = []
                        first_row = None
                        for i,line in enumerate(lines):
                            if 'species' in normalize_line(line):
                                first_row = i+1
                                break
                        if first_row is not None:
                            for line in lines[first_row:]:
                                t = line.split()
                                if len(t)<4:
                                    continue
                                values = line_numbers(' '.join(t[1:]))
                                if len(values)<3:
                                    continue
                                spec.append(t[0])
                                pos.append(values[:3])
                            if len(pos)>0:
                                setup_info.ion_positions = obj(
                                    units     = punits,
                                    atoms     = np.array(spec,dtype=object),
                                    positions = np.array(pos,dtype=float),
                                    )
        have_structure_data = (
            'lattice_setup' in setup_info and
            'axes' in setup_info.lattice_setup and
            'ion_positions' in setup_info and
            setup_info.ion_positions.get('units',None) in {'A','B'}
            )
        if have_structure_data:
            aunits       = setup_info.lattice_setup.get('units','B')
            axes         = np.asarray(setup_info.lattice_setup.axes,dtype=float)
            elem         = np.asarray(setup_info.ion_positions.atoms,dtype=object)
            pos          = np.asarray(setup_info.ion_positions.positions,dtype=float)
            punits       = setup_info.ion_positions.units
            valid_shapes = (
                axes.shape==(3,3) and
                pos.ndim==2 and pos.shape[1:]==(3,) and
                len(elem)==len(pos) and len(elem)>0
                )
            if valid_shapes:
                if aunits in {'a0','B'}:
                    aunits = 'B'
                else:
                    aunits = 'A' # assume for now
                units = 'B'
                axes  = convert(axes,aunits,units)
                pos   = convert(pos,punits,units)
                s     = generate_structure(
                    units = units,
                    axes  = axes,
                    elem  = elem,
                    pos   = pos,
                    )
                if 'k_points' in setup_info and 'kpoints_crystal' in setup_info.k_points:
                    kpu = setup_info.k_points.kpoints_crystal
                    if len(kpu)>0:
                        kw  = setup_info.k_points.kweights
                        kp  = np.dot(kpu,s.kaxes)
                        s.add_kpoints(kpoints=kp,kweights=kw)
                setup_info.structure = s
        if 'files' in setup_info and 'control_input_file' in setup_info.files:
            control_file = setup_info.files.control_input_file
            filepaths    = [
                os.path.join(self.path,control_file),
                os.path.join(self.path,os.path.basename(control_file)),
                os.path.join(os.path.dirname(self.path),control_file),
                ]
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
        """Parse final total-energy histories from the RMG log.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``energies`` and ``energy_units_history`` as NumPy arrays and
        ``energy`` and ``energy_units`` to their final entries. When direct
        energies are present, ``direct_energies`` and ``direct_energy_units``
        are also bound as NumPy arrays.
        """
        energies            = []
        direct_energies     = []
        energy_units        = []
        direct_energy_units = []
        # Match a final total energy and classify its eigenvalue-sum or direct source.
        # Example: final total energy from eigenvalue sum : -1.2345 Ha
        energy_pattern = re.compile(
            r'final\s+total\s+energy\s+from\s+'
            r'(?:(?P<eigenvalue>eig(?:envalue)?\s+sum)|(?P<direct>direct))'
            r'\s*[:=]\s*(?P<value>'+self.number_pattern+r')'
            r'(?:\s+(?P<units>[^\s,;]+))?',
            re.IGNORECASE,
            )
        for line in lines:
            match = energy_pattern.search(line)
            if match is not None:
                value = float(
                    match['value'].replace('D','E').replace('d','e'))
                if match['eigenvalue'] is not None:
                    energies.append(value)
                    energy_units.append(match['units'])
                else:
                    direct_energies.append(value)
                    direct_energy_units.append(match['units'])
        energies        = np.array(energies,dtype=float)
        direct_energies = np.array(direct_energies,dtype=float)
        if len(energies)>0:
            self.energy               = energies[-1]
            self.energy_units         = energy_units[-1] or 'Ha'
            self.energies             = energies
            self.energy_units_history = np.array(energy_units,dtype=object)
        if len(direct_energies)>0:
            self.direct_energies     = direct_energies
            self.direct_energy_units = np.array(direct_energy_units,dtype=object)
    #end def read_energies


    def read_electronic(self,lines):
        """Parse general electronic observables from the RMG log.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``electronic`` to an ``obj`` containing NumPy arrays for Fermi
        energies, band edges, band gaps, total charges, magnetizations, summed
        forces, volumes, and per-atom energies, together with their unit labels.
        The final complete Kohn--Sham table additionally supplies k-point-major
        eigenvalue and occupation arrays and crystal and Cartesian k-points.
        """
        def line_numbers(line):
            """Return all RMG-formatted numbers in a line as a float array."""
            # Find every RMG-formatted number occurring anywhere in a line.
            # Example: SUM FORCE = 0.1 0.2 0.3
            values = re.findall(
                self.number_pattern,line.replace('D','E').replace('d','e'))
            return np.array(values,dtype=float)
        #end def line_numbers

        def as_float(value):
            """Convert an RMG-formatted numeric token to a float."""
            return float(value.replace('D','E').replace('d','e'))
        #end def as_float

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
                        float(match.group(i).replace('D','E').replace('d','e'))
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
                    eigs.append(as_float(eigenvalue))
                    occs.append(as_float(occupation))
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
                    channel is None or len(channel[0])==0 or
                    len(channel[0])!=len(channel[1])
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
                return dotdict(
                    kpoints_crystal = kpoints,
                    eigenvalues     = eigenvalues,
                    occupations     = occupations,
                    )
            return None
        #end def read_eigenvalue_data

        data = obj(
            fermi_energies          = [],
            valence_band_maxima     = [],
            conduction_band_minima  = [],
            band_gaps               = [],
            total_charges           = [],
            total_magnetizations    = [],
            absolute_magnetizations = [],
            )

        sum_forces      = []
        volume_per_atom = []
        energy_per_atom = []
        npat            = self.number_pattern
        # Match and classify scalar electronic records, forces, and atom data.
        # Example: FERMI ENERGY : 5.25 eV
        electronic_pattern = re.compile(
            r'(?P<fermi>fermi\s+energy\s*[:=]\s*'
            r'(?P<fermi_value>'+npat+r'))|'
            r'(?P<band_edges>^(?=.*valence\s+band\s+maximum\s*[:=]\s*'
            r'(?P<vbm>'+npat+r'))(?=.*conduction\s+band\s+'
            r'(?:minimum|minumm)\s*[:=]\s*(?P<cbm>'+npat+r')).*$)|'
            r'(?P<band_gap>band\s+gap\s*[:=]\s*'
            r'(?P<band_gap_value>'+npat+r'))|'
            r'(?P<total_charge>total\s+charge\s+in\s+supercell\s*[:=]\s*'
            r'(?P<total_charge_value>'+npat+r'))|'
            r'(?P<total_magnetization>@@\s*total\s+magnetization\s*[:=]\s*'
            r'(?P<total_magnetization_value>'+npat+r'))|'
            r'(?P<absolute_magnetization>'
            r'@@\s*absolute\s+magnetization\s*[:=]\s*'
            r'(?P<absolute_magnetization_value>'+npat+r'))|'
            r'(?P<sum_force>^\s*sum\s+force\b.*$)|'
            r'(?P<atom_data>volume\s+and\s+energy\s+per\s+atom\b.*$)',
            re.IGNORECASE,
            )
        for line in lines:
            match = electronic_pattern.search(line)
            if match is None:
                continue
            if match['fermi'] is not None:
                data.fermi_energies.append(as_float(match['fermi_value']))
            elif match['band_edges'] is not None:
                data.valence_band_maxima.append(as_float(match['vbm']))
                data.conduction_band_minima.append(as_float(match['cbm']))
            elif match['band_gap'] is not None:
                data.band_gaps.append(as_float(match['band_gap_value']))
            elif match['total_charge'] is not None:
                data.total_charges.append(as_float(match['total_charge_value']))
            elif match['total_magnetization'] is not None:
                data.total_magnetizations.append(
                    as_float(match['total_magnetization_value']))
            elif match['absolute_magnetization'] is not None:
                data.absolute_magnetizations.append(
                    as_float(match['absolute_magnetization_value']))
            elif match['sum_force'] is not None:
                values = line_numbers(line.split('=',1)[-1])
                if len(values)>=3:
                    sum_forces.append(values[:3])
            elif match['atom_data'] is not None:
                values = line_numbers(line.split('=',1)[-1])
                if len(values)>=2:
                    volume_per_atom.append(values[0])
                    energy_per_atom.append(values[1])
        for name,values in data.items():
            data[name] = np.array(values,dtype=float)
        data.energy_units          = 'eV'
        data.magnetization_units   = 'Bohr mag/cell'
        data.sum_forces            = np.array(sum_forces,dtype=float)
        data.sum_force_units       = 'Ha/a0'
        data.volume_per_atom       = np.array(volume_per_atom,dtype=float)
        data.energy_per_atom       = np.array(energy_per_atom,dtype=float)
        data.energy_per_atom_units = 'eV'

        eigenvalue_data = read_eigenvalue_data()
        if eigenvalue_data is not None:
            data.kpoints_crystal = eigenvalue_data.kpoints_crystal
            data.eigenvalues     = eigenvalue_data.eigenvalues
            data.occupations     = eigenvalue_data.occupations
            if 'structure' in self.setup_info:
                data.kpoints = np.dot(
                    data.kpoints_crystal,self.setup_info.structure.kaxes)

        nfound = sum(len(v) for v in data.values() if isinstance(v,np.ndarray))
        if nfound>0:
            self.electronic = data
    #end def read_electronic


    def read_scf(self,lines):
        """Parse SCF energy components, iteration indices, residuals, and times.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``scf`` to an ``obj`` whose numerical members are NumPy arrays.
        Energy-component arrays use Hartree and timing arrays use seconds.
        """
        component_names = {
            'eigenvalue_sum'  : 'EIGENVALUE SUM',
            'ion_ion'         : 'ION_ION',
            'electrostatic'   : 'ELECTROSTATIC',
            'vxc'             : 'VXC',
            'exc'             : 'EXC',
            'total_energy'    : 'TOTAL ENERGY',
            'estimated_error' : 'estimated error',
            }

        # Match each @@ energy component and capture its value or overflow stars.
        # Example: @@ TOTAL ENERGY : -1.250000
        component_patterns = {
            name : re.compile(
                r'^\s*@@\s*'+re.escape(label)+r'\s*[:=]\s*'
                r'('+self.number_pattern+r'|\*+)',
                re.IGNORECASE,
                )
            for name,label in component_names.items()
            }
        # Match and classify fields in a detailed SCF-iteration summary.
        # Example: quench: [ RMS [ dV ] : 2.0D-5 scf: 3/20 ]
        detail_pattern = re.compile(
            r'\bmd\s*:\s*(?P<md>\d+)\s*/|'
            r'\bscf\s*:\s*(?P<scf>\d+)\s*/|'
            r'\bstep\s+time\s*:\s*(?P<step>'+self.number_pattern+r')|'
            r'\bscf\s+time\s*:\s*(?P<time>'+self.number_pattern+r')|'
            r'\brms\s*\[\s*dv\s*\]\s*:\s*(?P<rms>[^\]\s]+)',
            re.IGNORECASE,
            )
        values     = {name: [] for name in component_names}
        md_steps   = []
        scf_steps  = []
        step_times = []
        scf_times  = []
        rms_dv     = []
        for line in lines:
            for name,pattern in component_patterns.items():
                match = pattern.search(line)
                if match is not None:
                    token = match.group(1)
                    if '*' in token:
                        value = np.nan
                    else:
                        value = float(token.replace('D','E').replace('d','e'))
                    values[name].append(value)
                    break
            # Match the detailed SCF-iteration summary line.
            # Example: quench: [ RMS [ dV ] : 2.0D-5 scf: 3/20 ]
            if re.search(r'\bquench\s*:',line,re.IGNORECASE):
                details = {
                    match.lastgroup : match.group(match.lastgroup)
                    for match in detail_pattern.finditer(line)
                    }
                md_steps.append(int(details.get('md',-1)))
                scf_steps.append(int(details.get('scf',-1)))
                step_times.append(float(
                    details.get('step','nan').replace('D','E').replace('d','e')))
                scf_times.append(float(
                    details.get('time','nan').replace('D','E').replace('d','e')))
                rms = details.get('rms','*')
                try:
                    rms_dv.append(float(rms.replace('D','E').replace('d','e')))
                except ValueError:
                    rms_dv.append(np.nan)
        scf = obj()
        for name,array in values.items():
            scf[name] = np.array(array,dtype=float)
        scf.update(
            md_steps     = np.array(md_steps,dtype=int),
            scf_steps    = np.array(scf_steps,dtype=int),
            step_times   = np.array(step_times,dtype=float),
            scf_times    = np.array(scf_times,dtype=float),
            rms_dv       = np.array(rms_dv,dtype=float),
            energy_units = 'Ha',
            time_units   = 's',
            )
        nfound = len(scf.total_energy)
        if nfound>0:
            self.scf = scf
    #end def read_scf


    def read_ions(self,lines):
        """Parse ionic-step tables and construct trajectory-level data.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``ionic_steps`` to an ``obj`` mapping step indices to records that
        contain atom-label, position, charge, magnetization, force, and movable-
        flag arrays. For a consistent atom count, it also binds trajectory NumPy
        arrays through ``positions``, ``forces``, ``charges``,
        ``magnetizations``, and ``max_forces`` and maps steps to ``Structure``
        instances through ``structures``.
        """
        records = []
        i       = 0
        while i<len(lines):
            header_tokens = lines[i].split()
            is_header     = (
                len(header_tokens)>=3 and
                header_tokens[0].upper()=='@ION' and
                header_tokens[1].lower()=='ion' and
                header_tokens[2].lower()=='species'
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
                numeric_tokens = tokens[3:14]
                # Match a complete numeric token before converting an ionic row.
                # Example: 2.1
                valid = all(
                    re.fullmatch(self.number_pattern,v) is not None
                    for v in numeric_tokens)
                if not valid:
                    continue
                values      = [
                    float(v.replace('D','E').replace('d','e'))
                    for v in numeric_tokens]
                move_values = values[8:11]
                if not all(np.isfinite(v) and v.is_integer() for v in move_values):
                    continue
                atoms.append(tokens[2])
                positions.append(values[:3])
                charges.append(values[3])
                magnetizations.append(values[4])
                forces.append(values[5:8])
                movable.append([int(v) for v in move_values])
            if len(atoms)>0:
                records.append(obj(
                    atoms          = np.array(atoms,dtype=object),
                    positions      = np.array(positions,dtype=float),
                    position_units = 'a0',
                    charges        = np.array(charges,dtype=float),
                    magnetizations = np.array(magnetizations,dtype=float),
                    forces         = np.array(forces,dtype=float),
                    force_units    = 'Ha/a0',
                    movable        = np.array(movable,dtype=int),
                    ))
        if len(records)>0:
            self.ionic_steps    = obj(dict(enumerate(records)))
            self.position_units = 'a0'
            self.force_units    = 'Ha/a0'
            if len({len(r.atoms) for r in records})==1:
                self.positions      = np.array(
                    [r.positions for r in records],dtype=float)
                self.forces         = np.array(
                    [r.forces for r in records],dtype=float)
                self.charges        = np.array(
                    [r.charges for r in records],dtype=float)
                self.magnetizations = np.array(
                    [r.magnetizations for r in records],dtype=float)
                self.max_forces     = np.array([
                    np.sqrt((r.forces**2).sum(axis=1)).max() for r in records
                    ],dtype=float)
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
        """Collect cell volume and k-point data derived from the setup report.

        Notes
        -----
        Binds ``geometry`` to an ``obj`` containing a scalar volume and NumPy
        arrays for crystal k-points, Cartesian k-points, and k-point weights.
        """
        geometry = obj()
        if 'structure' in self.setup_info:
            structure             = self.setup_info.structure
            geometry.volume       = abs(np.linalg.det(structure.axes))
            geometry.volume_units = 'a0^3'
        if 'k_points' in self.setup_info:
            kpoints                  = self.setup_info.k_points
            geometry.kpoints_crystal = kpoints.kpoints_crystal
            geometry.kweights        = kpoints.kweights
            if 'structure' in self.setup_info and len(kpoints.kpoints_crystal)>0:
                geometry.kpoints_cart = np.dot(
                    kpoints.kpoints_crystal,self.setup_info.structure.kaxes)
        if len(geometry)>0:
            self.geometry = geometry
    #end def read_geometry


    def read_stress(self,lines):
        """Parse stress tensors and derive hydrostatic pressures.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``stress`` to a NumPy array of shape ``(nsteps, 3, 3)``,
        ``pressures`` to a one-dimensional NumPy array, and ``pressure`` to the
        final pressure. Stress and pressure values are reported in kbar.
        """
        def leading_numbers(line):
            """Return the leading whitespace-separated numbers as a float array."""
            # Match a whitespace-separated numeric sequence at the start of a line.
            # Example: 1 1.0 0.1 0.2 trailing
            npat    = self.number_pattern
            pattern = r'^\s*((?:'+npat+r')(?:\s+(?:'+npat+r'))*)'
            match   = re.match(pattern,line)
            if match is None:
                return np.array([],dtype=float)
            # Find every RMG-formatted number in the matched numeric sequence.
            # Example: 1 1.0 0.1 0.2
            values = re.findall(
                self.number_pattern,
                match.group(1).replace('D','E').replace('d','e'))
            return np.array(values,dtype=float)
        #end def leading_numbers

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
                values = leading_numbers(row_line)
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
            stress            = np.array(tensors,dtype=float)
            pressures         = -np.trace(stress,axis1=1,axis2=2)/3.0
            self.stress       = stress
            self.stress_units = 'kbar'
            self.pressures    = pressures
            self.pressure     = pressures[-1]
    #end def read_stress


    def read_convergence(self,lines):
        """Parse electronic and ionic convergence messages.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``convergence`` to an ``obj`` containing nullable Boolean
        convergence indicators and integer success and failure counts.
        """
        electronic_successes = 0
        electronic_failures  = 0
        ionic_converged      = None
        # Match and classify electronic or ionic convergence messages.
        # Example: Potential convergence has been achieved. stopping ...
        convergence_pattern = re.compile(
            r'(?P<electronic_failure>'
            r'potential\s+convergence.*\bnot\s+(?:been\s+)?achieved\b|'
            r'convergence\s+criterion.*\bnot\s+met\b)|'
            r'(?P<electronic_success>'
            r'potential\s+convergence.*\b(?:has\s+been\s+)?achieved\b)|'
            r'(?P<ionic_failure>'
            r'force\s+convergence.*\bnot\s+(?:been\s+)?achieved\b)|'
            r'(?P<ionic_success>'
            r'force\s+convergence.*\b(?:has\s+been\s+)?achieved\b)',
            re.IGNORECASE,
            )
        for line in lines:
            match = convergence_pattern.search(line)
            if match is None:
                continue
            if match['electronic_failure'] is not None:
                electronic_failures += 1
            elif match['electronic_success'] is not None:
                electronic_successes += 1
            elif match['ionic_failure'] is not None:
                ionic_converged = False
            elif match['ionic_success'] is not None:
                ionic_converged = True
        electronic_converged = None
        if electronic_successes+electronic_failures>0:
            electronic_converged = electronic_successes>0 and electronic_failures==0
        if electronic_successes+electronic_failures+(ionic_converged is not None)>0:
            self.convergence = obj(
                electronic_converged = electronic_converged,
                electronic_successes = electronic_successes,
                electronic_failures  = electronic_failures,
                ionic_converged      = ionic_converged,
                )
    #end def read_convergence


    def read_timing(self,lines):
        """Parse the RMG timing summary and its named sections.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``timing`` to an ``obj`` containing total and per-step seconds.
        Its ``sections`` member is an ``obj`` mapping normalized section names to
        nested ``obj`` instances with ``total`` and ``per_step`` values.
        """
        timing   = None
        sections = obj()
        # Match a finite, infinite, or NaN numeric token in the timing table.
        # Example: 3.00
        time_value = r'(?:'+self.number_pattern+r'|inf|nan)'
        # Match a numbered timing section followed by total and per-step times.
        # Example: 1-TOTAL 3.00 0.50
        pattern    = re.compile(
            r'^\s*(\d+\s*-\s*.*?)\s+('+time_value+r')\s+('+time_value+r')'
            r'(?:\s+.*)?$',re.IGNORECASE)
        for line in lines:
            match = pattern.match(line)
            if match is None:
                continue
            # Match spacing around the first dash so the section name can be normalized.
            # Example: 1 - TOTAL
            name     = re.sub(r'\s*-\s*','-',match.group(1).strip(),count=1)
            total    = float(match.group(2))
            per_step = float(match.group(3))
            # Match non-alphanumeric runs so a timing name can become an object key.
            # Example: 1-TOTAL
            key             = re.sub(r'[^a-z0-9]+','_',name.lower()).strip('_')
            sections[key]   = obj(total=total,per_step=per_step)
            if name.lower()=='1-total':
                timing = obj(total=total,per_step=per_step,units='s')
        if timing is not None:
            timing.sections = sections
            self.timing     = timing
    #end def read_timing


    def read_produced_files(self):
        """Locate files produced by QMCPACK-interface SCF runs.

        Notes
        -----
        Binds ``produced_files`` to an ``obj`` containing the QMCPACK restart
        filepath when present.
        """
        produced_files = obj()
        if 'files' in self.setup_info:
            data_output = self.setup_info.files.get('data_output_file',None)
            if data_output is not None:
                qmcpack_file = os.path.join(self.path,data_output+'.h5')
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
    arg0 : Simulation or str or None, optional
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
        RMG input reconstructed from the output when it can be obtained.
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
        and ``'atm'``.
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
        If ``arg0`` is neither a ``Simulation``, a string, nor ``None``.
    FileNotFoundError
        If a supplied output path does not exist.
    IsADirectoryError
        If a supplied output path does not identify a regular file.
    """

    all_modes         = frozenset({'scf','nscf','relax'})
    relaxation_modes  = frozenset({'relax'})
    pressure_units    = MappingProxyType({
        'Pa'   : 1.0,
        'bar'  : 1e5,
        'kbar' : 1e8,
        'Mbar' : 1e11,
        'GPa'  : 1e9,
        'atm'  : 1.01325e5,
        })


    def _require_supported(self,quantity,modes):
        """Require analyzed output and a run mode supporting the quantity."""
        if self.results is None:
            raise RuntimeError(
                f'RMG quantity "{quantity}" is unavailable because output has not been analyzed'
                )
        if self.run_mode not in modes:
            raise RuntimeError(
                f'RMG quantity "{quantity}" is not supported for run mode "{self.run_mode}"'
                )
    #end def _require_supported


    def initial_structure(self,units='A'):
        """Return the input ``Structure`` in Angstrom or bohr."""
        self._require_supported('initial_structure',self.all_modes)
        if units not in {'A','B'}:
            raise ValueError('initial_structure units must be one of: A, B')
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
            raise ValueError('energy units must be one of: eV, Ha, Ry')
        value = self.results.energy
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
            raise ValueError('kpoints units must be one of: A, B')
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
            raise ValueError('eigenvalues units must be one of: eV, Ha, Ry')
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
            raise ValueError('Ef units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.fermi_energies)==0:
            return None
        return convert(electronic.fermi_energies[-1],'eV',units)
    #end def Ef


    def Evbm(self,units='eV'):
        """Return the final valence-band maximum in selected energy units."""
        self._require_supported('Evbm',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Evbm units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.valence_band_maxima)==0:
            return None
        return convert(electronic.valence_band_maxima[-1],'eV',units)
    #end def Evbm


    def Ecbm(self,units='eV'):
        """Return the final conduction-band minimum in selected energy units."""
        self._require_supported('Ecbm',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Ecbm units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.conduction_band_minima)==0:
            return None
        return convert(electronic.conduction_band_minima[-1],'eV',units)
    #end def Ecbm


    def band_gap(self,units='eV'):
        """Return the final band gap in eV, Hartree, or Rydberg."""
        self._require_supported('band_gap',self.all_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('band_gap units must be one of: eV, Ha, Ry')
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
            occupations,full_occupation,rtol=0.0,atol=tolerance)
        return bool(np.any(~(empty|full)))
    #end def fractional_occs


    def relaxed_structure(self,units='A'):
        """Return the final relaxed ``Structure`` in Angstrom or bohr."""
        self._require_supported('relaxed_structure',self.relaxation_modes)
        if units not in {'A','B'}:
            raise ValueError('relaxed_structure units must be one of: A, B')
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
            raise ValueError('forces units must be one of: eV/A, Ry/B, Ha/B')
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
            raise ValueError(f'stress units must be one of: {supported}')
        stress = self.results.stress
        if stress is None:
            return None
        return stress*1e8/self.pressure_units[units]
    #end def stress


    def pressure(self,units='GPa'):
        """Return the final hydrostatic pressure in selected pressure units."""
        self._require_supported('pressure',self.all_modes)
        if units not in self.pressure_units:
            supported = ', '.join(sorted(self.pressure_units))
            raise ValueError(f'pressure units must be one of: {supported}')
        pressure = self.results.pressure
        if pressure is None:
            return None
        return pressure*1e8/self.pressure_units[units]
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
            if not isinstance(arg0,str):
                raise TypeError(
                    'invalid type provided for log_file\n'
                    'Type expected: str\n'
                    f'Type provided: {arg0.__class__.__name__}'
                    )
            elif not os.path.exists(arg0):
                raise FileNotFoundError(
                    'RMG log output file does not exist.\n'
                    f'Path provided: {arg0}'
                    )
            elif not os.path.isfile(arg0):
                raise IsADirectoryError(
                    'Path provided for RMG log output is not a file.\n'
                    f'Path provided: {arg0}'
                    )
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
        self.input    = results.input
    #end def analyze

#end class RmgAnalyzer
