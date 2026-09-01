##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os
import re
from copy import deepcopy
from glob import glob
from types import MappingProxyType

import numpy as np

from .developer import DevBase,obj
from .generic import NexusError
from .unit_converter import convert
from .numerics import simstats
from .simulation import SimulationAnalyzer, Simulation
from .structure import generate_structure
from .rmg_input import RmgInput, rmg_modes


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
        Short RMG calculation mode, such as ``"scf"``, ``"band"``, or
        ``"md_VE"``.
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
    bands : obj or None
        Mapping from spin index to an ``obj`` containing band-path distances,
        a two-dimensional energy array, units, and the source file path.
    md : obj or None
        Molecular-dynamics step, energy, temperature, and displacement arrays.
    md_stats : obj or None
        Mapping from molecular-dynamics quantity to an ``obj`` containing its
        mean, variance, statistical error, and autocorrelation factor.
    tddft : obj or None
        Time-dependent energy and spin-resolved dipole series represented by
        nested ``obj`` instances and NumPy arrays.
    neb : obj or None
        NEB controller settings, ordered input structures, path energies,
        relative energies, barriers, and image-local results. Path energies
        and barriers are in Hartree and the reaction coordinate is in bohr.
    produced_files : obj or None
        Lists or paths for mode-specific output files, such as exact-exchange
        integrals, STM data, or a QMCPACK restart file.

    Notes
    -----
    ``geometry``, ``convergence``, and ``timing`` apply to the supported modes
    ``scf``, ``nscf``, ``band``, ``exx``, ``relax``, ``md_VE``, ``md_TE``,
    ``tddft``, ``stm``, and ``neb``. Energy, SCF, ionic, and stress members
    apply to ``scf``, ``nscf``, ``relax``, ``md_VE``, ``md_TE``, ``tddft``,
    and ``neb``. ``electronic`` also applies to ``band``. ``bands`` applies
    only to ``band``; ``md`` and ``md_stats`` apply to ``md_VE`` and
    ``md_TE``; ``tddft`` applies only to ``tddft``; ``neb`` applies only to
    ``neb``; and ``produced_files`` applies to ``scf``, ``exx``, and
    ``stm``. A mode-applicable member is initialized to ``None`` and remains
    ``None`` when its data cannot be obtained. Members that do not apply to
    the detected mode are not bound.

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
            msg = (
                'invalid type provided for filepath\n'
                'Type expected: str\n'
                f'Type provided: {filepath.__class__.__name__}'
                )
            raise TypeError(msg)
        elif not os.path.exists(filepath):
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
        self.input         = None
        self.setup_info    = None

        with open(filepath,'r') as fobj:
            lines = fobj.read().splitlines()
        self.read_setup_info(lines)

        electronic_modes = {'scf','nscf','relax','md_VE','md_TE','tddft','neb'}
        eigenvalue_modes = electronic_modes|{'band'}
        supported_modes  = eigenvalue_modes|{'exx','stm'}

        # modes: scf, nscf, band, exx, relax, md_VE, md_TE, tddft, stm, neb
        if self.run_mode in supported_modes:
            self.geometry    = None
            self.convergence = None
            self.timing      = None

            self.read_geometry()
            self.read_convergence(lines)
            self.read_timing(lines)

        # modes: scf, nscf, relax, md_VE, md_TE, tddft, neb
        if self.run_mode in electronic_modes:
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

            self.read_energies(lines)
            self.read_scf(lines)
            self.read_ions(lines)
            self.read_stress(lines)

        # modes: scf, nscf, band, relax, md_VE, md_TE, tddft, neb
        if self.run_mode in eigenvalue_modes:
            self.electronic = None

            self.read_electronic(lines)

        # modes: md_VE, md_TE
        if self.run_mode in {'md_VE','md_TE'}:
            self.md       = None
            self.md_stats = None

            self.read_md(lines)

        # modes: band
        if self.run_mode=='band':
            self.bands = None

            self.read_band()

        # modes: tddft
        if self.run_mode=='tddft':
            self.tddft = None

            self.read_tddft()

        # modes: neb
        if self.run_mode=='neb':
            self.neb = None

            self.read_neb(lines)

        # modes: scf, exx, stm
        if self.run_mode in {'scf','exx','stm'}:
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
            return ' '.join(line.expandtabs().strip().split()).lower()
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
            # Example: FERMI ENERGY : 5.25 eV
            match = re.search(pattern,line,re.I)
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
                # Match the band-structure calculation label.
                # Example: Band structure
                (r'band\s+structure','band'),
                # Match a calculation of exact-exchange integrals.
                # Example: Calculate Exx integral's from saved wave functions
                (r'exx\s+integral','exx'),
                # Match either supported structural-relaxation label.
                # Example: Structure optimization
                (r'(?:structure\s+optimization|relax\s+structure)','relax'),
                # Match constant-volume, constant-energy molecular dynamics.
                # Example: Molecular dynamics - CVE
                (r'(?:molecular\s+dynamics\s*[-:]?\s*cve|constant\s+volume\s+and\s+energy)','md_VE'),
                # Match constant-volume, constant-temperature molecular dynamics.
                # Example: Molecular dynamics - CVT
                (r'(?:molecular\s+dynamics\s*[-:]?\s*cvt|constant\s+temperature\s+and\s+energy)','md_TE'),
                # Match either supported nudged-elastic-band label.
                # Example: Molecular dynamics using Nudged Elastic Band.
                (r'(?:nudged\s+elastic\s+band|neb\s+relax)','neb'),
                # Match either full or abbreviated time-dependent DFT labels.
                # Example: Time dependent DFT
                (r'(?:time\s+dependent\s+dft|\btddft\b)','tddft'),
                # Match either full or abbreviated scanning-tunneling labels.
                # Example: Calculate STM charge density
                (r'stm\s+charge\s+density|\bstm\b','stm'),
                )
            for pattern,run_mode in mode_patterns:
                if re.search(pattern,text,re.I):
                    return run_mode
            return rmg_modes.mode_match(text,short=True)
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
            s      = re.sub(r'\([^)]*\)','',s)
            tokens = s.strip().lower().split()
            name   = ''
            for t in tokens:
                if not t.startswith('('):
                    name += t+'_'
            name = name[:-1].replace('/','_').replace('-','_')
            return name
        #end def process_name
        def process_value(v,*,list=False):
            """Convert setup text to a scalar, array, list, or string and units."""
            v     = v.strip()
            units = None
            try:
                v = int(v)
            except (ValueError,TypeError,OverflowError):
                try:
                    v = float(v)
                except (ValueError,TypeError,OverflowError):
                    if ' ' in v or ',' in v:
                        vt = v.replace(',',' ')
                        if len(vt)>0:
                            tokens = vt.split()
                            if tokens[-1] in unit_set:
                                units  = tokens[-1]
                                tokens = tokens[:-1]
                            try:
                                if not list:
                                    v = np.array(tokens,dtype=float)
                                else:
                                    v = [process_value(t,list=True)[0] for t in tokens]
                            except (ValueError,TypeError,OverflowError):
                                units = None
                    elif v in on_off:
                        v = on_off[v]
            return v,units
        #end def process_value
        if mode is not None:
            # Match the standalone Files setup-section heading.
            # Example: Files:
            start_match = re.search(r'(?im)^\s*files\s*(?::\s*)?$',log_text)
            istart      = start_match.start() if start_match is not None else -1
            if istart!=-1:
                # Match the normal heading that terminates the setup section.
                # Example: Initial Ionic Positions And Displacements (Angstrom)
                end_pattern = (
                    r'(?im)^\s*initial\s+ionic\s+positions\s+and\s+'
                    r'displacements\s*\(\s*angstroms?\s*\)'
                    )
                end_match   = re.search(end_pattern,log_text[istart:])
                if end_match is not None:
                    iend = istart+end_match.start()
                else:
                    # Match alternate post-setup headings when ionic positions are absent.
                    # Example: Davidson converged
                    fallback_pattern = (
                        r'(?im)^\s*(?:diagonalization\s+using|davidson\s+'
                        r'(?:converged|incomplete)|-+\s*timing\s+information)'
                        )
                    end_match        = re.search(fallback_pattern,log_text[istart:])
                    iend             = istart+end_match.start() if end_match is not None else len(log_text)
                if iend!=-1:
                    text        = log_text[istart:iend].expandtabs()
                    blocks      = []
                    b           = ''
                    last_header = False
                    for line in text.splitlines():
                        if len(line)>0:
                            if line[0]!=' ':
                                if last_header:
                                    b += '\n'+line
                                else:
                                    if len(b)>0:
                                        blocks.append(b)
                                    b           = line
                                    last_header = True
                            else:
                                b += '\n'+line
                                last_header = False
                    if len(b)>0:
                        blocks.append(b)
                    other_blocks = obj()
                    for b in blocks:
                        if '\n' not in b:
                            continue
                        header,body  = b.split('\n',1)
                        bname         = process_name(header)
                        lines         = body.splitlines()
                        value_lines   = [line for line in lines if ':' in line]
                        table_blocks  = {
                            'k_points',
                            'initial_ionic_positions_and_displacements',
                            }
                        simple_values = bname not in table_blocks and len(value_lines)>0
                        if simple_values:
                            bvalues = obj()
                            for line in value_lines:
                                name,value    = line.split(':',1)
                                name          = process_name(name)
                                value,units   = process_value(value)
                                bvalues[name] = value
                                if units is not None:
                                    bvalues.units = units
                            setup_info[bname] = bvalues
                        else:
                            other_blocks[bname] = header,body,lines
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
                                r'\s*([^\s,;]+)',s,re.I)
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
                        header,body,lines = other_blocks.k_points
                        del other_blocks.k_points
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
                        header,body,lines = other_blocks[k]
                        del other_blocks[k]
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
                                ls = line.strip()
                                if len(ls)>0:
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
            for filepath in filepaths:
                if os.path.isfile(filepath):
                    break
            if os.path.isfile(filepath):
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
        # Match the final total energy obtained from the eigenvalue sum.
        # Example: final total energy from eigenvalue sum : -1.2345 Ha
        epat = re.compile(
            r'final\s+total\s+energy\s+from\s+eig(?:envalue)?\s+sum\s*[:=]\s*'
            r'('+self.number_pattern+r')(?:\s+([^\s,;]+))?',re.I)
        # Match the final total energy reported by the direct calculation.
        # Example: final total energy from direct = -1.2300 Ha
        dpat = re.compile(
            r'final\s+total\s+energy\s+from\s+direct\s*[:=]\s*'
            r'('+self.number_pattern+r')(?:\s+([^\s,;]+))?',re.I)
        for line in lines:
            match = epat.search(line)
            if match is not None:
                energies.append(float(match.group(1).replace('D','E').replace('d','e')))
                energy_units.append(match.group(2))
            match = dpat.search(line)
            if match is not None:
                direct_energies.append(float(match.group(1).replace('D','E').replace('d','e')))
                direct_energy_units.append(match.group(2))
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

        def match_float(pattern,line):
            """Return the first floating-point value captured by a regex."""
            # Apply a caller-provided expression that captures one physical value.
            # Example: FERMI ENERGY : 5.25 eV
            match = re.search(pattern,line,re.I)
            if match is None:
                return None
            return float(match.group(1).replace('D','E').replace('d','e'))
        #end def match_float

        def read_eigenvalue_data():
            """Return the final complete k-point eigenvalue and occupation table."""
            npat = self.number_pattern
            # Match an eigenvalue-table heading with its k-point index and coordinates.
            # Example: KOHN SHAM EIGENVALUES [eV] AT K-POINT [ 1]: 0.0 0.5 0.0
            header_pattern = re.compile(
                r'KOHN\s+SHAM\s+EIGENVALUES\s*\[\s*eV\s*\]\s+AT\s+'
                r'K-POINT\s*\[\s*(\d+)\s*\]\s*:\s*'
                r'('+npat+r')\s+('+npat+r')\s+('+npat+r')',re.I)
            # Match an eigenvalue-table row and capture its k-point index and values.
            # Example: [kpt 1 0 0] -6.4 [2.000] -1.9 [0.000]
            row_pattern = re.compile(
                r'^\s*\[\s*kpt\s+(\d+)\s+[-+]?\d+\s+\d+\s*\]\s*(.*)$',re.I)
            # Match one eigenvalue followed by its bracketed occupation.
            # Example: -6.4238 [2.000]
            pair_pattern = re.compile(
                r'('+npat+r')\s*\[\s*('+npat+r')\s*\]',re.I)

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
                        for i in sorted({2,3,4})
                        ],dtype=float)
                    dataset[index] = dict(
                        kpoint   = coordinates,
                        channels = {},
                        )
                    kpoint = index
                    spin   = 'none'
                    continue
                if kpoint is None:
                    continue
                normalized = ' '.join(line.expandtabs().strip().split()).lower()
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
                if spin not in channels:
                    channels[spin] = [[],[]]
                eigs,occs = channels[spin]
                for eigenvalue,occupation in pairs:
                    eigs.append(float(
                        eigenvalue.replace('D','E').replace('d','e')))
                    occs.append(float(
                        occupation.replace('D','E').replace('d','e')))
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
                nstates = None
                valid   = True
                for index in indices:
                    channels = dataset[index]['channels']
                    if set(channels)!=set(spins):
                        valid = False
                        break
                    for spin in spins:
                        eigs,occs = channels[spin]
                        if len(eigs)==0 or len(eigs)!=len(occs):
                            valid = False
                            break
                        if nstates is None:
                            nstates = len(eigs)
                        elif len(eigs)!=nstates:
                            valid = False
                            break
                    if not valid:
                        break
                if not valid:
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
                return obj(
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
        for line in lines:
            lower = ' '.join(line.expandtabs().strip().split()).lower()
            if 'fermi energy' in lower:
                # Match and capture the Fermi energy.
                # Example: FERMI ENERGY : 5.25 eV
                value = match_float(
                    r'fermi\s+energy\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.fermi_energies.append(value)
            elif 'valence band maximum' in lower and 'conduction band' in lower:
                # Match and capture the valence-band maximum.
                # Example: spin0: valence band maximum = 4.0 eV
                vbm = match_float(
                    r'valence\s+band\s+maximum\s*[:=]\s*('+self.number_pattern+r')',line)
                # Match and capture the conduction-band minimum.
                # Example: conduction band minimum = 6.0 eV
                cbm = match_float(
                    r'conduction\s+band\s+(?:minimum|minumm)\s*[:=]\s*'
                    r'('+self.number_pattern+r')',line)
                if vbm is not None and cbm is not None:
                    data.valence_band_maxima.append(vbm)
                    data.conduction_band_minima.append(cbm)
            elif 'band gap' in lower:
                # Match and capture the electronic band gap.
                # Example: spin0: Band gap : 2.0 eV
                value = match_float(
                    r'band\s+gap\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.band_gaps.append(value)
            elif 'total charge in supercell' in lower:
                # Match and capture the total charge in the simulation cell.
                # Example: Total charge in supercell = 0.0
                value = match_float(
                    r'total\s+charge\s+in\s+supercell\s*[:=]\s*'
                    r'('+self.number_pattern+r')',line)
                if value is not None:
                    data.total_charges.append(value)
            elif '@@ total magnetization' in lower:
                # Match and capture the signed total magnetization.
                # Example: @@ TOTAL MAGNETIZATION = 1.0
                value = match_float(
                    r'total\s+magnetization\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.total_magnetizations.append(value)
            elif '@@ absolute magnetization' in lower:
                # Match and capture the absolute total magnetization.
                # Example: @@ ABSOLUTE MAGNETIZATION = 1.5
                value = match_float(
                    r'absolute\s+magnetization\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.absolute_magnetizations.append(value)
            elif lower.lstrip().startswith('sum force'):
                values = line_numbers(line.split('=',1)[-1])
                if len(values)>=3:
                    sum_forces.append(values[:3])
            elif 'volume and energy per atom' in lower:
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

        nfound                     = sum(
            len(v) for v in data.values() if isinstance(v,np.ndarray))
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
        component_names = obj(
            eigenvalue_sum  = 'EIGENVALUE SUM',
            ion_ion         = 'ION_ION',
            electrostatic   = 'ELECTROSTATIC',
            vxc             = 'VXC',
            exc             = 'EXC',
            total_energy    = 'TOTAL ENERGY',
            estimated_error = 'estimated error',
            )

        values = obj()
        for name in component_names.keys():
            values[name] = []
        md_steps   = []
        scf_steps  = []
        step_times = []
        scf_times  = []
        rms_dv     = []
        for line in lines:
            for name,label in component_names.items():
                # Match an @@ energy component and capture its value or overflow stars.
                # Example: @@ TOTAL ENERGY : -1.250000
                pattern = (
                    r'^\s*@@\s*'+re.escape(label)+r'\s*[:=]\s*'
                    r'('+self.number_pattern+r'|\*+)'
                    )
                match   = re.search(pattern,line,re.I)
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
            if re.search(r'\bquench\s*:',line,re.I):
                # Match and capture the current molecular-dynamics step.
                # Example: md: 0/2
                md_match   = re.search(r'\bmd\s*:\s*(\d+)\s*/',line,re.I)
                # Match and capture the current SCF iteration.
                # Example: scf: 3/20
                scf_match  = re.search(r'\bscf\s*:\s*(\d+)\s*/',line,re.I)
                # Match and capture the elapsed time for the current step.
                # Example: step time: 0.20
                step_match = re.search(
                    r'\bstep\s+time\s*:\s*('+self.number_pattern+r')',line,re.I)
                # Match and capture the elapsed SCF time.
                # Example: scf time: 0.80
                time_match = re.search(
                    r'\bscf\s+time\s*:\s*('+self.number_pattern+r')',line,re.I)
                # Match and capture the RMS potential change inside brackets.
                # Example: RMS [ dV ] : 2.0D-5
                rms_match  = re.search(r'\brms\s*\[\s*dv\s*\]\s*:\s*([^\]\s]+)',line,re.I)
                md_steps.append(int(md_match.group(1)) if md_match is not None else -1)
                scf_steps.append(int(scf_match.group(1)) if scf_match is not None else -1)
                step_times.append(
                    float(step_match.group(1)) if step_match is not None else np.nan)
                scf_times.append(
                    float(time_match.group(1)) if time_match is not None else np.nan)
                rms = rms_match.group(1) if rms_match is not None else '*'
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
                atom          = tokens[2]
                position      = values[:3]
                charge        = values[3]
                magnetization = values[4]
                force         = values[5:8]
                move          = [int(v) for v in move_values]
                atoms.append(atom)
                positions.append(position)
                charges.append(charge)
                magnetizations.append(magnetization)
                forces.append(force)
                movable.append(move)
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
        ionic_steps = obj()
        for n,record in enumerate(records):
            ionic_steps[n] = record
        if len(records)>0:
            self.ionic_steps    = ionic_steps
            self.position_units = 'a0'
            self.force_units    = 'Ha/a0'
        if len(records)>0 and len(set(len(r.atoms) for r in records))==1:
            positions           = np.array([r.positions for r in records],dtype=float)
            forces              = np.array([r.forces for r in records],dtype=float)
            charges             = np.array([r.charges for r in records],dtype=float)
            magnetizations      = np.array([r.magnetizations for r in records],dtype=float)
            max_forces          = np.array([
                np.sqrt((r.forces**2).sum(axis=1)).max() for r in records
                ],dtype=float)

            self.positions      = positions
            self.forces         = forces
            self.charges        = charges
            self.magnetizations = magnetizations
            self.max_forces     = max_forces
            if 'structure' in self.setup_info:
                structures = obj()
                initial    = self.setup_info.structure
                for n,record in enumerate(records):
                    structure = generate_structure(
                        units = 'B',
                        axes  = initial.axes,
                        elem  = record.atoms,
                        pos   = record.positions,
                        )
                    structures[n] = structure
                self.structures = structures
    #end def read_ions


    def read_md(self,lines):
        """Parse molecular-dynamics records and compute summary statistics.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines.

        Notes
        -----
        Binds ``md`` to an ``obj`` of NumPy arrays for step, energy,
        temperature, and displacement data. ``md_stats`` is an ``obj`` mapping
        each numerical quantity to an ``obj`` containing ``mean``, ``var``,
        ``error``, and ``kappa``.
        """
        records = []
        for line in lines:
            stripped = line.strip()
            # Match constant-energy or constant-temperature MD record markers.
            # Example: @CVE 1 -1.0 0.1 -0.9 300.0 2.5e-4
            marker_match = re.match(r'^@(CVE|CVT)(?:\b|-)',stripped,re.I)
            if marker_match is None:
                continue
            fields = stripped.split(None,1)
            if len(fields)<2:
                continue
            # Match a whitespace-separated numeric sequence at the start of a line.
            # Example: 1 -1.0 0.1 -0.9 300.0 2.5e-4
            npat        = self.number_pattern
            value_match = re.match(
                r'^\s*((?:'+npat+r')(?:\s+(?:'+npat+r'))*)',fields[1])
            if value_match is None:
                continue
            # Find every RMG-formatted number in the matched numeric sequence.
            # Example: 1 -1.0 0.1 -0.9 300.0 2.5e-4
            values = re.findall(
                self.number_pattern,
                value_match.group(1).replace('D','E').replace('d','e'))
            values = np.array(values,dtype=float)
            if len(values)>=6:
                records.append(obj(
                    step             = int(values[0]),
                    potential_energy = values[1],
                    kinetic_energy   = values[2],
                    total_energy     = values[3],
                    temperature      = values[4],
                    displacement     = values[5],
                    ))
        if len(records)>0:
            md = obj()
            for name in records[0].keys():
                dtype    = int if name=='step' else float
                md[name] = np.array([r[name] for r in records],dtype=dtype)
            md.energy_units      = 'Ha'
            md.temperature_units = 'K'
            self.md              = md
            statistics           = obj()
            for name,values in md.items():
                if not isinstance(values,np.ndarray):
                    continue
                mean,var,error,kappa = simstats(values)
                statistics[name] = obj(
                    mean  = mean,
                    var   = var,
                    error = error,
                    kappa = kappa,
                    )
            self.md_stats = statistics
    #end def read_md


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
        elif (
            self.run_mode=='band' and
            isinstance(self.input,RmgInput) and
            'kpoints_bandstructure' in self.input
            ):
            band_path = self.input.kpoints_bandstructure
            endpoints = np.asarray(band_path.kpoints,dtype=float)
            counts    = np.asarray(band_path.counts,dtype=int)
            if (
                endpoints.ndim==2 and endpoints.shape[1:]==(3,) and
                len(endpoints)==len(counts) and len(endpoints)>0 and
                np.all(counts>=0)
                ):
                kpoints = [endpoints[0]]
                for i in range(1,len(endpoints)):
                    count = counts[i]
                    if count==0:
                        continue
                    start = endpoints[i-1]
                    stop  = endpoints[i]
                    for n in range(1,count+1):
                        kpoints.append(start+(stop-start)*n/count)
                geometry.kpoints_crystal = np.array(kpoints,dtype=float)
                if 'structure' in self.setup_info:
                    geometry.kpoints_cart = np.dot(
                        geometry.kpoints_crystal,self.setup_info.structure.kaxes)
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
            normalized = ' '.join(line.expandtabs().strip().split()).lower()
            # Match the heading for a total stress tensor reported in kbar.
            # Example: stress total in unit of kbar
            if not re.search(r'\bstress\s+total\b',normalized) or 'kbar' not in normalized:
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
                tensors.append(np.array(rows,dtype=float))
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
        for line in lines:
            lower = ' '.join(line.expandtabs().strip().split()).lower()
            # Match an explicit electronic-convergence failure message.
            # Example: Potential convergence has not been achieved
            if re.search(r'potential\s+convergence.*\bnot\s+(?:been\s+)?achieved\b',lower):
                electronic_failures += 1
            # Match an explicit electronic-convergence success message.
            # Example: Potential convergence has been achieved. stopping ...
            elif re.search(r'potential\s+convergence.*\b(?:has\s+been\s+)?achieved\b',lower):
                electronic_successes += 1
            # Match the alternate electronic-convergence failure wording.
            # Example: Convergence criterion not met
            elif re.search(r'convergence\s+criterion.*\bnot\s+met\b',lower):
                electronic_failures += 1
            # Match an explicit ionic force-convergence failure message.
            # Example: Force convergence has not been achieved
            elif re.search(r'force\s+convergence.*\bnot\s+(?:been\s+)?achieved\b',lower):
                ionic_converged = False
            # Match an explicit ionic force-convergence success message.
            # Example: Force convergence has been achieved
            elif re.search(r'force\s+convergence.*\b(?:has\s+been\s+)?achieved\b',lower):
                ionic_converged = True
        electronic_converged = None
        if electronic_successes+electronic_failures>0:
            electronic_converged = electronic_successes>0 and electronic_failures==0
        convergence = obj(
            electronic_converged = electronic_converged,
            electronic_successes = electronic_successes,
            electronic_failures  = electronic_failures,
            ionic_converged      = ionic_converged,
            )

        nfound = electronic_successes+electronic_failures+(ionic_converged is not None)
        if nfound>0:
            self.convergence = convergence
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
            r'(?:\s+.*)?$',re.I)
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
            normalized_name = ' '.join(name.expandtabs().strip().split()).lower()
            if normalized_name=='1-total':
                timing = obj(total=total,per_step=per_step,units='s')
        if timing is not None:
            timing.sections = sections
            self.timing     = timing
    #end def read_timing


    def read_band(self):
        """Read spin-resolved band structures from companion data files.

        Notes
        -----
        Binds ``bands`` to an ``obj`` mapping integer spin indices to nested
        ``obj`` instances. Each contains a path-distance array, a two-dimensional
        band-energy array in eV, and its source filepath.
        """
        def leading_numbers(line):
            """Return the leading whitespace-separated numbers as a float array."""
            # Match a whitespace-separated numeric sequence at the start of a line.
            # Example: 0.0 -1.0
            npat    = self.number_pattern
            pattern = r'^\s*((?:'+npat+r')(?:\s+(?:'+npat+r'))*)'
            match   = re.match(pattern,line)
            if match is None:
                return np.array([],dtype=float)
            # Find every RMG-formatted number in the matched numeric sequence.
            # Example: 0.0 -1.0
            values = re.findall(
                self.number_pattern,
                match.group(1).replace('D','E').replace('d','e'))
            return np.array(values,dtype=float)
        #end def leading_numbers

        prefix = self.outfile_name[:-4] if self.outfile_name.endswith('.log') else self.outfile_name
        files  = sorted(glob(os.path.join(self.path,prefix+'_spin*.bandstructure.dat')))
        bands  = obj()
        for filepath in files:
            # Match and capture the spin index in a band-structure filename.
            # Example: input_spin0.bandstructure.dat
            match = re.search(r'_spin(\d+)\.bandstructure\.dat$',filepath)
            if match is None:
                continue
            groups = []
            group  = []
            with open(filepath,'r') as fobj:
                for line in fobj:
                    if '&&' in line:
                        if len(group)>0:
                            groups.append(np.array(group,dtype=float))
                            group = []
                    else:
                        values = leading_numbers(line)
                        if len(values)>=2:
                            group.append(values[:2])
            if len(group)>0:
                groups.append(np.array(group,dtype=float))
            if len(groups)>0 and len(set(len(g) for g in groups))==1:
                spin = int(match.group(1))
                bands[spin] = obj(
                    distance     = groups[0][:,0],
                    energies     = np.array([g[:,1] for g in groups],dtype=float),
                    energy_units = 'eV',
                    filepath     = filepath,
                    )
        if len(bands)>0:
            self.bands = bands
    #end def read_band


    def read_tddft(self):
        """Read TDDFT energy and spin-resolved dipole time series.

        Notes
        -----
        Binds ``tddft`` to an ``obj``. Its optional ``energy`` member contains
        NumPy arrays for time and energy-component changes. Its optional
        ``dipoles`` member maps spin indices to nested ``obj`` instances holding
        time and dipole arrays, electric-field and ground-state vectors, and the
        source filepath.
        """
        def line_numbers(line):
            """Return all RMG-formatted numbers in a line as a float array."""
            # Find every RMG-formatted number occurring anywhere in a line.
            # Example: &&electric field: 0.0 0.0 0.1
            values = re.findall(
                self.number_pattern,line.replace('D','E').replace('d','e'))
            return np.array(values,dtype=float)
        #end def line_numbers

        def leading_numbers(line):
            """Return the leading whitespace-separated numbers as a float array."""
            # Match a whitespace-separated numeric sequence at the start of a line.
            # Example: 0.0 1.1, 2.1, 3.1
            line    = line.replace(',',' ')
            npat    = self.number_pattern
            pattern = r'^\s*((?:'+npat+r')(?:\s+(?:'+npat+r'))*)'
            match   = re.match(pattern,line)
            if match is None:
                return np.array([],dtype=float)
            return line_numbers(match.group(1))
        #end def leading_numbers

        prefix       = self.outfile_name[:-4] if self.outfile_name.endswith('.log') else self.outfile_name
        energy_file  = os.path.join(self.path,prefix+'_totalE')
        dipole_files = sorted(glob(os.path.join(self.path,prefix+'_spin*_dipole.dat')))
        tddft        = obj()
        if os.path.isfile(energy_file):
            rows = []
            with open(energy_file,'r') as fobj:
                for line in fobj:
                    if line.lstrip().startswith('&&'):
                        continue
                    values = leading_numbers(line)
                    if len(values)>=5:
                        rows.append(values[:5])
            if len(rows)>0:
                values = np.array(rows,dtype=float)

                tddft.energy = obj(
                    time                  = values[:,0],
                    kinetic_pseudo_change = values[:,1],
                    hartree_change        = values[:,2],
                    xc_change             = values[:,3],
                    total_energy_change   = values[:,4],
                    energy_units          = 'Ha',
                    filepath              = energy_file,
                    )
        dipoles = obj()
        for filepath in dipole_files:
            rows         = []
            field        = None
            ground_state = None
            with open(filepath,'r') as fobj:
                for line in fobj:
                    lower  = ' '.join(line.expandtabs().strip().split()).lower()
                    values = line_numbers(line)
                    if 'electric field' in lower and len(values)>=3:
                        values = line_numbers(line.split(':',1)[-1])
                        if len(values)>=3:
                            field = values[:3]
                    elif 'dipole at' in lower and len(values)>=3:
                        values = line_numbers(line.split(':',1)[-1])
                        if len(values)>=3:
                            ground_state = values[:3]
                    elif not line.lstrip().startswith('&&'):
                        values = leading_numbers(line)
                        if len(values)>=4:
                            rows.append(values[:4])
            if len(rows)>0:
                # Match and capture the spin index in a dipole filename.
                # Example: input_spin0_dipole.dat
                match  = re.search(r'_spin(\d+)_dipole\.dat$',filepath)
                spin   = int(match.group(1)) if match is not None else len(dipoles)
                values = np.array(rows,dtype=float)
                dipoles[spin] = obj(
                    time           = values[:,0],
                    dipole         = values[:,1:4],
                    electric_field = field,
                    ground_state   = ground_state,
                    filepath       = filepath,
                    )
        if len(dipoles)>0:
            tddft.dipoles = dipoles
        if len(tddft)>0:
            self.tddft = tddft
    #end def read_tddft


    def read_neb(self,lines):
        """Read NEB controller, path, energy-profile, and local-image data.

        Parameters
        ----------
        lines : list of str
            Complete RMG log split into lines. The NEB controller and sibling
            image inputs and logs are read relative to the current image log.

        Notes
        -----
        Binds ``neb`` to an ``obj`` containing controller settings and file
        paths, ordered input ``Structure`` objects, a cumulative input-path
        reaction coordinate, image energies and relative energies, forward and
        reverse barriers, the highest-energy image index, parallel-layout data,
        and an ``obj`` named ``local_image`` containing the current image index,
        energy histories, positions, structures, and constrained forces. Energy
        quantities are in Hartree, coordinates are in bohr, forces are in
        Hartree per bohr, and the spring constant is in Hartree per bohr squared.
        Missing or malformed companion data leave the affected fields as
        ``None`` without preventing extraction of the remaining fields.
        """
        def read_assignments(filepath):
            """Read quoted, possibly multiline controller assignments."""
            values = {}
            with open(filepath,'r') as fobj:
                control_lines = fobj.read().splitlines()
            i = 0
            while i<len(control_lines):
                line = control_lines[i]
                i += 1
                if '=' not in line:
                    continue
                key,value = line.split('=',1)
                key       = key.strip()
                value     = value.strip()
                if len(key)==0 or len(value)==0:
                    continue
                quote = value[0] if value[0] in {'"',"'"} else None
                if quote is None:
                    values[key] = value
                    continue
                value = value[1:]
                if value.endswith(quote):
                    values[key] = value[:-1]
                    continue
                parts = [value]
                while i<len(control_lines):
                    part = control_lines[i]
                    i += 1
                    if quote in part:
                        parts.append(part.split(quote,1)[0])
                        break
                    parts.append(part)
                values[key] = '\n'.join(parts)
            return values
        #end def read_assignments

        def as_int(value):
            """Convert a controller value to an integer when possible."""
            try:
                return int(value)
            except (TypeError,ValueError):
                return None
        #end def as_int

        def as_float(value):
            """Convert a controller value to a finite float when possible."""
            try:
                result = float(value.replace('D','E').replace('d','e'))
            except (AttributeError,TypeError,ValueError):
                return None
            return result if np.isfinite(result) else None
        #end def as_float

        def final_energy(filepath):
            """Return the last final eigenvalue-sum energy in Hartree."""
            value = None
            units = None
            try:
                with open(filepath,'r') as fobj:
                    energy_lines = fobj.read().splitlines()
            except (OSError,UnicodeError):
                return None
            for line in energy_lines:
                lower = ' '.join(line.expandtabs().strip().split()).lower()
                if (
                    'final total energy from eig sum' not in lower and
                    'final total energy from eigenvalue sum' not in lower
                    ):
                    continue
                # Find the first RMG-formatted number following the energy label.
                # Example: final total energy from eig sum = -1.73498899 Ha
                matches = re.findall(
                    self.number_pattern,line.split('=',1)[-1])
                if len(matches)==0:
                    continue
                candidate = as_float(matches[0])
                if candidate is None:
                    continue
                value = candidate
                for unit in {'eV','Ha','Ry'}:
                    if unit.lower() in lower.split():
                        units = unit
                        break
            if value is None or units is None:
                return None
            try:
                return convert(value,units,'Ha')
            except (KeyError,TypeError,ValueError):
                return None
        #end def final_energy

        neb = obj(
            controller_file           = None,
            calculation_mode          = None,
            spin_polarized            = None,
            num_intermediate_images   = None,
            num_images                = None,
            images_per_node           = None,
            max_steps                 = None,
            spring_constant           = None,
            spring_constant_units     = 'Ha/B^2',
            initial_input_file        = None,
            final_input_file          = None,
            image_directories         = None,
            image_input_files         = None,
            image_log_files           = None,
            image_mpi_processes       = None,
            input_structures          = None,
            reaction_coordinate       = None,
            reaction_coordinate_units = 'B',
            energies                  = None,
            relative_energies         = None,
            energy_units              = 'Ha',
            forward_barrier           = None,
            reverse_barrier           = None,
            barrier_image_index       = None,
            parallel                  = None,
            local_image               = None,
            )

        # Match the NEB image and MPI layout printed at RMG initialization.
        # Example: RMG initialization ... 2 image(s) total, 1 per node. 2 MPI processes/image.
        layout_pattern = (
            r'RMG\s+initialization.*?(\d+)\s+image\(s\)\s+total\s*,?\s*'
            r'(\d+)\s+per\s+node\s*\.?\s*(\d+)\s+MPI\s+process(?:es)?/image')
        layout = None
        for line in lines:
            match = re.search(layout_pattern,line,re.IGNORECASE)
            if match is not None:
                layout = obj(
                    num_intermediate_images = int(match.group(1)),
                    images_per_node         = int(match.group(2)),
                    mpi_processes_per_image = int(match.group(3)),
                    )
                break
        neb.parallel = layout

        constrained_images = []
        neb_calls          = 0
        for line in lines:
            lower = ' '.join(line.expandtabs().strip().split()).lower()
            if 'neb call' in lower:
                neb_calls += 1
            # Match the one-based intermediate-image index for constrained forces.
            # Example: Entering constrained forces for image 2
            match = re.search(
                r'entering\s+constrained\s+forces\s+for\s+image\s+(\d+)',
                line,re.IGNORECASE)
            if match is not None:
                constrained_images.append(int(match.group(1)))

        image_index = constrained_images[-1] if len(constrained_images)>0 else None
        neb.local_image = obj(
            index                   = image_index,
            neb_calls               = neb_calls,
            constrained_force_calls = np.array(constrained_images,dtype=int),
            energy                  = self.energy,
            energy_units            = self.energy_units,
            energies                = self.energies,
            positions               = self.positions,
            position_units          = self.position_units,
            forces                  = self.forces,
            force_units             = self.force_units,
            max_forces              = self.max_forces,
            structures              = self.structures,
            )

        controller_candidates = [
            os.path.join(self.path,'ctrl_init.dat'),
            os.path.join(os.path.dirname(self.path),'ctrl_init.dat'),
            ]
        controller_file = None
        for filepath in controller_candidates:
            if os.path.isfile(filepath):
                controller_file = os.path.abspath(filepath)
                break
        if controller_file is not None:
            neb.controller_file = controller_file
            try:
                control = read_assignments(controller_file)
            except (OSError,UnicodeError):
                control = {}
            controller_path     = os.path.dirname(controller_file)

            nintermediate               = as_int(control.get('num_images'))
            neb.num_intermediate_images = nintermediate
            if nintermediate is not None:
                neb.num_images = nintermediate+2
            neb.calculation_mode = control.get('calculation_mode')
            neb.images_per_node  = as_int(control.get('image_per_node'))
            neb.max_steps        = as_int(control.get('max_neb_steps'))
            neb.spring_constant  = as_float(control.get('neb_spring_constant'))
            spin_polarized       = control.get('spin_polarization')
            if spin_polarized is not None:
                spin_polarized = spin_polarized.strip().lower()
                if spin_polarized in {'true','yes','1'}:
                    neb.spin_polarized = True
                elif spin_polarized in {'false','no','0'}:
                    neb.spin_polarized = False

            initial_file = control.get('input_file_initial_image')
            final_file   = control.get('input_file_final_image')
            if initial_file is not None:
                neb.initial_input_file = os.path.abspath(
                    os.path.join(controller_path,initial_file))
            if final_file is not None:
                neb.final_input_file = os.path.abspath(
                    os.path.join(controller_path,final_file))

            image_directories   = []
            image_input_files   = []
            image_mpi_processes = []
            image_info          = control.get('image_infos')
            if image_info is not None:
                for line in image_info.splitlines():
                    tokens = line.split()
                    if len(tokens)<2:
                        continue
                    directory = os.path.abspath(
                        os.path.join(controller_path,tokens[0]))
                    image_directories.append(directory)
                    image_input_files.append(os.path.join(directory,tokens[1]))
                    processes = as_int(tokens[2]) if len(tokens)>=3 else None
                    image_mpi_processes.append(processes)
            if len(image_directories)>0:
                neb.image_directories   = image_directories
                neb.image_mpi_processes = np.array(image_mpi_processes,dtype=object)

            ordered_input_files = []
            if neb.initial_input_file is not None:
                ordered_input_files.append(neb.initial_input_file)
            ordered_input_files.extend(image_input_files)
            if neb.final_input_file is not None:
                ordered_input_files.append(neb.final_input_file)
            if len(ordered_input_files)>0:
                neb.image_input_files = ordered_input_files

            input_structures = []
            for filepath in ordered_input_files:
                structure = None
                if os.path.isfile(filepath):
                    try:
                        structure = RmgInput(filepath).return_structure('B')
                    except (NexusError,KeyError,TypeError,ValueError):
                        structure = None
                input_structures.append(structure)
            if any(structure is not None for structure in input_structures):
                neb.input_structures = input_structures
            if (
                len(input_structures)>0 and
                all(structure is not None for structure in input_structures)
                ):
                coordinate           = [0.0]
                valid_coordinate     = True
                for previous,current in zip(
                    input_structures[:-1],
                    input_structures[1:],
                    strict = True,
                    ):
                    if previous.pos.shape!=current.pos.shape:
                        valid_coordinate = False
                        break
                    displacement = current.pos-previous.pos
                    coordinate.append(
                        coordinate[-1]+np.sqrt((displacement**2).sum()))
                if valid_coordinate:
                    neb.reaction_coordinate = np.array(coordinate,dtype=float)

            image_logs = []
            for filepath in ordered_input_files:
                directory = os.path.dirname(filepath)
                prefix    = os.path.basename(filepath)
                logs      = sorted(glob(os.path.join(directory,prefix+'.*.log')))
                image_logs.append(logs[-1] if len(logs)>0 else None)
            if len(image_logs)>0:
                neb.image_log_files = image_logs

            energies = np.full(len(ordered_input_files),np.nan,dtype=float)
            if len(energies)>0:
                endpoint_initial_energy = as_float(
                    control.get('totale_initial_image'))
                endpoint_final_energy   = as_float(
                    control.get('totale_final_image'))
                if endpoint_initial_energy is not None:
                    energies[0] = endpoint_initial_energy
                if endpoint_final_energy is not None:
                    energies[-1] = endpoint_final_energy
                for i,filepath in enumerate(image_logs[1:-1],start=1):
                    if filepath is not None:
                        value = final_energy(filepath)
                        if value is not None:
                            energies[i] = value
                if np.any(np.isfinite(energies)):
                    neb.energies = energies
                if np.all(np.isfinite(energies)):
                    relative                   = energies-energies[0]
                    barrier_index              = int(np.argmax(energies))
                    neb.relative_energies      = relative
                    neb.forward_barrier        = energies[barrier_index]-energies[0]
                    neb.reverse_barrier        = energies[barrier_index]-energies[-1]
                    neb.barrier_image_index    = barrier_index

        self.neb = neb
    #end def read_neb


    def read_produced_files(self):
        """Locate files produced by EXX, STM, and QMCPACK-interface runs.

        Notes
        -----
        Binds ``produced_files`` to an ``obj`` whose members are lists of paths
        for EXX or STM products or a string path for a QMCPACK restart file.
        """
        produced_files = obj()
        mode           = self.run_mode
        if mode=='exx':
            exx_files = sorted(glob(os.path.join(self.path,'*exx*integral*.h5')))
            if len(exx_files)>0:
                produced_files.exx_integrals = exx_files
        if mode=='stm':
            stm_files      = sorted(glob(os.path.join(self.path,'STM','*.stm')))
            stm_cube_files = sorted(glob(os.path.join(self.path,'STM','*.cube')))
            if len(stm_files)>0:
                produced_files.stm = stm_files
            if len(stm_cube_files)>0:
                produced_files.stm_cube = stm_cube_files
        elif mode=='scf' and 'files' in self.setup_info:
            data_output = self.setup_info.files.get('data_output_file',None)
            if data_output is not None:
                qmcpack_file = os.path.join(self.path,data_output+'.h5')
                if os.path.isfile(qmcpack_file):
                    produced_files.qmcpack_restart = qmcpack_file
        if len(produced_files)>0:
            self.produced_files = produced_files
    #end def read_produced_files


    def return_initial_structure(self):
        """Return a detached copy of the initial structure, when available."""
        s = None
        if 'setup_info' in self and 'structure' in self.setup_info:
            s = deepcopy(self.setup_info.structure)
        return s
    #end def return_initial_structure

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
        Short RMG calculation mode determined during analysis, such as
        ``"scf"``, ``"band"``, ``"relax"``, or ``"md_VE"``.
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

    all_modes         = frozenset({
        'scf','nscf','band','exx','relax','md_VE','md_TE','tddft','stm','neb',
        })
    electronic_modes  = frozenset({
        'scf','nscf','relax','md_VE','md_TE','tddft','neb',
        })
    eigenvalue_modes  = electronic_modes|frozenset({'band'})
    relaxation_modes  = frozenset({'relax','neb'})
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
        structure = self.results.return_initial_structure()
        if structure is not None:
            structure.change_units(units)
        return structure
    #end def initial_structure


    def energy(self,units='Ha'):
        """Return the final total energy in eV, Hartree, or Rydberg."""
        self._require_supported('energy',self.electronic_modes)
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
        self._require_supported('eigenvalues',self.eigenvalue_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('eigenvalues units must be one of: eV, Ha, Ry')
        if self.run_mode=='band':
            bands = self.results.bands
            if bands is None:
                return None
            spin_values = []
            for spin in sorted(bands.keys()):
                values = np.asarray(bands[spin].energies,dtype=float).T
                if values.ndim!=2:
                    return None
                spin_values.append(values)
            if len(spin_values)==0:
                return None
            shape = spin_values[0].shape
            if any(values.shape!=shape for values in spin_values):
                return None
            if len(spin_values)==1:
                eigenvalues = spin_values[0]
            else:
                eigenvalues = np.stack(spin_values,axis=1)
            return convert(eigenvalues,'eV',units)
        electronic = self.results.electronic
        if electronic is None or 'eigenvalues' not in electronic:
            return None
        return convert(electronic.eigenvalues,'eV',units)
    #end def eigenvalues


    def occupations(self):
        """Return the dimensionless K-point-major occupation array."""
        self._require_supported('occupations',self.electronic_modes)
        electronic = self.results.electronic
        if electronic is None or 'occupations' not in electronic:
            return None
        return electronic.occupations
    #end def occupations


    def Ef(self,units='eV'):
        """Return the final Fermi energy in eV, Hartree, or Rydberg."""
        self._require_supported('Ef',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Ef units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.fermi_energies)==0:
            return None
        return convert(electronic.fermi_energies[-1],'eV',units)
    #end def Ef


    def Evbm(self,units='eV'):
        """Return the final valence-band maximum in selected energy units."""
        self._require_supported('Evbm',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Evbm units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.valence_band_maxima)==0:
            return None
        return convert(electronic.valence_band_maxima[-1],'eV',units)
    #end def Evbm


    def Ecbm(self,units='eV'):
        """Return the final conduction-band minimum in selected energy units."""
        self._require_supported('Ecbm',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('Ecbm units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.conduction_band_minima)==0:
            return None
        return convert(electronic.conduction_band_minima[-1],'eV',units)
    #end def Ecbm


    def band_gap(self,units='eV'):
        """Return the final band gap in eV, Hartree, or Rydberg."""
        self._require_supported('band_gap',self.electronic_modes)
        if units not in {'eV','Ha','Ry'}:
            raise ValueError('band_gap units must be one of: eV, Ha, Ry')
        electronic = self.results.electronic
        if electronic is None or len(electronic.band_gaps)==0:
            return None
        return convert(electronic.band_gaps[-1],'eV',units)
    #end def band_gap


    def fractional_occs(self):
        """Whether any occupation differs from empty or full by over ``1e-3``."""
        self._require_supported('fractional_occs',self.electronic_modes)
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
        self._require_supported('forces',self.electronic_modes)
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
        self._require_supported('stress',self.electronic_modes)
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
        self._require_supported('pressure',self.electronic_modes)
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
        elif isinstance(arg0,Simulation):
            sim = arg0

            path     = sim.locdir
            filename = sim.outfile
        else:
            filepath = arg0
            if not isinstance(filepath,str):
                msg = (
                    'invalid type provided for log_file\n'
                    'Type expected: str\n'
                    f'Type provided: {filepath.__class__.__name__}'
                    )
                raise TypeError(msg)
            elif not os.path.exists(filepath):
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
            path,filename = os.path.split(filepath)

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


    def return_initial_structure(self):
        """Return the initial structure recovered from the analyzed output."""
        return self.results.return_initial_structure()
    #end def return_initial_structure

#end class RmgAnalyzer
