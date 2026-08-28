##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os
import re
from glob import glob
from copy import deepcopy
import numpy as np
from .developer import obj
from .generic import NexusError
from .utilities import to_str
from .fileio import TextFile
from .unit_converter import convert
from .numerics import simstats
from .simulation import SimulationAnalyzer, Simulation
from .structure import generate_structure
from .rmg_input import RmgInput, rmg_modes


class RmgAnalyzer(SimulationAnalyzer):

    common_result_fields = (
        'geometry',
        'convergence',
        'timing',
        )

    electronic_result_fields = (
        'energy',
        'energy_units',
        'energies',
        'energy_units_history',
        'direct_energies',
        'direct_energy_units',
        'electronic',
        'scf',
        )

    ionic_result_fields = (
        'ionic_steps',
        'position_units',
        'force_units',
        'positions',
        'forces',
        'charges',
        'magnetizations',
        'max_forces',
        'structures',
        'stress',
        'stress_units',
        'pressures',
        'pressure',
        )

    result_fields_by_mode = dict(
        scf = common_result_fields+electronic_result_fields+ionic_result_fields+(
            'produced_files',
            ),
        nscf = common_result_fields+electronic_result_fields+ionic_result_fields,
        band = common_result_fields+(
            'electronic',
            'bands',
            ),
        exx = common_result_fields+(
            'produced_files',
            ),
        relax = common_result_fields+electronic_result_fields+ionic_result_fields,
        md_VE = common_result_fields+electronic_result_fields+ionic_result_fields+(
            'md',
            'md_stats',
            ),
        md_TE = common_result_fields+electronic_result_fields+ionic_result_fields+(
            'md',
            'md_stats',
            ),
        tddft = common_result_fields+electronic_result_fields+ionic_result_fields+(
            'tddft',
            ),
        stm = common_result_fields+(
            'produced_files',
            ),
        neb = common_result_fields+electronic_result_fields+ionic_result_fields,
        )

    @property
    def initialized(self):
        return self.path is not None and self.outfile_name is not None
    #end def initialized

    @property
    def analyzed(self):
        return 'parse_status' in self.info and self.info.parse_status.get('log',False)
    #end def analyzed

    @property
    def analysis_succeeded(self):
        succeeded = self.analyzed and len(self.setup_info)>0 and len(self.results)>0
        if succeeded and 'info' in self and 'parse_status' in self.info:
            status = self.info.parse_status
            succeeded = (
                status.get('setup',False) and
                status.get('results',False) and
                len(status.errors)==0
                )
        #end if
        return succeeded
    #end def analysis_succeeded

    @property
    def run_completed(self):
        return self.analyzed and self.results.get('timing',None) is not None
    #end def run_completed

    @property
    def calculation_mode(self):
        mode = self.calculation_shortmode
        if mode is not None:
            mode = rmg_modes.full_mode(mode)
        #end if
        return mode
    #end def calculation_mode

    @property
    def calculation_shortmode(self):
        mode = None
        if 'setup_info' in self and 'run_mode' in self.setup_info:
            mode = self.setup_info.run_mode
        #end if
        return mode
    #end def calculation_shortmode


    @staticmethod
    def normalize_line(line):
        return ' '.join(line.expandtabs().strip().split()).lower()
    #end def normalize_line


    @staticmethod
    def identify_mode(text):
        text = RmgAnalyzer.normalize_line(text)
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
        for pattern,mode in mode_patterns:
            if re.search(pattern,text,re.I):
                return mode
            #end if
        #end for
        return rmg_modes.mode_match(text,short=True)
    #end def identify_mode


    # Match a signed integer or decimal with an optional E- or D-exponent.
    # Example: -1.2345D+02
    number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


    @staticmethod
    def line_numbers(line):
        # Find every RMG-formatted number occurring anywhere in a line.
        # Example: SUM FORCE = 0.1 0.2 0.3
        values = re.findall(
            RmgAnalyzer.number_pattern,line.replace('D','E').replace('d','e'))
        return np.array(values,dtype=float)
    #end def line_numbers


    @staticmethod
    def leading_numbers(line):
        # Match a whitespace-separated numeric sequence at the start of a line.
        # Example: 1 1.0 0.1 0.2 trailing
        npat = RmgAnalyzer.number_pattern
        pattern = r'^\s*((?:'+npat+r')(?:\s+(?:'+npat+r'))*)'
        match = re.match(pattern,line)
        if match is None:
            return np.array([],dtype=float)
        #end if
        return RmgAnalyzer.line_numbers(match.group(1))
    #end def leading_numbers


    @staticmethod
    def match_float(pattern,line):
        # Apply a caller-provided expression that captures one physical value.
        # Example: FERMI ENERGY : 5.25 eV
        match = re.search(pattern,line,re.I)
        if match is None:
            return None
        #end if
        return float(match.group(1).replace('D','E').replace('d','e'))
    #end def match_float


    def __init__(self,arg0=None,*,analyze=False):
        self.path         = None
        self.abspath      = None
        self.outfile_name = None
        self._reset_analysis_data()

        if arg0 is None:
            return
        elif isinstance(arg0,Simulation):
            sim = arg0
            path     = sim.locdir
            filename = sim.outfile
        else:
            log_file = arg0
            if not isinstance(log_file,str):
                msg = (
                    'invalid type provided for log_file\n'
                    'Type expected: str\n'
                    'Type provided: {}'.format(log_file.__class__.__name__)
                    )
                raise TypeError(msg)
            elif not os.path.exists(log_file):
                msg = (
                    'RMG log output file does not exist.\n'
                    'Path provided: {}'.format(log_file)
                    )
                raise FileNotFoundError(msg)
            elif not os.path.isfile(log_file):
                msg = (
                    'Path provided for RMG log output is not a file.\n'
                    'Path provided: {}'.format(log_file)
                    )
                raise IsADirectoryError(msg)
            #end if
            path,filename = os.path.split(log_file)
        #end if
        
        self.path         = path
        self.abspath      = os.path.abspath(path)
        self.outfile_name = filename

        if analyze:
            self.analyze()
        #end if
    #end def __init__


    def _reset_analysis_data(self):
        self.info       = obj()
        self.input      = None
        self.setup_info = obj()
        self.results    = obj()
    #end def _reset_analysis_data


    def _initialize_results(self,mode):
        fields = self.result_fields_by_mode.get(mode,())
        for name in fields:
            self.results[name] = None
        #end for
    #end def _initialize_results


    def analyze(self,*,guard=True):
        if not self.initialized:
            return
        #end if
        self._reset_analysis_data()
        log_filepath = os.path.join(self.path,self.outfile_name)
        if not os.path.exists(log_filepath):
            msg = (
                'RMG analysis cannot be completed.\n'
                'Log file does not exist at path provided.\n'
                'Path provided: {}'.format(log_filepath)
                )
            raise FileNotFoundError(msg)
        #end if
        try:
            with open(log_filepath,'r') as fobj:
                lines = fobj.read().splitlines()
            #end with
        except (OSError,UnicodeError):
            if not guard:
                raise
            #end if
            return
        #end try
        logfile = TextFile(log_filepath)
        parse_status = obj(log=True,errors=obj())
        self.info.parse_status = parse_status
        parse_status.setup = self.read_setup_info(logfile)
        parse_status.results = self.read_results(logfile,lines)
    #end def analyze


    def read_setup_info(self,logfile):
        setup_info = obj()
        f = logfile
        log_text = to_str(f.mm[:])
        log_lines = log_text.splitlines()
        mode = None
        for line in log_lines:
            normalized = self.normalize_line(line)
            # Match the calculation-type heading and capture its value.
            # Example: Calculation type: Quench electrons
            match = re.match(r'^calculation\s*type\s*[:=]\s*(.*)$',normalized)
            if match is not None:
                mode = self.identify_mode(match.group(1))
                break
            #end if
        #end if
        setup_info.run_mode = mode
        self._initialize_results(mode)
        setup_start = None
        setup_end   = None
        if mode is not None:
            setup_start = 'Files'
            setup_end   = 'Initial Ionic Positions And Displacements (Angstrom)'
        else:
            # don't know how to handle other cases yet
            None
        #end if
        unit_set = set(['a0'])
        on_off = dict(ON=True,OFF=False)
        def process_name(s):
            # Match parenthetical annotations so they can be removed from names.
            # Example: Grid spacing (a0)
            s = re.sub(r'\([^)]*\)','',s)
            tokens = s.strip().lower().split()
            name = ''
            for t in tokens:
                if not t.startswith('('):
                    name += t+'_'
                #end if
            #end for
            name = name[:-1].replace('/','_').replace('-','_')
            return name
        #end def process_name
        def process_value(v,*,list=False):
            v = v.strip()
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
                                units = tokens[-1]
                                tokens = tokens[:-1]
                            #end if
                            try:
                                if not list:
                                    v = np.array(tokens,dtype=float)
                                else:
                                    v = [process_value(t,list=True)[0] for t in tokens]
                                #end if
                            except (ValueError,TypeError,OverflowError):
                                units = None
                            #end try
                        #end if
                    elif v in on_off:
                        v = on_off[v]
                    #end if
                #end try
            #end try
            return v,units
        #end def process_value
        if setup_start is not None:
            # Match the standalone Files setup-section heading.
            # Example: Files:
            start_match = re.search(r'(?im)^\s*files\s*(?::\s*)?$',log_text)
            istart = start_match.start() if start_match is not None else -1
            if istart!=-1:
                # Match the normal heading that terminates the setup section.
                # Example: Initial Ionic Positions And Displacements (Angstrom)
                end_pattern = (
                    r'(?im)^\s*initial\s+ionic\s+positions\s+and\s+'
                    r'displacements\s*\(\s*angstroms?\s*\)'
                    )
                end_match = re.search(end_pattern,log_text[istart:])
                if end_match is not None:
                    iend = istart+end_match.start()
                else:
                    # Match alternate post-setup headings when ionic positions are absent.
                    # Example: Davidson converged
                    fallback_pattern = (
                        r'(?im)^\s*(?:diagonalization\s+using|davidson\s+'
                        r'(?:converged|incomplete)|-+\s*timing\s+information)'
                        )
                    end_match = re.search(fallback_pattern,log_text[istart:])
                    iend = istart+end_match.start() if end_match is not None else len(log_text)
                #end if
                if iend!=-1:
                    text = log_text[istart:iend].expandtabs()
                    blocks = []
                    b = ''
                    last_header = False
                    for line in text.splitlines():
                        if len(line)>0:
                            if line[0]!=' ':
                                if last_header:
                                    b += '\n'+line
                                else:
                                    if len(b)>0:
                                        blocks.append(b)
                                    #end if
                                    b = line
                                    last_header = True
                                #end if
                            else:
                                b += '\n'+line
                                last_header = False
                            #end if
                        #end if
                    #end for
                    if len(b)>0:
                        blocks.append(b)
                    #end if
                    other_blocks = obj()
                    for b in blocks:
                        if '\n' not in b:
                            continue
                        #end if
                        header,body = b.split('\n',1)
                        bname = process_name(header)
                        lines = body.splitlines()
                        value_lines = [line for line in lines if ':' in line]
                        table_blocks = set([
                            'k_points',
                            'initial_ionic_positions_and_displacements',
                            ])
                        simple_values = bname not in table_blocks and len(value_lines)>0
                        if simple_values:
                            bvalues = obj()
                            for line in value_lines:
                                name,value = line.split(':',1)
                                name  = process_name(name)
                                value,units = process_value(value)
                                bvalues[name] = value
                                if units is not None:
                                    bvalues.units = units
                                #end if
                            #end for
                            setup_info[bname] = bvalues
                        else:
                            other_blocks[bname] = header,body,lines
                        #end if
                    #end for
                    # additional processing for specific blocks
                    if 'grid_points' in setup_info:
                        b = setup_info.grid_points
                        grid = []
                        grid_pe = []
                        spacing = []
                        grid_units = None
                        for c in 'xyz':
                            if c not in b:
                                continue
                            #end if
                            s = str(b[c])
                            # Match and capture the total grid-point count.
                            # Example: Total: 48 Per PE: 12 Spacing: 0.30 a0
                            total = self.match_float(
                                r'\btotal\s*:\s*('+self.number_pattern+r')',s)
                            # Match and capture the per-process grid-point count.
                            # Example: Total: 48 Per PE: 12 Spacing: 0.30 a0
                            per_pe = self.match_float(
                                r'\bper\s*pe\s*:\s*('+self.number_pattern+r')',s)
                            # Match and capture the real-space grid spacing.
                            # Example: Total: 48 Per PE: 12 Spacing: 0.30 a0
                            space = self.match_float(
                                r'\bspacing\s*:\s*('+self.number_pattern+r')',s)
                            if total is None or per_pe is None or space is None:
                                continue
                            #end if
                            # Match the unit following a grid-spacing value.
                            # Example: Spacing: 0.30 a0
                            unit_match = re.search(
                                r'\bspacing\s*:\s*'+self.number_pattern+
                                r'\s*([^\s,;]+)',s,re.I)
                            if unit_match is not None:
                                grid_units = unit_match.group(1)
                            #end if
                            grid.append(total)
                            grid_pe.append(per_pe)
                            spacing.append(space)
                        #end for
                        if len(grid)==3:
                            grid = np.array(grid,dtype=int)
                            grid_pe = np.array(grid_pe,dtype=int)
                            spacing = np.array(spacing,dtype=float)
                            b.update(
                                grid         = grid,
                                grid_pe      = grid_pe,
                                grid_spacing = spacing,
                                grid_units   = grid_units,
                                )
                        #end if
                        if 'equivalent_energy_cutoffs' in b:
                            cutoff_text = str(b.equivalent_energy_cutoffs)
                            cutoff_values = self.line_numbers(cutoff_text)
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
                            #end if
                        #end if
                    #end if
                    if 'lattice_setup' in setup_info:
                        b = setup_info.lattice_setup
                        axes = []
                        for name in ('x_basis_vector','y_basis_vector','z_basis_vector'):
                            if name not in b:
                                axes = []
                                break
                            #end if
                            values = self.line_numbers(str(b[name]))
                            if len(values)<3:
                                axes = []
                                break
                            #end if
                            axes.append(values[:3])
                        #end for
                        if len(axes)==3:
                            b.axes = np.array(axes,dtype=float)
                        #end if
                    #end if
                    if 'k_points' in other_blocks:
                        header,body,lines = other_blocks.k_points
                        del other_blocks.k_points
                        first_row = None
                        for i,line in enumerate(lines):
                            normalized = self.normalize_line(line)
                            if 'weight' in normalized and 'crystal' in normalized:
                                first_row = i+1
                                break
                            #end if
                        #end for
                        if first_row is not None:
                            kp = []
                            kw = []
                            for line in lines[first_row:]:
                                normalized = self.normalize_line(line)
                                if 'weight' in normalized and len(kp)>0:
                                    break
                                #end if
                                values = self.leading_numbers(line)
                                if len(values)<4:
                                    continue
                                #end if
                                kp.append(values[:3])
                                kw.append(values[3])
                            #end for
                            if len(kp)>0:
                                setup_info.k_points = obj(
                                    kpoints_crystal = np.array(kp,dtype=float),
                                    kweights        = np.array(kw,dtype=float),
                                    )
                            #end if
                        #end if
                    #end if
                    k = 'initial_ionic_positions_and_displacements'
                    if k in other_blocks:
                        header,body,lines = other_blocks[k]
                        del other_blocks[k]
                        h = header.lower()
                        punits = None
                        if 'bohr' in h:
                            punits = 'B'
                        elif 'angstrom' in h:
                            punits = 'A'
                        #end if
                        pos = []
                        spec = []
                        first_row = None
                        for i,line in enumerate(lines):
                            if 'species' in self.normalize_line(line):
                                first_row = i+1
                                break
                            #end if
                        #end for
                        if first_row is not None:
                            for line in lines[first_row:]:
                                ls = line.strip()
                                if len(ls)>0:
                                    t = line.split()
                                    if len(t)<4:
                                        continue
                                    #end if
                                    values = self.line_numbers(' '.join(t[1:]))
                                    if len(values)<3:
                                        continue
                                    #end if
                                    spec.append(t[0])
                                    pos.append(values[:3])
                                #end if
                            #end for
                            if len(pos)>0:
                                setup_info.ion_positions = obj(
                                    units     = punits,
                                    atoms     = np.array(spec,dtype=object),
                                    positions = np.array(pos,dtype=float),
                                    )
                            #end if
                        #end if
                    #end if
                #end if
            #end if
        #end if
        have_structure_data = (
            'lattice_setup' in setup_info and
            'axes' in setup_info.lattice_setup and
            'ion_positions' in setup_info and
            setup_info.ion_positions.get('units',None) in ('A','B')
            )
        if have_structure_data:
            aunits = setup_info.lattice_setup.get('units','B')
            axes   = np.asarray(setup_info.lattice_setup.axes)
            elem   = np.asarray(setup_info.ion_positions.atoms)
            pos    = np.asarray(setup_info.ion_positions.positions)
            punits = setup_info.ion_positions.units
            valid_shapes = (
                axes.shape==(3,3) and
                pos.ndim==2 and pos.shape[1:]==(3,) and
                len(elem)==len(pos) and len(elem)>0
                )
            if valid_shapes:
                if aunits in ('a0','B'):
                    aunits = 'B'
                else:
                    aunits = 'A' # assume for now
                #end if
                units = 'B'
                axes = convert(axes,aunits,units)
                pos  = convert(pos,punits,units)
                s = generate_structure(
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
                    #end if
                #end if
                setup_info.structure = s
            #end if
        #end if
        if 'files' in setup_info and 'control_input_file' in setup_info.files:
            control_file = setup_info.files.control_input_file
            filepaths = [
                os.path.join(self.path,control_file),
                os.path.join(self.path,os.path.basename(control_file)),
                os.path.join(os.path.dirname(self.path),control_file),
                ]
            for filepath in filepaths:
                if os.path.isfile(filepath):
                    break
                #end if
            #end for
            if os.path.isfile(filepath):
                try:
                    self.input = RmgInput(filepath)
                except (NexusError,OSError,TypeError,ValueError):
                    pass
                #end try
            #end if
        #end if
        self.setup_info = setup_info
        return len(setup_info)>0
    #end def read_setup_info


    def read_results(self,logfile,lines=None):
        if lines is None:
            lines = to_str(logfile.mm[:]).splitlines()
        #end if
        if len(self.results)==0 and 'run_mode' in self.setup_info:
            self._initialize_results(self.setup_info.run_mode)
        #end if
        parsers = (
            ('energies',self.analyze_energies),
            ('electronic',self.analyze_electronic),
            ('scf',self.analyze_scf),
            ('ions',self.analyze_ions),
            ('md',self.analyze_md),
            ('geometry',self.analyze_geometry),
            ('stress',self.analyze_stress),
            ('convergence',self.analyze_convergence),
            ('timing',self.analyze_timing),
            ('band',self.analyze_band),
            ('tddft',self.analyze_tddft),
            ('produced_files',self.analyze_produced_files),
            )
        parse_status = self.info.parse_status
        for name,parser in parsers:
            parse_status[name] = parser(lines)
        #end for
        return any(value is not None for value in self.results.values())
    #end def read_results


    def analyze_energies(self,lines):
        energies = []
        direct_energies = []
        energy_units = []
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
            #end if
            match = dpat.search(line)
            if match is not None:
                direct_energies.append(float(match.group(1).replace('D','E').replace('d','e')))
                direct_energy_units.append(match.group(2))
            #end if
        #end for
        energies = np.array(energies,dtype=float)
        direct_energies = np.array(direct_energies,dtype=float)
        if len(energies)>0:
            self.results.energy = energies[-1]
            self.results.energy_units = energy_units[-1] or 'Ha'
            self.results.energies = energies
            self.results.energy_units_history = np.array(energy_units,dtype=object)
        #end if
        if len(direct_energies)>0:
            self.results.direct_energies = direct_energies
            self.results.direct_energy_units = np.array(direct_energy_units,dtype=object)
        #end if
        return len(energies)
    #end def analyze_energies


    def analyze_electronic(self,lines):
        data = obj(
            fermi_energies          = [],
            valence_band_maxima     = [],
            conduction_band_minima = [],
            band_gaps               = [],
            total_charges           = [],
            total_magnetizations    = [],
            absolute_magnetizations = [],
            )
        sum_forces = []
        volume_per_atom = []
        energy_per_atom = []
        for line in lines:
            lower = self.normalize_line(line)
            if 'fermi energy' in lower:
                # Match and capture the Fermi energy.
                # Example: FERMI ENERGY : 5.25 eV
                value = self.match_float(
                    r'fermi\s+energy\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.fermi_energies.append(value)
                #end if
            elif 'valence band maximum' in lower and 'conduction band' in lower:
                # Match and capture the valence-band maximum.
                # Example: spin0: valence band maximum = 4.0 eV
                vbm = self.match_float(
                    r'valence\s+band\s+maximum\s*[:=]\s*('+self.number_pattern+r')',line)
                # Match and capture the conduction-band minimum.
                # Example: conduction band minimum = 6.0 eV
                cbm = self.match_float(
                    r'conduction\s+band\s+(?:minimum|minumm)\s*[:=]\s*'
                    r'('+self.number_pattern+r')',line)
                if vbm is not None and cbm is not None:
                    data.valence_band_maxima.append(vbm)
                    data.conduction_band_minima.append(cbm)
                #end if
            elif 'band gap' in lower:
                # Match and capture the electronic band gap.
                # Example: spin0: Band gap : 2.0 eV
                value = self.match_float(
                    r'band\s+gap\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.band_gaps.append(value)
                #end if
            elif 'total charge in supercell' in lower:
                # Match and capture the total charge in the simulation cell.
                # Example: Total charge in supercell = 0.0
                value = self.match_float(
                    r'total\s+charge\s+in\s+supercell\s*[:=]\s*'
                    r'('+self.number_pattern+r')',line)
                if value is not None:
                    data.total_charges.append(value)
                #end if
            elif '@@ total magnetization' in lower:
                # Match and capture the signed total magnetization.
                # Example: @@ TOTAL MAGNETIZATION = 1.0
                value = self.match_float(
                    r'total\s+magnetization\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.total_magnetizations.append(value)
                #end if
            elif '@@ absolute magnetization' in lower:
                # Match and capture the absolute total magnetization.
                # Example: @@ ABSOLUTE MAGNETIZATION = 1.5
                value = self.match_float(
                    r'absolute\s+magnetization\s*[:=]\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.absolute_magnetizations.append(value)
                #end if
            elif lower.lstrip().startswith('sum force'):
                values = self.line_numbers(line.split('=',1)[-1])
                if len(values)>=3:
                    sum_forces.append(values[:3])
                #end if
            elif 'volume and energy per atom' in lower:
                values = self.line_numbers(line.split('=',1)[-1])
                if len(values)>=2:
                    volume_per_atom.append(values[0])
                    energy_per_atom.append(values[1])
                #end if
            #end if
        #end for
        for name,values in data.items():
            data[name] = np.array(values,dtype=float)
        #end for
        data.energy_units = 'eV'
        data.magnetization_units = 'Bohr mag/cell'
        data.sum_forces = np.array(sum_forces,dtype=float)
        data.sum_force_units = 'Ha/a0'
        data.volume_per_atom = np.array(volume_per_atom,dtype=float)
        data.energy_per_atom = np.array(energy_per_atom,dtype=float)
        data.energy_per_atom_units = 'eV'
        nfound = sum(len(v) for v in data.values() if isinstance(v,np.ndarray))
        if nfound>0:
            self.results.electronic = data
        #end if
        return nfound
    #end def analyze_electronic


    def analyze_scf(self,lines):
        component_names = obj(
            eigenvalue_sum = 'EIGENVALUE SUM',
            ion_ion        = 'ION_ION',
            electrostatic  = 'ELECTROSTATIC',
            vxc            = 'VXC',
            exc            = 'EXC',
            total_energy   = 'TOTAL ENERGY',
            estimated_error = 'estimated error',
            )
        values = obj()
        for name in component_names.keys():
            values[name] = []
        #end for
        md_steps = []
        scf_steps = []
        step_times = []
        scf_times = []
        rms_dv = []
        for line in lines:
            for name,label in component_names.items():
                # Match an @@ energy component and capture its value or overflow stars.
                # Example: @@ TOTAL ENERGY : -1.250000
                pattern = (
                    r'^\s*@@\s*'+re.escape(label)+r'\s*[:=]\s*'
                    r'('+self.number_pattern+r'|\*+)'
                    )
                match = re.search(pattern,line,re.I)
                if match is not None:
                    token = match.group(1)
                    if '*' in token:
                        value = np.nan
                    else:
                        value = float(token.replace('D','E').replace('d','e'))
                    #end if
                    values[name].append(value)
                    break
                #end if
            #end for
            # Match the detailed SCF-iteration summary line.
            # Example: quench: [ RMS [ dV ] : 2.0D-5 scf: 3/20 ]
            if re.search(r'\bquench\s*:',line,re.I):
                # Match and capture the current molecular-dynamics step.
                # Example: md: 0/2
                md_match = re.search(r'\bmd\s*:\s*(\d+)\s*/',line,re.I)
                # Match and capture the current SCF iteration.
                # Example: scf: 3/20
                scf_match = re.search(r'\bscf\s*:\s*(\d+)\s*/',line,re.I)
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
                rms_match = re.search(r'\brms\s*\[\s*dv\s*\]\s*:\s*([^\]\s]+)',line,re.I)
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
                #end try
            #end if
        #end for
        scf = obj()
        for name,array in values.items():
            scf[name] = np.array(array,dtype=float)
        #end for
        scf.update(
            md_steps       = np.array(md_steps,dtype=int),
            scf_steps      = np.array(scf_steps,dtype=int),
            step_times     = np.array(step_times,dtype=float),
            scf_times      = np.array(scf_times,dtype=float),
            rms_dv         = np.array(rms_dv,dtype=float),
            energy_units   = 'Ha',
            time_units     = 's',
            )
        nfound = len(scf.total_energy)
        if nfound>0:
            self.results.scf = scf
        #end if
        return nfound
    #end def analyze_scf


    def analyze_ions(self,lines):
        records = []
        i = 0
        while i<len(lines):
            header_tokens = lines[i].split()
            is_header = (
                len(header_tokens)>=3 and
                header_tokens[0].upper()=='@ION' and
                header_tokens[1].lower()=='ion' and
                header_tokens[2].lower()=='species'
                )
            if not is_header:
                i += 1
                continue
            #end if
            atoms = []
            positions = []
            charges = []
            magnetizations = []
            forces = []
            movable = []
            i += 1
            while i<len(lines):
                tokens = lines[i].split()
                if len(tokens)==0 or tokens[0].upper()!='@ION':
                    break
                #end if
                i += 1
                if len(tokens)<14:
                    continue
                #end if
                numeric_tokens = tokens[3:14]
                # Match a complete numeric token before converting an ionic row.
                # Example: 2.1
                valid = all(
                    re.fullmatch(self.number_pattern,v) is not None
                    for v in numeric_tokens)
                if not valid:
                    continue
                #end if
                values = [
                    float(v.replace('D','E').replace('d','e'))
                    for v in numeric_tokens]
                move_values = values[8:11]
                if not all(np.isfinite(v) and v.is_integer() for v in move_values):
                    continue
                #end if
                atom = tokens[2]
                position = values[:3]
                charge = values[3]
                magnetization = values[4]
                force = values[5:8]
                move = [int(v) for v in move_values]
                atoms.append(atom)
                positions.append(position)
                charges.append(charge)
                magnetizations.append(magnetization)
                forces.append(force)
                movable.append(move)
            #end while
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
            #end if
        #end while
        ionic_steps = obj()
        for n,record in enumerate(records):
            ionic_steps[n] = record
        #end for
        if len(records)>0:
            self.results.ionic_steps = ionic_steps
            self.results.position_units = 'a0'
            self.results.force_units = 'Ha/a0'
        #end if
        if len(records)>0 and len(set(len(r.atoms) for r in records))==1:
            positions = np.array([r.positions for r in records],dtype=float)
            forces = np.array([r.forces for r in records],dtype=float)
            charges = np.array([r.charges for r in records],dtype=float)
            magnetizations = np.array([r.magnetizations for r in records],dtype=float)
            max_forces = np.array([
                np.sqrt((r.forces**2).sum(axis=1)).max() for r in records
                ],dtype=float)
            self.results.positions = positions
            self.results.forces = forces
            self.results.charges = charges
            self.results.magnetizations = magnetizations
            self.results.max_forces = max_forces
            if 'structure' in self.setup_info:
                structures = obj()
                initial = self.setup_info.structure
                for n,record in enumerate(records):
                    structure = generate_structure(
                        units = 'B',
                        axes  = initial.axes,
                        elem  = record.atoms,
                        pos   = record.positions,
                        )
                    structures[n] = structure
                #end for
                self.results.structures = structures
            #end if
        #end if
        return len(records)
    #end def analyze_ions


    def analyze_md(self,lines):
        records = []
        for line in lines:
            stripped = line.strip()
            # Match constant-energy or constant-temperature MD record markers.
            # Example: @CVE 1 -1.0 0.1 -0.9 300.0 2.5e-4
            marker_match = re.match(r'^@(CVE|CVT)(?:\b|-)',stripped,re.I)
            if marker_match is None:
                continue
            #end if
            fields = stripped.split(None,1)
            if len(fields)<2:
                continue
            #end if
            values = self.leading_numbers(fields[1])
            if len(values)>=6:
                records.append(obj(
                    step             = int(values[0]),
                    potential_energy = values[1],
                    kinetic_energy   = values[2],
                    total_energy     = values[3],
                    temperature      = values[4],
                    displacement     = values[5],
                    ))
            #end if
        #end for
        if len(records)>0:
            md = obj()
            for name in records[0].keys():
                dtype = int if name=='step' else float
                md[name] = np.array([r[name] for r in records],dtype=dtype)
            #end for
            md.energy_units = 'Ha'
            md.temperature_units = 'K'
            self.results.md = md
            self.results.md_stats = self.md_statistics()
        #end if
        return len(records)
    #end def analyze_md


    def md_statistics(self,equil=None):
        statistics = obj()
        if 'md' not in self.results:
            return statistics
        #end if
        for name,values in self.results.md.items():
            if not isinstance(values,np.ndarray):
                continue
            #end if
            if equil is not None:
                values = values[equil:]
            #end if
            mean,var,error,kappa = simstats(values)
            statistics[name] = obj(
                mean  = mean,
                var   = var,
                error = error,
                kappa = kappa,
                )
        #end for
        return statistics
    #end def md_statistics


    def analyze_geometry(self,lines):
        geometry = obj()
        if 'structure' in self.setup_info:
            structure = self.setup_info.structure
            geometry.volume = abs(np.linalg.det(structure.axes))
            geometry.volume_units = 'a0^3'
        #end if
        if 'k_points' in self.setup_info:
            kpoints = self.setup_info.k_points
            geometry.kpoints_crystal = kpoints.kpoints_crystal
            geometry.kweights = kpoints.kweights
            if 'structure' in self.setup_info and len(kpoints.kpoints_crystal)>0:
                geometry.kpoints_cart = np.dot(
                    kpoints.kpoints_crystal,self.setup_info.structure.kaxes)
            #end if
        #end if
        if len(geometry)>0:
            self.results.geometry = geometry
        #end if
        return len(geometry)>0
    #end def analyze_geometry


    def analyze_stress(self,lines):
        tensors = []
        for i,line in enumerate(lines):
            normalized = self.normalize_line(line)
            # Match the heading for a total stress tensor reported in kbar.
            # Example: stress total in unit of kbar
            if not re.search(r'\bstress\s+total\b',normalized) or 'kbar' not in normalized:
                continue
            #end if
            rows = []
            j = i+1
            while j<len(lines) and len(rows)<3:
                row_line = lines[j]
                j += 1
                if len(row_line.strip())==0:
                    continue
                #end if
                values = self.leading_numbers(row_line)
                if len(values)<3:
                    if len(rows)>0:
                        rows = []
                        break
                    #end if
                    continue
                #end if
                if len(values)>=4 and int(values[0])==len(rows)+1:
                    values = values[1:]
                #end if
                rows.append(values[:3])
            #end while
            if len(rows)==3:
                tensors.append(np.array(rows,dtype=float))
            #end if
        #end for
        if len(tensors)>0:
            stress = np.array(tensors,dtype=float)
            pressures = -np.trace(stress,axis1=1,axis2=2)/3.0
            self.results.stress = stress
            self.results.stress_units = 'kbar'
            self.results.pressures = pressures
            self.results.pressure = pressures[-1]
        #end if
        return len(tensors)
    #end def analyze_stress


    def analyze_convergence(self,lines):
        electronic_successes = 0
        electronic_failures = 0
        ionic_converged = None
        for line in lines:
            lower = self.normalize_line(line)
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
            #end if
        #end for
        electronic_converged = None
        if electronic_successes+electronic_failures>0:
            electronic_converged = electronic_successes>0 and electronic_failures==0
        #end if
        convergence = obj(
            electronic_converged = electronic_converged,
            electronic_successes = electronic_successes,
            electronic_failures  = electronic_failures,
            ionic_converged      = ionic_converged,
            )
        nfound = electronic_successes+electronic_failures+(ionic_converged is not None)
        if nfound>0:
            self.results.convergence = convergence
        #end if
        return nfound
    #end def analyze_convergence


    def analyze_timing(self,lines):
        timing = None
        sections = obj()
        # Match a finite, infinite, or NaN numeric token in the timing table.
        # Example: 3.00
        time_value = r'(?:'+self.number_pattern+r'|inf|nan)'
        # Match a numbered timing section followed by total and per-step times.
        # Example: 1-TOTAL 3.00 0.50
        pattern = re.compile(
            r'^\s*(\d+\s*-\s*.*?)\s+('+time_value+r')\s+('+time_value+r')'
            r'(?:\s+.*)?$',re.I)
        for line in lines:
            match = pattern.match(line)
            if match is None:
                continue
            #end if
            # Match spacing around the first dash so the section name can be normalized.
            # Example: 1 - TOTAL
            name = re.sub(r'\s*-\s*','-',match.group(1).strip(),count=1)
            total = float(match.group(2))
            per_step = float(match.group(3))
            # Match non-alphanumeric runs so a timing name can become an object key.
            # Example: 1-TOTAL
            key = re.sub(r'[^a-z0-9]+','_',name.lower()).strip('_')
            sections[key] = obj(total=total,per_step=per_step)
            if self.normalize_line(name)=='1-total':
                timing = obj(total=total,per_step=per_step,units='s')
            #end if
        #end for
        if timing is not None:
            timing.sections = sections
            self.results.timing = timing
            return True
        #end if
        return False
    #end def analyze_timing


    def analyze_band(self,lines):
        prefix = self.outfile_name[:-4] if self.outfile_name.endswith('.log') else self.outfile_name
        files = sorted(glob(os.path.join(self.path,prefix+'_spin*.bandstructure.dat')))
        bands = obj()
        for filepath in files:
            # Match and capture the spin index in a band-structure filename.
            # Example: input_spin0.bandstructure.dat
            match = re.search(r'_spin(\d+)\.bandstructure\.dat$',filepath)
            if match is None:
                continue
            #end if
            groups = []
            group = []
            with open(filepath,'r') as fobj:
                for line in fobj:
                    if '&&' in line:
                        if len(group)>0:
                            groups.append(np.array(group,dtype=float))
                            group = []
                        #end if
                    else:
                        values = self.leading_numbers(line)
                        if len(values)>=2:
                            group.append(values[:2])
                        #end if
                    #end if
                #end for
            #end with
            if len(group)>0:
                groups.append(np.array(group,dtype=float))
            #end if
            if len(groups)>0 and len(set(len(g) for g in groups))==1:
                spin = int(match.group(1))
                bands[spin] = obj(
                    distance = groups[0][:,0],
                    energies = np.array([g[:,1] for g in groups],dtype=float),
                    energy_units = 'eV',
                    filepath = filepath,
                    )
            #end if
        #end for
        if len(bands)>0:
            self.results.bands = bands
        #end if
        return len(bands)
    #end def analyze_band


    def analyze_tddft(self,lines):
        prefix = self.outfile_name[:-4] if self.outfile_name.endswith('.log') else self.outfile_name
        energy_file = os.path.join(self.path,prefix+'_totalE')
        dipole_files = sorted(glob(os.path.join(self.path,prefix+'_spin*_dipole.dat')))
        tddft = obj()
        if os.path.isfile(energy_file):
            rows = []
            with open(energy_file,'r') as fobj:
                for line in fobj:
                    if line.lstrip().startswith('&&'):
                        continue
                    #end if
                    values = self.leading_numbers(line)
                    if len(values)>=5:
                        rows.append(values[:5])
                    #end if
                #end for
            #end with
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
            #end if
        #end if
        dipoles = obj()
        for filepath in dipole_files:
            rows = []
            field = None
            ground_state = None
            with open(filepath,'r') as fobj:
                for line in fobj:
                    lower = self.normalize_line(line)
                    values = self.line_numbers(line)
                    if 'electric field' in lower and len(values)>=3:
                        values = self.line_numbers(line.split(':',1)[-1])
                        if len(values)>=3:
                            field = values[:3]
                        #end if
                    elif 'dipole at' in lower and len(values)>=3:
                        values = self.line_numbers(line.split(':',1)[-1])
                        if len(values)>=3:
                            ground_state = values[:3]
                        #end if
                    elif not line.lstrip().startswith('&&'):
                        values = self.leading_numbers(line)
                        if len(values)>=4:
                            rows.append(values[:4])
                        #end if
                    #end if
                #end for
            #end with
            if len(rows)>0:
                # Match and capture the spin index in a dipole filename.
                # Example: input_spin0_dipole.dat
                match = re.search(r'_spin(\d+)_dipole\.dat$',filepath)
                spin = int(match.group(1)) if match is not None else len(dipoles)
                values = np.array(rows,dtype=float)
                dipoles[spin] = obj(
                    time         = values[:,0],
                    dipole       = values[:,1:4],
                    electric_field = field,
                    ground_state = ground_state,
                    filepath     = filepath,
                    )
            #end if
        #end for
        if len(dipoles)>0:
            tddft.dipoles = dipoles
        #end if
        if len(tddft)>0:
            self.results.tddft = tddft
        #end if
        return len(tddft)
    #end def analyze_tddft


    def analyze_produced_files(self,lines):
        produced_files = obj()
        mode = self.calculation_shortmode
        if mode=='exx':
            exx_files = sorted(glob(os.path.join(self.path,'*exx*integral*.h5')))
            if len(exx_files)>0:
                produced_files.exx_integrals = exx_files
            #end if
        #end if
        if mode=='stm':
            stm_files = sorted(glob(os.path.join(self.path,'STM','*.stm')))
            stm_cube_files = sorted(glob(os.path.join(self.path,'STM','*.cube')))
            if len(stm_files)>0:
                produced_files.stm = stm_files
            #end if
            if len(stm_cube_files)>0:
                produced_files.stm_cube = stm_cube_files
            #end if
        elif mode=='scf' and 'files' in self.setup_info:
            data_output = self.setup_info.files.get('data_output_file',None)
            if data_output is not None:
                qmcpack_file = os.path.join(self.path,data_output+'.h5')
                if os.path.isfile(qmcpack_file):
                    produced_files.qmcpack_restart = qmcpack_file
                #end if
            #end if
        #end if
        if len(produced_files)>0:
            self.results.produced_files = produced_files
        #end if
        return sum(len(v) for v in produced_files.values())
    #end def analyze_produced_files


    def return_initial_structure(self):
        s = None
        if 'setup_info' in self and 'structure' in self.setup_info:
            s = deepcopy(self.setup_info.structure)
        #end if
        return s
    #end def return_initial_structure

#end class RmgAnalyzer
