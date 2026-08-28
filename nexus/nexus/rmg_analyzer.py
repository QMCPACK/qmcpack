##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################


import os
import re
from glob import glob
from copy import deepcopy
import numpy as np
from .developer import obj
from .utilities import to_str
from .fileio import TextFile
from .unit_converter import convert
from .numerics import simstats
from .simulation import SimulationAnalyzer, Simulation
from .structure import generate_structure
from .rmg_input import RmgInput, rmg_modes


class RmgAnalyzer(SimulationAnalyzer):

    @property
    def initialized(self):
        return 'path' in self
    #end def initialized

    @property
    def analyzed(self):
        return 'setup_info' in self and 'results' in self
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
        return self.analyzed and 'timing' in self.results
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


    def __init__(self,arg0=None,*,analyze=False):
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
        self.input        = None
        self.info         = obj()

        if analyze:
            self.analyze()
        #end if
    #end def __init__


    def analyze(self,*,guard=True):
        if not self.initialized:
            return
        #end if
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
        self.setup_info = obj()
        self.results    = obj()
        parse_status = obj(log=True,errors=obj())
        self.info.parse_status = parse_status

        def parse(name,method,*args):
            try:
                status = method(*args)
                parse_status[name] = True if status is None else status
            except Exception as e:
                parse_status[name] = False
                parse_status.errors[name] = '{}: {}'.format(e.__class__.__name__,e)
                if not guard:
                    raise
                #end if
            #end try
        #end def parse

        parse('setup',self.read_setup_info,logfile)
        parse('results',self.read_results,logfile,lines)
    #end def analyze


    def read_setup_info(self,logfile):
        setup_info = obj()
        f = logfile
        mode = None
        if f.seek('Calculation type',1) != -1:
            line = f.readline().lower()
            mode = self.identify_mode(line)
            setup_info.run_mode = mode
        #end if
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
            except:
                try:
                    v = float(v)
                except:
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
                            except:
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
            istart = f.seek(setup_start)
            if istart!=-1:
                istart = f.tell()
                iend = f.seek(setup_end,1)
                if iend!=-1:
                    iend = f.tell()
                    text = to_str(f.mm[istart:iend])
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
                        header,body = b.split('\n',1)
                        bname = process_name(header)
                        lines = body.splitlines()
                        simple_values = True
                        for line in lines:
                            simple_values &= ':' in line
                        #end for
                        if simple_values:
                            bvalues = obj()
                            for line in lines:
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
                        try:
                            grid = []
                            grid_pe = []
                            spacing = []
                            grid_units = None
                            for c in 'xyz':
                                if c in b:
                                    s = b[c].replace('Total:','')
                                    s = s.replace('Per PE:','')
                                    s = s.replace('Spacing:','')
                                    v,u = process_value(s,list=True)
                                    grid_units = u
                                    grid.append(v[0])
                                    grid_pe.append(v[1])
                                    spacing.append(v[2])
                                #end if
                            #end for
                            grid = np.array(grid,dtype=int)
                            grid_pe = np.array(grid_pe,dtype=int)
                            spacing = np.array(spacing,dtype=float)
                            ecut,ecut_charge,ecut_units = b.equivalent_energy_cutoffs.split()
                            b.update(
                                grid         = grid,
                                grid_pe      = grid_pe,
                                grid_spacing = spacing,
                                grid_units   = grid_units,
                                ecut         = float(ecut),
                                ecut_charge  = float(ecut_charge),
                                ecut_units   = ecut_units,
                                )
                        except:
                            pass
                        #end try
                    #end if
                    if 'lattice_setup' in setup_info:
                        b = setup_info.lattice_setup
                        try:
                            b.axes = np.array([b.x_basis_vector,b.y_basis_vector,b.z_basis_vector],dtype=float)
                        except:
                            pass
                        #end try
                    #end if
                    if 'k_points' in other_blocks:
                        try:
                            header,body,lines = other_blocks.k_points
                            del other_blocks.k_points
                            for i,line in enumerate(lines):  # noqa: B007
                                if 'Weight in crystal unit' in line:
                                    break
                                #end if
                            #end for
                            kp = []
                            kw = []
                            for line in lines[i+1:]:
                                if 'Weight in' in line:
                                    break
                                #end if
                                t = np.array(line.split(),dtype=float)
                                kp.append(t[:3])
                                kw.append(t[3])
                            #end for
                            setup_info.k_points = obj(
                                kpoints_crystal = np.array(kp,dtype=float),
                                kweights        = np.array(kw,dtype=float),
                                )
                        except:
                            pass
                        #end try
                    #end if
                    k = 'initial_ionic_positions_and_displacements'
                    if k in other_blocks:
                        try:
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
                            for i,line in enumerate(lines):  # noqa: B007
                                if 'Species' in line:
                                    break
                                #end if
                            #end for
                            for line in lines[i+1:]:
                                ls = line.strip()
                                if len(ls)>0:
                                    t = line.split()
                                    spec.append(t[0])
                                    pos.append(t[1:4])
                                #end if
                            #end for
                            setup_info.ion_positions = obj(
                                units     = punits,
                                atoms     = np.array(spec,dtype=object),
                                positions = np.array(pos,dtype=float),
                                )
                        except:
                            pass
                        #end try
                    #end if
                #end if
            #end if
        #end if
        if 'lattice_setup' in setup_info and 'ion_positions' in setup_info:
            try:
                aunits = setup_info.lattice_setup.get('units','B')
                axes   = setup_info.lattice_setup.axes
                elem   = setup_info.ion_positions.atoms
                pos    = setup_info.ion_positions.positions
                punits = setup_info.ion_positions.units
                kpu    = None
                kw     = None
                if aunits=='a0':
                    aunits = 'B'
                elif aunits!='B':
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
            except:
                pass
            #end try
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
                except:
                    pass
                #end try
            #end if
        #end if
        self.setup_info = setup_info
        return len(setup_info)>0
    #end def read_setup_info


    def identify_mode(self,text):
        text = text.lower()
        mode_patterns = (
            ('quench electrons','scf'),
            ('nscf','nscf'),
            ('band structure','band'),
            ('exx integral','exx'),
            ('structure optimization','relax'),
            ('molecular dynamics - cve','md_VE'),
            ('molecular dynamics - cvt','md_TE'),
            ('nudged elastic band','neb'),
            ('time dependent dft','tddft'),
            ('tddft','tddft'),
            ('stm charge density','stm'),
            )
        for pattern,mode in mode_patterns:
            if pattern in text:
                return mode
            #end if
        #end for
        return rmg_modes.mode_match(text,short=True)
    #end def identify_mode


    number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


    def line_numbers(self,line):
        values = re.findall(self.number_pattern,line.replace('D','E').replace('d','e'))
        return np.array(values,dtype=float)
    #end def line_numbers


    def match_float(self,pattern,line):
        match = re.search(pattern,line,re.I)
        if match is None:
            return None
        #end if
        return float(match.group(1).replace('D','E').replace('d','e'))
    #end def match_float


    def read_results(self,logfile,lines=None):
        if lines is None:
            lines = to_str(logfile.mm[:]).splitlines()
        #end if
        self.results = obj()
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
            ('artifacts',self.analyze_artifacts),
            )
        parse_status = self.info.parse_status
        for name,parser in parsers:
            status = parser(lines)
            parse_status[name] = status
        #end for
        return len(self.results)>0
    #end def read_results


    def analyze_energies(self,lines):
        energies = []
        direct_energies = []
        epat = r'final total energy from eig sum\s*=\s*('+self.number_pattern+r')'
        dpat = r'final total energy from direct\s*=\s*('+self.number_pattern+r')'
        for line in lines:
            value = self.match_float(epat,line)
            if value is not None:
                energies.append(value)
            #end if
            value = self.match_float(dpat,line)
            if value is not None:
                direct_energies.append(value)
            #end if
        #end for
        self.energies = np.array(energies,dtype=float)
        self.direct_energies = np.array(direct_energies,dtype=float)
        if len(energies)>0:
            self.E = energies[-1]
            self.results.energy = self.E
            self.results.energy_units = 'Ha'
            self.results.energies = self.energies
        #end if
        if len(direct_energies)>0:
            self.results.direct_energies = self.direct_energies
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
            lower = line.lower()
            if 'fermi energy' in lower:
                value = self.match_float(r'fermi energy\s*=\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.fermi_energies.append(value)
                #end if
            elif 'valence band maximum' in lower and 'conduction band' in lower:
                values = self.line_numbers(line)
                if len(values)>=2:
                    data.valence_band_maxima.append(values[-2])
                    data.conduction_band_minima.append(values[-1])
                #end if
            elif 'band gap' in lower:
                value = self.match_float(r'band gap\s*=\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.band_gaps.append(value)
                #end if
            elif 'total charge in supercell' in lower:
                value = self.match_float(r'=\s*('+self.number_pattern+r')',line)
                if value is not None:
                    data.total_charges.append(value)
                #end if
            elif '@@ total magnetization' in lower:
                values = self.line_numbers(line)
                if len(values)>0:
                    data.total_magnetizations.append(values[-1])
                #end if
            elif '@@ absolute magnetization' in lower:
                values = self.line_numbers(line)
                if len(values)>0:
                    data.absolute_magnetizations.append(values[-1])
                #end if
            elif lower.lstrip().startswith('sum force'):
                values = self.line_numbers(line)
                if len(values)>=3:
                    sum_forces.append(values[-3:])
                #end if
            elif 'volume and energy per atom' in lower:
                values = self.line_numbers(line)
                if len(values)>=2:
                    volume_per_atom.append(values[-2])
                    energy_per_atom.append(values[-1])
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
            self.electronic = data
            self.results.electronic = data
        #end if
        if len(data.fermi_energies)>0:
            self.Ef = data.fermi_energies[-1]
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
                if '@@ '+label.lower() in line.lower():
                    nums = self.line_numbers(line)
                    values[name].append(nums[-1] if len(nums)>0 else np.nan)
                    break
                #end if
            #end for
            if 'quench:' in line:
                match = re.search(
                    r'md:\s*(\d+)/\d+\s+scf:\s*(\d+)/\d+\s+'
                    r'step time:\s*('+self.number_pattern+r')\s+'
                    r'scf time:\s*('+self.number_pattern+r').*?RMS\[dV\]:\s*([^\]\s]+)',line,re.I)
                if match is not None:
                    md_steps.append(int(match.group(1)))
                    scf_steps.append(int(match.group(2)))
                    step_times.append(float(match.group(3)))
                    scf_times.append(float(match.group(4)))
                    rms = match.group(5)
                    rms_dv.append(float(rms) if '*' not in rms else np.nan)
                #end if
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
            self.scf_data = scf
            self.results.scf = scf
        #end if
        return nfound
    #end def analyze_scf


    def analyze_ions(self,lines):
        records = []
        i = 0
        while i<len(lines):
            if not lines[i].lstrip().startswith('@ION  Ion'):
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
            while i<len(lines) and lines[i].lstrip().startswith('@ION'):
                tokens = lines[i].split()
                if len(tokens)<13:
                    break
                #end if
                try:
                    atoms.append(tokens[2])
                    positions.append([float(v) for v in tokens[3:6]])
                    charges.append(float(tokens[6]))
                    magnetizations.append(float(tokens[7]))
                    forces.append([float(v) for v in tokens[8:11]])
                    movable.append([int(v) for v in tokens[11:14]])
                except (ValueError,IndexError):
                    break
                #end try
                i += 1
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
            self.ionic_steps = ionic_steps
            self.results.ionic_steps = ionic_steps
            self.results.position_units = 'a0'
            self.results.force_units = 'Ha/a0'
        #end if
        if len(records)>0 and len(set(len(r.atoms) for r in records))==1:
            self.positions = np.array([r.positions for r in records],dtype=float)
            self.forces = np.array([r.forces for r in records],dtype=float)
            self.charges = np.array([r.charges for r in records],dtype=float)
            self.magnetizations = np.array([r.magnetizations for r in records],dtype=float)
            self.max_forces = np.array([
                np.sqrt((r.forces**2).sum(axis=1)).max() for r in records
                ],dtype=float)
            self.results.positions = self.positions
            self.results.forces = self.forces
            self.results.charges = self.charges
            self.results.magnetizations = self.magnetizations
            self.results.max_forces = self.max_forces
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
                self.structures = structures
                self.results.structures = structures
            #end if
        #end if
        return len(records)
    #end def analyze_ions


    def analyze_md(self,lines):
        records = []
        for line in lines:
            stripped = line.strip()
            if not (stripped.startswith('@CVE') or stripped.startswith('@CVT')):
                continue
            #end if
            values = self.line_numbers(stripped)
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
            self.md_data = md
            self.md_stats = self.md_statistics()
            self.results.md = md
            self.results.md_stats = self.md_stats
        #end if
        return len(records)
    #end def analyze_md


    def md_statistics(self,equil=None):
        statistics = obj()
        if 'md_data' not in self:
            return statistics
        #end if
        for name,values in self.md_data.items():
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
            self.volume = geometry.volume
        #end if
        if 'k_points' in self.setup_info:
            kpoints = self.setup_info.k_points
            geometry.kpoints_crystal = kpoints.kpoints_crystal
            geometry.kweights = kpoints.kweights
            self.kpoints_unit = kpoints.kpoints_crystal
            self.kweights = kpoints.kweights
            if 'structure' in self.setup_info and len(kpoints.kpoints_crystal)>0:
                geometry.kpoints_cart = np.dot(
                    kpoints.kpoints_crystal,self.setup_info.structure.kaxes)
                self.kpoints_cart = geometry.kpoints_cart
            #end if
        #end if
        if len(geometry)>0:
            self.geometry = geometry
            self.results.geometry = geometry
        #end if
        return len(geometry)>0
    #end def analyze_geometry


    def analyze_stress(self,lines):
        tensors = []
        for i,line in enumerate(lines):
            if 'stress total' not in line.lower() or 'kbar' not in line.lower():
                continue
            #end if
            rows = []
            for row_line in lines[i+1:i+4]:
                values = self.line_numbers(row_line)
                if len(values)<3:
                    rows = []
                    break
                #end if
                rows.append(values[:3])
            #end for
            if len(rows)==3:
                tensors.append(np.array(rows,dtype=float))
            #end if
        #end for
        if len(tensors)>0:
            self.stress = np.array(tensors,dtype=float)
            self.pressures = -np.trace(self.stress,axis1=1,axis2=2)/3.0
            self.pressure = self.pressures[-1]
            self.results.stress = self.stress
            self.results.stress_units = 'kbar'
            self.results.pressures = self.pressures
        #end if
        return len(tensors)
    #end def analyze_stress


    def analyze_convergence(self,lines):
        electronic_successes = 0
        electronic_failures = 0
        ionic_converged = None
        for line in lines:
            lower = line.lower()
            if 'potential convergence has been achieved' in lower:
                electronic_successes += 1
            elif 'convergence criterion not met' in lower:
                electronic_failures += 1
            elif 'force convergence has been achieved' in lower:
                ionic_converged = True
            elif 'force convergence has not been achieved' in lower:
                ionic_converged = False
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
        self.convergence = convergence
        self.results.convergence = convergence
        return electronic_successes+electronic_failures+(ionic_converged is not None)
    #end def analyze_convergence


    def analyze_timing(self,lines):
        timing = None
        sections = obj()
        time_value = r'(?:'+self.number_pattern+r'|inf|nan)'
        pattern = re.compile(r'^\s*(.+?)\s{2,}('+time_value+r')\s+('+time_value+r')\s*$',re.I)
        for line in lines:
            match = pattern.match(line)
            if match is None:
                continue
            #end if
            name = match.group(1).strip()
            try:
                total = float(match.group(2))
                per_step = float(match.group(3))
            except ValueError:
                continue
            #end try
            key = re.sub(r'[^a-z0-9]+','_',name.lower()).strip('_')
            sections[key] = obj(total=total,per_step=per_step)
            if name=='1-TOTAL':
                timing = obj(total=total,per_step=per_step,units='s')
            #end if
        #end for
        if timing is not None:
            timing.sections = sections
            self.total_time = timing.total
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
                        values = self.line_numbers(line)
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
            self.bands = bands
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
                    values = self.line_numbers(line)
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
                    lower = line.lower()
                    values = self.line_numbers(line)
                    if 'electric field' in lower and len(values)>=3:
                        field = values[-3:]
                    elif 'dipole at' in lower and len(values)>=3:
                        ground_state = values[-3:]
                    elif not line.lstrip().startswith('&&') and len(values)>=4:
                        rows.append(values[:4])
                    #end if
                #end for
            #end with
            if len(rows)>0:
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
            self.tddft = tddft
            self.results.tddft = tddft
        #end if
        return len(tddft)
    #end def analyze_tddft


    def analyze_artifacts(self,lines):
        artifacts = obj()
        mode = self.calculation_shortmode
        if mode=='exx':
            exx_files = sorted(glob(os.path.join(self.path,'*exx*integral*.h5')))
            if len(exx_files)>0:
                artifacts.exx_integrals = exx_files
            #end if
        #end if
        if mode=='stm':
            stm_files = sorted(glob(os.path.join(self.path,'STM','*.stm')))
            stm_cube_files = sorted(glob(os.path.join(self.path,'STM','*.cube')))
            if len(stm_files)>0:
                artifacts.stm = stm_files
            #end if
            if len(stm_cube_files)>0:
                artifacts.stm_cube = stm_cube_files
            #end if
        elif mode=='scf' and 'files' in self.setup_info:
            data_output = self.setup_info.files.get('data_output_file',None)
            if data_output is not None:
                qmcpack_file = os.path.join(self.path,data_output+'.h5')
                if os.path.isfile(qmcpack_file):
                    artifacts.qmcpack_restart = qmcpack_file
                #end if
            #end if
        #end if
        if len(artifacts)>0:
            self.artifacts = artifacts
            self.results.artifacts = artifacts
        #end if
        return sum(len(v) for v in artifacts.values())
    #end def analyze_artifacts


    def return_initial_structure(self):
        s = None
        if 'setup_info' in self and 'structure' in self.setup_info:
            s = deepcopy(self.setup_info.structure)
        #end if
        return s
    #end def return_initial_structure

#end class RmgAnalyzer
