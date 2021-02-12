##################################################################
##  (c) Copyright 2020-  by Jaron T. Krogel                     ##
##################################################################

import os
import numpy as np

from generic import obj
from developer import to_str
from fileio import TextFile
from unit_converter import convert
from simulation import SimulationAnalyzer,Simulation
from structure import generate_structure
from rmg_input import RmgInput,rmg_modes



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
        return self.analyzed and len(self.setup_info)>0 and len(self.results)>0
    #end def analysis_succeeded

    @property
    def run_completed(self):
        return self.analyzed and 'timing' in self.results
    #end def run_completed

    @property
    def calculation_mode(self):
        mode = None
        return mode
    #end def calculation_mode

    @property
    def calculation_shortmode(self):
        mode = None
        return mode
    #end def calculation_shortmode


    def __init__(self,arg0=None,analyze=False):
        if arg0 is None:
            return
        elif isinstance(arg0,Simulation):
            sim = arg0
            path     = sim.locdir
            filename = sim.infile
        else:
            log_file = arg0
            if not isinstance(log_file,str):
                self.error('invalid type provided for log_file\nType expected: str\nType provided: {}'.format(log_file.__class__.__name__))
            elif not os.path.exists(log_file):
                self.error('RMG log output file does not exist.\nPath provided: {}'.format(log_file))
            elif not os.path.isfile(log_file):
                self.error('Path provided for RMG log output is not a file.\nPath provided: {}'.format(log_file))
            #end if
            path,filename = os.path.split(log_file)
        #end if
        
        self.path         = path
        self.outfile_name = filename
        self.input        = None

        if analyze:
            self.analyze()
        #end if
    #end def __init__


    def analyze(self,guard=True):
        if not self.initialized:
            return
        #end if
        log_filepath = os.path.join(self.path,self.outfile_name)
        if not os.path.exists(log_filepath):
            self.error('RMG analysis cannot be completed.\nLog file does not exist at path provided.\nPath provided: {}'.format(log_filepath))
        #end if
        logfile = TextFile(log_filepath)
        self.setup_info = obj()
        self.results    = obj()
        if guard:
            try:
                self.read_setup_info(logfile)
            except:
                None
            #end try
            try:
                self.read_results(logfile)
            except:
                None
            #end try
        else:
            self.read_setup_info(logfile)
            self.read_results(logfile)
        #end if
    #end def analyze


    def read_setup_info(self,logfile):
        setup_info = obj()
        f = logfile
        mode = None
        if f.seek('Calculation type',1) != -1:
            line = f.readline().lower()
            if 'quench electrons' in line:
                mode = 'scf'
            elif 'band structure' in line:
                mode = 'band'
            else:
                mode = rmg_modes.mode_match(line,short=True)
            #end if
            setup_info.run_mode = mode
        #end if
        setup_start = None
        setup_end   = None
        if mode=='scf':
            setup_start = 'Files'
            setup_end   = 'Diagonalization using'
        elif mode=='band':
            setup_start = 'Files'
            setup_end   = 'converged in'
        else:
            # don't know how to handle other cases yet
            None
        #end if
        unit_set = set(['a0'])
        on_off = dict(ON=True,OFF=False)
        def process_name(s):
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
        def process_value(v,list=False):
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
                            b.set(
                                grid         = grid,
                                grid_pe      = grid_pe,
                                grid_spacing = spacing,
                                grid_units   = grid_units,
                                ecut         = float(ecut),
                                ecut_charge  = float(ecut_charge),
                                ecut_units   = ecut_units,
                                )
                        except:
                            None
                        #end try
                    #end if
                    if 'lattice_setup' in setup_info:
                        b = setup_info.lattice_setup
                        try:
                            b.axes = np.array([b.x_basis_vector,b.y_basis_vector,b.z_basis_vector],dtype=float)
                        except:
                            None
                        #end try
                    #end if
                    if 'k_points' in other_blocks:
                        try:
                            header,body,lines = other_blocks.k_points
                            del other_blocks.k_points
                            for i,line in enumerate(lines):
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
                            None
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
                            for i,line in enumerate(lines):
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
                            None
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
                None
            #end try
        #end if
        if 'files' in setup_info and 'control_input_file' in setup_info.files:
            filepath = os.path.join(self.path,setup_info.files.control_input_file)
            if os.path.exists(filepath):
                try:
                    self.input = RmgInput(filepath)
                except:
                    None
                #end try
            #end if
        #end if
        self.setup_info = setup_info
    #end def read_setup_info


    def read_results(self,logfile):
        results = obj()
        if 'setup_info' in self and 'run_mode' in self.setup_info:
            mode = self.setup_info.run_mode
        else:
            return
        #end if
        f = logfile
        if mode=='scf':
            f.seek('final total energy',1)
            t = f.readtokens()
            results.energy = float(t[-2])
            results.energy_units = t[-1]
        elif mode=='band':
            None
        else:
            self.warn('Results not read.\nUnrecognized run mode: {}'.format(mode))
        #end if
        self.results = results
    #end def read_results


    def return_initial_structure(self):
        s = None
        if 'setup_info' in self and 'structure' in self.setup_info:
            s = self.setup_info.structure.copy()
        #end if
        return s
    #end def return_initial_structure

#end class RmgAnalyzer
