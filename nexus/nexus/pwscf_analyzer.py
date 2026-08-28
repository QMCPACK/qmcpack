##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf_analyzer.py                                                 #
#    Supports data analysis for PWSCF output.  Can handle log file   #
#    and XML output.                                                 #
#                                                                    #
#  Content summary:                                                  #
#    PwscfAnalyzer                                                   #
#      SimulationAnalyzer class for PWSCF.                           #
#      Reads log output and converts data to numeric form.           #
#      Can also read data-file.xml.  See pwscf_data_reader.py.       #
#                                                                    #
#====================================================================#


import os
import re
import numpy as np
from .developer import obj
from .unit_converter import convert
from .numerics import simstats, simplestats
from .simulation import SimulationAnalyzer, Simulation
from .structure import Structure, get_kpath
from .pwscf_input import PwscfInput
from .pwscf_data_reader import read_qexml
from .utilities import path_string
from . import numpy_extensions as npe


number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


def match_float(pattern,line):
    match = re.search(pattern,line)
    if match is None:
        return None
    return float(match.group(1).replace('D','E').replace('d','e'))
#end def match_float


def read_kpoint_tables(lines):
    kpoints_cart = None
    kpoints_unit = None
    kweights = None
    for i,line in enumerate(lines):
        if 'number of k points=' not in line:
            continue
        match = re.search(r'number of k points=\s*(\d+)',line)
        if match is None:
            continue
        nkpoints = int(match.group(1))
        if i+1>=len(lines) or 'cart. coord.' not in lines[i+1]:
            continue
        cart = []
        weights = []
        valid = len(lines[i+2:i+2+nkpoints])==nkpoints
        for kline in lines[i+2:i+2+nkpoints]:
            match = re.search(r'=\s*\((.*?)\),\s*wk\s*=\s*('+number_pattern+r')',kline)
            if match is None:
                valid = False
                break
            coordinates = np.array(
                re.findall(number_pattern,match.group(1).replace('D','E').replace('d','e')),
                dtype=float,
                )
            if len(coordinates)<3:
                valid = False
                break
            cart.append(coordinates[:3])
            weights.append(float(match.group(2).replace('D','E').replace('d','e')))
        if not valid:
            continue
        j = i+2+nkpoints
        while j<len(lines) and 'cryst. coord.' not in lines[j]:
            j+=1
        if j>=len(lines):
            continue
        unit = []
        valid = len(lines[j+1:j+1+nkpoints])==nkpoints
        for kline in lines[j+1:j+1+nkpoints]:
            match = re.search(r'=\s*\((.*?)\),\s*wk\s*=',kline)
            if match is None:
                valid = False
                break
            coordinates = np.array(
                re.findall(number_pattern,match.group(1).replace('D','E').replace('d','e')),
                dtype=float,
                )
            if len(coordinates)<3:
                valid = False
                break
            unit.append(coordinates[:3])
        if not valid:
            continue
        kpoints_cart = np.array(cart,dtype=float)
        kpoints_unit = np.array(unit,dtype=float)
        kweights = np.array(weights,dtype=float)
        break
    return kpoints_cart,kpoints_unit,kweights
#end def read_kpoint_tables



class PwscfAnalyzer(SimulationAnalyzer):
    def __init__(
        self,
        arg0              = None,
        infile_name       = None,
        outfile_name      = None,
        pw2c_outfile_name = None,
        *,
        analyze           = False,
        xml               = False,
        warn              = False,
        md_only           = False,
        ):
        self.results_out = None
        self.results_xml = None
        self.info = obj(
            warn        = warn,
            md_only     = md_only,
            data_status = obj(
                log                 = False,
                md                  = False,
                fermi_energies      = False,
                scf_conv_energy     = False,
                scf_conv_accuracy   = False,
                relax_energies      = False,
                bands               = False,
                relax_structures    = False,
                pressure            = False,
                volume              = False,
                stress              = False,
                forces              = False,
                total_forces        = False,
                timing              = False,
                kpoints             = False,
                pw2casino           = False,
                xml                 = False,
                ),
            )
        if isinstance(arg0,Simulation):
            sim = arg0
            path = sim.locdir
            infile_name = sim.infile
            outfile_name= sim.outfile
            self.input_structure = sim.system.structure
        elif arg0 is not None:
            path = path_string(arg0)
            if not os.path.exists(path):
                self.error(
                    'path to QE data does not exist\n'
                    f'path provided: {path}'
                    )
            if os.path.isfile(path):
                filepath = path
                path,filename = os.path.split(filepath)
                if filename.endswith('.in'):
                    infile_name = filename
                elif filename.endswith('.out'):
                    outfile_name = filename
                else:
                    self.error(
                        'could not determine whether file is QE input or output\n'
                        f'file provided: {filepath}'
                        )
            if outfile_name is None:
                outfile_name = f"{infile_name.rsplit('.',1)[0]}.out"
        else:
            return

        self.infile_name  = infile_name
        self.outfile_name = outfile_name

        self.path = path
        self.abspath = os.path.abspath(path)
        self.pw2c_outfile_name = pw2c_outfile_name

        self.input = None
        if self.infile_name is not None:
            self.input = PwscfInput(os.path.join(self.path,self.infile_name))
        self.initialize_results()
        if analyze:
            self.analyze(xml=xml)
    #end def __init__


    def initialize_results(self):
        self.results_out = obj()
        calculation = 'scf'
        if self.input is not None and 'control' in self.input and 'calculation' in self.input.control:
            calculation = self.input.control.calculation.lower()

        result_names = (
            'Ef','fermi_energies','bands','volume','cputime','walltime',
            'kpoints_cart','kpoints_unit','kweights','K',
            )
        if calculation not in ('nscf','bands'):
            result_names += (
                'E','relax_energies','scf_conv_energy','scf_conv_accuracy',
                'pressure','stress','forces','tot_forces','max_forces',
                )
        if calculation in ('relax','vc-relax','md','vc-md'):
            result_names += ('relax_structures',)
        if calculation in ('md','vc-md'):
            result_names += ('md_data','md_stats')
        for name in result_names:
            self.results_out[name] = None
    #end def initialize_results

    
    def analyze(self,*,xml=False):
        outfile = os.path.join(self.path,self.outfile_name)
        try:
            with open(outfile,'r') as fobj:
                lines = fobj.read().splitlines()
        except (OSError,UnicodeError) as e:
            if self.info.warn:
                self.warn(
                    'file read failed\n'
                    f'exception encountered: {e}'
                    )
            return
        self.info.data_status.log = True

        self.analyze_md(lines)
        if self.info.md_only:
            return
        self.analyze_fermi_energies(lines)
        self.analyze_scf_convergence(lines)
        self.analyze_energies(lines)
        self.analyze_bands(lines)
        self.analyze_structures(lines)
        self.analyze_pressure_volume(lines)
        self.analyze_stress(lines)
        self.analyze_forces(lines)
        self.analyze_timing(lines)
        self.analyze_kpoints(lines)
        self.analyze_pw2casino()
        if xml:
            self.analyze_xml()
            if self.results_xml is not None:
                self.info.xml_status_failed = bool(self.results_xml.failed)
    #end def analyze


    def analyze_md(self,lines):
        calculation = ''
        if self.input is not None and 'control' in self.input:
            calculation = self.input.control.calculation if 'calculation' in self.input.control else 'scf'
            calculation = calculation.lower()
        if calculation not in ('md','vc-md'):
            return 0

        records = []
        record = obj()
        required = ('total_energy','pressure','time','kinetic_energy','temperature')
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                energy = match_float(r'total energy\s*=\s*('+number_pattern+r')',line)
                record = obj()
                if energy is not None:
                    record.total_energy = energy
            elif 'total   stress' in line and 'P=' in line:
                pressure = match_float(r'P=\s*('+number_pattern+r')',line)
                if pressure is not None and 'total_energy' in record:
                    record.pressure = pressure
            if calculation=='md':
                if 'time      =' in line:
                    time = match_float(r'time\s*=\s*('+number_pattern+r')',line)
                    if time is not None and 'total_energy' in record:
                        record.time = time
                elif 'kinetic energy' in line and '=' in line:
                    kinetic = match_float(r'kinetic energy.*?=\s*('+number_pattern+r')',line)
                    if kinetic is not None and 'total_energy' in record:
                        record.kinetic_energy = kinetic
                elif line.strip().startswith('temperature') and '=' in line:
                    temperature = match_float(r'temperature\s*=\s*('+number_pattern+r')',line)
                    if temperature is not None and 'total_energy' in record:
                        record.temperature = temperature
            else:
                if 'Entering Dynamics;' in line and 'time' in line:
                    time = match_float(r'time\s*=\s*('+number_pattern+r')',line)
                    if time is not None and 'total_energy' in record:
                        record.time = time
                elif line.strip().startswith('Ekin'):
                    match = re.search(
                        r'Ekin\s*=\s*('+number_pattern+r')\s+Ry\s+T\s*=\s*('+number_pattern+r')',
                        line,
                        )
                    if match is not None and 'total_energy' in record:
                        record.kinetic_energy = float(match.group(1).replace('D','E').replace('d','e'))
                        record.temperature = float(match.group(2).replace('D','E').replace('d','e'))
            if all(name in record for name in required):
                records.append(record)
                record = obj()
        nsteps = len(records)
        if nsteps==0:
            return 0
        self.info.data_status.md = True
        md = obj(
            total_energy   = np.array([r.total_energy for r in records],dtype=float),
            pressure       = np.array([r.pressure for r in records],dtype=float),
            time           = np.array([r.time for r in records],dtype=float),
            kinetic_energy = np.array([r.kinetic_energy for r in records],dtype=float),
            temperature    = np.array([r.temperature for r in records],dtype=float),
            )
        md.potential_energy = md.total_energy - md.kinetic_energy
        self.results_out.md_data = md
        self.results_out.md_stats = self.md_statistics()
        return nsteps
    #end def analyze_md


    def analyze_fermi_energies(self,lines):
        fermi_energies = []
        for line in lines:
            if 'Fermi energ' in line:
                ev_index = line.lower().rfind('ev')
                segment = line[:ev_index] if ev_index!=-1 else line
                fermi_energies.extend(np.array(
                    re.findall(number_pattern,segment.replace('D','E').replace('d','e')),
                    dtype=float,
                    ))
        if len(fermi_energies)>0:
            self.info.data_status.fermi_energies = True
            self.results_out.Ef = fermi_energies[-1]
            self.results_out.fermi_energies = np.array(fermi_energies)
        return len(fermi_energies)
    #end def analyze_fermi_energies


    def analyze_energies(self,lines):
        relax_energies = []
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                energy = match_float(r'total energy\s*=\s*('+number_pattern+r')',line)
                if energy is not None:
                    relax_energies.append(energy)
        if len(relax_energies)>0:
            self.info.data_status.relax_energies = True
            self.results_out.E = relax_energies[-1]
            self.results_out.relax_energies = np.array(relax_energies)
        return len(relax_energies)
    #end def analyze_energies


    def analyze_scf_convergence(self,lines):
        scf_conv_energy = []
        scf_conv_accuracy = []
        capture_accuracy = False
        for line in lines:
            if 'total energy' in line and '=' in line:
                capture_accuracy = False
                if not line.lstrip().startswith('!'):
                    energy = match_float(r'total energy\s*=\s*('+number_pattern+r')',line)
                    if energy is not None:
                        scf_conv_energy.append(energy)
                        capture_accuracy = True
            elif capture_accuracy and 'estimated scf accuracy' in line:
                accuracy = match_float(r'estimated scf accuracy\s*[<=>]\s*('+number_pattern+r')',line)
                if accuracy is not None:
                    scf_conv_accuracy.append(accuracy)
                capture_accuracy = False
        if len(scf_conv_energy)>0:
            self.info.data_status.scf_conv_energy = True
            self.results_out.scf_conv_energy = np.array(scf_conv_energy)
        if len(scf_conv_accuracy)>0:
            self.info.data_status.scf_conv_accuracy = True
            self.results_out.scf_conv_accuracy = np.array(scf_conv_accuracy)
        return obj(energy=len(scf_conv_energy),accuracy=len(scf_conv_accuracy))
    #end def analyze_scf_convergence


    def analyze_bands(self,lines):
        def leading_numbers(line):
            pattern = r'^\s*((?:'+number_pattern+r')(?:\s*(?:'+number_pattern+r'))*)'
            match = re.match(pattern,line)
            if match is None:
                return np.array([],dtype=float)
            return np.array(
                re.findall(number_pattern,match.group(1).replace('D','E').replace('d','e')),
                dtype=float,
                )
        #end def leading_numbers

        table_cart,table_unit,_ = read_kpoint_tables(lines)
        nspin = 1
        if self.input is not None and 'system' in self.input and 'nspin' in self.input.system:
            nspin = self.input.system.nspin
        polarized = nspin>1
        bands = obj(up=obj(),down=obj())
        band_channel = bands.up
        up_spin = True
        nfound = 0
        for i,line in enumerate(lines):
            if 'End of self-consistent calculation' in line and nfound>0:
                bands = obj(up=obj(),down=obj())
                band_channel = bands.up
                up_spin = True
                nfound = 0
                continue
            elif '- SPIN UP -' in line:
                up_spin = True
                band_channel = bands.up
                continue
            elif '- SPIN DOWN -' in line:
                up_spin = False
                band_channel = bands.down
                continue
            elif 'bands (ev)' not in line:
                continue

            eigs = []
            j = i+1
            while j<len(lines):
                ls = lines[j].strip()
                if len(ls)==0:
                    if len(eigs)>0:
                        break
                elif 'occupation numbers' in ls or 'bands (ev)' in ls:
                    break
                else:
                    values = leading_numbers(ls)
                    if len(values)==0:
                        break
                    eigs.extend(values)
                j+=1
            while j<len(lines) and len(lines[j].strip())==0:
                j+=1
            occs = []
            if j<len(lines) and 'occupation numbers' in lines[j]:
                j+=1
                while j<len(lines):
                    ls = lines[j].strip()
                    if len(ls)==0:
                        if len(occs)>0:
                            break
                    else:
                        values = leading_numbers(ls)
                        if len(values)==0:
                            break
                        occs.extend(values)
                    j+=1

            match = re.search(r'k\s*=\s*(.*?)\s*\(',line)
            kpoint_cart = None
            if match is not None:
                values = np.array(
                    re.findall(number_pattern,match.group(1).replace('D','E').replace('d','e')),
                    dtype=float,
                    )
                if len(values)>=3:
                    kpoint_cart = values[:3]
            index = len(band_channel)
            kpoint_rel = table_unit[index] if table_unit is not None and index<len(table_unit) else kpoint_cart
            if table_cart is not None and index<len(table_cart):
                kpoint_cart = table_cart[index]
            if len(eigs)==0:
                continue
            pol = 'none'
            if polarized:
                pol = 'up' if up_spin else 'down'
            band_channel[index] = obj(
                index           = index,
                kpoint_2pi_alat = kpoint_cart,
                kpoint_rel      = kpoint_rel,
                eigs            = np.array(eigs,dtype=float),
                occs            = np.array(occs,dtype=float),
                pol             = pol,
                )
            nfound+=1
        if nfound==0:
            return 0
        self.info.data_status.bands = True
        self.analyze_band_edges(bands)
        self.results_out.bands = bands
        return nfound
    #end def analyze_bands


    def analyze_band_edges(self,bands):
        vbm = obj(energy=-1.0e6)
        cbm = obj(energy=1.0e6)
        direct_gap = obj(energy=1.0e6)
        found_occupations = False
        for band_channel in (bands.up,bands.down):
            for b in band_channel.values():
                if len(b.occs)!=len(b.eigs) or len(b.occs)==0:
                    continue
                occ = b.occs > 0.5
                unocc = b.occs < 0.5
                if not occ.any() or not unocc.any():
                    continue
                found_occupations = True
                e_val = np.max(b.eigs[occ])
                e_cond = np.min(b.eigs[unocc])
                if e_val > vbm.energy:
                    vbm = obj(
                        energy          = e_val,
                        kpoint_rel      = b.kpoint_rel,
                        kpoint_2pi_alat = b.kpoint_2pi_alat,
                        index           = b.index,
                        pol             = b.pol,
                        band_number     = np.max(np.where(occ)),
                        )
                if e_cond < cbm.energy:
                    cbm = obj(
                        energy          = e_cond,
                        kpoint_rel      = b.kpoint_rel,
                        kpoint_2pi_alat = b.kpoint_2pi_alat,
                        index           = b.index,
                        pol             = b.pol,
                        band_number     = np.min(np.where(unocc)),
                        )
                if e_cond-e_val < direct_gap.energy:
                    direct_gap = obj(
                        energy          = e_cond-e_val,
                        kpoint_rel      = b.kpoint_rel,
                        kpoint_2pi_alat = b.kpoint_2pi_alat,
                        index           = b.index,
                        pol             = [vbm.pol,cbm.pol],
                        )
        if not found_occupations:
            return
        if vbm.energy+0.025>=cbm.energy:
            electronic_structure = 'metallic' if vbm.band_number==cbm.band_number else 'semi-metal'
        else:
            electronic_structure = 'insulating'
            if (
                vbm.kpoint_rel is not None
                and cbm.kpoint_rel is not None
                and not np.equal(vbm.kpoint_rel,cbm.kpoint_rel).all()
                ):
                bands.indirect_gap = obj(
                    energy  = round(cbm.energy-vbm.energy,3),
                    kpoints = obj(vbm=vbm,cbm=cbm),
                    )
        bands.electronic_structure = electronic_structure
        bands.vbm = vbm
        bands.cbm = cbm
        bands.direct_gap = direct_gap
    #end def analyze_band_edges


    def analyze_structures(self,lines):
        structures = obj()
        conf = None
        i = 0
        while i<len(lines):
            line = lines[i]
            if line.strip().startswith('CELL_PARAMETERS'):
                axes = []
                if i+3<len(lines):
                    for d in range(3):
                        axis_line = lines[i+d+1].replace('D','E').replace('d','e')
                        row = np.array(re.findall(number_pattern,axis_line),dtype=float)
                        if len(row)<3:
                            axes = []
                            break
                        axes.append(row[:3])
                if len(axes)==3:
                    conf = obj()
                    axes = np.array(axes,dtype=float)
                    alat = match_float(r'alat\s*=\s*('+number_pattern+r')',line)
                    if alat is not None:
                        axes *= alat
                    conf.axes = axes
                    i+=3
                else:
                    conf = None
            elif 'ATOMIC_POSITIONS' in line:
                if conf is None:
                    conf = obj()
                atoms = []
                positions = []
                i+=1
                while i<len(lines):
                    tokens = lines[i].split()
                    if len(tokens)<4 or tokens[0].lower()=='end':
                        break
                    coordinates = tokens[1:4]
                    # Check whether every coordinate is a complete numeric value.
                    if not all(re.fullmatch(number_pattern,value.strip()) is not None for value in coordinates):
                        break
                    atoms.append(tokens[0])
                    positions.append(np.array(coordinates,dtype=float))
                    i+=1
                if len(positions)==0:
                    conf = None
                    continue
                conf.atoms = atoms
                conf.positions = np.array(positions)
                if 'crystal' in line.lower() and 'axes' in conf:
                    conf.positions = np.dot(conf.positions,conf.axes)
                structures[len(structures)] = conf
                conf = None
                continue
            i+=1
        if len(structures)>0:
            self.info.data_status.relax_structures = True
            self.results_out.relax_structures = structures
        return len(structures)
    #end def analyze_structures


    def analyze_pressure_volume(self,lines):
        pressure = None
        volume = None
        for line in lines:
            if 'unit-cell volume' in line:
                value = match_float(r'unit-cell volume\s*=\s*('+number_pattern+r')',line)
                if value is not None:
                    volume = value
            if 'total' in line and 'stress' in line and 'P=' in line:
                value = match_float(r'P=\s*('+number_pattern+r')',line)
                if value is not None:
                    pressure = value
        if pressure is not None:
            self.info.data_status.pressure = True
            self.results_out.pressure = pressure
        if volume is not None:
            self.info.data_status.volume = True
            self.results_out.volume = volume
        return obj(pressure=pressure is not None,volume=volume is not None)
    #end def analyze_pressure_volume


    def analyze_stress(self,lines):
        stress = []
        ntensors = 0
        for i,line in enumerate(lines):
            if 'total   stress' in line:
                rows = []
                if i+3<len(lines):
                    for stress_line in lines[i+1:i+4]:
                        stress_line = stress_line.replace('D','E').replace('d','e')
                        row = np.array(re.findall(number_pattern,stress_line),dtype=float)
                        if len(row)<6:
                            rows = []
                            break
                        rows.append(list(row[:6]))
                if len(rows)==3:
                    stress.extend(rows)
                    ntensors+=1
        if ntensors>0:
            self.info.data_status.stress = True
            self.results_out.stress = stress
        return ntensors
    #end def analyze_stress


    def analyze_forces(self,lines):
        forces = []
        tot_forces = []
        nat = None
        if self.input is not None and 'system' in self.input and 'nat' in self.input.system:
            nat = self.input.system.nat
        for i,line in enumerate(lines):
            if 'Forces acting on atoms' in line:
                aforces = []
                j = i+1
                while j<len(lines):
                    match = re.search(
                        r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*'
                        r'('+number_pattern+r')\s+('+number_pattern+r')\s+('+number_pattern+r')',
                        lines[j],
                        )
                    if match is not None:
                        values = [float(value.replace('D','E').replace('d','e')) for value in match.groups()]
                        aforces.append(np.array(values,dtype=float))
                    elif len(aforces)>0:
                        break
                    j+=1
                if len(aforces)>0 and (nat is None or len(aforces)==nat):
                    forces.append(aforces)
            if 'Total force' in line:
                match = re.search(r'Total force\s*=\s*('+number_pattern+r')',line)
                if match is not None:
                    tot_forces.append(float(match.group(1).replace('D','E').replace('d','e')))
        if len(forces)>0:
            self.info.data_status.forces = True
            self.results_out.forces = np.array(forces,dtype=float)
            self.results_out.tot_forces = np.array(tot_forces)
            self.results_out.max_forces = np.array([
                (np.sqrt((force**2).sum(1))).max()
                for force in self.results_out.forces
                ])
        if len(tot_forces)>0:
            self.info.data_status.total_forces = True
        return obj(forces=len(forces),total_forces=len(tot_forces))
    #end def analyze_forces


    def analyze_timing(self,lines):
        def pwscf_time(text):
            scales = dict(h=1.,m=60.,s=3600.)
            time = 0.
            for value,unit in re.findall(r'('+number_pattern+r')\s*([hms])',text):
                time += float(value.replace('D','E').replace('d','e'))/scales[unit]
            return time
        #end def pwscf_time

        for line in lines:
            if 'PWSCF        :' in line:
                match = re.search(r'PWSCF\s*:\s*(.*?)\s+CPU\s+(.*?)\s+WALL',line)
                if match is not None:
                    self.info.data_status.timing = True
                    self.results_out.cputime = pwscf_time(match.group(1))
                    self.results_out.walltime = pwscf_time(match.group(2))
                    return True
        return False
    #end def analyze_timing


    def analyze_kpoints(self,lines):
        kpoints_cart,kpoints_unit,kweights = read_kpoint_tables(lines)
        if kpoints_cart is not None:
            self.info.data_status.kpoints = True
            self.results_out.kpoints_cart = kpoints_cart
            self.results_out.kpoints_unit = kpoints_unit
            self.results_out.kweights = kweights
            return len(kpoints_cart)
        return 0
    #end def analyze_kpoints


    def analyze_pw2casino(self):
        if self.pw2c_outfile_name is None:
            return False
        filepath = os.path.join(self.path,self.pw2c_outfile_name)
        try:
            with open(filepath,'r') as fobj:
                lines = fobj.readlines()
        except (OSError,UnicodeError) as e:
            if self.info.warn:
                self.warn(
                    'pw2casino file read failed\n'
                    f'exception encountered: {e}'
                    )
            return False
        for line in lines:
            if 'Kinetic' in line:
                tokens = line.split()
                # Check whether the kinetic energy token is a complete numeric value.
                if len(tokens)>5 and re.fullmatch(number_pattern,tokens[5].strip()) is not None:
                    self.info.data_status.pw2casino = True
                    self.results_out.K = float(tokens[5].replace('D','E').replace('d','e'))
                    return True
        return False
    #end def analyze_pw2casino


    def analyze_xml(self):
        import xml.etree.ElementTree as ET

        self.results_xml = None
        if 'xml_status_failed' in self.info:
            del self.info.xml_status_failed
        if self.input is None or 'control' not in self.input:
            if self.info.warn:
                self.warn('xml data is not available\nreason: input control section is not available')
            return 0
        cont = self.input.control
        if 'outdir' not in cont or 'prefix' not in cont:
            if self.info.warn:
                self.warn('xml data is not available\nreason: input outdir/prefix is not available')
            return 0
        savedir = os.path.join(self.path,cont.outdir,f'{cont.prefix}.save')
        schema_file = os.path.join(savedir,'data-file-schema.xml')
        legacy_file = os.path.join(savedir,'data-file.xml')
        legacy_dir = savedir
        if os.path.exists(schema_file):
            self.results_xml = obj(data=None,kpoints=None,failed=False)
            try:
                root = ET.parse(schema_file).getroot()
            except (OSError,ET.ParseError) as e:
                self.xml_failed(
                    'encountered an exception during xml read, this data will not be available\n'
                    f'exception encountered: {e}'
                    )
                return 0
            npoints = self.analyze_schema_xml(root)
        else:
            if not os.path.exists(legacy_file):
                legacy_dir = os.path.join(self.path,cont.outdir)
                legacy_file = os.path.join(legacy_dir,f'{cont.prefix}.xml')
            if not os.path.exists(legacy_file):
                if self.info.warn:
                    self.warn(
                        'xml data is not available\n'
                        f'file not found: {legacy_file}'
                        )
                return 0
            self.results_xml = obj(data=None,kpoints=None,failed=False)
            try:
                data = read_qexml(legacy_file)
            except Exception as e:
                self.xml_failed(
                    'encountered an exception during xml read, this data will not be available\n'
                    f'exception encountered: {e}'
                    )
                return 0
            npoints = self.analyze_legacy_xml(data,legacy_dir)
        if self.results_xml.data is not None:
            self.info.data_status.xml = True
        return npoints
    #end def analyze_xml


    def xml_failed(self,message):
        self.results_xml.failed = True
        if self.info.warn:
            self.warn(message)
    #end def xml_failed


    def analyze_legacy_xml(self,data,datadir):
        def object_path(value,*names):
            for name in names:
                if value is None or name not in value:
                    return None
                value = value[name]
            return value
        #end def object_path

        kpdata = object_path(data,'root','eigenvalues','k_point')
        if kpdata is None:
            self.results_xml.update(data=data,kpoints=obj())
            self.xml_failed(
                'xml data is incomplete, some data will not be available\n'
                'reason: legacy eigenvalue k-points are not available'
                )
            return 0
        kpoints = obj()
        for ki,kpd in kpdata.items():
            if 'k_point_coords' not in kpd or 'weight' not in kpd or 'datafile' not in kpd:
                self.results_xml.failed = True
                continue
            kp = obj(kpoint=kpd.k_point_coords,weight=kpd.weight)
            kpoints[ki] = kp
            for si,dfile in kpd.datafile.items():
                efilepath = os.path.join(datadir,dfile.iotk_link)
                if not os.path.exists(efilepath):
                    self.results_xml.failed = True
                    continue
                try:
                    edata = read_qexml(efilepath)
                except Exception as e:
                    self.xml_failed(
                        'encountered an exception during xml read, this data will not be available\n'
                        f'exception encountered: {e}'
                        )
                    continue
                eunits = object_path(edata,'root','units_for_energies','units')
                eigenvalues = object_path(edata,'root','eigenvalues')
                occupations = object_path(edata,'root','occupations')
                if eunits is None or eigenvalues is None or occupations is None:
                    self.results_xml.failed = True
                    continue
                eunits = eunits.lower()
                if eunits.startswith('ha'):
                    units = 'Ha'
                elif eunits.startswith('ry'):
                    units = 'Ry'
                elif eunits.startswith('ev'):
                    units = 'eV'
                else:
                    units = 'Ha'
                spin = obj(units=units,eigenvalues=eigenvalues,occupations=occupations)
                if si==1:
                    kp.up = spin
                elif si==2:
                    kp.down = spin
        self.results_xml.update(data=data,kpoints=kpoints)
        return len(kpoints)
    #end def analyze_legacy_xml


    def analyze_schema_xml(self,root):
        def xml_value(text):
            text = text.strip()
            if len(text)==0:
                return ''
            elif text.lower() in ('true','false'):
                return text.lower()=='true'
            tokens = text.replace('D','E').replace('d','e').split()
            # Check whether every token is a complete numeric value.
            if all(re.fullmatch(number_pattern,token.strip()) is not None for token in tokens):
                values = np.array(tokens,dtype=float)
                if len(values)==1:
                    value = values[0]
                    if value.is_integer() and all(c not in text.lower() for c in ('.','e')):
                        return int(value)
                    return value
                return values
            return text
        #end def xml_value

        def xml_element(element):
            node = obj()
            for name,value in element.attrib.items():
                node[name.lower()] = xml_value(value)
            groups = obj()
            for child in element:
                name = child.tag.rsplit('}',1)[-1].lower()
                if name not in groups:
                    groups[name] = []
                groups[name].append(child)
            for name,group in groups.items():
                if len(group)==1:
                    node[name] = xml_element(group[0])
                else:
                    node[name] = obj({i+1:xml_element(child) for i,child in enumerate(group)})
            text = (element.text or '').strip()
            if len(text)>0:
                value = xml_value(text)
                if len(node)==0:
                    return value
                node.value = value
            return node
        #end def xml_element

        def xml_child(element,name):
            if element is not None:
                name = name.lower()
                for child in element:
                    if child.tag.rsplit('}',1)[-1].lower()==name:
                        return child
            return None
        #end def xml_child

        data = obj(root=xml_element(root))
        output = xml_child(root,'output')
        band_structure = xml_child(output,'band_structure')
        kpoints = obj()
        if output is None or band_structure is None:
            self.results_xml.update(data=data,kpoints=kpoints)
            self.xml_failed(
                'xml data is incomplete, some data will not be available\n'
                'reason: output band_structure is not available'
                )
            return 0
        lsda_element = xml_child(band_structure,'lsda')
        lsda = xml_value(lsda_element.text) if lsda_element is not None else False
        records = [child for child in band_structure if child.tag.rsplit('}',1)[-1].lower()=='ks_energies']
        if len(records)==0:
            self.results_xml.update(data=data,kpoints=kpoints)
            self.xml_failed(
                'xml data is incomplete, some data will not be available\n'
                'reason: ks_energies records are not available'
                )
            return 0
        coordinate_map = {}
        for record in records:
            k_element = xml_child(record,'k_point')
            e_element = xml_child(record,'eigenvalues')
            o_element = xml_child(record,'occupations')
            if k_element is None or e_element is None or o_element is None:
                self.results_xml.failed = True
                continue
            text = (k_element.text or '').replace('D','E').replace('d','e')
            coordinates = np.array(re.findall(number_pattern,text),dtype=float)
            text = (e_element.text or '').replace('D','E').replace('d','e')
            eigenvalues = np.array(re.findall(number_pattern,text),dtype=float)
            text = (o_element.text or '').replace('D','E').replace('d','e')
            occupations = np.array(re.findall(number_pattern,text),dtype=float)
            if (
                len(coordinates)<3
                or len(eigenvalues)==0
                or len(occupations)==0
                or len(eigenvalues)!=len(occupations)
                ):
                self.results_xml.failed = True
                continue
            coordinates = coordinates[:3]
            weight_text = k_element.attrib.get('weight','')
            # Check whether the weight text is a complete numeric value.
            if re.fullmatch(number_pattern,weight_text.strip()) is None:
                self.results_xml.failed = True
                continue
            weight = float(weight_text.replace('D','E').replace('d','e'))
            key = tuple(np.round(coordinates,12))
            if not lsda or key not in coordinate_map:
                ki = len(kpoints)+1
                coordinate_map[key] = ki
                kpoints[ki] = obj(kpoint=coordinates,weight=weight)
            else:
                ki = coordinate_map[key]
            kp = kpoints[ki]
            spin = obj(
                units       = 'Ha',
                eigenvalues = eigenvalues,
                occupations = occupations,
                )
            if not lsda or 'up' not in kp:
                kp.up = spin
            else:
                kp.down = spin
        self.results_xml.update(data=data,kpoints=kpoints)
        if lsda and any('down' not in kp for kp in kpoints.values()):
            self.results_xml.failed = True
        if len(records)>0 and len(kpoints)==0:
            self.xml_failed(
                'xml data is incomplete, some data will not be available\n'
                'reason: no complete ks_energies records are available'
                )
        return len(kpoints)
    #end def analyze_schema_xml


    def write_electron_counts(
        self,
        filepath = None,
        *,
        return_flag = False,
        ):
        results_xml = self.results_xml
        if not return_flag:
            if results_xml is None:
                self.error('xml data has not been processed\ncannot write electron counts')
            elif results_xml.failed:
                self.error('xml data processing failed\ncannot write electron counts')
        elif results_xml is None or results_xml.failed:
            return False
        kpoints = results_xml.kpoints
        if 'down' in kpoints[1]:
            spins = obj(up='up',down='down')
        else:
            spins = obj(up='up',down='up')
        tot = obj(up=0,down=0)
        for kp in kpoints.values():
            w = kp.weight
            for s,sl in spins.items():
                tot[s] += w*kp[sl].occupations.sum()
        text = 'total electron counts\n'
        total_count = tot.up+tot.down
        spin_count = tot.up-tot.down
        text += f'  {total_count: 3.2f}  {spin_count: 3.2f}  {tot.up: 3.2f}  {tot.down: 3.2f}\n'
        text += '\nkpoint electron counts\n'
        weights = np.array([kp.weight for kp in kpoints.values()],dtype=float)
        mult = (weights/weights.min()).sum()
        for ik in sorted(kpoints.keys()):
            kp = kpoints[ik]
            kpt = obj()
            for s,sl in spins.items():
                kpt[s] = kp.weight*kp[sl].occupations.sum()*mult
            total_count = kpt.up+kpt.down
            spin_count = kpt.up-kpt.down
            text += (
                f'  {ik:>3}  {kp.weight: 8.6f}    {total_count: 3.2f}  '
                f'{spin_count: 3.2f}  {kpt.up: 3.2f}  {kpt.down: 3.2f}\n'
                )
        if filepath is not None:
            with open(filepath,'w') as fobj:
                fobj.write(text)
        if not return_flag:
            return text
        else:
            return True
    #end def write_electron_counts


    def md_statistics(self,equil=None,autocorr=None):
        mds = obj()
        for q,v in self.results_out.md_data.items():
            if equil is not None:
                v = v[equil:]
            if autocorr is None:
                mean,_,error,_ = simstats(v)
            else:
                nv = len(v)
                nb = int(np.floor(float(nv)/autocorr))
                nexclude = nv-nb*autocorr
                v = v[nexclude:]
                npe.reshape_inplace(v, (nb, autocorr))
                mean,error = simplestats(v.mean(axis=1))
            mds[q] = mean,error
        return mds
    #end def md_statistics


    def md_plots(self,*,show=True):

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
        structures = self.results_out.get('relax_structures')
        if structures is None:
            return
        if filepath is None:
            filepath = self.abspath
        filepath = os.path.join(filepath,filename)

        movie = ''
        alat_angstrom = convert(self.input.system['celldm(1)'],'B','A')
        cell = self.input.cell_parameters.vectors
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
        show_vbm_cbm = True,
        k_labels      = None,
        ):
        import matplotlib.pyplot as plt
        bands = self.results_out.bands
        params = {
            'legend.fontsize'      : 14,
            'figure.facecolor'     : 'white',
            'figure.subplot.hspace': 0.,
            'axes.labelsize'       : 16,
            'xtick.labelsize'      : 14,
            'ytick.labelsize'      : 14,
            }
        plt.rcParams.update(params)
        if bands is not None:
            if filename is None:
                filename = 'band_structure.pdf'
            if filepath is None:
                filepath = os.path.join(self.abspath,filename)
            else:
                filepath = os.path.join(filepath,filename)

            fig    = plt.figure()
            ax     = plt.gca()
            nbands = self.input.system.nbnd

            if k_labels is None:
                kpath  = get_kpath(structure=self.input_structure, check_standard=False)
                x      = kpath['explicit_path_linearcoords']
                labels = kpath['explicit_kpoints_labels']
            else:
                labels = k_labels
                # Calculate linear coordinates from self.results_out.kpoints_cart
                x = []
                prev_label = ''
                ref_kpt = self.results_out.kpoints_cart[0]
                lincoord = 0.0
                for kpt_idx,kpt in enumerate(self.results_out.kpoints_cart):
                    curr_label = labels[kpt_idx]
                    if (curr_label != '' and prev_label == '') or curr_label == '':
                        lincoord+=np.linalg.norm(kpt-ref_kpt)
                        ref_kpt = kpt
                    else:
                        ref_kpt = kpt
                        lincoord+=np.linalg.norm(kpt-ref_kpt)
                    x.append(lincoord)
                    prev_label = curr_label
            for nb in range(nbands):
                y = np.array([band.eigs[nb] for band in bands.up]) - bands.vbm.energy
                plt.plot(x, y, 'k')
                if len(bands.down) > 0:
                    y = np.array([band.eigs[nb] for band in bands.down]) - bands.vbm.energy
                    plt.plot(x, y, 'r')
            for ln, li in enumerate(labels):
                if li != '':
                    plt.axvline(x[ln], ymin=-100, ymax=100, linewidth=3, color='k')
                    if li == 'GAMMA':
                        labels[ln] = r'$\Gamma$'
                    else:
                        labels[ln] = f'${li}$'
                    if ln>0 and labels[ln-1] != '':
                        labels[ln] = f'{labels[ln-1]}|{labels[ln]}'
                        labels[ln-1] = ''
            
            plt.xlim([np.min(x), np.max(x)])
            if max_min_e is None:
                plt.ylim(-5, +5)
            else:
                plt.ylim(max_min_e[0],max_min_e[1])
            plt.ylabel('Energy (eV)')
            plt.xticks(x, labels)
            ax.tick_params(axis='x', which='both', length=0)
            ax.tick_params(axis='x', which='both', pad=10)
        if show_vbm_cbm and bands is not None:
            vbm = bands.vbm
            cbm = bands.cbm
            for kn, ki in enumerate(bands.up):
                if (vbm.kpoint_rel == ki.kpoint_rel).all():
                    plt.scatter(x[kn], 0, c='green', s=100)
                if (cbm.kpoint_rel == ki.kpoint_rel).all():
                    plt.scatter(x[kn], cbm.energy-vbm.energy, c='r', s=100)
        if save:
            plt.savefig(filename, format='pdf',bbox_inches='tight')
        if show:
            plt.show()
    #end def plot_bandstructure

#end class PwscfAnalyzer
        
