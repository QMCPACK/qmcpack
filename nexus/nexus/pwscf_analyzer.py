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
from .developer import obj,dotdict,FileFormatError
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
    """Return the first captured floating-point value matching a pattern."""
    match = re.search(pattern,line)
    if match is None:
        return None
    return float(match.group(1).replace('D','E').replace('d','e'))
#end def match_float


def read_kpoint_tables(lines):
    """Read Cartesian and crystal k-point tables from PWscf output lines."""
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
        unit  = []
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
        cart    = np.array(cart,dtype=float)
        unit    = np.array(unit,dtype=float)
        weights = np.array(weights,dtype=float)
        return cart,unit,weights
    return None,None,None
#end def read_kpoint_tables



class PwscfAnalyzer(SimulationAnalyzer):
    """Analyze output produced by Quantum ESPRESSO PWscf calculations.

    The analyzer gathers physical results from text and XML output across
    common PWscf run modes. It also provides summaries and visualizations
    useful for inspecting molecular dynamics and electronic structure.
    """

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
        """Initialize a PWscf output analyzer.

        Parameters
        ----------
        arg0 : Simulation or path-like, optional
            PWscf simulation to analyze, or path to a calculation directory,
            input file, or output file. If omitted, only the basic analyzer
            state is initialized.
        infile_name : str, optional
            Name of the PWscf input file within the calculation directory.
        outfile_name : str, optional
            Name of the PWscf output file. It is inferred from ``infile_name``
            when possible.
        pw2c_outfile_name : str, optional
            Name of an accompanying PW2CASINO output file.
        analyze : bool, optional
            Analyze the available output during initialization.
        xml : bool, optional
            Include XML output when analysis is requested.
        warn : bool, optional
            Issue warnings when output is missing, malformed, or incomplete.
        md_only : bool, optional
            Limit analysis to molecular-dynamics data.
        """
        self.results_out = None
        self.results_xml = None
        status_names     = (
            'log','md','fermi_energies','scf_conv_energy',
            'scf_conv_accuracy','relax_energies','bands',
            'relax_structures','pressure','volume','stress','forces',
            'total_forces','timing','kpoints','pw2casino','xml',
            )
        self.info = obj(
            warn        = warn,
            md_only     = md_only,
            data_status = obj({name:False for name in status_names}),
            )
        if isinstance(arg0,Simulation):
            sim                  = arg0
            path                 = sim.locdir
            infile_name          = sim.infile
            outfile_name         = sim.outfile
            self.input_structure = sim.system.structure
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
            if outfile_name is None:
                outfile_name = f"{infile_name.rsplit('.',1)[0]}.out"
        else:
            return

        inp = None
        if infile_name is not None:
            inp = PwscfInput(os.path.join(path,infile_name))

        self.infile_name       = infile_name
        self.outfile_name      = outfile_name
        self.path              = path
        self.abspath           = os.path.abspath(path)
        self.pw2c_outfile_name = pw2c_outfile_name
        self.input             = inp
        self.initialize_results()
        if analyze:
            self.analyze(xml=xml)
    #end def __init__


    def initialize_results(self):
        """Initialize result fields appropriate to the PWscf run mode."""
        calculation = 'scf'
        if (
            self.input is not None
            and 'control' in self.input
            and 'calculation' in self.input.control
            ):
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
        self.results_out = obj({name:None for name in result_names})
    #end def initialize_results

    
    def analyze(self,*,xml=False):
        """Analyze the available PWscf text and optional XML output."""
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
        """Parse molecular-dynamics histories from output lines."""
        def record_value(name,pattern,line):
            """Add a matched value to the current dynamics record."""
            value = match_float(pattern,line)
            if value is not None and 'total_energy' in record:
                record[name] = value
        #end def record_value

        calculation = 'scf'
        if self.input is not None and 'control' in self.input:
            if 'calculation' in self.input.control:
                calculation = self.input.control.calculation.lower()
        if calculation not in ('md','vc-md'):
            return 0

        records  = []
        record   = dotdict()
        required = ('total_energy','pressure','time','kinetic_energy','temperature')
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                energy = match_float(r'total energy\s*=\s*('+number_pattern+r')',line)
                record = dotdict()
                if energy is not None:
                    record.total_energy = energy
            elif 'total   stress' in line and 'P=' in line:
                record_value('pressure',r'P=\s*('+number_pattern+r')',line)
            if calculation=='md':
                if 'time      =' in line:
                    record_value('time',r'time\s*=\s*('+number_pattern+r')',line)
                elif 'kinetic energy' in line and '=' in line:
                    record_value(
                        'kinetic_energy',
                        r'kinetic energy.*?=\s*('+number_pattern+r')',
                        line,
                        )
                elif line.strip().startswith('temperature') and '=' in line:
                    record_value(
                        'temperature',
                        r'temperature\s*=\s*('+number_pattern+r')',
                        line,
                        )
            else:
                if 'Entering Dynamics;' in line and 'time' in line:
                    record_value('time',r'time\s*=\s*('+number_pattern+r')',line)
                elif line.strip().startswith('Ekin'):
                    match = re.search(
                        r'Ekin\s*=\s*('+number_pattern+r')\s+Ry\s+T\s*=\s*('+number_pattern+r')',
                        line,
                        )
                    if match is not None and 'total_energy' in record:
                        record.kinetic_energy = float(
                            match.group(1).replace('D','E').replace('d','e')
                            )
                        record.temperature    = float(
                            match.group(2).replace('D','E').replace('d','e')
                            )
            if all(name in record for name in required):
                records.append(record)
                record = dotdict()
        if len(records)==0:
            return 0
        self.info.data_status.md = True
        md                       = obj({
            name:np.array([record[name] for record in records],dtype=float)
            for name in required
            })
        md.potential_energy       = md.total_energy - md.kinetic_energy
        self.results_out.md_data  = md
        self.results_out.md_stats = self.md_statistics()
        return len(records)
    #end def analyze_md


    def analyze_fermi_energies(self,lines):
        """Parse the sequence of reported Fermi energies."""
        fermi_energies = []
        for line in lines:
            if 'Fermi energ' in line:
                ev_index = line.lower().rfind('ev')
                segment  = line[:ev_index] if ev_index!=-1 else line
                values   = re.findall(
                    number_pattern,
                    segment.replace('D','E').replace('d','e'),
                    )
                fermi_energies.extend(float(value) for value in values)
        if len(fermi_energies)>0:
            self.info.data_status.fermi_energies = True
            self.results_out.Ef                  = fermi_energies[-1]
            self.results_out.fermi_energies      = np.array(fermi_energies,dtype=float)
        return len(fermi_energies)
    #end def analyze_fermi_energies


    def analyze_energies(self,lines):
        """Parse total energies reported for completed SCF calculations."""
        relax_energies = []
        for line in lines:
            if line.lstrip().startswith('!') and 'total energy' in line:
                energy = match_float(r'total energy\s*=\s*('+number_pattern+r')',line)
                if energy is not None:
                    relax_energies.append(energy)
        if len(relax_energies)>0:
            self.info.data_status.relax_energies = True
            self.results_out.E                   = relax_energies[-1]
            self.results_out.relax_energies      = np.array(relax_energies,dtype=float)
        return len(relax_energies)
    #end def analyze_energies


    def analyze_scf_convergence(self,lines):
        """Parse SCF iteration energies and estimated accuracies."""
        scf_conv_energy   = []
        scf_conv_accuracy = []
        capture_accuracy  = False
        for line in lines:
            if 'total energy' in line and '=' in line:
                capture_accuracy = False
                if not line.lstrip().startswith('!'):
                    energy = match_float(r'total energy\s*=\s*('+number_pattern+r')',line)
                    if energy is not None:
                        scf_conv_energy.append(energy)
                        capture_accuracy = True
            elif capture_accuracy and 'estimated scf accuracy' in line:
                accuracy = match_float(
                    r'estimated scf accuracy\s*[<=>]\s*('+number_pattern+r')',
                    line,
                    )
                if accuracy is not None:
                    scf_conv_accuracy.append(accuracy)
                capture_accuracy = False
        if len(scf_conv_energy)>0:
            self.info.data_status.scf_conv_energy = True
            self.results_out.scf_conv_energy      = np.array(scf_conv_energy,dtype=float)
        if len(scf_conv_accuracy)>0:
            self.info.data_status.scf_conv_accuracy = True
            self.results_out.scf_conv_accuracy      = np.array(scf_conv_accuracy,dtype=float)
        return dotdict(energy=len(scf_conv_energy),accuracy=len(scf_conv_accuracy))
    #end def analyze_scf_convergence


    def analyze_bands(self,lines):
        """Parse band energies and occupations for each reported k-point."""
        def leading_numbers(line):
            """Return the leading sequence of numeric values from a line."""
            pattern = r'^\s*((?:'+number_pattern+r')(?:\s*(?:'+number_pattern+r'))*)'
            match   = re.match(pattern,line)
            if match is None:
                return np.array([],dtype=float)
            return np.array(
                re.findall(number_pattern,match.group(1).replace('D','E').replace('d','e')),
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

        table_cart,table_unit,_ = read_kpoint_tables(lines)
        nspin     = 1
        if (
            self.input is not None
            and 'system' in self.input
            and 'nspin' in self.input.system
            ):
            nspin = self.input.system.nspin
        polarized    = nspin>1
        bands        = obj(up=obj(),down=obj())
        band_channel = bands.up
        up_spin      = True
        for i,line in enumerate(lines):
            if ('End of self-consistent calculation' in line
                and len(bands.up)+len(bands.down)>0
                ):
                bands        = obj(up=obj(),down=obj())
                band_channel = bands.up
                up_spin      = True
                continue
            elif '- SPIN UP -' in line:
                up_spin      = True
                band_channel = bands.up
                continue
            elif '- SPIN DOWN -' in line:
                up_spin      = False
                band_channel = bands.down
                continue
            elif 'bands (ev)' not in line:
                continue

            eigs,j = read_values(i+1,('occupation numbers','bands (ev)'))
            if len(eigs)==0:
                continue
            while j<len(lines) and len(lines[j].strip())==0:
                j+=1
            occs = []
            if j<len(lines) and 'occupation numbers' in lines[j]:
                occs,_ = read_values(j+1)

            match       = re.search(r'k\s*=\s*(.*?)\s*\(',line)
            kpoint_cart = None
            if match is not None:
                values = np.array(
                    re.findall(number_pattern,match.group(1).replace('D','E').replace('d','e')),
                    dtype=float,
                    )
                if len(values)>=3:
                    kpoint_cart = values[:3]
            index      = len(band_channel)
            kpoint_rel = kpoint_cart
            if table_unit is not None and index<len(table_unit):
                kpoint_rel = table_unit[index]
            if table_cart is not None and index<len(table_cart):
                kpoint_cart = table_cart[index]
            pol = ('up' if up_spin else 'down') if polarized else 'none'
            band_channel[index] = obj(
                index           = index,
                kpoint_2pi_alat = kpoint_cart,
                kpoint_rel      = kpoint_rel,
                eigs            = np.array(eigs,dtype=float),
                occs            = np.array(occs,dtype=float),
                pol             = pol,
                )
        nfound = len(bands.up)+len(bands.down)
        if nfound==0:
            return 0
        self.info.data_status.bands = True
        self.analyze_band_edges(bands)
        self.results_out.bands = bands
        return nfound
    #end def analyze_bands


    def analyze_band_edges(self,bands):
        """Identify band edges, gaps, and the electronic structure type."""
        def edge_data(band,energy,band_number):
            """Build a band-edge record for one k-point."""
            return obj(
                energy          = energy,
                kpoint_rel      = band.kpoint_rel,
                kpoint_2pi_alat = band.kpoint_2pi_alat,
                index           = band.index,
                pol             = band.pol,
                band_number     = band_number,
                )
        #end def edge_data

        vbm        = None
        cbm        = None
        direct_gap = None
        for band_channel in (bands.up,bands.down):
            for b in band_channel.values():
                if len(b.occs)!=len(b.eigs) or len(b.occs)==0:
                    continue
                occ   = b.occs > 0.5
                unocc = b.occs < 0.5
                if not occ.any() or not unocc.any():
                    continue
                e_val  = np.max(b.eigs[occ])
                e_cond = np.min(b.eigs[unocc])
                if vbm is None or e_val>vbm.energy:
                    vbm = edge_data(b,e_val,np.max(np.where(occ)))
                if cbm is None or e_cond<cbm.energy:
                    cbm = edge_data(b,e_cond,np.min(np.where(unocc)))
                if direct_gap is None or e_cond-e_val<direct_gap.energy:
                    direct_gap = obj(
                        energy          = e_cond-e_val,
                        kpoint_rel      = b.kpoint_rel,
                        kpoint_2pi_alat = b.kpoint_2pi_alat,
                        index           = b.index,
                        pol             = [vbm.pol,cbm.pol],
                        )
        if vbm is None:
            return
        if vbm.energy+0.025>=cbm.energy:
            electronic_structure = 'metallic' if vbm.band_number==cbm.band_number else 'semi-metal'
        else:
            electronic_structure = 'insulating'
            if (vbm.kpoint_rel is not None
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
    #end def analyze_band_edges


    def analyze_structures(self,lines):
        """Parse structures reported during relaxation or dynamics."""
        structures = obj()
        conf       = None
        i          = 0
        while i<len(lines):
            line = lines[i]
            if line.strip().startswith('CELL_PARAMETERS'):
                axes = []
                if i+3<len(lines):
                    for axis_line in lines[i+1:i+4]:
                        axis_line = axis_line.replace('D','E').replace('d','e')
                        row       = np.array(re.findall(number_pattern,axis_line),dtype=float)
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
                    positions.append(coordinates)
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
            self.info.data_status.relax_structures = True
            self.results_out.relax_structures      = structures
        return len(structures)
    #end def analyze_structures


    def analyze_pressure_volume(self,lines):
        """Parse the final reported pressure and unit-cell volume."""
        pressure = None
        volume   = None
        for line in lines:
            if 'unit-cell volume' in line:
                value = match_float(r'unit-cell volume\s*=\s*('+number_pattern+r')',line)
                if value is not None:
                    volume = value
            if 'total' in line and 'stress' in line and 'P=' in line:
                value = match_float(r'P=\s*('+number_pattern+r')',line)
                if value is not None:
                    pressure = value
        for name,value in (('pressure',pressure),('volume',volume)):
            if value is not None:
                self.info.data_status[name] = True
                self.results_out[name] = value
        return dotdict(pressure=pressure is not None,volume=volume is not None)
    #end def analyze_pressure_volume


    def analyze_stress(self,lines):
        """Parse the sequence of reported stress tensors."""
        stress = []
        for i,line in enumerate(lines):
            if 'total   stress' in line:
                rows = []
                if i+3<len(lines):
                    for stress_line in lines[i+1:i+4]:
                        stress_line = stress_line.replace('D','E').replace('d','e')
                        row         = np.array(re.findall(number_pattern,stress_line),dtype=float)
                        if len(row)<6:
                            rows = []
                            break
                        rows.append(list(row[:6]))
                if len(rows)==3:
                    stress.extend(rows)
        if len(stress)>0:
            self.info.data_status.stress = True
            self.results_out.stress      = stress
        return len(stress)//3
    #end def analyze_stress


    def analyze_forces(self,lines):
        """Parse atomic and total forces reported during the run."""
        forces     = []
        tot_forces = []
        nat        = None
        if (
            self.input is not None
            and 'system' in self.input
            and 'nat' in self.input.system
            ):
            nat = self.input.system.nat
        for i,line in enumerate(lines):
            if 'Forces acting on atoms' in line:
                aforces = []
                j       = i+1
                while j<len(lines):
                    match = re.search(
                        r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*'
                        r'('+number_pattern+r')\s+('+number_pattern+r')\s+('+number_pattern+r')',
                        lines[j],
                        )
                    if match is not None:
                        values = [
                            float(value.replace('D','E').replace('d','e'))
                            for value in match.groups()
                            ]
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
            self.info.data_status.forces   = True
            forces                         = np.array(forces,dtype=float)
            self.results_out.forces        = forces
            self.results_out.max_forces    = np.linalg.norm(forces,axis=2).max(axis=1)
        if len(tot_forces)>0:
            self.info.data_status.total_forces = True
            self.results_out.tot_forces        = np.array(tot_forces,dtype=float)
        return dotdict(forces=len(forces),total_forces=len(tot_forces))
    #end def analyze_forces


    def analyze_timing(self,lines):
        """Parse total CPU and wall-clock times."""
        def pwscf_time(text):
            """Convert a PWscf timing string to hours."""
            scales = dict(h=1.,m=60.,s=3600.)
            return sum(
                float(value.replace('D','E').replace('d','e'))/scales[unit]
                for value,unit in re.findall(r'('+number_pattern+r')\s*([hms])',text)
                )
        #end def pwscf_time

        for line in lines:
            if 'PWSCF        :' in line:
                match = re.search(r'PWSCF\s*:\s*(.*?)\s+CPU\s+(.*?)\s+WALL',line)
                if match is not None:
                    self.info.data_status.timing = True
                    self.results_out.cputime     = pwscf_time(match.group(1))
                    self.results_out.walltime    = pwscf_time(match.group(2))
                    return True
        return False
    #end def analyze_timing


    def analyze_kpoints(self,lines):
        """Parse Cartesian k-points, crystal k-points, and weights."""
        kpoints_cart,kpoints_unit,kweights = read_kpoint_tables(lines)
        if kpoints_cart is not None:
            self.info.data_status.kpoints   = True
            self.results_out.kpoints_cart   = kpoints_cart
            self.results_out.kpoints_unit   = kpoints_unit
            self.results_out.kweights       = kweights
            return len(kpoints_cart)
        return 0
    #end def analyze_kpoints


    def analyze_pw2casino(self):
        """Parse the kinetic energy from accompanying PW2CASINO output."""
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
                    self.results_out.K              = float(
                        tokens[5].replace('D','E').replace('d','e')
                        )
                    return True
        return False
    #end def analyze_pw2casino


    def analyze_xml(self):
        """Locate and parse schema or legacy PWscf XML output."""
        import xml.etree.ElementTree as ET

        self.results_xml = None
        if 'xml_status_failed' in self.info:
            del self.info.xml_status_failed
        if self.input is None or 'control' not in self.input:
            if self.info.warn:
                self.warn(
                    'xml data is not available\n'
                    'reason: input control section is not available'
                    )
            return 0
        cont = self.input.control
        if 'outdir' not in cont or 'prefix' not in cont:
            if self.info.warn:
                self.warn('xml data is not available\nreason: input outdir/prefix is not available')
            return 0
        savedir     = os.path.join(self.path,cont.outdir,f'{cont.prefix}.save')
        schema_file = os.path.join(savedir,'data-file-schema.xml')
        legacy_file = os.path.join(savedir,'data-file.xml')
        schema      = os.path.exists(schema_file)
        if not schema:
            legacy_dir = savedir
            if not os.path.exists(legacy_file):
                legacy_dir  = os.path.join(self.path,cont.outdir)
                legacy_file = os.path.join(legacy_dir,f'{cont.prefix}.xml')
            if not os.path.exists(legacy_file):
                if self.info.warn:
                    self.warn(
                        'xml data is not available\n'
                        f'file not found: {legacy_file}'
                        )
                return 0

        self.results_xml = obj(data=None,kpoints=None,failed=False)
        if schema:
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
        """Record an XML parsing failure and optionally issue a warning."""
        self.results_xml.failed = True
        if self.info.warn:
            self.warn(message)
    #end def xml_failed


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
                        'encountered an exception during xml read, '
                        'this data will not be available\n'
                        f'exception encountered: {e}'
                        )
                    continue
                eunits      = object_path(edata,'root','units_for_energies','units')
                eigenvalues = object_path(edata,'root','eigenvalues')
                occupations = object_path(edata,'root','occupations')
                if eunits is None or eigenvalues is None or occupations is None:
                    self.results_xml.failed = True
                    continue
                units = dict(ha='Ha',ry='Ry',ev='eV').get(eunits.lower()[:2],'Ha')
                spin  = obj(units=units,eigenvalues=eigenvalues,occupations=occupations)
                if si==1:
                    kp.up = spin
                elif si==2:
                    kp.down = spin
        self.results_xml.update(data=data,kpoints=kpoints)
        return len(kpoints)
    #end def analyze_legacy_xml


    def analyze_schema_xml(self,root):
        """Extract k-point and orbital data from schema-based PWscf XML."""
        def xml_value(text):
            """Convert XML text to its natural scalar or array value."""
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
            """Convert an XML element tree to an object hierarchy."""
            node = obj()
            for name,value in element.attrib.items():
                node[name.lower()] = xml_value(value)
            groups = dotdict()
            for child in element:
                name = child.tag.rsplit('}',1)[-1].lower()
                groups.setdefault(name,[]).append(child)
            for name,group in groups.items():
                node[name] = (
                    xml_element(group[0])
                    if len(group)==1
                    else obj({i+1:xml_element(child) for i,child in enumerate(group)})
                    )
            text = (element.text or '').strip()
            if len(text)>0:
                value = xml_value(text)
                if len(node)==0:
                    return value
                node.value = value
            return node
        #end def xml_element

        def xml_child(elem,name):
            """Return the named direct XML child when present."""
            if elem is None:
                return None
            name = name.lower()
            return next(
                (child for child in elem if child.tag.rsplit('}',1)[-1].lower()==name),
                None,
                )
        #end def xml_child

        data           = obj(root=xml_element(root))
        output         = xml_child(root,'output')
        band_structure = xml_child(output,'band_structure')
        kpoints        = obj()
        self.results_xml.update(data=data,kpoints=kpoints)
        if output is None or band_structure is None:
            self.xml_failed(
                'xml data is incomplete, some data will not be available\n'
                'reason: output band_structure is not available'
                )
            return 0
        lsda_element = xml_child(band_structure,'lsda')
        lsda         = xml_value(lsda_element.text) if lsda_element is not None else False
        records      = [child for child in band_structure
                        if child.tag.rsplit('}',1)[-1].lower()=='ks_energies']
        if len(records)==0:
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
            text        = (k_element.text or '').replace('D','E').replace('d','e')
            coordinates = np.array(re.findall(number_pattern,text),dtype=float)
            text        = (e_element.text or '').replace('D','E').replace('d','e')
            eigenvalues = np.array(re.findall(number_pattern,text),dtype=float)
            text        = (o_element.text or '').replace('D','E').replace('d','e')
            occupations = np.array(re.findall(number_pattern,text),dtype=float)
            if (len(coordinates)<3
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
            key    = tuple(np.round(coordinates,12))
            if not lsda or key not in coordinate_map:
                ki = len(kpoints)+1
                coordinate_map[key] = ki
                kpoints[ki] = obj(kpoint=coordinates,weight=weight)
            else:
                ki = coordinate_map[key]
            kp   = kpoints[ki]
            spin = obj(
                units       = 'Ha',
                eigenvalues = eigenvalues,
                occupations = occupations,
                )
            if not lsda or 'up' not in kp:
                kp.up = spin
            else:
                kp.down = spin
        if lsda and any('down' not in kp for kp in kpoints.values()):
            self.results_xml.failed = True
        if len(kpoints)==0:
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
        """Write or return spin-resolved electron counts from parsed XML."""
        results_xml = self.results_xml
        if not return_flag:
            if results_xml is None:
                msg = 'xml data has not been processed\ncannot write electron counts'
                raise RuntimeError(msg)
            elif results_xml.failed:
                msg = 'xml data processing failed\ncannot write electron counts'
                raise FileFormatError(msg)
        elif results_xml is None or results_xml.failed:
            return False
        kpoints      = results_xml.kpoints
        first_kpoint = next(iter(kpoints.values()))
        spins        = dotdict(
            up   = 'up',
            down = 'down' if 'down' in first_kpoint else 'up',
            )
        tot          = dotdict({
            spin:sum(kp.weight*kp[label].occupations.sum() for kp in kpoints.values())
            for spin,label in spins.items()
            })
        text        = 'total electron counts\n'
        total_count = tot.up+tot.down
        spin_count  = tot.up-tot.down
        text += f'  {total_count: 3.2f}  {spin_count: 3.2f}  {tot.up: 3.2f}  {tot.down: 3.2f}\n'
        text += '\nkpoint electron counts\n'
        weights = np.array([kp.weight for kp in kpoints.values()],dtype=float)
        mult    = (weights/weights.min()).sum()
        for ik in sorted(kpoints.keys()):
            kp  = kpoints[ik]
            kpt = dotdict({
                spin:kp.weight*kp[label].occupations.sum()*mult
                for spin,label in spins.items()
                })
            total_count = kpt.up+kpt.down
            spin_count  = kpt.up-kpt.down
            text += (
                f'  {ik:>3}  {kp.weight: 8.6f}    {total_count: 3.2f}  '
                f'{spin_count: 3.2f}  {kpt.up: 3.2f}  {kpt.down: 3.2f}\n'
                )
        if filepath is not None:
            with open(filepath,'w') as fobj:
                fobj.write(text)
        return True if return_flag else text
    #end def write_electron_counts


    def md_statistics(self,equil=None,autocorr=None):
        """Calculate summary statistics for molecular-dynamics histories."""
        mds = obj()
        for q,v in self.results_out.md_data.items():
            if equil is not None:
                v = v[equil:]
            if autocorr is None:
                mean,_,error,_ = simstats(v)
            else:
                nv = len(v)
                nb = int(np.floor(float(nv)/autocorr))
                v  = v[nv-nb*autocorr:]
                npe.reshape_inplace(v, (nb, autocorr))
                mean,error = simplestats(v.mean(axis=1))
            mds[q] = mean,error
        return mds
    #end def md_statistics


    def md_plots(self,*,show=True):
        """Plot energy, temperature, and pressure histories from dynamics."""

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
        """Write relaxed or dynamic structures as a tiled XYZ movie."""
        if ('relax_structures' not in self.results_out
            or self.results_out.relax_structures is None
            ):
            return
        structures = self.results_out.relax_structures
        filepath = self.abspath if filepath is None else filepath
        filepath = os.path.join(filepath,filename)

        movie         = ''
        alat_angstrom = convert(self.input.system['celldm(1)'],'B','A')
        cell          = self.input.cell_parameters.vectors
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
        show_vbm_cbm  = True,
        k_labels      = None,
        ):
        """Plot the analyzed band structure along its reciprocal-space path."""
        import matplotlib.pyplot as plt
        bands = self.results_out.bands
        if bands is None:
            return
        params = {
            'legend.fontsize'      : 14,
            'figure.facecolor'     : 'white',
            'figure.subplot.hspace': 0.,
            'axes.labelsize'       : 16,
            'xtick.labelsize'      : 14,
            'ytick.labelsize'      : 14,
            }
        plt.rcParams.update(params)
        filename = 'band_structure.pdf' if filename is None else filename
        filepath = os.path.join(self.abspath if filepath is None else filepath,filename)

        plt.figure()
        ax     = plt.gca()
        nbands = self.input.system.nbnd

        if k_labels is None:
            kpath  = get_kpath(structure=self.input_structure, check_standard=False)
            x      = kpath['explicit_path_linearcoords']
            labels = list(kpath['explicit_kpoints_labels'])
        else:
            labels = list(k_labels)
            # Calculate linear coordinates from self.results_out.kpoints_cart
            x          = []
            prev_label = ''
            ref_kpt    = self.results_out.kpoints_cart[0]
            lincoord   = 0.0
            for kpt_idx,kpt in enumerate(self.results_out.kpoints_cart):
                curr_label = labels[kpt_idx]
                if (curr_label != '' and prev_label == '') or curr_label == '':
                    lincoord += np.linalg.norm(kpt-ref_kpt)
                ref_kpt = kpt
                x.append(lincoord)
                prev_label = curr_label
        for band_channel,color in ((bands.up,'k'),(bands.down,'r')):
            if len(band_channel)==0:
                continue
            for nb in range(nbands):
                eigs = np.array([band.eigs[nb] for band in band_channel],dtype=float)
                y = eigs - bands.vbm.energy
                plt.plot(x,y,color)
        for ln,li in enumerate(labels):
            if li != '':
                plt.axvline(x[ln],ymin=-100,ymax=100,linewidth=3,color='k')
                labels[ln] = r'$\Gamma$' if li=='GAMMA' else f'${li}$'
                if ln>0 and labels[ln-1]!='':
                    labels[ln] = f'{labels[ln-1]}|{labels[ln]}'
                    labels[ln-1] = ''

        plt.xlim([np.min(x),np.max(x)])
        plt.ylim((-5,+5) if max_min_e is None else max_min_e)
        plt.ylabel('Energy (eV)')
        plt.xticks(x,labels)
        ax.tick_params(axis='x',which='both',length=0)
        ax.tick_params(axis='x',which='both',pad=10)
        if show_vbm_cbm:
            vbm = bands.vbm
            cbm = bands.cbm
            for kn,ki in enumerate(bands.up):
                if (vbm.kpoint_rel==ki.kpoint_rel).all():
                    plt.scatter(x[kn],0,c='green',s=100)
                if (cbm.kpoint_rel==ki.kpoint_rel).all():
                    plt.scatter(x[kn],cbm.energy-vbm.energy,c='r',s=100)
        if save:
            plt.savefig(filepath,format='pdf',bbox_inches='tight')
        if show:
            plt.show()
    #end def plot_bandstructure

#end class PwscfAnalyzer
        
