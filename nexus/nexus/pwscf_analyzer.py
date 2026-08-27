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
import xml.etree.ElementTree as ET
import numpy as np
from .developer import obj
from .unit_converter import convert
from .periodic_table import Elements
from .numerics import simstats, simplestats
from .simulation import SimulationAnalyzer, Simulation
from .structure import Structure, get_kpath
from .pwscf_input import PwscfInput
from .pwscf_data_reader import read_qexml
from .utilities import path_string
from . import numpy_extensions as npe


elements = set([e.symbol for e in Elements])


def is_number(s):
    number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
    return re.fullmatch(number_pattern,s.strip()) is not None
#end def is_number


def pwscf_time(tsin):
    ts = tsin
    h,m,s='','',''
    if ts!='' and ts.find('h')!=-1:
        sp = ts.split('h')
        h = sp[0]
        ts = sp[1]
    #end if
    if ts!='' and ts.find('m')!=-1:
        sp = ts.split('m')
        m = sp[0]
        ts = sp[1]
    #end if
    if ts!='' and ts.find('s')!=-1:
        sp = ts.split('s')
        s = sp[0]
        ts = sp[1]
    #end if

    times = [h,m,s]
    time = 0.
    for n in range(3):
        t = times[n]
        if is_number(t):
            t=float(t)
        else:
            t=0
        #end if
        time += t/(60.)**n
    #end for
    
    return time
#end def pwscf_time


class PwscfAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,infile_name=None,outfile_name=None,pw2c_outfile_name=None,*,analyze=False,xml=False,warn=False,md_only=False):
        if isinstance(arg0,Simulation):
            sim = arg0
            path = sim.locdir
            infile_name = sim.infile
            outfile_name= sim.outfile
            self.input_structure = sim.system.structure
        elif arg0 is not None:
            path = path_string(arg0)
            if not os.path.exists(path):
                self.error('path to QE data does not exist\npath provided: {}'.format(path))
            #end if
            if os.path.isfile(path):
                filepath = path
                path,filename = os.path.split(filepath)
                if filename.endswith('.in'):
                    infile_name = filename
                elif filename.endswith('.out'):
                    outfile_name = filename
                else:
                    self.error('could not determine whether file is QE input or output\nfile provided: {}'.format(filepath))
                #end if
            #end if
            if outfile_name is None:
                outfile_name = infile_name.rsplit('.',1)[0]+'.out'
            #end if
        else:
            return
        #end if

        self.infile_name  = infile_name
        self.outfile_name = outfile_name

        self.path = path
        self.abspath = os.path.abspath(path)
        self.pw2c_outfile_name = pw2c_outfile_name

        self.info = obj(xml=xml,warn=warn,md_only=md_only)

        self.input = None
        if self.infile_name is not None:
            self.input = PwscfInput(os.path.join(self.path,self.infile_name))
        #end if
        if analyze:
            self.analyze()
        #end if
    #end def __init__

    
    number_pattern = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


    def analyze(self):
        outfile = os.path.join(self.path,self.outfile_name)
        try:
            with open(outfile,'r') as fobj:
                lines = fobj.read().splitlines()
        except OSError as e:
            if self.info.warn:
                self.warn('file read failed\nexception encountered: '+str(e))
            #end if
            return
        #end try

        self.analyze_md(lines)
        if self.info.md_only:
            return
        #end if
        self.analyze_fermi_energies(lines)
        self.analyze_energies(lines)
        self.analyze_bands(lines)
        self.analyze_structures(lines)
        self.analyze_pressure_volume(lines)
        self.analyze_stress(lines)
        self.analyze_forces(lines)
        self.analyze_timing(lines)
        self.analyze_kpoints(lines)
        self.analyze_pw2casino()
        if self.info.xml:
            self.analyze_xml()
        #end if
    #end def analyze


    def line_numbers(self,line):
        values = re.findall(self.number_pattern,line.replace('D','E').replace('d','e'))
        return np.array(values,dtype=float)
    #end def line_numbers


    def analyze_md(self,lines):
        calculation = ''
        if self.input is not None and 'control' in self.input:
            calculation = self.input.control.calculation if 'calculation' in self.input.control else 'scf'
            calculation = calculation.lower()
        #end if
        if calculation not in ('md','vc-md'):
            return
        #end if

        energies = []
        pressures = []
        times = []
        kinetic_energies = []
        temperatures = []
        for l in lines:
            if l.lstrip().startswith('!') and 'total energy' in l:
                energies.append(float(l.split('=',1)[1].split()[0]))
            elif 'total   stress' in l and 'P=' in l:
                pressures.append(float(l.rsplit('P=',1)[1].split()[0]))
            #end if
            if calculation=='md':
                if 'time      =' in l:
                    times.append(float(l.split('=',1)[1].split()[0]))
                elif 'kinetic energy' in l and '=' in l:
                    kinetic_energies.append(float(l.split('=',1)[1].split()[0]))
                elif l.strip().startswith('temperature') and '=' in l:
                    temperatures.append(float(l.split('=',1)[1].split()[0]))
                #end if
            else:
                if 'Entering Dynamics;' in l and 'time' in l:
                    times.append(float(l.rsplit('=',1)[1].split()[0]))
                elif l.strip().startswith('Ekin'):
                    match = re.search(r'Ekin\s*=\s*('+self.number_pattern+r')\s+Ry\s+T\s*=\s*('+self.number_pattern+r')',l)
                    kinetic_energies.append(float(match.group(1)))
                    temperatures.append(float(match.group(2)))
                #end if
            #end if
        #end for
        nsteps = len(times)
        if nsteps==0:
            return
        #end if
        quantities = (energies,pressures,times,kinetic_energies,temperatures)
        if not all(len(q)==nsteps for q in quantities):
            self.error('inconsistent MD record lengths\ncalculation: {0}\nrecord lengths: {1}'.format(
                calculation,[len(q) for q in quantities]))
        #end if
        md = obj(
            total_energy   = np.array(energies,dtype=float),
            pressure       = np.array(pressures,dtype=float),
            time           = np.array(times,dtype=float),
            kinetic_energy = np.array(kinetic_energies,dtype=float),
            temperature    = np.array(temperatures,dtype=float),
            )
        md.potential_energy = md.total_energy - md.kinetic_energy
        self.md_data = md
        self.md_stats = self.md_statistics()
    #end def analyze_md


    def analyze_fermi_energies(self,lines):
        fermi_energies = []
        for l in lines:
            if 'Fermi energ' in l:
                tokens = l.lower().replace('ev','').split()
                values = []
                for token in tokens[::-1]:
                    if is_number(token):
                        values.append(float(token))
                    elif len(values)>0:
                        break
                    #end if
                #end for
                fermi_energies.extend(values[::-1])
            #end if
        #end for
        self.Ef = fermi_energies[-1] if len(fermi_energies)>0 else 0.0
        self.fermi_energies = np.array(fermi_energies)
    #end def analyze_fermi_energies


    def analyze_energies(self,lines):
        energies = []
        for l in lines:
            if l.lstrip().startswith('!') and 'total energy' in l:
                energies.append(float(l.split('=',1)[1].split()[0]))
            #end if
        #end for
        self.E = energies[-1] if len(energies)>0 else 0.0
        self.energies = np.array(energies)
    #end def analyze_energies


    def read_kpoint_tables(self,lines):
        kpoints_cart = None
        kpoints_unit = None
        kweights = None
        for i,l in enumerate(lines):
            if 'number of k points=' not in l:
                continue
            #end if
            match = re.search(r'number of k points=\s*(\d+)',l)
            nkpoints = int(match.group(1))
            if i+1>=len(lines) or 'cart. coord.' not in lines[i+1]:
                continue
            #end if
            cart = []
            weights = []
            for kl in lines[i+2:i+2+nkpoints]:
                match = re.search(r'=\s*\((.*?)\),\s*wk\s*=\s*('+self.number_pattern+r')',kl)
                cart.append(self.line_numbers(match.group(1)))
                weights.append(float(match.group(2)))
            #end for
            j = i+2+nkpoints
            while j<len(lines) and 'cryst. coord.' not in lines[j]:
                j+=1
            #end while
            if j>=len(lines):
                continue
            #end if
            unit = []
            for kl in lines[j+1:j+1+nkpoints]:
                match = re.search(r'=\s*\((.*?)\),\s*wk\s*=',kl)
                unit.append(self.line_numbers(match.group(1)))
            #end for
            kpoints_cart = np.array(cart,dtype=float)
            kpoints_unit = np.array(unit,dtype=float)
            kweights = np.array(weights,dtype=float)
            break
        #end for
        return kpoints_cart,kpoints_unit,kweights
    #end def read_kpoint_tables


    def analyze_bands(self,lines):
        table_cart,table_unit,_ = self.read_kpoint_tables(lines)
        nspin = self.input.system.nspin if self.input is not None and 'nspin' in self.input.system else 1
        polarized = nspin>1
        bands = obj(up=obj(),down=obj())
        band_channel = bands.up
        up_spin = True
        nfound = 0
        for i,l in enumerate(lines):
            if 'End of self-consistent calculation' in l and nfound>0:
                bands = obj(up=obj(),down=obj())
                band_channel = bands.up
                up_spin = True
                nfound = 0
                continue
            elif '- SPIN UP -' in l:
                up_spin = True
                band_channel = bands.up
                continue
            elif '- SPIN DOWN -' in l:
                up_spin = False
                band_channel = bands.down
                continue
            elif 'bands (ev)' not in l:
                continue
            #end if

            eigs = []
            j = i+1
            while j<len(lines):
                ls = lines[j].strip()
                if len(ls)==0:
                    if len(eigs)>0:
                        break
                    #end if
                elif 'occupation numbers' in ls or 'bands (ev)' in ls:
                    break
                else:
                    values = self.line_numbers(ls)
                    remainder = re.sub(self.number_pattern,'',ls)
                    if len(values)==0 or len(remainder.strip())>0:
                        break
                    #end if
                    eigs.extend(values)
                #end if
                j+=1
            #end while
            while j<len(lines) and len(lines[j].strip())==0:
                j+=1
            #end while
            occs = []
            if j<len(lines) and 'occupation numbers' in lines[j]:
                j+=1
                while j<len(lines):
                    ls = lines[j].strip()
                    if len(ls)==0:
                        if len(occs)>0:
                            break
                        #end if
                    else:
                        values = self.line_numbers(ls)
                        remainder = re.sub(self.number_pattern,'',ls)
                        if len(values)==0 or len(remainder.strip())>0:
                            break
                        #end if
                        occs.extend(values)
                    #end if
                    j+=1
                #end while
            #end if

            match = re.search(r'k\s*=\s*(.*?)\s*\(',l)
            kpoint_cart = self.line_numbers(match.group(1))
            index = len(band_channel)
            kpoint_rel = table_unit[index] if table_unit is not None and index<len(table_unit) else kpoint_cart.copy()
            if table_cart is not None and index<len(table_cart):
                kpoint_cart = table_cart[index]
            #end if
            pol = 'none'
            if polarized:
                pol = 'up' if up_spin else 'down'
            #end if
            band_channel[index] = obj(
                index           = index,
                kpoint_2pi_alat = kpoint_cart,
                kpoint_rel      = kpoint_rel,
                eigs            = np.array(eigs,dtype=float),
                occs            = np.array(occs,dtype=float),
                pol             = pol,
                )
            nfound+=1
        #end for
        if nfound==0:
            return
        #end if
        self.analyze_band_edges(bands)
        self.bands = bands
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
                #end if
                occ = b.occs > 0.5
                unocc = b.occs < 0.5
                if not occ.any() or not unocc.any():
                    continue
                #end if
                found_occupations = True
                e_val = np.max(b.eigs[occ])
                e_cond = np.min(b.eigs[unocc])
                if e_val > vbm.energy:
                    vbm = obj(energy=e_val,kpoint_rel=b.kpoint_rel,kpoint_2pi_alat=b.kpoint_2pi_alat,
                              index=b.index,pol=b.pol,band_number=np.max(np.where(occ)))
                #end if
                if e_cond < cbm.energy:
                    cbm = obj(energy=e_cond,kpoint_rel=b.kpoint_rel,kpoint_2pi_alat=b.kpoint_2pi_alat,
                              index=b.index,pol=b.pol,band_number=np.min(np.where(unocc)))
                #end if
                if e_cond-e_val < direct_gap.energy:
                    direct_gap = obj(energy=e_cond-e_val,kpoint_rel=b.kpoint_rel,kpoint_2pi_alat=b.kpoint_2pi_alat,
                                     index=b.index,pol=[vbm.pol,cbm.pol])
                #end if
            #end for
        #end for
        if not found_occupations:
            return
        #end if
        if vbm.energy+0.025>=cbm.energy:
            electronic_structure = 'metallic' if vbm.band_number==cbm.band_number else 'semi-metal'
        else:
            electronic_structure = 'insulating'
            if not np.equal(vbm.kpoint_rel,cbm.kpoint_rel).all():
                bands.indirect_gap = obj(energy=round(cbm.energy-vbm.energy,3),kpoints=obj(vbm=vbm,cbm=cbm))
            #end if
        #end if
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
            l = lines[i]
            if l.strip().startswith('CELL_PARAMETERS'):
                conf = obj()
                axes = np.array([lines[i+d+1].split() for d in range(3)],dtype=float)
                match = re.search(r'alat\s*=\s*('+self.number_pattern+r')',l)
                if match is not None:
                    axes *= float(match.group(1))
                #end if
                conf.axes = axes
                i+=3
            elif 'ATOMIC_POSITIONS' in l:
                if conf is None:
                    conf = obj()
                #end if
                atoms = []
                positions = []
                i+=1
                while i<len(lines):
                    tokens = lines[i].split()
                    if len(tokens) not in (4,7) or len(tokens)==0 or tokens[0].lower()=='end':
                        break
                    #end if
                    atoms.append(tokens[0])
                    positions.append(np.array(tokens[1:4],dtype=float))
                    i+=1
                #end while
                conf.atoms = atoms
                conf.positions = np.array(positions)
                if 'crystal' in l.lower() and 'axes' in conf:
                    conf.positions = np.dot(conf.positions,conf.axes)
                #end if
                structures[len(structures)] = conf
                conf = None
                continue
            #end if
            i+=1
        #end while
        if len(structures)>0:
            self.structures = structures
        #end if
    #end def analyze_structures


    def analyze_pressure_volume(self,lines):
        press = 0.
        vol = 0.
        for l in lines:
            if 'unit-cell volume' in l:
                match = re.search(r'unit-cell volume\s*=\s*('+self.number_pattern+r')',l)
                vol = float(match.group(1))
            #end if
            if 'total' in l and 'stress' in l and 'P=' in l:
                press = float(l.rsplit('P=',1)[1].split()[0])
            #end if
        #end for
        self.pressure = press
        self.volume = vol
    #end def analyze_pressure_volume


    def analyze_stress(self,lines):
        stress = []
        for i,l in enumerate(lines):
            if 'total   stress' in l:
                stress.extend([list(np.array(sl.split(),dtype=float)) for sl in lines[i+1:i+4]])
            #end if
        #end for
        self.stress = stress
    #end def analyze_stress


    def analyze_forces(self,lines):
        forces = []
        tot_forces = []
        for i,l in enumerate(lines):
            if 'Forces acting on atoms' in l:
                aforces = []
                j = i+1
                while j<len(lines):
                    match = re.search(r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*(.*)',lines[j])
                    if match is not None:
                        aforces.append(self.line_numbers(match.group(1)))
                    elif len(aforces)>0:
                        break
                    #end if
                    j+=1
                #end while
                forces.append(aforces)
            #end if
            if 'Total force' in l:
                match = re.search(r'Total force\s*=\s*('+self.number_pattern+r')',l)
                tot_forces.append(float(match.group(1)))
            #end if
        #end for
        if len(forces)>0:
            self.forces = np.array(forces,dtype=float)
            self.tot_forces = np.array(tot_forces)
            self.max_forces = np.array([(np.sqrt((f**2).sum(1))).max() for f in self.forces])
        #end if
    #end def analyze_forces


    def analyze_timing(self,lines):
        tc = 0.
        tw = 0.
        for l in lines:
            if 'PWSCF        :' in l:
                times = l.split(':',1)[1].split('CPU')
                tc = pwscf_time(times[0])
                tw = pwscf_time(times[1].replace('WALL',''))
                break
            #end if
        #end for
        self.cputime = tc
        self.walltime = tw
    #end def analyze_timing


    def analyze_kpoints(self,lines):
        kpoints_cart,kpoints_unit,kweights = self.read_kpoint_tables(lines)
        if kpoints_cart is not None:
            self.kpoints_cart = kpoints_cart
            self.kpoints_unit = kpoints_unit
            self.kweights = kweights
        #end if
    #end def analyze_kpoints


    def analyze_pw2casino(self):
        if self.pw2c_outfile_name is None:
            return
        #end if
        filepath = os.path.join(self.path,self.pw2c_outfile_name)
        try:
            with open(filepath,'r') as fobj:
                lines = fobj.readlines()
        except OSError as e:
            if self.info.warn:
                self.warn('pw2casino file read failed\nexception encountered: '+str(e))
            #end if
            return
        #end try
        for l in lines:
            if 'Kinetic' in l:
                self.K = float(l.split()[5])
                break
            #end if
        #end for
    #end def analyze_pw2casino


    def analyze_xml(self):
        self.xmldata = obj(data=None,kpoints=None,failed=False)
        cont = self.input.control
        savedir = os.path.join(self.path,cont.outdir,cont.prefix+'.save')
        schema_file = os.path.join(savedir,'data-file-schema.xml')
        legacy_file = os.path.join(savedir,'data-file.xml')
        legacy_dir = savedir
        if os.path.exists(schema_file):
            try:
                root = ET.parse(schema_file).getroot()
            except (OSError,ET.ParseError) as e:
                self.xml_read_failed(e)
                return
            #end try
            self.analyze_schema_xml(root)
        else:
            if not os.path.exists(legacy_file):
                legacy_dir = os.path.join(self.path,cont.outdir)
                legacy_file = os.path.join(legacy_dir,cont.prefix+'.xml')
            #end if
            try:
                data = read_qexml(legacy_file)
            except OSError as e:
                self.xml_read_failed(e)
                return
            #end try
            self.analyze_legacy_xml(data,legacy_dir)
        #end if
    #end def analyze_xml


    def xml_read_failed(self,exception):
        self.xmldata.failed = True
        if self.info.warn:
            self.warn('encountered an exception during xml read, this data will not be available\nexception encountered: '+str(exception))
        #end if
    #end def xml_read_failed


    def analyze_legacy_xml(self,data,datadir):
        kpdata = data.root.eigenvalues.k_point
        kpoints = obj()
        for ki,kpd in kpdata.items():
            kp = obj(kpoint=kpd.k_point_coords,weight=kpd.weight)
            kpoints[ki] = kp
            for si,dfile in kpd.datafile.items():
                efilepath = os.path.join(datadir,dfile.iotk_link)
                if not os.path.exists(efilepath):
                    self.xmldata.failed = True
                    continue
                #end if
                try:
                    edata = read_qexml(efilepath)
                except OSError as e:
                    self.xml_read_failed(e)
                    continue
                #end try
                eunits = edata.root.units_for_energies.units.lower()
                if eunits.startswith('ha'):
                    units = 'Ha'
                elif eunits.startswith('ry'):
                    units = 'Ry'
                elif eunits.startswith('ev'):
                    units = 'eV'
                else:
                    units = 'Ha'
                #end if
                spin = obj(units=units,eigenvalues=edata.root.eigenvalues,occupations=edata.root.occupations)
                if si==1:
                    kp.up = spin
                elif si==2:
                    kp.down = spin
                #end if
            #end for
        #end for
        self.xmldata.update(data=data,kpoints=kpoints)
    #end def analyze_legacy_xml


    def xml_local_name(self,element):
        return element.tag.rsplit('}',1)[-1].lower()
    #end def xml_local_name


    def xml_value(self,text):
        text = text.strip()
        if len(text)==0:
            return ''
        elif text.lower() in ('true','false'):
            return text.lower()=='true'
        #end if
        tokens = text.replace('D','E').replace('d','e').split()
        if all(is_number(token) for token in tokens):
            values = np.array(tokens,dtype=float)
            if len(values)==1:
                value = values[0]
                if value.is_integer() and all(c not in text.lower() for c in ('.','e')):
                    return int(value)
                #end if
                return value
            #end if
            return values
        #end if
        return text
    #end def xml_value


    def xml_element(self,element):
        node = obj()
        for name,value in element.attrib.items():
            node[name.lower()] = self.xml_value(value)
        #end for
        children = list(element)
        groups = obj()
        for child in children:
            name = self.xml_local_name(child)
            if name not in groups:
                groups[name] = []
            #end if
            groups[name].append(child)
        #end for
        for name,group in groups.items():
            if len(group)==1:
                node[name] = self.xml_element(group[0])
            else:
                node[name] = obj({i+1:self.xml_element(child) for i,child in enumerate(group)})
            #end if
        #end for
        text = (element.text or '').strip()
        if len(text)>0:
            value = self.xml_value(text)
            if len(node)==0:
                return value
            #end if
            node.value = value
        #end if
        return node
    #end def xml_element


    def xml_child(self,element,name):
        name = name.lower()
        for child in element:
            if self.xml_local_name(child)==name:
                return child
            #end if
        #end for
        return None
    #end def xml_child


    def analyze_schema_xml(self,root):
        data = obj(root=self.xml_element(root))
        output = self.xml_child(root,'output')
        band_structure = self.xml_child(output,'band_structure')
        lsda_element = self.xml_child(band_structure,'lsda')
        lsda = self.xml_value(lsda_element.text) if lsda_element is not None else False
        records = [child for child in band_structure if self.xml_local_name(child)=='ks_energies']
        kpoints = obj()
        coordinate_map = dict()
        for record in records:
            k_element = self.xml_child(record,'k_point')
            e_element = self.xml_child(record,'eigenvalues')
            o_element = self.xml_child(record,'occupations')
            coordinates = np.array(self.xml_value(k_element.text),ndmin=1,dtype=float)
            weight = float(k_element.attrib.get('weight',0.0))
            key = tuple(np.round(coordinates,12))
            if not lsda or key not in coordinate_map:
                ki = len(kpoints)+1
                coordinate_map[key] = ki
                kpoints[ki] = obj(kpoint=coordinates,weight=weight)
            else:
                ki = coordinate_map[key]
            #end if
            kp = kpoints[ki]
            spin = obj(
                units       = 'Ha',
                eigenvalues = np.array(self.xml_value(e_element.text),ndmin=1,dtype=float),
                occupations = np.array(self.xml_value(o_element.text),ndmin=1,dtype=float),
                )
            if not lsda or 'up' not in kp:
                kp.up = spin
            else:
                kp.down = spin
            #end if
        #end for
        self.xmldata.update(data=data,kpoints=kpoints)
    #end def analyze_schema_xml


    def write_electron_counts(self,filepath=None,*,return_flag=False):
        if not return_flag:
            if not self.info.xml:
                self.error('xml data has not been processed\ncannot write electron counts')
            elif self.xmldata.failed:
                self.error('xml data processing failed\ncannot write electron counts')
            #end if
        elif not self.info.xml or self.xmldata.failed:
            return False
        #end if
        kpoints = self.xmldata.kpoints
        if 'down' in kpoints[1]:
            spins = obj(up='up',down='down')
        else:
            spins = obj(up='up',down='up')
        #end if
        tot = obj(up=0,down=0)
        for kp in kpoints.values():
            w = kp.weight
            for s,sl in spins.items():
                tot[s] += w*kp[sl].occupations.sum()
            #end for
        #end for
        text = 'total electron counts\n'
        text += '  {0: 3.2f}  {1: 3.2f}  {2: 3.2f}  {3: 3.2f}\n'.format(tot.up+tot.down,tot.up-tot.down,tot.up,tot.down)
        text += '\nkpoint electron counts\n'
        weights = []
        for kp in kpoints.values():
            weights.append(kp.weight)
        #end for
        weights = np.array(weights,dtype=float)
        mult = (weights/weights.min()).sum()
        for ik in sorted(kpoints.keys()):
            kp = kpoints[ik]
            kpt = obj()
            for s,sl in spins.items():
                kpt[s] = w*kp[sl].occupations.sum()*mult
            #end for
            #text+='  {0:>3}  {1: 8.6f}    {2: 3.2f}  {3: 3.2f}  {4: 3.2f}    {5}\n'.format(ik,kp.weight,kpt.up+kpt.down,kpt.up,kpt.down,kp.kpoint[0])
            text+='  {0:>3}  {1: 8.6f}    {2: 3.2f}  {3: 3.2f}  {4: 3.2f}  {5: 3.2f}\n'.format(ik,kp.weight,kpt.up+kpt.down,kpt.up-kpt.down,kpt.up,kpt.down)
        #end for
        if filepath is not None:
            with open(filepath,'w') as fobj:
                fobj.write(text)
        #end if
        if not return_flag:
            return text
        else:
            return True
        #end if
    #end def write_electron_counts


    def md_statistics(self,equil=None,autocorr=None):
        mds = obj()
        for q,v in self.md_data.items():
            if equil is not None:
                v = v[equil:]
            #end if
            if autocorr is None:
                mean,var,error,kappa = simstats(v)
            else:
                nv = len(v)
                nb = int(np.floor(float(nv)/autocorr))
                nexclude = nv-nb*autocorr
                v = v[nexclude:]
                npe.reshape_inplace(v, (nb, autocorr))
                mean,error = simplestats(v.mean(axis=1))
            #end if
            mds[q] = mean,error
        #end for
        return mds
    #end def md_statistics


    def md_plots(self,*,show=True):

        md = self.md_data

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
        #end if

        return fig
    #end def md_plots


    def make_movie(self,filename,filepath=None):
        if 'structures' in self:
            if filepath is None:
                filepath = os.path.join(self.abspath,filename)
            else:
                filepath = os.path.join(filepath,filename)
            #end if
            movie = ''
            structures = self.structures

            aA = convert(self.input.system['celldm(1)'],'B','A')
            cell = self.input.cell_parameters.vectors
            for i in range(len(structures)):
                s = structures[i]
                struct = Structure(elem=s.atoms,pos=s.positions,axes=cell,scale=aA,units='A')
                struct=struct.tile(2,2,2)
                ss=struct.write_xyz()
                movie += ss
                with open(filepath,'w') as fobj:
                    fobj.write(movie)
            #end for
        #end for
    #end def make_movie


    def plot_bandstructure(self, filename=None, filepath=None, max_min_e = None, *, show=False, save=True, show_vbm_cbm=True,k_labels=None):
        import matplotlib.pyplot as plt
        params = {
            'legend.fontsize'      : 14,
            'figure.facecolor'     : 'white',
            'figure.subplot.hspace': 0.,
            'axes.labelsize'       : 16,
            'xtick.labelsize'      : 14,
            'ytick.labelsize'      : 14,
            }
        plt.rcParams.update(params)
        if 'bands' in self:
            success = True
            if filename is None:
                filename = 'band_structure.pdf'
            if filepath is None:
                filepath = os.path.join(self.abspath,filename)
            else:
                filepath = os.path.join(filepath,filename)
            #end if

            fig    = plt.figure()
            ax     = plt.gca()
            nbands = self.input.system.nbnd

            if k_labels is None:
                kpath  = get_kpath(structure=self.input_structure, check_standard=False)
                x      = kpath['explicit_path_linearcoords']
                labels = kpath['explicit_kpoints_labels']
            else:
                labels = k_labels
                # Calculate linear coordinates from self.kpoints_cart
                x = []
                prev_label = ''
                ref_kpt = self.kpoints_cart[0]
                lincoord = 0.0
                for kpt_idx,kpt in enumerate(self.kpoints_cart):
                    curr_label = labels[kpt_idx]
                    if (curr_label != '' and prev_label == '') or curr_label == '':
                        lincoord+=np.linalg.norm(kpt-ref_kpt)
                        ref_kpt = kpt
                    else:
                        ref_kpt = kpt
                        lincoord+=np.linalg.norm(kpt-ref_kpt)
                    #end if
                    x.append(lincoord)
                    prev_label = curr_label
                #end for
            #end if
            for nb in range(nbands):
                y = []
                for bi in self.bands.up:
                    y.append(bi['eigs'][nb])
                #end for
                y = np.array(y) - self.bands.vbm.energy
                plt.plot(x, y, 'k')
                if len(self.bands.down) > 0:
                    y = []
                    for bi in self.bands.down:
                        y.append(bi['eigs'][nb])
                    #end for
                    y = np.array(y) - self.bands.vbm.energy
                    plt.plot(x, y, 'r')
                #end if              
            #end for
            for ln, li in enumerate(labels):
                if li != '':
                    plt.axvline(x[ln], ymin=-100, ymax=100, linewidth=3, color='k')
                    if li == 'GAMMA':
                        labels[ln] = r'$\Gamma$'
                    elif li != '':
                        labels[ln] = '${0}$'.format(li)
                    #end if
                    if labels[ln-1] != '' and ln > 0:
                        labels[ln] = labels[ln-1]+'|'+labels[ln]
                        labels[ln-1] = ''
                    #end if
                #end if
            #end for
            
            plt.xlim([np.min(x), np.max(x)])
            if max_min_e is None:
                plt.ylim(-5, +5)
            else:
                plt.ylim(max_min_e[0],max_min_e[1])
            #end if
            plt.ylabel('Energy (eV)')
            plt.xticks(x, labels)
            ax.tick_params(axis='x', which='both', length=0)
            ax.tick_params(axis='x', which='both', pad=10)
        #end if
        if show_vbm_cbm:
            vbm = self.bands.vbm
            cbm = self.bands.cbm
            for kn, ki in enumerate(self.bands.up):
                if (vbm.kpoint_rel == ki['kpoint_rel']).all():
                    plt.scatter(x[kn], 0, c='green', s=100)
                #end if
                if (cbm.kpoint_rel == ki['kpoint_rel']).all():
                    plt.scatter(x[kn], cbm.energy-vbm.energy, c='r', s=100)
                #end if
            #end for
        #end if
        if save:
            plt.savefig(filename, format='pdf',bbox_inches='tight')
        #end if
        if show:
            plt.show()
        #end if
    #end def plot_bandstructure

#end class PwscfAnalyzer
        
