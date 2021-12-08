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
import numpy as np
from numpy import array,fromstring,sqrt,dot,max,equal,zeros,min,where
from generic import obj
from unit_converter import convert
from periodic_table import PeriodicTable
from simulation import SimulationAnalyzer,Simulation
from pwscf_input import PwscfInput
from pwscf_data_reader import read_qexml
from fileio import TextFile
from debug import *
import code
import pdb

pt = PeriodicTable()
elements = set(pt.elements.keys())




def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    #end try
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
    def __init__(self,arg0=None,infile_name=None,outfile_name=None,pw2c_outfile_name=None,analyze=False,xml=False,warn=False,md_only=False):
        if isinstance(arg0,Simulation):
            sim = arg0
            path = sim.locdir
            infile_name = sim.infile
            outfile_name= sim.outfile
            self.input_structure = sim.system.structure
        elif arg0 is not None:
            path = arg0
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

    
    def analyze(self):
        path = self.path
        infile_name = self.infile_name
        outfile_name = self.outfile_name
        pw2c_outfile_name = self.pw2c_outfile_name

        nx=0

        outfile = os.path.join(path,outfile_name)

        try:
            # perform MD analysis
            f = TextFile(outfile)
            n = 0
            md_res = []
            while f.seek('!',1)!=-1:
                E = float(f.readtokens()[-2])
                f.seek('P=',1)
                P = float(f.readtokens()[-1])
                f.seek('time      =',1)
                t = float(f.readtokens()[-2])
                f.seek('kinetic energy',1)
                K = float(f.readtokens()[-2])
                f.seek('temperature',1)
                T = float(f.readtokens()[-2])
                md_res.append((E,P,t,K,T))
                n+=1
            #end while
            md_res = array(md_res,dtype=float).T
            quantities = ('total_energy','pressure','time','kinetic_energy',
                          'temperature')
            md = obj()
            for i,q in enumerate(quantities):
                md[q] = md_res[i]
            #end for
            md.potential_energy = md.total_energy - md.kinetic_energy
            self.md_data = md
            self.md_stats = self.md_statistics()
            if self.info.md_only:
                return
            #end if
        except:
            nx+=1
            if self.info.warn:
                self.warn('MD analysis failed')
            #end if
        #end try

        try:
            lines = open(outfile,'r').read().splitlines()
        except:
            nx+=1
            if self.info.warn:
                self.warn('file read failed')
            #end if
        #end try

        try:
            fermi_energies = []
            for l in lines:
                if l.find('Fermi energ')!=-1:
                    toks = l.split()[::-1]
                    assert toks[0] == 'ev'
                    for tok in toks[1:]:
                      try:
                        ef1 = float(tok)
                        fermi_energies.append(ef1)
                      except ValueError:
                        fermi_energies = fermi_energies[::-1]
                        break
                      #end try
                    #end for
                #end if
            #end for
            if len(fermi_energies)==0:
                self.Ef = 0.0
            else:
                self.Ef = fermi_energies[-1]
            #end if
            self.fermi_energies = array(fermi_energies)
        except:
            nx+=1
            if self.info.warn:
                self.warn('fermi energy read failed')
            #end if
        #end try
        try:
            energies = []
            for l in lines:
                if l.find('!  ')!=-1:
                    energies.append( float( l.split('=')[1].split()[0] ) )
                #end if
            #end for
            if len(energies)==0:
                self.E = 0.0
            else:
                self.E = energies[-1]
            #end if
            self.energies = array(energies)
        except:
            nx+=1
            if self.info.warn:
                self.warn('total energy read failed')
            #end if
        #end try
        try:
            # get bands and occupations
            nfound = 0
            index = -1
            bands = obj()
            bands.up = obj()
            bands.down = obj()
            polarized = False
            if self.input.system.nspin > 1:
                polarized = True
            #end if
            read_kpoints  = False
            read_2pi_alat = False
            read_rel      = False
            for i in range(len(lines)):
                l = lines[i]
                if 'End of self-consistent calculation' in l:
                    # Initialize each time in case a hybrid functional was used
                    if nfound > 0:
                        nfound = 0
                        index = -1
                        bands = obj()
                        bands.up = obj()
                        bands.down = obj()
                    #end if
                #end if

                if '- SPIN UP -' in l:
                    up_spin   = True
                elif '- SPIN DOWN -' in l:
                    up_spin   = False
                    index = -1
                #end if
                              
                if 'number of k points=' in l:
                    try:
                        num_kpoints      = int(l.strip().split()[4])
                    except:
                        print("Number of k-points {0} is not an integer".format(num_kpoints))
                    #end try

                    kpoints_2pi_alat = lines[i+2:i+2+num_kpoints]
                    kpoints_rel      = lines[i+4+num_kpoints:i+4+2*num_kpoints]
                    kpoints_2pi_alat = array([k.strip().split()[4:6] + [k.strip().split()[6][0:-2]] for k in kpoints_2pi_alat], dtype=float)
                    kpoints_rel      = array([k.strip().split()[4:6] + [k.strip().split()[6][0:-2]] for k in kpoints_rel], dtype=float)
                #end if
                if 'bands (ev)' in l:
                    index+=1
                    nfound+=1
                    i_occ = -1
                    j = i
                    while i_occ==-1:
                        j+=1
                        if 'occupation numbers' in lines[j]:
                            i_occ = j
                        #end if
                    #end while
                    seigs = ''
                    for j in range(i+1,i_occ):
                        seigs+=lines[j]
                    #end for
                    seigs = seigs.replace('-',' -') # For cases where the eigenvalues "touch", e.g. -144.9938-144.9938 -84.3023
                    seigs = seigs.strip()
                    eigs = array(seigs.split(),dtype=float)

                    soccs = ''
                    for j in range(i_occ+1,i_occ+1+(i_occ-i)-2):
                        soccs+= lines[j]
                    #end for
                    occs   = array(soccs.split(),dtype=float)
                    bk = obj(
                        index           = index,
                        kpoint_2pi_alat = kpoints_2pi_alat[index],
                        kpoint_rel      = kpoints_rel[index],
                        eigs            = eigs,
                        occs            = occs,
                        pol             = 'none',
                        )
                    band_channel = bands.up
                    if polarized:                  
                        if up_spin:
                            bk.pol = 'up'
                        elif not up_spin:
                            bk.pol = 'down'
                            band_channel = bands.down
                        #end if
                    else:
                        index = nfound -1 
                    #end if
                    band_channel.append(bk)
                    #if nfound==1:
                    #    bands.up = obj(
                    #        eigs = eigs,
                    #        occs = occs
                    #        )
                    #elif nfound==2:
                    #    bands.down = obj(
                    #        eigs = eigs,
                    #        occs = occs
                    #        )
                    #end if
                #end if
            #end for
            vbm        = obj(energy=-1.0e6)
            cbm        = obj(energy=1.0e6)
            direct_gap = obj(energy=1.0e6)
            for band_channel in bands:
                for b in band_channel:
                    e_val  = max(b.eigs[b.occs > 0.5])
                    e_cond = min(b.eigs[b.occs < 0.5])

                    if e_val > vbm.energy:
                        vbm.energy          = e_val
                        vbm.kpoint_rel      = b.kpoint_rel
                        vbm.kpoint_2pi_alat = b.kpoint_2pi_alat
                        vbm.index           = b.index
                        vbm.pol             = b.pol
                        vbm.band_number     = max(where(b.occs > 0.5))
                    #end if
                    if e_cond < cbm.energy:
                        cbm.energy          = e_cond
                        cbm.kpoint_rel      = b.kpoint_rel
                        cbm.kpoint_2pi_alat = b.kpoint_2pi_alat
                        cbm.index           = b.index
                        cbm.pol             = b.pol
                        cbm.band_number     = min(where(b.occs < 0.5))
                    #end if
                    if (e_cond - e_val) < direct_gap.energy:
                        direct_gap.energy          = e_cond - e_val
                        direct_gap.kpoint_rel      = b.kpoint_rel
                        direct_gap.kpoint_2pi_alat = b.kpoint_2pi_alat
                        direct_gap.index           = b.index
                        direct_gap.pol             = [vbm.pol, cbm.pol]
                    #end if
                #end for
            #end for
            electronic_structure = ''
            if (vbm.energy +0.025) >= cbm.energy:
                if vbm.band_number == cbm.band_number:
                    electronic_structure = 'metallic'
                else:
                    electronic_structure = 'semi-metal'
                #end if
            else:
                electronic_structure = 'insulating'
                if not equal(vbm.kpoint_rel, cbm.kpoint_rel).all():
                    indirect_gap = obj(energy=round(cbm.energy-vbm.energy, 3), kpoints=obj(vbm=vbm, cbm=cbm))
                    bands.indirect_gap = indirect_gap
            #end if
            bands.electronic_structure = electronic_structure
            bands.vbm = vbm
            bands.cbm = cbm
            bands.direct_gap = direct_gap
            if nfound>0:
                self.bands = bands
            #end if
            # Kayahan edited --end
        except:
            nx+=1
            if self.info.warn:
                self.warn('band read failed')
            #end if
        #end try

        try:
            # read structures
            structures = obj()
            i=0
            found = False
            cont  = False
            while i<len(lines):
                l = lines[i]
                if l.find('CELL_PARAMETERS')!=-1 and l.strip().startswith('CELL'):
                    conf = obj()
                    axes = []
                    cont = True
                    for d in (0,1,2):
                        i+=1
                        axes.append(array(lines[i].split(),dtype=float))
                    #end for
                    conf.axes = array(axes)
                #end if
                if l.find('ATOMIC_POSITIONS')!=-1:
                    found = True
                    if not cont:
                        conf = obj()
                    #end if
                    atoms = []
                    positions = []
                    i+=1
                    tokens = lines[i].split()

                    while len(tokens)>0 and tokens[0].lower()!='end' and (len(tokens)==4 or (len(tokens)==7 and tokens[-1] in '01')):
                        atoms.append(tokens[0])
                        positions.append(array(tokens[1:4],dtype=float))
                        i+=1
                        tokens = lines[i].split()
                    #end while
                    conf.atoms = atoms
                    conf.positions = array(positions)
                    if 'crystal' in l.lower() and 'axes' in conf:
                        conf.positions = dot(conf.positions,conf.axes)
                    #end if
                    nconf = len(structures)
                    structures[nconf]=conf
                    cont = False
                #end if
                i+=1
            #end while
            if found:
                self.structures = structures
            #end if
        except:
            nx+=1
            if self.info.warn:
                self.warn('structure read failed')
            #end if
        #end try

        #begin added by Yubo "Paul" Yang: cell, stress, pressure and volume
        # 01/21/2016: grab cells in a vc-relax run, one at each optimization step
        # 09/29/2016: obselete after rev7131

        # grab stress, pressure and volume
        try:
            press= 0.
            vol=   0.
            for l in lines:
                if l.find('unit-cell volume')!=-1:
                    #vol = float( l.split('=')[-1].split()[-2] )
                    vol = l.split('=')[-1].split()[-2]
                # end if
                if (l.find('total')!=-1) and (l.find('stress')!=-1):
                    press= l.split('=')[-1]
                # end if
            # end for
            self.pressure = float(press)
            self.volume   = float(vol)
        except:
            nx+=1
            if self.info.warn:
                self.warn('pressure/volume read failed')
            #end if
        #end try

        try:
            stress = []
            nlines = len(lines)
            i=0
            while i<nlines:
                l = lines[i]
                if l.find('total   stress')!=-1:
                    for j in range(3):
                        i+=1
                        stress.append(list(np.array(lines[i].split(),dtype=float)))
                    #end for
                #end if
                i += 1
            #end while
            self.stress=stress
        except:
            nx+=1
            if self.info.warn:
                self.warn('stress read failed')
            #end if
        #end try
        #end added by Yubo "Paul" Yang

        try:
            forces = []
            tot_forces = []
            i=0
            found = False
            nlines = len(lines)
            while i<nlines:
                l = lines[i]
                if l.find('Forces acting on atoms')!=-1:
                    found = True
                    conf = obj()
                    aforces = []
                    found_atom = False
                    for j in range(10):
                        i+=1
                        if i<nlines and 'atom' in lines[i]:
                            found_atom = True
                            break
                        #end if
                    #end for
                    if found_atom:
                        tokens = lines[i].split()
                        while len(tokens)==9 and tokens[4]=='force':
                            aforces.append(tokens[6:])
                            i+=1
                            tokens = lines[i].split()
                        #end while
                    #end if
                    forces.append(aforces)
                    i+=1
                elif 'Total force' in l:
                    tokens = l.split()
                    if len(tokens)==9 and tokens[1]=='force':
                        tot_forces.append(float(tokens[3]))
                    #end if
                #end if
                i+=1
            #end while
            if found:
                self.forces = array(forces,dtype=float)
                self.tot_forces = array(tot_forces)
                max_forces = []
                for f in self.forces:
                    if len(f.shape)==2:
                        max_forces.append((sqrt((f**2).sum(1))).max())
                    #end if
                #end for
                self.max_forces = array(max_forces)
            #end if
        except:
            nx+=1
            if self.info.warn:
                self.warn('force read failed')
            #end if
        #end try

        try:
            tc= 0.
            tw= 0.
            for l in lines:
                if l.find('PWSCF        :')!=-1:
                    t1 = l.split(':')[1].split('CPU')
                    tc = pwscf_time(t1[0])
                    tw = pwscf_time(t1[1].replace('WALL',''))
                    break
                #end if
            #end for
            self.cputime = tc
            self.walltime= tw
        except:
            nx+=1
            if self.info.warn:
                self.warn('time read failed')
            #end if
        #end try


        try:
            # read symmetrized k-points
            nkpoints = None
            i=0
            for l in lines:
                if 'number of k points' in l:
                    tokens = l.replace('=',' ').split()
                    nkpoints = int(tokens[4])
                    break
                #end if
                i+=1
            #end for
            if nkpoints is not None:
                i+=2
                klines_cart = lines[i:i+nkpoints]
                i+=nkpoints+2
                klines_unit = lines[i:i+nkpoints]
                kpoints_cart = []
                for l in klines_cart:
                    tokens = l.replace('= (',':').replace('), wk =',':').split(':')
                    kpoints_cart.append(tokens[1].split())
                #end for
                kpoints_unit = []
                kweights = []
                for l in klines_unit:
                    tokens = l.replace('= (',':').replace('), wk =',':').split(':')
                    kpoints_unit.append(tokens[1].split())
                    kweights.append(tokens[2])
                #end for
                self.kpoints_cart = array(kpoints_cart,dtype=float)
                self.kpoints_unit = array(kpoints_unit,dtype=float)
                self.kweights     = array(kweights,dtype=float)
            #end if
        except:
            nx+=1
            if self.info.warn:
                self.warn('symmetrized kpoint read failed')
            #end if
        #end try
            

        try:
            if pw2c_outfile_name!=None:
                lines = open(os.path.join(path,pw2c_outfile_name),'r').readlines()
                for l in lines:
                    if l.find('Kinetic')!=-1:
                        tokens = l.split()
                        self.K = eval(tokens[5])
                        break
                    #end if
                #end for
            #end if        
        except:
            nx+=1
            if self.info.warn:
                self.warn('pw2casino read failed')
            #end if
        #end try

        if nx>0 and self.info.warn:
            self.warn('encountered an exception, some quantities will not be available')
        #end try

        if self.info.xml:
            self.xmldata = obj(
                data    = None,
                kpoints = None,
                failed  = False
                )
            try:
                cont = self.input.control
                datadir = os.path.join(self.path,cont.outdir,cont.prefix+'.save')
                data_file = os.path.join(datadir,'data-file.xml')
                if not os.path.exists(data_file):
                    datadir = os.path.join(self.path,cont.outdir)
                    data_file = os.path.join(datadir,cont.prefix+'.xml')
                #end if
                data = read_qexml(data_file)
                kpdata = data.root.eigenvalues.k_point
                kpoints = obj()
                for ki,kpd in kpdata.items():
                    kp = obj(
                        kpoint = kpd.k_point_coords,
                        weight = kpd.weight
                        )
                    kpoints[ki]=kp
                    for si,dfile in kpd.datafile.items():
                        efilepath = os.path.join(datadir,dfile.iotk_link)
                        if os.path.exists(efilepath):
                            edata = read_qexml(efilepath)
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
                            spin = obj(
                                units       = units,
                                eigenvalues = edata.root.eigenvalues,
                                occupations = edata.root.occupations
                                )
                            if si==1:
                                kp.up = spin
                            elif si==2:
                                kp.down = spin
                            #end if
                        else:
                            self.xmldata.failed = True
                        #end if
                    #end for
                #end for
                self.xmldata.set(
                    data    = data,
                    kpoints = kpoints
                    )
            except Exception as e:
                if self.info.warn:
                    self.warn('encountered an exception during xml read, this data will not be available\nexception encountered: '+str(e))
                #end if
                self.xmldata.failed = True
            #end try
        #end if
    #end def analyze


    def write_electron_counts(self,filepath=None,return_flag=False):
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
        for kp in kpoints:
            w = kp.weight
            for s,sl in spins.items():
                tot[s] += w*kp[sl].occupations.sum()
            #end for
        #end for
        text = 'total electron counts\n'
        text += '  {0: 3.2f}  {1: 3.2f}  {2: 3.2f}  {3: 3.2f}\n'.format(tot.up+tot.down,tot.up-tot.down,tot.up,tot.down)
        text += '\nkpoint electron counts\n'
        weights = []
        for kp in kpoints:
            weights.append(kp.weight)
        #end for
        weights = array(weights,dtype=float)
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
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        if not return_flag:
            return text
        else:
            return True
        #end if
    #end def write_electron_counts


    def md_statistics(self,equil=None,autocorr=None):
        import numpy as np
        from numerics import simstats,simplestats
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
                v.shape = nb,autocorr
                mean,error = simplestats(v.mean(axis=1))
            #end if
            mds[q] = mean,error
        #end for
        return mds
    #end def md_statistics


    def md_plots(self,show=True):

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
            from structure import Structure
            if filepath==None:
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
                open(filepath,'w').write(movie)
            #end for
        #end for
    #end def make_movie


    def plot_bandstructure(self, filename=None, filepath=None, max_min_e = None, show=False, save=True, show_vbm_cbm=True,k_labels=None):
        if 'bands' in self:
            success = True
            from structure import get_kpath
            if filename==None:
                filename = 'band_structure.pdf'
            if filepath==None:
                filepath = os.path.join(self.abspath,filename)
            else:
                filepath = os.path.join(filepath,filename)
            #end if
            try:
                import matplotlib
                gui_envs = ['GTKAgg','TKAgg','agg','Qt4Agg','WXAgg']
                for gui in gui_envs:
                    try:
                        matplotlib.use(gui,warn=False, force=True)
                        from matplotlib import pyplot
                        success = True
                        break
                    except:
                        continue
                    #end try
                #end for
                from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,rc,savefig,gca,xticks,axvline, scatter
                params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
                              'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
                rcParams.update(params)
            except(ImportError, RuntimeError):
                success = False
            if not success:
                figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text = unavailable('matplotlib.pyplot','figure','plot','xlabel','ylabel','title','show','ylim','legend','xlim','rcParams','savefig','bar','xticks','subplot','grid','setp','errorbar','loglog','semilogx','semilogy','text')
            #end if
            fig    = figure()
            ax     = gca()
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
                y = array(y) - self.bands.vbm.energy
                plot(x, y, 'k')
                if len(self.bands.down) > 0:
                    y = []
                    for bi in self.bands.down:
                        y.append(bi['eigs'][nb])
                    #end for
                    y = array(y) - self.bands.vbm.energy
                    plot(x, y, 'r')
                #end if              
            #end for
            for ln, li in enumerate(labels):
                if li != '':
                    axvline(x[ln], ymin=-100, ymax=100, linewidth=3, color='k')
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
            
            xlim([min(x), max(x)])
            if max_min_e is None:
                ylim(-5, +5)
            else:
                ylim(max_min_e[0],max_min_e[1])
            #end if
            ylabel('Energy (eV)')
            xticks(x, labels)
            ax.tick_params(axis='x', which='both', length=0)
            ax.tick_params(axis='x', which='both', pad=10)
        #end if
        if show_vbm_cbm:
            vbm = self.bands.vbm
            cbm = self.bands.cbm
            for kn, ki in enumerate(self.bands.up):
                if (vbm.kpoint_rel == ki['kpoint_rel']).all():
                    scatter(x[kn], 0, c='green', s=100)
                #end if
                if (cbm.kpoint_rel == ki['kpoint_rel']).all():
                    scatter(x[kn], cbm.energy-vbm.energy, c='r', s=100)
                #end if
            #end for
        #end if
        if save:
            savefig(filename, format='pdf',bbox_inches='tight')
        #end if
        if show:
            show()
        #end if
    #end def plot_bandstructure

#end class PwscfAnalyzer
        
