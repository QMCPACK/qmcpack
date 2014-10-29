import os
from numpy import array,fromstring,sqrt
from generic import obj
from unit_converter import convert
from periodic_table import PeriodicTable
from simulation import SimulationAnalyzer,Simulation
from pwscf_input import PwscfInput
from pwscf_data_reader import read_qexml
from debug import *
import code

pt = PeriodicTable()
elements = set(pt.elements.keys())




def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
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
    def __init__(self,arg0=None,infile_name=None,outfile_name=None,pw2c_outfile_name=None,analyze=False,xml=False,warn=False):
        if isinstance(arg0,Simulation):
            sim = arg0
            path = sim.locdir
            self.infile_name = sim.infile
            self.outfile_name= sim.outfile
        elif arg0!=None:
            path = arg0
            self.infile_name = infile_name
            self.outfile_name = outfile_name
        else:
            return
        #end if
        self.path = path
        self.abspath = os.path.abspath(path)
        self.pw2c_outfile_name = pw2c_outfile_name

        self.info = obj(xml=xml,warn=warn)

        self.input = PwscfInput(os.path.join(self.path,self.infile_name))

        if analyze:
            self.analyze()
        #end if
    #end def __init__

    
    def analyze(self):
        try:
            path = self.path
            infile_name = self.infile_name
            outfile_name = self.outfile_name
            pw2c_outfile_name = self.pw2c_outfile_name

            lines = open(os.path.join(path,outfile_name),'r').read().splitlines()

            energies = []
            for l in lines:
                if l.find('!  ')!=-1:
                    energies.append( eval( l.split('=')[1].split()[0] ) )
                #end if
            #end for
            if len(energies)==0:
                self.E = 0.0
            else:
                self.E = energies[-1]
            #end if
            self.energies = array(energies)

            # get bands and occupations
            nfound = 0
            bands = obj()
            for i in range(len(lines)):
                l = lines[i]
                if 'bands (ev)' in l:
                    nfound+=1
                    i_occ = -1
                    j = i
                    while i_occ==-1:
                        j+=1
                        if 'occupation numbers' in lines[j]:
                            i_occ = j
                        #end if
                    #end while
                    try:
                        seigs = ''
                        for j in range(i+1,i_occ):
                            seigs+=lines[j]
                        #end for
                        seigs = seigs.strip()
                        eigs = array(seigs.split(),dtype=float)

                        soccs = ''
                        for j in range(i_occ+1,i_occ+1+(i_occ-i)-2):
                            soccs+= lines[j]
                        #end for
                        occs = array(soccs.split(),dtype=float)
                    except Exception:
                        eigs = array([])
                        occs = array([])
                        if self.info.warn:
                            self.warn('band read failed, line: '+l)
                        #end if
                    #end try
                    

                    if nfound==1:
                        bands.up = obj(
                            eigs = eigs,
                            occs = occs
                            )
                    elif nfound==2:
                        bands.down = obj(
                            eigs = eigs,
                            occs = occs
                            )
                    #end if
                #end if
                #try:
                #    seigs = ''
                #    for j in range(i+1,i_occ):
                #        seigs+=lines[j]
                #    #end for
                #    seigs = seigs.strip()
                #    eigs = array(seigs.split(),dtype=float)
                #
                #    soccs = ''
                #    for j in range(i_occ+1,i_occ+1+(i_occ-i)-2):
                #        soccs+= lines[j]
                #    #end for
                #    occs = array(soccs.split(),dtype=float)
                #except Exception:
                #    eigs = array([])
                #    occs = array([])
                ##end try

                if nfound==1:
                    bands.up = obj(
                        eigs = eigs,
                        occs = occs
                        )
                elif nfound==2:
                    bands.down = obj(
                        eigs = eigs,
                        occs = occs
                        )
                #end if
            #end for
            if nfound>0:
                self.bands = bands
            #end if


            structures = obj()
            i=0
            found = False
            while i<len(lines):
                l = lines[i]
                if l.find('ATOMIC_POSITIONS')!=-1:
                    found = True
                    conf = obj()
                    atoms = []
                    positions = []
                    i+=1
                    tokens = lines[i].split()

                    while len(tokens)>0 and tokens[0] in elements and (len(tokens)==4 or (len(tokens)==7 and tokens[-1] in '01')):
                        atoms.append(tokens[0])
                        positions.append(array(tokens[1:4],dtype=float))
                        i+=1
                        tokens = lines[i].split()
                    #end while
                    conf.atoms = atoms
                    conf.positions = array(positions)
                    nconf = len(structures)
                    structures[nconf]=conf
                #end if
                i+=1
            #end while
            if found:
                self.structures = structures
            #end if

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
                    if i<nlines:
                        tokens = lines[i].split()
                        if len(tokens)==9 and tokens[1]=='force':
                            tot_forces.append(float(tokens[3]))
                        #end if
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
        except Exception:
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
                data = read_qexml(os.path.join(datadir,'data-file.xml'))
                kpdata = data.root.eigenvalues.k_point
                kpoints = obj()
                for ki,kpd in kpdata.iteritems():
                    kp = obj(
                        kpoint = kpd.k_point_coords,
                        weight = kpd.weight
                        )
                    kpoints[ki]=kp
                    for si,dfile in kpd.datafile.iteritems():
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
            except Exception,e:
                self.warn('encountered an exception during xml read, this data will not be available\nexception encountered: '+str(e))
                self.xmldata.failed = True
            #end try
        #end if
    #end def analyze


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
                

#end class PwscfAnalyzer
        
