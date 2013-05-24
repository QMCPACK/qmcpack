import os
from numpy import array,fromstring,sqrt
from generic import obj
from unit_converter import convert
from periodic_table import PeriodicTable
from simulation import SimulationAnalyzer,Simulation
from pwscf_input import PwscfInput

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
    def __init__(self,arg0=None,infile_name=None,outfile_name=None,pw2c_outfile_name=None,analyze=False):
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

        self.input = PwscfInput(os.path.join(self.path,self.infile_name))

        if analyze:
            self.analyze()
        #end if
    #end def __init__

    
    def analyze(self):
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


        found = False
        for i in range(len(lines)):
            l = lines[i]
            if l.find(' bands')!=-1:
                lb = lines[i+2]
                ls = lb.split()
                all_number = True
                for sval in ls:
                    all_number = all_number and is_number(sval)
                #end for
                if all_number:
                    bands = array(ls,float)
                else:
                    bands = None
                #end if
                found = True
            #end if
        #end for
        if found:
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

            print structures


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
        
