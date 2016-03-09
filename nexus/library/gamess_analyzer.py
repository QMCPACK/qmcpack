##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  gamess_analyzer.py                                                #
#    Support for analysis of GAMESS output data.                     #
#                                                                    #
#  Content summary:                                                  #
#    GamessAnalyzer                                                  #
#      Analyzer class for the GAMESS code.                           #
#      Reads numerical data from GAMESS output logs.                 #
#                                                                    #                                        
#====================================================================#


import os
from numpy import array,ndarray,abs
from generic import obj
from developer import DevBase
from fileio import TextFile
from debug import *
from simulation import SimulationAnalyzer,Simulation
from gamess_input import GamessInput


def assign_value(host,dest,file,string):
    if file.seek(string)!=-1:
        host[dest] = float(file.readtokens()[-1])
    #end if
#end def assign_value


class GamessAnalyzer(SimulationAnalyzer):

    lset_full = 'spdfg'

    lxyz = obj(
        s = set('s'.split()),
        p = set('x y z'.split()),
        d = set('xx yy zz xy xz yz'.split()),
        f = set('xxx yyy zzz xxy xxz yyx yyz zzx zzy xyz'.split()),
        g = set('xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy'.split()),
        )

    lxyz_reverse = obj()
    for l,xyzs in lxyz.iteritems():
        for xyz in xyzs:
            lxyz_reverse[xyz] = l
        #end for
    #end for



    def __init__(self,arg0=None,prefix=None,analyze=False,exit=False,**outfilenames):
        self.info = obj(
            exit   = exit,
            path   = None,
            input  = None,
            prefix = None,
            files  = obj(),
            initialized = False
            )
        infile = None
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = os.path.join(sim.locdir,sim.infile)
        else:
            infile = arg0
        #end if
        if infile!=None:
            info = self.info
            info.path = os.path.dirname(infile)
            info.input = GamessInput(infile)
            infilename = os.path.split(infile)[1]
            if prefix is None:
                prefix = infilename.rsplit('.',1)[0]
            #end if
            info.prefix = prefix
            files = info.files
            for file,unit in GamessInput.file_units.iteritems():
                files[file.lower()] = '{0}.F{1}'.format(prefix,str(unit).zfill(2))
            #end for
            files.input  = infilename
            files.output = '{0}.out'.format(prefix)
            for name,filename in outfilenames:
                if name in files:
                    files[name] = filename
                else:
                    self.error('unknown GAMESS file: {0}'.format(name))
                #end if
            #end for
            info.initialized = True
            if analyze:
                self.analyze()
            #end if
        #end if
    #end def __init__


    def analyze(self):
        if not self.info.initialized:
            self.error('cannot perform analysis\nGamessAnalyzer has not been initialized')
        #end if
        self.analyze_log()
        self.analyze_punch()
    #end def analyze


    def get_output(self,filetag):
        filename = self.info.files[filetag]
        outfile = os.path.join(self.info.path,filename)
        if os.path.exists(outfile):
            filepath = outfile
        elif os.path.exists(filename):
            filepath = filename
        elif self.info.exit:
            self.error('output file does not exist at either of the locations below:\n  {0}\n  {1}'.format(outfile,filename))
        else:
            return None
        #end if
        try:
            return TextFile(filepath)
        except:
            return None
        #end try
    #end def get_output


    def analyze_log(self):
        # read the log file
        log = self.get_output('output')

        # try to get the energy components
        energy = obj()
        try:
            self.read_energy_components(log,energy)
        except:
            if self.info.exit:
                self.error('log file analysis failed (energy components)')
            #end if
        #end try
        if len(energy)>0:
            self.energy = energy
        #end if

        # try to get the orbital count from the log file
        counts = obj()
        try:
            if log.seek('TOTAL NUMBER OF BASIS SET SHELLS',0)!=-1:
                counts.shells = int(log.readtokens()[-1])
            #end if
            if log.seek('NUMBER OF CARTESIAN ATOMIC ORBITALS',0)!=-1:
                counts.cao = int(log.readtokens()[-1])
            #end if
            if log.seek('TOTAL NUMBER OF MOS IN VARIATION SPACE',1)!=-1:
                counts.mos = int(log.readtokens()[-1])
            #end if
        except:
            if self.info.exit:
                self.error('log file analysis failed (counts)')
            #end if
        #end try
        if len(counts)>0:
            self.counts = counts
        #end if

        # try to get the up/down orbitals
        if 'counts' in self and False: # don't read orbitals, large
            orbitals = obj()
            try:
                self.read_orbitals(log,orbitals,'up'  ,'-- ALPHA SET --')
                self.read_orbitals(log,orbitals,'down','-- BETA SET --' )
            except:
                if self.info.exit:
                    self.error('log file analysis failed (orbitals)')
                #end if
            #end try
            if len(orbitals)>0:
                self.orbitals = orbitals
            #end if
        #end if

        # try to get the mulliken/lowdin populations in each ao
        if 'counts' in self:
            ao_populations = obj()
            try:
                self.read_ao_populations(log,ao_populations)
            except:
                if self.info.exit:
                    self.error('log file analysis failed (ao populations)')
                #end if
            #end try
            if len(ao_populations)>0:
                self.ao_populations = ao_populations
            #end if
        #end if
    #end def analyze_log


    def read_energy_components(self,log,energy):
        if log!=None and log.seek('ENERGY COMPONENTS',0)!=-1:
            for n in xrange(18):
                line = log.readline()
                if '=' in line and 'ENERGY' in line:
                    nameline,value = line.split('=')
                    tokens = nameline.lower().split()
                    name = ''
                    for token in tokens:
                        if token!='energy':
                            name += token.replace('-','_')+'_'
                        #end if
                    #end for
                    name = name[:-1]
                    value = float(value.strip())
                    energy[name]=value
                #end if
            #end for
        #end if
        if log!=None and log.seek('COUPLED-CLUSTER ENERGY',0)!=-1:
            line = log.readline()
            energy['ccsd(t)'] = float( line.split()[-1] )
        # end if

    #end def read_energy_components


    def read_orbitals(self,log,orbs,spec,header):
        success = True
        cao_tot   = self.counts.cao
        mos_tot   = self.counts.mos
        mos_found = 0
        eigenvalue   = []
        symmetry     = []
        coefficients = []
        have_basis   = False
        element      = []
        spec_index   = []
        angular      = []
        if log.seek(header,0)!=-1:
            imax = mos_tot
            jmax = 100
            i=0
            while mos_found<mos_tot and i<imax:
                j=0
                norbs = -1
                while j<jmax:
                    line = log.readline()
                    if line.strip().replace(' ','').isdigit():
                        found_orbs = True
                        norbs = len(line.split())
                        break
                    #end if
                    j+=1
                #end while
                if j==jmax:
                    success=False
                    if self.info.exit:
                        self.error('could not find start of orbitals for {0}\nnumber of orbitals read successfully: {1}\nnumber of orbitals not read: {2}'.format(header,mos_found,mos_tot-mos_found))
                    else:
                        break
                    #end if
                #end if
                eigenvalue.extend(log.readtokens())
                symmetry.extend(log.readtokens())
                coeff = []
                for icao in xrange(cao_tot):
                    tokens = log.readtokens()
                    if not have_basis:
                        e = tokens[1]
                        element.append(e[0].upper()+e[1:].lower())
                        spec_index.append(tokens[2])
                        angular.append(tokens[3])
                    #end if
                    coeff.append(tokens[4:])
                #end for
                coefficients.extend(array(coeff,dtype=float).T)
                if not have_basis:
                    stype = []
                    ptype = []
                    dtype = []
                    ftype = []
                    for ia in xrange(len(angular)):
                        a = angular[ia].lower()
                        if a in GamessAnalyzer.stypes:
                            stype.append(ia)
                        elif a in GamessAnalyzer.ptypes:
                            ptype.append(ia)
                        elif a in GamessAnalyzer.dtypes:
                            dtype.append(ia)
                        elif a in GamessAnalyzer.ftypes:
                            ftype.append(ia)
                        elif self.info.exit:
                            self.error('unrecognized angular type: {0}'.format(angular[ia]))
                        #end if
                    #end for
                #end if
                have_basis = True
                mos_found = len(eigenvalue)
                i+=1
            #end while
            if i==imax:
                success=False
                if self.info.exit:
                    self.error('orbital read failed for {0}\nnumber of orbitals read successfully: {1}\nnumber of orbitals not read: {2}'.format(header,mos_found,mos_tot-mos_found))
            #end if
            if success:
                orbs[spec] = obj(
                    eigenvalue   = array(eigenvalue  ,dtype=float),
                    symmetry     = array(symmetry    ,dtype=str  ),
                    coefficients = array(coefficients,dtype=float),
                    basis = obj(
                        element    = array(element   ,dtype=str),
                        spec_index = array(spec_index,dtype=int),
                        angular    = array(angular   ,dtype=str),
                        stype      = array(stype     ,dtype=int),
                        ptype      = array(ptype     ,dtype=int),
                        dtype      = array(dtype     ,dtype=int),
                        ftype      = array(ftype     ,dtype=int),
                        )
                    )
            #end if
            return success
        else:
            return False
        #end if
    #end def read_orbitals


    def read_ao_populations(self,log,ao_populations):
        cao_tot   = self.counts.cao
        if log.seek('-- POPULATIONS IN EACH AO --',0)!=-1:
            log.readline()
            log.readline()
            mulliken   = []
            lowdin     = []
            element    = []
            spec_index = []
            angular    = []
            linds = obj()
            for l in GamessAnalyzer.lset_full:
                linds[l]=[]
            #end for
            for icao in xrange(cao_tot):
                tokens = log.readtokens()
                e = tokens[1]
                element.append(e[0].upper()+e[1:].lower())
                if len(tokens)==6:
                    spec_index.append(tokens[2])
                    angular.append(tokens[3])
                    mulliken.append(tokens[4])
                    lowdin.append(tokens[5])
                elif len(tokens)==5:
                    spec_index.append(tokens[2][:-4])
                    angular.append(tokens[2][-4:])
                    mulliken.append(tokens[3])
                    lowdin.append(tokens[4])
                #end if
            #end for
            for ia in xrange(len(angular)):
                a = angular[ia].lower()
                if a in GamessAnalyzer.lxyz_reverse:
                    l = GamessAnalyzer.lxyz_reverse[a]
                    linds[l].append(ia)
                elif self.info.exit:
                    self.error('unrecognized angular type: {0}'.format(angular[ia]))
                #end if
            #end for
            mulliken = array(mulliken,dtype=float)
            lowdin   = array(lowdin  ,dtype=float)
            for l,lind in linds.iteritems():
                linds[l] = array(lind,dtype=int)
            #end for

            mulliken_shell = []
            lowdin_shell = []
            shell = obj()
            n=0
            for l in GamessAnalyzer.lset_full:
                if l in GamessAnalyzer.lxyz:
                    lind = linds[l]
                    nxyz = len(GamessAnalyzer.lxyz[l])
                    nshell = len(lind)/nxyz
                    shellinds = []
                    for ns in range(nshell):
                        inds = lind[ns*nxyz:(ns+1)*nxyz]
                        mulliken_shell.append(mulliken[inds].sum())
                        lowdin_shell.append(lowdin[inds].sum())
                        shellinds.append(n)
                        n+=1
                    #end for
                    shell[l] = array(shellinds,dtype=int)
                #end if
            #end for
            mulliken_shell = array(mulliken_shell)
            lowdin_shell   = array(lowdin_shell)
            mulliken_angular = obj()
            lowdin_angular   = obj()
            for l,shellinds in shell.iteritems():
                mulliken_angular[l] = mulliken_shell[shellinds].sum()
                lowdin_angular[l]   = lowdin_shell[shellinds].sum()
            #end for
            basis = obj(
                element    = array(element   ,dtype=str),
                spec_index = array(spec_index,dtype=int),
                angular    = array(angular   ,dtype=str),
                )
            basis.transfer_from(linds)
            ao_populations.set(
                mulliken         = mulliken,
                lowdin           = lowdin,
                mulliken_shell   = mulliken_shell,
                lowdin_shell     = lowdin_shell,
                mulliken_angular = mulliken_angular,
                lowdin_angular   = lowdin_angular,
                basis            = basis,
                shell            = shell,
                )
        #end if
    #end def read_ao_populations


    def analyze_punch(self):
        # read the punch file
        try:
            text = self.get_output('punch')
            if text!=None:
                #text = text.read()
                punch = obj(norbitals=0)
                group_name = None
                group_text = ''
                new = False
                #for line in text.splitlines():
                for line in text:
                    ls = line.strip()
                    if len(ls)>0 and ls[0]=='$':
                        if ls=='$END':
                            punch[group_name] = group_text
                            group_name = None
                            group_text = ''
                        else:
                            group_name = ls[1:].lower()
                            group_text = ''
                            new        = True
                        #end if
                    #end if
                    if group_name!=None and not new:
                        group_text += line#+'\n'
                    #end if
                    new = False
                #end for
                if len(punch)>0:
                    self.punch=punch
                #end if
                if 'vec' in punch:
                    # count the orbitals in vec
                    norbs = 0
                    ind_counts = dict()
                    nbasis = 0
                    for line in punch.vec.splitlines()[1:-1]:
                        ind = int(line[:2])
                        iln = int(line[2:5])
                        nln = max(nbasis,iln)
                        if ind not in ind_counts:
                            ind_counts[ind] = 1
                        else:
                            ind_counts[ind] += 1
                        #end if
                    #end for
                    all_double = True
                    norbs = 0
                    for ind,cnt in ind_counts.iteritems():
                        count = cnt/nln
                        norbs+=count
                        all_double = all_double and count%2==0
                    #end for
                    if all_double:
                        norbs/=2
                    #end if
                    punch.norbitals = norbs
                #end if
            #end if
        except:
            if self.info.exit:
                self.error('punch file analysis failed')
            #end if
        #end try
    #end def analyze_punch
#end class GamessAnalyzer
