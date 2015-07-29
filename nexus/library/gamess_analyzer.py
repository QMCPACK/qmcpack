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
            if log.seek('TOTAL NUMBER OF MOS IN VARIATION SPACE',0)!=-1:
                counts.mos = int(log.readtokens()[-1])
            #end if
        except:
            if self.info.exit:
                self.error('log file analysis failed (mo count)')
            #end if
        #end try
        if len(counts)>0:
            self.counts = counts
        #end if
    #end def analyze_log


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
