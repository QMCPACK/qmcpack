
import os
from numpy import array,ndarray,abs
from generic import obj
from developer import DevBase
from debug import *
from simulation import SimulationAnalyzer,Simulation
from gamess_input import GamessInput


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
        info = self.info
        if not info.initialized:
            self.error('cannot perform analysis\nGamessAnalyzer has not been initialized')
        #end if
        files = info.files

        # read the log file
        try:
            output = self.get_output(files.output)
            if output!=None:
                sections = obj()
                insection = False
                sec = None
                prev=''
                lm2 = False
                lm1 = False
                for line in output.splitlines():
                    if insection:
                        sec.append(line)
                    #end if
                    ls = line.strip()
                    chars = list(set(ls))
                    lm0 = len(chars)==1 and chars[0]=='-'
                    if lm0 and lm2:
                        secname = ''
                        tokens = prev.split()
                        for token in tokens:
                            secname+=token+'_'
                        #end for
                        secname = secname.replace('-','_').lower().split(',')[0].split('=')[0].split('(')[0].rstrip('_').replace('2','two')
                        insection = secname!=''
                        if insection:
                            sec = []
                            sections[secname] = sec
                        #end if
                    #end if
                    lm2=lm1
                    lm1=lm0
                    prev = ls
                #end for
                #for secname in list(sections.keys()):
                #    sec = ''
                #    for line in sections[secname]:
                #        sec+=line+'\n'
                #    #end for
                #    sections[secname]=sec
                ##end for
                if 'energy_components' in sections:
                    energies = obj()
                    for line in sections.energy_components:
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
                            energies[name]=value
                        #end if
                    #end for
                    self.energy_components = energies
                #end if
            #end if
        except:
            if self.info.exit:
                self.error('log file analysis failed')
            #end if
        #end try

        # read the punch file
        try:
            text = self.get_output(files.punch)
            if text!=None:
                punch = obj(norbitals=0)
                group_name = None
                group_text = ''
                new = False
                for line in text.splitlines():
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
                        group_text += line+'\n'
                    #end if
                    new = False
                #end for
                if len(punch)>0:
                    self.punch=punch
                #end if
                if 'vec' in punch:
                    norbs = 0
                    for line in punch.vec.splitlines()[1:-1]:
                        norbs = max(norbs,int(line[:2]))
                    #end for
                    punch.norbitals = norbs
                #end if
            #end if
        except:
            if self.info.exit:
                self.error('punch file analysis failed')
            #end if
        #end try
    #end def analyze


    def get_output(self,filename):
        outfile = os.path.join(self.info.path,filename)
        if os.path.exists(outfile):
            return open(outfile,'r').read()
        elif os.path.exists(filename):
            return open(filename,'r').read()
        elif self.info.exit:
            self.error('output file does not exist at either of the locations below:\n  {0}\n  {1}'.format(outfile,filename))
        else:
            return None
        #end if
    #end def get_output
#end class GamessAnalyzer
