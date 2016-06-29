##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  gamess_input.py                                                   #
#    Support for GAMESS input file I/O                               #
#                                                                    #
#  Content summary:                                                  #
#    GamessInput                                                     #
#      Input class for the GAMESS code.                              #
#      Capable of reading/writing arbitrary GAMESS input files.      #
#                                                                    #
#    generate_gamess_input                                           #
#      User function to create arbitrary GAMESS input.               #
#                                                                    #
#    KeywordGroup                                                    #
#      Represents an arbitary keyword group in the input file.       #
#                                                                    #
#    KeywordSpecGroup                                                #
#      Base class for specialized keyword groups.                    #
#      Derived classes enforce the keyword specification.            #
#      See ContrlGroup, SystemGroup, GuessGroup, ScfGroup,           #
#        McscfGroup, DftGroup, GugdiaGroup, DrtGroup, CidrtGroup,    #
#        and DetGroup                                                #
#                                                                    #
#    FormattedGroup                                                  #
#      Represents strict machine-formatted input groups.             #
#                                                                    #
#====================================================================#



import os
from numpy import array,ndarray,abs
from generic import obj
from periodic_table import pt
from developer import DevBase
from nexus_base import nexus_noncore
from simulation import SimulationInput
from debug import *



class GIbase(DevBase):
    def message(self,msg,**kwargs):
        self.error(msg,**kwargs)
    #end def message
#end class GIbase



class GIarray(GIbase):
    def __init__(self,d):
        for n,v in d.iteritems():
            if not isinstance(n,int):
                self.error("keys must be integers\nattempted to initialize array from input provided: {0}\nnote that dict's are used only for arrays".format(d))
            #end if
            if isinstance(v,(tuple,list,ndarray)):
                nv = array(v,type(v[0]))
            else:
                nv = array([v],type[v])
            #end if
            self[n]=nv
        #end for
    #end def __init__
#end class GIarray



class Group(GIbase):
    def __init__(self,text=None,**kwargs):
        if text!=None:
            self.read(text)
        #end if
        self.set(**kwargs)
    #end def __init__

    def read(self,text):
        self.not_implemented()
    #end def read

    def write(self,text):
        self.not_implemented()
    #end def read
#end class Group



class KeywordGroup(Group):

    linewrap   = 70
    lineindent = '         '


    booldict = {'.true.':True,'.TRUE.':True,'.t.':True,'.T.':True,
                '.false.':False,'.FALSE.':False,'.f.':False,'.F.':False}

    def readval(self,val):
        fail = False
        if val in self.booldict:
            v = self.booldict[val]
        else:
            try:
                v = int(val)
            except:
                try:
                    v = float(val.replace('d','e'))
                except:
                    #val = val.replace(',',' ')
                    if ' ' in val:
                        val = val.split()
                        try:
                            v = array(val,dtype=int)
                        except:
                            try:
                                v = array(val,dtype=float)
                            except:
                                try:
                                    v = array(val,dtype=str)
                                except:
                                    fail = True
                                #end try
                            #end try
                        #end try
                    else:
                        v = val
                    #end if
                #end try
            #end try
        #end if
        if fail:
            self.error('failed to read value: "{0}"'.format(val))
        #end if
        return v
    #end def readval
            

    def read(self,text):
        tokens = text.replace(',',' ').split()
        for token in tokens:
            if '=' in token:
                var,val = token.split('=')
                var = var.lower()
                val = val.lower()
                self[var]=val
            else:
                self[var]+=' '+token.lower()
            #end if
        #end for
        vars = list(self.keys())
        for var in vars:
            val = self.readval(self[var])
            if not '(' in var:
                self[var] = val
            else:
                del self[var]
                var,index = var.replace('(',' ').replace(')','').split()
                index = int(index)
                if not var in self:
                    arr = GIarray()
                    arr[index] = val
                    self[var] = arr
                else:
                    self[var][index]=val
                #end if
            #end if
        #end for
    #end def read


    def writeval(self,val):
        if isinstance(val,bool):
            if val:
                sval = '.true.'
            else:
                sval = '.false.'
            #end if
        elif isinstance(val,str):
            sval = val
        elif isinstance(val,int):
            sval = str(val)
        elif isinstance(val,float):
            sval = str(val).replace('e','d')
        elif isinstance(val,(ndarray,list)):
            sval = ''
            for v in val:
                vs = str(v)+','
                if len(sval)+len(vs)<self.linewrap:
                    sval+=vs
                else:
                    sval+='\n'+self.lineindent+vs
                #end if
            #end for
            sval = sval[0:-1]
        else:
            self.error('unknown type encountered on write: {0}'.format(val))
        #end if
        return sval
    #end def writeval


    def write(self,name):
        text = ''
        line = ' ${0:<6} '.format(name)
        for var in sorted(self.keys()):
            val = self[var]
            if not isinstance(val,GIarray):
                vtext='{0}={1} '.format(var,self.writeval(val))
                if len(line)+len(vtext) < self.linewrap:
                    line+=vtext
                else:
                    text+=line+'\n'
                    line = self.lineindent+vtext
                #end if
            else:
                for n in sorted(val.keys()):
                    vtext = '{0}({1})={2} '.format(var,n,self.writeval(val[n]))
                    if len(line)+len(vtext) < self.linewrap:
                        line+=vtext
                    else:
                        text+=line+'\n'
                        line = self.lineindent+vtext
                    #end if
                #end for
            #end if
        #end for
        text += line+' $end\n'
        return text
    #end def write
#end class KeywordGroup



class CardGroup(Group):

    #input spec page numbers
    # ecp 287
    # data 37

    def readval(self,val):
        try:
            v = int(val)
        except:
            try:
                v = float(val.replace('d','e'))
            except:
                v = val
            #end try
        #end try
        return v
    #end def readval


    def read_tokens(self,line):
        tokens = []
        for token in line.split():
            tokens.append(self.readval(token))
        #end for
        return tokens
    #end def read_tokens


    def read_line_tokens(self,text):
        line_tokens = []
        for line in text.splitlines():
            line_tokens.append(self.read_tokens(line))
        #end for
        return line_tokens
    #end def read_line_tokens


    def append_text(self,text):
        for tokens in self.read_line_tokens(text):
            self.append(tokens)
        #end for
    #end def append_text

            
    def append_list(self,lst):
        for tokens in lst:
            self.append(tokens)
        #end for
    #end def append_list


    def read(self,inp):
        self.clear()
        if isinstance(inp,str):
            self.append_text(inp)
        elif isinstance(inp,list):
            self.append_list(inp)
        #end if
    #end def read


    def writeval(self,val):
        if isinstance(val,float):
            sval = str(val).replace('e','d')
            if len(sval)>8 and abs(val)>=10.0:
                sval = '{0:16.8e}'.format(val).replace('e','d')
            #end if
        else:
            sval = str(val)
        #end if
        return sval
    #end def writeval


    def write(self,name):
        text = ' ${0}\n'.format(name)
        contents = ''
        for n in range(len(self)):
            for token in self[n]:
                contents += self.writeval(token)+' '
            #end for
            contents += '\n'
        #end for
        text+= contents.lstrip()
        text+=' $end\n'
        return text
    #end def write


    def list(self):
        lst = []
        for n in range(len(self)):
            lst.append(self[n])
        #end for
        return lst
    #end def list
#end class CardGroup



class FormattedGroup(Group):
    def read(self,text):
        self.text = str(text)
    #end def read

    def write(self,name):
        #return ' ${0}\n{1} $END\n'.format(name.upper(),self.text.lstrip())
        return ' ${0}\n{1} $END\n'.format(name.upper(),self.text)
    #end def write
#end class FormattedGroup





# detailed keyword specification groups to check names and types of keyword inputs

class KeywordSpecGroup(KeywordGroup):
    keywords = set()
    integers = set()
    reals    = set()
    bools    = set()
    strings  = set()
    arrays   = set()
    allowed_values = obj()

    def is_consistent(self):
        return len(set(self.keys())-self.keywords)==0
    #end def is_consistent

    def is_valid(self):
        valid = self.is_consistent()
        for name,val in self.iteritems():
            if name in self.allowed_values:
                if isinstance(val,str):
                    val = val.lower()
                #end if
                valid &= val in self.allowed_values[name]
            #end if
        #end for
        return valid
    #end def is_valid
#end class KeywordSpecGroup
    


class ContrlGroup(KeywordSpecGroup):
    keywords = set([
            'scftyp','dfttyp','tddft' ,'vbtyp' ,'mplevl','cityp' ,'cctyp' ,
            'cimtyp','relwfn','runtyp','numgrd','exetyp','icharg','mult'  ,
            'coord' ,'units' ,'nzvar' ,'pp'    ,'local' ,'ispher','qmttol',
            'maxit' ,'molplt','pltorb','aimpac','friend','nfflvl','nprint',
            'nosym' ,'etollz','inttyp','grdtyp','normf' ,'normp' ,'itol'  ,
            'icut'  ,'iskprp','irest' ,'geom'  ,'ecp'   ,'casino'
            ])
    integers = set([
            'mplevl','icharg','mult' ,'nzvar'  ,'ispher','maxit' ,'nfflvl',
            'nprint','nosym' ,'normf','normp'  ,'itol'  ,'icut'  ,'iskprp',
            'irest'
            ])
    reals    = set(['qmttol' ,'etollz'])
    bools    = set(['numgrd' ,'molplt','pltorb','aimpac','casino'])
    strings  = set([
            'scftyp','dfttyp','tddft' ,'vbtyp' ,'cityp' ,'cctyp' ,'cimtyp',
            'relwfn','runtyp','exetyp','coord' ,'units' ,'pp'    ,'local' ,
            'friend','inttyp','grdtyp','geom'  ,'ecp'
            ])

    allowed_values = obj(
        scftyp = set(['rhf','uhf','rohf','gvb','mcscf','none']),
        dfttyp = set(['none','slater','becke','gill','optx','pw91x','pbex',
                      'vwn','vwn3','vwn1rpa','pz81','p86','lyp','pw91c','pbec',
                      'op','svwn','wvwn1rpa','blyp','bop','bp86','gvwn','gpw91',
                      'pbevwn','pbeop','olyp','pw91','pbe','edf1','pbe','revpbe',
                      'rpbe','pbesol','hcth93','hcth120','hcth147','hcth407',
                      'sogga','mohlyp','b97-d','sogga11','bhhlyp','b3pw91',
                      'b3lyp','b3lypv1r','b3lypv3','b3p86','b3p86v1r','b3p86v5',
                      'b97','b97-1','b97-2','b97-3','b97-k','b98','pbe0','x3lyp',
                      'sogga11x','camb3lyp','wb97','wb97x','wb97x-d','b2plyp',
                      'wb97x-2','wb97x-2l','vs98','pkzb','thcth','thcthhyb','bmk',
                      'tpss','tpssh','tpssm','revtpss','dldf','m05','m05-2x',
                      'm06','m06-l','m06-2x','m06-hf','m08-hx','m08-s0','m11','m11-l',
                      'none','xalpha','slater','becke','depristo','cama','half',
                      'vwn','pwloc','lyp','bvwn','blyp','bpwloc','b3lyp','camb',
                      'xvwn','xpwloc','svwn','spwloc','wigner','ws','wigexp']),
        tddft  = set(['none','excite','spnflp']),
        vbtyp  = set(['none','vb2000']),
        mplevl = set([0,2]),
        cityp  = set(['none','cis','sfcis','aldet','ormas','fsoci','genci','guga']),
        cctyp  = set(['none','lccd','ccd','ccsd','ccsd(t)','r-cc','cr-cc','cr-ccl',
                      'ccsd(tq)','cr-cc(q)','eom-ccsd','cr-eom','cr-eoml','ip-eom2',
                      'ip-eom2','ip-eom3a','ea-eom2','ea-eom3a']),
        cimtyp = set(['none','secim','decim','gsecim']),
        relwfn = set(['none','iotc','dk','resc','nesc']),
        runtyp = set(['energy','gradient','hessian','gamma','optimize','trudge',
                      'sadpoint','mex','conical','irc','vscf','drc','md','globop',
                      'optfmo','gradextr','surface','comp','g3mp2','prop','raman',
                      'nacme','nmr','eda','qmefpea','transitn','ffield','tdhf',
                      'tdhfx','makefp','fmo0']),
        exetyp = set(['run','check']),
        coord  = set(['unique','hint','prinaxis','zmt','zmtmpc','fragonly']),
        units  = set(['angs','bohr']),
        pp     = set(['none','read','sbkjc','hw','mcp']),
        local  = set(['none','boys','ruednbrg','pop','svd']),
        ispher = set([-1,0,1]),
        friend = set(['hondo','meldf','gamessuk','gaussian','all']),
        nfflvl = set([2,3]),
        nprint = set([-7,-6,-5,-4,-3,-2,1,2,3,4,5,6,7,8,9]),
        nosym  = set([0,1]),
        inttyp = set(['best','rotaxis','eric','rysquad']),
        grdtyp = set(['best rsyquad']),
        normf  = set([0,1]),
        normp  = set([0,1]),
        iskprp = set([0,1]),
        irest  = set([-1,0,1,2,3,4]),
        geom   = set(['input','daf']),
        )
#end class ContrlGroup



class SystemGroup(KeywordSpecGroup):
    keywords = set(['mwords','memddi','timlim','parall','kdiag','corefl',
                    'baltyp','mxseq2','mxseq3','nodext','iosmp','modio' ,
                    'memory'])

    integers = set(['mwords','memddi','kdiag','mxseq2','mxseq3','modio','memory'])
    reals    = set(['timlim'])
    bools    = set(['parall','corefl'])
    strings  = set(['baltyp'])
    arrays   = set(['nodext','iosmp'])

    allowed_values = obj(
        kdiag  = set([0,1,2,3]),
        baltyp = set(['slb','dlb','loop','nxtval']),
        modio  = set([1,2,4,8,15]),
        )
#end class SystemGroup



class GuessGroup(KeywordSpecGroup):
    keywords = set(['guess' ,'prtmo' ,'punmo' ,'mix' ,'norb','norder','iorder',
                    'jorder','insorb','purify','tolz','tole','symden'])

    integers = set(['norb','norder','insorb'])
    reals    = set(['tolz','tole'])
    bools    = set(['prtmo','punmo','mix','purify','symden'])
    strings  = set(['guess'])
    arrays   = set(['iorder','jorder'])

    allowed_values = obj(
        guess  = set(['huckel','hcore','moread','rdmini','mosaved','skip','fmo','hucsub','dmread']),
        norder = set([0,1]),
        )
#end class GuessGroup



class ScfGroup(KeywordSpecGroup):
    keywords = set([
            'dirscf','fdiff' ,'noconv','diis'  ,'soscf' ,'extrap','damp'  ,
            'shift' ,'rstrct','dem'   ,'cuhf'  ,'conv'  ,'sogtol','ethrsh',
            'maxdii','swdiis','locopt','demcut','dmpcut','uhfnos','vvos'  ,
            'mvoq'  ,'acavo' ,'pacavo','uhfchk','nhomo' ,'nlumo' ,'mom'   ,
            'kproj' ,'nco'   ,'nseto' ,'no'    ,'npair' ,'cicoef','couple',
            'f'     ,'alpha' ,'beta'  ,'npunch','npreo' ,'vtscal','scalf' ,
            'maxvt' ,'vtconv'
            ])
    integers = set([
            'maxdii','mvoq'  ,'nhomo'  ,'nlumo' ,'kproj','nco','nseto',
            'npair' ,'npunch','maxvt'
            ])
    reals    = set([
            'conv'  ,'sogtol','ethrsh' ,'swdiis','demcut','dmpcut',
            'scalf' ,'vtconv'
            ])
    bools    = set([
            'dirscf','fdiff' ,'noconv' ,'diis'  ,'soscf' ,'extrap',
            'damp'  ,'shift' ,'rstrct' ,'dem'   ,'cuhf'  ,'locopt',
            'uhfnos','vvos'  ,'acavo'  ,'uhfchk','mom'   ,'couple',
            'vtscal'
            ])
    arrays   = set([
            'pacavo','no'    ,'cicoef','f'     ,'alpha' ,'beta'  ,
            'npreo'
            ])

    allowed_values = obj(
        kproj = set([0,1,2]),
        )
#end class ScfGroup



class McscfGroup(KeywordSpecGroup):
    keywords = set([
            'cistep','focas' ,'soscf' ,'fullnr','quad'  ,'jacobi','acurcy',
            'engtol','maxit' ,'micit' ,'nword' ,'fors'  ,'canonc','finci' ,
            'diabat','ekt'   ,'npunch','npflg' ,'nofo'  ,'mcfmo' ,'casdii',
            'cashft','nrmcas','qudthr','damp'  ,'method','linser','fcore' ,
            'mofrz' ,'norb'  ,'norot' ,'dropc'
            ])
    integers = set(['maxit','micit','nword','npunch','nofo','mcfmo','nrmcas','norb'])
    reals    = set(['acurcy','engtol','casdii','cashft','qudthr','damp'])
    bools    = set(['focas','soscf','fullnr','quad','jacobi','fors','canonc',
                    'diabat','ekt','linser','fcore','dropc'])
    strings  = set(['cistep','finci','method'])
    arrays   = set(['npflg','mofrz','norot'])

    allowed_values = obj(
        cistep = set(['aldet','ormas','guga','genci','gmcci']),
        finci  = set(['none','mos','nos']),
        nrmcas = set([0,1]),
        method = set(['dm2','tei']),
        )
#end class McscfGroup



class DftGroup(KeywordSpecGroup):
    keywords = set([
            'method','dc'    ,'idcver','dcchg' ,'dcabc' ,'dcalp' ,'dcsr'  ,
            'dcs6'  ,'dcs8'  ,'lrdflg','mltint','lambda','kappa' ,'rzero' ,
            'prpol' ,'prcoef','prpair','lc'    ,'mu'    ,'chf'   ,'cmp2'  ,
            'nrad'  ,'nleb'  ,'sg1'   ,'jans'  ,'nrad'  ,'nthe'  ,'nphi'  ,
            'swoff' ,'switch','nrad0' ,'nleb0' ,'nthe0'  ,'nphi0' ,'thresh',
            'gthre' ,'auxfun','three'
            ])
    integers = set(['idcver','prcoef','prpair','nrad','nleb','jans','nthe',
                    'nphi','nrad0','nleb0','nthe0','nphi0','gthre'])
    reals    = set(['dcalp','dcsr','dcs6','dcs8','lambda','kappa','rzero',
                    'mu','chf','cmp2','swoff','switch','thresh'])
    bools    = set(['dc','dcchg','dcabc','lrdflg','mltint','prpol','lc','sg1',
                    'three'])
    strings  = set(['method','auxfun'])

    allowed_values = obj(
        method = set(['grid','gridfree']),
        idcver = set([1,2,3]),
        jans   = set([1,2]),
        auxfun = set(['aux0','aux3']),
        )
#end class DftGroup



class GugdiaGroup(KeywordSpecGroup):
    keywords = set([
            'nstate','prttol','mxxpan','itermx','cvgtol' ,'nword' ,'maxham',
            'maxdia','nimprv','nselct','selthr','nextra','kprint','nref','eref'
            ])

    integers = set(['nstate','mxxpan','itermx','nword','maxham','maxdia',
                    'nimprv','nselct','nextra','nref'])
    reals    = set(['prttol','cvgtol','selthr','eref'])
    arrays   = set(['kprint'])
#end class GugdiaGroup



class DrtGroup(KeywordSpecGroup):
    keywords = set([
            'group','fors'  ,'foci'  ,'soci','iexcit','intact','nmcc',
            'ndoc' ,'naos'  ,'nbos'  ,'nalp','nval'  ,'next'  ,'nfzv','stsym',
            'noirr','mxnint','mxneme','nprt'
            ])

    integers = set(['iexcit','nmcc','ndoc','naos','nbos','nalp','nval',
                    'next','nfzv','noirr','mxnint','mxneme','nprt'])
    bools    = set(['fors','foci','soci','intact'])
    strings  = set(['group','stsym'])

    allowed_values = obj(
        group = set(['c1','c2','ci','cs','c2v','c2h','d2','d2h','c4v','d4','d4h']),
        stsym = set(['a','ag','au','ap','app','a','b','a1','a2','b1','b2','ag',
                     'bu','bg','au','a','b1','b2','b3','ag','b1g','b2g','b3g',
                     'au','b1u','b2u','b3u']),
        nprt = set([0,1,2,3]),
        )
#end class DrtGroup



class CidrtGroup(KeywordSpecGroup):
    keywords = set([
            'group','fors'  ,'foci'  ,'soci','iexcit','intact','nfzc' ,
            'ndoc' ,'naos'  ,'nbos'  ,'nalp','nval'  ,'next'  ,'nfzv' ,'stsym',
            'noirr','mxnint','mxneme','nprt'
            ])

    integers = set(['iexcit','nfzc','ndoc','naos','nbos','nalp','nval',
                    'next','nfzv','noirr','mxnint','mxneme','nprt'])
    bools    = set(['fors','foci','soci','intact'])
    strings  = set(['group','stsym'])

    allowed_values = obj(
        group = set(['c1','c2','ci','cs','c2v','c2h','d2','d2h','c4v','d4','d4h']),
        stsym = set(['a','ag','au','ap','app','a','b','a1','a2','b1','b2','ag',
                     'bu','bg','au','a','b1','b2','b3','ag','b1g','b2g','b3g',
                     'au','b1u','b2u','b3u']),
        nprt = set([0,1,2,3]),
        )
#end class CidrtGroup
 


class DetGroup(KeywordSpecGroup):
    keywords = set([
            'ncore' ,'nact'  ,'nels'  ,'sz'    ,'group' ,'stsym' ,'irreps',
            'nstate','prttol','analys','itermx','cvgtol','nhgss' ,'nstgss',
            'mxxpan','clobbr','pures' ,'iroot' ,'nflgdm','saflg' ,'wstate',
            'idwref','dwparm'
            ])

    integers = set(['ncore','nact','nels','nstate','itermx','nhgss','nstgss',
                    'mxxpan','iroot','idwref'])
    reals    = set(['sz','prttol','cvgtol','dwparm'])
    bools    = set(['analys','clobbr','pures','saflg'])
    strings  = set(['group','stsym'])
    arrays   = set(['irreps','nflgdm','wstate'])

    allowed_values = obj(
        group = set(['c1','c2','ci','cs','c2v','c2h','d2','d2h','c4v','d4','d4h']),
        stsym = set(['a','ag','au','ap','app','a','b','a1','a2','b1','b2','ag',
                     'bu','bg','au','a','b1','b2','b3','ag','b1g','b2g','b3g',
                     'au','b1u','b2u','b3u']),
        )
#end class DetGroup



class BasisGroup(KeywordSpecGroup):
    keywords = set([
            'gbasis','ngauss','ndfunc','npfunc','diffsp','diffs',
            'polar' ,'split2','split3','basnam','extfil'
            ])

    integers = set(['ngauss','ndfunc','nffunc'])
    bools    = set(['diffsp','diffs','extfil'])
    strings  = set(['gbasis','polar'])
    arrays   = set(['split2','split3','basname'])

    allowed_values = obj(
        #gbasis = set(['sto','n21','n31','n311','g3l','g3lx','mini','midi','dzv',
        #              'dh','tzv','mc']) # many others
        ndfunc = set([0,1,2,3]),
        nffunc = set([0,1]),
        polar  = set(['common','popn31','popn311','dunning','huzinaga','hondo7']),
        )
#end class BasisGroup



#class XGroup(KeywordSpecGroup):
#    keywords = set([''])
#    integers = set([''])
#    reals    = set([''])
#    bools    = set([''])
#    strings  = set([''])
#    arrays   = set([''])
#    allowed_values = obj(
#         = set([]),
#        )
##end class XGroup






class GamessInput(SimulationInput,GIbase):
    group_order = '''
        contrl   system   basis    ecp      data     zmat     libe     
        scf      scfmi    dft      tddft    cis      cisvec   mp2      
        rimp2    auxbas   ccinp    eominp   mopac    guess    vec      
        mofrz    statpt   trudge   trurst   force    cphf     cpmchf   
        mass     hess     grad     dipdr    vib      vib2     vscf     
        vibscf   gamma    eqgeom   hlowt    glowt    irc      drc      
        mex      conicl   md       rdf      globop   gradex   surf     
        local    truncn   elmom    elpot    eldens   elfldg   points   
        grid     pdc      mgc      radial   molgrf   stone    raman    
        alpdr    comp     nmr      morokm   lmoeda   qmefp    ffcalc   
        tdhf     tdhfx    efrag    fragname frgrpl   ewald    makefp   
        prtefp   damp     dampgs   pcm      pcmgrd   mcpcav   tescav   
        newcav   iefpcm   pcmitr   disbs    disrep   svp      svpirf   
        cosgms   scrf     mcp      relwfn   efield   intgrl   fmm      
        trans    fmo      fmoprp   fmoxyz   optfmo   fmohyb   fmobnd   
        fmoenm   fmoend   optrst   gddi     elg      dandc    dccorr   
        subscf   subcor   mp2res   ccres    ciminp   cimatm   cimfrg   
        ffdata   ffpdb    ciinp    det      cidet    gen      cigen    
        ormas    ceeis    cedata   gcilst   gmcpt    pdet     adddet   
        remdet   sodet    drt      cidrt    mcscf    mrmp     detpt    
        mcqdpt   excorr   casci    ivoorb   cisort   gugem    gugdia   
        gugdm    gugdm2   lagran   trfdm2   diabat   transt
        drt1     drt2     vec1     vec2     det1     det2     hess2
        '''.split()

    all_groups = set(group_order)

    key_groups  = set(['contrl','system','guess','scf','mcscf','dft',
                       'gugdia','drt','cidrt','det','basis'])

    card_groups = set()
    #card_groups = set(['ecp','data','mcp','gcilst','points','stone','efrag',
    #                   'fragname','frgrpl','dampgs'])#,'fmoxyz'])

    formatted_groups = set()


    # detailed specifications for certain groups
    keyspec_groups = obj(
        contrl = ContrlGroup,
        system = SystemGroup,
        guess  = GuessGroup,
        scf    = ScfGroup,
        mcscf  = McscfGroup,
        dft    = DftGroup,
        gugdia = GugdiaGroup,
        drt    = DrtGroup,
        cidrt  = CidrtGroup,
        det    = DetGroup,
        basis  = BasisGroup
        )
    keyspec_group_order = []
    for gname in group_order:
        if gname in keyspec_groups:
            keyspec_group_order.append(gname)
        #end if
    #end for
    all_keywords = set()
    for g in keyspec_groups:
        all_keywords |= g.keywords
    #end for
    group_keyword_overlap = all_groups & all_keywords
    all_names = all_groups | all_keywords
    
    #cardspec_groups = obj()

    # aliases for generate_gamess_input
    group_aliases = obj()
    for gname in group_order:
        group_aliases['gamess_'+gname]=gname
    #end for
    all_group_aliases = all_groups | set(group_aliases.keys())
    all_name_aliases = all_group_aliases | all_keywords


    # gamess file I/O
    file_units = obj(
        #MCPPATH = -5,BASPATH = -4,EXTCAB  = -3,
        #MAKEFP  =  1, ERICFMT =  2, EXTBAS  =  3, 
        TRAJECT =  4, INPUT   =  5, 
        OUTPUT  =  6, PUNCH   =  7, AOINTS  =  8, MOINTS  =  9, DICTNRY = 10, 
        DRTFILE = 11, CIVECTR = 12, CASINTS = 13, CIINTS  = 14, WORK15  = 15, 
        WORK16  = 16, CSFSAVE = 17, FOCKDER = 18, WORK19  = 19, DASORT  = 20, 
        DFTINTS = 21, DFTGRID = 22, JKFILE  = 23, ORDINT  = 24, EFPIND  = 25, 
        PCMDATA = 26, PCMINTS = 27, MLTPL   = 28, MLTPLT  = 29, DAFL30  = 30, 
        RESTART = 35, HESSIAN = 38, SOCCDAT = 40, AABB41  = 41, BBAA42  = 42, 
        BBBB43  = 43, REMD    = 44, MCQD50  = 50, MCQD51  = 51, MCQD52  = 52, 
        MCQD53  = 53, MCQD54  = 54, MCQD55  = 55, MCQD56  = 56, MCQD57  = 57, 
        MCQD58  = 58, MCQD59  = 59, MCQD60  = 60, MCQD61  = 61, MCQD62  = 62, 
        MCQD63  = 63, MCQD64  = 64, DCPHFH2 = 67, NMRINT1 = 61, CCREST  = 70, 
        CCDIIS  = 71, CCINTS  = 72, CCT1AMP = 73, CCT2AMP = 74, CCT3AMP = 75, 
        CCVM    = 76, CCVE    = 77, CCQUADS = 78, QUADSVO = 79, EOMSTAR = 80, 
        EOMVEC1 = 81, EOMVEC2 = 82, EOMHC1  = 83, EOMHC2  = 84, EOMHHHH = 85, 
        EOMPPPP = 86, EOMRAMP = 87, EOMRTMP = 88, EOMDG12 = 89, MMPP    = 90, 
        MMHPP   = 91, MMCIVEC = 92, MMCIVC1 = 93, MMCIITR = 94, EOMVL1  = 95, 
        EOMVL2  = 96, EOMLVEC = 97, EOMHL1  = 98, EOMHL2  = 99, EFMOI   = 102, 
        EFMOF   = 103 
        )

    def __init__(self,filepath=None):
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def read_text(self,contents,filepath=None):
        #print 8*'\n'
        #print contents

        groups = obj()
        lines = contents.splitlines()
        ingroup = False
        incard  = False
        group_name = None
        group_text = ''
        gname = ''
        gtext = ''
        n=0
        for line in lines:
            ended = False
            ls = line.strip()
            # specialized parsing for unknown card groups
            if ingroup and ls!='$END' and ls!='$end':
                gtext+=line+'\n'
            #end if
            if incard:
                ended = ls=='$END' or ls=='$end'
                ingroup = not ended
                incard  = not ended
                if ended:
                    groups[group_name] = group_text
                    group_name = None
                    group_text = ''
                else:
                    group_text+=line+'\n'
                #end if
            elif len(line)>0 and line[0]==' ' and ls!='':
                if len(line)>1 and line[1]=='$' and not ingroup:
                    if not ' ' in ls:
                        group_name = ls.replace('$','').lower()
                        gname = group_name
                        ingroup = True
                    else:
                        group_name,ls = ls.split(' ',1)
                        group_name = group_name.replace('$','').lower()
                        gname = group_name
                        text,ended = self.process_line(ls)
                        group_text += text
                        ingroup = not ended
                        if ended:
                            groups[group_name] = group_text
                            group_name = None
                            group_text = ''
                        #end if
                    #end if
                    incard = group_name in self.card_groups
                elif ingroup:
                    text,ended = self.process_line(ls)
                    group_text += text
                    ingroup = not ended
                    if ended:
                        groups[group_name] = group_text
                        group_name = None
                        group_text = ''
                    #end if
                elif not ingroup:
                    None
                else:
                    self.error('invalid text encountered during read of line number {0}:\n{1}'.format(n,line))
                #end if
            elif ls=='' or line[0]!=' ' or not ingroup:
                None
            else:
                self.error('invalid text encountered during read of line number {0}:\n{1}'.format(n,line))
            #end if                    
            # specialized parsing for unknown card groups
            if ended:
                if not '=' in groups[gname]:
                    groups[gname]=gtext
                #end if
                gtext = ''
                gname = ''
            #end if
        #end for

        #print groups

        for group_name,group_text in groups.iteritems():
            failed = False
            if group_name in self.keyspec_groups:
                self[group_name] = self.keyspec_groups[group_name](group_text)
            #elif group_name in self.cardspec_groups:
            #    self[group_name] = self.cardspec_groups[group_name](group_text)
            elif group_name in self.key_groups:
                self[group_name] = KeywordGroup(group_text)
            elif group_name in self.card_groups:
                self[group_name] = CardGroup(group_text)
            elif '=' in group_text:
                try:
                    self[group_name] = KeywordGroup(group_text)
                except:
                    try:
                        self[group_name] = FormattedGroup(group_text)
                    except:
                        failed = True
                    #end try
                #end try
            else:
                try:
                    self[group_name] = FormattedGroup(group_text)
                except:
                    failed = True
                #end try
            #end if
            if failed:
                self.message('Read failure: group "{0}" does not appear to be a keyword group\nand a generic read of card data failed\ndata for this group will not be available'.format(group_name))
            #end if
        #end for
                    
        #print self

        #exit()
    #end def read_text

        
    def process_line(self,ls):
        ended = True
        if ls.endswith('$END'):
            text = ls.replace('$END','')
        elif ls.endswith('$end'):
            text = ls.replace('$end','')
        else:
            text = ls
            ended = False
        #end if
        cloc = text.find('!')
        if cloc!=-1:
            text = text[0:cloc]
        #end if
        text +='\n'
        return text,ended
    #end def process_line


    def write_text(self,filepath=None):
        contents = ''
        extra_groups = set(self.keys())-set(self.group_order)
        if len(extra_groups)>0:
            self.error('write failed\nthe following groups are unknown: {0}'.format(sorted(extra_groups)))
        #end if
        for group in self.group_order:
            if group in self and isinstance(self[group],KeywordGroup):
                contents += self[group].write(group)
            #end if
        #end for
        for group in self.group_order:
            if group in self and isinstance(self[group],(CardGroup,FormattedGroup)):
                contents += self[group].write(group)
            #end if
        #end for
        return contents
    #end def write_text


    def incorporate_system(self,system):
        self.not_implemented()
    #end def incorporate_system
#end class GamessInput





def generate_gamess_input(**kwargs):
    if 'input_type' in kwargs:
        input_type = kwargs['input_type']
        del kwargs['input_type']
    else:
        input_type = 'general'
    #end if
    if input_type=='general':
        gi = generate_any_gamess_input(**kwargs)
    else:
        GamessInput.class_error('input_type {0} is unrecognized\nvalid options are: general'.format(input_type))
    #end if
    return gi
#end def generate_gamess_input



ps_inputs = set('descriptor symmetry system pseudos pseudo_bases bases'.split())
ps_defaults = obj()
for var in ps_inputs:
    ps_defaults[var]=None
#end for
ps_defaults.set(
    descriptor = 'A molecule.',
    symmetry   = 'C1'
    )
kw_defaults = obj()
for var in GamessInput.all_keywords:
    kw_defaults[var]=None
#end for


def generate_any_gamess_input(**kwargs):
    kwset = set(kwargs.keys())
    pskw = ps_defaults.copy()
    ps_overlap = ps_inputs & kwset
    if len(ps_overlap)>0:
        pskw.move_from(kwargs,ps_overlap)
        kwset = set(kwargs.keys())
    #end if
    for name in kwargs.keys():
        val = kwargs[name]
        if isinstance(val,dict):
            kwargs[name] = GIarray(val)
        #end if
    #end for
    kw = kw_defaults.copy()
    kw.set(**kwargs)
    kwrem = obj(**kwargs)

    invalid_names = kwset-GamessInput.all_name_aliases
    if len(invalid_names)>0:
        GamessInput.class_error('invalid group names or keywords encountered\ninvalid names/keywords provided: {0}\nplease check if these group names or keywords are actually valid GAMESS inputs\nif so, unsupported groups can be generated by providing the keywords as a single argument:\n  generate_gamess_input(\n    ...,\n    group_name = obj(assign keywords),\n    ...,\n    )'.format(sorted(invalid_names)),'generate_gamess_input')
    #end if

    gi = GamessInput()
    
    # handle groups provided directly by the user
    #   use aliases to guard against namespace collisions w/ nexus (e.g. system)
    group_names = kwset & GamessInput.all_group_aliases
    for name in group_names:
        group_info = kw[name]
        vname = name
        if name in GamessInput.group_aliases:
            name = GamessInput.group_aliases[name]
        #end if
        if isinstance(group_info,obj):
            for n in group_info.keys():
                v = group_info[n]
                if isinstance(v,dict):
                    group_info[n] = GIarray(v)
                #end if
            #end for
            if isinstance(group_info,Group):
                gi[name] = group_info
            elif name in GamessInput.keyspec_groups:
                gi[name] = GamessInput.keyspec_groups[name](**group_info)
            #elif name in GamessInput.cardspec_groups:
            #    gi[name] = GamessInput.cardspec_groups[name](**group_info)
            elif name in GamessInput.key_groups:
                gi[name] = KeywordGroup(**group_info)
            elif name in GamessInput.card_groups:
                GamessInput.class_error('card group {0} cannot be generated from a keyword list\n  keyword list provided:\n  {1}'.format(name,group_info),'generate_gamess_input')
            elif name in GamessInput.formatted_groups:
                GamessInput.class_error('formatted group {0} cannot be generated from a keyword list\n  keyword list provided:\n  {1}'.format(name,group_info),'generate_gamess_input')
            else:
                gi[name] = KeywordGroup(**group_info) # assume keyword group
            #end if
            del kw[vname]
            del kwrem[vname]
        elif name in GamessInput.group_keyword_overlap:
            None
        else:
            GamessInput.class_error('invalid information provided to initialize group {0}\nyou must provide a dict, obj, or Group\nyou provided {1}'.format(vname,group_info),'generate_gamess_input')
        #end if
    #end for

    # load keywords into groups by group order
    #   this may not be correct for overlapping keywords between groups!
    #   user will have to supply explicit keyword subsets by group in obj's as above
    for name in GamessInput.keyspec_group_order:
        group_type = GamessInput.keyspec_groups[name]
        keywords = group_type.keywords & set(kwrem.keys())
        if len(keywords)>0:
            group_info = obj()
            group_info.move_from(kwrem,keywords)
            gi[name] = group_type(**group_info)
        #end if
    #end for
    if len(kwrem)>0:
        GamessInput.class_error('encountered unrecognized keywords\n  unrecognized keywords: {0}\n  these keywords may belong to groups not fully implemented here\n  fully supported groups: {1}\n  unsupported groups can be generated by providing the keywords as a single argument: group_name = obj(assign keywords)'.format(sorted(kwrem),GamessInput.keyspec_group_order))
    #end if

    # handle nexus specific input generation keywords
    #  ecp 287
    #  data 37
    if pskw.system!=None and not 'data' in gi:
        system = pskw.system
        if not 'contrl' in gi:
            gi.contrl = ContrlGroup()
        #end if
        gi.contrl.set(
            icharg = system.net_charge,
            mult   = system.net_spin+1
            )
        elem = system.structure.elem
        pos  = system.structure.pos
        pskw.symmetry = pskw.symmetry.strip()
        data = '{0}\n{1}\n'.format(pskw.descriptor,pskw.symmetry)
        if pskw.symmetry!='C1':
            data+='\n'
        #end if
        if pskw.pseudos is None:
            if pskw.bases!=None:
                bss = nexus_noncore.basissets.bases_by_atom(*pskw.bases)
            else:
                bss = obj()
                if 'coord' not in gi.contrl:
                    gi.contrl.coord = 'unique'
                #end if
            #end if
            for i in range(len(elem)):
                a = elem[i]
                Z = pt[a].atomic_number
                data+='{0} {1:3.2f} {2:16.8f} {3:16.8f} {4:16.8f}\n'.format(a,Z,*pos[i])
                if a in bss:
                    data+=bss[a].text+'\n\n'
                #end if
            #end for
        else:
            gi.contrl.set(
                coord = 'unique',
                ecp   = 'read'
                )
            ecp = ''
            pps = nexus_noncore.pseudopotentials.pseudos_by_atom(*pskw.pseudos)
            atoms = set()
            for i in range(len(elem)):
                a = elem[i]
                Z = pt[a].atomic_number
                data+='{0} {1} {2:16.8f} {3:16.8f} {4:16.8f}\n'.format(a,Z,*pos[i])
                if a in pps:
                    pp = pps[a]
                    data += pp.basis_text+'\n\n'
                    if a in atoms:
                        ecp += pp.pp_name+'\n'
                    else:
                        ecp += pp.pp_text+'\n'
                    #end if
                #end if
                atoms.add(a)
            #end for
            gi.ecp = FormattedGroup(ecp)
        #end if
        gi.data = FormattedGroup(data)
    #end if

    return gi
#end def generate_any_gamess_input









def check_keyspec_groups():
    print 'checking GamessInput KeywordSpecGroups'    
    groups      = GamessInput.keyspec_groups
    group_order = GamessInput.group_order
    glist = []
    for group_name in group_order:
        if group_name in groups:
            glist.append(group_name)
        #end if
    #end for

    failed = False

    #check for unrecognized groups
    extra_groups = set(groups.keys())-set(group_order)
    if len(extra_groups)>0:
        failed = True
        print '  encountered unrecognized keyspec groups: {0}'.format(sorted(extra_groups))
    #end if

    #check that integers, reals, bools, strings, and arrays are non-overlapping subsets of keywords
    #check that allowed_values are a subset of keywords and values specified are of the correct type
    for group_name in glist:
        g = groups[group_name]
        go = obj(
            integers = g.integers,
            reals    = g.reals,
            bools    = g.bools,
            strings  = g.strings,
            arrays   = g.arrays
            )
        overlaps = obj()
        for tname1,tset1 in go.iteritems():
            for tname2,tset2 in go.iteritems():
                if tname1!=tname2:
                    overlap = tset1 & tset2
                    if len(overlap)>0:
                        overlaps[tname1,tname2] = sorted(overlap)
                    #end if
                #end if
            #end for
        #end for
        if len(overlaps)>0:
            failed = True
            msg = '  keyspec group {0} has overlapping keywords'.format(g.__name__)
            for tname1,tname2 in sorted(overlaps.keys()):
                msg += '    \n {0} {1} overlap: {2}'.format(tname1,tname2,overlaps[tname1,tname2])
            #end for
            print msg
        #end if
        for tname in sorted(go.keys()):
            extra_keys = go[tname]-g.keywords
            if len(extra_keys)>0:
                failed = True
                print '  keyspec group {0} has unrecognized {1} keywords:\n    {2}'.format(g.__name__,tname,sorted(extra_keys))
            #end if
        #end for
        extra_keys = set(g.allowed_values.keys())-g.keywords
        if len(extra_keys)>0:
            failed = True
            print '  keyspec group {0} has unrecognized allowed_value keywords:\n    {1}'.format(g.__name__,sorted(extra_keys))
        #end if
        type_keys = set()
        for keys in go:
            type_keys |= keys
        #end for
        undefined = g.keywords-type_keys
        if len(undefined)>0:
            print '  keyspec group {0} has keywords w/o type assignment:\n    {1}'.format(g.__name__,sorted(undefined))
        #end if

        #check that allowed values for each keyword have the right type
        to = obj(
            integers = int,
            reals    = float,
            bools    = bool,
            strings  = str,
            arrays   = ndarray
            )
        for tname in sorted(go.keys()):
            type = to[tname]
            for kw in sorted(go[tname]):
                if kw in g.allowed_values:
                    for val in g.allowed_values[kw]:
                        if not isinstance(val,type):
                            failed = True
                            print '  allowed values of {0} keyword {1} are not all {2}: {3}'.format(g.__name__,kw,tname,sorted(g.allowed_values[kw]))
                            break
                        #end if
                    #end for
                #end if
            #end for
        #end for
    #end for

    #note any overlapping keywords between groups (this is a feature, not an error)
    overlaps = obj()
    for gname1 in glist:
        kw1 = groups[gname1].keywords
        for gname2 in glist:
            kw2 = groups[gname2].keywords
            if gname1!=gname2:
                overlap = kw1 & kw2
                if len(overlap)>0:
                    tup = tuple(sorted((gname1,gname2)))
                    overlaps[tup] = sorted(overlap)
                #end if
            #end if
        #end for
    #end for
    if len(overlaps)>0:
        print '\n  Note: some groups have overlapping keywords'
        for gname1,gname2 in sorted(overlaps.keys()):
            print '    groups {0} and {1} have overlapping keywords:\n      {2}'.format(gname1,gname2,overlaps[gname1,gname2])
        #end for
    #end if

    #note any overlapping keyword and group names (also a feature)
    overlap = GamessInput.all_keywords & set(GamessInput.group_order)
    if len(overlap)>0:
        print '\n  Note: some group names overlap with keywords:\n    {0}'.format(sorted(overlap))
    #end if

    if failed:
        print '\ncheck failed, see messages above and fix implementation'
    else:
        print '\ncheck passed'
    #end if
    exit()
#end def check_keyspec_groups

#check_keyspec_groups()  # uncomment this to check keyword spec group self-consistency





