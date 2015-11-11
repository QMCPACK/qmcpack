##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_pp_validation.py                                          #
#    Tool to test pseudopotentials with PWSCF and SQD.  Input format #
#    is complicated and the tool will be superseded in the future.   #
#    A usage example can be found at the end of the file.            #
#                                                                    #
#====================================================================#


#! /usr/bin/env python


import os
from copy import deepcopy
from numpy import array,empty,identity

from generic import obj
from project import generate_physical_system
from project import generate_pwscf
from project import generate_pw2casino,generate_pw2qmcpack
from project import generate_sqd
from project import generate_qmcpack,generate_jastrow
from project import loop,vmc,dmc,linear,cslinear
from project import settings,ProjectManager,Job,Settings

from qmcpack_calculations import basic_qmc,standard_qmc

from periodic_table import periodic_table
from unit_converter import convert
from developer import DevBase
from plotting import *
from debug import *




def remove(o,*names):
    m = o.copy()
    keys = list(m.keys())
    for key in keys:
        for name in names:
            if name in key:
                del m[key]
                break
            #end if
        #end for
    #end for
    return m
#end def remove



def sort_pseudos(pps,dext=set(['ncpp','upf']),qext=set(['xml'])):
    dpp = []
    qpp = []
    for pp in pps:
        n,ext = pp.lower().rsplit('.',1)
        if ext in dext:
            dpp.append(pp)
        elif ext in qext:
            qpp.append(pp)
        else:
            print 'Error: pseudopotential format could not be determined\n  filename: {0}\n  extensions checked: {1}'.format(pp,list(dext)+list(qext))
            exit()
        #end if
    #end for
    return dpp,qpp
#end def sort_pseudos




class ColorWheel(DevBase):
    def __init__(self):
        colors = 'Black Maroon DarkOrange Green DarkBlue Purple Gray Firebrick Orange MediumSeaGreen DodgerBlue MediumOrchid'.split()
        styles = '- -- -. :'.split()
        c = []
        for style in styles:
            for color in colors:
                c.append((color,style))
            #end for
        #end for
        self.colors = c
        self.icolor = -1
    #end def __init__


    def next(self):
        self.icolor = (self.icolor+1)%len(self.colors)
        return self.colors[self.icolor]
    #end def next


    def reset(self):
        self.icolor = -1
    #end def reset
#end class ColorWheel
color_wheel = ColorWheel()
    



class ValidatePPBase(DevBase):
    settings = None

    def performed_runs(self):
        s = self.settings
        return not s.generate_only and not s.status_only
    #end def performed_runs
#end class ValidatePPBase




class ValidationStage(ValidatePPBase):
    stage_inputs = set()
    stage_dependencies = set()
    stage_result = None
    final_result = None


    stage_type = None
    data_type  = None
    title      = None
    xlabel     = None
    ylabel     = None
    xunit      = None
    yunit      = None


    def __init__(self):
        self.set(
            name         = None,
            ready        = False,
            sims         = obj(),
            saved_inputs = obj(),
            simlist      = [],
            inputs       = obj(),
            results      = obj(),
            result_sims  = None
            )
    #end def __init__


    def get_tags(self):
        s=self
        return s.title,s.xlabel,s.ylabel,s.xunit,s.yunit,s.stage_type,s.data_type
    #end def get_tags


    def set_name(self,name):
        self.name = name
    #end def set_name


    def set_result_sims(self,result_sims):
        self.result_sims = result_sims
    #end def set_result_sims    


    def check_dependencies(self):
        inputs = self.inputs
        results = self.results
        have_inputs = True
        for inp in self.stage_inputs:
            have_inputs = have_inputs and inp in inputs and inputs[inp]!=None
        #end for
        have_deps = True
        for dep in self.stage_dependencies:
            have_deps = have_deps and dep in results and results[dep]!=None
        #end for
        self.ready = have_inputs and have_deps
    #end def check_dependencies


    def set_stage_inputs(self,**inputs):
        self.inputs.set(**inputs)
    #end def set_stage_inputs


    def set_stage_results(self,**results):
        self.results.set(**results)
    #end def set_stage_results


    def check_inputs(self,inputs):
        inkeys = set(inputs.keys())
        if not self.stage_inputs <= inkeys:
            self.error('stage inputs are missing\n  inputs required: {0}\n  inputs provided: {1}'.format(list(self.stage_inputs),list(inkeys)))
        #end if
        if not self.stage_dependencies <= inkeys:
            self.error('stage dependencies are missing\n  dependencies required: {0}\n  dependencies provided: {1}'.format(list(self.stage_dependencies),list(inkeys)))
        #end if
    #end def check_inputs

        
    def save_sims(self,sims,key,inputs,inputs_save=None):
        if inputs_save is None:
            inputs_save = inputs
        #end if
        if 'result' in sims:
            if self.stage_result is None:
                self.error('encountered result sim for a stage with no results')
            #end if
            self.result_sims[self.stage_result,key] = sims.result
            del sims.result
        #end if
        self.sims[key] = sims
        self.saved_inputs[key] = inputs_save
        block = False
        for dep in self.stage_dependencies:
            if inputs[dep] is None:
                block = True
            #end if
        #end for
        if not block:
            self.simlist.extend(sims.list())
        #end if
    #end def save_sims

            
    def get_result_sim(self,key):
        if not key in self.result_sims:
            self.error('attempted to retrieve non-existent stage result\n  key requested: {0}\n  keys present: {1}'.format(key,self.result_sims.keys()))
        #end if
        return self.result_sims[key]
    #end def get_result_sim


    def sim_list(self):
        return self.simlist
    #end def sim_list


    def status(self,pad=' ',n=0,*args):
        if self.ready:
            status = 'ready'
        else:
            status = 'waiting'
        #end if
        print n*pad+'{0:<20}  {1}'.format(self.name,status)
        if not self.ready or 'verbose' in args:
            p = (n+1)*pad
            presence = {True:'present',False:'absent'}
            for name in self.stage_inputs:
                print p+'{0:<20}  {1}'.format(name,presence[name in self.inputs])
            #end for
            for name in self.stage_dependencies:
                print p+'{0:<20}  {1}'.format(name,presence[name in self.results])
            #end for
        #end if        
    #end def status


    def make_all_sims(self,basepath):
        self.not_implemented()
    #end def make_all_sims


    def get_results(self):
        self.not_implemented()
    #end def get_results
#end class ValidationStage




class ValidationProcess(ValidatePPBase):

    allowed_stage_inputs  = set()
    allowed_stage_results = set()

    def __init__(self,**kwargs):
        self.name = None
        self.stages = obj()
        self.stage_order = None
        self.results = obj()
        si = obj()
        for name in self.allowed_stage_inputs:
            si[name] = None
        #end for
        self.stage_inputs = si
        sr = obj()
        for name in self.allowed_stage_results:
            sr[name] = None
        #end for
        self.stage_results = sr
        self.result_sims = obj()
        self.initialize(**kwargs)
        for stage in self.stages:
            stage.set_result_sims(self.result_sims)
        #end for    
    #end def __init__

        
    def initialize(self,**kwargs):
        self.not_implemented()
    #end def initialize


    def set_name(self,name):
        self.name = name
    #end def set_name

        
    def get_stage_order(self):
        if self.stage_order!=None:
            stage_order = self.stage_order
        else:
            stage_order = self.stages.keys()
            stage_order.sort()
        #end if
        return stage_order
    #end def get_stage_order


    def status(self,pad='  ',n=0,*args):
        print n*pad+self.name
        stage_order = self.get_stage_order()
        for sname in stage_order:
            self.stages[sname].status(pad,n+1,*args)
        #end for
        print n*pad+'end '+self.name
    #end def status


    known_functionals = set('lda pbe pw91 hse pbe0'.split())
    pp_extensions = [['upf','ncpp'],['xml']]

    def check_funcs_pseudos(self,**func_pp):
        ppmap = obj()
        ppdir = self.settings.pseudo_dir
        if not os.path.exists(ppdir):
            self.error('pseudopotential directory does not exist\n  directory provided: '+ppdir)
        #end if
        ppfiles = os.listdir(ppdir)
        ppfiles.sort()
        unknown_functionals = []
        missing_pp_files = []
        for functional,pseudos in func_pp.iteritems():
            if not functional in self.known_functionals:
                unknown_functionals.append(functional)
            else:
                for pseudo_key in pseudos:
                    if isinstance(pseudo_key,str):
                        pseudo_prefix = pseudo_key
                        pp_found,pp_missing = self.findpp(pseudo_prefix,ppfiles)
                    elif isinstance(pseudo_key,tuple):
                        pseudo_tuple = pseudo_key
                        pp_found = []
                        pp_missing = []
                        for pseudo_prefix in pseudo_tuple:
                            ppf,ppm = self.findpp(pseudo_prefix,ppfiles)
                            pp_found.extend(ppf)
                            pp_missing.extend(ppm)
                        #end for
                    #end if
                    ppmap[pseudo_key] = pp_found
                    #end for
                #end for
            #end if
        #end for
        msg=''
        if len(unknown_functionals)>0:
            msg+='inputted functionals are unrecognized\n  functionals: '+str(unknown_functionals)+'\n'
        #end if
        if len(missing_pp_files)>0:
            msg+='pseudopotential prefixes had no corresponding files\n  prefixes: '+str(missing_pp_files)+'\n  extensions required: '+str(self.pp_extensions)+'\n  directory searched: '+ppdir+'\n  directory contents: '+str(ppfiles)+'\n'
        #end if
        if msg!='':
            self.error(msg)
        #end if
        self.functionals = list(func_pp.keys())
        self.func_pp     = obj(**func_pp)
        self.ppmap       = ppmap
    #end def check_funcs_pseudos


    def findpp(self,pseudo_prefix,ppfiles):
        pp_found   = []
        pp_missing = []
        for pp_ext in self.pp_extensions:
            pp_prefix = pseudo_prefix.lower()
            for ppfile in ppfiles:
                prefix,ext = ppfile.rsplit('.',1)
                prefix = prefix.lower()
                ext = ext.lower()
                if prefix==pp_prefix and ext in pp_ext:
                    pp_found.append(ppfile)
                    break
                #end if
            #end for
        #end for
        if len(pp_found)!=len(self.pp_extensions):
            pp_missing.append(pseudo_prefix)
        #end if
        return pp_found,pp_missing
    #end def findpp

    def pseudo_list(self):
        return list(self.ppmap.keys())
    #end def pseudo_list


    def add_stage(self,name,stage):
        if not isinstance(stage,ValidationStage):
            self.error(name+' is not a ValidationStage')
        #end if
        self.stages[name] = stage
        stage.set_name(name)
    #end def add_stage


    def add_stages(self,**stages):
        for name,stage in stages.iteritems():
            self.add_stage(name,stage)
        #end for
    #end def add_stages


    def set_stage_order(self,stage_order):
        stages = set(self.stages.keys())
        order = set(stage_order)
        if stages!=order:
            self.error('stage order must account for all stages\n  stages: {0}  stage_order: {1}'.format(list(stages),stage_order))
        #end if
        self.stage_order = stage_order
    #end def set_stage_order


    def set_stage_inputs(self,**inputs):
        inkeys = set(inputs.keys())
        invalid = inkeys-self.allowed_stage_inputs
        if len(invalid)>0:
            self.error('invalid stage inputs encountered\n  invalid inputs: {0}\n valid options are: {1}'.format(invalid,self.allowed_stage_inputs))
        #end if
        self.stage_inputs.set(**inputs)
        for stage in self.stages:
            stage.set_stage_inputs(**inputs)
        #end for
    #end def set_stage_inputs


    def check_stage_results(self,inkeys):
        invalid = set(inkeys)-self.allowed_stage_results
        if len(invalid)>0:
            self.error('invalid stage results encountered\n  invalid results: {0}\n valid options are: {1}'.format(invalid,self.allowed_stage_results))
        #end if
    #end def check_stage_results


    def set_stage_results(self,**results):
        self.stage_results.set(**results)
        for stage in self.stages:
            stage.set_stage_results(**results)
        #end for
    #end def set_stage_results


    def check_dependencies(self):
        for stage in self.stages:
            stage.check_dependencies()
        #end for
    #end def check_dependencies


    def pre_make_sims(self):
        None
    #end def pre_make_sims


    def make_sims(self,basepath):
        print 'make_sims',self.__class__.__name__
        self.pre_make_sims()
        self.check_dependencies()
        stage_order = self.get_stage_order()
        for name in stage_order:
            stage = self.stages[name]
            path = os.path.join(basepath,name)
            print ' ',name,stage.__class__.__name__,stage.ready
            if stage.ready:
                stage.make_all_sims(path)
            else:
                stage.status('  ',3,'verbose')
            #end if
        #end for
    #end def make_sims


    def sim_list(self):
        sims = []
        for stage in self.stages:
            sims.extend(stage.sim_list())
        #end for
        return sims
    #end def sim_list
#end class ValidationProcess



class ValidationProcesses(ValidatePPBase):

    process_type = None

    @classmethod
    def set_process_type(cls,ptype):
        cls.process_type = ptype
        cls.allowed_stage_inputs  = ptype.allowed_stage_inputs
        cls.allowed_stage_results = ptype.allowed_stage_results
    #end def set_process_type


    def __init__(self,**processes):
        self.name = None
        self.processes = obj()
        for pname,pinfo in processes.iteritems():
            pinfo['name'] = pname
            process = self.process_type(**pinfo)
            #process.set_name(pname)
            self.add_process(pname,process)
        #end for
        self.initialize()
    #end def __init__


    def set_name(self,name):
        self.name = name
    #end def set_name


    def status(self,pad='  ',n=0,*args):
        print n*pad+self.name
        pnames = self.processes.keys()
        pnames.sort()
        for pname in pnames:
            self.processes[pname].status(pad,n+1,*args)
        #end for
        print n*pad+'end '+self.name
    #end def status


    def initialize(self):
        None
    #end def initialize


    def add_process(self,name,process):
        if not isinstance(process,ValidationProcess):
            self.error('Validation input improperly constructed\n  {0} is not a ValidationProcess'.format(name))
        #end if
        self[name] = process
        self.processes[name] = process
    #end def add_process

        
    def sim_list(self):
        sims = []
        for process in self.processes:
            sims.extend(process.sim_list())
        #end for
        return sims
    #end def sim_list

    
    def pre_make_sims(self):
        None
    #end def pre_make_sims


    def make_sims(self,basepath):
        print 'make_sims',self.__class__.__name__
        self.pre_make_sims()
        for name,process in self.processes.iteritems():
            path = os.path.join(basepath,name)
            process.make_sims(path)
        #end for
    #end def make_sims


    def show_results(self,
                     stages     = None,
                     quantities = None,
                     processes  = None,
                     mode       = 'plot',
                     titles     = None,
                     xlabels    = None,
                     ylabels    = None,
                     xlims      = None,
                     ylims      = None,
                     lw         = 2,
                     units      = None,
                     xunits     = None,
                     figlevel   = 'quantity',
                     suppress   = False):
        if not self.performed_runs():
            return
        #end if
        modes = ['plot','print','all']
        figlevels = [None,'top','stage','quantity','process']    
        if not mode in modes:
            self.error('invalid mode encountered\n  you provided: {0}\n  valid options are: {1}'.format(mode,modes))
        #end if    
        if not figlevel in figlevels:
            self.error('invalid figlevel encountered\n  you provided: {0}\n  valid options are: {1}'.format(figlevel,figlevels))
        #end if    
        if stages is None:
            self.error('results cannot be shown unless the stages variable is specified')
        #end if    
        if isinstance(stages,str):
            stages = [stages]
        #end if    
        if quantities is None or isinstance(quantities,str):
            quantities = [quantities]
        #end if    
        if processes is None:
            processes = self.processes.keys()
        #end if
        plotting = mode=='plot'  or mode=='all'
        printing = mode=='print' or mode=='all'
        color_wheel.reset()
        if printing:
            print
            print '{0} results'.format(self.name)
        #end if
        if figlevel=='top' and plotting:
            figure()
        #end if    
        for sname in stages:
            if printing:
                print '  {0} stage results'.format(sname)
            #end if
            if figlevel=='stage' and plotting:
                figure()
            #end if    
            for qname in quantities:
                if qname!=None and printing:
                    print '    quantity {0} results'.format(qname)
                    pad = '      '
                else:
                    pad = '    '
                #end if
                if figlevel=='quantity' and plotting:
                    figure()
                #end if    
                for pname in processes:
                    if printing:
                        print pad+pname
                    #end if
                    if figlevel=='process' and plotting:
                        figure()
                    #end if    
                    if not pname in self.processes:
                        self.error('process {0} is unrecognized\n  valid options are: {1}'.format(pname,self.processes.keys()))
                    #end if
                    process = self.processes[pname]
                    if not sname in process.stages:
                        self.error('stage {0} is not part of process {1}\n  valid options are: {2}'.format(sname,process.name,process.stages.keys()))
                    #end if    
                    stage = process.stages[sname] 

                    res = stage.get_results(quantity=qname)
                    rkeys = res.keys()
                    rkeys.sort()

                    titlet,xlabelt,ylabelt,xunitt,yunitt,stype,dtype = stage.get_tags()
                    if titles!=None:
                        ptitle = titles
                    else:
                        ptitle = titlet
                    #end if
                    if xlabels!=None:
                        pxlabel = xlabels
                    else:
                        pxlabel = xlabelt
                    #end if
                    if ylabels!=None:
                        pylabel = ylabels
                    else:
                        pylabel = ylabelt
                    #end if

                    rtypes = {
                        ('value','scalar'):'y',
                        ('value','stat'  ):'ye',
                        ('scan' ,'scalar'):'xy',
                        ('scan' ,'stat'  ):'xye',
                        }
                    rtype = rtypes[stype,dtype]
                    convertx = xunits!=None and xunitt!=None
                    converty = units!=None and yunitt!=None
                    xunitb = ''
                    yunitb = ''
                    if convertx:
                        xunitb = '('+xunits+')'
                    elif xunitt!=None:
                        xunitb = '('+xunitt+')'
                    #end if
                    if converty:
                        yunitb = '('+units+')'
                    elif yunitt!=None:
                        yunitb = '('+yunitt+')'
                    #end if
                    if rtype=='y':
                        phead = '{0} {1:<4}'.format(ylabelt,yunitb)
                        pfmt  = '{0}     '
                        plotter = plot
                    elif rtype=='ye':
                        phead = '{0} {1:<4}  {2} {3:<4}'.format(ylabelt,yunitb,'error',yunitb)
                        pfmt  = '{0}       {1}'
                        plotter = errorbar
                    elif rtype=='xy':
                        phead = '{0} {1:<4}  {2} {3:<4}'.format(xlabelt,xunitb,ylabelt,yunitb)
                        pfmt  = '{0}       {1}'
                        plotter = plot
                    elif rtype=='xye':
                        phead = '{0} {1:<4}  {2} {3:<4}  {4} {5:<4}'.format(xlabelt,xunitb,ylabelt,yunitb,'error',yunitb)
                        pfmt  = '{0}       {1}       {2}'
                        plotter = errorbar
                    #end if
                    keylen = 0
                    for key in rkeys:
                        keylen = max(keylen,len(str(key)))
                    #end for
                    keylen+=4

                    if printing:
                        p = '  '
                        print pad+p+titlet
                        print pad+p+keylen*' '+phead
                    #end if
                    if stype=='value':
                        color,style = color_wheel.next()
                    #end if
                    for key in rkeys:
                        r = res[key]
                        x,y,yerr = None,None,None
                        xlen = 0
                        if stype=='value':
                            if isinstance(key,tuple):
                                x = key[-1] # hope this works
                            else:
                                x = key
                            #end if
                            if dtype=='scalar':
                                y    = r
                                yerr = 0
                            elif dtype=='stat':
                                y,yerr = r
                            #end if
                            x = array([x])
                            y = array([y])
                            yerr = array([yerr])
                            ufmt = pfmt
                        elif stype=='scan':
                            color,style = color_wheel.next()
                            if dtype=='scalar':
                                x,y = r
                                yerr = empty([])
                            elif dtype=='stat':
                                x,yp = r
                                yp = array(yp)
                                y    = yp[:,0]
                                yerr = yp[:,1]
                            #end if
                            x    = array(x)
                            y    = array(y)
                            yerr = array(yerr)
                            for xv in x:
                                xlen = max(xlen,len(str(xv)))
                            #end for
                            ufmt = pfmt.replace('{0}','{0:<'+str(xlen)+'}')
                        #end if
                        if convertx:
                            x = convert(x,xunitt,xunits)
                        #end if
                        if converty:
                            y    = convert(y,   yunitt,units)
                            yerr = convert(yerr,yunitt,units)
                        #end if
                        if rtype=='y':
                            data = array([y])
                            pdata = array([x,y])
                        elif rtype=='ye':
                            data = array([y,yerr])
                            pdata = array([x,y,yerr])
                        elif rtype=='xy':
                            data = array([x,y])
                            pdata = data
                        elif rtype=='xye':
                            data = array([x,y,yerr])
                            pdata = data
                        #end if

                        if printing:
                            skey = str(key)
                            if len(skey)<keylen:
                                skey+=(keylen-len(skey))*' '
                            #end if
                            skey = pad+2*p+skey
                            nw,nl = data.shape
                            for i in range(nl):
                                print skey+ufmt.format(*tuple(data[:,i]))
                                skey = pad+2*p+keylen*' '
                            #end for
                        #end if
                        if plotting:
                            plotter(*pdata,color=color,fmt=style,lw=lw,label=str(tuple(pname,*key)))
                            title(ptitle)
                            xlabel(pxlabel)
                            ylabel(pylabel)
                            xticks(x)
                        #end if
                        del x,y,yerr
                    #end for

                    if printing:
                        print pad+'end '+pname
                    #end if
                #end for
                if qname!=None and printing:
                    print '    quantity {0} results'.format(qname)
                #end if
            #end for
            if printing:
                print '  end {0} stage results'.format(sname)
            #end if
        #end for        
        if printing:
            print 'end {0} results'.format(self.name)
        #end if
        if mode=='plot' and not suppress:
            show()
        #end if    
        return
    #end def show_results

#end class ValidationProcesses




class AtomicValidationStage(ValidationStage):
    systypes = ['ae','pp']
    systype = None

    stage_type = 'scan'
    data_type  = 'stat'
    yunit      = 'Ha'

    def __init__(self,**vars):
        if not self.systype in self.systypes:
            self.error('system type {0} is unrecognized\n  valid system types are: {1}'.format(systype,self.systypes))
        #end if
        self.set(**vars)
        ValidationStage.__init__(self)
        svars = obj(
            stage_inputs = self.stage_inputs,
            stage_dependencies = self.stage_dependencies
            )
        prefix = self.systype+'_'
        for svar,varlist in svars.iteritems():
            prefixed = True
            for var in varlist:
                prefixed = prefixed and var.startswith(prefix)
            #end for
            if not prefixed:
                vars = []
                for var in varlist:
                    if var.startswith(prefix):
                        vars.append(var)
                    else:
                        vars.append(prefix+var)
                    #end if
                #end for
                self[svar] = vars
            #end if
        #end for
    #end def __init__


    def get_ae_inputs(self,iq):
        inputs = obj()
        pinputs = obj()
        for name,value in self.inputs.iteritems():
            if name.startswith('ae_'):
                inputs[name.replace('ae_','')] = value[iq]
                pinputs[name] = value[iq]
            #end if
        #end for
        for name,value in self.results.iteritems():
            if name.startswith('ae_'):
                inputs[name.replace('ae_','')] = value[iq]
                pinputs[name] = value[iq]
            #end if
        #end for 
        return inputs,pinputs
    #end def get_ae_inputs


    def get_pp_inputs(self,pp,iq):
        inputs = obj()
        pinputs = obj()
        for name,value in self.inputs.iteritems():
            if name.startswith('pp_'):
                inputs[name.replace('pp_','')] = value[iq]
                pinputs[name] = value[iq]
            #end if
        #end for
        for name,value in self.results.iteritems():
            if name.startswith('pp_'):
                inputs[name.replace('pp_','')] = value[pp][iq]
                pinputs[name] = value[pp][iq]
            #end if
        #end for 
        return inputs,pinputs
    #end def get_pp_inputs


    def make_all_sims(self,basepath):
        print '    make_all_sims',self.__class__.__name__,self.systype
        sinfo = self.system_info
        if self.systype=='ae':
            for iq in range(len(sinfo.q)):
                q = sinfo.q[iq]
                path = os.path.join(basepath,'q'+str(q))
                v,vp = self.get_ae_inputs(iq)
                self.check_inputs(vp)
                v.path = path
                v.q    = q
                sims = self.make_sims(v)
                self.save_sims(sims,q,vp,v)
            #end for
        elif self.systype=='pp':
            for functional,pplist in sinfo.func_pp.iteritems():
                for pp in pplist:
                    pseudos = sinfo.ppmap[pp]
                    for iq in range(len(sinfo.q)):
                        q = sinfo.q[iq]
                        path = os.path.join(basepath,functional,pp,'q'+str(q))
                        v,vp = self.get_pp_inputs(pp,iq)
                        self.check_inputs(vp)
                        v.path       = path
                        v.q          = q
                        v.pseudos    = pseudos
                        v.functional = functional
                        v.pp         = pp
                        sims = self.make_sims(v)
                        self.save_sims(sims,(functional,pp,q),vp,v)
                    #end for
                #end for
            #end for
        #end if
    #end def make_all_sims


    def make_system(self,v,L=None,spin=None):
        atom = self.system_info.atom
        if not 'q' in v:
            self.error('charge is not present, physical system cannot be made')
        #end if
        if self.systype=='ae':
            system = generate_physical_system(
                type       = 'atom',
                atom       =  atom,
                net_charge = v.q,
                net_spin   = 'low'
                )
        elif self.systype=='pp':
            if L is None:
                if not 'L' in v:
                    self.error('box size is not present, physical system cannot be made')
                #end if
                L = v.L
            #end if
            if spin is None:
                if not 'spin' in v:
                    self.error('total spin is not present, physical system cannot be made')
                #end if
                spin = v.spin
            #end if
            valency = {atom:v.Zeff}
            system = generate_physical_system(
                lattice    = 'orthorhombic',
                cell       = 'primitive',
                centering  = 'P',
                constants  = (L,1.0000001*L,1.0000002*L),
                units      = 'A',
                atoms      = atom,
                net_charge = v.q,
                net_spin   = spin,
                kgrid      = (1,1,1),
                kshift     = (0,0,0),
                **valency
                )
            s = system.structure
            s.slide(s.axes.sum(0)/2)
        #end if
        return system
    #end def make_system


    def get_result_sim(self,v,name):
        if self.systype=='ae':
            key = v.q
            name = 'ae_'+name
        elif self.systype=='pp':
            key = v.tuple('functional','pp','q')
            name = 'pp_'+name
        #end if
        sim = ValidationStage.get_result_sim(self,(name,key))
        return sim
    #end def get_result_sim


    def get_results(self,quantity=None):
        res = obj()
        sinfo = self.system_info
        if self.systype=='ae':
            for q in sinfo.q:
                key = q
                v = self.saved_inputs[key]
                sims = self.sims[key]
                res[key] = self.get_result(v,sims,quantity)
            #end for
        elif self.systype=='pp':
            for functional,pplist in sinfo.func_pp.iteritems():
                for pp in pplist:
                    for q in sinfo.q:
                        key = (functional,pp,q)
                        v = self.saved_inputs[key]
                        sims = self.sims[key]
                        res[key] = self.get_result(v,sims,quantity)
                    #end for
                #end for
            #end for
        #end if
        return res
    #end def get_results


    def make_sims(self,v):
        self.not_implemented()
    #end def make_sims


    def get_result(self,v,sims,name):
        self.not_implemented()
    #end def get_result
#end class AtomicValidationStage



class AtomicHFOccupationScan(AtomicValidationStage):
    systype = 'ae'
    stage_dependencies = set(['ae_hfjob'])

    data_type = 'scalar'
    title  = 'Hartree-Fock energies vs. orbital occupation'
    xlabel = 'Orbital occupations'
    ylabel = 'HF energy'
    yunit  = 'Ha'

    def make_sims(self,v):
        sims = obj()
        atom = self.make_system(v)
        for occ in v.occupations:
            up,down = occ
            path = os.path.join(v.path,'{0}__{1}'.format(up,down))
            path = path.replace(',','').replace('(','').replace(')','').replace(' ','')
            hf = generate_sqd(
                identifier = 'hf',
                path       = path,
                job        = v.hfjob,
                system     = atom,
                up         = up,
                down       = down,
                grid_type  = 'log',
                ri         = 1e-6 ,
                rf         = 400  ,
                npts       = 10001,
                max_iter   = 1000 ,
                etot_tol   = 1e-8 ,
                eig_tol    = 1e-12,
                mix_ratio  = 0.7  
                )
            sims[occ] = hf
        #end for
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        occs = []
        energies = []
        for occ in v.occupations:
            a = sims[occ].load_analyzer_image()
            energies.append(a.E_tot)
            occs.append(occ[0]+' '+occ[1])
        #end for
        return occs,energies
    #end def get_result
#end class AtomicHFOccupationScan




class AtomicHFCalc(AtomicValidationStage):
    systype = 'ae'
    stage_dependencies = set(['ae_hfjob','ae_occupation'])
    stage_result = 'ae_orbitals'
    final_result = 'Ehf_ae'

    stage_type = 'value'
    data_type  = 'scalar'
    title      = 'Hartree-Fock energy vs. Ion charge'
    ylabel     = 'HF energy'
    xlabel     = 'Ion charge'
    yunit      = 'Ha'

    def make_sims(self,v):
        sims = obj()
        atom = self.make_system(v)
        up,down = v.occupation
        hf = generate_sqd(
            identifier = 'hf',
            path       = v.path,
            job        = v.hfjob,
            system     = atom,
            up         = up,
            down       = down,
            grid_type  = 'log',
            ri         = 1e-6 ,
            rf         = 400  ,
            npts       = 10001,
            max_iter   = 1000 ,
            etot_tol   = 1e-8 ,
            eig_tol    = 1e-12,
            mix_ratio  = 0.7  
            )
        sims.orb = hf
        sims.result = hf
        return sims
    #end def make_sims


    def get_result(self,v,sims,name='E'):
        ha = sims.hf.load_analyzer_image()
        if name=='E':
            res = ha.E
        elif name=='B':
            res = 1./ha.moment(n=1)
        else:
            self.error(name+' is not a valid quantity\n valid options are E, B')
        #end if
        return res
    #end def get_result
#end class AtomicHFCalc




class AtomicDFTBoxScan(AtomicValidationStage):
    systype = 'pp'
    stage_inputs       = set(['pp_Ls'])
    stage_dependencies = set(['pp_dftjob','pp_Ecut0','pp_spin0','pp_assume_isolated0','pp_Zeff'])

    data_type = 'scalar'
    title  = 'DFT energy vs. box size'
    ylabel = 'DFT energy'
    xlabel = 'box size'
    yunit  = 'Ry'
    xunit  = 'A'

    def make_sims(self,v):
        sims = obj()
        dftpp,qmcpp = sort_pseudos(v.pseudos)
        for L in v.Ls:
            atom = self.make_system(v,L=L,spin=v.spin0)
            s = atom.structure
            s.slide(s.axes.sum(0)/2)
            path = os.path.join(v.path,'L_'+str(L))
            scf = generate_pwscf(
                identifier   = 'scf',
                path         = path,
                job          = v.dftjob,
                input_type   = 'scf',
                input_dft    = v.functional,
                ecut         = v.Ecut0,
                conv_thr     = 1e-8,
                mixing_beta  = .7,
                nosym        = True,
                assume_isolated = v.assume_isolated0,
                pseudos      = dftpp,
                system       = atom,
                kgrid        = (1,1,1),
                kshift       = (0,0,0),
                use_folded   = False
                )
            sims[L] = scf
        #end for
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        Ls = []
        energies = []
        for L in v.Ls:
            pa = sims[L].load_analyzer_image()
            Ls.append(L)
            energies.append(pa.E)
        #end for
        return Ls,energies
    #end def get_result
#end class AtomicDFTBoxScan




class AtomicDFTEcutScan(AtomicValidationStage):
    systype = 'pp'
    stage_inputs       = set(['pp_Ecuts'])
    stage_dependencies = set(['pp_dftjob','pp_p2cjob','pp_L','pp_assume_isolated','pp_spin0','pp_Zeff'])

    data_type = 'scalar'
    title  = 'DFT energy vs. planewave energy cutoff'
    ylabel = 'DFT energy'
    xlabel = 'Ecut'
    yunit  = 'Ry'
    xunit  = 'Ry'

    def make_sims(self,v):
        sims = obj()
        dftpp,qmcpp = sort_pseudos(v.pseudos)
        atom = self.make_system(v,spin=v.spin0)
        for ecut in v.Ecuts:
            path = os.path.join(v.path,'Ecut_'+str(ecut))
            scf = generate_pwscf(
                identifier   = 'scf',
                path         = path,
                job          = v.dftjob,
                input_type   = 'scf',
                input_dft    = v.functional,
                ecut         = ecut,
                conv_thr     = 1e-8,
                mixing_beta  = .7,
                nosym        = True,
                assume_isolated = v.assume_isolated,
                pseudos      = dftpp,
                system       = atom,
                kgrid        = (1,1,1),
                kshift       = (0,0,0),
                use_folded   = False
                )
            p2c = generate_pw2casino(
                identifier   = 'p2c',
                path         = path,
                job          = v.p2cjob
                )
            p2c.depends(scf,'orbitals')
            sims[ecut,'E'] = scf
            sims[ecut,'KE'] = p2c
        #end for
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        ecuts = []
        energies = []
        if not name in ['E','KE']:
            self.error('quantity {0} is not computed by {1}'.format(name,self.name))
        #end if
        for ecut in v.Ecuts:
            ecuts.append(ecut)
            pa = sims[ecut,name].load_analyzer_image()
            if name=='E':
                energies.append(pa.E)
            elif name=='KE':
                energies.append(pa.energies.Kinetic)
            #end if
        #end for
        return ecuts,energies
    #end def get_result
#end class AtomicDFTEcutScan




class AtomicDFTSpinScan(AtomicValidationStage):
    systype = 'pp'
    stage_inputs       = set(['pp_spins'])
    stage_dependencies = set(['pp_dftjob','pp_L','pp_Ecut','pp_assume_isolated','pp_spin0','pp_Zeff'])

    data_type = 'scalar'
    title  = 'DFT energy vs. spin state'
    ylabel = 'DFT energy'
    xlabel = 'Total spin (in units of 1/2)'
    yunit  = 'Ry'

    def make_sims(self,v):
        sims = obj()
        dftpp,qmcpp = sort_pseudos(v.pseudos)
        for spin in v.spins:
            atom = self.make_system(v,spin=spin)
            path = os.path.join(v.path,'spin_'+str(spin))
            scf = generate_pwscf(
                identifier   = 'scf',
                path         = path,
                job          = v.dftjob,
                input_type   = 'scf',
                input_dft    = v.functional,
                ecut         = v.Ecut,
                conv_thr     = 1e-8,
                mixing_beta  = .7,
                nosym        = True,
                assume_isolated = v.assume_isolated,
                pseudos      = dftpp,
                system       = atom,
                kgrid        = (1,1,1),
                kshift       = (0,0,0),
                use_folded   = False
                )
            sims[spin] = scf
        #end for
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        spins = []
        energies = []
        for spin in v.spins:
            spin.append(spin)
            pa = sims[spin].load_analyzer_image()
            energies.append(pa.E)
        #end for
        return spins,energies
    #end def get_result
#end class AtomicDFTSpinScan




class AtomicDFTCalc(AtomicValidationStage):
    systype = 'pp'
    stage_dependencies = set(['pp_dftjob','pp_p2qjob','pp_Zeff','pp_L','pp_assume_isolated','pp_Ecut','pp_spin'])
    stage_result = 'pp_orbitals'
    final_result = 'Edft_pp'

    stage_type = 'value'
    data_type  = 'scalar'
    title      = 'DFT energy vs Ion charge'
    ylabel     = 'DFT energy'
    xlabel     = 'Ion charge'
    yunit      = 'Ry'

    def make_sims(self,v):
        dftpp,qmcpp = sort_pseudos(v.pseudos)
        atom = self.make_system(v)
        scf = generate_pwscf(
            identifier   = 'scf',
            path         = v.path,
            job          = v.dftjob,
            input_type   = 'scf',
            input_dft    = v.functional,
            ecut         = v.Ecut,
            conv_thr     = 1e-8,
            mixing_beta  = .7,
            nosym        = True,
            assume_isolated = v.assume_isolated,
            pseudos      = dftpp,
            system       = atom,
            kgrid        = (1,1,1),
            kshift       = (0,0,0),
            use_folded   = False
            )
        p2q = generate_pw2qmcpack(
            identifier   = 'p2q',
            path         = v.path,
            job          = v.p2qjob,
            write_psir   = False
            )
        p2q.depends(scf,'orbitals')
        sims = obj(
            scf = scf,
            p2q = p2q,
            result = p2q
            )
        return sims
    #end def make_sims


    def get_result(self,v,sims,name='E'):
        pa = sims.scf.load_analyzer_image()
        return pa.E
    #end def get_result
#end class AtomicDFTCalc




class AtomicOptJ1RcutScan(AtomicValidationStage):
    systype = 'pp'
    stage_inputs       = set(['pp_J1_rcuts','pp_opt_calcs'])
    stage_dependencies = set(['pp_Zeff','pp_L','pp_spin','pp_optjob','pp_orbitals','pp_pade_b'])
    stage_result       = 'pp_J1_jastrow'

    title  = 'Optimal VMC Energy vs. J1 rcut'
    ylabel = 'Opt. Energy'
    xlabel = 'J1 rcut'
    yunit  = 'Ha'
    xunit  = 'B'

    def make_sims(self,v):
        sims = obj()
        res = obj()
        dftpp,qmcpp = sort_pseudos(v.pseudos)
        atom = self.make_system(v)
        orb  = self.get_result_sim(v,'orbitals')
        bu,bd = v.pade_b
        for rcut in v.J1_rcuts:
            path = os.path.join(v.path,'rcut_'+str(rcut))
            jastrows = [
                generate_jastrow('J1','bspline',8,rcut,system=atom),
                generate_jastrow('J2','pade',bu,bd,system=atom)
                ]
            opt = generate_qmcpack(
                identifier = 'opt',
                path       = path,
                job        = v.optjob,
                input_type = 'opt_jastrow',
                system     = atom,
                bconds     = 'nnn',
                pseudos    = qmcpp,
                jastrows   = jastrows,
                corrections = [],
                opt_calcs  = v.opt_calcs
                )
            opt.depends(orb,'orbitals')
            sims[rcut] = opt
            res[rcut] = opt
        #end for
        sims.result = res
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name=None):
        rcuts    = []
        energies = []
        for rcut in v.J1_rcuts:
            qa = sims[rcut].load_analyzer_image()
            if 'vmc' in qa and len(qa.vmc)>0:
                series = array(qa.vmc.keys())
                qmc = qa.vmc[series.max()]
            else:
                series = array(qa.opt.keys())
                qmc = qa.opt[series.max()]
            #end if
            le = qmc.scalar.LocalEnergy
            rcuts.append(rcut)
            energies.append((le.mean,le.error))
        #end for
        return rcuts,energies
    #end def get_result
#end class AtomicOptJ1RcutScan




class AtomicOptCalc(AtomicValidationStage):

    title  = 'Opt. Energy vs. iteration #'
    ylabel = 'Opt. Energy'
    xlabel = 'iteration #'
    yunit  = 'Ha'

    def make_sims(self,v):
        atom = self.make_system(v)
        orb = self.get_result_sim(v,'orbitals')
        if self.systype=='ae':
            qmcpp = None
            bu,bd = v.pade_b
            J3 = generate_jastrow('J3','polynomial',4,4,5.0,system=atom)
            J3.source = 'atom'
            jastrows = [
                generate_jastrow('J2','pade',bu,bd,system=atom),
                J3
                ]
        elif self.systype=='pp':
            dftpp,qmcpp = sort_pseudos(v.pseudos)
            jastrows = []
        #end if
        opt = generate_qmcpack(
            identifier = 'opt',
            path       = v.path,
            job        = v.optjob,
            input_type = 'opt_jastrow',
            system     = atom,
            bconds     = 'nnn',
            pseudos    = qmcpp,
            jastrows   = jastrows,
            corrections = [],
            opt_calcs  = v.opt_calcs
            )
        opt.depends(orb,'orbitals')
        if self.systype=='pp':
            J1_jastrows = self.get_result_sim(v,'J1_jastrow')
            if not v.J1_rcut in J1_jastrows:
                self.error('requested rcut {0} could not be found in J1 Jastrow rcuts\n  possible rcuts: {1}'.format(v.J1_rcut,J1_jastrows.keys()))
            #end if
            jsim = J1_jastrows[v.J1_rcut]
            opt.depends(jsim,'jastrow')
        #end if
        sims = obj(
            opt    = opt,
            result = opt
            )
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        qa = sims.opt.load_analyzer_image()
        series = []
        energies = []
        for s,opt in qa.opt.iteritems():
            series.append(s)
            le = opt.scalar.LocalEnergy
            energies.append((le.mean,le.error))
        #end for
        return series,energies
    #end def get_result
#end class AtomicOptCalc


class AtomicAEOptCalc(AtomicOptCalc):
    systype = 'ae'
    stage_inputs = set(['ae_opt_calcs'])
    stage_dependencies = set(['ae_optjob','ae_occupation','ae_orbitals','ae_pade_b'])
    stage_result = 'ae_jastrow'
#end class AtomicAEOptCalc


class AtomicPPOptCalc(AtomicOptCalc):
    systype = 'pp'
    stage_inputs = set(['pp_opt_calcs'])
    stage_dependencies = set(['pp_Zeff','pp_L','pp_spin','pp_optjob','pp_orbitals','pp_J1_rcut'])
    stage_result = 'pp_jastrow'
#end class AtomicPPOptCalc




class AtomicVMCCalc(AtomicValidationStage):

    stage_type = 'value'
    title  = 'VMC Energy vs. Ion charge'
    ylabel = 'VMC Energy'
    xlabel = 'Ion charge'
    yunit  = 'Ha'

    def make_sims(self,v):
        atom = self.make_system(v)
        orb = self.get_result_sim(v,'orbitals')
        jastrow = self.get_result_sim(v,'jastrow')
        if self.systype=='ae':
            qmcpp = None
        elif self.systype=='pp':
            dftpp,qmcpp = sort_pseudos(v.pseudos)
        #end if
        vmc = generate_qmcpack(
            identifier   = 'vmc',
            path         = v.path,
            job          = v.vmcjob,
            input_type   = 'basic',
            system       = atom,
            bconds       = 'nnn',
            pseudos      = qmcpp,
            jastrows     = [],
            corrections  = [],
            calculations = v.vmc_calcs
            )
        vmc.depends(orb,'orbitals')
        vmc.depends(jastrow,'jastrow')
        sims = obj(vmc=vmc)
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        qa = sims.vmc.load_analyzer_image()
        qmc = qa.qmc
        vmc = qmc[len(qmc)-1]
        le = vmc.scalar.LocalEnergy
        return (le.mean,le.error)
    #end def get_result
#end class AtomicVMCCalc


class AtomicAEVMCCalc(AtomicVMCCalc):
    systype = 'ae'
    stage_inputs = set(['ae_vmc_calcs'])
    stage_dependencies = set(['ae_vmcjob','ae_occupation','ae_orbitals','ae_jastrow'])
    final_result = 'Evmc_ae'
#end class AtomicAEVMCCalc


class AtomicPPVMCCalc(AtomicVMCCalc):
    systype = 'pp'
    stage_inputs = set(['pp_vmc_calcs'])
    stage_dependencies = set(['pp_Zeff','pp_L','pp_spin','pp_vmcjob','pp_orbitals','pp_jastrow'])
    final_result = 'Evmc_pp'
#end class AtomicPPVMCCalc




class AtomicDMCPopulationScan(AtomicValidationStage):

    title  = 'DMC Energy vs. Walker population'
    ylabel = 'DMC Energy'
    xlabel = 'Population'
    yunit  = 'Ha'

    def make_sims(self,v):
        sims = obj()
        atom = self.make_system(v)
        orb = self.get_result_sim(v,'orbitals')
        jastrow = self.get_result_sim(v,'jastrow')
        if self.systype=='ae':
            qmcpp = None
        elif self.systype=='pp':
            dftpp,qmcpp = sort_pseudos(v.pseudos)
        #end if
        for population in v.populations:
            calcs = deepcopy(v.dmc_calcs)
            found_vmc = False
            found_dmc = False
            for calc in calcs:
                if isinstance(calc,vmc):
                    if 'samplesperthread' in calc:
                        del calc.samplesperthread
                    #end if
                    calc.samples = population
                    found_vmc = True
                elif isinstance(calc,dmc):
                    found_dmc = True
                #end if
            #end for
            if not found_vmc or not found_dmc:
                self.error('vmc and dmc blocks must be present in dmc_calcs')
            #end if
            path = os.path.join(v.path,'pop_'+str(population))
            qmc = generate_qmcpack(
                identifier   = 'dmc',
                path         = path,
                job          = v.dmcjob,
                input_type   = 'basic',
                system       = atom,
                bconds       = 'nnn',
                pseudos      = qmcpp,
                jastrows     = [],
                corrections  = [],
                calculations = calcs
                )
            qmc.depends(orb,'orbitals')
            qmc.depends(jastrow,'jastrow')
            sims[population] = qmc
        #end if
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        pops = []
        energies = []
        for pop in v.populations:
            qa = sims[pop].load_analyzer_image()
            qmc = qa.qmc
            le = qmc[len(qmc)-1].dmc.LocalEnergy
            pops.append(pop)
            energies.append((le.mean,le.error))
        #end for
    #end def get_result
#end class AtomicDMCPopulationScan


class AtomicAEDMCPopulationScan(AtomicDMCPopulationScan):
    systype = 'ae'
    stage_inputs = set(['ae_dmc_calcs','ae_populations'])
    stage_dependencies = set(['ae_dmcjob','ae_occupation','ae_orbitals','ae_jastrow'])
#end class AtomicAEDMCPopulationScan


class AtomicPPDMCPopulationScan(AtomicDMCPopulationScan):
    systype = 'pp'
    stage_inputs = set(['pp_dmc_calcs','pp_populations'])
    stage_dependencies = set(['pp_Zeff','pp_L','pp_spin','pp_dmcjob','pp_orbitals','pp_jastrow'])
#end class AtomicPPDMCPopulationScan




class AtomicDMCTimestepScan(AtomicValidationStage):

    title  = 'DMC Energy vs. Timestep'
    ylabel = 'DMC Energy'
    xlabel = 'Timestep'
    yunit  = 'Ha'
    xunit  = '1/Ha'

    def make_sims(self,v):
        sims = obj()
        atom = self.make_system(v)
        orb = self.get_result_sim(v,'orbitals')
        jastrow = self.get_result_sim(v,'jastrow')
        if self.systype=='ae':
            qmcpp = None
        elif self.systype=='pp':
            dftpp,qmcpp = sort_pseudos(v.pseudos)
        #end if
        calcs = deepcopy(v.dmc_calcs)
        found_vmc = False
        found_dmc = False
        vmc_calc = None
        dmc_calc = None
        for calc in calcs:
            if isinstance(calc,vmc):
                if 'samplesperthread' in calc:
                    del calc.samplesperthread
                #end if
                calc.samples = v.population
                vmc_calc = calc
                found_vmc = True
            elif isinstance(calc,dmc):
                dmc_calc = calc
                found_dmc = True
            #end if
        #end for
        if not found_vmc or not found_dmc:
            self.error('vmc and dmc blocks must be present in dmc_calcs')
        #end if
        calcs = [vmc_calc]
        for timestep in v.timesteps:
            dcalc = dmc_calc.copy()
            tfac = dcalc.timestep/timestep
            dcalc.warmupsteps = int(round(tfac*dcalc.warmupsteps))
            dcalc.steps       = int(round(tfac*dcalc.steps))
            dcalc.timestep    = timestep
            calcs.append(dcalc)
        #end for
        qmc = generate_qmcpack(
            identifier   = 'dmc',
            path         = v.path,
            job          = v.dmcjob,
            input_type   = 'basic',
            system       = atom,
            bconds       = 'nnn',
            pseudos      = qmcpp,
            jastrows     = [],
            corrections  = [],
            calculations = calcs
            )
        qmc.depends(orb,'orbitals')
        qmc.depends(jastrow,'jastrow')
        sims.dmc = qmc
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        s=0
        timesteps = []
        energies  = []
        qa = sims.dmc.load_analyzer_image()
        dmc = qa.dmc
        for timestep in v.timesteps:
            s+=1
            le = dmc[s].dmc.LocalEnergy
            timesteps.append(timestep)
            energies.append((le.mean,le.error))
        #end for
        return timesteps,energies
    #end def get_result
#end class AtomicDMCTimestepScan


class AtomicAEDMCTimestepScan(AtomicDMCTimestepScan):
    systype = 'ae'
    stage_inputs = set(['ae_dmc_calcs'])
    stage_dependencies = set(['ae_dmcjob','ae_occupation','ae_orbitals','ae_jastrow','ae_population'])
#end class AtomicAEDMCTimestepScan


class AtomicPPDMCTimestepScan(AtomicDMCTimestepScan):
    systype = 'pp'
    stage_inputs = set(['pp_dmc_calcs'])
    stage_dependencies = set(['pp_Zeff','pp_L','pp_spin','pp_dmcjob','pp_orbitals','pp_jastrow','pp_population'])
#end class AtomicPPDMCTimestepScan




class AtomicDMCCalc(AtomicValidationStage):

    stage_type = 'value'
    title  = 'DMC Energy vs. Ion charge'
    ylabel = 'DMC Energy'
    xlabel = 'Ion charge'
    yunit  = 'Ha'

    def make_sims(self,v):
        atom = self.make_system(v)
        orb     = self.get_result_sim(v,'orbitals')
        jastrow = self.get_result_sim(v,'jastrow')
        if self.systype=='ae':
            qmcpp = None
        elif self.systype=='pp':
            dftpp,qmcpp = sort_pseudos(v.pseudos)
        #end if
        calcs = deepcopy(v.dmc_calcs)
        found_vmc = False
        found_dmc = False
        for calc in calcs:
            if isinstance(calc,vmc):
                if 'samplesperthread' in calc:
                    del calc.samplesperthread
                #end if
                calc.samples = v.population
                found_vmc = True
            elif isinstance(calc,dmc):
                dmc_calc = calc
                found_dmc = True
            #end if
        #end for
        if found_dmc:
            stepfac = calc.timestep/v.timestep
            dmc_calc.warmupsteps = int(round(dmc_calc.warmupsteps*stepfac))
            dmc_calc.steps = int(round(dmc_calc.steps*stepfac))
            dmc_calc.timestep = v.timestep
        #end if

        if not found_vmc or not found_dmc:
            self.error('vmc and dmc blocks must be present in dmc_calcs')
        #end if
        qmc = generate_qmcpack(
            identifier   = 'dmc',
            path         = v.path,
            job          = v.dmcjob,
            input_type   = 'basic',
            system       = atom,
            pseudos      = qmcpp,
            jastrows     = [],
            corrections  = [],
            calculations = calcs
            )
        qmc.depends(orb,'orbitals')
        qmc.depends(jastrow,'jastrow')
        sims = obj(dmc=qmc)
        return sims
    #end def make_sims

        
    def get_result(self,v,sims,name='E'):
        qa = sims.dmc.load_analyzer_image()
        qmc = qa.qmc
        dmc = qmc[len(qmc)-1]
        le = dmc.dmc.LocalEnergy
        return (le.mean,le.error)
    #end def get_result
#end class AtomicDMCCalc


class AtomicAEDMCCalc(AtomicDMCCalc):
    systype = 'ae'
    stage_inputs = set(['ae_dmc_calcs'])
    stage_dependencies = set(['ae_dmcjob','ae_orbitals','ae_jastrow','ae_population','ae_timestep'])
    final_result = 'Edmc_ae'
#end class AtomicAEDMCCalc


class AtomicPPDMCCalc(AtomicDMCCalc):
    systype = 'pp'
    stage_inputs = set(['pp_dmc_calcs'])
    stage_dependencies = set(['pp_Zeff','pp_L','pp_spin','pp_dmcjob','pp_orbitals','pp_jastrow','pp_population','pp_timestep'])
    final_result = 'Edmc_pp'
#end class AtomicPPDMCCalc





class ValidateAtomPP(ValidationProcess):

    allowed_stage_inputs = set(
        ['ae_occupations', 'ae_opt_calcs', 'ae_vmc_calcs', 'ae_dmc_calcs', 
         'ae_populations', 'ae_timesteps', 
         'pp_Ls', 'pp_Ecuts','pp_spins', 'pp_opt_calcs', 'pp_J1_rcuts', 
         'pp_vmc_calcs', 'pp_populations', 'pp_timesteps', 'pp_dmc_calcs'])

    allowed_stage_results = set(
        ['ae_hfjob', 'ae_optjob', 'ae_vmcjob', 'ae_dmcjob', 'ae_pade_b',
         'ae_occupation','ae_orbitals','ae_jastrow','ae_population','ae_timestep',
         'pp_dftjob', 'pp_optjob', 'pp_vmcjob', 'pp_dmcjob', 'pp_pade_b',
         'pp_p2cjob','pp_p2qjob',
         'pp_Ecut0','pp_spin0','pp_Zeff','pp_L','pp_Ecut','pp_spin','pp_orbitals',
         'pp_assume_isolated0','pp_assume_isolated','pp_J1_rcut','pp_jastrow',
         'pp_population','pp_timestep'])
        

    def initialize(self,**func_pp):
        if not 'q' in func_pp:
            self.error('variable q must be provided')
        #end if
        self.q = func_pp['q']
        del func_pp['q']
        ref = None
        if 'ref' in func_pp:
            ref = func_pp['ref']
            del func_pp['ref']
        #end if
        if 'name' in func_pp:
            self.name = func_pp['name']
            del func_pp['name']
        #end if
        self.check_funcs_pseudos(**func_pp)
        self.ref = ref

        self.add_stages(
            ae_occupation_scan = AtomicHFOccupationScan(),
            ae_hf_calc         = AtomicHFCalc(),
            ae_opt_calc        = AtomicAEOptCalc(),
            ae_vmc_calc        = AtomicAEVMCCalc(),
            ae_dmc_pop_scan    = AtomicAEDMCPopulationScan(),
            ae_dmc_tau_scan    = AtomicAEDMCTimestepScan(),
            ae_dmc_calc        = AtomicAEDMCCalc(),
            pp_box_scan        = AtomicDFTBoxScan(),
            pp_ecut_scan       = AtomicDFTEcutScan(),
            pp_spin_scan       = AtomicDFTSpinScan(),
            pp_dft_calc        = AtomicDFTCalc(),
            pp_J1_rcut_scan    = AtomicOptJ1RcutScan(),
            pp_opt_calc        = AtomicPPOptCalc(),
            pp_vmc_calc        = AtomicPPVMCCalc(),
            pp_dmc_pop_scan    = AtomicPPDMCPopulationScan(),
            pp_dmc_tau_scan    = AtomicPPDMCTimestepScan(),
            pp_dmc_calc        = AtomicPPDMCCalc()
            )

        self.set_stage_order(
            ['ae_occupation_scan','ae_hf_calc','ae_opt_calc','ae_vmc_calc',
             'ae_dmc_pop_scan','ae_dmc_tau_scan','ae_dmc_calc',
             'pp_box_scan','pp_ecut_scan','pp_spin_scan','pp_dft_calc',
             'pp_J1_rcut_scan','pp_opt_calc','pp_vmc_calc','pp_dmc_pop_scan',
             'pp_dmc_tau_scan','pp_dmc_calc']
            )

        info = obj(
            q = self.q,
            ppmap = self.ppmap,
            func_pp = self.func_pp,
            functionals = self.functionals,
            atom = self.name
            )
        for stage in self.stages:
            stage.system_info = info.copy()
        #end for

    #end def initialize
#end class ValidateAtomPP




class ValidateAtomPPs(ValidationProcesses):
    def initialize(self):
        atoms = list(self.processes.keys())
        pseudos = obj()
        for atom in atoms:
            pseudos[atom] = self.processes[atom].pseudo_list()
        #end for
        pp_to_atom = obj()
        for atom,pplist in pseudos.iteritems():
            for pp in pplist:
                pp_to_atom[pp] = atom
            #end for
        #end for
        self.atoms = atoms
        self.pseudos = pseudos
        self.pp_to_atom = pp_to_atom
        nqset = set()
        nq = obj()
        charges = obj()
        for atom,process in self.processes.iteritems():
            nq[atom] = len(process.q)
            nqset.add(len(process.q))
            charges[atom] = list(process.q)
        #end for
        self.charges = charges
        self.nq = nq
        self.nq_same = len(nqset)==1
    #end def initialize


    def set_stage_inputs(self,**kwargs):
        inputs = obj()
        atoms = set(self.processes.keys())
        errors = False
        for name,inp in kwargs.iteritems():
            input = self.expand_stage_input(inp,name)
            if not atoms<=set(input.keys()):
                self.error('all atoms must be specified in input variable {0}\n  atoms required: {1}\n  atoms specified: {2}'.format(name,list(atoms),input.keys()),exit=False,trace=False)
                errors = True
            #end if
            inputs[name] = input
        #end for
        if errors:
            self.error('errors encountered')
        #end if
        for atom,process in self.processes.iteritems():
            pinputs = obj()
            for name,input in inputs.iteritems():
                pinputs[name] = inputs[name][atom]
            #end for
            pi = pinputs.copy()
            pkeys = list(pi.keys())
            for key in pkeys:
                if key.endswith('job') or key.endswith('calcs'):
                    del pi[key]
                #end if
            #end for
            process.set_stage_inputs(**pinputs)
        #end for
    #end def set_stage_inputs


    def set_stage_results(self,**kwargs):
        ae_res = obj()
        pp_res = obj()
        for name,result in kwargs.iteritems():
            if name.startswith('ae_'):
                ae_res[name] = self.expand_ae_stage_result(result)
            elif name.startswith('pp_'):
                pp_res[name] = self.expand_pp_stage_result(result)
            else:
                self.error('invalid stage result variable encountered\n  variable encountered {0}\n  variable must be prefixed with ae_ or pp_'.format(name))
            #end if
        #end for

        for atom,process in self.processes.iteritems():
            res = obj()
            for name,rcoll in ae_res.iteritems():
                res[name] = deepcopy(rcoll[atom])
            #end for
            for name,rcoll in pp_res.iteritems():
                res[name] = deepcopy(rcoll[atom])
            #end for
            process.check_stage_results(res.keys())
            process.set_stage_results(**res)
        #end for
    #end def set_stage_results


    def expand_stage_input(self,iin,name=None):
        i = obj()
        if isinstance(iin,obj):
            atoms = set(iin.keys())
            if len(set(self.atoms)-atoms)>0:
                self.error('cannot expand stage input {0}\n  information is required for all requested atoms\n  atoms requested: {1}\n  information provided:\n{2}'.format(name,self.atoms,iin))
            #end if
            for atom in self.atoms:
                ainfo = iin[atom]
                nq = len(self.charges[atom])
                if isinstance(ainfo,list):
                    nlists = 0
                    for elem in ainfo:
                        if isinstance(elem,list):
                            nlists+=1
                        #end if
                    #end if
                    nelem = len(ainfo)
                    if nlists==0:
                        i[atom] = nq*[ainfo]
                    elif nlists==nelem:
                        if nelem==nq:
                            i[atom] = ainfo
                        else:
                            self.error('cannot expand stage input {0}\n  number of elements provided for atom {1} does not match the number of charge states\n  number of elements provided: {2}\n  number of charge states: {3}'.format(name,atom,nelem,nq))
                        #end if
                    else:
                        self.error('cannot expand stage input {0}\n  some elements in the list provided for atom {1} are lists and some are not\n  if all are lists, it is assumed that each list corresponds to a different charge state\n  if none are lists, the elments are taken to be scan parameters applied the same to every charge state'.format(name,element))
                    #end if
                else:
                    self.error('invalid type for stage input {0} atom {1}\n  type encountered: {2}\n  valid types are: list'.format(name,atom,ainfo.__class__.__name__))
                #end if
            #end for
        elif isinstance(iin,list):
            for atom in self.atoms:
                nq = len(self.charges[atom])
                i[atom] = nq*[iin]
            #end for
        else:
            self.error('invalid type for stage input {0}\n  type encountered: {1}\n  valid types are: obj,list'.format(name,iin.__class__.__name__))
        #end if
        return i
    #end def expand_stage_input


    def expand_ae_stage_result(self,rin):
        r = obj()
        single_types = (str,int,float,tuple,type(None),Job)
        wrong_type = False
        if isinstance(rin,single_types):
            for atom in self.atoms:
                r[atom] = self.nq[atom]*[rin]
            #end for
        elif isinstance(rin,obj):
            atoms = set(rin.keys())
            if len(set(self.atoms)-atoms)>0:
                self.error('cannot expand stage result\n  information is required for all requested atoms\n  atoms requested: {0}\n  information provided:\n{1}'.format(self.atoms,rin))
            #end if
            for atom,ainfo in rin.iteritems():
                if atom in self.atoms:
                    ares = None
                    if isinstance(ainfo,list):
                        if not len(ainfo)==self.nq[atom]:
                            self.error('stage result list for atom {0} is not the same length as the charge list\n  result list: {1}\n  charge list: {2}'.format(atom,ainfo,self.processes[atom].q))
                        #end if
                        ares = list(ainfo)
                    elif isinstance(ainfo,single_types):
                        ares = self.nq[atom]*[ainfo]
                    else:
                        wrong_type = True
                    #end if
                    r[atom] = ares
                #end if
            #end for
        elif isinstance(rin,list):
            if not self.nq_same:
                self.error('charge lists are not the same length for each atom\n  stage results cannot be set by list provided\n  list provided: '+str(rin))
            #end if
            for atom in self.atoms:
                r[atom] = list(rin)
            #end for
        else:
            wrong_type = True
        #end if

        if wrong_type:
            self.error('cannot expand stage result\n  received type: {0}\n  allowed_types: {1}'.format(rin.__class__.__name__,single_types))
        #end if
        return r
    #end def expand_ae_stage_result


    def expand_pp_stage_result(self,rin):
        r = obj()
        single_types = (str,int,float,tuple,type(None),Job)
        wrong_type = False
        if isinstance(rin,single_types):
            for atom in self.atoms:
                ares = obj()
                for pp in self.pseudos[atom]:
                    ares[pp] = self.nq[atom]*[rin]
                #end for
                r[atom] = ares
            #end for
        elif isinstance(rin,obj):
            atoms = set(rin.keys())
            if len(set(self.atoms)-atoms)>0:
                self.error('cannot expand stage result\n  information is required for all requested atoms\n  atoms requested: {0}\n  information provided:\n{1}'.format(self.atoms,rin))
            #end if
            for atom,ainfo in rin.iteritems():
                if atom in self.atoms:
                    ares = obj()
                    if isinstance(ainfo,dict):
                        pps = set(ainfo.keys())
                        if len(set(self.pseudos[atom])-pps)>0:
                            self.error('cannot expand stage result\n  information is required for all requested pseudopotentials\n  pseudopotentials requested: {0}\n  information provided:\n{1}'.format(self.pseudos[atom],rin))
                        #end if
                        for pp,ppinfo in ainfo.iteritems():
                            if isinstance(ppinfo,list):
                                ares[pp] = list(ppinfo)
                            elif isinstance(ppinfo,single_types):
                                ares[pp] = self.nq[atom]*[ppinfo]
                            else:
                                wrong_type = True
                            #end if
                        #end for
                    elif isinstance(ainfo,list):
                        if not len(ainfo)==self.nq[atom]:
                            self.error('stage result list for atom {0} is not the same length as the charge list\n  result list: {1}\n  charge list: {2}'.format(atom,ainfo,self.processes[atom].q))
                        #end if
                        for pp in self.pseudos[atom]:
                            ares[pp] = list(ainfo)
                        #end for
                    elif isinstance(ainfo,single_types):
                        for pp in self.pseudos[atom]:
                            ares[pp] = self.nq[atom]*[ainfo]
                        #end for
                    else:
                        wrong_type = True
                    #end if
                    r[atom] = ares
                #end if
            #end for
        elif isinstance(rin,dict):
            pps = set(rin.keys())
            if len(set(self.pp_to_atom.keys())-pps)>0:
                self.error('cannot expand stage result\n  information is required for all requested pseudopotentials\n  pseudopotentials requested: {0}\n  information provided:\n{1}'.format(self.pp_to_atom.keys(),rin))
            #end if
            for pp,atom in self.pp_to_atom.iteritems():
                if pp in rin:
                    if not atom in r:
                        r[atom] = obj()
                    #end if
                    ppinfo = rin[pp]
                    ppres = obj()
                    if isinstance(ppinfo,list):
                        if not len(ppinfo)==self.nq[atom]:
                            self.error('stage result list for atom {0} is not the same length as the charge list\n  result list: {1}\n  charge list: {2}'.format(atom,ppinfo,self.processes[atom].q))
                        #end if
                        ppres = list(ppinfo)
                    elif isinstance(ppinfo,single_types):
                        ppres = self.nq[atom]*[ppinfo]
                    else:
                        wrong_type = True
                    #end if
                    r[atom][pp] = ppres
                #end if
            #end for
        elif isinstance(rin,list):
            if not self.nq_same:
                self.error('charge lists are not the same length for each atom\n  stage results cannot be set by list provided\n  list provided: '+str(rin))
            #end if
            for atom,pseudos in self.pseudos.iteritems():
                ares = obj()
                for pp in pseudos:
                    ares[pp] = list(rin)
                #end for
                r[atom] = ares
            #end for
        else:
            wrong_type = True
        #end if

        if wrong_type:
            self.error('cannot expand stage result\n  received type: {0}\n  allowed_types: {1}'.format(rin.__class__.__name__,single_types))
        #end if

        return r
    #end def expand_pp_stage_result
        
#end class ValidateAtomPPs
ValidateAtomPPs.set_process_type(ValidateAtomPP)




class DimerValidationStage(ValidationStage):
    stage_type = 'scan'
    data_type  = 'stat'


    def get_pp_inputs(self,pp):
        inputs = obj()
        for name,value in self.results.iteritems():
            inputs[name] = value[pp]
        #end for
        return inputs
    #end def get_pp_inputs


    def make_all_sims(self,basepath):
        print '    make_all_sims',self.__class__.__name__
        sinfo = self.system_info
        for functional,pplist in sinfo.func_pp.iteritems():
            for pp in pplist:
                if not isinstance(pp,tuple) or not len(pp)==2:
                    self.error('expected length 2 tuple key for pseudopotentials\n  received: {0}'.format(pp))
                #end if
                pseudos = sinfo.ppmap[pp]
                path = os.path.join(basepath,functional,pp[0]+'_'+pp[1])
                v = self.get_pp_inputs(pp)
                self.check_inputs(v)
                v.path       = path
                v.pseudos    = pseudos
                v.functional = functional
                v.pp         = pp
                sims = self.make_sims(v)
                self.save_sims(sims,(functional,pp[0],pp[1]),v)
            #end for
        #end for
    #end def make_all_sims


    def make_system(self,v,bond_length=None):
        if bond_length is None:
            if not 'bond_length' in v:
                self.error('bond length is not present, physical system cannot be made')
            #end if
            bond_length = v.bond_length
        #end if
        if not 'L' in v:
            self.error('box size is not present, physical system cannot be made')
        #end if
        L = v.L
        self.error('varied pseudopotential Zeff has not yet been accounted for\n  please contact the developer')
        system = generate_physical_system(
            lattice       = 'orthorhombic',
            cell          = 'primitive',
            centering     = 'P',
            constants     = (L,1.0000001*L,1.0000002*L),
            units         = 'A',
            atoms         = self.system_info.dimer_atoms,
            basis         = [[0,0,0],[bond_length,0,0]],
            basis_vectors = identity(3),
            net_charge    = 0,
            net_spin      = 'low',
            kgrid         = (1,1,1),
            kshift        = (0,0,0)
            )
        s = system.structure
        s.slide(s.axes.sum(0)/2-s.pos.mean(0))
        return system
    #end def make_system


    def get_result_sim(self,v,name):
        key = v.functional,v.pp[0],v.pp[1]
        sim = ValidationStage.get_result_sim(self,(name,key))
        return sim
    #end def get_result_sim


    def get_results(self,quantity=None):
        res = obj()
        sinfo = self.system_info
        for functional,pplist in sinfo.func_pp.iteritems():
            for pp in pplist:
                key = (functional,pp[0],pp[1])
                v = self.saved_inputs[key]
                sims = self.sims[key]
                res[key] = self.get_result(v,sims,quantity)
            #end for
        #end for
        return res
    #end def get_results


    def insert_pop_tau_dmc(self,v):
        calcs = deepcopy(v.dmc_calcs)
        found_vmc = False
        found_dmc = False
        for calc in calcs:
            if isinstance(calc,vmc):
                if 'samplesperthread' in calc:
                    del calc.samplesperthread
                #end if
                calc.samples = v.population
                found_vmc = True
            elif isinstance(calc,dmc):
                dmc_calc = calc
                found_dmc = True
            #end if
        #end for
        if found_dmc:
            stepfac = calc.timestep/v.timestep
            dmc_calc.warmupsteps = int(round(dmc_calc.warmupsteps*stepfac))
            dmc_calc.steps = int(round(dmc_calc.steps*stepfac))
            dmc_calc.timestep = v.timestep
        #end if
        return calcs
    #end def insert_pop_tau_dmc


    def make_sims(self,v):
        self.not_implemented()
    #end def make_sims


    def get_result(self,v,sims,name):
        self.not_implemented()
    #end def get_result
        
#end class DimerValidationStage




class DimerBondLengthScan(DimerValidationStage):

    stage_dependencies = set(['bond_lengths','vmc_calcs','opt_calcs','dmc_calcs',
                              'dftjob','p2cjob','p2qjob','optjob','vmcjob','dmcjob',
                              'L','Ecut','population','timestep','J1_rcut','pade_b'])

    stage_type = 'scan'
    xlabel     = 'Bond length'
    xunit      = 'A'

    tag_sets = obj(
        Edft = obj(
            data_type = 'scalar',
            title     = 'DFT energy vs bond length',
            ylabel    = 'DFT energy',
            yunit     = 'Ry'
            ),
        Ekin_pw = obj(
            data_type = 'scalar',
            title     = 'PW kinetic energy vs bond length',
            ylabel    = 'PW kinetic energy',
            yunit     = 'Ha'
            ),
        Ekin_vmc = obj(
            data_type = 'stat',
            title     = 'VMC kinetic energy vs bond length',
            ylabel    = 'VMC kinetic energy',
            yunit     = 'Ha'
            ),
        Evmc_noJ = obj(
            data_type = 'stat',
            title     = 'Jastrowless VMC energy vs bond length',
            ylabel    = 'Jastrowless VMC energy',
            yunit     = 'Ha'
            ),
        Evmc = obj(
            data_type = 'stat',
            title     = 'VMC energy vs bond length',
            ylabel    = 'VMC energy',
            yunit     = 'Ha'
            ),
        Edmc = obj(
            data_type = 'stat',
            title     = 'DMC energy vs bond length',
            ylabel    = 'DMC energy',
            yunit     = 'Ha'
            )
        )


    def make_sims(self,v):
        sims = obj()
        dftpp,qmcpp = sort_pseudos(v.pseudos)
        rcut  = v.J1_rcut
        bu,bd = v.pade_b
        dmc_calcs = self.insert_pop_tau_dmc(v)
        for bond_length in v.bond_lengths:
            bpath = os.path.join(v.path,'bond_length_{0}'.format(bond_length))
            dftpath    = os.path.join(bpath,'dft')
            vmcnoJpath = os.path.join(bpath,'vmc_noJ')
            optJ12path = os.path.join(bpath,'opt_J12')
            optJ3path  = os.path.join(bpath,'opt_J3')
            vmcpath    = os.path.join(bpath,'vmc')
            dmcpath    = os.path.join(bpath,'dmc')

            dimer = self.make_system(v,bond_length)

            scf = generate_pwscf(
                identifier   = 'scf',
                path         = dftpath,
                job          = v.dftjob,
                input_type   = 'scf',
                input_dft    = v.functional,
                ecut         = v.Ecut,
                conv_thr     = 1e-8,
                mixing_beta  = .7,
                nosym        = True,
                pseudos      = dftpp,
                system       = dimer,
                kgrid        = (1,1,1),
                kshift       = (0,0,0),
                use_folded   = False
                )
            p2q = generate_pw2qmcpack(
                identifier   = 'p2q',
                path         = dftpath,
                job          = v.p2qjob,
                write_psir   = False
                )
            p2c = generate_pw2casino(
                identifier   = 'p2c',
                path         = dftpath,
                job          = v.p2cjob
                )
            vmc_noJ = generate_qmcpack(
                identifier   = 'vmc',
                path         = vmcnoJpath,
                job          = v.vmcjob,
                input_type   = 'basic',
                system       = dimer,
                bconds       = 'nnn',
                pseudos      = qmcpp,
                jastrows     = [],
                corrections  = [],
                calculations = v.vmc_calcs
                )
            optJ12 = generate_qmcpack(
                identifier  = 'opt',
                path        = optJ12path,
                job         = v.optjob,
                input_type  = 'opt_jastrow',
                system      = dimer,
                bconds      = 'nnn',
                pseudos     = qmcpp,
                jastrows    = [
                    generate_jastrow('J1','bspline',8,rcut,system=dimer),
                    generate_jastrow('J2','pade',bu,bd,system=dimer)
                    ],
                corrections = [],
                opt_calcs   = v.opt_calcs
                )
            optJ3 = generate_qmcpack(
                identifier  = 'opt',
                path        = optJ3path,
                job         = v.optjob,
                input_type  = 'opt_jastrow',
                system      = dimer,
                bconds      = 'nnn',
                pseudos     = qmcpp,
                jastrows    = [
                    generate_jastrow('J3','polynomial',4,4,5.0,system=dimer)
                    ],
                corrections = [],
                opt_calcs   = v.opt_calcs
                )
            vmc = generate_qmcpack(
                identifier   = 'vmc',
                path         = vmcpath,
                job          = v.vmcjob,
                input_type   = 'basic',
                system       = dimer,
                bconds       = 'nnn',
                pseudos      = qmcpp,
                corrections  = [],
                calculations = v.vmc_calcs
                )
            qmc = generate_qmcpack(
                identifier   = 'dmc',
                path         = dmcpath,
                job          = v.dmcjob,
                input_type   = 'basic',
                system       = dimer,
                bconds       = 'nnn',
                pseudos      = qmcpp,
                corrections  = [],
                calculations = dmc_calcs
                )
            p2q.depends(    scf,'orbitals')
            p2c.depends(    scf,'orbitals')
            vmc_noJ.depends(scf,'orbitals')
            optJ12.depends( scf,'orbitals')
            optJ3.depends( (scf,'orbitals'),(optJ12,'jastrow'))
            vmc.depends(   (scf,'orbitals'),(optJ3, 'jastrow'))
            qmc.depends(   (scf,'orbitals'),(optJ3, 'jastrow'))
            
            sims[bond_length,'scf'    ] = scf
            sims[bond_length,'p2q'    ] = p2q
            sims[bond_length,'p2c'    ] = p2c
            sims[bond_length,'vmc_noJ'] = vmc_noJ
            sims[bond_length,'optJ12' ] = optJ12
            sims[bond_length,'optJ3'  ] = optJ3
            sims[bond_length,'vmc'    ] = vmc
            sims[bond_length,'dmc'    ] = qmc
        #end for
        return sims
    #end def make_sims


    def get_result(self,v,sims,name='Edmc'):
        if not name in self.tag_sets:
            self.error('requested unrecognized result\n  you requested: {0}\n  valid options are: {1}'.format(name,self.tag_sets.keys()))
        #end if
        tags = self.tag_sets[name]
        self.transfer_from(tags)
        bond_lengths = []
        energies = []
        for bond_length in v.bond_lengths:
            if name=='Edft':
                a = sims[bond_length,'scf'].load_analyzer_image()
                energy = a.E
            elif name=='Ekin_pw':
                a = sims[bond_length,'p2c'].load_analyzer_image()
                energy = a.energies.Kinetic
            elif name=='Ekin_vmc':
                a = sims[bond_length,'vmc_noJ'].load_analyzer_image()
                ke = a.vmc.last().scalar.Kinetic
                energy = ke.mean,ke.error
            elif name=='Evmc_noJ':
                a = sims[bond_length,'vmc_noJ'].load_analyzer_image()
                le = a.vmc.last().scalar.LocalEnergy
                energy = le.mean,le.error
            elif name=='Evmc':
                a = sims[bond_length,'vmc'].load_analyzer_image()
                le = a.vmc.last().scalar.LocalEnergy
                energy = le.mean,le.error
            elif name=='Edmc':
                a = sims[bond_length,'dmc'].load_analyzer_image()
                le = a.dmc.last().dmc.LocalEnergy
                energy = le.mean,le.error
            #end if
            bond_lengths.append(bond_length)
            energies.append(energy)
        #end for
        return bond_lengths,energies
    #end def get_result
#end class DimerBondLengthScan




class ValidateDimerPP(ValidationProcess):

    allowed_stage_results = set(
        ['bond_lengths','vmc_calcs','opt_calcs','dmc_calcs',
         'dftjob','optjob','vmcjob','dmcjob','L','Ecut',
         'population','timestep','J1_rcut','pade_b']
        )

    def initialize(self,**func_pp):
        if 'name' in func_pp:
            self.name = func_pp['name']
            del func_pp['name']
        #end if
        self.check_funcs_pseudos(**func_pp)

        self.add_stages(
            bond_length_scan = DimerBondLengthScan()
            )

        self.set_stage_order(
            ['bond_length_scan']
            )

        info = obj(
            ppmap = self.ppmap,
            func_pp = self.func_pp,
            functionals = self.functionals,
            dimer = self.name,
            dimer_atoms = self.read_dimer_atoms()
            )

        for stage in self.stages:
            stage.system_info = info.copy()
        #end for
    #end def initialize


    def read_dimer_atoms(self):
        dimer = self.name
        atoms = []
        i=0
        while i<len(dimer):
            if dimer[i].isupper():
                if i+1<len(dimer) and dimer[i+1].islower():
                    atoms.append(dimer[i:i+2])
                    i+=2
                else:
                    atoms.append(dimer[i:i+1])
                    i+=1
                #end if
            else:
                self.error('cannot read dimer string: '+dimer)
            #end if
        #end for
        if len(atoms)!=2:
            self.error('dimer string {0} contains {1} atoms ({2}), but should contain exactly 2'.format(dimer,len(atoms),atoms))
        #end if
        return tuple(atoms)
    #end def read_dimer_atoms

#end class ValidateDimerPP




class ValidateDimerPPs(ValidationProcesses):
    def initialize(self):
        dimers = list(self.processes.keys())
        pseudos = obj()
        for dimer,process in self.processes.iteritems():
            pseudos[dimer] = process.pseudo_list()
        #end for
        self.dimers = dimers
        self.pseudos = pseudos
    #end def initialize


    def set_stage_inputs(self,**kwargs):
        self.set_stage_results(**kwargs)
    #end def set_stage_inputs


    def set_stage_results(self,**kwargs):
        results = obj()
        for name,result in kwargs.iteritems():
            results[name] = self.expand_stage_result(result)
        #end for

        for dimer,process in self.processes.iteritems():
            res = obj()
            for name,rcoll in results.iteritems():
                res[name] = deepcopy(rcoll[dimer])
            #end for
            process.check_stage_results(res.keys())
            process.set_stage_results(**res)
        #end for
    #end def set_stage_results


    def expand_stage_result(self,rin):
        r = obj()
        single_types = (str,int,float,tuple,list,dict,type(None),Job)
        wrong_type = False
        wrong_class = None
        if isinstance(rin,single_types):
            for dimer in self.dimers:
                dres = obj()
                for pp in self.pseudos[dimer]:
                    dres[pp] = deepcopy(rin)
                #end for
                r[dimer] = dres
            #end for
        elif isinstance(rin,obj):
            dimers = set(rin.keys())
            if len(set(self.dimers)-dimers)>0:
                self.error('cannot expand stage result\n  information is required for all requested dimers\n  dimers requested: {0}\n  information provided:\n{1}'.format(self.dimers,rin))
            #end if
            for dimer,dinfo in rin.iteritems():
                if dimer in self.dimers:
                    dres = obj()
                    if isinstance(dinfo,dict):
                        pps = set(dinfo.keys())
                        if len(set(self.pseudos[dimer])-pps)>0:
                            self.error('cannot expand stage result\n  information is required for all requested pseudopotentials\n  pseudopotentials requested: {0}\n  information provided:\n{1}'.format(self.pseudos[dimer],rin))
                        #end if
                        for pp,ppinfo in dinfo.iteritems():
                            if isinstance(ppinfo,single_types):
                                dres[pp] = deepcopy(ppinfo)
                            else:
                                wrong_type = True
                                wrong_class = ppinfo.__class__.__name__
                            #end if
                        #end for
                    elif isinstance(dinfo,single_types) or isinstance(dinfo,obj):
                        for pp in self.pseudos[dimer]:
                            dres[pp] = deepcopy(dinfo)
                        #end for
                    else:
                        wrong_type  = True
                        wrong_class = dinfo.__class__.__name__ 
                    #end if
                    r[dimer] = dres
                #end if
            #end for
        else:
            wrong_type = True
            wrong_class =rin.__class__.__name__
        #end if
        if wrong_type:
            self.error('cannot expand stage result\n  received type: {0}\n  allowed_types: {1}'.format(wrong_class,single_types))
        #end if
        return r
    #end def expand_stage_result

#end class ValidateDimerPPs
ValidateDimerPPs.set_process_type(ValidateDimerPP)




class ValidatePPs(ValidatePPBase):
    test_order = ['atoms','dimers']

    def __init__(self,settings=None,atoms=None,dimers=None):
        if not isinstance(settings,Settings):
            self.error('input variable settings must be the settings object from Nexus')
        #end if
        ValidatePPBase.settings = settings
        if atoms!=None:
            if not isinstance(atoms,(obj,dict)):
                self.error('input variable atoms must be of type dict or obj\n  you provided: '+atoms.__class__.__name__)
            #end if
            self.atoms = ValidateAtomPPs(**atoms)
        #end if
        if dimers!=None:
            if not isinstance(dimers,(obj,dict)):
                self.error('input variable dimers must be of type dict or obj\n  you provided: '+dimers.__class__.__name__)
            #end if
            self.dimers = ValidateDimerPPs(**dimers)
        #end if
        tests = obj()
        for name,test in self.iteritems():
            if test!=None:
                test.set_name(name)
                tests[name] = test
            #end if
        #end for
        self.tests = tests
    #end def __init__

            
    def sim_list(self):
        sims = []
        for test in self.tests:
            sims.extend(test.sim_list())
        #end for
        return sims 
    #end def sim_list


    def make_sims(self):
        print 'make_sims',self.__class__.__name__
        for name,test in self.tests.iteritems():
            test.make_sims(name)
        #end for
        return self.sim_list()
    #end def make_sims


    def status(self,*args):
        pad='  '
        n=0
        print 
        print n*pad+'ValidatePPs status'
        for name in self.test_order:
            if name in self.tests:
                self.tests[name].status(pad,n+1,*args)
            #end if
        #end for
        print n*pad+'end ValidatePPs status'
        print
        if not 'continue' in args:
            exit()
        #end if
    #end def status
#end class ValidatePPs




def validate_qmc_pp(**kwargs):
    return ValidateQMCPPs(**kwargs)
#end def validate_qmc_pp




if __name__=='__main__':

    settings(
        pseudo_dir    = './pseudopotentials',
        status_only   = 0,
        generate_only = 1,
        sleep         = .3,
        machine       = 'node32'
        )
 


    v = ValidatePPs(
        settings = settings,
        atoms = obj(
            #C = obj(
            #    q   = [0,1,2],
            #    lda = ['C.BFD'],
            #    pbe = ['C.BFD']
            #    ),
            O = obj(
                q   = [0,1],
                lda = ['O.6_lda']
                ),
            ),
        dimers = obj(
            CO = obj(
                lda = [('C.BFD','O.6_lda')],
                pbe = [('C.BFD','O.6_lda')]
                ),
            CuO = obj(
                lda = [('Cu.17_lda'    ,'O.6_lda'),
                       ('Cu.19_lda_opt','O.6_lda')]
                )
            )
        )



    v.atoms.set_stage_inputs(
        ae_occupations = obj(     #varies w/ a,q
            C  = [[('He2s2p(-1,0)','He2s'),  # 0
                   ('He2s2p(-1,1)','He2s'),
                   ('He2s2p( 0,1)','He2s'),
                   ('He2s2p(-1)'  ,'He2s2p(-1)'),
                   ('He2s2p( 0)'  ,'He2s2p( 0)'),
                   ('He2s2p( 1)'  ,'He2s2p( 1)'),
                   ('He2s2p(-1)'  ,'He2s2p( 0)'),
                   ('He2s2p(-1)'  ,'He2s2p( 1)')
                   ],
                  [('He2s2p(-1)'  ,'He2s'),  # 1
                   ('He2s2p( 0)'  ,'He2s'),
                   ('He2s2p( 1)'  ,'He2s'),
                   ('He2s2p(-1,0)','He'),
                   ('He2s2p(-1)'  ,'He2p(-1)')
                   ],
                  [('He2s'      ,'He2s'),    # 2
                   ('He2s2p(-1)','He'),
                   ('He2s'      ,'He2p(-1)')
                   ]
                  ],
            O  = [#[('He2s2p','He2s2p'),      # -2
                  # ('He2s2p','He2s2p(-1,0)3s'),
                  # ],
                  #[('He2s2p','He2s2p(-1,0)'),# -1
                  # ('He2s2p','He2s2p(-1,1)'),
                  # ('He2s2p','He2s2p( 0,1)'),
                  # ('He2s2p','He2p')
                  # ],
                  [('He2s2p','He2s2p(-1)'),  # 0
                   ('He2s2p','He2s2p( 0)'),
                   ('He2s2p','He2s2p( 1)'),
                   ('He2s2p(-1,0)','He2s2p(-1,1)'),
                   ('He2s2p(-1,1)','He2s2p(-1,0)'),
                   ('He2s2p(-1,0)','He2s2p( 0,1)'),
                   ('He2s2p( 0,1)','He2s2p(-1,0)'),
                   ('He2s2p(-1,1)','He2s2p( 0,1)'),
                   ('He2s2p( 0,1)','He2s2p(-1,1)')
                   ],
                  [('He2s2p','He2s'),        # 1
                   ('He2s2p(-1,1)','He2s2p(0)'   ),
                   ('He2s2p(-1,0)','He2s2p(1)'   ),
                   ('He2s2p(-1)'  ,'He2s2p(0,1)' ),
                   ('He2s2p(-1,0)','He2s2p(-1)'  ),
                   ('He2s2p(-1,1)','He2s2p(-1)'  ),
                   ('He2s2p(-1)'  ,'He2s2p(-1,0)'),
                   ('He2s2p(-1)'  ,'He2s2p(-1,1)'),
                   ]
                  #[('He2s2p(-1,0)','He2s'),  # 2
                  # ('He2s2p(-1,1)','He2s'),
                  # ('He2s2p( 0,1)','He2s'),
                  # ('He2s2p(-1)'  ,'He2s2p(-1)'),
                  # ('He2s2p( 0)'  ,'He2s2p( 0)'),
                  # ('He2s2p( 1)'  ,'He2s2p( 1)'),
                  # ('He2s2p(-1)'  ,'He2s2p( 0)'),
                  # ('He2s2p(-1)'  ,'He2s2p( 1)')
                  # ]
                  ]
            ),
        #pp_Ls          = [10,15,20],             #same for all a,q,pp
        #pp_Ecuts       = [100,150,200,250,300],  #same for all a,q,pp
        pp_spins       = obj(
            C = [[0,2],[1,3],[0,2]],
            O = [[0,2],[1,3]]
            ),
        pp_J1_rcuts    = [3,4,5,6,7],            #same for all a,q,pp
        ae_populations = [500,1000,2000],        #same for all a,q
        pp_populations = [500,1000,2000],        #same for all a,q,pp
        ae_timesteps   = [.008,.004,.002,.001],  #same for all a,q
        pp_timesteps   = [.08,.04,.02,.01,.005], #same for all a,q,pp
        ae_opt_calcs   = [
            loop(max=4,
                 qmc = linear(
                    energy               = 0.0,
                    unreweightedvariance = 1.0,
                    reweightedvariance   = 0.0,
                    walkers              = 1,
                    warmupsteps          = 200,
                    blocks               = 500,
                    timestep             = 0.1,
                    usedrift             = True,
                    samples              = 16000,
                    stepsbetweensamples  = 100,
                    minwalkers           = 0.0,
                    bigchange            = 15.0,
                    alloweddifference    = 1e-4
                    )
                 ),
            loop(max=4,
                 qmc = linear(
                    energy               = 0.5,
                    unreweightedvariance = 0.0,
                    reweightedvariance   = 0.5,
                    walkers              = 1,
                    warmupsteps          = 200,
                    blocks               = 500,
                    timestep             = 0.1,
                    usedrift             = True,
                    samples              = 32000,
                    stepsbetweensamples  = 100,
                    minwalkers           = 0.0,
                    bigchange            = 15.0,
                    alloweddifference    = 1e-4
                    )
                 )
            ],
        ae_vmc_calcs   = [
            vmc(
                walkers     = 1,   
                warmupsteps = 50,
                blocks      = 1000,
                steps       = 10,
                substeps    = 10,
                timestep    = 0.1,
                usedrift    = True
                )
            ],
        ae_dmc_calcs   = [
            vmc(
                walkers     = 1,   
                warmupsteps = 20,
                blocks      = 100,
                steps       = 10,
                substeps    = 10,
                timestep    = 0.1,
                usedrift    = True,
                samples     = 1000
                ),
            dmc(
                warmupsteps = 24,
                blocks      = 1000,
                steps       = 12,
                timestep    = 0.004
                ) 
            ],
        pp_opt_calcs   = [
            loop(max=4,
                 qmc = linear(
                    energy               = 0.0,
                    unreweightedvariance = 1.0,
                    reweightedvariance   = 0.0,
                    walkers              = 1,
                    warmupsteps          =  50,
                    blocks               = 500,
                    timestep             = 0.4,
                    usedrift             = True,
                    samples              = 16000,
                    stepsbetweensamples  = 10,
                    minmethod            = 'quartic',
                    minwalkers           = 0.5,
                    bigchange            = 2.0,
                    alloweddifference    = 1e-4,
                    maxweight            = 1e9,
                    beta                 = 0.025,
                    exp0                 = -16,
                    stepsize             = 0.2,
                    stabilizerscale      = 1.0,
                    nstabilizers         = 3,
                    nonlocalpp           = True,
                    usebuffer            = True
                    )
                 ),
            loop(max=4,
                 qmc = linear(
                    energy               = 0.5,
                    unreweightedvariance = 0.0,
                    reweightedvariance   = 0.5,
                    walkers              = 1,
                    warmupsteps          =  50,
                    blocks               = 500,
                    timestep             = 0.4,
                    usedrift             = True,
                    samples              = 32000,
                    stepsbetweensamples  = 10,
                    minmethod            = 'quartic',
                    minwalkers           = 0.5,
                    bigchange            = 2.0,
                    alloweddifference    = 1e-4,
                    maxweight            = 1e9,
                    beta                 = 0.025,
                    exp0                 = -16,
                    stepsize             = 0.2,
                    stabilizerscale      = 1.0,
                    nstabilizers         = 3,
                    nonlocalpp           = True,
                    usebuffer            = True
                    )
                 )
            ],
        pp_vmc_calcs   = [
            vmc(
                walkers     = 1,   
                warmupsteps = 50,
                blocks      = 1000,
                steps       = 20,
                substeps    =  3,
                timestep    = 0.4,
                usedrift    = True
                )
            ],
        pp_dmc_calcs   = [
            vmc(
                walkers     = 1,   
                warmupsteps = 50,
                blocks      = 1000,
                steps       = 20,
                substeps    =  3,
                timestep    = 0.4,
                usedrift    = True,
                samples     = 1000,
                ),
            dmc(
                warmupsteps = 12,
                blocks      = 1000,
                steps       = 4,
                timestep    = 0.04,
                nonlocalmoves = True
                ) 
            ] 
        )


    v.atoms.set_stage_results(
        ae_hfjob      = Job(cores=1),
        ae_optjob     = Job(cores=16),
        ae_vmcjob     = Job(cores=16),
        ae_dmcjob     = Job(cores=16),
        pp_dftjob     = Job(cores=16), 
        pp_p2cjob     = Job(cores=1), 
        pp_p2qjob     = Job(cores=1), 
        pp_optjob     = Job(cores=16),
        pp_vmcjob     = Job(cores=16),
        pp_dmcjob     = Job(cores=16),
        pp_Ecut0      = 150,   
        pp_spin0      = obj(
            C = 2,
            O = 2
            ),
        #ae_occupation = obj(   
        #    C         = [('He2s2p(-1,0)','He2s'        ),
        #                 ('He2s2p(-1)'  ,'He2s'        ),
        #                 ('He2s'        ,'He2s'        )],
        #    O         = [('He2s2p'      ,'He2s2p(-1,0)'),
        #                 ('He2s2p'      ,'He2s2p(-1)'  ),
        #                 ('He2s2p'      ,'He2s'        ),
        #                 ('He2s2p(-1,0)','He2s'        )]
        #    ),
        pp_L          = 10,    
        pp_Ecut = obj(
            C         = 150,
            O         = 150
            ),
        #pp_spin = obj(
        #    C = 0,
        #    O = 2
        #    ),
        #ae_orbitals   = 'finished',
        #pp_orbitals   = 'finished',
        #ae_pade_b     = (4.,4.),
        #pp_pade_b     = (4.,4.),
        #pp_J1_rcut = obj( 
        #    C         = {'C.BFD' : [5, 4, 3]},
        #    O         = 4
        #    ),
        #ae_jastrow    = 'finished',
        #pp_jastrow    = 'finished',
        #ae_population = obj(
        #    C         = 1000,  
        #    O         = 2000
        #    ),
        #pp_population = obj(
        #    C         = [2000,3000,4000], 
        #    O         = [1500,2500,3500,4500]
        #    ),
        #ae_timestep   = 1e-3, 
        #pp_timestep = {     
        #    'C.BFD'   : .01,
        #    'O.6_lda' : .02
        #    }
        )



    v.dimers.set_stage_inputs(
        #bond_lengths = obj(
        #    CO  = [2,3,4,5,6],
        #    #CuO = [2.5,3.5,4.5,5.5] 
        #    ),
        dftjob = Job(cores=16),
        optjob = Job(cores=16),
        vmcjob = Job(cores=16),
        dmcjob = Job(cores=16),
        L = obj(
            CO  = 20,
            CuO = 20,
            ),
        Ecut = obj(
            CO  = 200,
            CuO = {
                ('Cu.17_lda'    ,'O.6_lda') : 200,
                ('Cu.19_lda_opt','O.6_lda') : 300
                }
            ),
        population = obj(
            CO  = 1000,
            CuO = 1500
            ),
        timestep = obj(
            CO  = .02,
            CuO = {
                ('Cu.17_lda'    ,'O.6_lda') : .01,
                ('Cu.19_lda_opt','O.6_lda') : .005
                }
            ),
        J1_rcut = obj(
            CO  = obj(C =4,O=4.5),
            CuO = obj(Cu=5,O=4.5)
            ),
        pade_b = obj(
            CO  = (3,4),
            CuO = (5,6)
            ),
        opt_calcs   = [
            loop(max=4,
                 qmc = linear(
                    energy               = 0.0,
                    unreweightedvariance = 1.0,
                    reweightedvariance   = 0.0,
                    walkers              = 1,
                    warmupsteps          =  50,
                    blocks               = 500,
                    timestep             = 0.4,
                    usedrift             = True,
                    samples              = 16000,
                    stepsbetweensamples  = 10,
                    minmethod            = 'quartic',
                    minwalkers           = 0.5,
                    bigchange            = 2.0,
                    alloweddifference    = 1e-4,
                    maxweight            = 1e9,
                    beta                 = 0.025,
                    exp0                 = -16,
                    stepsize             = 0.2,
                    stabilizerscale      = 1.0,
                    nstabilizers         = 3,
                    nonlocalpp           = True,
                    usebuffer            = True
                    )
                 ),
            loop(max=4,
                 qmc = linear(
                    energy               = 0.5,
                    unreweightedvariance = 0.0,
                    reweightedvariance   = 0.5,
                    walkers              = 1,
                    warmupsteps          =  50,
                    blocks               = 500,
                    timestep             = 0.4,
                    usedrift             = True,
                    samples              = 32000,
                    stepsbetweensamples  = 10,
                    minmethod            = 'quartic',
                    minwalkers           = 0.5,
                    bigchange            = 2.0,
                    alloweddifference    = 1e-4,
                    maxweight            = 1e9,
                    beta                 = 0.025,
                    exp0                 = -16,
                    stepsize             = 0.2,
                    stabilizerscale      = 1.0,
                    nstabilizers         = 3,
                    nonlocalpp           = True,
                    usebuffer            = True
                    )
                 )
            ],
        vmc_calcs   = [
            vmc(
                walkers     = 1,   
                warmupsteps = 50,
                blocks      = 1000,
                steps       = 20,
                substeps    =  3,
                timestep    = 0.4,
                usedrift    = True
                )
            ],
        dmc_calcs   = [
            vmc(
                walkers     = 1,   
                warmupsteps = 50,
                blocks      = 1000,
                steps       = 20,
                substeps    =  3,
                timestep    = 0.4,
                usedrift    = True,
                samples     = 1000,
                ),
            dmc(
                warmupsteps = 12,
                blocks      = 1000,
                steps       = 4,
                timestep    = 0.04,
                nonlocalmoves = True
                ) 
            ]
        )


    # make the simulations
    sims = v.make_sims()



    #v.status()


    pm = ProjectManager()

    pm.add_simulations(sims)

    pm.run_project()


    #v.atoms.show_results(
    #    stages = 'pp_box_scan',
    #    mode   = 'print',
    #    units  = 'eV',
    #    xunits = 'B'
    #    )


    v.atoms.show_results(
        stages = 'ae_occupation_scan',
        mode   = 'print'
        )










    #v = ValidatePPs(
    #    settings = settings,
    #    atoms = obj(
    #        Cu = obj(
    #            ref = (0,'lda','Cu.pz'),
    #            q   = [0,1,2],
    #            lda = ['Cu.17_lda','Cu.17_lda_f','Cu.19_lda_opt'],
    #            pbe = ['Cu.17_lda']
    #            ),
    #        Zn = obj(
    #            q   = [0,1,2],
    #            lda = ['Zn.lda.pp5','Zn.lda.pp6'],
    #            pbe = ['Zn.lda.pp6']
    #            ),
    #        #Co = obj(
    #        #    q = [0,1,2],
    #        #    pbe = ['Co.pbe'],
    #        #    hse = ['Co.pbe','Co.pbe2']
    #        #    )
    #        ),
    #    #dimers = obj(
    #    #    CuO = obj(
    #    #        lda = ['Cu.pz','O.pz'],
    #    #        pbe = ['Cu.pbe','O.pbe']
    #    #        ),
    #    #    ZnO = obj(
    #    #        lda = ['Zn.pz','Zn.pz']
    #    #        ),
    #    #    CoO = obj(
    #    #        lda = ['Co.pz','O.pz'],
    #    #        pbe = ['Co.pbe','O.pbe']
    #    #        )
    #    #    )
    #    )
    #
    #v.atoms.set_stage_inputs(
    #    hfjob  = Job(cores=1),
    #    dftjob = Job(cores=16), 
    #    optjob = Job(cores=16),
    #    vmcjob = Job(cores=16),
    #    dmcjob = Job(cores=16),
    #    ae_occupations = obj(     #varies w/ a,q
    #        Cu = [('abc1','def')],
    #        Zn = [('abc2','def')],
    #        Co = [('abc3','def')]
    #        ),
    #    ae_opt_calcs       = [],    #same for all a,q
    #    ae_vmc_calcs       = [],    #same for all a,q
    #    ae_populations     = [],    #same for all a,q
    #    ae_timesteps       = [],    #same for all a,q
    #    ae_dmc_calcs       = [],    #same for all a,q
    #    pp_Ls              = [],    #same for all a,q,pp
    #    pp_Ecuts           = [],    #same for all a,q,pp
    #    pp_opt_calcs       = [],    #same for all a,q,pp
    #    pp_J1_rcuts        = [],    #same for all a,q,pp
    #    pp_vmc_calcs       = [],    #same for all a,q,pp
    #    pp_dmc_populations = [],    #same for all a,q,pp
    #    pp_dmc_timesteps   = [],    #same for all a,q,pp
    #    pp_dmc_calcs       = []     #same for all a,q,pp
    #    )
    #
    ##v.atoms.set_stage_results(
    ##    ae_occupation = obj(      #varies w/ a,q
    ##        Cu = [],
    ##        Zn = [],
    ##        Co = []
    ##        ),
    ##    ae_population = 1000,     #can vary or be same
    ##    ae_timestep   = .0001,    #can vary or be same
    ##    pp_L          = 20,       #can vary or be same
    ##    pp_Ecut       = 300,      #will vary w/ a,pp
    ##    pp_J1_rcut    = obj(      #will vary w/ a,q,pp
    ##        Cu = {'Cu.pz' : [3.4, 4.5, 3.6],
    ##              'Cu.pbe': [3.5, 3.4, 4.2]},
    ##        Co = {}
    ##        ),
    ##    pp_population = 2000,     #can vary or be same
    ##    pp_timestep   = .01       #can vary or be same
    ##    )
    #
    #v.atoms.set_stage_results(
    #    ae_occupation = obj(      #varies w/ a,q
    #        Cu = [('up','down'),('up','down'),('up','down')],
    #        Zn = [('up','down'),('up','down'),('up','down')],
    #        Co = []
    #        ),
    #    ae_orbitals = 'finished',
    #    ae_jastrow  = 'finished',
    #    ae_population = obj(
    #        Cu = 1000,     #can vary or be same
    #        Zn = 2000
    #        ),
    #    ae_timestep   = [1e-4,2e-4,3e-4],    #can vary or be same
    #    pp_Ecut0      = 200,   #same for all a,q,pp
    #    pp_L          = 20,       #can vary or be same
    #    pp_Ecut       = obj(
    #        Cu = 150,
    #        Zn = 200
    #        ),      #will vary w/ a,pp
    #    pp_orbitals = 'finished',
    #    pp_jastrow  = 'finished',
    #    pp_J1_rcut    = obj(      #will vary w/ a,q,pp
    #        Cu = {'Cu.17_lda'    : [3.4, 4.5, 3.6],
    #              'Cu.17_lda_f'  : [3.5, 3.4, 4.2],
    #              'Cu.19_lda_opt': 5.0  },
    #        Zn = 2.7,
    #        Co = {}
    #        ),
    #    pp_population = obj(
    #        Cu = [2000,3000,4000],     #can vary or be same
    #        Zn = [2500,3500,4500]
    #        ),
    #    pp_timestep   = {       #can vary or be same
    #        'Cu.17_lda'    :.01,
    #        'Cu.17_lda_f'  :.02,
    #        'Cu.19_lda_opt':.001,
    #        'Zn.lda.pp5':.005,
    #        'Zn.lda.pp6':.006
    #        }
    #    )
    #
#end if

