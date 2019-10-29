##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_method_analyzers.py                                       #
#    Analyzer classes at the level of QMC methods.  Instances        #
#    contain all data outputted by VMC, Opt, DMC, etc. sub-runs      #
#    carried out by QMCPACK.                                         #
#                                                                    #
#  Content summary:                                                  #
#    MethodAnalyzer                                                  #
#      Base class for specific method analyzers.                     #
#      Derived classes are OptAnalyzer, VmcAnalyzer, DmcAnalyzer     #
#                                                                    #
#====================================================================#


import os
import re
from generic import obj
from hdfreader import HDFreader
from qmcpack_analyzer_base import Checks,QAanalyzer,QAdata,QAHDFdata
from qmcpack_property_analyzers import WavefunctionAnalyzer
from qmcpack_quantity_analyzers import HDFAnalyzer
from debug import *


class MethodAnalyzer(QAanalyzer):
    def __init__(self,series=None,calc=None,input=None,nindent=0):
        QAanalyzer.__init__(self,nindent=nindent)
        if series!=None and calc!=None and input!=None:
            self.init_sub_analyzers(series,calc,input)
        #end if
    #end def __init__


    def init_sub_analyzers(self,series,calc,input):
        request  = QAanalyzer.request
        run_info = QAanalyzer.run_info 
        
        source_path = run_info.source_path
        file_prefix = run_info.file_prefix+'.s'+str(series).zfill(3)
        method = calc.method

        files  = obj()
        outfiles = os.listdir(source_path)
        self.vlog('looking for file prefix: '+file_prefix,n=2)
        matched = False
        for file in outfiles:
            if file.startswith(file_prefix):
                local_match = True
                if file.endswith('scalar.dat'):
                    files.scalar = file
                elif file.endswith('stat.h5'):
                    files.stat   = file
                elif file.endswith('storeConfig.h5'):
                    files.storeconfig = file
                elif file.endswith('opt.xml'):
                    files.opt    = file
                elif file.endswith('dmc.dat') and method=='dmc':
                    files.dmc    = file
                elif '.traces.' in file:
                    if not 'traces' in files:
                        files.traces = []
                    #end if
                    files.traces.append(file)
                else:
                    local_match = False
                #end if
                matched = matched or local_match
                self.vlog('match found: '+file,n=3)
            #end if
        #end for
        complete = matched
        complete &= 'scalar' in files
        if 'linear' in method or method=='opt':
            complete &= 'opt' in files
        #end if
        equil = request.equilibration
        nblocks_exclude = -1
        if isinstance(equil,int):
            nblocks_exclude = equil
        elif isinstance(equil,(dict,obj)):
            if series in equil:
                nblocks_exclude = equil[series]
            #end if
        elif equil!=None:
            self.error('invalid input for equilibration which must be an int, dict, or obj\n  you provided: {0}\n  with type {1}'.format(equil,equil.__class__.__name__))
        #end if
        data_sources     = request.data_sources & set(files.keys())
        method_info = obj(
            method       = method,
            series       = series,
            file_prefix  = file_prefix,
            files        = files,
            data_sources = data_sources,
            method_input = calc.copy(),
            nblocks_exclude = nblocks_exclude,
            complete     = complete,
            )
        self.info.transfer_from(method_info)

        self.vlog('requested sources = '+str(list(request.data_sources)),n=2)
        self.vlog('files available   = '+str(files.keys()),n=2)
        self.vlog('available sources = '+str(list(data_sources)),n=2)

        if not matched:
            msg = 'no data files found\n  file prefix used for matching: {0}\n  checked all files in directory: {1}'.format(file_prefix,source_path)
            #self.error(msg,trace=False)
            #self.warn(msg)
            return
        #end if

        self.set_global_info()

        try:
            analyzers = self.capabilities.analyzers
            if 'scalar' in data_sources:
                filepath = os.path.join(source_path,files.scalar)
                self.scalars = analyzers.scalars_dat(filepath,equilibration='LocalEnergy',nindent=self.subindent())
            #end if
            if 'stat' in data_sources:
                #determine scalars and analyzer quantities
                analyzer_quantities = self.capabilities.analyzer_quantities
                analyzer_quants = obj()
                ignored_quantities = set()
                ham = input.get('hamiltonian')
                ham = ham.get_single('h0')
                ham_est  = ham.get('estimator')
                calc_est = calc.get('estimator')
                estimators = obj()
                if ham_est!=None:
                    estimators.transfer_from(ham_est)
                #end if
                if calc_est!=None:
                    estimators.transfer_from(calc_est)
                #end if
                for estname,est in estimators.iteritems():
                    if est==None:
                        self.error('estimators have not been read properly by QmcpackInput',trace=False)
                    #end if
                    has_type = 'type' in est
                    has_name = 'name' in est
                    if has_type and has_name:
                        type = est.type
                        name = est.name
                    elif has_name:
                        type = est.name
                        name = est.name
                    elif has_type:
                        type = est.type
                        name = est.type
                    else:
                        self.error('estimator '+estname+' has no type or name')
                    #end if
                    cname = self.condense_name(name)
                    ctype = self.condense_name(type)

                    if ctype=='density' and not has_name:
                        name = 'any'
                    #end if

                    if ctype in analyzer_quantities:
                        analyzer_quants[name] = self.condense_name(type)
                    #end if
                #end for
                not_scalars = set(analyzer_quants.keys())

                self.scalars_hdf = analyzers.scalars_hdf(not_scalars,nindent=self.subindent())

                analyzer_quantities = analyzer_quantities & request.quantities
                for name,type in analyzer_quants.iteritems():
                    if type in analyzer_quantities:
                        if type in analyzers:
                            qqa = analyzers[type](name,nindent=self.subindent())
                            qqa.init_sub_analyzers()
                            self[name] = qqa
                        else:
                            ignored_quantities.add(name)
                        #end if
                    #end if
                #end for
                self.info.ignored_quantities = ignored_quantities
            #end if
            if 'dmc' in data_sources:
                filepath = os.path.join(source_path,files.dmc)
                self.dmc = analyzers.dmc_dat(filepath,nindent=self.subindent())
            #end if
            if 'traces' in data_sources and 'traces' in files:
                self.traces = analyzers.traces(source_path,files.traces,nindent=self.subindent())
            #end if
        except:
            self.info.complete = False
        #end try

        self.unset_global_info()

        return
    #end def init_sub_analyzers


    def load_data_local(self):
        source_path = QAanalyzer.run_info.source_path
        data_sources = self.info.data_sources
        files  = self.info.files
        if 'stat' in data_sources:
            filepath = os.path.join(source_path,files.stat)
            hr = HDFreader(filepath)
            if not hr._success:
                self.warn('  hdf file seems to be corrupted, skipping contents:\n    '+filepath)
            #end if
            hdf = hr.obj
            self.data = QAHDFdata()
            self.data.transfer_from(hdf)
        #end if
        remove = []
        for name,value in self.iteritems():
            if isinstance(value,HDFAnalyzer):
                value.load_data_local(self.data)
                value.info.data_loaded = True
                if value.info.should_remove:
                    remove.append(name)
                #end if
            #end if
        #end for       
        for name in remove:
            del self[name]
        #end for
    #end def load_data_local
        

    
    def set_global_info(self):
        QAanalyzer.method_info = self.info
    #end def set_global_info

    def unset_global_info(self):
        QAanalyzer.method_info = None
    #end def unset_global_info


    def check_traces(self,pad=None):
        verbose = pad!=None
        method = self.info.method
        series = self.info.series
        if verbose:
            desc = 'method {0} series {1}'.format(method,series)
        #end if
        if 'traces' in self:
            check = {None:True,False:False,True:True}
            if verbose:
                self.log(pad+'Checking traces in '+desc)
            #end if
            scalars     = None
            scalars_hdf = None
            dmc         = None
            if 'scalars' in self and 'data' in self.scalars:
                scalars = self.scalars.data
            #end if
            if 'scalars_hdf' in self and 'data' in self.scalars_hdf:
                scalars_hdf = self.scalars_hdf.data
            #end if
            if 'dmc' in self and 'data' in self.dmc:
                dmc = self.dmc.data
            #end if
            checks = Checks('traces')
            checks.exclude(None)
            traces = self.traces
            traces.form_diagnostic_data()
            checks.psums   = traces.check_particle_sums()
            if method=='dmc':
                checks.dmc = traces.check_dmc(dmc)
            else:
                svalid,shvalid = traces.check_scalars(scalars,scalars_hdf)
                checks.scalars     = svalid
                checks.scalars_hdf = shvalid
            #end if
            valid = checks.valid()
            if verbose:
                checks.write(pad+'  ')
            #end if
            return valid
        else:
            if verbose:
                self.log(pad+'No traces in '+desc)
            #end if
            return None
        #end if
    #end def check_traces

#end class MethodAnalyzer




class OptAnalyzer(MethodAnalyzer):
    def init_sub_analyzers(self,series,calc,input):
        MethodAnalyzer.init_sub_analyzers(self,series,calc,input)

        source_path = QAanalyzer.run_info.source_path
        files = self.info.files

        if 'opt' in files:
            opt_out_xml = os.path.join(source_path,files.opt)
            self.wavefunction = WavefunctionAnalyzer(opt_out_xml)
        #end if
    #ed def init_sub_analyzers
#end class OptAnalyzer

class VmcAnalyzer(MethodAnalyzer):
    None
#end class OptAnalyzer

class DmcAnalyzer(MethodAnalyzer):
    None
#end class OptAnalyzer
