##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_analyzer.py                                               #
#    Supports data analysis for QMCPACK output.                      #
#                                                                    #
#  Content summary:                                                  #
#    QmcpackAnalyzer                                                 #
#      SimulationAnalyzer class for QMCPACK.                         #
#                                                                    #
#    QmcpackCapabilities                                             #
#      Class to pair QMCPACK output data with analyzer classes.      #
#                                                                    #
#    QmcpackAnalysisRequest                                          #
#      Offers detailed control over analysis.                        #
#      Serves as a record of requested analysis.                     #
#                                                                    #
#====================================================================#



from time import time

#python standard library imports
import os
import re
import sys
import traceback
from numpy import arange,array
#custom library imports
from generic import obj
from developer import unavailable
from xmlreader import XMLreader
from plotting import *
#QmcpackAnalyzer classes imports
from qmcpack_analyzer_base import QAobject,QAanalyzer,QAanalyzerCollection
from qmcpack_property_analyzers \
    import WavefunctionAnalyzer
from qmcpack_quantity_analyzers \
    import ScalarsDatAnalyzer,ScalarsHDFAnalyzer,DmcDatAnalyzer,\
    EnergyDensityAnalyzer,TracesAnalyzer,DensityMatricesAnalyzer,\
    SpinDensityAnalyzer,StructureFactorAnalyzer,DensityAnalyzer
from qmcpack_method_analyzers \
    import OptAnalyzer,VmcAnalyzer,DmcAnalyzer
from qmcpack_result_analyzers \
    import OptimizationAnalyzer,TimestepStudyAnalyzer
from simulation import SimulationAnalyzer,Simulation
from qmcpack_input import QmcpackInput
from debug import *

try:
    import h5py
    h5py_unavailable = False
except ImportError:
    h5py = unavailable('h5py')
    h5py_unavailable = True
#end try



class QmcpackAnalyzerCapabilities(QAobject):

    def __init__(self):

        self.methods=set(['opt','vmc','dmc','rmc'])
        self.data_sources = set(['scalar','stat','dmc','storeconfig','opt','traces'])
        self.scalars=set(['localenergy','localpotential','kinetic','elecelec','localecp','nonlocalecp','ionion','localenergy_sq','acceptratio','blockcpu','blockweight','mpc','kecorr'])
        self.fields=set(['energydensity','density','dm1b','spindensity','structurefactor'])

        hdf_data_sources = set(['stat','storeconfig','traces'])
        if h5py_unavailable:
            self.data_sources -= hdf_data_sources
        #end if

        self.analyzer_quantities = set(self.fields)

        self.analyzers = obj(
            scalars_dat     = ScalarsDatAnalyzer,
            scalars_hdf     = ScalarsHDFAnalyzer,
            dmc_dat         = DmcDatAnalyzer,
            traces          = TracesAnalyzer,
            energydensity   = EnergyDensityAnalyzer,
            dm1b            = DensityMatricesAnalyzer,
            spindensity     = SpinDensityAnalyzer,
            structurefactor = StructureFactorAnalyzer,
            density         = DensityAnalyzer
            )

        self.quantities = self.scalars | self.fields

        self.ignorable_estimators=set(['LocalEnergy'])

        self.quantity_aliases=dict()
        for q in self.analyzer_quantities:
            self.quantity_aliases[q]=q
        #end for

        self.future_quantities=set(['StructureFactor','MomentumDistribution'])
        return
    #end def __init__
#end class QmcpackCapabilities


QAanalyzer.capabilities = QmcpackAnalyzerCapabilities()



class QmcpackAnalysisRequest(QAobject):
    def __init__(self,source=None,destination=None,savefile='',
                 methods=None,calculations=None,data_sources=None,quantities=None,
                 warmup_calculations=None,
                 output=set(['averages','samples']),
                 ndmc_blocks=1000,equilibration=None,group_num=None,
                 traces=False,dm_settings=None):
        self.source          = source          
        self.destination     = destination     
        self.savefile        = str(savefile)
        self.output          = set(output)
        self.ndmc_blocks     = int(ndmc_blocks)
        self.group_num       = group_num
        self.traces          = traces
        self.dm_settings     = dm_settings

        cap = QAanalyzer.capabilities

        if methods is None:
            self.methods = set(cap.methods)
        else:
            self.methods = set(methods) & cap.methods
        #end if
        if calculations is None:
            self.calculations = set()
        else:
            self.calculations = set(calculations)
        #end if
        if data_sources is None:
            self.data_sources = set(cap.data_sources)
        else:
            self.data_sources = set(data_sources) & cap.data_sources
        #end if
        if quantities is None:
            self.quantities   = set(cap.quantities)
        else:
            quants = set()
            for q in quantities:
                qc = self.condense_name(q)
                quants.add(qc)
            #end for
            self.quantities  = quants & cap.quantities
        #end if
        if warmup_calculations is None:
            self.warmup_calculations = set()
        else:
            self.warmup_calculations = set(warmup_calculations)
        #end if
        if isinstance(equilibration,(dict,obj)):
            eq = obj()
            eq.transfer_from(equilibration)
        else:
            eq = equilibration
        #end if
        self.equilibration = eq

        return
    #end def __init__

    def complete(self):
        spath,sfile = os.path.split(self.source)
        if spath=='':
            self.source = os.path.join('./',self.source)
        #end if
        if self.destination==None:
            self.destination = os.path.split(self.source)[0]
        #end if
        return True
    #end def complete
#end class QmcpackAnalysisRequest:



"""
class QmcpackAnalyzer
  used to analyze all data produced by QMCPACK

  usage:
     results = QmcpackAnalyzer("qmcpack.in.xml")
       |  QMC methods used and observables estimated are determined
       \  Each observable is calculated by an object contained in results
"""
class QmcpackAnalyzer(SimulationAnalyzer,QAanalyzer):
    def __init__(self,arg0=None,**kwargs):

        verbose = False
        if 'verbose' in kwargs:
            verbose=kwargs['verbose']
            del kwargs['verbose']
        #end if
        QAanalyzer.verbose_vlog = verbose or QAanalyzer.verbose_vlog

        nindent = 0
        if 'nindent' in kwargs:
            nindent=kwargs['nindent']
            del kwargs['nindent']
        #end if
        QAanalyzer.__init__(self,nindent=nindent)

        analyze = False
        if 'analyze' in kwargs:
            analyze=kwargs['analyze']
            del kwargs['analyze']
        #end if

        if isinstance(arg0,Simulation):
            sim = arg0
            if 'analysis_request' in sim:
                request = sim.analysis_request.copy()
            else:
                request = QmcpackAnalysisRequest(
                    source = os.path.join(sim.resdir,sim.infile),
                    destination = sim.resdir
                    )
                if 'stat' in request.data_sources:
                    request.data_sources.remove('stat')
                #end if
                if 'storeconfig' in request.data_sources:
                    request.data_sources.remove('storeconfig')
                #end if
                if 'traces' in request.data_sources:
                    request.data_sources.remove('traces')
                #end if
            #end if
        elif isinstance(arg0,QmcpackAnalysisRequest):
            request = arg0
        elif isinstance(arg0,str):
            kwargs['source']=arg0
            request = QmcpackAnalysisRequest(**kwargs)
        else:
            if 'source' not in kwargs:
                kwargs['source']='./qmcpack.in.xml'
            #end if
            request = QmcpackAnalysisRequest(**kwargs)
        #end if

        self.change_request(request)

        if request!=None and os.path.exists(request.source):
            self.init_sub_analyzers(request)
        #end if

        savefile = request.savefile
        savefilepath = os.path.join(request.destination,request.savefile)
        self.info.savefile = savefile
        self.info.savefilepath = savefilepath
        self.info.error = None
        if os.path.exists(savefilepath) and savefile!='':
            self.load()
        elif analyze:
            self.analyze()
        #end if

        return
    #end def __init__


    def change_request(self,request):
        if not isinstance(request,QmcpackAnalysisRequest):
            self.error('input request must be a QmcpackAnalysisRequest',exit=False)
            self.error('  type provided: '+str(type(request)))
        #end if
        request.complete()
        self.info.request = request
    #end def change_request



    def init_sub_analyzers(self,request=None):        
        own_request = request==None
        if request==None:
            request = self.info.request
        #end if
        group_num = request.group_num
        
        #determine if the run was bundled
        if request.source.endswith('.xml'):
            self.info.type = 'single'
        else:
            self.info.type = 'bundled'
            self.bundle(request.source)
            return
        #end if

        self.vlog('reading input file: '+request.source,n=1)
        input = QmcpackInput(request.source)
        input.pluralize()
        input.unroll_calculations()
        calculations = input.simulation.calculations
        self.info.set(
            input = input,
            ordered_input = input.read_xml(request.source)
            )

        project,wavefunction = input.get('project','wavefunction')
        wavefunction = wavefunction.get_single('psi0')

        subindent = self.subindent()

        self.wavefunction = WavefunctionAnalyzer(wavefunction,nindent=subindent)

        self.vlog('project id: '+project.id,n=1)
        file_prefix  = project.id
        if group_num!=None:
            group_ext = '.g'+str(group_num).zfill(3)
            if not file_prefix.endswith(group_ext):
                file_prefix += group_ext
            #end if
        elif self.info.type=='single':
            resdir,infile = os.path.split(request.source)
            #ifprefix = infile.replace('.xml','')
            ifprefix = infile.replace('.xml','.')
            ls = os.listdir(resdir)
            for filename in ls:
                if filename.startswith(ifprefix) and filename.endswith('.qmc'):
                    group_tag = filename.split('.')[-2]
                    #file_prefix = 'qmc.'+group_tag
                    file_prefix = project.id+'.'+group_tag
                    break
                #end if
            #end for
        #end if
        if 'series' in project:
            series_start = int(project.series)
        else:
            series_start = 0
        #end if

        self.vlog('data file prefix: '+file_prefix,n=1)

        run_info = obj(
            file_prefix  = file_prefix,
            series_start = series_start,
            source_path  = os.path.split(request.source)[0],
            group_num    = group_num,
            system       = input.return_system()
            )
        self.info.transfer_from(run_info)

        self.set_global_info()        

        if len(request.calculations)==0:
            request.calculations = set(series_start+arange(len(calculations)))
        #end if

        method_aliases = dict()
        for method in self.opt_methods:
            method_aliases[method]='opt'
        #end for
        for method in self.vmc_methods:
            method_aliases[method]='vmc'
        #end for
        for method in self.dmc_methods:
            method_aliases[method]='dmc'
        #end for
        
        method_objs = ['qmc','opt','vmc','dmc']
        for method in method_objs:
            self[method] = QAanalyzerCollection()
        #end for
        for index,calc in calculations.iteritems():
            method = calc.method
            if method in method_aliases:
                method_type = method_aliases[method]
            else:
                self.error('method '+method+' is unrecognized')
            #end if
            if method_type in request.methods:
                series = series_start + index
                if series in request.calculations:
                    if method in self.opt_methods:
                        qma = OptAnalyzer(series,calc,input,nindent=subindent)
                        primary = self.opt
                    elif method in self.vmc_methods:
                        qma = VmcAnalyzer(series,calc,input,nindent=subindent)
                        primary = self.vmc
                    elif method in self.dmc_methods:
                        qma = DmcAnalyzer(series,calc,input,nindent=subindent)
                        primary = self.dmc
                    #end if
                    primary[series]  = qma
                    self.qmc[series] = qma
                #end if
            #end if
        #end for            
        for method in method_objs:
            if len(self[method])==0:
                del self[method]
            #end if
        #end for

        #Check for multi-qmc results such as
        # optimization or timestep studies
        results = QAanalyzerCollection()
        if 'opt' in self and len(self.opt)>0:
            optres = OptimizationAnalyzer(input,self.opt,nindent=subindent)
            results.optimization = optres
        #end if
        if 'dmc' in self and len(self.dmc)>1:
            maxtime = 0
            times = dict()
            for series,dmc in self.dmc.iteritems():
                blocks,steps,timestep = dmc.info.method_input.list('blocks','steps','timestep')
                times[series] = blocks*steps*timestep
                maxtime = max(times[series],maxtime)
            #end for
            dmc = QAanalyzerCollection()            
            for series,time in times.iteritems():
                if abs(time-maxtime)/maxtime<.5:
                    dmc[series] = self.dmc[series]
                #end if
            #end for
            if len(dmc)>1:
                results.timestep_study = TimestepStudyAnalyzer(dmc,nindent=subindent)
            #end if
        #end if

        if len(results)>0:
            self.results = results
        #end if

        self.unset_global_info()

    #end def init_sub_analyzers


    def set_global_info(self):
        QAanalyzer.request  = self.info.request
        QAanalyzer.run_info = self.info
    #end def set_global_info

    def unset_global_info(self):
        QAanalyzer.request  = None
        QAanalyzer.run_info = None
    #end def unset_global_info


    def load_data(self):
        request = self.info.request
        if not os.path.exists(request.source):
            self.error('path to source\n  '+request.source+'\n  does not exist\n ensure that request.source points to a valid qmcpack input file')
        #end if
        self.set_global_info()
        self.propagate_indicators(data_loaded=False)
        if self.info.type=='bundled' and self.info.perform_bundle_average:
            self.prevent_average_load()
        #end if        
        QAanalyzer.load_data(self)
        if self.info.type=='bundled' and self.info.perform_bundle_average:
            self.average_bundle_data()
        #end if
        self.unset_global_info()
    #end def load_data


    def analyze(self,force=False):
        if not self.info.analyzed or force:
            if not self.info.data_loaded:
                self.load_data()
            #end if
            self.vlog('main analysis of QmcpackAnalyzer data',n=1)
            try:
                self.set_global_info()
                self.propagate_indicators(analyzed=False)
                if self.info.type!='bundled':
                    QAanalyzer.analyze(self,force=force)
                else:
                    for analyzer in self.bundled_analyzers:
                        analyzer.analyze()
                    #end for
                    QAanalyzer.analyze(self,force=force)
                #end if
                self.unset_global_info()
            except:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
                msg = ''
                for line in lines:
                    msg+=line
                #end for
                self.info.error = exc_type
                self.warn('runtime exception encountered\n'+msg)
            #end try
            self.vlog('end main analysis of QmcpackAnalyzer data',n=1)
            if self.info.request.savefile!='':
                self.save()
            #end if
        #end if
    #end def analyze



    def bundle(self,source):
        self.vlog('bundled run detected',n=1)
        if os.path.exists(source):
            fobj = open(source,'r')
            lines = fobj.read().split('\n')
            fobj.close()
        else:
            self.error('source file '+source+' does not exist')
        #end if
        infiles = []
        for line in lines:
            ls = line.strip()
            if ls!='':
                infiles.append(ls)
            #end if
        #end for
        self.info.input_infiles = list(infiles)
        analyzers = QAanalyzerCollection()
        request = self.info.request
        path = os.path.split(request.source)[0]
        files = os.listdir(path)
        outfiles = []
        for file in files:
            if file.endswith('qmc'):
                outfiles.append(file)
            #end if
        #end for
        del files
        for i in range(len(infiles)):
            infile = infiles[i]
            prefix = infile.replace('.xml','')
            gn = i
            for outfile in outfiles:
                if outfile.startswith(prefix):
                    gn = int(outfile.split('.')[-2][1:])
                    break
                #end if
            #end for
            req = request.copy()
            req.source = os.path.join(path,infile)
            req.group_num = gn
            qa = QmcpackAnalyzer(req,nindent=self.subindent())
            #qa.init_sub_analyzers(group_num=gn)
            analyzers[gn] = qa
        #end for
        self.bundled_analyzers = analyzers
        self.info.perform_bundle_average = False
        #check to see if twist averaging
        #  indicated by distinct twistnums
        #  or twist in all ids
        twistnums = set()
        twist_ids = True
        for analyzer in analyzers:
            input = analyzer.info.input
            twistnum = input.get('twistnum')
            project = input.get('project')
            if twistnum!=None:
                twistnums.add(twistnum)
            #end if
            twist_ids = twist_ids and 'twist' in project.id
        #end for
        distinct_twistnums = len(twistnums)==len(analyzers)
        twist_averaging = distinct_twistnums or twist_ids
        if twist_averaging:
            self.info.perform_bundle_average = True
        #end if
        example = analyzers.list()[0]
        input,system = example.info.tuple('input','system')
        self.info.set(
            input  = input.copy(),
            system = system.copy()
            )
        self.vlog('average over bundled runs?  {0}'.format(self.info.perform_bundle_average),n=1)
    #end def bundle


    def prevent_average_load(self):
        for method_type in self.capabilities.methods:
            if method_type in self:
                self[method_type].propagate_indicators(data_loaded=True)
            #end if
        #end for
    #end def prevent_average_load


    def average_bundle_data(self):
        analyzers = self.bundled_analyzers
        if len(analyzers)>0:
            self.vlog('performing bundle (e.g. twist) averaging',n=1)
            #create local data structures to match those in the bundle
            example = analyzers.list()[0].copy()
            for method_type in self.capabilities.methods:
                if method_type in self:
                    del self[method_type]
                #end if
                if method_type in example:
                    self.vlog('copying {0} methods from analyzer 0'.format(method_type),n=2)
                    self[method_type] = example[method_type]
                #end if            
            #end if
            if 'qmc' in self:
                del self.qmc
            #end if            
            if 'qmc' in example:
                self.vlog('copying qmc methods from analyzer 0',n=2)
                self.qmc = example.qmc
            #end if
            if 'wavefunction' in self:
                del self.wavefunction
            #end if            
            if 'wavefunction' in example:
                self.vlog('copying wavefunction from analyzer 0',n=2)
                self.wavefunction = example.wavefunction
            #end if
            del example

            if 'qmc' in self:
                #zero out the average data
                self.vlog('zeroing own qmc data',n=2)
                for qmc in self.qmc:
                    qmc.zero_data()
                #end for

                #resize the average data
                self.vlog('finding minimum data size (for incomplete runs)',n=2)
                for analyzer in analyzers:
                    for series,qmc in self.qmc.iteritems():
                        qmc.minsize_data(analyzer.qmc[series])
                    #end for
                #end for
    
                #accumulate the average data
                self.vlog('accumulating data from bundled runs',n=2)
                for analyzer in analyzers:
                    for series,qmc in self.qmc.iteritems():
                        qmc.accumulate_data(analyzer.qmc[series])
                    #end for
                #end for
    
                #normalize the average data
                norm_factor = len(analyzers)
                self.vlog('normalizing bundle average (factor={0})'.format(norm_factor),n=2)
                for qmc in self.qmc:
                    qmc.normalize_data(norm_factor)
                #end for
            #end if
        #end if
    #end def average_bundle_data




    def save(self,filepath=None,overwrite=True):
        if filepath==None:
            filepath = self.info.savefilepath
        #end if
        self.vlog('saving QmcpackAnalyzer in file {0}'.format(filepath),n=1)
        if not overwrite and os.path.exists(filepath):
            return
        #end if
        self._unlink_dynamic_methods()
        self.saved_global = QAobject._global
        self._save(filepath)
        self._relink_dynamic_methods()
        return
    #end def save

    def load(self,filepath=None):
        if filepath==None:
            filepath = self.info.savefilepath
        #end if
        self.vlog('loading QmcpackAnalyzer from file {0}'.format(filepath),n=1)
        self._load(filepath)
        QAobject._global = self.saved_global
        del self.saved_global
        self._relink_dynamic_methods()
        return
    #end def load




    def check_traces(self,verbose=False,pad=None,header=None):
        if pad is None:
            pad = ''
        #end if
        if header is None:
            header = '\nChecking traces'
        #end if
        if 'qmc' in self:
            if verbose:
                self.log(pad+header)
                pad += '  '
            #end if
            for method in self.qmc:
                method.check_traces(pad)
            #end for
        else:
            if verbose:
                self.log(pad+'\nNo traces to check')
            #end if
            return None
        #end if
    #end def check_traces


    def plot_trace(self,quantity,style='b-',offset=0,source='scalar',mlabels=True,
                   mlines=True,show=True,alloff=False):
        mlabels &= not alloff
        mlines  &= not alloff
        show    &= not alloff
        shw = show
        offset = int(offset)
        id = self.info.input.get('project').id
        sdata = obj()
        series = sorted(self.qmc.keys())
        q = []
        soffset = offset
        for s in series:
            qmc = self.qmc[s]
            method = qmc.info.method
            if source=='scalar' or method=='vmc':
                src = qmc.scalars.data
            elif source=='dmc':
                src = qmc.dmc.data
            else:
                self.error('invalid source: '+source)
            #end if
            if quantity in src:
                qn = list(src[quantity])
            else:
                qn = len(src.LocalEnergy)*[0]
            #end if
            q.extend(qn)
            sdata[s] = obj(
                mlab = method+' '+str(s),
                mloc = soffset + len(qn)/2,
                line_loc = soffset + len(qn)-1
                )
            soffset += len(qn)
        #end for
        q = array(q)
        qmin = q.min()
        qmax = q.max()
        mlabel_height = qmin + .8*(qmax-qmin)
        if shw:
            figure()
        #end if
        plot(offset+arange(len(q)),q,style,label=id)
        for s in series:
            sd = sdata[s]
            if mlabels:
                text(sd.mloc,mlabel_height,sd.mlab)
            #end if
            if mlines:
                plot([sd.line_loc,sd.line_loc],[qmin,qmax],'k-')
            #end if
        #end for
        if shw:
            title('{0} vs series for {1}'.format(quantity,id))
            xlabel('blocks')
            ylabel(quantity)
            legend()
            show()
        #end if
    #end def plot_trace
          
#end class QmcpackAnalyzer


