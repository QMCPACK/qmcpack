##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_result_analyzers.py                                       #
#    Analyzer classes for results of multi-step processes carried    #
#    out by QMCPACK.  This includes basic analysis of wavefunction   #
#    optimization and DMC timestep studies.                          #
#                                                                    #
#  Content summary:                                                  #
#    ResultAnalyzer                                                  #
#      Empty base class to distinguish result analyzers from other   #
#      types.                                                        #
#                                                                    #
#    OptimizationAnalyzer                                            #
#      Supports analysis of optimization convergence, including      #
#      plots of energy and variance vs. series and plots detailing   #
#      the convergence of bspline Jastrows.                          #
#                                                                    #
#    TimestepStudyAnalyzer                                           #
#      Supports basic plotting and reporting of DMC timestep study   #
#      data.                                                         #
#                                                                    #
#====================================================================#


from numpy import array,empty,zeros,sqrt,arange
from generic import obj
from unit_converter import convert
from qmcpack_input import QmcpackInput
from qmcpack_analyzer_base import QAobject,QAanalyzer
from debug import *



class ResultAnalyzer(QAanalyzer):
    None
#end class ResultAnalyzer


class OptimizationAnalyzer(ResultAnalyzer):
    def __init__(self,input,opts,energy_weight=None,variance_weight=None,nindent=0):
        QAanalyzer.__init__(self,nindent=nindent)

        self.opts  = opts
        self.energy_weight = energy_weight
        self.variance_weight = variance_weight


        ew,vw = energy_weight,variance_weight
        if ew==None or vw==None:
            opts_in = []            
            for qmc in input.simulation.calculations:
                if qmc.method in self.opt_methods:
                    opts_in.append(qmc)
                #end if
            #end for
            optin = opts_in[-1] #take cost info from the last optimization section
            curv,crv,ce = optin.get(['unreweightedvariance','reweightedvariance','energy'])
            tol=1e-4
            if crv>tol:
                cv = crv
            elif curv>tol:
                cv = curv
            else:
                cv = 0
            #end if
            if ew==None:
                ew = ce
            #end if
            if vw==None:
                vw = cv
            #end if
        #end if

        if self.optimize=='lastcost':
            self.optimize = ew,vw
        #end if

    #end def __init__


    def init_sub_analyzers(self):
        None
    #end def init_sub_analyzers

    def analyze_local(self):
        input = QAanalyzer.run_info.input
        self.info.system = QAanalyzer.run_info.system
        opts  = self.opts
        ew    = self.energy_weight
        vw    = self.variance_weight

        Efail        = 1e6
        Vfail        = 1e3
        EVratio_fail = 0.30
        EVratio_soft_fail = 0.15
        
        #save the energies and variances of opt iterations
        res = obj()
        variance_present = False
        any_complete = False
        all_complete = True
        unstable = False
        any_stable = False
        for s,opt in opts.iteritems():
            complete = opt.info.complete
            any_complete |= complete
            all_complete &= complete
            if complete:
                fail = False
                le = opt.scalars.LocalEnergy
                en     = le.mean
                enerr  = le.error
                fail |= abs(en)>Efail
                if 'LocalEnergyVariance' in opt.scalars:
                    variance_present = True
                    lev = opt.scalars.LocalEnergyVariance
                    va    = lev.mean
                    vaerr = lev.error
                    fail |= abs(va)>Vfail or abs(va/en)>EVratio_fail
                #end if
                if not fail:
                    any_stable = True
                    sres = obj()
                    sres.en    = en
                    sres.enerr = enerr
                    if variance_present:
                        sres.va    = va
                        sres.vaerr = vaerr
                    #end if
                    res[s] = sres
                #end if
                unstable|=fail
            #end if
        #end for
        unstable |= not any_complete

        nseries = len(res)
        en    = zeros((nseries,),dtype=float)
        enerr = zeros((nseries,),dtype=float)
        va    = zeros((nseries,),dtype=float)
        vaerr = zeros((nseries,),dtype=float)

        series = array(sorted(res.keys()),dtype=int)
        i = 0
        for s in series:
            sres = res[s]
            en[i]    = sres.en
            enerr[i] = sres.enerr
            if variance_present:
                va[i]    = sres.va
                vaerr[i] = sres.vaerr
            #end if
            i+=1
        #end for
                

        self.set(
            any_complete   = any_complete,
            all_complete   = all_complete,
            unstable       = unstable,
            series         = series,
            energy         = en,
            energy_error   = enerr,
            variance       = va,
            variance_error = vaerr,
            )


        #find the optimal coefficients
        optimize = self.optimize
        if variance_present and optimize=='variance':
            ew = 0.0
            vw = 1.0
        elif optimize=='energy':
            ew = 1.0
            vw = 0.0
        elif optimize=='energy_within_variance_tol' or optimize=='ewvt':
            None
        elif optimize=='last':
            None
        elif isinstance(optimize,(tuple,list)) and len(optimize)==2:
            ew,vw = optimize
        else:
            self.error('selection for optimization is invalid\noptimize setting: {0}\nvalid options are: energy, variance, energy_within_variance_tol, or a length 2 tuple containing the cost of energy and variance, e.g. (.5,.5)'.format(optimize))
        #end if

        self.failed = True
        self.optimal_series = None
        self.optimal_file   = None
        self.optimal_wavefunction = None
        if any_stable:
            if optimize=='energy_within_variance_tol' or optimize=='ewvt':
                indices = arange(len(series),dtype=int)
                vartol  = 0.2
                vmin    = va.min()
                vind    = indices[abs(va-vmin)/vmin<vartol]
                index   = vind[en[vind].argmin()]
                opt_series = series[index]
            elif optimize=='last':
                index = len(en)-1
                opt_series = series[index]
            else:
                cost = en*ew+va*vw
                index = cost.argmin()
                opt_series = series[index]
            #end if
            failed = abs(en[index])>Efail or abs(va[index])>Vfail or abs(va[index]/en[index])>EVratio_soft_fail 

            self.failed = failed
            self.optimal_series = opt_series
            self.optimal_file = opts[opt_series].info.files.opt
            self.optimal_wavefunction = opts[opt_series].wavefunction.info.wfn_xml.copy()
        #end if
    #end def analyze_local


    def summarize(self,units='eV',norm=1.,energy=True,variance=True,header=True):
        if isinstance(norm,str):
            norm = norm.replace('_',' ').replace('-',' ')
            if norm=='per atom':
                norm = len(self.info.system.structure.elem)
            else:
                self.error('norm must be a number or "per atom"\n you provided '+norm)
            #end if
        #end if
        econv = convert(1.0,'Ha',units)/norm
        en    = econv*self.energy
        enerr = econv*self.energy_error
        va    = econv**2*self.variance
        vaerr = econv**2*self.variance_error
        emax = en.max()
        vmax = va.max()
        if header:
            print 'Optimization summary:'
            print '===================='
        #end if
        if energy:
            if header:
                print '  Energies ({0}):'.format(units)
            #end if
            for i in range(len(en)):
                print '    {0:>2}    {1:9.6f} +/-{2:9.6f}'.format(i,en[i]-emax,enerr[i])
            #end for
            print '    ref {0:9.6f}'.format(emax)
        #end if
        if variance:
            if header:
                print '  Variances ({0}^2):'.format(units)
            #end if
            for i in range(len(en)):
                print '    {0:>2}    {1:9.6f} +/- {2:9.6f}'.format(i,va[i],vaerr[i])
            #end for
        #end if
    #end def summarize


    def plot_opt_convergence(self,title=None,saveonly=False):
        if title is None:
            ts = 'Optimization: Energy/Variance Convergence'
        else:
            ts = title
        #end if
        from matplotlib.pyplot import figure,subplot,xlabel,ylabel,plot,errorbar,title,xticks,xlim

        opt = self.opts
        nopt = len(opt)
        if nopt==0:
            return
        #end if

        en    = self.energy
        enerr = self.energy_error
        va    = self.variance
        vaerr = self.variance_error

        #plot energy and variance
        figure()
        r = range(nopt)
        subplot(3,1,1)
        errorbar(r,en,enerr,fmt='b')
        ylabel('Energy (Ha)')
        title(ts)
        xticks([])
        xlim([r[0]-.5,r[-1]+.5])
        subplot(3,1,2)
        errorbar(r,va,vaerr,fmt='r')
        ylabel('Var. ($Ha^2$)')
        xticks([])
        xlim([r[0]-.5,r[-1]+.5])
        subplot(3,1,3)
        plot(r,abs(sqrt(va)/en),'k')
        ylabel('Var.^(1/2)/|En.|')
        xlabel('Optimization attempts')
        xticks(r)
        xlim([r[0]-.5,r[-1]+.5])
    #end def plot_opt_convergence
    

    def plot_jastrow_convergence(self,title=None,saveonly=False,optconv=True):
        if title is None:
            tsin = None
        else:
            tsin = title
        #end if
        from matplotlib.pyplot import figure,subplot,xlabel,ylabel,plot,errorbar,title,xticks,xlim

        opt = self.opts
        nopt = len(opt)
        if nopt==0:
            return
        #end if

        if optconv:
            self.plot_opt_convergence(saveonly=saveonly)
        #end if

        #plot Jastrow functions
        w = opt[0].wavefunction
        jtypes = w.jastrow_types
        order = QAobject()
        for jt in jtypes:
            if jt in w:
                order[jt] = list(w[jt].__dict__.keys())
                order[jt].sort()
            #end if
        #end for
        cs = array([1.,0,0])
        ce = array([0,0,1.])
        for jt in jtypes:
            if jt in w:
                figure()
                nsubplots = len(order[jt])
                n=0
                for o in order[jt]:
                    n+=1
                    subplot(nsubplots,1,n)
                    if n==1:
                        if tsin is None:
                            ts = 'Optimization: '+jt+' Convergence'
                        else:
                            ts = tsin
                        #end if
                        title(ts)
                    #end if
                    for i in range(len(opt)):
                        f = float(i)/len(opt)
                        c = f*ce + (1-f)*cs
                        J = opt[i].wavefunction[jt][o]
                        J.plot(color=c)
                    #end for
                    ylabel(o)
                #end for
                xlabel('r (Bohr)')
            #end if
        #end for
    #end def plot_jastrow_convergence


#end class OptimizationAnalyzer






class TimestepStudyAnalyzer(ResultAnalyzer):
    def __init__(self,dmc,nindent=0):
        QAanalyzer.__init__(self,nindent=nindent)
        self.set(
            dmc = dmc,
            timesteps = [],
            energies  = [],
            errors    = []
            )
    #end def __init__

    def init_sub_analyzers(self):
        None
    #end def init_sub_analyzers

    def analyze_local(self):
        timesteps = []
        energies  = []
        errors    = []
        for dmc in self.dmc:
            timesteps.append(dmc.info.method_input.timestep)
            energies.append(dmc.scalars.LocalEnergy.mean)
            errors.append(dmc.scalars.LocalEnergy.error)
        #end for
        timesteps = array(timesteps)
        energies  = array(energies)
        errors    = array(errors)
        order = timesteps.argsort()
        self.timesteps = timesteps[order]
        self.energies  = energies[order]
        self.errors    = errors[order]
    #end def analyze_local

    def summarize(self,units='eV',header=True):
        timesteps = self.timesteps
        energies  = convert(self.energies.copy(),'Ha',units)
        errors    = convert(self.errors.copy(),'Ha',units)
        Esmall = energies[0]
        if header:
            print 'Timestep study summary:'
            print '======================'
        #end if
        for i in range(len(timesteps)):
            ts,E,Eerr = timesteps[i],energies[i],errors[i]
            print '    {0:>6.4f}   {1:>6.4f} +/- {2:>6.4f}'.format(ts,E-Esmall,Eerr)
        #end for
    #end def summarize

    def plot_timestep_convergence(self):
        from matplotlib.pyplot import figure,subplot,xlabel,ylabel,plot,errorbar,title,text,xticks,rcParams,savefig,xlim

        params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
        rcParams.update(params) 


        timesteps = self.timesteps
        energies  = convert(self.energies.copy(),'Ha','eV')
        errors    = convert(self.errors.copy(),'Ha','eV')
        Esmall = energies[0]

        figure()
        tsrange = [0,1.1*timesteps[-1]]
        plot(tsrange,[0,0],'k-')
        errorbar(timesteps,energies-Esmall,errors,fmt='k.')
        text(array(tsrange).mean(),0,'{0:6.4f} eV'.format(Esmall))
        xticks(timesteps)
        xlim(tsrange)
        xlabel('Timestep (Ha)')
        ylabel('Total Energy (eV)')
        title('DMC Timestep Convergence')

        savefig('TimestepConvergence.png',format='png',bbox_inches ='tight',pad_inches=1)
    #end def plot_timestep_convergence
#end class TimestepStudyAnalyzer
