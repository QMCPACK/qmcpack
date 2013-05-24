import re
import copy
from numpy import array,floor,sqrt,zeros,prod,dot,ones,empty,min,max
from numpy import pi,sin,cos,arccos as acos,arctan2 as atan2
from numpy.linalg import inv,det
from extended_numpy import simplestats,ndgrid,ogrid,arange,simstats
from hdfreader import HDFgroup
from qaobject import QAobject

#simple constants
o2pi = 1./(2.*pi)

#simple functions
def is_integer(i):
    return abs(floor(i)-i)<1e-6
#end def is_integer


class SpaceGridInitializer(QAobject):
    def __init__(self):
        self.coord              = None # string
        return
    #end def __init__

    def check_complete(self,exit_on_fail=True):
        succeeded = True
        for k,v in self._iteritems():
            if v==None:
                succeeded=False
                if exit_on_fail:
                    self.error('  SpaceGridInitializer.'+k+' must be provided',exit=False)
                #end if
            #end if
        #end if
        if not succeeded and exit_on_fail:
            self.error('  SpaceGridInitializer is incomplete')
        #end if
        return succeeded
    #end def check_complete
#end class SpaceGridInitializer


class SpaceGridBase(QAobject):
    cnames=['cartesian','cylindrical','spherical','voronoi']
    coord_s2n = dict()
    coord_n2s = dict()
    i=0
    for name in cnames:
        exec name+'='+str(i)
        coord_s2n[name]=i
        coord_n2s[i]=name
        i+=1
    #end for

    xlabel = 0
    ylabel = 1
    zlabel = 2
    rlabel = 3
    plabel = 4
    tlabel = 5
    axlabel_s2n = {'x':xlabel,'y':ylabel,'z':zlabel,'r':rlabel,'phi':plabel,'theta':tlabel}
    axlabel_n2s = {xlabel:'x',ylabel:'y',zlabel:'z',rlabel:'r',plabel:'phi',tlabel:'theta'}

    axindex = {'x':0,'y':1,'z':2,'r':0,'phi':1,'theta':2}

    quantities=['D','T','V','E','P']

    def __init__(self,initobj,options):
        if options==None:
            options = QAobject()
            options.wasNone = True
            options.points       = None
            options.exit_on_fail = True
            options.nblocks_exclude = 0
        else:
            if 'points' not in options: 
                options.points = None
            if 'exit_on_fail' not in options:
                options.exit_on_fail = True
            if 'nblocks_exclude' not in options:
                options.nblocks_exclude = 0
        #end if

        self.points          = options.points
        self.init_exit_fail  = options.exit_on_fail
        self.nblocks_exclude = options.nblocks_exclude
        self.keep_data = True
        delvars = ['init_exit_fail','keep_data']
            
        self.coord          = None # string
        self.coordinate     = None
        self.ndomains       = None
        self.domain_volumes = None
        self.domain_centers = None
        self.nvalues_per_domain = -1
        self.nblocks            = -1
        self.D  = QAobject() #Number Density
        self.T  = QAobject() #Kinetic Energy Density
        self.V  = QAobject() #Potential Energy Density
        self.E  = QAobject() #Energy Density, T+V
        self.P  = QAobject() #Local Pressure, (Volume)*P=(2*T+V)/3


        self.init_special()

        if initobj==None:
            return
        #end if

        self.DIM=3

        iname = initobj.__class__.__name__
        self.iname=iname
        if iname==self.__class__.__name__+'Initializer':
            self.init_from_initializer(initobj)
        elif iname==self.__class__.__name__:
            self.init_from_spacegrid(initobj)
        elif iname=='HDFgroup':
            self.init_from_hdfgroup(initobj)
        elif iname=='XMLelement':
            self.init_from_xmlelement(initobj)
        else:
            self.error('Spacegrid cannot be initialized from '+iname)
        #end if
        delvars.append('iname')

        self.check_complete()

        for dv in delvars:
            del self[dv]
        #end for

        self._reset_dynamic_methods()
        self._register_dynamic_methods()
        return
    #end def __init__

    def copy(self,other):
        None
    #end def copy

    def init_special(self):
        None
    #end def init_special

    def init_from_initializer(self,init):
        None
    #end def init_from_initializer

    def init_from_spacegrid(self,init):
        None
    #end def init_from_spacegrid

    def init_from_hdfgroup(self,init):
        #copy all datasets from hdf group
        value_pattern = re.compile('value')
        gmap_pattern = re.compile(r'gmap\d*')
        for k,v in init._iteritems():
            exclude = k[0]=='_' or gmap_pattern.match(k) or value_pattern.match(k)
            if not exclude:
                self.__dict__[k]=v                
            #end if
        #end for

        #convert 1x and 1x1 numpy arrays to just numbers
        #convert Nx1 and 1xN numpy arrays to Nx arrays
        array_type = type(array([]))
        exclude = set(['value','value_squared'])
        for k,v in self._iteritems():
            if k[0]!='_' and type(v)==array_type and k not in exclude:
                sh=v.shape
                ndim = len(sh)
                if ndim==1 and sh[0]==1:
                    self.__dict__[k]=v[0]
                elif ndim==2:
                    if sh[0]==1 and sh[1]==1:
                        self.__dict__[k]=v[0,0]
                    elif sh[0]==1 or sh[1]==1:
                        self.__dict__[k]=v.reshape((sh[0]*sh[1],))
                    #end if
                #end if
            #end if
        #end for

        #set coord string
        self.coord = SpaceGridBase.coord_n2s[self.coordinate]

        #determine if chempot grid
        chempot = 'min_part' in init
        self.chempot = chempot
        if chempot:
            npvalues = self.max_part-self.min_part+1
            self.npvalues = npvalues
        #end if

        #process the data in hdf value,value_squared
        nbe = self.nblocks_exclude
        nquant   = self.nvalues_per_domain
        ndomains = self.ndomains
        nblocks,ntmp = init.value.shape
        self.nblocks = nblocks

        if not chempot:
            value = init.value.reshape(nblocks,ndomains,nquant).transpose(2,1,0)
        else:
            value = init.value.reshape(nblocks,ndomains,npvalues,nquant).transpose(3,2,1,0)
        #end if
        value = value[...,nbe:]

        #(mean,error)=simplestats(value)
        (mean,var,error,kappa)=simstats(value)
        quants = ['D','T','V']
        for i in range(len(quants)):
            q=quants[i]
            self[q].mean  =  mean[i,...]
            self[q].error = error[i,...]
            exec 'i'+q+'='+str(i)
        #end for
        
        E = value[iT,...]+value[iV,...]
#        (mean,error)=simplestats(E)
        (mean,var,error,kappa)=simstats(E)
        self.E.mean  =  mean
        self.E.error = error
        
        P = 2./3.*value[iT,...]+1./3.*value[iV,...]
        #(mean,error)=simplestats(P)
        (mean,var,error,kappa)=simstats(P)
        self.P.mean  =  mean
        self.P.error = error


        #convert all quantities into true densities
        ovol = 1./self.domain_volumes
        sqovol = sqrt(ovol)
        for q in SpaceGridBase.quantities:
            self[q].mean  *= ovol
            self[q].error *= sqovol
        #end for

        #keep original data, if requested
        if self.keep_data:
            self.data = QAobject()
            for i in range(len(quants)):
                q=quants[i]
                self.data[q] = value[i,...]
            #end for
            self.data.E = E
            self.data.P = P
        #end if

        #print 'sg'
        #import code
        #code.interact(local=locals())
            
        return
    #end def init_from_hdfgroup

    def init_from_xmlelement(self,init):
        None
    #end def init_from_xmlelement

    def check_complete(self,exit_on_fail=True):
        succeeded = True
        for k,v in self._iteritems():
            if k[0]!='_' and v==None:
                succeeded=False
                if exit_on_fail:
                    self.error('SpaceGridBase.'+k+' must be provided',exit=False)
                #end if
            #end if
        #end if
        if not succeeded:
            self.error('SpaceGrid attempted initialization from '+self.iname,exit=False)
            self.error('SpaceGrid is incomplete',exit=False)
            if exit_on_fail:
                exit()
            #end if
        #end if
        return succeeded
    #end def check_complete

    def _reset_dynamic_methods(self):
        None
    #end def _reset_dynamic_methods

    def _unset_dynamic_methods(self):
        None
    #end def _unset_dynamic_methods

    def add_all_attributes(self,o):
        for k,v in o.__dict__.iteritems():
            if not k.startswith('_'):
                vc = copy.deepcopy(v)
                self._add_attribute(k,vc)
            #end if
        #end for
        return
    #end def add_all_attributes


    def reorder_atomic_data(self,imap):
        None
    #end if


    def integrate(self,quantity,domain=None):
        if quantity not in SpaceGridBase.quantities:
            msg = 'requested integration of quantity '+quantity+'\n'
            msg +='  '+quantity+' is not a valid SpaceGrid quantity\n'
            msg +='  valid quantities are:\n'
            msg +='  '+str(SpaceGridBase.quantities)
            self.error(msg)
        #end if
        dv = self.domain_volumes
        if domain==None:
            mean = (self[quantity].mean*dv).sum()
            error = sqrt((self[quantity].error**2*dv).sum())
        else:
            mean = (self[quantity].mean[domain]*dv[domain]).sum()
            error = sqrt((self[quantity].error[domain]**2*dv[domain]).sum())
        #end if
        return mean,error
    #end def integrate

    def integrate_data(self,quantity,*domains,**kwargs):
        return_list = False
        if 'domains' in kwargs:
            domains = kwargs['domains']
            return_list = True
        #end if
        if 'return_list' in kwargs:
            return_list = kwargs['return_list']
        #end if
        if quantity not in SpaceGridBase.quantities:
            msg = 'requested integration of quantity '+quantity+'\n'
            msg +='  '+quantity+' is not a valid SpaceGrid quantity\n'
            msg +='  valid quantities are:\n'
            msg +='  '+str(SpaceGridBase.quantities)
            self.error(msg)
        #end if
        q = self.data[quantity]
        results = list()
        nblocks = q.shape[-1]
        qi = zeros((nblocks,))
        if len(domains)==0:
            for b in xrange(nblocks):
                qi[b] = q[...,b].sum()
            #end for
            (mean,var,error,kappa)=simstats(qi)
        else:
            for domain in domains:
                for b in xrange(nblocks):
                    qb = q[...,b]
                    qi[b] = qb[domain].sum()
                #end for                
                (mean,var,error,kappa)=simstats(qi)
                res = QAobject()
                res.mean  = mean
                res.error = error
                res.data  = qi.copy()
                results.append(res)
            #end for
        #end for
        if len(domains)<2:
            return mean,error
        else:
            if not return_list:
                return tuple(results)
            else:
                means = list()
                errors = list()
                for res in results:
                    means.append(res.mean)
                    errors.append(res.error)
                #end for
                return means,errors
            #end if
        #end if
    #end def integrate_data

#end class SpaceGridBase




class RectilinearGridInitializer(SpaceGridInitializer):
    def __init__(self):
        SpaceGridInitializer.__init__(self)
        self.origin             = None # 3x1 array
        self.axes               = None # 3x3 array
        self.axlabel            = None # 3x1 string list
        self.axgrid             = None # 3x1 string list
    #end def __init__
#end class RectilinearGridInitializer


class RectilinearGrid(SpaceGridBase):
    def __init__(self,initobj=None,options=None):
        SpaceGridBase.__init__(self,initobj,options)
        return
    #end def __init__

    def init_special(self):
        self.origin         = None # 3x1 array
        self.axes           = None # 3x3 array
        self.axlabel        = None # 3x1 string list
        self.axinv          = None
        self.volume         = None
        self.dimensions     = None
        self.gmap           = None
        self.umin           = None
        self.umax           = None
        self.odu            = None
        self.dm             = None
        self.domain_uwidths = None
        return
    #end def init_special

    def copy(self):
        return RectilinearGrid(self)
    #end def copy

    def _reset_dynamic_methods(self):
        p2d=[self.points2domains_cartesian,   \
             self.points2domains_cylindrical, \
             self.points2domains_spherical]
        self.points2domains = p2d[self.coordinate]

        p2u=[self.point2unit_cartesian,   \
             self.point2unit_cylindrical, \
             self.point2unit_spherical]
        self.point2unit = p2u[self.coordinate]
        return
    #end def _reset_dynamic_methods

    def _unset_dynamic_methods(self):
        self.points2domains = None
        self.point2unit     = None
        return
    #end def _unset_dynamic_methods

    def init_from_initializer(self,init):
        init.check_complete()
        for k,v in init._iteritems():
            if k[0]!='_':
                self.__dict__[k]=v
            #end if
        #end for
        self.initialize()
        return
    #end def init_from_initializer

    def init_from_spacegrid(self,init):
        for q in SpaceGridBase.quantities:
            self[q].mean = init[q].mean.copy()
            self[q].error = init[q].error.copy()
        #end for
        array_type = type(array([1]))
        exclude = set(['point2unit','points2domains','points'])
        for k,v in init._iteritems():
            if k[0]!='_':
                vtype = type(v)
                if k in SpaceGridBase.quantities:
                    self[k].mean  = v.mean.copy()
                    self[k].error = v.error.copy()
                elif vtype==array_type:
                    self[k] = v.copy()
                elif vtype==HDFgroup:
                    self[k] = v
                elif k in exclude:
                    None
                else:
                    self[k] = vtype(v)
                #end if
            #end for            
        #end for
        self.points = init.points
        return
    #end def init_from_spacegrid

    def init_from_hdfgroup(self,init):
        SpaceGridBase.init_from_hdfgroup(self,init)
        self.gmap=[init.gmap1,init.gmap2,init.gmap3]
        #set axlabel strings
        self.axlabel=list()
        for d in range(self.DIM):
            label = SpaceGridBase.axlabel_n2s[self.axtypes[d]]
            self.axlabel.append(label)
        #end for
        del self.axtypes
        for i in range(len(self.gmap)):
            self.gmap[i]=self.gmap[i].reshape((len(self.gmap[i]),))
        #end for
        return
    #end def init_from_hdfgroup


    def init_from_xmlelement(self,init):
        DIM=self.DIM
        self.axlabel=list()
        self.axgrid =list()
        #coord
        self.coord = init.coord
        #origin
        p1 = self.points[init.origin.p1]
        if 'p2' in init.origin:
            p2 = self.points[init.origin.p2]
        else:
            p2 = self.points['zero']
        #end if
        if 'fraction' in init.origin:
            frac = eval(init.origin.fraction)
        else:
            frac = 0.0
        self.origin = p1 + frac*(p2-p1)
        #axes
        self.axes = zeros((DIM,DIM))
        for d in range(DIM):
            exec 'axis=init.axis'+str(d+1)
            p1 = self.points[axis.p1]
            if 'p2' in axis:
                p2 = self.points[axis.p2]
            else:
                p2 = self.points['zero']
            #end if
            if 'scale' in axis:
                scale = eval(axis.scale)
            else:
                scale = 1.0
            #end if
            for dd in range(DIM):
                self.axes[dd,d] = scale*(p1[dd]-p2[dd])
            #end for
            self.axlabel.append(axis.label)
            self.axgrid.append(axis.grid)
        #end for
        self.initialize()
        return
    #end def init_from_xmlelement

    def initialize(self): #like qmcpack SpaceGridBase.initialize
        write=False
        succeeded=True
    
        ndomains=-1
    
        DIM = self.DIM

        coord   = self.coord
        origin  = self.origin
        axes    = self.axes
        axlabel = self.axlabel
        axgrid  = self.axgrid
        del self.axgrid


    
        ax_cartesian   = ["x" , "y"   , "z"    ] 
        ax_cylindrical = ["r" , "phi" , "z"    ] 
        ax_spherical   = ["r" , "phi" , "theta"] 
    
        cmap = dict()
        if(coord=="cartesian"):
            for d in range(DIM):
                cmap[ax_cartesian[d]]=d
                axlabel[d]=ax_cartesian[d]
            #end 
        elif(coord=="cylindrical"):
            for d in range(DIM):
                cmap[ax_cylindrical[d]]=d
                axlabel[d]=ax_cylindrical[d]
            #end 
        elif(coord=="spherical"):
            for d in range(DIM):
                cmap[ax_spherical[d]]=d
                axlabel[d]=ax_spherical[d]
            #end 
        else:
            self.error("  Coordinate supplied to spacegrid must be cartesian, cylindrical, or spherical\n  You provided "+coord,exit=False)
            succeeded=False
        #end 
        self.coordinate = SpaceGridBase.coord_s2n[self.coord]
        coordinate = self.coordinate    
    
    
        #loop over spacegrid xml elements
        naxes =DIM
        # variables for loop
        utol = 1e-5
        dimensions=zeros((DIM,),dtype=int)
        umin=zeros((DIM,))
        umax=zeros((DIM,))
        odu=zeros((DIM,))
        ndu_per_interval=[None,None,None]
        gmap=[None,None,None]
        for dd in range(DIM):
            iaxis = cmap[axlabel[dd]]
            grid = axgrid[dd]
            #read in the grid contents
            #  remove spaces inside of parentheses
            inparen=False
            gtmp=''
            for gc in grid:
                if(gc=='('):
                    inparen=True
                    gtmp+=' '
                #end 
                if(not(inparen and gc==' ')):
                    gtmp+=gc
                if(gc==')'):
                    inparen=False
                    gtmp+=' '
                #end 
            #end 
            grid=gtmp
            #  break into tokens
            tokens = grid.split()
            if(write):
                print "      grid   = ",grid
                print "      tokens = ",tokens
            #end 
            #  count the number of intervals
            nintervals=0
            for t in tokens:
                if t[0]!='(':
                    nintervals+=1
                #end 
            #end 
            nintervals-=1
            if(write):
                print "      nintervals = ",nintervals
            #  allocate temporary interval variables
            ndom_int = zeros((nintervals,),dtype=int)
            du_int = zeros((nintervals,))
            ndu_int = zeros((nintervals,),dtype=int)
            #  determine number of domains in each interval and the width of each domain
            u1=1.0*eval(tokens[0])
            umin[iaxis]=u1
            if(abs(u1)>1.0000001):
                self.error("  interval endpoints cannot be greater than 1\n  endpoint provided: "+str(u1),exit=False)
                succeeded=False
            #end 
            is_int=False
            has_paren_val=False
            interval=-1
            for i in range(1,len(tokens)):
                if not tokens[i].startswith('('):
                    u2=1.0*eval(tokens[i])
                    umax[iaxis]=u2
                    if(not has_paren_val):
                        du_i=u2-u1
                    #end
                    has_paren_val=False
                    interval+=1
                    if(write):
                        print "      parsing interval ",interval," of ",nintervals
                        print "      u1,u2 = ",u1,",",u2
                    #end 
                    if(u2<u1):
                        self.error("  interval ("+str(u1)+","+str(u2)+") is negative",exit=False)
                        succeeded=False
                    #end 
                    if(abs(u2)>1.0000001):
                        self.error("  interval endpoints cannot be greater than 1\n  endpoint provided: "+str(u2),exit=False)
                        succeeded=False
                    #end 
                    if(is_int):
                        du_int[interval]=(u2-u1)/ndom_i
                        ndom_int[interval]=ndom_i
                    else:
                        du_int[interval]=du_i
                        ndom_int[interval]=floor((u2-u1)/du_i+.5)
                        if(abs(u2-u1-du_i*ndom_int[interval])>utol):
                            self.error("  interval ("+str(u1)+","+str(u2)+") not divisible by du="+str(du_i),exit=False)
                            succeeded=False
                        #end 
                    #end 
                    u1=u2
                else:
                    has_paren_val=True
                    paren_val=tokens[i][1:len(tokens[i])-1]
                    if(write):
                        print "      interval spacer = ",paren_val
                    #end if
                    is_int=tokens[i].find(".")==-1
                    if(is_int):
                        ndom_i = eval(paren_val)
                        du_i = -1.0
                    else:
                        ndom_i = 0
                        du_i = eval(paren_val)
                    #end 
                #end 
            #end 
            # find the smallest domain width
            du_min=min(du_int)
            odu[iaxis]=1.0/du_min
            # make sure it divides into all other domain widths
            for i in range(len(du_int)):
                ndu_int[i]=floor(du_int[i]/du_min+.5)
                if(abs(du_int[i]-ndu_int[i]*du_min)>utol):
                    self.error("interval {0} of axis {1} is not divisible by smallest subinterval {2}".format(i+1,iaxis+1,du_min),exit=False)
                    succeeded=False
                #end 
            #end      
    
            if(write):
                print "      interval breakdown"
                print "        interval,ndomains,nsubdomains_per_domain"
                for i in range(len(ndom_int)):
                    print "      ",i,",",ndom_int[i],",",ndu_int[i]
                #end 
            #end 
       
            # set up the interval map such that gmap[u/du]==domain index
            gmap[iaxis] = zeros((floor((umax[iaxis]-umin[iaxis])*odu[iaxis]+.5),),dtype=int)
            n=0
            nd=-1
            if(write):
                print "        i,j,k    ax,n,nd  "
            #end if
            for i in range(len(ndom_int)):
                for j in range(ndom_int[i]):
                    nd+=1
                    for k in range(ndu_int[i]):
                        gmap[iaxis][n]=nd
                        if(write):
                            print "      ",i,",",j,",",k,"    ",iaxis,",",n,",",nd
                        #end
                        n+=1
                    #end 
                #end 
            #end 
            dimensions[iaxis]=nd+1
            #end read in the grid contents
            
            #save interval width information
            ndom_tot=sum(ndom_int)
            ndu_per_interval[iaxis] = zeros((ndom_tot,),dtype=int)
            idom=0
            for i in range(len(ndom_int)):
                for ii in range(ndom_int[i]):
                    ndu_per_interval[iaxis][idom] = ndu_int[i]
                    idom+=1
                #end 
          #end       
        #end 

        axinv = inv(axes)
    
        #check that all axis grid values fall in the allowed intervals
        cartmap = dict()
        for d in range(DIM):
            cartmap[ax_cartesian[d]]=d
        #end for
        for d in range(DIM):
            if axlabel[d] in cartmap:
                if(umin[d]<-1.0 or umax[d]>1.0):
                    self.error("  grid values for {0} must fall in [-1,1]\n".format(axlabel[d])+"  interval provided: [{0},{1}]".format(umin[d],umax[d]),exit=False)
                    succeeded=False
                #end if
            elif(axlabel[d]=="phi"):
                if(abs(umin[d])+abs(umax[d])>1.0):
                    self.error("  phi interval cannot be longer than 1\n  interval length provided: {0}".format(abs(umin[d])+abs(umax[d])),exit=False)
                    succeeded=False
                #end if
            else:
                if(umin[d]<0.0 or umax[d]>1.0):
                    self.error("  grid values for {0} must fall in [0,1]\n".format(axlabel[d])+"  interval provided: [{0},{1}]".format(umin[d],umax[d]),exit=False)
                    succeeded=False
                #end if
            #end if
        #end for
    
    
        #set grid dimensions
        # C/Python style indexing
        dm=array([0,0,0],dtype=int)
        dm[0] = dimensions[1]*dimensions[2]
        dm[1] = dimensions[2]
        dm[2] = 1
    
        ndomains=prod(dimensions)
    
        volume = abs(det(axes))*8.0#axes span only one octant
    
        #compute domain volumes, centers, and widths
        domain_volumes = zeros((ndomains,))
        domain_centers = zeros((ndomains,DIM))
        domain_uwidths = zeros((ndomains,DIM))
        interval_centers = [None,None,None]
        interval_widths  = [None,None,None]
        for d in range(DIM):
            nintervals = len(ndu_per_interval[d])
            interval_centers[d] = zeros((nintervals))
            interval_widths[d] = zeros((nintervals))
            interval_widths[d][0]=ndu_per_interval[d][0]/odu[d]
            interval_centers[d][0]=interval_widths[d][0]/2.0+umin[d]
            for i in range(1,nintervals):
                interval_widths[d][i] = ndu_per_interval[d][i]/odu[d]
                interval_centers[d][i] = interval_centers[d][i-1] \
                    +.5*(interval_widths[d][i]+interval_widths[d][i-1])
            #end for
        #end for
        du,uc,ubc,rc = zeros((DIM,)),zeros((DIM,)),zeros((DIM,)),zeros((DIM,))
        vol = -1e99
        vol_tot=0.0
        vscale = abs(det(axes))
        
        for i in range(dimensions[0]):                           
            for j in range(dimensions[1]):                           
                for k in range(dimensions[2]):                           
                    idomain = dm[0]*i + dm[1]*j + dm[2]*k
                    du[0] = interval_widths[0][i]
                    du[1] = interval_widths[1][j]
                    du[2] = interval_widths[2][k]
                    uc[0] = interval_centers[0][i]
                    uc[1] = interval_centers[1][j]
                    uc[2] = interval_centers[2][k]
    
                    if(coordinate==SpaceGridBase.cartesian):
                        vol=du[0]*du[1]*du[2]
                        ubc=uc
                    elif(coordinate==SpaceGridBase.cylindrical):
                        uc[1]=2.0*pi*uc[1]-pi
                        du[1]=2.0*pi*du[1]
                        vol=uc[0]*du[0]*du[1]*du[2]
                        ubc[0]=uc[0]*cos(uc[1])
                        ubc[1]=uc[0]*sin(uc[1])
                        ubc[2]=uc[2]
                    elif(coordinate==SpaceGridBase.spherical):
                        uc[1]=2.0*pi*uc[1]-pi
                        du[1]=2.0*pi*du[1]
                        uc[2]=    pi*uc[2]
                        du[2]=    pi*du[2]
                        vol=(uc[0]*uc[0]+du[0]*du[0]/12.0)*du[0] \
                           *du[1]                                \
                           *2.0*sin(uc[2])*sin(.5*du[2])          
                        ubc[0]=uc[0]*sin(uc[2])*cos(uc[1])
                        ubc[1]=uc[0]*sin(uc[2])*sin(uc[1])
                        ubc[2]=uc[0]*cos(uc[2])
                    #end if
                    vol*=vscale
    
                    vol_tot+=vol
    
                    rc = dot(axes,ubc) + origin
    
                    domain_volumes[idomain] = vol
                    for d in range(DIM):
                        domain_uwidths[idomain,d] = du[d]
                        domain_centers[idomain,d] = rc[d]
                    #end for
                #end for
            #end for
        #end for
    
        #find the actual volume of the grid
        du = umax-umin
        uc = .5*(umax+umin)
        if coordinate==SpaceGridBase.cartesian:
            vol=du[0]*du[1]*du[2]
        elif coordinate==SpaceGridBase.cylindrical:
            uc[1]=2.0*pi*uc[1]-pi
            du[1]=2.0*pi*du[1]
            vol=uc[0]*du[0]*du[1]*du[2]
        elif coordinate==SpaceGridBase.spherical:
            uc[1]=2.0*pi*uc[1]-pi
            du[1]=2.0*pi*du[1]
            uc[2]=    pi*uc[2]
            du[2]=    pi*du[2]
            vol=(uc[0]*uc[0]+du[0]*du[0]/12.0)*du[0]*du[1]*2.0*sin(uc[2])*sin(.5*du[2])
        #end if
        volume = vol*abs(det(axes))

        for q in SpaceGridBase.quantities:
            self[q].mean  = zeros((ndomains,))
            self[q].error = zeros((ndomains,))
        #end for

        #save the results
        self.axinv              = axinv         
        self.volume             = volume        
        self.gmap               = gmap          
        self.umin               = umin          
        self.umax               = umax      
        self.odu                = odu
        self.dm                 = dm            
        self.dimensions         = dimensions
        self.ndomains           = ndomains      
        self.domain_volumes     = domain_volumes
        self.domain_centers     = domain_centers
        self.domain_uwidths     = domain_uwidths


        #succeeded = succeeded and check_grid()
    
        if(self.init_exit_fail and not succeeded):
            self.error(" in def initialize")
        #end 

        return succeeded
    #end def initialize

    def point2unit_cartesian(point):
        u = dot(self.axinv,(point-self.origin)) 
        return u
    #end def point2unit_cartesian

    def point2unit_cylindrical(point):
        ub = dot(self.axinv,(point-self.origin)) 
        u=zeros((self.DIM,))
        u[0] = sqrt(ub[0]*ub[0]+ub[1]*ub[1]) 
        u[1] = atan2(ub[1],ub[0])*o2pi+.5 
        u[2] = ub[2] 
        return u
    #end def point2unit_cylindrical

    def point2unit_spherical(point):
        ub = dot(self.axinv,(point-self.origin)) 
        u=zeros((self.DIM,))
        u[0] = sqrt(ub[0]*ub[0]+ub[1]*ub[1]+ub[2]*ub[2]) 
        u[1] = atan2(ub[1],ub[0])*o2pi+.5 
        u[2] = acos(ub[2]/u[0])*o2pi*2.0 
        return u
    #end def point2unit_spherical

    def points2domains_cartesian(self,points,domains,points_outside):        
        u  = zeros((self.DIM,))
        iu = zeros((self.DIM,),dtype=int)
        ndomains=-1
        npoints,ndim = points.shape
        for p in xrange(npoints):
            u = dot(self.axinv,(points[p]-self.origin)) 
            if (u>self.umin).all() and (u<self.umax).all():
                points_outside[p]=False 
                iu=floor( (u-self.umin)*self.odu )
                iu[0] = self.gmap[0][iu[0]]
                iu[1] = self.gmap[1][iu[1]]
                iu[2] = self.gmap[2][iu[2]]
                ndomains+=1
                domains[ndomains,0] = p
                domains[ndomains,1] = dot(self.dm,iu)
            #end 
        #end 
        ndomains+=1
        return ndomains 
    #end def points2domains_cartesian

    def points2domains_cylindrical(self,points,domains,points_outside):        
        u  = zeros((self.DIM,))
        iu = zeros((self.DIM,),dtype=int)
        ndomains=-1
        npoints,ndim = points.shape
        for p in xrange(npoints):
            ub = dot(self.axinv,(points[p]-self.origin)) 
            u[0] = sqrt(ub[0]*ub[0]+ub[1]*ub[1]) 
            u[1] = atan2(ub[1],ub[0])*o2pi+.5 
            u[2] = ub[2] 
            if (u>self.umin).all() and (u<self.umax).all():
                points_outside[p]=False 
                iu=floor( (u-self.umin)*self.odu )
                iu[0] = self.gmap[0][iu[0]]
                iu[1] = self.gmap[1][iu[1]]
                iu[2] = self.gmap[2][iu[2]]
                ndomains+=1
                domains[ndomains,0] = p
                domains[ndomains,1] = dot(self.dm,iu)
            #end 
        #end 
        ndomains+=1
        return ndomains 
    #end def points2domains_cylindrical

    def points2domains_spherical(self,points,domains,points_outside):        
        u  = zeros((self.DIM,))
        iu = zeros((self.DIM,),dtype=int)
        ndomains=-1
        npoints,ndim = points.shape
        for p in xrange(npoints):
            ub = dot(self.axinv,(points[p]-self.origin)) 
            u[0] = sqrt(ub[0]*ub[0]+ub[1]*ub[1]+ub[2]*ub[2]) 
            u[1] = atan2(ub[1],ub[0])*o2pi+.5 
            u[2] = acos(ub[2]/u[0])*o2pi*2.0 
            if (u>self.umin).all() and (u<self.umax).all():
                points_outside[p]=False 
                iu=floor( (u-self.umin)*self.odu )
                iu[0] = self.gmap[0][iu[0]]
                iu[1] = self.gmap[1][iu[1]]
                iu[2] = self.gmap[2][iu[2]]
                ndomains+=1
                domains[ndomains,0] = p
                domains[ndomains,1] = dot(self.dm,iu)
            #end 
        #end 
        ndomains+=1
        return ndomains 
    #end def points2domains_spherical

    
    def shift_origin(self,shift):
        self.origin += shift
        for i in range(self.domain_centers.shape[0]):
            self.domain_centers[i,:] += shift
        #end for
        return
    #end def shift_origin


    def set_origin(self,origin):
        self.shift_origin(origin-self.origin)
        return
    #end def set_origin


    def interpolate_across(self,quantities,spacegrids,outside,integration=False,warn=False):
        #if the grid is to be used for integration confirm that domains 
        #  of this spacegrid subdivide source spacegrid domains
        if integration:
            #setup checking variables
            am_cartesian   = self.coordinate==Spacegrid.cartesian
            am_cylindrical = self.coordinate==Spacegrid.cylindrical
            am_spherical   = self.coordinate==Spacegrid.spherical
            fine_interval_centers = [None,None,None]
            fine_interval_domains = [None,None,None]
            for d in range(self.DIM):
                ndu = round( (self.umax[d]-self.umin[d])*self.odu[d] )
                if len(self.gmap[d])!=ndu:
                    self.error('ndu is different than len(gmap)')
                #end if
                du = 1./self.odu[d]
                fine_interval_centers[d] = self.umin + .5*du + du*array(range(ndu))
                find_interval_domains[d] = zeros((ndu,))
            #end for
            #checks are done on each source spacegrid to determine interpolation compatibility 
            for s in spacegrids:
                # all the spacegrids must have coordinate system to satisfy this
                if s.coordinate!=self.coordinate:
                    if warn:
                        self.warn('SpaceGrids must have same coordinate for interpolation')
                    #end if
                    return False
                #end if
                # each spacegrids' axes must be int mult of this spacegrid's axes
                #   (this ensures that isosurface shapes conform)
                tile = dot(self.axinv,s.axes)
                for d in range(self.DIM):
                    if not is_integer(tile[d,d]):
                        if warn:
                            self.warn("source axes must be multiples of interpolant's axes")
                        #end if
                        return False
                    #end if
                #end for
                # origin must be at r=0 for cylindrical or spherical
                uo = self.point2unit(s.origin)
                if am_cylindrical or am_spherical:
                    if uo[0]>1e-6:
                        if warn:
                            self.warn('source origin must lie at interpolant r=0')
                        #end if
                        return False
                    #end if
                #end if
                # fine meshes must align
                #  origin must be an integer multiple of smallest dom width
                if am_cylindrical:
                    mdims=[2]
                elif am_cartesian:
                    mdims=[0,1,2]
                else:
                    mdims=[]
                #end if
                for d in mdims:
                    if not is_integer(uo[d]*self.odu[d]):
                        if warn:
                            self.warn('source origin does not lie on interpolant fine mesh')
                        #end if
                        return False
                    #end if
                #end for
                #  smallest dom width must be multiple of this smallest dom width 
                for d in range(self.DIM):
                    if not is_integer(self.odu[d]/s.odu[d]):
                        if warn:
                            self.warn('smallest source domain width must be a multiple of interpolants smallest domain width')
                        #end if
                        return False
                    #end if
                #end for
                #  each interval along each direction for interpolant must map to only one source interval
                #    construct points at each fine interval center of interpolant, run them through source gmap to get interval indices
                for d in range(self.DIM):
                    fine_interval_domains[d][:]=-2
                    gmlen = len(s.gmap[d])
                    for i in range(len(fine_interval_centers[d])):
                        uc = fine_interval_centers[d][i]
                        ind = floor((uc-s.umin[d])*s.odu[d])
                        if ind < gmlen:
                            idom=s.gmap[d][ind]
                        else:
                            idom=-1
                        #end if
                        fine_interval_domains[d][i]=idom
                    #end for
                    cind   = self.gmap[d][0]
                    istart = 0
                    iend   = 0
                    for i in range(len(self.gmap[d])):
                        if self.gmap[d][i]==cind:
                            iend+=1
                        else:
                            source_ind = fine_interval_domains[istart]
                            for j in range(istart+1,iend):
                                if fine_interval_domains[j]!=source_ind:
                                    if warn:
                                        self.warn('an interpolant domain must not fall on multiple source domains')
                                    #end if
                                    return False
                                #end if
                            #end for
                            istart=iend
                        #end if
                    #end for
                #end for                    
            #end for
        #end if


        #get the list of domains points from this grid fall in
        #  and interpolate requested quantities on them
        domain_centers  = self.domain_centers
        domind = zeros((self.ndomains,2),dtype=int)
        domout = ones((self.ndomains,) ,dtype=int)
        for s in spacegrids:
            domind[:,:] = -1
            ndomin = s.points2domains(domain_centers,domind,domout)
            for q in quantities:
                self[q].mean[domind[0:ndomin,0]]  = s[q].mean[domind[0:ndomin,1]].copy()
                self[q].error[domind[0:ndomin,0]] = s[q].error[domind[0:ndomin,1]].copy()
            #end for
        #end for
        for d in xrange(self.ndomains):
            if domout[d]:
                for q in quantities:
                    self[q].mean[d]  = outside[q].mean
                    self[q].error[d] = outside[q].error
                #end for
            #end if
        #end for
        return True
    #end def interpolate_across


    def interpolate(self,points,quantities=None):
        if quantities==None:
            quantities=SpaceGridBase.quantities
        #end if
        npoints,ndim = points.shape
        ind = empty((npoints,2),dtype=int)
        out = ones((npoints,) ,dtype=int)
        nin = self.points2domains(points,ind,out)
        result = QAobject()
        for q in quantities:
            result._add_attribute(q,QAobject())
            result[q].mean  = zeros((npoints,))
            result[q].error = zeros((npoints,))
            result[q].mean[ind[0:nin,0]]  = self[q].mean[ind[0:nin,1]].copy()
            result[q].error[ind[0:nin,0]] = self[q].error[ind[0:nin,1]].copy()
        #end for
        return result
    #end def interpolate


    def isosurface(self,quantity,contours=5,origin=None):
        if quantity not in SpaceGridBase.quantities:
            self.error()
        #end if
        dimensions = self.dimensions
        if origin==None:
            points     = self.domain_centers
        else:
            npoints,ndim = self.domain_centers.shape
            points = empty((npoints,ndim))
            for i in range(npoints):
                points[i,:] = origin + self.domain_centers[i,:]
            #end for
        #end if
        scalars    = self[quantity].mean
        name       = quantity
        self.plotter.isosurface(points,scalars,contours,dimensions,name)
        return 
    #end def isosurface

    
    def surface_slice(self,quantity,x,y,z,options=None):
        if quantity not in SpaceGridBase.quantities:
            self.error()
        #end if
        points = empty( (x.size,self.DIM) )
        points[:,0] = x.ravel()
        points[:,1] = y.ravel()
        points[:,2] = z.ravel()
        val = self.interpolate(points,[quantity])
        scalars = val[quantity].mean
        scalars.shape = x.shape
        self.plotter.surface_slice(x,y,z,scalars,options)
        return
    #end def surface_slice


    def plot_axes(self,color=None,radius=.025,origin=None):
        if color is None:
            color = (0.,0,0)
        #end if
        if origin is None:
            origin = array([0.,0,0])
        #end if
        colors=array([[1.,0,0],[0,1.,0],[0,0,1.]])
        for d in range(self.DIM):
            a=self.axes[:,d]+origin
            ax=array([-a[0],a[0]])
            ay=array([-a[1],a[1]])
            az=array([-a[2],a[2]])
            self.plotter.plot3d(ax,ay,az,tube_radius=radius,color=tuple(colors[:,d]))
        #end for
        return
    #end def plot_axes

    def plot_box(self,color=None,radius=.025,origin=None):
        if color is None:
            color = (0.,0,0)
        #end if
        if origin is None:
            origin = array([0.,0,0])
        #end if
        p = self.points
        p1=p.cmmm+origin
        p2=p.cmpm+origin
        p3=p.cpmm+origin
        p4=p.cppm+origin
        p5=p.cmmp+origin
        p6=p.cmpp+origin
        p7=p.cpmp+origin
        p8=p.cppp+origin
        bline = array([p1,p2,p4,p3,p1,p5,p6,p8,p7,p5,p7,p3,p4,p8,p6,p2])
        self.plotter.plot3d(bline[:,0],bline[:,1],bline[:,2],color=color)
        return
    #end def plot_box
#end class RectilinearGrid





class VoronoiGridInitializer(SpaceGridInitializer):
    def __init__(self):
        SpaceGridInitializer.__init__(self)
    #end def __init__
#end class VoronoiGridInitializer


class VoronoiGrid(SpaceGridBase):
    def __init__(self,initobj=None,options=None):
        SpaceGridBase.__init__(self,initobj,options)
        return
    #end def __init__

    def copy(self,other):
        return VoronoiGrid(other)
    #end def copy


    def reorder_atomic_data(self,imap):
        for q in self.quantities:
            qv = self[q]
            qv.mean  = qv.mean[...,imap]
            qv.error = qv.error[...,imap]
        #end for
        if 'data' in self:
            data = self.data
            for q in self.quantities:
                data[q] = data[q][...,imap,:]
            #end for
        #end if
    #end def reorder_atomic_data
#end class VoronoiGrid






def SpaceGrid(init,opts=None):
    SpaceGrid.count+=1
    
    iname = init.__class__.__name__
    if iname=='HDFgroup':
        coordinate = init.coordinate[0]
    #end if
    coord = SpaceGrid.coord_n2s[coordinate]

    if coord in SpaceGrid.rect:
        return RectilinearGrid(init,opts)
    elif coord=='voronoi':
        return VoronoiGrid(init,opts)
    else:
        print 'SpaceGrid '+coord+' has not been implemented, exiting...'
        exit()
    #end if

#end def SpaceGrid
SpaceGrid.count = 0
SpaceGrid.coord_n2s = SpaceGridBase.coord_n2s
SpaceGrid.rect = set(['cartesian','cylindrical','spherical'])
