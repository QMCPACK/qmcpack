
import os
import re
from numpy import array,zeros,dot,loadtxt,floor,empty,sqrt
from extended_numpy import ndgrid,simstats,simplestats,equilibration_length
from generic import obj
from qaobject import QAobject
from qmcpack_analyzer_base import QAanalyzer,QAdata,QAHDFdata


class QuantityAnalyzer(QAanalyzer):
    def __init__(self):
        QAanalyzer.__init__(self)
        self.method_info = QAanalyzer.method_info
    #end def __init__

    def plot_trace(self,quantity,*args,**kwargs):
        from matplotlib.pyplot import plot,xlabel,ylabel,title,ylim
        if 'data' in self:
            if not quantity in self.data:
                self.error('quantity '+quantity+' is not present in the data')
            #end if
            nbe = self.get_nblocks_exclude()
            q = self.data[quantity]
            middle = int(len(q)/2)
            qmean = q[middle:].mean()
            qmax = q[middle:].max()
            qmin = q[middle:].min()
            ylims = [qmean-2*(qmean-qmin),qmean+2*(qmax-qmean)]
            smean,svar = self[quantity].tuple('mean','sample_variance')
            sstd = sqrt(svar)
            plot(q,*args,**kwargs)
            plot([nbe,nbe],ylims,'k-.',lw=2)
            plot([0,len(q)],[smean,smean],'r-')
            plot([0,len(q)],[smean+sstd,smean+sstd],'r-.')
            plot([0,len(q)],[smean-sstd,smean-sstd],'r-.')
            ylim(ylims)
            ylabel(quantity)
            xlabel('samples')
            title('Trace of '+quantity)
        #end if
    #end def QuantityAnalyzer

    def init_sub_analyzers(self):
        None
    #end def init_sub_analyzers
    
    def get_nblocks_exclude(self):
        return self.info.nblocks_exclude
    #end def get_nblocks_exclude
#end class QuantityAnalyzer


class DatAnalyzer(QuantityAnalyzer):
    def __init__(self,filepath=None,equilibration=None):
        QuantityAnalyzer.__init__(self)
        self.info.filepath = filepath
        if equilibration!=None:
            self.load_data()
            nbe = equilibration_length(self.data[equilibration])
            self.method_info.nblocks_exclude = nbe
        #end if
    #end def __init__

    def analyze_local(self):
        self.not_implemented()
    #end def load_data_local
#end class ScalarsDatAnalyzer


class ScalarsDatAnalyzer(DatAnalyzer):
    def load_data_local(self):
        filepath = self.info.filepath
        quantities = QAanalyzer.request.quantities

        lt = loadtxt(filepath)
        if len(lt.shape)==1:
            lt.shape = (1,len(lt))
        #end if

        data = lt[:,1:].transpose()

        fobj = open(filepath,'r')
        variables = fobj.readline().split()[2:]
        fobj.close()

        self.data = QAdata()
        for i in range(len(variables)):
            var = variables[i]
            cvar = self.condense_name(var)
            if cvar in quantities:
                self.data[var]=data[i,:]
            #end if
        #end for
    #end def load_data_local


    def analyze_local(self):
        nbe = QAanalyzer.method_info.nblocks_exclude
        self.info.nblocks_exclude = nbe
        data = self.data
        for varname,samples in data.iteritems():
            (mean,var,error,kappa)=simstats(samples[nbe:])
            self[varname] = obj(
                mean            = mean,
                sample_variance = var,
                error           = error,
                kappa           = kappa
                )
        #end for
        
        if 'LocalEnergy_sq' in data:
            v = data.LocalEnergy_sq - data.LocalEnergy**2
            (mean,var,error,kappa)=simstats(v[nbe:])
            self.LocalEnergyVariance = obj(
                mean            = mean,
                sample_variance = var,
                error           = error,
                kappa           = kappa
                )
        #end if            
    #end def load_data_local
#end class ScalarsDatAnalyzer


class DmcDatAnalyzer(DatAnalyzer):
    def load_data_local(self):
        filepath = self.info.filepath

        lt = loadtxt(filepath)
        if len(lt.shape)==1:
            lt.shape = (1,len(lt))
        #end if

        data = lt[:,1:].transpose()

        fobj = open(filepath,'r')
        variables = fobj.readline().split()[2:]
        fobj.close()

        self.data = QAdata()
        for i in range(len(variables)):
            var = variables[i]
            self.data[var]=data[i,:]
        #end for
    #end def load_data_local


    def analyze_local(self):
        nbe = QAanalyzer.method_info.nblocks_exclude
        self.info.nblocks_exclude = nbe
        data = self.data

        input       = self.run_info.input
        series      = self.method_info.series
        ndmc_blocks = self.run_info.request.ndmc_blocks

        qmc    = input.simulation.calculations[series]
        blocks = qmc.blocks
        steps  = qmc.steps
        nse    = nbe*steps

        self.info.nsteps_exclude = nse

        nsteps = len(data.list()[0])-nse

        #nsteps = blocks*steps-nse
        block_avg = nsteps > 2*ndmc_blocks
        if block_avg:
            block_size  = int(floor(float(nsteps)/ndmc_blocks))
            ndmc_blocks = int(floor(float(nsteps)/block_size))
            nse += nsteps-ndmc_blocks*block_size
            nsteps      = ndmc_blocks*block_size
        #end if

        for varname,samples in data.iteritems():
            samp = samples[nse:]
            if block_avg:
                samp.shape = ndmc_blocks,block_size
                samp = samp.mean(axis=1)
            #end if
            (mean,var,error,kappa)=simstats(samp)
            self[varname] = obj(
                mean            = mean,
                sample_variance = var,
                error           = error,
                kappa           = kappa
                )
        #end for
    #end def load_data_local

    
    def get_nblocks_exclude(self):
        return self.info.nsteps_exclude
    #end def get_nblocks_exclude
#end class DmcDatAnalyzer


class HDFAnalyzer(QuantityAnalyzer):
    def __init__(self):
        QuantityAnalyzer.__init__(self)
        self.info.should_remove = False
    #end def __init__
#end class HDFAnalyzer


class ScalarsHDFAnalyzer(HDFAnalyzer):
    corrections = obj(
        mpc = obj(ElecElec=-1,MPC=1),
        kc  = obj(KEcorr=1)
        )
    
    def __init__(self,exclude):
        HDFAnalyzer.__init__(self)
        self.info.exclude = exclude
    #end def


    def load_data_local(self,data=None):
        if data==None:
            self.error('attempted load without data')
        #end if
        exclude = self.info.exclude
        self.data = QAHDFdata()
        for var in data.keys():
            if not var in exclude and not str(var)[0]=='_':
                self.data[var] = data[var]
                del data[var]
            #end if
        #end for
        corrvars = ['LocalEnergy','ElecElec','MPC','KEcorr']
        if set(corrvars)<set(self.data.keys()):
            Ed,Ved,Vmd,Kcd = self.data.tuple(*corrvars)
            E,E2 = Ed.value,Ed.value_squared
            Ve,Ve2 = Ved.value,Ved.value_squared
            Vm,Vm2 = Vmd.value,Vmd.value_squared
            Kc,Kc2 = Kcd.value,Kcd.value_squared
            self.data.LocalEnergy_mpc_kc = obj(
                value = E-Ve+Vm+Kc,
                value_squared = E2+Ve2+Vm2+Kc2 + 2*(E*(-Ve+Vm+Kc)-Ve*(Vm+Kc)+Vm*Kc)
                )
        #end if
    #end def load_data_local


    def analyze_local(self):
        nbe = QAanalyzer.method_info.nblocks_exclude
        self.info.nblocks_exclude = nbe
        for varname,val in self.data.iteritems():
            (mean,var,error,kappa)=simstats(val.value[nbe:,...].ravel())
            self[varname] = obj(
                mean            = mean,
                variance        = val.value_squared[nbe:,...].mean()-mean**2,
                sample_variance = var,
                error           = error,
                kappa           = kappa
                )
        #end for
        self.correct('mpc','kc')
    #end def analyze_local


    def correct(self,*corrections):
        corrkey=''
        for corr in corrections:
            corrkey+=corr+'_'
        #end for
        corrkey=corrkey[:-1]
        if set(corrections)>set(self.corrections.keys()):
            self.warn('correction '+corrkey+' is unknown and cannot be applied')
            return
        #end if
        if not 'data' in self:
            self.warn('correction '+corrkey+' cannot be applied because data is not present')
            return
        #end if
        varname = 'LocalEnergy_'+corrkey
        if varname in self and varname in self.data:
            return
        #end if
        corrvars = ['LocalEnergy']
        signs    = [1]
        for corr in corrections:
            for var,sign in self.corrections[corr].iteritems():
                corrvars.append(var)
                signs.append(sign)
            #end for
        #end for
        missing = list(set(corrvars)-set(self.data.keys()))
        if len(missing)>0:            
            self.warn('correction '+corrkey+' cannot be applied because '+str(missing)+' are missing')
            return
        #end if

        le = self.data.LocalEnergy
        E,E2 = 0*le.value,0*le.value_squared
        n = len(corrvars)
        for i in range(n):
            ed = self.data[corrvars[i]]
            e,e2 = ed.value,ed.value_squared
            s = signs[i]
            E += s*e
            E2 += e2
            for j in range(i+1,n):
                eo = self.data[corrvars[j]].value
                so = signs[j]
                E2 += 2*s*e*so*eo
            #end for
        #end for
        val = obj(value=E,value_squared=E2)
        self.data[varname] = val
        nbe = self.info.nblocks_exclude
        (mean,var,error,kappa)=simstats(val.value[nbe:,...].ravel())
        self[varname] = obj(
            mean            = mean,
            variance        = val.value_squared[nbe:,...].mean()-mean**2,
            sample_variance = var,
            error           = error,
            kappa           = kappa
            )
    #end def correct
#end class ScalarsHDFAnalyzer



class EnergyDensityAnalyzer(HDFAnalyzer):
    def __init__(self,name):
        HDFAnalyzer.__init__(self)
        self.info.set(
            name = name,
            reordered = False
            )
    #end def __init__


    def load_data_local(self,data=None):
        if data==None:
            self.error('attempted load without data')
        #end if
        name = self.info.name
        self.data = QAdata()
        if name in data:
            self.data = data[name]
            del data[name]
        else:
            self.info.should_remove = True
        #end if
    #end def load_data_local


    def analyze_local(self):
        from spacegrid import SpaceGrid

        nbe = QAanalyzer.method_info.nblocks_exclude
        self.info.nblocks_exclude = nbe 
        data = self.data

        #why is this called 3 times?
        print nbe

        #transfer hdf data
        sg_pattern = re.compile(r'spacegrid\d*')
        nspacegrids=0
        #  add simple data first
        for k,v in data._iteritems():
            if not sg_pattern.match(k):
                self._add_attribute(k,v)
            else:
                nspacegrids+=1
            #end if
        #end for
        #  add spacegrids second
        opts = QAobject()
        opts.points = self.reference_points
        opts.nblocks_exclude = nbe
        self.spacegrids=[]
        if nspacegrids==0:
            self.spacegrids.append(SpaceGrid(data.spacegrid,opts))
        else:
            for ig in range(nspacegrids):
                sg=SpaceGrid(data['spacegrid'+str(ig+1)],opts)
                self.spacegrids.append(sg)
            #end for
        #end if

        #reorder atomic data to match input file for Voronoi grids
        if self.run_info.type=='bundled':
            self.info.reordered=True
        #end if
        if not self.info.reordered:
            self.reorder_atomic_data()
        #end if

        #convert quantities outside all spacegrids
        outside = QAobject()
        iD,iT,iV = tuple(range(3))        
        outside.D  = QAobject()
        outside.T  = QAobject()
        outside.V  = QAobject()
        outside.E  = QAobject()
        outside.P  = QAobject()

        value = self.outside.value.transpose()[...,nbe:]

        #mean,error = simplestats(value)
        mean,var,error,kappa = simstats(value)
        outside.D.mean   = mean[iD]
        outside.D.error  = error[iD]
        outside.T.mean   = mean[iT]
        outside.T.error  = error[iT]
        outside.V.mean   = mean[iV]
        outside.V.error  = error[iV]

        E  = value[iT,:]+value[iV,:]
        #mean,error = simplestats(E)
        mean,var,error,kappa = simstats(E)
        outside.E.mean  = mean
        outside.E.error = error

        P  = 2./3.*value[iT,:]+1./3.*value[iV,:]
        #mean,error = simplestats(P)
        mean,var,error,kappa = simstats(P)
        outside.P.mean  = mean
        outside.P.error = error

        self.outside = outside

        self.outside.data = obj(
            D = value[iD,:],
            T = value[iT,:],
            V = value[iV,:],
            E = E,
            P = P
            )

        return
    #end def analyze_local

    
    def reorder_atomic_data(self):
        input = self.run_info.input
        xml   = self.run_info.ordered_input
        ps = input.get('particlesets')
        if 'ion0' in ps and len(ps.ion0.groups)>1:
            qsx = xml.simulation.qmcsystem
            if len(ps)==1:
                psx = qsx.particleset
            else:
                psx=None
                for pst in qsx.particleset:
                    if pst.name=='ion0':
                        psx=pst
                    #end if
                #end for
                if psx==None:
                    self.error('ion0 particleset not found in qmcpack xml file for atomic reordering of Voronoi energy density')
                #end if
            #end if

            #ordered ion names
            # xml groups are ordered the same as in qmcpack's input file
            ion_names = []
            for gx in psx.group:
                ion_names.append(gx.name)
            #end for

            #create the mapping to restore proper ordering
            nions = ps.ion0.size
            ions = ps.ion0.ionid
            imap=empty((nions,),dtype=int)
            icurr = 0
            for ion_name in ion_names:
                for i in range(len(ions)):
                    if ions[i]==ion_name:
                        imap[i]=icurr
                        icurr+=1
                    #end if
                #end for
            #end for

            #reorder the atomic data
            for sg in self.spacegrids:

                #print 'qqa'
                #import code
                #code.interact(local=locals())

                sg.reorder_atomic_data(imap)
            #end for
        #end if
        self.info.reordered=True
        return
    #end def reorder_atomic_data


    def remove_data(self):
        QAanalyzer.remove_data(self)
        if 'spacegrids' in self:
            for sg in self.spacegrids:
                if 'data' in sg:
                    del sg.data
                #end if
            #end for
        #end if
        if 'outside' in self and 'data' in self.outside:
            del self.outside.data
        #end if
    #end def remove_data

    
    def prev_init(self):
        if data._contains_group("spacegrid1"):
            self.points = data.spacegrid1.domain_centers
            self.axinv  = data.spacegrid1.axinv
            val = data.spacegrid1.value
            npoints,ndim = self.points.shape
            self.E = zeros((npoints,))
            print 'p shape ',self.points.shape
            print 'v shape ',val.shape
            nblocks,nvpoints = val.shape
            for b in range(nblocks):
                for i in range(npoints):
                    ind = 6*i
                    self.E[i] += val[b,ind+1] + val[b,ind+2]
                #end for
            #end for
        #end if
    #end def prev_init

    

    def isosurface(self):
        from enthought.mayavi import mlab

        npoints,ndim = self.points.shape
        dimensions = array([20,20,20])

        x = zeros(dimensions)
        y = zeros(dimensions)
        z = zeros(dimensions)
        s = zeros(dimensions)

        ipoint = 0
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                for k in range(dimensions[2]):
                    r = self.points[ipoint,:]
                    u = dot(self.axinv,r)
                    #u=r
                    x[i,j,k] = u[0]
                    y[i,j,k] = u[1]
                    z[i,j,k] = u[2]
                    s[i,j,k] = self.E[ipoint]
                    ipoint+=1
                #end for
            #end for
        #end for

        mlab.contour3d(x,y,z,s)
        mlab.show()

        return
    #end def isosurface

    def mesh(self):
        return
    #end def mesh

    def etest(self):
        from enthought.mayavi import mlab
        from numpy import pi, sin, cos, exp, arange, array
        ni=10
        dr, dphi, dtheta = 1.0/ni, 2*pi/ni, pi/ni

        rlin = arange(0.0,1.0+dr,dr)
        plin = arange(0.0,2*pi+dphi,dphi)
        tlin = arange(0.0,pi+dtheta,dtheta)
        r,phi,theta = ndgrid(rlin,plin,tlin)

        a=1

        fr = .5*exp(-r/a)*(cos(2*pi*r/a)+1.0)
        fp = (1.0/6.0)*(cos(3.0*phi)+5.0)
        ft = (1.0/6.0)*(cos(10.0*theta)+5.0)

        f = fr*fp*ft

        x = r*sin(theta)*cos(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(theta)


        #mayavi
        #mlab.contour3d(x,y,z,f)
        #mlab.contour3d(r,phi,theta,f)
        i=7
        #mlab.mesh(x[i],y[i],z[i],scalars=f[i])
        mlab.mesh(f[i]*x[i],f[i]*y[i],f[i]*z[i],scalars=f[i])
        mlab.show()

        return
    #end def test


    def mtest(self):
        from enthought.mayavi import mlab
        # Create the data.
        from numpy import pi, sin, cos, mgrid, arange, array
        ni = 100.0
        dtheta, dphi = pi/ni, pi/ni

        #[theta,phi] = mgrid[0:pi+dtheta:dtheta,0:2*pi+dphi:dphi]

        #tlin = arange(0,pi+dtheta,dtheta)
        #plin = arange(0,2*pi+dphi,dphi)
        tlin = pi*array([0,.12,.2,.31,.43,.56,.63,.75,.87,.92,1])
        plin = 2*pi*array([0,.11,.22,.34,.42,.58,.66,.74,.85,.97,1])
        theta,phi = ndgrid(tlin,plin)

        fp = (1.0/6.0)*(cos(3.0*phi)+5.0)
        ft = (1.0/6.0)*(cos(10.0*theta)+5.0)

        r = fp*ft

        x = r*sin(theta)*cos(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(theta)

        # View it.
        s = mlab.mesh(x, y, z, scalars=r)
        mlab.show()
        return
    #end def


    def test(self):
        from enthought.mayavi import mlab
        from numpy import array,dot,arange,sin,ogrid,mgrid,zeros

        n=10
        n2=2*n
        s = '-'+str(n)+':'+str(n)+':'+str(n2)+'j'
        exec 'x, y, z = ogrid['+s+','+s+','+s+']'
        del s

        #x, y, z = ogrid[-10:10:20j, -10:10:20j, -10:10:20j]
        #x, y, z = mgrid[-10:11:1, -10:11:1, -10:11:1]

        s = sin(x*y*z)/(x*y*z)


        #xl = [-5.0,-4.2,-3.5,-2.1,-1.7,-0.4,0.7,1.8,2.6,3.7,4.3,5.0]
        #yl = [-5.0,-4.3,-3.6,-2.2,-1.8,-0.3,0.8,1.7,2.7,3.6,4.4,5.0]
        #zl = [-5.0,-4.4,-3.7,-2.3,-1.9,-0.4,0.9,1.6,2.8,3.5,4.5,5.0]
        dx = 2.0*n/(2.0*n-1.0)
        xl = arange(-n,n+dx,dx)
        yl = xl
        zl = xl

        x,y,z = ndgrid(xl,yl,zl)

        s2 = sin(x*y*z)/(x*y*z)

        #shear the grid
        nx,ny,nz = x.shape
        A = array([[1,1,-1],[1,-1,1],[-1,1,1]])
        #A = array([[3,2,1],[0,2,1],[0,0,1]])
        #A = array([[4,7,2],[8,4,3],[2,5,3]])
        #A = 1.0*array([[1,2,3],[4,5,6],[7,8,9]]).transpose()
        r = zeros((3,))
        np=0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    r[0] = x[i,j,k]
                    r[1] = y[i,j,k]
                    r[2] = z[i,j,k]
                    
                    #print np,r[0],r[1],r[2]
                    np+=1

                    r = dot(A,r)
                    x[i,j,k] = r[0]  
                    y[i,j,k] = r[1]  
                    z[i,j,k] = r[2]  
                #end for
            #end for
        #end for
        s2 = sin(x*y*z)/(x*y*z)

        mlab.contour3d(x,y,z,s2)
        mlab.show()

        out = QAobject()
        out.x=x
        out.y=y
        out.z=z
        out.s=s2
        out.A=A

        return out
    #end def


    def test_structured(self):

        import numpy as np
        from numpy import cos, sin, pi
        from enthought.tvtk.api import tvtk
        from enthought.mayavi import mlab

        def generate_annulus(r=None, theta=None, z=None):
            """ Generate points for structured grid for a cylindrical annular
                volume.  This method is useful for generating a unstructured
                cylindrical mesh for VTK (and perhaps other tools).

                Parameters
                ----------
                r : array : The radial values of the grid points.
                            It defaults to linspace(1.0, 2.0, 11).

                theta : array : The angular values of the x axis for the grid
                                points. It defaults to linspace(0,2*pi,11).

                z: array : The values along the z axis of the grid points.
                           It defaults to linspace(0,0,1.0, 11).

                Return
                ------
                points : array
                    Nx3 array of points that make up the volume of the annulus.
                    They are organized in planes starting with the first value
                    of z and with the inside "ring" of the plane as the first 
                    set of points.  The default point array will be 1331x3.
            """
            # Default values for the annular grid.
            if r is None: r = np.linspace(1.0, 2.0, 11)
            if theta is None: theta = np.linspace(0, 2*pi, 11)
            if z is None: z = np.linspace(0.0, 1.0, 11)

            # Find the x values and y values for each plane.
            x_plane = (cos(theta)*r[:,None]).ravel()
            y_plane = (sin(theta)*r[:,None]).ravel()

            # Allocate an array for all the points.  We'll have len(x_plane)
            # points on each plane, and we have a plane for each z value, so
            # we need len(x_plane)*len(z) points.
            points = np.empty([len(x_plane)*len(z),3])

            # Loop through the points for each plane and fill them with the
            # correct x,y,z values.
            start = 0
            for z_plane in z:
                end = start + len(x_plane)
                # slice out a plane of the output points and fill it
                # with the x,y, and z values for this plane.  The x,y
                # values are the same for every plane.  The z value
                # is set to the current z 
                plane_points = points[start:end]    
                plane_points[:,0] = x_plane
                plane_points[:,1] = y_plane    
                plane_points[:,2] = z_plane
                start = end

            return points

        # Make the data.
        dims = (51, 25, 25)
        # Note here that the 'x' axis corresponds to 'theta'
        theta = np.linspace(0, 2*np.pi, dims[0])
        # 'y' corresponds to varying 'r'
        r = np.linspace(1, 10, dims[1])
        z = np.linspace(0, 5, dims[2])
        pts = generate_annulus(r, theta, z)
        # Uncomment the following if you want to add some noise to the data.
        #pts += np.random.randn(dims[0]*dims[1]*dims[2], 3)*0.04
        sgrid = tvtk.StructuredGrid(dimensions=dims)
        sgrid.points = pts
        s = np.sqrt(pts[:,0]**2 + pts[:,1]**2 + pts[:,2]**2)
        sgrid.point_data.scalars = np.ravel(s.copy())
        sgrid.point_data.scalars.name = 'scalars'

        contour = mlab.pipeline.contour(sgrid)
        mlab.pipeline.surface(contour)


        return
    #end def test_structured



#end class EnergyDensityAnalyzer



