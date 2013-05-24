
from abilities import AllAbilities

import numpy as np

class Plotter(AllAbilities):
    def __init__(self):
        self.initialized = False
        return
    #end def __init__

    def ensure_init(self):
        if not self.initialized:
            from enthought.mayavi import mlab
            from enthought.tvtk.api import tvtk
            self.mlab = mlab
            self.tvtk = tvtk

            self.show   = mlab.show
            self.plot3d = mlab.plot3d
            self.mesh   = mlab.mesh

            self.initialized = True
        #end if
    #end def ensure_init

    def isosurface(self,points,scalars,contours,dimensions,name='val'):
        self.ensure_init()
        mlab = self.mlab
        tvtk = self.tvtk
        sg=tvtk.StructuredGrid(dimensions=dimensions,points=points)
        sg.point_data.scalars = scalars
        sg.point_data.scalars.name = name
        d = mlab.pipeline.add_dataset(sg)
        iso = mlab.pipeline.iso_surface(d)
        if isinstance(contours,int):
            iso.contour.number_of_contours = contours
        elif isinstance(contours,list):
            iso.contour.auto_contours = False
            iso.contour.contours = contours
        else:
            self.error('isosurface contours must be an int or list\n  a '+str(type(contours))+' was provided instead')
        #end if
        return
    #end def isosurface

    def surface_slice(self,x,y,z,scalars,options=None):
        scale = 1.0
        opacity= 1.0
        if options!=None:
            if 'norm_height' in options:
                scale = options.norm_height/abs(scalars.max())
            if 'scale' in options:
                scale = options.scale
            if 'opacity' in options:
                opacity = options.opacity
        #end if
        self.ensure_init()
        from extended_numpy import surface_normals
        self.mesh(x,y,z,opacity=.2)
        surfnorm = scale*surface_normals(x,y,z)
        xs=x.copy()
        ys=y.copy()
        zs=z.copy()
        xs[...] = x[...] + surfnorm[...,0]*scalars[...]
        ys[...] = y[...] + surfnorm[...,1]*scalars[...]
        zs[...] = z[...] + surfnorm[...,2]*scalars[...]
        self.mesh(xs,ys,zs,scalars=scalars,opacity=opacity)
        return
    #end def surface_slice
#end class Plotter
