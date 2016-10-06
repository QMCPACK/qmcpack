#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
#//
#// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign 
#//////////////////////////////////////////////////////////////////////////////////////


import pygtk
import numpy
pygtk.require('2.0')
from trackball import *


haveGTK_GL = True
try:
    import gtk.gtkgl
    from OpenGL.GL import *
    from OpenGL.GLU import *
except ImportError:
    haveGTK_GL = False


def ViewerOkay():
    return haveGTK_GL


#(file, pathname, desc) = imp.find_module("gtk.gtkgl.apputils")
#if (file != None):
#    print "Found gtk.gtkgl.apputils.  Description = " + desc
#    imp.load_module("gtk.gtkgl.apputils", file, pathname, desc)
    

class StructureViewer(gtk.gtkgl.DrawingArea):
    def __init__(self):
        try:
            # try double-buffered
            glconfig = gtk.gdkgl.Config(mode=(gtk.gdkgl.MODE_RGB    |
                                              gtk.gdkgl.MODE_DOUBLE |
                                              gtk.gdkgl.MODE_DEPTH))
        except gtk.gdkgl.NoMatches:
            # try single-buffered
            glconfig = gtk.gdkgl.Config(mode=(gtk.gdkgl.MODE_RGB    |
                                              gtk.gdkgl.MODE_DEPTH))

        gtk.gtkgl.DrawingArea.__init__(self, glconfig)
        
        self.connect_after("realize"      , self.init        )
        self.connect("configure_event"    , self.reshape     )
        self.connect("expose_event"       , self.display     )
        self.connect("map_event"          , self.map         )
        self.connect("button_press_event" , self.button_press)
        self.connect("motion_notify_event", self.button_motion)
        self.atom_pos = 8.0*numpy.array([[0.00, 0.00, 0.00],\
                                         [0.25, 0.25, 0.25]])
        self.lattice = 8.0*numpy.array([[0.50, 0.50, 0.00],\
                                        [0.50, 0.00, 0.50],\
                                        [0.00, 0.50, 0.00]])
        self.set_events(gtk.gdk.BUTTON1_MOTION_MASK    |
                        gtk.gdk.BUTTON2_MOTION_MASK    |
                        gtk.gdk.BUTTON3_MOTION_MASK    |
                        gtk.gdk.BUTTON_PRESS_MASK      |
                        gtk.gdk.BUTTON_RELEASE_MASK    |
                        gtk.gdk.VISIBILITY_NOTIFY_MASK |
                        gtk.gdk.SCROLL_MASK)
        self.Scale = 1.0
        self.Distance = 3.0
        self.tb = Trackball()
    

    def button_press(self, glDrawArea, event):
        if (event.button == 1):
            self.StartX = event.x
            self.StartY = event.y
            self.Button1Pressed=True
            self.Button2Pressed=False
            self.Button3Pressed=False
        elif (event.button == 2):
            self.StartX = event.x
            self.StartY = event.y
            self.Button1Pressed=False
            self.Button2Pressed=True
            self.Button3Pressed=False
            self.OldScale = self.Scale
        elif (event.button == 3):
            self.StartX = event.x
            self.StartY = event.y
            self.Button1Pressed=False
            self.Button2Pressed=False
            self.Button3Pressed=True
        print 'button pressed at (%3d,%3d)' % (event.x, event.y)

    def button_motion(self, glDrawArea, event):
        if (self.Button3Pressed):
            print 'translate'
        elif (self.Button1Pressed):
            w = self.get_allocation()[2]
            h = self.get_allocation()[3]
            x = event.x
            y = event.y
            self.tb.update(self.StartX, self.StartY, x, y, w, h)
            self.display(self, None)
            #print self.tb.matrix
         
    def set_lattice(self, lattice):
        self.lattice = lattice
        self.BoxList = glGenLists(1)
        glNewList (self.BoxList, GL_COMPILE)
        a = []
        ma = []
        a.append (numpy.array([lattice[0,0], lattice[0,1], lattice[0,2]]))
        a.append (numpy.array([lattice[1,0], lattice[1,1], lattice[1,2]]))
        a.append (numpy.array([lattice[2,0], lattice[2,1], lattice[2,2]]))
        ma.append(-1.0*a[0])
        ma.append(-1.0*a[1])
        ma.append(-1.0*a[2])
        r = []
        r.append(0.5*(ma[0] + ma[1] + ma[2]));
        r.append(0.5*(ma[0] + ma[1] +  a[2]));
        r.append(0.5*(ma[0] +  a[1] +  a[2]));
        r.append(0.5*(ma[0] +  a[1] + ma[2]));
        r.append(0.5*( a[0] + ma[1] + ma[2]));
        r.append(0.5*( a[0] + ma[1] +  a[2]));
        r.append(0.5*( a[0] +  a[1] +  a[2]));
        r.append(0.5*( a[0] +  a[1] + ma[2]));
        glColor3d (0.3, 0.3, 0.3)
        c = (1.0, 1.0, 1.0, 1.0)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, c)
        glLineWidth (2.0);


        p01 = numpy.cross (self.lattice[0], self.lattice[1])
        p12 = numpy.cross (self.lattice[1], self.lattice[2]);
        p20 = numpy.cross (self.lattice[2], self.lattice[0]);
        p01 = 1.0/numpy.sqrt(numpy.dot(p01, p01)) * p01;
        p12 = 1.0/numpy.sqrt(numpy.dot(p12, p12)) * p12;
        p20 = 1.0/numpy.sqrt(numpy.dot(p20, p20)) * p20;
        d01 = abs(0.5001*numpy.dot(self.lattice[2], p01));
        d12 = abs(0.5001*numpy.dot(self.lattice[0], p12));
        d20 = abs(0.5001*numpy.dot(self.lattice[1], p20));
        
        eqn0 = ( p01[0], p01[1], p01[2], d01);
        eqn1 = (-p01[0],-p01[1],-p01[2], d01);
        eqn2 = ( p12[0], p12[1], p12[2], d12);
        eqn3 = (-p12[0],-p12[1],-p12[2], d12);
        eqn4 = ( p20[0], p20[1], p20[2], d20);
        eqn5 = (-p20[0],-p20[1],-p20[2], d20);
        glClipPlane(GL_CLIP_PLANE0, eqn0);
        glClipPlane(GL_CLIP_PLANE1, eqn1);
        glClipPlane(GL_CLIP_PLANE2, eqn2);
        glClipPlane(GL_CLIP_PLANE3, eqn3);
        glClipPlane(GL_CLIP_PLANE4, eqn4);
        glClipPlane(GL_CLIP_PLANE5, eqn5);
#        glEnable(GL_CLIP_PLANE0);
#        glEnable(GL_CLIP_PLANE1);
#        glEnable(GL_CLIP_PLANE2);
#        glEnable(GL_CLIP_PLANE3);
#        glEnable(GL_CLIP_PLANE4);
#        glEnable(GL_CLIP_PLANE5);
        


        glBegin(GL_LINES);
        glNormal3d(1.0, 1.0, 1.0)
        glVertex3d(r[0][0],r[0][1],r[0][2]); 
        glVertex3d(r[1][0],r[1][1],r[1][2]);
        glVertex3d(r[1][0],r[1][1],r[1][2]); 
        glVertex3d(r[2][0],r[2][1],r[2][2]);
        glVertex3d(r[2][0],r[2][1],r[2][2]); 
        glVertex3d(r[3][0],r[3][1],r[3][2]);
        glVertex3d(r[3][0],r[3][1],r[3][2]); 
        glVertex3d(r[0][0],r[0][1],r[0][2]);
        glVertex3d(r[4][0],r[4][1],r[4][2]); 
        glVertex3d(r[5][0],r[5][1],r[5][2]);
        glVertex3d(r[5][0],r[5][1],r[5][2]); 
        glVertex3d(r[6][0],r[6][1],r[6][2]);
        glVertex3d(r[6][0],r[6][1],r[6][2]); 
        glVertex3d(r[7][0],r[7][1],r[7][2]);
        glVertex3d(r[7][0],r[7][1],r[7][2]); 
        glVertex3d(r[4][0],r[4][1],r[4][2]);
        glVertex3d(r[0][0],r[0][1],r[0][2]); 
        glVertex3d(r[4][0],r[4][1],r[4][2]);
        glVertex3d(r[1][0],r[1][1],r[1][2]); 
        glVertex3d(r[5][0],r[5][1],r[5][2]);
        glVertex3d(r[2][0],r[2][1],r[2][2]); 
        glVertex3d(r[6][0],r[6][1],r[6][2]);
        glVertex3d(r[3][0],r[3][1],r[3][2]); 
        glVertex3d(r[7][0],r[7][1],r[7][2]);
        glEnd()
        glEndList()

    def init (self, glDrawingArea):
        glcontext  = self.get_gl_context()
        gldrawable = self.get_gl_drawable()
        
        self.DisplayLists = []
        self.SphereList = glGenLists(1)
        glNewList(self.SphereList, GL_COMPILE)
        gtk.gdkgl.draw_sphere(True, 1.0, 30, 30)
        glEndList()
        self.set_lattice(self.lattice)
        self.redraw()


        glShadeModel(GL_SMOOTH);
        glEnable (GL_LIGHTING);
        glEnable (GL_LINE_SMOOTH);
        glEnable (GL_POLYGON_SMOOTH);
        glDisable (GL_POLYGON_SMOOTH);
#        glEnable (GL_MULTISAMPLE);
        glEnable (GL_COLOR_MATERIAL);
        glEnable(GL_LIGHT0)

        diffuse  = (1.0, 1.0, 1.0, 1.0)
        ambient  = (0.001, 0.001, 0.001, 1.0)
        specular = (1.0, 1.0, 1.0, 1.0)
        position = (1.0, 1.0, 1.0, 0.0)
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse)
        glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient)
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular)
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, specular)
        glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 92)
        glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE)

        (width, height) = self.window.get_size()
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(40.0, 1.0, 1.0, 10.0)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
#        gluLookAt(0.0, 0.0, 3.0,\
#                  0.0, 0.0, 0.0,\
#                  0.0, 1.0, 0.0)
        glEnable(GL_AUTO_NORMAL)
        glEnable(GL_NORMALIZE)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_COLOR_MATERIAL)
        self.reshape (None, None)

    def reshape(self, glDrawArea, event):
        # get GLContext and GLDrawable
        glcontext = self.get_gl_context()
        gldrawable = self.get_gl_drawable()
        
        # GL calls
        if not gldrawable.gl_begin(glcontext): return
        
        x, y, width, height = self.get_allocation()
        
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if width > height:
            w = float(width) / float(height)
            glFrustum(-w, w, -1.0, 1.0, 5.0, 60.0)
        else:
            h = float(height) / float(width)
            glFrustum(-1.0, 1.0, -h, h, 5.0, 60.0)
            
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -40.0)
            
        gldrawable.gl_end()
    
        return True

    def redraw(self):
#        glMatrixMode(GL_PROJECTION)
#        glLoadIdentity();
#        width  = self.get_allocation()[2]
#        height = self.get_allocation()[3]
#        ratio = width/height
#        gluPerspective(40.0, ratio, 1.0, 8.0*self.Distance/Scale)
#        glMatrixMode(GL_MODELVIEW);
#        glLoadIdentity();
#        gluLookAt(0.0, 0.0, self.Distance/Scale,\
#                  0.0, 0.0,            0.0,\
#                  0.0, 1.0, 0.0);
#        glTranslatef(0.0, 0.0, -self.Distance/Scale);
#        glScaled(Scale, Scale, Scale);

        for l in self.DisplayLists:
            glDeleteLists(l,1)
        for r in self.atom_pos:
            rtran = r - 0.5*(self.lattice[0] + self.lattice[1] + self.lattice[2])
            list = glGenLists(1)
            self.DisplayLists.append(list)
            glNewList(list, GL_COMPILE)
            glPushMatrix();
            glTranslated (rtran[0], rtran[1], rtran[2])
            glCallList(self.SphereList)
            glPopMatrix();
            glEndList()

    def display(self, glDrawArea, event):
        # get GLContext and GLDrawable
        glcontext  = self.get_gl_context()
        gldrawable = glDrawArea.get_gl_drawable()
        
        # GL calls
        if not gldrawable.gl_begin(glcontext): return
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glPushMatrix()
        glMultMatrixd(self.tb.matrix);
        # Do stuff here
#       glMatrixMode(GL_MODELVIEW)
#       glLoadIdentity()
        glTranslatef(0.0, 0.0, -5.0)
        
        #gtk.gdkgl.draw_sphere(True, 1.0, 30, 30)
        #glCallList(self.SphereList)
        glCallList(self.BoxList)
        glColor ([1.0, 0.0, 0.0, 1.0])
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        for list in self.DisplayLists:
            glCallList(list)

        glPopMatrix()
        if gldrawable.is_double_buffered():
            gldrawable.swap_buffers()
        else:
            glFlush()
        gldrawable.gl_end()
        # Invalidate whole window.
        # self.window.invalidate_rect(self.allocation, False)
        # Update window synchronously (fast).
        # self.window.process_updates(False)

    def map(self, event, dummy):
        print "map_event"

#    def set_lattice (self, lattice):
#        print "set_lattice"

    def set_atoms (self, atomList):
        print "set_atoms"
        
