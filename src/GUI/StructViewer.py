import pygtk
pygtk.require('2.0')


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
        
        self.connect_after("realize"  , self.init   )
        self.connect("configure_event", self.reshape)
        self.connect("expose_event"   , self.display)
        self.connect("map_event"      , self.map    )

    def init (self, glDrawingArea):
        glcontext  = self.get_gl_context()
        gldrawable = self.get_gl_drawable()

        glShadeModel(GL_SMOOTH);
        glEnable (GL_LIGHTING);
        glEnable (GL_LINE_SMOOTH);
        glEnable (GL_POLYGON_SMOOTH);
        glDisable (GL_POLYGON_SMOOTH);
#        glEnable (GL_MULTISAMPLE);
        glEnable (GL_COLOR_MATERIAL);
        glEnable(GL_LIGHT0)

        diffuse  = (1.0, 1.0, 1.0, 1.0)
        ambient  = (0.2, 0.2, 0.2, 1.0)
        specular = (1.0, 1.0, 1.0, 1.0)
        position = (1.0, 1.0, 2.0, 0.0)
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse)
        glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient)
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular)
        glLightfv(GL_LIGHT0, GL_POSITION, diffuse)

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

    def display(self, glDrawArea, event):
        # get GLContext and GLDrawable
        glcontext  = self.get_gl_context()
        gldrawable = glDrawArea.get_gl_drawable()
        
        # GL calls
        if not gldrawable.gl_begin(glcontext): return
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glPushMatrix()
        # Do stuff here
#       glMatrixMode(GL_MODELVIEW)
#       glLoadIdentity()
        glTranslatef(0.0, 0.0, -5.0)
        
        glColor ([1.0, 0.0, 0.0, 1.0])
        gtk.gdkgl.draw_sphere(True, 1.0, 30, 30)

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

    def set_lattice (self, lattice):
        print "set_lattice"

    def set_atoms (self, atomList):
        print "set_atoms"
        
