#
# try to find GTK (and glib) and GTKGLArea
#

# GTK_INCLUDE_DIR   - Directories to include to use GTK
# GTK_LIBRARIES     - Files to link against to use GTK
# GTK_FOUND         - If false, don't try to use GTK
# GTK_GL_FOUND      - If false, don't try to use GTK's GL features

# modified CMake/Modules/FindGTK.cmake
# - add gtkmm search
# don't even bother under WIN32
IF(UNIX)

  # common search paths for header files
  SET(GTK_INCLUDE_SEARCH_PATH
    $ENV{GTK_HOME}/include
    /usr/include
    /sw/include
    /usr/local/include
    /opt/gnome/include
    /usr/openwin/share/include
    /usr/openwin/include
  )

  # common search paths for libraries
  SET(GTK_LIBRARY_SEARCH_PATH
    $ENV{GTK_HOME}/lib
    /usr/lib
    /sw/lib
    /usr/openwin/lib
    /usr/X11R6/lib
    /opt/gnome/lib
  )

  FIND_PATH( GTK_pango_INCLUDE_PATH pango/pango.h
    /sw/include/pango-${PANGO_VERSION}
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  FIND_PATH( GTK_gtk_INCLUDE_PATH gtk/gtk.h
    /sw/include/gtk-${GTK_VERSION}
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  # Some Linux distributions (e.g. Red Hat) have glibconfig.h
  # and glib.h in different directories, so we need to look
  # for both.
  #  - Atanas Georgiev <atanas@cs.columbia.edu>

  FIND_PATH( GTK_glibconfig_INCLUDE_PATH glibconfig.h
    /sw/lib/glib-${GTK_VERSION}/include
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  FIND_PATH( GTK_glib_INCLUDE_PATH glib.h
    /sw/include/glib-${GTK_VERSION}
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  FIND_PATH( GTK_gtkgl_INCLUDE_PATH gtkgl/gtkglarea.h
    /sw/include/gtkgl-${GTK_VERSION}
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  FIND_PATH( GTK_gdkconfig_INCLUDE_PATH gdkconfig.h
    /sw/lib/gtk-${GTK_VERSION}/include
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  FIND_PATH( GTK_atk_INCLUDE_PATH atk/atk.h
    /sw/include/atk-1.0
    ${GTK_INCLUDE_SEARCH_PATH}
  )

  FIND_LIBRARY( GTK_gtkgl_LIBRARY gtkgl-${GTK_VERSION}
    ${GTK_LIBRARY_SEARCH_PATH}
  )

  #
  # The 12 suffix is thanks to the FreeBSD ports collection
  #

  FIND_LIBRARY( GTK_gtk_LIBRARY
    NAMES  gtk-x11-${GTK_VERSION} gtk 
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
  )

  FIND_LIBRARY( GTK_gdk_LIBRARY
    NAMES  gdk-x11-${GTK_VERSION} gdk
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
  )

  FIND_LIBRARY( GTK_gmodule_LIBRARY
    NAMES  gmodule-${GTK_VERSION} gmodule
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
  )

  FIND_LIBRARY( GTK_glib_LIBRARY
    NAMES  glib-${GTK_VERSION} glib
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
  )

  FIND_LIBRARY( GTK_Xi_LIBRARY 
    NAMES Xi 
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
    ) 

  FIND_LIBRARY( GTK_gthread_LIBRARY
    NAMES  gthread-${GTK_VERSION} gthread
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
  )

  FIND_LIBRARY( GTK_gobject_LIBRARY
    NAMES  gobject-${GTK_VERSION} gobject
    PATHS ${GTK_LIBRARY_SEARCH_PATH}
  )

  IF(GTK_pango_INCLUDE_PATH)
  IF(GTK_atk_INCLUDE_PATH)
  IF(GTK_gtk_INCLUDE_PATH)
  IF(GTK_glibconfig_INCLUDE_PATH)
  IF(GTK_gdkconfig_INCLUDE_PATH)
  IF(GTK_glib_INCLUDE_PATH)
  IF(GTK_gtk_LIBRARY)
  IF(GTK_glib_LIBRARY)

    # Assume that if gtk and glib were found, the other
    # supporting libraries have also been found.

    SET( GTK_FOUND "YES" )
    SET( GTK_INCLUDE_DIR  ${GTK_pango_INCLUDE_PATH}
                           ${GTK_atk_INCLUDE_PATH}
                           ${GTK_gtk_INCLUDE_PATH}
                           ${GTK_glibconfig_INCLUDE_PATH}
                           ${GTK_gdkconfig_INCLUDE_PATH}
                           ${GTK_glib_INCLUDE_PATH} )
    SET( GTK_LIBRARIES  ${GTK_gtk_LIBRARY}
                        ${GTK_gdk_LIBRARY}
                        ${GTK_glib_LIBRARY} )

    IF(GTK_gmodule_LIBRARY)
      SET(GTK_LIBRARIES ${GTK_LIBRARIES} ${GTK_gmodule_LIBRARY})
    ENDIF(GTK_gmodule_LIBRARY)
    IF(GTK_gthread_LIBRARY)
      SET(GTK_LIBRARIES ${GTK_LIBRARIES} ${GTK_gthread_LIBRARY})
    ENDIF(GTK_gthread_LIBRARY)
    IF(GTK_Xi_LIBRARY)
      SET(GTK_LIBRARIES ${GTK_LIBRARIES} ${GTK_Xi_LIBRARY})
    ENDIF(GTK_Xi_LIBRARY)
    IF(GTK_gobject_LIBRARY)
      SET(GTK_LIBRARIES ${GTK_LIBRARIES} ${GTK_gobject_LIBRARY})
    ENDIF(GTK_gobject_LIBRARY)

    IF(GTK_gtkgl_INCLUDE_PATH)
      IF(GTK_gtkgl_LIBRARY)
        SET( GTK_GL_FOUND "YES" )
        SET( GTK_INCLUDE_DIR  ${GTK_INCLUDE_DIR}
                               ${GTK_gtkgl_INCLUDE_PATH} )
        SET( GTK_LIBRARIES  ${GTK_gtkgl_LIBRARY} ${GTK_LIBRARIES} )
        MARK_AS_ADVANCED(
          GTK_gtkgl_LIBRARY
          GTK_gtkgl_INCLUDE_PATH
          )
      ENDIF(GTK_gtkgl_LIBRARY)
    ENDIF(GTK_gtkgl_INCLUDE_PATH)

  ENDIF(GTK_glib_LIBRARY)
  ENDIF(GTK_gtk_LIBRARY)
  ENDIF(GTK_glib_INCLUDE_PATH) 
  ENDIF(GTK_gdkconfig_INCLUDE_PATH)
  ENDIF(GTK_glibconfig_INCLUDE_PATH)
  ENDIF(GTK_gtk_INCLUDE_PATH)
  ENDIF(GTK_atk_INCLUDE_PATH)
  ENDIF(GTK_pango_INCLUDE_PATH)


  MARK_AS_ADVANCED(
    GTK_gdk_LIBRARY
    GTK_glib_INCLUDE_PATH
    GTK_glib_LIBRARY
    GTK_glibconfig_INCLUDE_PATH
    GTK_gmodule_LIBRARY
    GTK_gthread_LIBRARY
    GTK_Xi_LIBRARY
    GTK_gtk_INCLUDE_PATH
    GTK_gtk_LIBRARY
    GTK_gtkgl_INCLUDE_PATH
    GTK_gtkgl_LIBRARY
    GTK_pango_INCLUDE_PATH
    GTK_gdkconfig_INCLUDE_PATH
    GTK_atk_INCLUDE_PATH
  )

  # gtkmm
  FIND_PATH( GTK_gtkmm_INCLUDE_PATH gtkmm.h
    /usr/include
    /sw/include/gtkmm-${GTKMM_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_gtkmmconfig_INCLUDE_PATH gtkmmconfig.h
    /usr/include
    /sw/lib/gtkmm-${GTKMM_VERSION}/include
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_glibmm_INCLUDE_PATH glibmm.h
    /usr/include
    /sw/include/glibmm-${GTKMM_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_glibmmconfig_INCLUDE_PATH glibmmconfig.h
    /usr/include
    /sw/lib/glibmm-${GTKMM_VERSION}/include
    /sw/include/glibmm-${GTKMM_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_atkmm_INCLUDE_PATH atkmm.h
    /usr/include
    /sw/include/atkmm-1.6
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_gdkmm_INCLUDE_PATH gdkmm.h
    /usr/include
    /sw/include/gdkmm-${GTKMM_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_gdkmmconfig_INCLUDE_PATH gdkmmconfig.h
    /usr/include
    /sw/lib/gdkmm-${GTKMM_VERSION}/include
    /sw/include/glibmm-${GTKMM_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_pangomm_INCLUDE_PATH pangomm.h
    /usr/include
    /sw/include/pangomm-${PANGOMM_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_sigc_INCLUDE_PATH sigc++/sigc++.h
    /usr/include
    /sw/include/sigc++-${SIGCPP_VERSION}
    /sw/include
    /opt/gnome/include
  )

  FIND_PATH( GTK_sigcconfig_INCLUDE_PATH sigc++config.h
    /usr/include
    /sw/lib/sigc++-${SIGCPP_VERSION}/include
    /sw/include
    /opt/gnome/include
  )

  FIND_LIBRARY( GTK_gtkmm_LIBRARY
    NAMES  gtkmm-${GTKMM_VERSION} gtkmm
    PATHS  /usr/lib
           /sw/lib
           /usr/openwin/lib
           /usr/X11R6/lib
           /opt/gnome/lib
  )

  FIND_LIBRARY( GTK_atkmm_LIBRARY
    NAMES  atkmm-1.6 atkmm
    PATHS  /usr/lib
           /sw/lib
           /usr/openwin/lib
           /usr/X11R6/lib
           /opt/gnome/lib
  )

  FIND_LIBRARY( GTK_gdkmm_LIBRARY
    NAMES  gdkmm-${GTKMM_VERSION} gdkmm
    PATHS  /usr/lib
           /sw/lib
           /usr/openwin/lib
           /usr/X11R6/lib
           /opt/gnome/lib
  )

  FIND_LIBRARY( GTK_glibmm_LIBRARY
    NAMES  glibmm-${GTKMM_VERSION} glibmm
    PATHS  /usr/lib
           /sw/lib
           /usr/openwin/lib
           /usr/X11R6/lib
           /opt/gnome/lib
  )

  FIND_LIBRARY( GTK_pangomm_LIBRARY
    NAMES  pangomm-${PANGOMM_VERSION} pangomm
    PATHS  /usr/lib
           /sw/lib
           /usr/openwin/lib
           /usr/X11R6/lib
           /opt/gnome/lib
  )

  FIND_LIBRARY( GTK_sigc_LIBRARY
    NAMES  sigc-${SIGCPP_VERSION} sigc
    PATHS  /usr/lib
           /sw/lib
           /usr/openwin/lib
           /usr/X11R6/lib
           /opt/gnome/lib
  )

 IF(GTK_gtkmm_INCLUDE_PATH)
 IF(GTK_glibmmconfig_INCLUDE_PATH)
 IF(GTK_sigc_INCLUDE_PATH)
 IF(GTK_gtkmm_INCLUDE_PATH)
 IF(GTK_atkmm_INCLUDE_PATH)
 IF(GTK_gdkmm_INCLUDE_PATH)
 IF(GTK_glibmm_INCLUDE_PATH)
 IF(GTK_pangomm_INCLUDE_PATH)
   SET( GTKMM_FOUND "YES" )
   SET( GTKMM_INCLUDE_DIR  ${GTK_INCLUDE_DIR} 
                           ${GTK_gtkmm_INCLUDE_PATH} 
                           ${GTK_gtkmmconfig_INCLUDE_PATH} 
                           ${GTK_sigc_INCLUDE_PATH}
                           ${GTK_sigcconfig_INCLUDE_PATH}
                           ${GTK_glibmm_INCLUDE_PATH}
                           ${GTK_glibmmconfig_INCLUDE_PATH}
                           ${GTK_gdkmm_INCLUDE_PATH}
                           ${GTK_gdkmmconfig_INCLUDE_PATH}
                           ${GTK_pangomm_INCLUDE_PATH}
                           ${GTK_atkmm_INCLUDE_PATH}
   )
   #link everything
   IF(GTK_gtkmm_LIBRARY)
   IF(GTK_atkmm_LIBRARY)
   IF(GTK_gdkmm_LIBRARY)
   IF(GTK_glibmm_LIBRARY)
   IF(GTK_pangomm_LIBRARY)
   IF(GTK_sigc_LIBRARY)
     SET( GTKMM_LIBRARIES  ${GTK_gtkmm_LIBRARY} 
                           ${GTK_atkmm_LIBRARY} 
                           ${GTK_gdkmm_LIBRARY} 
                           ${GTK_glibmm_LIBRARY} 
                           ${GTK_pangomm_LIBRARY} 
                           ${GTK_sigc_LIBRARY} 
                           ${GTK_LIBRARIES} )
   ENDIF(GTK_sigc_LIBRARY)
   ENDIF(GTK_pangomm_LIBRARY)
   ENDIF(GTK_glibmm_LIBRARY)
   ENDIF(GTK_gdkmm_LIBRARY)
   ENDIF(GTK_atkmm_LIBRARY)
   ENDIF(GTK_gtkmm_LIBRARY)

 ENDIF(GTK_pangomm_INCLUDE_PATH)
 ENDIF(GTK_glibmm_INCLUDE_PATH)
 ENDIF(GTK_gdkmm_INCLUDE_PATH)
 ENDIF(GTK_atkmm_INCLUDE_PATH)
 ENDIF(GTK_gtkmm_INCLUDE_PATH)
 ENDIF(GTK_sigc_INCLUDE_PATH)
 ENDIF(GTK_glibmmconfig_INCLUDE_PATH)
 ENDIF(GTK_gtkmm_INCLUDE_PATH)

 MARK_AS_ADVANCED(
   GTK_INCLUDE_DIR  
   GTK_gtkmm_INCLUDE_PATH  
   GTK_gtkmmconfig_INCLUDE_PATH  
   GTK_sigc_INCLUDE_PATH 
   GTK_sigcconfig_INCLUDE_PATH 
   GTK_glibmm_INCLUDE_PATH 
   GTK_glibmmconfig_INCLUDE_PATH 
   GTK_gdkmm_INCLUDE_PATH 
   GTK_gdkmmconfig_INCLUDE_PATH 
   GTK_pangomm_INCLUDE_PATH 
   GTK_atkmm_INCLUDE_PATH 
   GTK_sigc_LIBRARY 
   GTK_pangomm_LIBRARY 
   GTK_glibmm_LIBRARY 
   GTK_gdkmm_LIBRARY 
   GTK_atkmm_LIBRARY 
   GTK_gtkmm_LIBRARY 
  )
ENDIF(UNIX)

