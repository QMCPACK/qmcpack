QMC_APP_SRCS =\
  ParticleSetPool.cpp\
  WaveFunctionPool.cpp\
  HamiltonianPool.cpp\
  InitMolecularSystem.cpp\
  QMCAppBase.cpp\
  QMCDriverFactory.cpp\
  QMCMain.cpp

lib_LIBRARIES = libqmc.a
libqmc_a_SOURCES = \
$(QMC_APP_SRCS)

bin_PROGRAMS = qmcapp
qmcapp_SOURCES = \
qmcapp.cpp

if NEEDTOBUILDOOMPI
OOMPI_INC = -I$(top_srcdir)/src/OOMPI
OOMPI_LIB = ${top_builddir}/src/OOMPI/liboompi.a
endif

INCLUDES = -I$(top_srcdir)/src $(OOMPI_INC)

QMC_LIB=$(top_builddir)/src/QMCApp/libqmc.a 

qmcapp_LDADD = $(QMC_LIB) @QMCBASE_LIBS@  $(OOMPI_LIB)
qmcapp_LDFLAGS = @LDFLAGS@ 
qmcapp_DEPENDENCIES = ${top_builddir}/src/libqmcbase.a ${top_builddir}/src/libqmcutil.a $(OOMPI_LIB) $(QMC_LIB)
