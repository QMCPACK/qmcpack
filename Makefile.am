SUBDIRS = src

# Configuration for config directory
AUX_DIST          = \
                  $(ac_aux_dir)/config.guess \
                  $(ac_aux_dir)/config.sub \
                  $(ac_aux_dir)/install.sh \
                  $(ac_aux_dir)/ltconfig \
                  $(ac_aux_dir)/ltmain.sh \
                  $(ac_aux_dir)/mdate-sh \
                  $(ac_aux_dir)/missing \
                  $(ac_aux_dir)/mkinstalldirs
AUX_DIST_EXTRA    = \
                  $(ac_aux_dir)/ac_cxx_have_complex_math1.m4 \
                  $(ac_aux_dir)/ac_cxx_have_complex_math2.m4 \
                  $(ac_aux_dir)/ac_cxx_namespaces.m4 \
                  $(ac_aux_dir)/ac_path_generic.m4 \
                  $(ac_aux_dir)/acx_blas.m4 \
                  $(ac_aux_dir)/acx_lapack.m4 \
                  $(ac_aux_dir)/acx_mpi.m4 \
                  $(ac_aux_dir)/acx_pthread.m4 \
                  $(ac_aux_dir)/acx_hdf5.m4 \
                  $(ac_aux_dir)/acx_restrict.m4
EXTRA_DIST        = bootstrap

#MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure stamp.h.in $(AUX_DIST)
MAINTAINERCLEANFILES = Makefile.in configure stamp.h.in $(AUX_DIST)

dist-hook:
	(cd $(distdir) && mkdir $(ac_aux_dir))
	for file in $(AUX_DIST) $(AUX_DIST_EXTRA); do \
		cp $$file $(distdir)/$$file; \
	done
