#! /usr/bin/env python

import sys

def error(mark):
    print 'label_error;mark_'+mark
    sys.exit(1)
#end def error


# universally broken tests (irrespective of build)
universally_broken = set('''
    developer-heg_54_J2rpa

    short-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi
    short-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell
    short-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell

    short-diamondC_1x1x1_pp-dmc-estimator-spindensity
    short-diamondC_1x1x1_pp-dmc-estimator-density

    short-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm-4-4

    short-diamondC-afqmc_hyb_nn2
    short-diamondC-afqmc_incmf_nn2
    short-diamondC-afqmc_nn2
    '''.split())


label_sets = dict(
    unstable = universally_broken,
    )



# extract test name and build flags from args
try:
    full_test,qmc_cuda,enable_soa,qmc_complex,qmc_mixed = sys.argv[1:]
    qmc_cuda    = qmc_cuda=='1'
    enable_soa  = enable_soa=='1'
    qmc_complex = qmc_complex=='1'
    qmc_mixed   = qmc_mixed=='1'
    cpu   = not qmc_cuda
    gpu   = qmc_cuda
    aos   = not enable_soa
    soa   = enable_soa
    real  = not qmc_complex
    comp  = qmc_complex
    full  = not qmc_mixed
    mixed = qmc_mixed
except:
    error('extract')
#end try


# final set of labels for test
labels = []


# get a shortened name for the test without #mpi and #omp parts
test = full_test
tokens = test.replace('-','_').split('_')
if len(tokens)>1:
    mpi = tokens[-2]
    omp = tokens[-1]
    if mpi.isdigit() and omp.isdigit():
        test = full_test.rsplit('-',2)[0]
    #end if
#end if


# get labels from known sets
for label,label_set in label_sets.iteritems():
    if test in label_set or full_test in label_set:
        labels.append(label)
    #end if
#end for


# make a ctest list of the labels
ctest_labels = ''
for label in labels:
    ctest_labels += label+';'
#end for
ctest_labels = ctest_labels.rstrip(';')


# print out the label list for ctest to read
print ctest_labels
