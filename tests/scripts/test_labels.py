#! /usr/bin/env python

import sys

def error(mark):
    sys.stdout.write('label_error;mark_'+mark)
    sys.exit(1)
#end def error


# label descriptions
#
#   main classification
#     stable   - test is expected to pass
#     unstable - test may fail (definitely or with some probability)
#     (stable/unstable are exclusive)
#
#   feature classification of unstable tests
#     broken      - test fails because tested feature is broken (bug)
#     unsupported - test fails because tested feature is not supported and never will be
#     (broken/unsupported are exclusive but not comprehensive)
#   
#   failure type classification of unstable tests
#     hard_fail              - ungraceful crash (segfault, etc) (bug)
#     abort                  - controlled abort (bug)
#     definite_stat_fail     - statistical failure by > 5 sigma (bug)
#     intermittent_stat_fail - statistical failure by < 5 sigma
#
#     code_fail              - hard fail or abort
#     stat_fail              - definite or intermittent statistical failure
#
#   other classifications
#     bug - broken or code_fail or definite_stat_fail
#


def create_label_sets():

    # primary sets, all others derived from these
    broken        = set()
    unsupported   = set()
    hard_fail     = set()
    abort         = set()
    def_stat_fail = set()
    int_stat_fail = set()

    # universally broken tests (irrespective of build)
    broken |= set([
        # https://github.com/QMCPACK/qmcpack/issues/848
        'developer-heg_54_J2rpa',

        # https://github.com/QMCPACK/qmcpack/issues/934
        'short-diamondC_1x1x1_pp-dmc-estimator-density',
        'short-diamondC_1x1x1_pp-dmc-estimator-spindensity',
        
        # https://github.com/QMCPACK/qmcpack/issues/995
        # https://github.com/QMCPACK/qmcpack/issues/1052
        'short-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm-4-4',

        # https://github.com/QMCPACK/qmcpack/issues/982
        'short-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi',
        'short-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell',
        'short-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell',

        # https://github.com/QMCPACK/qmcpack/issues/1040
        'short-diamondC-afqmc_hyb_nn2',
        'short-diamondC-afqmc_incmf_nn2',
        'short-diamondC-afqmc_nn2',
        ])

    # aos specific (but general otherwise)
    if aos:
        unsupported |= set()
    #end if

    # soa specific (but general otherwise)
    if soa:
        broken |= set()
    #end if

    # gpu specific (but general otherwise)
    if gpu:
        broken      |= set()
        unsupported |= set()
    #end if

    # real specific (but general otherwise)
    if real:
        unsupported |= set()
    #end if

    # complex specific (but general otherwise)
    if comp:
        unsupported |= set()
    #end if
    
    # mixed precision specific (but general otherwise)
    if mixed:
        broken      |= set()
        unsupported |= set()
    #end if

    # finer details based on build
    # AoS issues
    if aos and cpu and real and full:
        None
    elif aos and cpu and real and mixed:
        None
    elif aos and cpu and comp and full:
        None
    elif aos and cpu and comp and mixed:
        None
    elif aos and gpu and real:
        None
    elif aos and gpu and comp:
        None
    # SoA issues
    elif soa and cpu and real and full:
        None
    elif soa and cpu and real and mixed:
        None
    elif soa and cpu and comp and full:
        None
    elif soa and cpu and comp and mixed:
        None
    elif soa and gpu and real:
        None
    elif soa and gpu and comp:
        None
    #end if

    code_fail = hard_fail | abort
    stat_fail = def_stat_fail | int_stat_fail

    bug = broken | code_fail | def_stat_fail

    unstable = broken | unsupported | code_fail | stat_fail

    label_sets = dict(
        unstable               = unstable,
        broken                 = broken,
        unsupported            = unsupported,
        hard_fail              = hard_fail,
        abort                  = abort,
        definite_stat_fail     = def_stat_fail,
        intermittent_stat_fail = int_stat_fail,
        code_fail              = code_fail,
        stat_fail              = stat_fail,
        bug                    = bug,
        )

    return label_sets
#end def create_label_sets


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


# get a shortened name for the test without #mpi and #omp parts
try:
    test = full_test
    tokens = test.replace('-','_').split('_')
    if len(tokens)>1:
        mpi = tokens[-2]
        omp = tokens[-1]
        if mpi.isdigit() and omp.isdigit():
            test = full_test.rsplit('-',2)[0]
        #end if
    #end if
except:
    error('shorten')
#end try


# get known label sets
try:
    label_sets = create_label_sets()
except:
    error('create_labels')
#end try

# get labels from known sets
try:
    labels = []
    for label,label_set in label_sets.iteritems():
        if test in label_set or full_test in label_set:
            labels.append(label)
        #end if
    #end for
except:
    error('find_labels')
#end try

# make a ctest list of the labels
try:
    ctest_labels = ''
    for label in labels:
        ctest_labels += label+';'
    #end for
    ctest_labels = ctest_labels.rstrip(';')
except:
    error('ctest_list')
#end try


# print out the label list for ctest to read
sys.stdout.write(ctest_labels)
sys.exit(0)

