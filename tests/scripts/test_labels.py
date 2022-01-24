#! /usr/bin/env python3

import sys
import traceback


def error(msg=None):
    if msg is None:
        sys.stdout.write(traceback.format_exc())
    else:
        sys.stdout.write('error: '+msg+'\n')
    #end if
    sys.exit(1)
#end def error



# label descriptions
#
#   main classification
#     stable   - test is expected to pass
#     unstable - test may fail (definitely or with some probability)
#     (stable/unstable are exclusive and comprehensive for all tests)
#
#   reason for instability
#     bug           - test failure definitely a bug and QMCPACK needs to be fixed
#     unsupported   - test fails because feature is not supported and was not before
#     poor_test     - test fails because test is insufficient (false positive)
#     cause_unknown - cause/nature of failure is not known and needs investigation
#     (bug/unsupported/poor/unknown are comprehensive for failures)
#
#   failure type classification of unstable tests
#     hard_fail              - ungraceful crash (segfault, etc) 
#     abort                  - controlled abort 
#     check_fail             - fails rigorous non-statistical check (e.g. checksum)
#     reference_stat_fail    - fails separate reference level check by > 5 sigma
#     definite_stat_fail     - statistical failure by > 5 sigma
#     intermittent_stat_fail - statistical failure by < 5 sigma
#     deterministic_fail     - fails non-rigorous deterministic check
#     (failure types are comprehensive for failures)
#  
#   test quality categories
#     good_test       - test well designed: test fail means there is a bug  
#                       any failure means: action needed to fix QMCPACK
#     poor_test       - test poorly designed and failure likely reflects this
#                       regardless of pass/fail: action needed to improve test 
#     quality_unknown - test quality has not yet been assessed
#                       code failure means: action needed to fix QMCPACK
#                       statistical failure means: prioritize quality assessment
#     (good/poor/unknown are exclusive and comprehensive for all tests)
#
#   put failing tests into the following categories (these determine the others)
#     unsupported
#     hard_fail
#     abort
#     check_fail
#     reference_stat_fail
#     definite_stat_fail
#     intermittent_stat_fail
#     deterministic_fail
#     poor_test
#
# main actions for developers and corresponding categories
#   fix known bugs in QMCPACK  - identified by label "bug"
#   fix poor tests             - identified by label "poor_test"
#   assess test quality        - identified by label "quality_unknown"
#   investigate failure causes - identified by label "cause_unknown"
#                                run reference-level statistical test in this case

def create_label_sets():

    # primary sets, all others derived from these
    unsupported            = set() # failure categories may depend on build
    hard_fail              = set()
    abort                  = set()
    check_fail             = set()
    reference_stat_fail    = set()
    definite_stat_fail     = set()
    intermittent_stat_fail = set()
    deterministic_fail     = set()
    good_test              = set() # test quality does not depend on build
    poor_test              = set()

    # universal failures (irrespective of build)
    hard_fail |= set([
        # https://github.com/QMCPACK/qmcpack/issues/848
        'developer-heg_54_J2rpa',
        ])

    abort |= set([
        # https://github.com/QMCPACK/qmcpack/issues/998
        'long-LiH_dimer_pp-vmc_hf_sdj_estimator_spindensity',
        'short-LiH_dimer_pp-vmc_hf_sdj_estimator_spindensity',

        # https://github.com/QMCPACK/qmcpack/issues/1040
        'short-diamondC-afqmc_hyb_nn2',
        'short-diamondC-afqmc_incmf_nn2',
        'short-diamondC-afqmc_nn2',
        ])

    check_fail |= set([
        # https://github.com/QMCPACK/qmcpack/issues/934
        'short-diamondC_1x1x1_pp-dmc-estimator-density',
        'short-diamondC_1x1x1_pp-dmc-estimator-spindensity',
        'long-diamondC_1x1x1_pp-dmc-estimator-density',
        'long-diamondC_1x1x1_pp-dmc-estimator-spindensity',
        ])

    definite_stat_fail |= set([
        # https://github.com/QMCPACK/qmcpack/issues/995
        # https://github.com/QMCPACK/qmcpack/issues/1052
        'short-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm-4-4',
        'short-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm-4-4',
        'long-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm-4-4',
        'long-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm-4-4',

        # https://github.com/QMCPACK/qmcpack/issues/982
        'long-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi',
        'long-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell',
        'long-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell',
        'long-LiH_dimer_ae-vmc_hf_noj_estimator_energydensity_voronoi',
        'long-LiH_dimer_ae-vmc_hf_sdj_estimator_energydensity_voronoi',
        'short-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi',
        'short-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell',
        'short-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell',
        'short-LiH_dimer_ae-vmc_hf_noj_estimator_energydensity_voronoi',
        'short-LiH_dimer_ae-vmc_hf_sdj_estimator_energydensity_voronoi',
        ])

    intermittent_stat_fail |= set([
        'long-diamondC_1x1x1_pp-vmc-dmc-allp_sdj',
        'long-diamondC_2x1x1_pp-dmc-reconf_sdj',
        'long-diamondC_2x1x1_pp-vmc_sdj',
        'long-H2O_dimer_sep_pp-j3_dmc_la',
        'short-bccH_1x1x1_ae-csvmc-all_sdj',
        'short-bccH_1x1x1_ae-csvmc-all-nodrift_sdj',
        'short-bccH_1x1x1_ae-csvmc-pbyp-nodrift_sdj',
        'short-bccH_1x1x1_ae-dmc-all_sdj',
        'short-bccH_1x1x1_ae-vmc-all-nodrift_sdj',
        'short-C2_pp-msdj_vmc',
        'short-diamondC_1x1x1_pp-vmc-dmc-allp_sdj',
        'short-diamondC_2x1x1_pp-dmc-reconf_sdj',
        'short-H4-opt-adaptive',
        'long-LiH_dimer_ae_qp-vmc_hf_noj',
        'short-LiH_dimer_pp-vmc_hf_sdj_hdf5',
        'short-LiH_dimer_pp-vmc_hf_sdj_xml',
        'short-LiH_solid_1x1x1_pp-x-dmcnl-hf_noj',
        'short-LiH_solid_1x1x1_hybridrep_pp-x-vmc_hf_noj',
        'short-diamondC_1x1x1_pp-vmc_sdj_kspace-1-16',
        'short-diamondC_2x1x1_pp-dmc_sdj-1-16',
        'short-diamondC_2x1x1_hybridrep_pp-vmc_sdj',
        'short-bccH_1x1x1_ae-dmc_sdj',
        'short-NiO_a4_e48_pp-dmc-TMv1v3_sdj',
        'long-heg_14_gamma-sj-1-16',            
        'short-chn_ae_cuspCorrection-vmc',
        'short-li2_sto-sj_dmc',
        'short-LiH_dimer_ae_qp-vmc_hf_noj',
        'short-LiH_dimer_ae_pyscf-vmc_hf_noj',
        'short-LiH_ae-vmc_msdj-1-16',
        'short-LiH_ae-vmc_msdj_noj-1-16',        
        'vmc_short_C2_pp_msdj-H5',    
        ])

    poor_test |= set([
        'short-bccH_2x2x2_ae-deriv',
        'short-bccH_2x2x2_ae-gamma-deriv',
        'short-bccH_2x2x2_ae-grad_lap',
        'short-bccH_3x3x3_ae-deriv',
        'short-bccH_3x3x3_ae-gamma-deriv',
        'short-bccH_3x3x3_ae-grad_lap',
        'short-bccH_3x3x3_ae-not_orth-deriv',
        ])

    # aos specific (but general otherwise)
    if aos:
        intermittent_stat_fail |= set([
            'short-H4-orb-opt-dmc',
            ])
        check_fail |= set([
            'short-bccH_2x2x2_ae-deriv',
            'short-bccH_2x2x2_ae-gamma-deriv',
            'short-bccH_2x2x2_ae-grad_lap',
            'short-bccH_3x3x3_ae-deriv',
            'short-bccH_3x3x3_ae-gamma-deriv',
            'short-bccH_3x3x3_ae-grad_lap',
            'short-bccH_3x3x3_ae-not_orth-deriv',
            ])
    #end if

    # soa specific (but general otherwise)
    if soa:
        hard_fail |= set([
            'long-heg_14_gamma-sjb',
            'short-heg_14_gamma-sjb',
            ])
        abort |= set([
            'long-c_no-hf_vmc',
            'long-c_no-sj_dmc',
            'short-bccH_2x2x2_ae-deriv',
            'short-bccH_3x3x3_ae-deriv',
            'short-bccH_3x3x3_ae-grad_lap',
            'short-c_no-hf_vmc',
            'short-c_no-sj_dmc',
            'short-H2-FDLR',
            'short-H2-orb-opt',
            'short-H4-FDLR',
            'short-H4-orb-opt',
            'short-H4-orb-opt-dmc',
            'short-bccH_2x2x2_ae-gamma-deriv',
            'short-bccH_2x2x2_ae-grad_lap',
            'short-bccH_3x3x3_ae-gamma-deriv',
            'short-bccH_3x3x3_ae-not_orth-deriv',
            ])
        check_fail |= set([
            'short-bccH_2x2x2_ae-deriv',
            'short-bccH_3x3x3_ae-deriv',
            'short-bccH_3x3x3_ae-grad_lap',
            ])
        intermittent_stat_fail |= set([
            'short-H4-opt-cslinear-rescale',
            ])
    #end if

    # gpu specific (but general otherwise)
    if gpu:
        check_fail |= set([
            'restart-1-16',
            'short-diamondC_1x1x1_pp-dmc-estimator-density',
            'short-diamondC_1x1x1_pp-dmc-estimator-spindensity',
            ])
        poor_test |= set([
            'restart-1-16',
            ])
        unsupported |= set([
            'estimator-skinetic',
            'short-afqmc-N2_vdz',
            'short-diamondC-afqmc',
            'short-diamondC-afqmc_hyb',
            'short-diamondC-afqmc_hyb_nn2',
            'short-diamondC-afqmc_incmf',
            'short-diamondC-afqmc_incmf_nn2',
            'short-diamondC-afqmc_nn2',
            'short-diamondC_1x1x1_hybridrep_pp-vmc_sdj',
            'short-diamondC_1x1x1_pp-vmc_sdj_kspace',
            'short-diamondC_2x1x1_hybridrep_pp-vmc_sdj',
            'short-LiH_solid_1x1x1_hybridrep_pp-x-vmc_hf_noj',
            'short-monoO_1x1x1_pp-vmc_sdj3',
            'short-NiO_a4_e48-batched_pp-vmc_sdj3',
            'short-NiO_a4_e48-hybridrep-batched_pp-vmc_sdj3',
            'short-NiO_a4_e48-hybridrep-pp-vmc_sdj3',
            'short-NiO_a4_e48_pp-vmc_sdj3',
            ])
    #end if

    # real specific (but general otherwise)
    if real:
        None
    #end if

    # complex specific (but general otherwise)
    if comp:
        None
    #end if
    
    # mixed precision specific (but general otherwise)
    if mixed:
        abort |= set([
            'short-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell',
            'short-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi',
            'short-LiH_dimer_ae-vmc_hf_noj_estimator_energydensity_voronoi',
            'short-LiH_dimer_pp-vmc_hf_sdj_estimator_spindensity',
            ])
    #end if


    # finer details based on build

    # AoS issues
    if aos and cpu and real and full:
        abort |= set([
            'long-NiO_afm-afqmc',
            'short-c_no-hf_vmc-4-4',
            ])
        definite_stat_fail |= set([
            'long-monoO_1x1x1_pp-vmc_sdj3',     # flux fails
            'short-H4-opt-cslinear-rescale',    # energy fails
            'short-LiH_dimer_ae_qp-vmc_hf_noj', # energy fails
            ])
    elif aos and cpu and real and mixed:
        definite_stat_fail |= set([
            'short-c_no-sj_dmc', # variance fails
            ])
    elif aos and cpu and comp and full:
        None
    elif aos and cpu and comp and mixed:
        None
    elif aos and gpu and real and full:
        None
    elif aos and gpu and real and mixed:
        None
    elif aos and gpu and comp and full:
        None
    elif aos and gpu and comp and mixed:
        None
    #end if

    # SoA issues
    if soa and comp and cpu:
        abort |= set([
            'short-diamondC_1x1x1_pp-vmc_gaussian_sdj',
            'short-diamondC_2x1x1_pp-vmc_gaussian_sdj',
            'long-diamondC_2x1x1_pp-dmc_gaussian_sdj',
            ])
    #end if
    if soa and cpu and real and full:
        abort |= set([
            'long-NiO_afm-afqmc',
            ])
        check_fail |= set([
            'short-diamondC_2x1x1_pp-vmc_gaussian_sdj', # ionion fails
            ])
        definite_stat_fail |= set([
            'long-C2_pp-msdj_vmc',
            'long-diamondC_1x1x1_pp-vmc-dmc-allp_sdj',
            ])
    elif soa and cpu and real and mixed:
        check_fail |= set([
            'short-diamondC_2x1x1_pp-vmc_gaussian_sdj', # ionion fails
            ])
    elif soa and cpu and comp and full:
        None
    elif soa and cpu and comp and mixed:
        None
    elif soa and gpu and real and full:
        None
    elif soa and gpu and real and mixed:
        None
    elif soa and gpu and comp and full:
        None
    elif soa and gpu and comp and mixed:
        None
    #end if

    # derive all other sets from the primary sets

    # don't require failure type to be enumerated for unsupported
    #   just assume abort is used for now
    abort |= unsupported - abort - hard_fail

    # weak failures: insufficient to tell if there is a bug
    weak_fail   = intermittent_stat_fail \
                | deterministic_fail

    # strong failures: sufficient to tell if there is a bug
    strong_fail = hard_fail           \
                | abort               \
                | check_fail          \
                | reference_stat_fail \
                | definite_stat_fail 

    fail = weak_fail | strong_fail

    # a test is unstable if it fails for any reason
    unstable = fail | poor_test | unsupported

    # a failure is a bug if it is a strong failure (strong fail => bug)
    #   intermittent failures imply bugs only with verified good test data
    #   unsupported features are not bugs
    bug = (strong_fail | (good_test & intermittent_stat_fail)) - unsupported

    # a failing test needs to be followed up with a reference-level statistical test
    #   if it is insufficient on its own to imply a bug
    # currently this includes intermittent failures with test data of unverified quality
    #   and deterministic failures  
    cause_unknown = unstable - bug - unsupported - poor_test


    positive_label_sets = dict(
        # main classification
        unstable               = unstable,
        # reason for failure
        bug                    = bug,
        unsupported            = unsupported,
        poor_test              = poor_test,
        cause_unknown          = cause_unknown,
        # failure type classification
        hard_fail              = hard_fail,
        abort                  = abort,
        check_fail             = check_fail,
        reference_stat_fail    = reference_stat_fail,
        definite_stat_fail     = definite_stat_fail,
        intermittent_stat_fail = intermittent_stat_fail,
        deterministic_fail     = deterministic_fail,
        # test quality categories (also poor_test)
        good_test              = good_test,
        )

    negative_label_sets = dict(
        # main classification
        #stable          = unstable, # can't add this because ctest regex matches "unstable"
        # test quality categories
        quality_unknown = good_test | poor_test,
        )

    return positive_label_sets,negative_label_sets
#end def create_label_sets



def check_exclusive(*labels):
    # fast check of mutual exclusivity
    exclusive = True
    lsets = [positive_label_sets[l] for l in labels]
    combined = set(lsets[0])
    for lset in lsets[1:]:
        exclusive &= len(lset & combined)==0
        if not exclusive:
            break
        #end if
        combined |= lset
    #end for
    # slower process to report errors
    if not exclusive:
        for i in range(len(labels)):
            for j in range(i):
                overlap = lsets[i] & lsets[j]
                if(len(overlap)>0):
                    error('label sets {0} and {1} are not mutually exclusive\n  overlap: {2}'.format(labels[i],labels[j],sorted(overlap)))
                #end if
            #end for
        #end for
    #end if
#end def check_exclusive



def check_comprehensive(full_label,*other_labels):
    full = set(positive_label_sets[full_label])
    for label in other_labels:
        full -= positive_label_sets[label]
    #end for
    if len(full)>0:
        error('the following sets are not comprehensive:\n  {0}\nthese should equate to set: {1}\nbut the following tests are not labeled: {2}'.format(other_labels,full_label,sorted(full)))
    #end if
#end def check_comprehensive



def check_positive_label_sets(positive_label_sets):

    check_exclusive('good_test','poor_test')

    check_comprehensive('unstable',
                        'bug',
                        'unsupported',
                        'poor_test',
                        'cause_unknown'
                        )

    check_comprehensive('unstable',
                        'hard_fail',
                        'abort',
                        'check_fail',
                        'reference_stat_fail',
                        'definite_stat_fail',
                        'intermittent_stat_fail',
                        'deterministic_fail'
                        )
#end def check_positive_label_sets



#================#
# main execution #
#================#

# extract test name and build flags from args
try:
    full_test,qmc_cuda,qmc_complex,qmc_mixed = sys.argv[1:]
    qmc_cuda    = qmc_cuda=='1'
    qmc_complex = qmc_complex=='1'
    qmc_mixed   = qmc_mixed=='1'
    cpu   = not qmc_cuda
    gpu   = qmc_cuda
    aos   = False
    soa   = True
    real  = not qmc_complex
    comp  = qmc_complex
    full  = not qmc_mixed
    mixed = qmc_mixed
except:
    error("command line args: " + str(sys.argv))
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
    error()
#end try


# get known label sets
try:
    positive_label_sets,negative_label_sets = create_label_sets()
except:
    error()
#end try


# check positive label sets
try:
    check_positive_label_sets(positive_label_sets)
except:
    error()
#end try


# get labels from known sets
try:
    labels = []
    for label,label_set in positive_label_sets.items():
        if test in label_set or full_test in label_set:
            labels.append(label)
        #end if
    #end for
    for label,label_set in negative_label_sets.items():
        if test not in label_set and full_test not in label_set:
            labels.append(label)
        #end if
    #end for
except:
    error()
#end try


# directly apply labels based on test names/pattern matching
try:
    if test.startswith('deterministic'):
        labels.append('deterministic')
    #end if
except:
    error()
#end try


try:
    if test.startswith('qe-'):
        labels.append('converter')
    #end if
except:
    error()
#end try


# mark all statistical tests as unstable
#   some of these work well, and others do not
#   cause of intermittent statistical failures needs further investiagion
try: 
    if test.startswith('short-') or test.startswith('long-') or test.startswith('estimator-'):
        if 'unstable' not in labels:
            labels.append('unstable')
        #end if
    #end if
except:
    error()
#end try


# make a ctest list of the labels
try:
    ctest_labels = ';'.join(labels)
except:
    error()
#end try


# print out the label list for ctest to read
sys.stdout.write(ctest_labels)
sys.exit(0)

