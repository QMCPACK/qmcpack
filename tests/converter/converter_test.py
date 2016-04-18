from __future__ import print_function

import argparse
import filecmp
import glob
import os
import subprocess
import sys

# Test convert4qmc for GAMESS output files

# Each test directory contains
# *.inp - Gamess input file.  Not strictly necessary for this test, but useful for reproducing the run
# *.out - Gamess output file.  Serves as input to convert4qmc
# One of the following
#  gold.Gaussian-G2.xml - expected version of output from convert4qmc
#  'expect_fail.txt' - if present, converter should fail

# Only the wavefunction conversion is tested currently.
# Structure conversion (gold.Gaussian-G2.ptcl.xml) is not tested.


def run_test(test_name, c4q_exe, conv_inp, gold_file, expect_fail):
    okay = True

    # Example invocation of converter
    #convert4qmc -nojastrow -prefix gold -gamessAscii be.out


    print("LD_LIBRARY_PATH = ",os.environ.get("LD_LIBRARY_PATH", "not set"))
    print("PATH = ",os.environ.get("PATH", "not set"))
    cmd = c4q_exe.split()
    print("first command (likely mpirun) = ",cmd[0])
    which_mpirun = subprocess.check_output(["which",cmd[0]])
    print("which mpirun = ",which_mpirun)
    if len(cmd) > 3:
        convert_exe = cmd[3]
        convert_libs = subprocess.check_output(["ldd",convert_exe])
        print("converter libs:")
        print(convert_libs)
    
    cmd.extend(['-nojastrow', '-prefix', 'test', '-gamessAscii', conv_inp])
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    ret = p.returncode

    if expect_fail:
        if ret == 0:
            print("Return code zero, but expected failure")
    else:
        if ret != 0:
            print("Return code nonzero: ", ret)
            okay = False

        if len(stderr.strip()) != 0:
            # some MPI output on stderr is okay
            # TODO - more general way of checking acceptable stderr strings
            if not stderr.startswith('Rank'):
                print("Stderr not emptry ")
                print(stderr)
                okay = False

        if not os.path.exists(gold_file):
            print("Gold file missing")
            okay = False

        test_file = gold_file.replace('gold', 'test')
        if not filecmp.cmp(gold_file, test_file):
            print("Gold file comparison failed")
            okay = False
        # TODO: use difflib module to print out the diff

    if okay:
        print("  pass")
    else:
        print("  FAIL")

    return okay


def run_one_converter_test(c4q_exe):
    test_name = os.path.split(os.getcwd())[-1]

    conv_input_files = glob.glob('*.out')
    if len(conv_input_files) != 1:
        print("Unexpected number of inputs files (should be 1): ",
              len(conv_input_files))
        return False
    conv_input_file = conv_input_files[0]

    expect_fail = os.path.exists('expect_fail.txt')
    gold_file = 'gold.Gaussian-G2.xml'
    if expect_fail:
        gold_file = ''
    else:
        if not os.path.exists(gold_file):
            print("Gold file missing")
            return False

    return run_test(test_name, c4q_exe, conv_input_file, gold_file,
                    expect_fail)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test convert4qmc')
    parser.add_argument('test_name',
                        help='Name of test to run (name of directory)')
    parser.add_argument('--exe',
                        default='conver4qmc',
                        help='Location of convert4qmc executable')
    args = parser.parse_args()

    test_dir = args.test_name
    if not os.path.exists(test_dir):
        print("Test not found: ", test_dir)
        sys.exit(1)

    curr_dir = os.getcwd()
    os.chdir(test_dir)

    ret = run_one_converter_test(args.exe)

    os.chdir(curr_dir)

    if ret:
        sys.exit(0)
    sys.exit(1)
