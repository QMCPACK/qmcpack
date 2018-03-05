from __future__ import print_function

import argparse
import difflib
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
#  gold.wfnoj.xml - expected version of output from convert4qmc
#  'expect_fail.txt' - if present, converter should fail

# Only the wavefunction conversion is tested currently.
# Structure conversion (gold.structure.xml) is not tested.

def compare(gold_file,test_file):
        if not filecmp.cmp(gold_file, test_file):
            print("Gold file comparison failed")
            with open(gold_file, 'r') as f_gold:
                gold_lines = f_gold.readlines()
            with open(test_file, 'r') as f_test:
                test_lines = f_test.readlines()

            diff = difflib.unified_diff(gold_lines, test_lines, fromfile=gold_file, tofile=test_file)
            diff_line_limit = 200
            for i,diff_line in enumerate(diff):
                print(diff_line,end="")
                if i > diff_line_limit:
                    print('< diff truncated due to line limit >')
                    break

	    return False
	else:
	    return True


def run_test(test_name, c4q_exe, h5diff_exe, conv_inp, gold_file, expect_fail, extra_cmd_args,code):
    okay = True

    # Example invocation of converter
    #convert4qmc -nojastrow -prefix gold -gamessAscii be.out

    cmd = c4q_exe.split()
    if code=='pyscf':
        cmd.extend(['-nojastrow', '-prefix', 'test', '-pyscf', conv_inp])
    if code=='qp':
        cmd.extend(['-nojastrow', '-prefix', 'test', '-QP',  conv_inp])
    if code=='gamess':
        cmd.extend(['-nojastrow', '-prefix', 'test', '-gamessAscii', conv_inp])

    for ex_arg in extra_cmd_args:
        if ex_arg == '-ci':
            cmd.extend(['-ci', conv_inp])
        else:
            cmd.append(ex_arg)

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
	else:
            if (code != 'pyscf'): 
                if '-hdf5' in extra_cmd_args:
                   ret = os.system(h5diff_exe + ' gold.orbs.h5 test.orbs.h5')
                   if ret==0:
                      print("  pass")
                      return True
                   else:
                      print("  FAIL")
                      return False
            test_file = gold_file.replace('gold', 'test')
            okay = compare(gold_file, test_file)

    if okay:
        print("  pass")
    else:
        print("  FAIL")

    return okay

def read_extra_args():
    extra_cmd_args = []
    if os.path.exists('cmd_args.txt'):
        with open('cmd_args.txt', 'r') as f_cmd_args:
            for line in f_cmd_args:
                line = line.strip()
                if line.startswith('#'):
                    continue
                extra_cmd_args.append(line)
    return extra_cmd_args


def run_one_converter_test(c4q_exe, h5diff_exe):
    code='gamess'
    if os.path.exists('pyscf'):
       code='pyscf'
    if os.path.exists('quantum_package'):
       code='qp'
    
    test_name = os.path.split(os.getcwd())[-1]
    
    if code=='gamess': 
       conv_input_files = glob.glob('*.out')

    if code=='pyscf': 
       conv_input_files = glob.glob('*.h5')

    if code=='qp': 
       conv_input_files = glob.glob('*.dump')

    if len(conv_input_files) != 1:
        print("Unexpected number of inputs files (should be 1): ",
              len(conv_input_files))
        return False
    conv_input_file = conv_input_files[0]

    extra_cmd_args = read_extra_args()

    expect_fail = os.path.exists('expect_fail.txt')
    gold_file = 'gold.wfnoj.xml'
    if expect_fail:
        gold_file = ''
    else:
        if not os.path.exists(gold_file):
            print("Gold file missing")
            return False
    return run_test(test_name, c4q_exe, h5diff_exe, conv_input_file, gold_file,
                    expect_fail, extra_cmd_args,code)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test convert4qmc')
    parser.add_argument('test_name',
                        help='Name of test to run (name of directory)')
    parser.add_argument('--exe',
                        default='convert4qmc',
                        help='Location of convert4qmc executable')
    parser.add_argument('--h5diff',
                        default='h5diff',
                        help='Location of h5diff executable')
    args = parser.parse_args()

    test_dir = args.test_name
    if not os.path.exists(test_dir):
        print("Test not found: ", test_dir)
        sys.exit(1)

    curr_dir = os.getcwd()
    os.chdir(test_dir)

    ret = run_one_converter_test(args.exe, args.h5diff)

    os.chdir(curr_dir)

    if ret:
        sys.exit(0)
    sys.exit(1)
