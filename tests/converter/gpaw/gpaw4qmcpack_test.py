#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

# unit test GPAW4QMCPACK

# TODO: implement and test converter generated QMCPACK input

def run_gpaw4qmcpack_test(test_name, g4q_exe, h5diff_exe):
    okay = True

    # Example invocation of converter
    #gpaw4qmcpack infile.gpw outfile.h5 --density

    infile = 'restart.gpw'
    outfile = 'test.orbs.h5'
    cmd = [g4q_exe,infile,outfile]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = p.communicate()

    file_out = open('stdout.txt', 'w')
    file_out.write(stdout)
    file_out.close()
    if len(stderr) > 0 :
        file_err = open('stderr.txt', 'w')
        file_err.write(stderr)
        file_err.close()

    ret = p.returncode

    if ret != 0:
        print("Return code nonzero: ", ret)
        okay = False

    if len(stderr.strip()) != 0:
        # some MPI output on stderr is okay
        # TODO - more general way of checking acceptable stderr strings
        if not stderr.startswith('Rank'):
            print("Stderr not empty")
            print(stderr)
            okay = False

    ret = os.system(h5diff_exe + ' -d 0.000001 gold.orbs.h5 test.orbs.h5')
    # if it's okay up to this point
    if ret==0 and okay:
        print("  pass")
        return True
    else:
        print("h5diff reported a difference")
        print("  FAIL")
        return False

    return okay


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test gpaw4qmcpck')
    parser.add_argument('test_name',
                        default='test_Si_diamond',
                        help='Name of test to run (name of directory)')
    parser.add_argument('--exe',
                        default='gpaw4qmcpack',
                        help='Location of gpaw4qmcpack executable')
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

    ret = run_gpaw4qmcpack_test(test_dir, args.exe, args.h5diff)

    os.chdir(curr_dir)

    if ret:
        sys.exit(0)
    sys.exit(1)

