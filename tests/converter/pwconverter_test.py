#! /usr/bin/env python3

import argparse
import filecmp
import glob
import os
import subprocess
import sys

def run_test(cpw4q_exe, h5diff_exe, gold_file, conv_inp):
    okay = True
    
    # Example invocation of converter for qbox sample file:
    # convertpw4qmc foo.sample -o eshdf.h5

    cmd = cpw4q_exe.split()
    cmd.extend([conv_inp, '-o', 'eshdf.h5'])
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    ret = p.returncode

    if ret != 0:
        print ("Return code from conversion nonzero: ", ret)
        okay = False
    
    if len(stderr.strip()) != 0:
        print ("Stderr not empty ")
        print (stderr)
        san = stderr.decode("utf-8").find("Sanitizer")
        sup = stderr.decode("utf-8").find("Suppressions")
        if sup >= 0 or san >= 0 :
            print("Ignoring sanitizer suppressions")
            okay = True
        else:
            okay = False
        
    if not os.path.exists(gold_file):
        print ("Gold file missing")
        okay = False
    
    if okay:
        ret = os.system(h5diff_exe + ' -d 1e-9 ' + gold_file + ' eshdf.h5')
        if ret!=0:
            print("h5diff reported a difference")
            okay = False
        
    if okay:
        print("   pass")
    else:
        print("   FAIL")
    
    return okay
        

def run_one_converter_test(cpw4q_exe, h5diff_exe):
    code = 'qbox'
    if 'pwscf' in os.getcwd():
        code = 'pwscf'

    if code == 'qbox':
        conv_input_files = glob.glob('*.sample')

    if code == 'pwscf':
        conv_input_files = glob.glob('*.xml')
    
    if len(conv_input_files) != 1:
        print ("Unexpected number of input files (should be 1): ",
               len(conv_input_files))
        return False
    conv_input_file = conv_input_files[0]
    gold_file = 'gold.h5'
    return run_test(cpw4q_exe, h5diff_exe, gold_file, conv_input_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test convertpw4qmc')
    parser.add_argument('test_name',
                        help='Name of test to run (name of directory)')
    parser.add_argument('--exe',
                        default='convertpw4qmc',
                        help='Location of convertpw4qmc executable')
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
