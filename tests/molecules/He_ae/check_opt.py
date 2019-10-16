#!/usr/bin/env python

from __future__ import print_function

from optparse import OptionParser
import xml.etree.ElementTree as ET
import sys
import os


def exit_fail(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: fail')
    exit(1)
#end def exit_fail


def exit_pass(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: pass')
    exit(0)
#end def exit_pass


# Open the XML file and return the jud_b value
def get_jud_b(file):
    tree = ET.parse(file)
    root = tree.getroot()
    for elem in root.iterfind('wavefunction/jastrow/correlation/'):
        return [float(i) for i in elem.text.split()]
#end def get_jud_b


passfail = {True:'pass',False:'fail'}
def run_opt_test(options):


    prefix_file = options.prefix+'.s'+str(options.series).zfill(3)+'.opt.xml'    

    # Set default for the reference
    ref_file = './qmc-ref/{0}'.format(prefix_file) if options.ref is None else options.ref

    if not os.path.exists(prefix_file):
       exit_fail("Test not found:" + prefix_file)       

    if not os.path.exists(ref_file):
       exit_fail("Reference not found:" + ref_file)

    output = get_jud_b(prefix_file)
    reference = get_jud_b(ref_file)

    if len(output) != len(reference):
       exit_fail('Number of coefficient in test({0}) does not match with the reference({1})'.format(len(output),len(reference)))

    success = True
    tolerance = 1e-06
    deviation = []
   
    for i in range(len(output)):
       deviation.append(abs(float(output[i])-float(reference[i])))
       quant_success = deviation[i] <= tolerance 
       if quant_success is False:
           success &= quant_success
       #end if
    #end for
  
    msg='\n  Testing Series: {0}\n'.format(options.series)
    msg+='   reference coefficients   : {0}\n'.format(reference)
    msg+='   computed  coefficients   : {0}\n'.format(output)
    msg+='   pass tolerance           : {0: 12.6f}\n'.format(tolerance) 
    msg+='   deviation from reference : {0}\n'.format(deviation)
    msg+='   status of this test      :   {0}\n'.format(passfail[success])

    return success, msg

#end def run_opt_test


if __name__ == '__main__':

    parser = OptionParser(
       usage='usage: %prog [options]',
       add_help_option=False,
       )
    parser.add_option('-h','--help',dest='help',
                     action='store_true',default=False,
                     help='Print help information and exit (default=%default).'
                     )
    parser.add_option('-p','--prefix',dest='prefix',
                     default='qmc',
                     help='Prefix for output files (default=%default).'
                     )
    parser.add_option('-s','--series',dest='series',
                     default='0',
                     help='Output series to analyze (default=%default).'
                     )
    parser.add_option('-r','--ref',dest='ref',
                     help='Reference to check output files (default=./qmc-ref/$PREFIX).'
                     )

  
    options,files_in = parser.parse_args()

    if options.help:
     print('\n'+parser.format_help().strip())
     exit()
   #end if

    success,msg = run_opt_test(options)

    if success:
       exit_pass(msg)
    else:
       exit_fail(msg)

#end if
