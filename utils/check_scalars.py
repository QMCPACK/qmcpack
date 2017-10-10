#!/usr/bin/env python
from __future__ import print_function

# Statical error checking code for use by testing framework
# Jaron Krogel/ORNL

# To maximize portability, only standard Python modules should
# be used (that is, no numpy)

import os
from optparse import OptionParser
import math


# standalone definition of error function from Abramowitz & Stegun
# credit: http://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/
def erf(x):
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of x
    sign = 1
    if x < 0:
        sign = -1
    x = abs(x)

    # A & S 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)

    return sign*y
#end def erf


# Returns failure error code to OS.
# Explicitly prints 'fail' after an optional message.
def exit_fail(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: fail')
    exit(1)
#end def exit_fail


# Returns success error code to OS.
# Explicitly prints 'pass' after an optional message.
def exit_pass(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: pass')
    exit(0)
#end def exit_pass

def compute_mean(v):
    if len(v) == 0:
        return 0.0
    return sum(v)/len(v)


def compute_variance(v):
    if len(v) == 0:
        return 0.0
    mean = compute_mean(v)
    return sum([(x-mean)**2 for x in v])/len(v)



# Calculates the mean, variance, errorbar, and autocorrelation time
# for a 1-d array of statistical data values.
# If 'exclude' is provided, the first 'exclude' values will be 
# excluded from the analysis.
def simstats(x,exclude=None):
    if exclude!=None:
        x = x[exclude:]
    #end if

    N    = len(x)
    mean = compute_mean(x)
    var = compute_variance(x)

    i=0          
    tempC=0.5
    kappa=0.0
    if abs(var)<1e-15:
        kappa = 1.0
    else:
        ovar=1.0/var
        while (tempC>0 and i<(N-1)):
            kappa=kappa+2.0*tempC
            i=i+1
            tempC = ovar/(N-i)*sum([ (x[idx]-mean)*(x[idx+i]-mean) for idx in range(N-i) ])
        #end while
        if kappa == 0.0:
            kappa = 1.0
        #end if
    #end if
    Neff=(N+0.0)/(kappa+0.0)
    if (Neff == 0.0):
        Neff = 1.0
    #end if
    error=math.sqrt(var/Neff)

    return (mean,var,error,kappa)
#end def simstats



# Reads command line options.
# For example:
#   check_scalars.py  --ns 2  -p Li  -s '1 3' -e 10  --le '-7.478011 0.000035  -7.478059 0.000035'
# This invocation expects two scalar files to be present:
#   Li.s001.scalar.dat  Li.s003.scalar.dat
# The local energy ('le') will be computed excluding ('e') 10 blocks each.
# the test passes if both of the following are true:
#   |local_energy_mean_of_series_1 - (-7.478011)| < 2*0.000035
#   |local_energy_mean_of_series_3 - (-7.478059)| < 2*0.000035
# Note that the factor of 2 is specified by 'ns', the allowed number of sigma deviations for the test.
# The inputted errorbar (sigma) to check against (0.000035) should satisfy the following:
#   Let err_ref be the error bar of the reference solution.
#   Let err_comp be the expected errorbar of the completed test run.
#   Then the provided errorbar/sigma should be sqrt( err_ref**2 + err_comp**2 ).
def read_command_line():
    try:

        parser = OptionParser(
            usage='usage: %prog [options]',
            add_help_option=False,
            version='%prog 0.1'
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
        parser.add_option('-e','--equilibration',dest='equilibration',
                          default='0',
                          help='Equilibration length in blocks (default=%default).'
                          )
        parser.add_option('-n','--nsigma',dest='nsigma',
                          default='3',
                          help='Sigma requirement for pass/fail (default=%default).'
                          )

        quantities = dict(
            ar   = 'AcceptRatio',
            le   = 'LocalEnergy',
            va   = 'Variance',
            ke   = 'Kinetic',
            lp   = 'LocalPotential',
            ee   = 'ElecElec',
            cl   = 'Coulomb',
            ii   = 'IonIon',
            lpp  = 'LocalECP',
            nlpp = 'NonLocalECP',
            mpc  = 'MPC',
            kec  = 'KEcorr',
            bw   = 'BlockWeight',
            ts   = 'TotalSamples',
            fl   = 'Flux',
#now for some RMC estimators
            ke_m = "Kinetic_m",
            ke_p = "Kinetic_p",
            ee_m = "ElecElec_m",
            ee_p = "ElecElec_p",
            lp_p = "LocalPotential_pure",
#and some CSVMC estimators
            le_A = "LocEne_0",
            le_B = "LocEne_1",
            dle_AB = "dLocEne_0_1",
            ii_A = "IonIon_0",
            ii_B = "IonIon_1",
            dii_AB = "dIonIon_0_1",
            ee_A = "ElecElec_0",
            ee_B = "ElecElec_1",
            dee_AB = "dElecElec_0_1"
            )

        for qshort in sorted(quantities.keys()):
            qlong = quantities[qshort]
            parser.add_option('--'+qshort,'--'+qlong,dest=qlong,
                              default=None,
                              help='Reference value and errorbar for '+qlong+' (one value/error pair per series).'
                              )
        #end for

        options,files_in = parser.parse_args()

        if options.help:
            print('\n'+parser.format_help().strip())
            exit()
        #end if

        options.series         = [int(a) for a in options.series.split()]
        options.equilibration  = [int(a) for a in options.equilibration.split()]
        if len(options.series)>0 and len(options.equilibration)==1:
            options.equilibration = len(options.series)*[options.equilibration[0]]
        #end if
        options.nsigma         = float(options.nsigma)

        quants_check = []
        for q in quantities.values():
            v = options.__dict__[q]
            if v!=None:
                vref = [float(a) for a in v.split()]
                if len(vref)!=2*len(options.series):
                    exit_fail('must provide one reference value and errorbar for '+q)
                #end if
                options.__dict__[q] = vref
                quants_check.append(q)
            #end if
        #end for
    except Exception as e:
        exit_fail('error during command line read:\n'+str(e))
    #end try

    return options,quants_check
#end def read_command_line



# Reads scalar.dat files and performs statistical analysis.
def process_scalar_files(options,quants_check):
    values = dict()

    try:
        ns = 0
        for s in options.series:
            svals = dict()
            scalar_file = options.prefix+'.s'+str(int(s)).zfill(3)+'.scalar.dat'

            if os.path.exists(scalar_file):
                fobj = open(scalar_file,'r')
                quantities = fobj.readline().split()[2:]
                rawdata_list = []
                for line in fobj:
                    vals = line.strip().split()[1:]
                    fp_vals = [float(v) for v in vals]
                    rawdata_list.append(fp_vals)
                fobj.close()

                rawdata = [[]]
                if len(rawdata_list) == 0:
                    exit_fail('scalar file has no data: '+scalar_file)
                # end if

                # transpose
                ncols = len(rawdata_list[0])
                rawdata = [[] for i in range(ncols)]
                for line in rawdata_list:
                    for i in range(ncols):
                        rawdata[i].append(line[i])
                    # end for
                # end for


                equil = options.equilibration[ns]

                data  = dict()
                stats = dict()
                for i in range(len(quantities)):
                    q = quantities[i]
                    d = rawdata[i]
                    data[q]  = d
                    stats[q] = simstats(d,equil)
                #end for

                if 'LocalEnergy_sq' in data and 'LocalEnergy' in data:
                    v = [(le_sq-le**2) for le_sq,le in zip(data['LocalEnergy_sq'],data['LocalEnergy'])]
                    stats['Variance'] = simstats(v,equil)
                #end if

                if 'BlockWeight' in data:
                    ts = sum(data['BlockWeight'])
                    stats['TotalSamples'] = (ts,0.0,0.0,1.0) # mean, var, error, kappa
                #end if

                for q in quants_check:
                    if q in stats:
                        mean,var,error,kappa = stats[q]
                        svals[q] = mean,error
                    else:
                        exit_fail('{0} is not present in file {1}'.format(q,scalar_file))
                    #end if
                #end for
            else:
                exit_fail('scalar file does not exist: '+scalar_file)
            #end if

            values[s] = svals
            ns += 1
        #end for        
    except Exception as e:
        exit_fail('error during scalar file processing:\n'+str(e))
    #end try

    return values
#end def process_scalar_files


# Checks computed values from scalar.dat files 
# against specified reference values.
passfail = {True:'pass',False:'fail'}
def check_values(options,quants_check,values):
    success = True
    msg = ''

    try:
        ns = 0
        for s in options.series:
            msg+='Tests for series {0}\n'.format(s)
            for q in quants_check:
                msg+='  Testing quantity: {0}\n'.format(q)

                ref = options.__dict__[q]
                mean_ref  = ref[2*ns]
                error_ref = ref[2*ns+1]
                mean_comp,error_comp = values[s][q]

                quant_success = abs(mean_comp-mean_ref) <= options.nsigma*error_ref

                success &= quant_success

                delta = mean_comp-mean_ref
                delta_err = math.sqrt(error_comp**2+error_ref**2)

                msg+='    reference mean value     : {0: 12.8f}\n'.format(mean_ref)
                msg+='    reference error bar      : {0: 12.8f}\n'.format(error_ref)
                msg+='    computed  mean value     : {0: 12.8f}\n'.format(mean_comp)
                msg+='    computed  error bar      : {0: 12.8f}\n'.format(error_comp)
                msg+='    pass tolerance           : {0: 12.8f}  ({1: 12.8f} sigma)\n'.format(options.nsigma*error_ref,options.nsigma)
                if error_ref > 0.0:
                    msg+='    deviation from reference : {0: 12.8f}  ({1: 12.8f} sigma)\n'.format(delta,delta/error_ref)
                msg+='    error bar of deviation   : {0: 12.8f}\n'.format(delta_err)
                if error_ref > 0.0:
                    msg+='    significance probability : {0: 12.8f}  (gaussian statistics)\n'.format(erf(abs(delta/error_ref)/math.sqrt(2.0)))
                msg+='    status of this test      :   {0}\n'.format(passfail[quant_success])
            #end for
            ns+=1
        #end for
    except Exception as e:
        exit_fail('error during value check:\n'+str(e))
    #end try

    return success,msg
#end def check_values



# Main execution
if __name__=='__main__':
    # Read and interpret command line options.
    options,quants_check = read_command_line()

    # Compute means of desired quantities from scalar.dat files.
    values = process_scalar_files(options,quants_check)

    # Check computed means agains reference solutions.
    success,msg = check_values(options,quants_check,values)

    # Pass success/failure exit codes and strings to the OS.
    if success:
        exit_pass(msg)
    else:
        exit_fail(msg)
    #end if
#end if
