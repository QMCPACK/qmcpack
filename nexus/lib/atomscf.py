##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  atomscf.py                                                        #
#    Nexus interface to Pitzer's ATOMSCF atomic HF code.             #
#                                                                    #
#    ATOMSCF is capable of HF calculations of specific atomic states #
#    in gaussian or slater type orbitals and with or without         #
#    semilocal pseudopotentials.  It has native support for exponent #
#    optimization.                                                   #
#                                                                    #
#    See Pitzer, Comp. Phys. Comm., 170 239 (2005)                   #
# http://www.sciencedirect.com/science/article/pii/S0010465505002638 #
#                                                                    #
#  Content summary:                                                  #
#    Atomscf                                                         #
#      Simulation class for ATOMSCF.                                 #
#                                                                    #
#    AtomscfInput                                                    #
#      Input class for the ATOMSCF code.                             #
#      Capable of reading/writing arbitrary ATOMSCF input files.     #
#                                                                    #
#    AtomscfAnalyzer                                                 #
#      SimulationAnalyzer class for ATOMSCF.  Essentially vacuous.   #
#                                                                    #
#    generate_atomscf_input                                          #
#      User function to create arbitrary ATOMSCF input.              #
#                                                                    #
#    generate_atomscf                                                #
#      User-facing function to create Atomscf simulation objects.    #
#                                                                    #
#    KJTables                                                        #
#      Class for parse and lookup of K and J energy coefficients     #
#      from Pitzer CPC 170 239 (2005).  Also with corrections from   #
#      2012 update.                                                  #
#                                                                    #
#====================================================================#


# a couple of usage examples
#
##! /usr/bin/env python
#
#from nexus import *
#from atomscf import generate_atomscf
#
#settings(
#    pseudo_dir    = './pseudopotentials',
#    status_only   = 0,
#    generate_only = 0,
#    sleep         = 3,
#    machine       = 'ws16',
#    )
#
#sims = []
#
#atomscf = generate_atomscf(  # example 2 from Pitzer CPC 2005
#    # nexus inputs
#    identifier = 'atomscf',
#    path       = 'O_RECP1',
#    job        = job(cores=1,serial=True),
#    # input file specification
#    input_type = 'basic',
#    #  title
#    aname      = 'Oxygen RECP pseudo with VTZ basis set',
#    #  general controls
#    nflag1     = 1, # gaussian basis                         (from default)
#    nprint     = 0, # do not print integrals                 (from default)
#    nsym       = 2, # s,p                                    (from state or basis?)
#    nvar       = 0, # no exponents being varied              (from default)
#    mxvar      = 0, #   for exp optimization                 (set automatically if not varied)
#    nextra     = 0, #   for exp optimization
#    iscale     = 0, #   for exp optimization
#    iamp       = 0, # do not print radial functions          (from default)
#    ird        = 0, # use default convergence parameters     (from default)
#    nzet       = 0, # normal default for no of exponents     (from default!)
#    nj         = 0, # no. of J-type open shell energy coeff  (from state!)
#    nk         = 2, # no. of K-type open shell energy coeff  (from state!)
#    ipnch      = 11, # if !=0, unit no. to write expon and coeff
#    nlab       = 1, # do not read labels for orb exponents
#    lcpu       = 2, # s,p nonlocal pseudopotential           (from pseudo!)
#    ileg       = 0, # unknown
#    #  basis set
#    nbas       = (5,5), # no. of s,p basis functions           (from basis!)
#    ncsh       = (1,0), # no. of closed s,p shells             (from state!)
#    nosh       = (0,1), # no. of open s,p shells               (from state!)
#    nccup      = (0,4), # open shell occ number for s,p        (from state!)
#    nlbl       = (2,2), # we have a 2s shell and a 2p shell    (from state!)
#    ipqn       = 5*[1]+5*[2], # s,p labels for basis functions (from basis!)
#    zn         = 8.0,   # nuclear charge                       (from element)
#    zscale     = 0.002, # fractional change of expon during opt(from default)
#                        # exponents for basis functions        (from basis)
#    zeta       = '''37.70  6.840  1.053  0.4163  0.1706
#                    34.57  7.760  2.282  0.7160  0.2140'''.split(),
#    kcoeff     = [(1,1,0,-8,3), # K-type energy coeff          (from state)
#                  (1,1,2,-4,3)],
#    coeff      = [[-0.02,-0.14, 0.38, 0.54, 0.21], # basis coeff (from basis)
#                  [ 0.02, 0.10, 0.31, 0.49, 0.34]],
#    #  pseudopotential
#    pseudo     = './pseudopotentials/O.RECP.VTZ.gms',
#    ppformat   = 'gamess',
#    use_none   = True,
#    )
#sims.append(atomscf)
#
#atomscf = generate_atomscf(  # example 2 again, but with all the redundancy handled by Nexus
#    # nexus inputs
#    identifier = 'atomscf',
#    path       = 'O_RECP2',
#    job        = job(cores=1,serial=True),
#    # input file specification
#    input_type = 'basic',
#    aname      = 'Oxygen RECP pseudo with VTZ basis set',
#    state      = '2s2 2p4',
#    guess      = 'random',
#    pseudo     = './pseudopotentials/O.RECP.VTZ.gms',
#    ppformat   = 'gamess',
#    ipnch      = 11,
#    )
#sims.append(atomscf)
#
#for n in range(100):
#    atomscf = generate_atomscf(  # similar to example 2, but for O BFD PP
#        # nexus inputs           #   run it many times with random inputs
#        identifier = 'atomscf',  #   to get convergence
#        path       = 'O_BFD_'+str(n),
#        job        = job(cores=1,serial=True),
#        # input file specification
#        input_type = 'basic',
#        aname      = 'Oxygen BFD pseudo with VQZ basis set',
#        state      = '2s2 2p4',
#        guess      = 'random',
#        pseudo     = './pseudopotentials/O.BFD.VQZ.gms',
#        ppformat   = 'gamess',
#        ipnch      = 11,
#        )
#    sims.append(atomscf)
##end for
#
#
#run_project(sims)
#
#


# contents of O.RECP.VTZ.gms:
#O 8 0. 0. 0.
#s 5 1.00
#1  37.70   1.00
#2  6.840   1.00
#3  1.053   1.00
#4  0.4163  1.00
#5  0.1706  1.00
#p 5 1.00
#1  34.57   1.00
#2  7.760   1.00
#3  2.282   1.00
#4  0.7160  1.00
#5  0.2140  1.00
#
#
#O-RECP GEN 2 1
# 3
#   -1.486456  1 100.0039
#   -5.766847  2  34.1980
#   -0.798420  2  10.0286
# 4
#    2.193891  0   2.1892
#    1.042944  1   4.3740
#  -16.344477  2   2.4049
#   11.216304  2   2.2479


# contents of O.BFD.VQZ.gms
#O 8    0. 0. 0.
#s 9 1.00
#1 0.125346     0.055741
#2 0.268022     0.304848
#3 0.573098     0.453752
#4 1.225429     0.295926
#5 2.620277     0.019567
#6 5.602818     -0.128627
#7 11.980245     0.012024
#8 25.616801     0.000407
#9 54.775216     -0.000076
#s 1 1.00
#1 0.224380     1.000000
#s 1 1.00
#1 0.843157     1.000000
#s 1 1.00
#1 1.351771     1.000000
#p 9 1.00
#1 0.083598     0.044958
#2 0.167017     0.150175
#3 0.333673     0.255999
#4 0.666627     0.281879
#5 1.331816     0.242835
#6 2.660761     0.161134
#7 5.315785     0.082308
#8 10.620108     0.039899
#9 21.217318     0.004679
#p 1 1.00
#1 0.148562     1.000000
#p 1 1.00
#1 0.452364     1.000000
#p 1 1.00
#1 1.106737     1.000000
#d 1 1.00
#1 0.455711     1.000000
#d 1 1.00
#1 1.344331     1.000000
#d 1 1.00
#1 4.008867     1.000000
#f 1 1.00
#1 0.876289     1.000000
#f 1 1.00
#1 2.763115     1.000000
#g 1 1.00
#1 1.759081     1.000000
#
#
#O-QMC GEN 2 1
#3
#6.00000000 1 9.29793903
#55.78763416 3 8.86492204
#-38.81978498 2 8.62925665
#1
#38.41914135 2 8.71924452




from simulation import Simulation,SimulationInput,SimulationAnalyzer


import os
from numpy import array,arange,zeros
from numpy.random import random
from generic import obj
from periodic_table import ptable
from pseudopotential import GaussianPP,GaussianBasisSet

from debug import *

class KJTables:
    kj_table_data = '''

        # s intrashell coefficients
        s              K000
        s1       2S    -1


        # p intrashell coefficients
        p              K110    K112  
        p1,p5    2P    -5/3     2/3   
        p2,p4    3P    -8/3    -4/3  
        p2,p4    1D    -8/3    52/15 
        p2,p4    1S    -8/3    32/3  
        p3       4S    -3      -6    
        p3       2D    -3       6/5   
        p3       2P    -3       6     


        # d intrashell coefficients
        d              K220    K222      K224
        d1,d9    2D    -9/5    2/7       18/35    
        d2,d8    3F    -16/5   -104/49   324/245  
        d2,d8    3P    -16/5   4         -24/5    
        d2,d8    1G    -16/5   136/49    524/245  
        d2,d8    1D    -16/5   -4/49     1224/245 
        d2,d8    1S    -16/5   48/7      432/35   
        d3,d7    4F    -21/5   -174/49   -306/245 
        d3,d7    4P    -21/5   18/7      -258/35  
        d3,d7    2H    -21/5   6/49      894/245  
        d3,d7    2G    -21/5   -94/49    1394/245 
        d3,d7    2F    -21/5   306/49    -606/245 
        d3,d7    2P    -21/5   6/49      894/245  
        d4,d6    5D    -24/5   -4        -36/5    
        d4,d6    3H    -24/5   -116/49   636/245  
        d4,d6    3G    -24/5   -16/49    136/245  
        d4,d6    3D    -24/5   124/49    -564/245 
        d4,d6    1I    -24/5   -76/49    1836/245 
        d4,d6    1F    -24/5   32/7      48/35    
        d5       6S    -5      -50/7     -90/7    
        d5       4G    -5      -150/49   -130/49  
        d5       4F    -5      90/49     -90/49   
        d5       4D    -5      -10/49    -270/49  
        d5       4P    -5      -30/7     30/7     
        d5       2I    -5      -130/49   270/49   
        d5       2H    -5      -90/49    510/49   
        d5       2P    -5      750/49    -330/49  
        d5       2S    -5      290/49    -150/49  


        # f intrashell coefficients
        f              K330    K332       K334        K336
        f1,f13   2F    -13/7   4/21       18/77       100/231      
        f2,f12   3H    -24/7   -92/63     -636/847    12500/7623   
        f2,f12   3F    -24/7   -8/63      -12/77      -200/693     
        f2,f12   3P    -24/7   100/21     156/77      -1700/231    
        f2,f12   1I    -24/7   188/63     1044/847    172300/99099 
        f2,f12   1G    -24/7   -40/21     3508/847    5800/2541    
        f2,f12   1D    -24/7   772/315    -180/77     4700/693     
        f2,f12   1S    -24/7   128/21     576/77      3200/231     
        f3,f11   4I    -33/7   -256/63    -2166/847   17800/7623   
        f3,f11   4G    -33/7   52/63      -318/847    -36100/7623  
        f3,f11   4F    -33/7   -20/21     -90/77      -500/231     
        f3,f11   4D    -33/7   248/63     78/77       -6400/693    
        f3,f11   4S    -33/7   -20/21     -90/77      -500/231     
        f3,f11   2L    -33/7   12/7       18/847      41500/11011  
        f3,f11   2K    -33/7   -116/63    1810/847    412700/99099 
        f3,f11   2I    -33/7   80/63      1614/847    172600/99099 
        f3,f11   2P    -33/7   -32/63     50/77       3400/693     
        f4,f10   5I    -40/7   -340/63    -3552/847   -5300/7623   
        f4,f10   5G    -40/7   -32/63     -1704/847   -59200/7623  
        f4,f10   5F    -40/7   -16/7      -216/77     -400/77      
        f4,f10   5D    -40/7   164/63     -48/77      -8500/693      
        f4,f10   5S    -40/7   -16/7      -216/77     -400/77        
        f4,f10   3M    -40/7   -116/63    -1032/847   538700/99099   
        f4,f10   3L    -40/7   -200/63    228/847     465200/99099   
        f4,f10   1N    -40/7   52/63      760/847     685700/99099   
        f4,f10   1K    -40/7   -328/315   2524/847    582800/99099   
        f4,f10   1F    -40/7   304/63     5604/847    -29800/7623    
        f5,f9    6H    -45/7   -344/63    -4794/847   -56800/7623    
        f5,f9    6F    -45/7   -260/63    -390/77     -6500/693      
        f5,f9    6P    -45/7   16/21      -222/77     -3800/231      
        f5,f9    4M    -45/7   -32/7      -138/77     9400/3003      
        f5,f9    4L    -45/7   -68/21     -2778/847   127900/33033   
        f5,f9    4S    -45/7   44/21      18/7        100/21         
        f5,f9    2O    -45/7   -148/63    -90/847     873700/99099   
        f5,f9    2N    -45/7   -1048/315  -1322/847   1143200/99099  
        f6,f8    7F    -48/7   -136/21    -612/77     -3400/231      
        f6,f8    5L    -48/7   -352/63    -2280/847   -263000/99099  
        f6,f8    5K    -48/7   -128/63    -4072/847   -302200/99099  
        f6,f8    5P    -48/7   -44/63     -164/77     -5300/693      
        f6,f8    5S    -48/7   -16/63     -24/77      -400/693       
        f6,f8    3O    -48/7   -296/63    -1356/847   815000/99099   
        f6,f8    3N    -48/7   -1172/315  -124/847    545500/99099   
        f6,f8    1Q    -48/7   -152/35    -516/847    145000/11011   
        f6,f8    1P    -48/7   2132/315   1724/847    500/7623       
        f7       8S    -7      -28/3      -126/11     -700/33        
        f7       6I    -7      -56/9      -630/121    -128800/14157  
        f7       6H    -7      -16/9      -390/121    -9800/1089     
        f7       6G    -7      -4/3       -982/121    -3500/363      
        f7       6F    -7      -28/9      -42/11      -700/99        
        f7       6D    -7      -256/45    -18/11      -1400/99       
        f7       6P    -7      -8         -6          0              
        f7       4N    -7      -16/3      -274/121    4200/1573      
        f7       4M    -7      -8/3       -18/121     19600/4719     
        f7       2Q    -7      -236/45    -174/121    160300/14157   
        f7       2O    -7      -44/9      -54/121     230300/14157   


        # g intrashell coefficients
        g              K440    K442      K444         K446      K448
        g1,g17   2G    -17/9   100/693   162/1001     20/99     490/1287   
        g8,g10   9G    -80/9  -6200/693  -10044/1001  -1240/99  -30380/1287
        g9       10S   -9      -900/77   -13122/1001  -180/11   -4410/143  



        # p,s intershell coefficients for states with open shells of different symmetry
        p,s                    K101
        p1,p5(2P) s1(2S)  3P   -1
        p1,p5(2P) s1(2S)  1P   3 
        p2,p4(3P) s1(2S)  4P   -2
        p2,p4(1D) s1(2S)  2D   0 
        p2,p4(3P) s1(2S)  2P   4 
        p2,p4(1S) s1(2S)  2S   0 
        p3(4S)    s1(2S)  5S   -3
        p3(2D)    s1(2S)  3D   -1
        p3(2P)    s1(2S)  3P   -1
        p3(4S)    s1(2S)  3S   5 
        p3(2D)    s1(2S)  1D   3 
        p3(2P)    s1(2S)  1P   3 


        # d,s intershell coefficients for states with open shells of different symmetry
        d,s                    K202
        d1,d9(2D) s1(2S)  3D   -1
        d1,d9(2D) s1(2S)  1D   3 
        d2,d8(3F) s1(2S)  4F   -2
        d2,d8(3P) s1(2S)  4P   -2
        d2,d8(1G) s1(2S)  2G   0 
        d2,d8(3F) s1(2S)  2F   4 
        d2,d8(1D) s1(2S)  2D   0 
        d2,d8(3P) s1(2S)  2P   4 
        d2,d8(1S) s1(2S)  2S   0 
        d3,d7(4F) s1(2S)  5F   -3
        d3,d7(4P) s1(2S)  5P   -3
        d3,d7(2H) s1(2S)  3H   -1
        d3,d7(2G) s1(2S)  3G   -1
        d3,d7(2H) s1(2S)  1H   3 
        d3,d7(2G) s1(2S)  1G   3 
        d3,d7(2F) s1(2S)  1F   3 
        d3,d7(2P) s1(2S)  1P   3 
        d4,d6(5D) s1(2S)  6D   -4
        d4,d6(3H) s1(2S)  4H   -2
        d4,d6(3G) s1(2S)  4G   -2
        d4,d6(1I) s1(2S)  2I   0 
        d4,d6(3H) s1(2S)  2H   4 
        d5(6S)    s1(2S)  7S   -5
        d5(4G)    s1(2S)  5G   -3
        d5(4F)    s1(2S)  5F   -3
        d5(4D)    s1(2S)  5D   -3
        d5(4P)    s1(2S)  5P   -3
        d5(6S)    s1(2S)  5S   7 
        d5(2I)    s1(2S)  3I   -1
        d5(2H)    s1(2S)  3H   -1
        d5(2S)    s1(2S)  3S   -1
        d5(2I)    s1(2S)  1I   3 
        d5(2H)    s1(2S)  1H   3 
        d5(2P)    s1(2S)  1P   3 
        d5(2S)    s1(2S)  1S   3 



        # f,s intershell coefficients for states with open shells of different symmetry
        f,s                    K303
        f1,f13(2F) s1(2S)  3F  -1
        f1,f13(2F) s1(2S)  1F  3 
        f2,f12(3H) s1(2S)  4H  -2
        f2,f12(3F) s1(2S)  4F  -2
        f2,f12(3P) s1(2S)  4P  -2
        f2,f12(1I) s1(2S)  2I  0 
        f2,f12(3H) s1(2S)  2H  4 
        f2,f12(1G) s1(2S)  2G  0 
        f2,f12(3F) s1(2S)  2F  4 
        f2,f12(1D) s1(2S)  2D  0 
        f2,f12(3P) s1(2S)  2P  4 
        f2,f12(1S) s1(2S)  2S  0 
        f3,f11(4I) s1(2S)  5I  -3
        f3,f11(4G) s1(2S)  5G  -3
        f3,f11(4F) s1(2S)  5F  -3
        f3,f11(4D) s1(2S)  5D  -3
        f3,f11(4S) s1(2S)  5S  -3
        f3,f11(2L) s1(2S)  3L  -1
        f3,f11(2K) s1(2S)  3K  -1
        f3,f11(2P) s1(2S)  3P  -1
        f3,f11(4S) s1(2S)  3S  5 
        f3,f11(2L) s1(2S)  1L  3 
        f3,f11(2K) s1(2S)  1K  3 
        f3,f11(2I) s1(2S)  1I  3 
        f3,f11(2P) s1(2S)  1P  3 
        f4,f10(5I) s1(2S)  6I  -4
        f4,f10(5G) s1(2S)  6G  -4
        f4,f10(5F) s1(2S)  6F  -4
        f4,f10(5D) s1(2S)  6D  -4
        f4,f10(5S) s1(2S)  6S  -4
        f4,f10(3M) s1(2S)  4M  -2
        f4,f10(3L) s1(2S)  4L  -2
        f4,f10(5S) s1(2S)  4S  6 
        f4,f10(1N) s1(2S)  2N  0 
        f4,f10(3M) s1(2S)  2M  4 
        f5,f9(6H)  s1(2S)  7H  -5
        f5,f9(6F)  s1(2S)  7F  -5
        f5,f9(6P)  s1(2S)  7P  -5
        f5,f9(4M)  s1(2S)  5M  -3
        f5,f9(4L)  s1(2S)  5L  -3
        f5,f9(4S)  s1(2S)  5S  -3
        f5,f9(2O)  s1(2S)  3O  -1
        f5,f9(2N)  s1(2S)  3N  -1
        f5,f9(4S)  s1(2S)  3S  5 
        f5,f9(2O)  s1(2S)  1O  3 
        f5,f9(2N)  s1(2S)  1N  3 
        f6,f8(7F)  s1(2S)  8F  -6
        f6,f8(5L)  s1(2S)  6L  -4
        f6,f8(5K)  s1(2S)  6K  -4
        f6,f8(5P)  s1(2S)  6P  -4
        f6,f8(5S)  s1(2S)  6S  -4
        f6,f8(3O)  s1(2S)  4O  -2
        f6,f8(3N)  s1(2S)  4N  -2
        f6,f8(5S)  s1(2S)  4S  6 
        f6,f8(1Q)  s1(2S)  2Q  0 
        f6,f8(3O)  s1(2S)  2O  4 
        f7(8S)     s1(2S)  9S  -7
        f7(6I)     s1(2S)  7I  -5
        f7(6H)     s1(2S)  7H  -5
        f7(6G)     s1(2S)  7G  -5
        f7(6F)     s1(2S)  7F  -5
        f7(6D)     s1(2S)  7D  -5
        f7(6P)     s1(2S)  7P  -5
        f7(8S)     s1(2S)  7S  9 
        f7(4N)     s1(2S)  5N  -3
        f7(4M)     s1(2S)  5M  -3
        f7(2Q)     s1(2S)  3Q  -1
        f7(2O)     s1(2S)  3O  -1
        f7(2Q)     s1(2S)  1Q  3 
        f7(2O)     s1(2S)  1O  3 



        # g,s intershell coefficients for states with open shells of different symmetry
        g,s                     K404
        g1,g17(2G) s1(2S)  3G   -1
        g1,g17(2G) s1(2S)  1G   3 
        g8,g10(9G) s1(2S)  10G  -8
        g9(10S)    s1(2S)  11S  -9
        g9(10S)    s1(2S)  9S   11


        # d,p(2,14) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d1p1;d9p5;2D,2P   3F   -2/7   -2     3/7  
        d1p1;d9p5;2D,2P   3D   1      8/5    -3/5 
        d1p1;d9p5;2D,2P   3P   -1     0      -3   
        d1p1;d9p5;2D,2P   1F   -2/7   14/5   27/35
        d1p1;d9p5;2D,2P   1D   1      -4/5   9/5  
        d1p1;d9p5;2D,2P   1P   -1     4/5    21/5 


        # d,p(3,13) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d1p2;d9p4;2D,3P   4F   2/7    -8/5   6/35    
        d1p2;d9p4;2D,3P   4D   -1     -8/5   -12/5   
        d1p2;d9p4;2D,3P   4P   1      12/5   -12/5   
        d1p2;d9p4;2D,1D   2G   -4/7   -8/5   36/35   
        d1p2;d9p4;2D,1D   2S   -2     2/5    -12/5   
        d2p1;d8p5;3F,2P   4G   -1/7   -14/5  18/35   
        d2p1;d8p5;3F,2P   4F   3/7    2/5    -54/35  
        d2p1;d8p5;3F,2P   4D   -12/35 28/25  -666/175
        d2p1;d8p5;3P,2P   4D   1/5    -28/25 -12/25  
        d2p1;d8p5;3P,2P   4P   -1     -8/5   -12/5   
        d2p1;d8p5;3P,2P   4S   2      16/5   -6/5    
        d2p1;d8p5;1G,2P   2H   -4/7   -8/5   36/35   
        d2p1;d8p5;3P,2P   2S   2      -2/5   12/5    



        # d,p(4,12) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d1p3;d9p3;2D,4S   5D   0      -6/5   -9/5     
        d1p3;d9p3;2D,2D   3G   0      -6/5   27/35    
        d1p3;d9p3;2D,2D   3S   0      14/5   -9/5     
        d1p3;d9p3;2D,2D   1G   0      -6/5   87/35    
        d1p3;d9p3;2D,2D   1S   0      -6/5   -9/5     
        d2p2;d8p4;3F,3P   5G   1/7    -16/5  -18/35   
        d2p2;d8p4;3F,3P   5F   -3/7   -8/5   -114/35  
        d2p2;d8p4;3F,3P   5D   12/35  32/25  -804/175 
        d2p2;d8p4;3P,3P   5D   -1/5   -62/25 -48/25   
        d2p2;d8p4;3P,3P   5P   1      2/5    -12/5    
        d2p2;d8p4;3P,3P   5S   -2     -16/5  -24/5    
        d2p2;d8p4;3P,3P   3S   -2     8/5    12/5     
        d2p2;d8p4;1G,1D   1I   -8/7   -16/5  72/35    
        d2p2;d8p4;1G,1D   1H   10/7   2/5    36/35    
        d3p1;d7p5;4F,2P   5G   1/7    -14/5  3/35     
        d3p1;d7p5;4F,2P   5F   -3/7   -6/5   -93/35   
        d3p1;d7p5;4F,2P   5D   12/35  42/25  -699/175 
        d3p1;d7p5;4P,2P   5D   -1/5   -52/25 -33/25   
        d3p1;d7p5;4P,2P   5P   1      4/5    -9/5     
        d3p1;d7p5;4P,2P   5S   -2     -14/5  -21/5    
        d3p1;d7p5;2H,2P   3I   -3/7   -12/5  39/35    
        d3p1;d7p5;2H,2P   1I   -3/7   0      15/7     
        d3p1;d7p5;2P,2P   1S   4/7    -14/5  33/35    



        # d,p(5,11) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d1p4;d9p2;2D,3P   4F   -2/7   -4/5   -6/5      
        d1p4;d9p2;2D,3P   4D   1      -4/5   -6/5      
        d1p4;d9p2;2D,3P   4P   -1     -4/5   -6/5      
        d1p4;d9p2;2D,1D   2G   4/7    -4/5   48/35     
        d1p4;d9p2;2D,1D   2S   2      16/5   -6/5      
        d2p3;d8p3;3F,4S   6F   0      -12/5  -18/5     
        d2p3;d8p3;3P,4S   6P   0      -12/5  -18/5     
        d2p3;d8p3;3F,2D   4H   0      -12/5  24/35     
        d2p3;d8p3;1G,2D   2I   0      -12/5  12/5      
        d3p2;d7p4;4F,3P   6G   -1/7   -22/5  -66/35    
        d3p2;d7p4;4F,3P   6F   3/7    -18/5  -138/35   
        d3p2;d7p4;4F,3P   6D   -12/35 -12/5  -1086/175 
        d3p2;d7p4;4P,3P   6D   1/5    -68/25 -72/25    
        d3p2;d7p4;4P,3P   6P   -1     -16/5  -24/5     
        d3p2;d7p4;4P,3P   6S   2      8/5    -18/5     
        d3p2;d7p4;2H,3P   4I   3/7    -12/5  24/35     
        d3p2;d7p4;2H,1D   2K   -6/7   -18/5  96/35     
        d4p1;d6p5;5D,2P   6F   2/7    -12/5  -36/35    
        d4p1;d6p5;5D,2P   6D   -1     -12/5  -18/5     
        d4p1;d6p5;5D,2P   6P   1      8/5    -18/5     
        d4p1;d6p5;3H,2P   4I   -1/7   -12/5  24/35     
        d4p1;d6p5;1I,2P   2K   -2/7   -2     12/7      



        # d,p(6,10) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d1p5;d9p1;2D,2P   3F   2/7    -2/5   -3/5   
        d1p5;d9p1;2D,2P   3D   -1     -2/5   -3/5   
        d1p5;d9p1;2D,2P   3P   1      -2/5   -3/5   
        d1p5;d9p1;2D,2P   1F   2/7    -2/5   159/35 
        d1p5;d9p1;2D,2P   1D   -1     -2/5   -3/5   
        d1p5;d9p1;2D,2P   1P   1      38/5   -3/5   
        d2p4;d8p2;3F,3P   5G   -1/7   -8/5   -12/5  
        d2p4;d8p2;3F,3P   5F   3/7    -8/5   -12/5  
        d2p4;d8p2;3F,3P   5D   -12/35 -8/5   -12/5  
        d2p4;d8p2;3P,3P   5D   1/5    -8/5   -12/5  
        d2p4;d8p2;3P,3P   5P   -1     -8/5   -12/5  
        d2p4;d8p2;3P,3P   5S   2      -8/5   -12/5  
        d2p4;d8p2;3P,3P   3S   2      28/5   12/5   
        d2p4;d8p2;1G,1D   1I   8/7    -8/5   96/35  
        d2p4;d8p2;1G,1D   1H   -10/7  -8/5   6/35   
        d3p3;d7p3;4F,4S   7F   0      -18/5  -27/5  
        d3p3;d7p3;4P,4S   7P   0      -18/5  -27/5  
        d3p3;d7p3;4P,2P   5S   0      -2     -3     
        d3p3;d7p3;2H,2D   3K   0      -18/5  81/35  
        d3p3;d7p3;2H,2D   1K   0      -6/5   177/35 
        d4p2;d6p4;5D,3P   7F   -2/7   -24/5  -132/35
        d4p2;d6p4;5D,3P   7D   1      -6/5   -24/5  
        d4p2;d6p4;5D,3P   7P   -1     -14/5  -36/5  
        d4p2;d6p4;3H,3P   5I   1/7    -18/5  -24/35 
        d4p2;d6p4;1I,1D   1L   -4/7   -4     24/7   
        d4p2;d6p4;1I,1D   1K   4/7    -4/5   48/35  
        d5p1;d5p5;6S,2P   7P   0      -2     -3     
        d5p1;d5p5;4G,2P   5H   0      -2     -3/7   
        d5p1;d5p5;4P,2P   5S   0      -2     -3     
        d5p1;d5p5;2I,2P   3K   0      -2     9/7    
        d5p1;d5p5;2I,2P   1K   0      -6/5   117/35 
        d5p1;d5p5;2P,2P   1S   0      -2     -3/7   



        # d,p(7,9) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d2p5;d8p1;3F,2P   4G   1/7    -4/5   -6/5 
        d2p5;d8p1;3F,2P   4F   -3/7   -4/5   -6/5 
        d2p5;d8p1;3F,2P   4D   12/35  -4/5   -6/5 
        d2p5;d8p1;3P,2P   4D   -1/5   -4/5   -6/5 
        d2p5;d8p1;3P,2P   4P   1      -4/5   -6/5 
        d2p5;d8p1;3P,2P   4S   -2     -4/5   -6/5 
        d2p5;d8p1;1G,2P   2H   4/7    -4/5   48/35
        d2p5;d8p1;3P,2P   2S   -2     -4/5   -6/5 
        d3p4;d7p2;4F,3P   6G   1/7    -12/5  -18/5
        d3p4;d7p2;4F,3P   6F   -3/7   -12/5  -18/5
        d3p4;d7p2;4F,3P   6D   12/35  -12/5  -18/5
        d3p4;d7p2;4P,3P   6D   -1/5   -12/5  -18/5
        d3p4;d7p2;4P,3P   6P   1      -12/5  -18/5
        d3p4;d7p2;4P,3P   6S   -2     -12/5  -18/5
        d3p4;d7p2;2H,3P   4I   -3/7   -12/5  -6/35 
        d3p4;d7p2;2H,1D   2K   6/7    -12/5  114/35
        d4p3;d6p3;5D,4S   8D   0      -24/5  -36/5 
        d4p3;d6p3;3H,4S   6H   0      -12/5  -18/5 
        d4p3;d6p3;5D,2D   6S   0      -6/5   -24/5 
        d4p3;d6p3;3H,2D   4K   0      -22/5  54/35 
        d4p3;d6p3;1I,2D   2L   0      -18/5  18/5  
        d5p2;d5p4;6S,3P   8P   0      -4     -6    
        d5p2;d5p4;4G,3P   6H   0      -4     -18/7 
        d5p2;d5p4;2I,3P   4K   0      -14/5  18/35 
        d5p2;d5p4;2I,1D   2L   0      -18/5  18/5  
        d6p1;d4p5;5D,2P   6F   -2/7   -8/5   -12/5 
        d6p1;d4p5;5D,2P   6D   1      -8/5   -12/5 
        d6p1;d4p5;5D,2P   6P   -1     -8/5   -12/5 
        d6p1;d4p5;3H,2P   4I   1/7    -8/5   6/35  
        d6p1;d4p5;1I,2P   2K   2/7    -8/5   66/35 




        # d,p(8) intershell coefficients for states with open shells of different symmetry
        d,p                    J212   K211   K213
        d3p5;d7p1;4F,2P   5G   -1/7   -6/5   -9/5  
        d3p5;d7p1;4F,2P   5F   3/7    -6/5   -9/5  
        d3p5;d7p1;4F,2P   5D   -12/35 -6/5   -9/5  
        d3p5;d7p1;4P,2P   5D   1/5    -6/5   -9/5  
        d3p5;d7p1;4P,2P   5P   -1     -6/5   -9/5  
        d3p5;d7p1;4P,2P   5S   2      -6/5   -9/5  
        d3p5;d7p1;2H,2P   3I   3/7    -6/5   27/35 
        d3p5;d7p1;2H,2P   1I   3/7    -6/5   21/5  
        d3p5;d7p1;2P,2P   1S   -4/7   -6/5   -9/5  
        d4p4;d6p2;5D,3P   7F   2/7    -16/5  -24/5 
        d4p4;d6p2;5D,3P   7D   -1     -16/5  -24/5 
        d4p4;d6p2;5D,3P   7P   1      -16/5  -24/5 
        d4p4;d6p2;3H,3P   5I   -1/7   -16/5  -48/35
        d4p4;d6p2;1I,1D   1L   4/7    -16/5  132/35
        d4p4;d6p2;1I,1D   1K   -4/7   -8/5   36/35 
        d5p3;6S,4S        9S   0      -6     -9    
        d5p3;4G,4S        7G   0      -18/5  -27/5 
        d5p3;4F,4S        7F   0      -18/5  -27/5 
        d5p3;6S,4S        7S   0      2/5    3/5   
        d5p3;2I,2D        3L   0      -22/5  99/35 
        d5p3;2I,2D        1L   0      -6/5   207/35



        # f,p intershell coefficients for states with open shells of different symmetry
        f,p                       J312   K312   K314
        f1(2F) p1(2P)    3G     -1/3   -15/7  10/21 
        f1(2F) p1(2P)    3F     1      9/7    -2/7 
        f1(2F) p1(2P)    3D     -4/5   9/35   -20/7
        f1(2F) p1(2P)    1G     -1/3   3      2/3
        f1(2F) p1(2P)    1F     1      -3/7   10/7 
        f1(2F) p1(2P)    1D     -4/5   3/5    4    
        f7(8S) p1(2P)    9P     0      -3     -4   
        f1(2F) p3(4S)    5F     0      -9/7   -12/7




        # f,d intershell coefficients for states with open shells of different symmetry
        f,d                     J322    J324   K321    K323    K325
        f1(2F) d1(2D)    3H     -10/21  -3/77  -81/35  -8/45   320/693 
        f1(2F) d1(2D)    3G     5/7     2/7    69/35   -58/45  20/63   
        f1(2F) d1(2D)    3F     11/21   -6/7   -27/35  -26/45  -20/63  
        f1(2F) d1(2D)    3D     -2/7    9/7    27/35   32/15   -40/21  
        f1(2F) d1(2D)    3P     -8/7    -6/7   3/35    -4/5    -30/7   
        f1(2F) d1(2D)    1H     -10/21  -3/77  99/35   32/45   340/693 
        f1(2F) d1(2D)    1G     5/7     2/7    -51/35  82/45   40/63   
        f1(2F) d1(2D)    1F     11/21   -6/7   9/7     10/9    80/63   
        f1(2F) d1(2D)    1D     -2/7    9/7    -9/35   -8/5    20/7    
        f1(2F) d1(2D)    1P     -8/7    -6/7   3/7     4/3     110/21  
        f7(8S) d1(2D)    9D     0       0      -9/5    -28/15  -10/3   
        f1(2F) d5(6S)    7F     0       0      -9/7    -4/3    -50/21  


        # s,s intershell and intrashell coefficients for states with two open shells of the same symmetry
        s,s           K000
        s1 s1    3S   -1



        # p,p intershell and intrashell coefficients for states with two open shells of the same symmetry
        p,p           J112   K110   K112
        p1 p1    3D   -1/5   -5/3   4/15  
        p1 p1    3S   -2     -5/3   -10/3 
        p1 p1    1P   1      -5/3   8/3   
        p2 p2    5D   -1/5   -8/3   -26/15
        p2 p2    5S   -2     -8/3   -16/3 
        p2 p2    1G   -4/5   -8/3   28/15 
        p3 p3    7S   0      -3     -6    
        p3 p3    3G   0      -3     -6/5  



        # d,d intershell and intrashell coefficients for states with two open shells of the same symmetry
        d,d           J222     J224     K220    K222      K224
        d1 d1    3G   -20/49   -1/49    -9/5    -26/49    116/245  
        d1 d1    3D   15/49    -36/49   -9/5    44/49     -234/245 
        d1 d1    3S   -10/7    -18/7    -9/5    -18/7     -162/35  
        d1 d1    1F   40/49    9/49     -9/5    94/49     216/245  
        d1 d1    1P   -5/7     12/7     -9/5    -8/7      138/35   
        d2 d2    5I   -5/49    -9/49    -16/5   -114/49   234/245  
        d2 d2    1L   -80/49   -4/49    -16/5   -24/49    484/245  
        d3 d3    7I   -5/49    -9/49    -21/5   -184/49   -396/245 
        d3 d3    3N   -45/49   -4/49    -21/5   -12/7     122/35   
        d4 d4    9G   -20/49   -1/49    -24/5   -236/49   -1774/245
        d4 d4    9D   -15/49   -36/49   -24/5   -166/49   -2124/245
        d4 d4    9S   -10/7    -18/7    -24/5   -48/7     -432/35  
        d4 d4    5N   -5/49    -16/49   -24/5   -18/7     68/35    
        d4 d4    1Q   -20/49   -20/49   -24/5   -58/7     1476/245 
        d5 d5   11S   0        0        -5      -50/7     -90/7    
        d5 d5    7L   0        0        -5      -150/49   -130/49  
        d5 d5    3Q   0        0        -5      -130/49   270/49   


        # f,f intershell and intrashell coefficients for states with two open shells of the same symmetry
        f,f           J332     J334     J336        K330    K332       K334        K336        
        f1 f1    3I   -5/9     -9/121   -25/14157   -13/7   -58/63     72/847      42550/99099 
        f1 f1    3G   2/3      -97/121  -50/363     -13/7   32/21      -1160/847   400/2541    
        f1 f1    3D   -19/45   9/11     -125/99     -13/7   -206/315   144/77      -1450/693   
        f1 f1    3S   -4/3     -18/11   -100/33     -13/7   -52/21     -234/77     -1300/231   
        f1 f1    1H   5/9      51/121   25/1089     -13/7   82/63      912/847     3650/7623   
        f1 f1    1F   2/9      3/11     50/99       -13/7   40/63      60/77       1000/693    
        f1 f1    1P   -1       -3/11    25/11       -13/7   -38/21     -24/77      1150/231    

        # keep this comment here
    '''

    lines = kj_table_data.splitlines()
    tables = []
    curtable = []
    for line in lines:
        loc = line.find('#')
        if loc!=-1:
            line = line[:loc]
        #end if
        line = line.strip()
        if len(line)>0:
            curtable.append(line)
        elif len(curtable)>0:
            tables.append(curtable)
            curtable = []
        #end if
    #end for
    table_lines = tables

    tables = obj()
    for tlines in table_lines:
        tokens = tlines[0].split()
        id_token = tokens[0]
        clabels = tokens[1:]
        coeff_labels = []
        for clabel in clabels:
            etype = clabel[0]
            ind   = tuple(array(tuple(clabel[1:]),dtype=int))
            coeff_labels.append((etype,ind))
        #end for
        table = obj()
        if ',' not in id_token:
            nl = 1
            key = id_token
        else:
            nl = 2
            key = tuple(id_token.split(','))
        #end if
        if key not in tables:
            tables[key] = obj()
        #end if
        table = tables[key]

        for line in tlines[1:]:
            tokens   = line.split()
            term_pair = None
            term_pairs = None
            if nl==1:
                id_token = tokens[0]
                term     = tokens[1]
                ctokens  = tokens[2:]
                if ',' in id_token:
                    lorbs = id_token.split(',')
                else:
                    lorbs = [id_token]
                #end if
            elif nl==2:
                if len(tokens)==len(coeff_labels)+3:
                    term_pairs = []
                    id_tokens = tokens[0:2]
                    term      = tokens[2]
                    ctokens   = tokens[3:]
                    lorb_factors = []
                    lorb_terms   = []
                    for id_token in id_tokens:
                        loc = id_token.find('(')
                        if loc!=-1:
                            id_token,iterm = id_token.replace('(',' ').replace(')',' ').split()
                        else:
                            iterm = None
                        #end if
                        if ',' in id_token:
                            ids = id_token.split(',')
                        else:
                            ids = [id_token]
                        #end if
                        lorb_fac = ids
                        lorb_factors.append(lorb_fac)
                        lorb_terms.append(len(ids)*[iterm])
                    #end for
                    lorbs = []
                    n1 = 0
                    for lf in lorb_factors[0]:
                        n2=0
                        for rf in lorb_factors[1]:
                            lorbs.append(lf+rf)
                            lorbs.append(rf+lf)
                            lt = lorb_terms[0][n1]
                            rt = lorb_terms[1][n2]
                            term_pairs.append((lt,rt))
                            term_pairs.append((rt,lt))
                            n1+=1
                        #end for
                        n2+=1
                    #end for
                    lorbs = list(set(lorbs))
                elif len(tokens)==len(coeff_labels)+2:
                    id_token = tokens[0]
                    term     = tokens[1]
                    ctokens  = tokens[2:]
                    id_tokens = id_token.split(';')
                    lorb_pairs = id_tokens[0:-1]
                    term_pair  = id_tokens[-1].split(',')
                    lorbs = []
                    for lorb_pair in lorb_pairs:
                        lf = lorb_pair[0:2]
                        rf = lorb_pair[2:4]
                        lorbs.append(lf+rf)
                        lorbs.append(rf+lf)
                    #end for
                    lorbs = list(set(lorbs))
                else:
                    KJTables.class_error('table misformatted')
                #end if
            #end if
            jcoeff   = []
            kcoeff   = []
            n=0
            for etype,ind in coeff_labels:
                ctext = ctokens[n]
                if '/' in ctext:
                    num,denom = ctext.split('/')
                else:
                    num = ctext
                    denom = 1
                #end if
                ecoeff = ind[0],ind[1],ind[2],int(num),int(denom)
                if etype=='J':
                    jcoeff.append(ecoeff)
                elif etype=='K':
                    kcoeff.append(ecoeff)
                else:
                    KJTables.class_error('ecoeff')
                #end if
                n+=1
            #end for
            n=0
            for lorb in lorbs:
                if term_pairs!=None:
                    term_pair = term_pairs[n]
                    if term_pair==(None,None):
                        term_pair = None
                    #end if
                #end if
                if lorb not in table:
                    table[lorb] = obj()
                #end if
                ltable = table[lorb]
                ecoeff = obj(term=term,jcoeff=jcoeff,kcoeff=kcoeff)
                if term_pair!=None:
                    ecoeff.term_pair = term_pair
                #end if
                ltable.append( ecoeff )
                n+=1
            #end for
        #end for
    #end for
    for key in list(tables.keys()):
        if isinstance(key,tuple):
            l1,l2 = key
            table_entry = tables[key]
            del tables[key]
            tables[l1+l2] = table_entry
            tables[l2+l1] = table_entry
        #end if
    #end for

    ordered_tables = tables
    tables = obj()
    for k1,otable in ordered_tables.iteritems():
        table = obj()
        for k2,oentries in otable.iteritems():
            entries = obj()
            for v in oentries:
                entries[v.term] = v
                if 'term_pair' in v:
                    tp = v.term_pair
                    term_triple = tp[0],tp[1],v.term
                    entries[term_triple] = v
                #end if
            #end for
            table[k2] = entries
        #end for
        tables[k1] = table
    #end for

    #print tables.dp.d8p1.keys()
    #print tables.ds.d8s1.keys()
    #exit()


    ne_filled = obj(s=2,p=6,d=10,f=14,g=18,h=22)
    l2sym = obj({
            0 :'s',
            1 :'p',
            2 :'d',
            3 :'f',
            4 :'g',
            5 :'h',
            6 :'i',
            7 :'k',
            8 :'l',
            9 :'m',
            10:'n',
            11:'o',
            12:'q',
            13:'r',
            14:'t',
            15:'u',
            16:'v',
            17:'w',
            18:'x',
            19:'y',
            20:'z',
            })
    sym2l = l2sym.inverse()

    @classmethod
    def get_kjcoeff(cls,shells,openshell_pairs=None,verbose=False,info=False):
        if isinstance(shells,str):
            shells = shells.split()
        #end if
        open_shells = []
        closed_shells = []
        terms = []
        for shell in shells:
            if '(' in shell:
                shell,term = shell.replace('(',' ').replace(')',' ').split()
            else:
                term = None
            #end if
            nshell = -1
            if shell[0].isdigit():
                nshell = int(shell[0])
                shell  = shell[1:]
            #end if
            shell = shell.lower()
            l = shell[0]
            ne = int(shell[1:])
            if ne<cls.ne_filled[l]:
                open_shells.append((l,shell,term,nshell))
            else:
                closed_shells.append((l,shell,term,nshell))
            #end if
        #end for
        nopen = len(open_shells)
        if nopen==0:
            jcoeff = []
            kcoeff = []
        elif nopen==1:
            l,shell,term,nshell = open_shells[0]
            ecoeff = cls.get_ecoeff(l,shell,term)
            if ecoeff!=None:
                jcoeff = sorted(ecoeff.jcoeff)
                kcoeff = sorted(ecoeff.kcoeff)
            else:
                jcoeff = None
                kcoeff = None
            #end if
        else:
            jcoeff = []
            kcoeff = []
            for l,shell,term,nshell in open_shells:
                ecoeff = cls.get_ecoeff(l,shell,term)
                #print shell,ecoeff
                if ecoeff!=None:
                    jcoeff.extend(ecoeff.jcoeff)
                    kcoeff.extend(ecoeff.kcoeff)
                else:
                    jcoeff = None
                    kcoeff = None
                    break
                #end if
            #end for
            if jcoeff!=None:
                if openshell_pairs!=None:
                    if isinstance(openshell_pairs,str):
                        openshell_pairs = openshell_pairs.split()
                    #end if
                    for opair in openshell_pairs:
                        if '(' in opair:
                            otokens = opair.replace('(',' ').replace(')',' ').split()
                            if len(otokens)==2:
                                opair,term = otokens
                            elif len(otokens)==5:
                                o1,t1,o2,t2,t3 = otokens
                                opair = o1+o2
                                term  = t1,t2,t3
                            else:
                                cls.class_error('openshell pairs')
                            #end if
                        else:
                            term = None
                        #end if
                        shellkey = opair
                        lkey = opair[0]+opair[2]
                        ecoeff = cls.get_ecoeff(lkey,shellkey,term)
                        #print lkey,shellkey,ecoeff
                        if ecoeff!=None:
                            jcoeff.extend(ecoeff.jcoeff)
                            kcoeff.extend(ecoeff.kcoeff)
                        else:
                            jcoeff = None
                            kcoeff = None
                            break
                        #end if
                    #end for
                else:
                    for i in xrange(len(open_shells)):
                        for j in xrange(i+1,len(open_shells)):
                            l1,shell1,term1,nshell1 = open_shells[i]
                            l2,shell2,term2,nshell2 = open_shells[j]
                            lkey = l1+l2
                            shellkey = shell1+shell2
                            term = None # ground state only
                            ecoeff = cls.get_ecoeff(lkey,shellkey,term)
                            #print lkey,shellkey,ecoeff
                            if ecoeff!=None:
                                jcoeff.extend(ecoeff.jcoeff)
                                kcoeff.extend(ecoeff.kcoeff)
                            else:
                                jcoeff = None
                                kcoeff = None
                                break
                            #end if
                        #end for
                    #end for
                #end if
            #end if
            if jcoeff!=None:
                jcoeff = sorted(jcoeff)
                kcoeff = sorted(kcoeff)
            #end if
        #end if
        if info==False:
            return jcoeff,kcoeff
        elif jcoeff is None:
            return None,None,None
        else:
            shells = closed_shells+open_shells
            lmax = 0
            for lsym,shell,term,nshell in shells:
                lmax = max(lmax,cls.sym2l[lsym])
            #end for
            nsym  = lmax+1
            ncsh  = zeros((nsym,),dtype=int)
            nosh  = zeros((nsym,),dtype=int)
            nccup = zeros((nsym,),dtype=int)
            nlbl  = []
            for lsym,shell,term,nshell in closed_shells:
                lval = cls.sym2l[lsym]
                ncsh[lval] += 1
                nlbl.append(nshell)
            #end for
            for lsym,shell,term,nshell in open_shells:
                lval = cls.sym2l[lsym]
                locc = int(shell[1:])
                nosh[lval]  += 1
                nccup[lval] += locc
                nlbl.append(nshell)
            #end for
            info = obj(
                nsym  = nsym,
                nj    = len(jcoeff),
                nk    = len(kcoeff),
                ncsh  = ncsh,
                nosh  = nosh,
                nccup = nccup,
                nlbl  = nlbl,
                )
            return jcoeff,kcoeff,info
        #end if
    #end def get_kjcoeff


    @classmethod
    def get_ecoeff(cls,lkey,skey,term):
        ecoeff = None
        if term is None:
            term = 0
            tabs = cls.ordered_tables
        else:
            tabs = cls.tables
        #end if
        if lkey in tabs:
            table = tabs[lkey]
            if skey in table:
                stable = table[skey]
                if term in stable:
                    ecoeff = stable[term]
                #end if
            #end if
        #end if
        return ecoeff
    #end def get_ecoeff

#end class KJTables


#jc,kc = KJTables.get_kjcoeff('5d2(3P)')
#jc,kc = KJTables.get_kjcoeff('2s2 2p4')
#jc,kc = KJTables.get_kjcoeff('4s1 3d8 4p1')
#jc,kc = KJTables.get_kjcoeff('4s1 3d8 4p1','p1s1(3P) d8s1(4F) d8p1(4D)')
#jc,kc = KJTables.get_kjcoeff('4s1 3d8(3F) 4p1',
#                             'p1s1(3P) d8(3F)s1(2S)4F d8(3F)p1(2P)4D')
#if jc!=None:
#    for j in jc:
#        print j
#    #end for
#    for k in kc:
#        print k
#    #end for
#else:
#    print 'mistake!'
##end if
#exit()


#Lsymbols_dict = {
#    0 :'S',
#    1 :'P',
#    2 :'D',
#    3 :'F',
#    4 :'G',
#    5 :'H',
#    6 :'I',
#    7 :'K',
#    8 :'L',
#    9 :'M',
#    10:'N',
#    11:'O',
#    12:'Q',
#    13:'R',
#    14:'T',
#    15:'U',
#    16:'V',
#    17:'W',
#    18:'X',
#    19:'Y',
#    20:'Z',
#    }
#Lsymbols = obj()
#Lsymbols.transfer_from(Lsymbols_dict)
#Lsymbols_inv = Lsymbols.inverse()
#
#
#
#import itertools
#def all_symbols(term):
#    print 'this function is wrong, do not use'
#    exit()
#    ne = int(term[1:])
#    if ne%2==0:
#        ne_o2 = float(ne)/2
#    else:
#        ne_o2 = float(ne)/2+1
#    #end if
#    lsym = term[0].upper()
#    l = Lsymbols_inv[lsym]
#    nmax = 2*l+1
#    all_ml = list(arange(nmax)-l)
#    all_ml.reverse()
#    all_ml = array(all_ml,dtype=int)
#    states = []
#    nu = min(ne,nmax)+1
#    while nu>ne_o2:
#        nu-=1
#        nd = ne-nu
#        S = nu-nd+1
#        uocc_base = nu*[1]+(nmax-nu)*[0]
#        uoccs = sorted(list(set(itertools.permutations(uocc_base))))
#        uoccs.reverse()
#        if nd>0:
#            docc_base = nd*[1]+(nmax-nd)*[0]
#            doccs = sorted(list(set(itertools.permutations(docc_base))))
#            doccs.reverse()
#        #end if
#        for uocc in uoccs:
#            L = 0
#            energy = 0
#            uml = all_ml[array(uocc)==1]
#            L+=uml.sum()
#            energy+=(l-uml).sum()
#            if nd==0:
#                states.append((nu,nd,L,S,energy,uml,[]))
#            else:
#                for docc in doccs:
#                    dml = all_ml[array(docc)==1]
#                    L+=dml.sum()
#                    energy+=(l-dml).sum()
#                    states.append((nu,nd,L,S,energy,uml,dml))
#                #end for
#            #end if
#        #end for
#    #end while
#
#    terms = obj()
#    for nu,nd,L,S,energy,uml,dml in states:
#        if L in Lsymbols:
#            Lsym = Lsymbols[L]
#            Tsym = str(S)+Lsym
#            if Tsym not in terms:
#                terms[Tsym] = obj()
#            #end if
#            terms[Tsym][energy] = (list(uml),list(dml))
#        #end if
#    #end for
#    
#    return terms
##end def all_symbols
#
##print
##print sorted(all_symbols('d3').keys())
##print all_symbols('s1')
##print all_symbols('p1')
##print all_symbols('d8')
#print 'all_symbols'
#exit()




class AtomscfInput(SimulationInput):

    keywords = set([
        'lim','aname','nflag1','nprint','nsym','nvar',
        'mxvar','nextra','iscale','iamp','ird','nzet',
        'nj','nk','ipnch','nlab','lcpu','ileg','nbas',
        'ncsh','nosh','nccup','nlbl','nbvar','ipqn',
        'zn','zscale','nxzet','zeta','ndiag','nxtrp',
        'mxxtrp','bias','dgath','scfat','expmn','dmnxp',
        'rmned','rmxed'
        ])

    input_keywords = set([
        'jcoeff','kcoeff','coeff','pseudo','ppformat' # not atomscf keywords
        ]) | keywords


    defaults = obj(
        nflag1 = 1,
        nprint = 0,
        nvar   = 0,
        iamp   = 0,
        ird    = 0,
        nzet   = 0,
        ipnch  = 0,
        nlab   = 1,
        ileg   = 0,
        zscale = 0.002,
        nextra = 0,
        )


    def __init__(self,filepath=None):
        self.lim      = 0
        self.sections = obj()
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def read_text(self,text,filepath=None):
        raw_lines = text.splitlines()
        lines = []
        for rline in raw_lines:
            loc = rline.find('/') # strip comments
            if loc==-1:
                lines.append(rline)
            else:
                lines.append(rline[:loc])
            #end if
        #end for
        self.lines = lines
        self.n=1
        self.cur_section = None
        self.lim = int(self.lines[0].strip())
        for ilim in range(self.lim):
            cs = obj()
            self.cur_section = cs
            self.read_val('aname',str)
            self.read_vals(['nflag1','nprint','nsym','nvar','mxvar','nextra','iscale','iamp','ird','nzet','nj','nk','ipnch','nlab','lcpu','ileg'],int)
            self.read_array('nbas',int)
            self.read_array('ncsh',int)
            self.read_array('nosh',int)
            self.read_array('nccup',int)
            if self.nonzero('nlab'):
                self.read_array('nlbl',int)
                self.set_count('nsht','nlbl')
            #end if
            if self.nonzero('nvar'):
                self.read_array('nbvar',int)
            #end if
            self.read_array('ipqn',int)
            self.set_count('nbast','ipqn')
            self.read_vals(['zn','zscale'],float)
            nzeta = cs.nbast
            if self.nonzero('nzet'):
                self.read_array('nxzet',float)
                nzeta = cs.nzet
            #end if
            self.read_array('zeta',float,nzeta)
            if self.nonzero('nj'):
                jcoeff = []
                for j in range(cs.nj):
                    jc = tuple(array(self.curtokens(),dtype=int))
                    jcoeff.append(jc)
                #end for
            #end if
            if self.nonzero('nk'):
                kcoeff = []
                for k in range(cs.nk):
                    kc = tuple(array(self.curtokens(),dtype=int))
                    kcoeff.append(kc)
                #end for
            #end if
            coeff = []
            nsh = cs.ncsh+cs.nosh
            for isym in range(cs.nsym):
                for ish in range(nsh[isym]):
                    coeff.append(array(self.curtokens(),dtype=float))
                #end for
            #end for
            cs.coeff = coeff
            #if self.nonzero('nj'):
            #    jcoeff = obj()
            #    for j in range(cs.nj):
            #        jcont = obj()
            #        self.read_vals(['lambda','mu','nu','num','nden'],int,jcont)
            #        jcoeff.append(jcont)
            #    #end for
            #    cs.jcoeff = jcoeff
            ##end if
            #if self.nonzero('nk'):
            #    kcoeff = obj()
            #    for k in range(cs.nk):
            #        kcont = obj()
            #        self.read_vals(['lambda','mu','nu','num','nden'],int,kcont)
            #        kcoeff.append(kcont)
            #    #end for
            #    cs.kcoeff = kcoeff
            ##end if
            #coeff = obj()
            #nsh = cs.ncsh+cs.nosh
            #for isym in range(cs.nsym):
            #    lcoeff = obj()
            #    for ish in range(nsh[isym]):
            #        self.read_array(ish,float,cs.nbas[isym],lcoeff)
            #    #end for
            #    coeff.append(lcoeff)
            ##end for
            cs.coeff = coeff
            if self.nonzero('ird'):
                self.read_vals(['ndiag','nxtrp','mxxtrp'],int)
                self.read_vals(['bias','dgath','scfat','expmn','dmnxp','rmned','rmxed'],int)
            #end if
            if self.nonzero('lcpu'):
                pptext = ''
                pptext += self.curline()+'\n' # cname
                pptext += self.curline()+'\n' # nzcor
                for l in range(cs.lcpu):
                    line = self.curline() # icpu
                    pptext += line+'\n'
                    icpu = int(line)
                    for i in range(icpu):
                        pptext += self.curline()+'\n' # channel primitives
                    #end for
                #end for
                pp = GaussianPP()
                pp.read_text(pptext,format='atomscf')
                Z = int(cs.zn)
                pp.Zval = Z-pp.Zcore
                pp.element = ptable.simple_elements[Z].symbol
                cs.pseudo = pp
            #end if
            self.sections.append(self.cur_section)
        #end for
        del self.lines
        del self.n
        del self.cur_section
    #end def read_text


    def write_text(self,filepath=None):
        self.text = ''
        self.cur_section = None
        self.text += '{0}\n'.format(self.lim)
        for ilim in range(self.lim):
            cs = self.sections[ilim]
            self.cur_section = cs
            self.write_val('aname',str)
            #self.write_vals(['nflag1','nprint','nsym','nvar','mxvar','nextra','iscale','iamp','ird','nzet','nj','nk','ipnch','nlab','lcpu','ileg'],int,cut=True)
            self.write_vals(['nflag1','nprint','nsym','nvar','mxvar','nextra','iscale','iamp','ird','nzet','nj','nk','ipnch','nlab','lcpu','ileg'],int,cut=False)
            self.write_array('nbas',int)
            self.write_array('ncsh',int)
            self.write_array('nosh',int)
            self.write_array('nccup',int)
            if self.nonzero('nlab'):
                self.write_array('nlbl',int)
            #end if
            if self.nonzero('nvar'):
                self.write_array('nbvar',int)
            #end if
            self.write_array('ipqn',int)
            self.write_vals(['zn','zscale'],float)
            if self.nonzero('nzet'):
                self.write_array('nxzet',float)
            #end if
            self.write_array('zeta',float)
            if self.nonzero('nj'):
                self.write_kjcoeff('jcoeff',cs.nj)
            #end if
            if self.nonzero('nk'):
                self.write_kjcoeff('kcoeff',cs.nk)
            #end if
            if 'coeff' not in cs:
                self.error('coeff is missing')
            #end if
            coeff = cs.coeff
            nsh = array(cs.ncsh)+array(cs.nosh)
            n = 0
            s = ''
            for isym in range(cs.nsym):
                for ish in range(nsh[isym]):
                    c = coeff[n]
                    for v in c:
                        s += str(v)+' '
                    #end for
                    s += '\n'
                    n += 1
                #end for
            #end for
            #print s
            self.text += s
            #if self.nonzero('nj'):
            #    if 'jcoeff' not in cs:
            #        self.error('jcoeff is missing')
            #    #end if
            #    jcoeff = cs.joeff
            #    for j in range(cs.nj):
            #        jcont = jcoeff[j]
            #        self.write_vals(['lambda','mu','nu','num','nden'],int,jcont)
            #    #end for
            ##end if
            #if self.nonzero('nk'):
            #    if 'kcoeff' not in cs:
            #        self.error('kcoeff is missing')
            #    #end if
            #    kcoeff = cs.kcoeff
            #    for k in range(cs.nk):
            #        kcont = kcoeff[k]
            #        self.write_vals(['lambda','mu','nu','num','nden'],int,kcont)
            #    #end for
            ##end if
            #coeff = cs.coeff
            #nsh = cs.ncsh+cs.nosh
            #for isym in range(cs.nsym):
            #    lcoeff = coeff[isym]
            #    for ish in range(nsh[isym]):
            #        self.write_array(ish,float,cs.nbas[isym],lcoeff)
            #    #end for
            ##end for
            if self.nonzero('ird'):
                self.write_vals(['ndiag','nxtrp','mxxtrp'],int)
                self.write_vals(['bias','dgath','scfat','expmn','dmnxp','rmned','rmxed'],int)
            #end if
            if self.nonzero('lcpu'):
                if 'pseudo' not in cs:
                    self.error('pseudo is missing')
                #end if
                pptext = cs.pseudo.write_text(format='atomscf')
                self.text += pptext
            #end if
            self.sections.append(self.cur_section)
        #end for
        text = self.text
        del self.text
        del self.cur_section
        return text
    #end def write_text


    def nonzero(self,name):
        return name in self.cur_section and self.cur_section[name]!=0
    #end def nonzero


    def set_count(self,count_name,array_name):
        self.cur_section[count_name] = len(self.cur_section[array_name])
    #end def set_count


    def curline(self):
        ls = self.lines[self.n].strip()
        self.n+=1
        return ls
    #end def curline


    def curtokens(self):
        tokens = self.lines[self.n].split()
        self.n+=1
        return tokens
    #end def curline
    

    def read_val(self,name,type):
        self.cur_section[name] = type(self.curline())
    #end def read_val


    def read_vals(self,names,type,container=None):
        if container is None:
            container = self.cur_section
        #end if
        tokens = self.curtokens()
        for n in range(len(tokens)):
            container[names[n]] = type(tokens[n])
        #end for
    #end def read_vals


    def read_array(self,name,type,nread=None,container=None):
        tokens = self.curtokens()
        if nread is not None:
            nmax = 100
            n = 0
            while len(tokens)<nread and n<nmax:
                tokens.extend(self.curtokens())
                n+=1
            #end while
            if len(tokens)!=nread:
                self.error('array read of {0} failed\nexpected {1} values\nfound {2} values'.format(name,nread,len(tokens)))
            #end if
            if n==nmax:
                self.error('array read of {0} failed, too many lines read'.format(name))
            #end if
        #end if
        if container is None:
            container = self.cur_section
        #end if
        container[name] = array(tokens,dtype=type)
    #end def read_array
    

    def write_val(self,name,type):
        if name not in self.cur_section:
            self.error('name or keyword is missing: {0}'.format(name))
        #end if
        s = str(self.cur_section[name])
        #print s
        self.text += s+'\n'
    #end def write_val


    def write_vals(self,names,type,container=None,cut=False):
        if container is None:
            container = self.cur_section
        #end if
        s = ''
        if not cut:
            for name in names:
                if name not in container:
                    self.error('name or keyword is missing: {0}'.format(name))
                #end if
                s += str(container[name])+' '
            #end for
        else:
            for name in names:
                if name in container:
                    s += str(container[name])+' '
                else:
                    break
                #end if
            #end for
        #end if
        #print s
        self.text += s+ '\n'
    #end def write_vals


    def write_array(self,name,type,nwrite=None,container=None):
        if container is None:
            container = self.cur_section
        #end if
        if name not in container:
            self.error('name or keyword is missing: {0}'.format(name))
        #end if
        s = ''
        for v in container[name]:
            s += str(v)+' '
        #end for
        #print s
        self.text += s + '\n'
    #end def write_array


    def write_kjcoeff(self,name,count):
        cs = self.cur_section
        if name not in cs:
            self.error(name+' is missing')
        #end if
        coeff = cs[name]
        s = ''
        for j in range(count):
            c = coeff[j]
            for v in c:
                s += str(v)+' '
            #end for
            s += '\n'
        #end for
        #print s
        self.text += s
    #end def write_kjcoeff


    def new_section(self,**kwargs):
        section = obj()
        if len(kwargs)>0:
            kw = obj(**kwargs)

            element  = kw.delete_optional('element')
            state    = kw.delete_optional('state')
            jcoeff   = kw.delete_optional('jcoeff')
            kcoeff   = kw.delete_optional('kcoeff')
            coeff    = kw.delete_optional('coeff')
            pseudo   = kw.delete_optional('pseudo')
            ppformat = kw.delete_optional('ppformat')
            basis    = kw.delete_optional('basis')
            guess    = kw.delete_optional('guess')

            use_none     = kw.delete_optional('use_none'    ,False)
            use_defaults = kw.delete_optional('use_defaults',True)
            use_state    = kw.delete_optional('use_state'   ,True)
            use_pseudo   = kw.delete_optional('use_pseudo'  ,True)
            use_basis    = kw.delete_optional('use_basis'   ,True)
            if use_none:
                use_defaults = False
                use_state    = False
                use_pseudo   = False
                use_basis    = False
            #end if

            if use_defaults:
                section.transfer_from(AtomscfInput.defaults)
            #end if

            if use_state and state is not None:
                jcoeff,kcoeff,info = KJTables.get_kjcoeff(state,info=True)
                section.transfer_from(info)
                section.jcoeff = jcoeff
                section.kcoeff = kcoeff
            #end if

            if element!=None:
                if element not in ptable:
                    self.error('element is not in the periodic table: {0}'.format(element))
                #end if
                section.zn = ptable[element].atomic_number
            #end if

            if pseudo!=None:
                if isinstance(pseudo,str):
                    if not os.path.exists(pseudo):
                        self.error('pseudopotential file does not exist: {0}'.format(pseudo))
                    #end if
                    if ppformat is None:
                        self.error('must provide ppformat to read pseudopotential file')
                    #end if
                    pseudo = GaussianPP(pseudo,ppformat)
                #end if
                if not isinstance(pseudo,GaussianPP):
                    self.error('pseudo must be a GaussianPP object\nyou provided: {0}'.format(pseudo.__class__.__name__))
                #end if
                if basis is None:
                    basis = pseudo.get_basis()
                #end if
                section.pseudo = pseudo
                section.lcpu = pseudo.lmax+1
                if not 'zn' in section and use_pseudo:
                    section.zn = ptable[pseudo.element].atomic_number
                #end if
            #end if

            if use_basis and basis!=None:
                if not isinstance(basis,GaussianBasisSet):
                    self.error('basis must be a GaussianBasisSet object\nyou provided: {0}'.format(basis.__class__.__name__))
                #end if
                if 'nsym' in kw:
                    nsym = kw.nsym
                elif 'nsym' in section:
                    nsym = section.nsym
                else:
                    self.error('must provide nsym to read basis')
                #end if
                basis = basis.copy()
                basis.uncontract()
                lbasis = basis.lbasis()
                lstr = 'spdfghiklmnoqrtuvwxyz'
                nbas = zeros((nsym,),dtype=int)
                ipqn = []
                zeta = []
                for l in range(nsym):
                    lsym = lstr[l]
                    lprims = lbasis[lsym]
                    nbas[l] = len(lprims)
                    ipqn.extend(len(lprims)*[l+1])
                    le = []
                    for i in xrange(len(lprims)):
                        le.append(lprims[i].terms[0].expon)
                    #end for
                    le.sort()
                    le.reverse()
                    zeta.extend(le)
                #end for
                section.nbas = nbas
                section.ipqn = ipqn
                section.zeta = zeta
            #end if



            # pull in all the user given keywords
            keys  = set(kw.keys())
            junk  = keys-AtomscfInput.keywords
            if len(junk)>0:
                self.error('received arguments that are not atomscf inputs: {0}'.format(sorted(junk)))
            #end if
            section.transfer_from(kw)

            if jcoeff!=None:
                section.jcoeff = jcoeff
            #end if
            if kcoeff!=None:
                section.kcoeff = kcoeff
            #end if
            if coeff!=None:
                section.coeff = coeff
            #end if

            if use_defaults and 'nvar' in section and section.nvar==0:
                for name in ('mxvar','nextra','iscale'):
                    if name not in kw:
                        section[name] = 0
                    #end if
                #end for
            #end if

            # apply a guess for the coefficients, if desired
            if guess!=None and 'coeff' not in section:
                nsym = section.nsym
                ncsh = section.ncsh
                nosh = section.nosh
                nbas = section.nbas
                nsh = ncsh+nosh
                coeff = []
                if guess=='random':
                    for isym in range(nsym):
                        nb = nbas[isym]
                        for ish in range(nsh[isym]):
                            coeff.append(random(nb)-.5)
                        #end for
                    #end for
                else:
                    self.error('invalid value for guess: {0}\nvalid options are: random'.format(guess))
                #end if
            #end if
            section.coeff = coeff

            # handle optimization
            if 'mxvar' in section and section.mxvar>0:
                if 'nvar' not in section or section.nvar==0:
                    nvar = 0
                    for c in coeff:
                        nvar += len(c)
                    #end for
                    section.nvar = nvar
                #end if
                if 'nvar' in section and 'nbvar' not in section:
                    section.nbvar = arange(section.nvar,dtype=int)+1
                #end if
                if 'iscale' not in section:
                    section.iscale = 0
                #end if
            #end if
        #end if
        self.lim += 1
        self.sections.append(section)
        return section
    #end def new_section

#end class AtomscfInput



def generate_atomscf_input(selector,**kwargs):
    if selector=='basic':
        return generate_basic_input(**kwargs)
    else:
        AtomscfInput.class_error('selection '+str(selector)+' has not been implemented for atomscf input generation')
    #end if
#end def generate_atomscf_input


def generate_basic_input(**kwargs):
    input = AtomscfInput()
    input.new_section(**kwargs)
    return input
#end def generate_basic_input



class AtomscfAnalyzer(SimulationAnalyzer):
    def __init__(self,filepath):
        None
    #end def __init__

    def analyze(self):
        None
    #end def analyze
#end class AtomscfAnalyzer



class Atomscf(Simulation):
    input_type    = AtomscfInput
    analyzer_type = AtomscfAnalyzer
    generic_identifier = 'atomscf'
    application   = 'atomscf'
    application_properties = set(['serial'])

    def check_sim_status(self):
        self.finished = True
    #end def check_sim_status

    def get_output_files(self):
        return []
    #end def get_output_files

    def app_command(self):
        return self.app_name+'<'+self.infile
    #end def app_command
#end class Atomscf



def generate_atomscf(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    if not 'input' in sim_args:
        input_type = inp_args.input_type
        del inp_args.input_type
        sim_args.input = generate_atomscf_input(input_type,**inp_args)
    #end if
    atomscf = Atomscf(**sim_args)

    return atomscf
#end def generate_atomscf
