##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  differences.py                                                    #
#    Support for detailed comparisons (full diffs) of collections    #
#    of files.  Useful for rapid integration of branch files into    #
#    main trunk of Nexus, QMCPACK, etc.                              #
#                                                                    #
#  Content summary:                                                  #
#    Differences                                                     #
#      Class to find and summarize differences of a file collection. #
#                                                                    #                                        
#====================================================================#


#! /usr/bin/env python

import os
import sys
from numpy import array
from difflib import Differ
from generic import obj
from developer import DevBase 


class Differences(DevBase):
    def __init__(self,old,new,old_files,filemap=None):
        if filemap is None:
            filemap = dict()
        #end if
        if not os.path.exists(old):
            self.error('old path '+old+' does not exist')
        #end if
        if not os.path.exists(new):
            self.error('new path '+new+' does not exist')
        #end if
        files = old_files
        if isinstance(files,str):
            files = files.split()
        #end if
        remove = []
        extensions = set('h cpp'.split())
        for file in files:
            tokens = file.split('.')
            if len(tokens)==1 or not tokens[-1] in extensions:
                remove.append(file)
            #end if
        #end for
        for file in remove:
            files.remove(file)
        #end for
        for file in files:
            if not file in filemap:
                filemap[file] = file
            #end if
        #end for
        absent = []
        paths = [old,new]
        for path in paths:
            if path==old:
                filelist = filemap.keys()
            elif path==new:
                filelist = filemap.values()
            #end if
            for file in filelist:
                filepath = os.path.join(path,file)
                if not os.path.exists(filepath):
                    absent.append(filepath+'\n')
                #end if
            #end for
        #end for
        if len(absent)>0:
            self.log('\nFiles not found:')
            sys.stdout.writelines(absent)
            self.error('some files could not be found, see list above.')
        #end if
        diffs = obj()
        d = Differ()
        for oldfile,newfile in filemap.iteritems():
            oldc = open(os.path.join(old,oldfile),'r').read().splitlines()
            newc = open(os.path.join(new,newfile),'r').read().splitlines()
            diffs[(oldfile,newfile)] = list(d.compare(oldc,newc))
        #end for
        files.sort()
        self.set(
            old       = old,
            new       = new,
            oldfiles  = files,
            filemap   = filemap,
            diffs     = diffs
            )
    #end def __init__


    def summarize(self,order=None):
        oldfiles = self.oldfiles
        filemap = self.filemap
        diffs = self.diffs
        if order is None:
            ofile_list = oldfiles
        else:
            ofile_list = order
        #end if
        print
        print 'Difference summary:'
        print '==================='
        print '  old_path =',self.old
        print '  new_path =',self.new
        print '  differences: '
        scores = []
        slines = []
        ofiles = []
        for oldfile in ofile_list:
            newfile = filemap[oldfile]
            difflist = diffs[(oldfile,newfile)]
            np,nm,nq = 0,0,0
            nlines = 0
            for line in difflist:
                if len(line)>0:
                    nlines+=1
                    c = line[0]
                    if c=='+':
                        np+=1
                    elif c=='-':
                        nm+=1
                    elif c=='?':
                        nq+=1
                    #end if
                #end if
            #end for
            nt = np+nm+nq
            scores.append(nt)
            slines.append('    {0:<30}   {1:>4}   #(+,-,?) = ({2:>3},{3:>3},{4:>3})  {5:>5}'.format(oldfile,nt,np,nm,nq,nlines))
            ofiles.append(oldfile)
        #end for
        scores = array(scores)
        slines = array(slines)
        ofiles = array(ofiles)
        if order is None:
            order = scores.argsort()
            slines = slines[order]
            ofiles = ofiles[order]
        #end if
        for sline in slines:
            print sline
        #end for
        return ofiles
    #end def summarize


    def write(self,order=None):
        oldfiles = self.oldfiles
        filemap = self.filemap
        diffs = self.diffs
        if order is None:
            ofile_list = oldfiles
        else:
            ofile_list = order
        #end if
        print
        print
        print 'Difference lists:'
        print '================='
        print '  old_path =',self.old
        print '  new_path =',self.new
        print '  differences: '
        fpad = 4*' '
        lpad = 6*' '
        lpad_orig = lpad
        changes = set('+-?')
        for oldfile in ofile_list:
            newfile = filemap[oldfile]
            difflist = diffs[(oldfile,newfile)]
            header = oldfile+' -> '+newfile
            print
            print
            print '`ch '+header
            print fpad+len(header)*'-'
            for line in difflist:
                if len(line)>0 and line[0] in changes:
                    lpad = '  `ch '
                else:
                    lpad = lpad_orig
                #end if
                print lpad+line
            #end for
        #end for
    #end def write
#end class Differences



if __name__=='__main__':
    dp = Differences(
        old = './r5784_pristine',
        new = './r5906_pristine',
        old_files = '''
    AsymmetricDistanceTableData.h  get_source                QMCDriverFactory.cpp
    BareKineticEnergy.h            get_source~               QMCDriver.h
    CloneManager.cpp               HamiltonianFactory.cpp    QMCHamiltonianBase.cpp
    CloneManager.h                 hdf_dataproxy.h           QMCHamiltonianBase.h
    CMakeLists                     HDFNumericAttrib.h        QMCHamiltonian.cpp
    CoulombPBCAATemp.cpp           LocalECPotential.cpp      QMCHamiltonian.h
    CoulombPBCAATemp.h             LocalECPotential.h        QMCMain.cpp
    CoulombPBCABTemp.cpp           NonLocalECPComponent.cpp  QMCMain.h
    CoulombPBCABTemp.h             NonLocalECPComponent.h    QMCUpdateBase.cpp
    CoulombPotential.h             NonLocalECPotential.cpp   QMCUpdateBase.h
    DistanceTableData.h            NonLocalECPotential.h     ScalarEstimatorBase.h
    DMCOMP.cpp                     OhmmsArray.h              scalar_traits.h
    DMCUpdatePbyPFast.cpp          ParticleSet.cpp           SymmetricDistanceTableData.h
    EstimatorManager.cpp           ParticleSet.h             VMCSingleOMP.cpp
    EstimatorManager.h             QMCDriver.cpp             VMCUpdatePbyP.cpp
    ''',
        filemap = {'CoulombPBCAATemp.h'  :'CoulombPBCAA.h',
                   'CoulombPBCAATemp.cpp':'CoulombPBCAA.cpp',
                   'CoulombPBCABTemp.h'  :'CoulombPBCAB.h',
                   'CoulombPBCABTemp.cpp':'CoulombPBCAB.cpp'}
        )


    dmo = Differences(
        old = './r5784_pristine',
        new = './r5784_trace',
        old_files = '''
    AsymmetricDistanceTableData.h  get_source                QMCDriverFactory.cpp
    BareKineticEnergy.h            get_source~               QMCDriver.h
    CloneManager.cpp               HamiltonianFactory.cpp    QMCHamiltonianBase.cpp
    CloneManager.h                 hdf_dataproxy.h           QMCHamiltonianBase.h
    CMakeLists                     HDFNumericAttrib.h        QMCHamiltonian.cpp
    CoulombPBCAATemp.cpp           LocalECPotential.cpp      QMCHamiltonian.h
    CoulombPBCAATemp.h             LocalECPotential.h        QMCMain.cpp
    CoulombPBCABTemp.cpp           NonLocalECPComponent.cpp  QMCMain.h
    CoulombPBCABTemp.h             NonLocalECPComponent.h    QMCUpdateBase.cpp
    CoulombPotential.h             NonLocalECPotential.cpp   QMCUpdateBase.h
    DistanceTableData.h            NonLocalECPotential.h     ScalarEstimatorBase.h
    DMCOMP.cpp                     OhmmsArray.h              scalar_traits.h
    DMCUpdatePbyPFast.cpp          ParticleSet.cpp           SymmetricDistanceTableData.h
    EstimatorManager.cpp           ParticleSet.h             VMCSingleOMP.cpp
    EstimatorManager.h             QMCDriver.cpp             VMCUpdatePbyP.cpp
    '''
        )






    order = dp.summarize()
    #dp.write()

    dmo.summarize(order)

    dmo.write(order)

    dp.write(order)
#end if
