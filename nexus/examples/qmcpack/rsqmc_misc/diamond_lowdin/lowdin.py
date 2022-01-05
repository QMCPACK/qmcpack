#!/usr/bin/env python3
import sys
import numpy as np

def collectValuesFromAtomicProj(xmlfile):

    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()
    
    header = root.find('.//HEADER')
    
    # Find number of bands
    nBands = int(header.attrib['NUMBER_OF_BANDS'])
    # Find number of kpoints
    nKpoints = int(header.attrib['NUMBER_OF_K-POINTS'])
    # Find number of atomic wave functions
    nAtomicWFC = int(header.attrib['NUMBER_OF_ATOMIC_WFC'])
    # Find number of spin components
    nSpin = int(header.attrib['NUMBER_OF_SPIN_COMPONENTS'])

    kWeights = np.empty((nKpoints),dtype=float)

    atomicProjections = np.empty((nKpoints,nSpin,nAtomicWFC,nBands),dtype=complex)
    # Find atomic projections
    for k in range(nKpoints):
        kWeights[k] = float(root.findall('EIGENSTATES/K-POINT')[k].attrib['Weight'])
        for s in range(nSpin):
            for awfc in range(nAtomicWFC):
                if nSpin==1:
                    for b, text in enumerate(root.findall('EIGENSTATES/PROJS')[k][awfc].text.strip().splitlines()):
                        proj = float(text.split()[0])
                        proj = proj+complex(0,float(text.split()[1]))
                        # zeroth element below is for spin-type. In this case there is only one
                        atomicProjections[k][0][awfc][b]=proj
                    #end for
                else:
                    for b, text in enumerate(root.findall('EIGENSTATES/PROJS')[s*nKpoints+k][awfc].text.strip().splitlines()):
                        proj = float(text.split()[0])
                        proj = proj+complex(0,float(text.split()[1]))
                        atomicProjections[k][s][awfc][b]=proj
                    #end for
                    #for b, text in enumerate(root.find('EIGENSTATES/PROJS')[k][s][awfc].text.strip().splitlines()):
                    #    proj = float(text.split()[0])
                    #    proj = proj+complex(0,float(text.split()[1]))
                    #    atomicProjections[k][s][awfc][b]=proj
                    ##end for
                #end if
            #end for
        #end for
    #end for

    atomicOverlaps = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)

    # Find atomic overlaps
    for k in range(nKpoints):
        for s in range(nSpin):
            if nSpin==1:
                for o, text in enumerate(root.findall('OVERLAPS/OVPS')[k].text.strip().splitlines()):
                    ovlp = float(text.split()[0])
                    ovlp = ovlp+complex(0,float(text.split()[1]))
                    atomicOverlaps[k][0][o//nAtomicWFC][o%nAtomicWFC]=ovlp
                #end for
            else:
                for o, text in enumerate(root.findall('OVERLAPS/OVPS')[s*nKpoints+k].text.strip().splitlines()):
                    ovlp = float(text.split()[0])
                    ovlp = ovlp+complex(0,float(text.split()[1]))
                    atomicOverlaps[k][s][o//nAtomicWFC][o%nAtomicWFC]=ovlp
                #end for
            #end if
        #end for
    #end for

    invAtomicOverlaps = np.copy(atomicOverlaps)
    tmp = np.copy(atomicOverlaps)
    # Store inverse of atomic overlaps
    for k in range(nKpoints):
        for s in range(nSpin):
            invAtomicOverlaps[k][s] = np.linalg.inv(tmp[k][s])
        #end for
    #end for

    return nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps

#end def

def collectValuesFromXML(xmlfile):

    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()

    totmag = int(float(root.find('.//magnetization/total').text))
    nElec = int(float(root.find('.//nelec').text))
    nAtom = int(float(root.find('.//atomic_structure').attrib['nat']))

    return nAtom,nElec,int((nElec+totmag)/2),int((nElec-totmag)/2)

#end def

def matprint(m):
    for row in m:
        for element in row:
            print("%0.5f" % element),
        #end for
        print("\n")
    #end for
#end def

if __name__ == '__main__':

    from developer import ci
    from qmcpack_analyzer import QmcpackAnalyzer
    from uncertainties import ufloat,unumpy

    # Exit if atomic_proj.xml and outdir locations not given
    if(len(sys.argv)<5):
        print("Usage: lowdin.py <pw_prefix> <pw_outdir> <qmc_directory> <qmc_identifier> <spin>")
        quit()
    #end if

    pw_prefix = sys.argv[1]
    pw_outdir = sys.argv[2]

    qmc_directory = sys.argv[3]
    qmc_identifier = sys.argv[4]

    # spin (up=0,down=1)
    sp = int(sys.argv[5])

    if not sp in (0,1):
        print('Invalid spin specfied: {}'.format(sp))
        print('Must be either 0 (up) or 1 (down)')
        quit()
    #end if

    # Collect parameters from atomic_proj.xml.
    # Note: if atomic_proj.xml was not generated from projwfc.x or if the <OVERLAPS> element is not present in atomic_proj.xml then
    #       you should try re-running projwfc.x on a single core and single thread with "-ndiag 1" -- this can sometimes help
    nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps = collectValuesFromAtomicProj(pw_outdir+"/"+pw_prefix+".save/atomic_proj.xml")

    # Collect parameters from <prefix>.xml
    nAtom,nElec,nOccUp,nOccDown = collectValuesFromXML(pw_outdir+"/"+pw_prefix+".xml")

    print('\nNumber of up electrons: {}'.format(nOccUp))
    print('Number of down electrons: {}'.format(nOccDown))

    # Analyze QMC data
    qa = [] # qmcpack_analyzer instance
    nm = [] # number matrix
    for tn in range(nKpoints):
        qa_tmp = QmcpackAnalyzer('{}/{}.g{:03d}.twistnum_{}.in.xml'.format(qmc_directory,qmc_identifier,tn,tn),verbose=False)
        qa_tmp.analyze()
        qa.append(qa_tmp)

        # get the density matrix (called the number matrix here)
        nm_tmp = []

        if sp==0:
            nm_tmp.append(qa[tn].qmc[0].DensityMatrices.number_matrix.u.data)
        else:
            nm_tmp.append(qa[tn].qmc[0].DensityMatrices.number_matrix.d.data)
        #end if

        nm.append(nm_tmp)
    #end for

    nm = np.array(nm)

    # Obtain dimensions of number matrices

    nblocks,nstates,nstates = nm[0][0].shape

    # Store stats of number matrix corresponding to single determinant with no jastrow, projected
    # on MO basis

    from numerics import simstats

    m_mo,v_mo,e_mo,k_mo = simstats(nm,dim=2) # stats over blocks

    # Perform "unitary" transform on each block's number matrix individually
    # and store in nmqmcu (i.e., up component of number matrix prime)
    # After the transformation, number matrix has been transformed from
    # the MO basis to the AO basis

    s=sp

    nmqmc = np.empty((nKpoints,nSpin,nblocks,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        for b in range(nblocks):
            nmqmc[k][s][b] = kWeights[k]*np.matmul(atomicProjections[k][s][:,:],np.matmul(nm[k][0][b][:,:],np.conj(atomicProjections[k][s][:,:].T)))
        #end for
    #end for
    m_ao,v_ao,e_ao,k_ao = simstats(nmqmc,dim=2)
    m_mo_avg = np.sum(unumpy.uarray(m_mo.real,e_mo.real),axis=0)
    m_ao_avg = np.sum(unumpy.uarray(m_ao.real,e_ao.real),axis=0)

    # Obtain exact number matrix corresponding to single determinant with no jastrow, projected
    # on AO basis.

    exct_nmqmc = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        exct_nmqmc[k][s] = kWeights[k]*np.matmul(atomicProjections[k][s][:,:nOccUp],np.conj(atomicProjections[k][s][:,:nOccUp].T))
    #end for
    exavg = np.sum(exct_nmqmc,axis=0)


    # Print real part of mean of number matrix in MO basis
    print('nElec',nElec)

    print("\n     Total Charge of system (QMCPACK): " + str(np.trace(m_ao_avg[s])) +"\n")
    for a in range(nAtomicWFC):
        print("          charge on AO "+str(a)+" = "+str(m_ao_avg[sp][a][a]))
    #end for

    print("\n     Total Charge of system (QE): " + str(np.trace(exavg[s].real)) +"\n")
    for a in range(nAtomicWFC):
        print("          charge on AO "+str(a)+" = "+str(exavg[sp][a][a].real))
    #end for

    print()

#end if
