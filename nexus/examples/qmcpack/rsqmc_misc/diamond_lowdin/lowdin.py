#!/usr/bin/env python3
import sys
import numpy as np

def collectValuesFromAtomicProj(xmlfile):

    import xml.etree.ElementTree as ET
    
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    
    # Find number of bands
    nBands = int(root.find('.//NUMBER_OF_BANDS').text)
    # Find number of kpoints
    nKpoints = int(root.find('.//NUMBER_OF_K-POINTS').text)
    # Find number of atomic wave functions
    nAtomicWFC = int(root.find('.//NUMBER_OF_ATOMIC_WFC').text)
    # Find number of spin components
    nSpin = int(root.find('.//NUMBER_OF_SPIN_COMPONENTS').text)

    kWeights = np.empty((nKpoints),dtype=float)
    # Find kpoint weights
    for ki, kw in enumerate(root.find('WEIGHT_OF_K-POINTS').text.strip().split()):
        kWeights[ki]=float(kw)

    atomicProjections = np.empty((nKpoints,nSpin,nAtomicWFC,nBands),dtype=complex)
    # Find atomic projections
    for k in range(nKpoints):
        for s in range(nSpin):
            for awfc in range(nAtomicWFC):
                if nSpin==1:
                    for b, text in enumerate(root.find('PROJECTIONS')[k][awfc].text.strip().splitlines()):
                        proj = float(text.split(',')[0])
                        proj = proj+complex(0,float(text.split(',')[1]))
                        atomicProjections[k][0][awfc][b]=proj
                else:
                    for b, text in enumerate(root.find('PROJECTIONS')[k][s][awfc].text.strip().splitlines()):
                        proj = float(text.split(',')[0])
                        proj = proj+complex(0,float(text.split(',')[1]))
                        atomicProjections[k][s][awfc][b]=proj
    
    atomicOverlaps = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)
    # Find atomic overlaps
    for k in range(nKpoints):
        for s in range(nSpin):
            for o, text in enumerate(root.find('OVERLAPS')[k][s].text.strip().splitlines()):
                ovlp = float(text.split(',')[0])
                ovlp = ovlp+complex(0,float(text.split(',')[1]))
                atomicOverlaps[k][s][o//nAtomicWFC][o%nAtomicWFC]=ovlp

    invAtomicOverlaps = np.copy(atomicOverlaps)
    tmp = np.copy(atomicOverlaps)
    # Store inverse of atomic overlaps
    for k in range(nKpoints):
        for s in range(nSpin):
            invAtomicOverlaps[k][s] = np.linalg.inv(tmp[k][s])

    return nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps

def collectValuesFromXML(xmlfile):

    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()

    totmag = int(float(root.find('.//magnetization/total').text))
    nElec = int(float(root.find('.//nelec').text))
    nAtom = int(float(root.find('.//atomic_structure').attrib['nat']))

    return nAtom,nElec,int((nElec+totmag)/2),int((nElec-totmag)/2)

def matprint(m):

    for row in m:
        for element in row:
            print("%0.8f" % element),
        print("\n")

if __name__ == '__main__':

    from developer import ci
    from qmcpack_analyzer import QmcpackAnalyzer
    import os
    
    qa = QmcpackAnalyzer('./runs/vmc_1rdm_noJ/vmc_1rdm_noJ.in.xml',verbose=True)
    qa.analyze()
   
    # get the density matrix (called the number matrix here)
    
    nm = qa.qmc[0].DensityMatrices.number_matrix
    
    # get the up and down data
    
    nmu = nm.u.data
    nmd = nm.d.data
    
    # shape of nmu: 200 x 8 x 8
    #   number of basis (bloch states) = 8
    #   number of samples (blocks) from qmcpack = 200

    # Exit if atomic_proj.xml location not given
    if(len(sys.argv)<2):
        print("Usage: spillage.py <pw.x prefix> <pw.x outdir>")
        quit()
    #end if
    pw_outdir = sys.argv[2]
    pw_prefix = sys.argv[1]

    # Collect values from atomic_proj.xml

    if not os.path.exists('{}/{}.save/atomic_proj.xml'.format(pw_outdir,pw_prefix)):
        print('{}/{}.save/atomic_proj.xml does not exist\nExiting...'.format(pw_outdir,pw_prefix))
        exit()
    #end if
    nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps = collectValuesFromAtomicProj(pw_outdir+"/"+pw_prefix+".save/atomic_proj.xml")

    # Collect values from <prefix>.xml

    nAtom,nElec,nOccUp,nOccDown = collectValuesFromXML(pw_outdir+"/"+pw_prefix+".xml") 

    # Obtain dimensions of number matrix

    nblocks,nstates,nstates = nmu.shape

    # Store stats of number matrix corresponding to single determinant with no jastrow, projected
    # on MO basis

    from numerics import simstats

    m_mo,v_mo,e_mo,k_mo = simstats(nmu,dim=0) # stats over blocks

    # Perform "unitary" transform on each block's number matrix individually
    # and store in nmqmcu (i.e., up component of number matrix prime)
    # After the transformation, number matrix has been transformed from
    # the MO basis to the AO basis

    nmqmcu = np.empty((nblocks,nstates,nstates),dtype=complex)
    for b in range(nblocks):
        nmqmcu[b] = np.matmul(atomicProjections[0][0],np.matmul(nmu[b],np.conj(atomicProjections[0][0].T)))
    m_ao,v_ao,e_ao,k_ao = simstats(nmqmcu,dim=0)

    # Obtain exact number matrix corresponding to single determinant with no jastrow, projected
    # on AO basis. 

    exct_nmqmcu = np.matmul(atomicProjections[0][0][:,:nOccUp],np.conj(atomicProjections[0][0][:,:nOccUp].T))

    # Print real part of mean of number matrix in MO basis
    print("\n     (*********************************************************************************)")
    print("     (******************** Density matrix in MO basis from QMCPACK ********************)")
    print("     (*********************************************************************************)\n")
    matprint(2*m_mo.real)
    print( "     Total Charge of system: " + str(np.trace(2*m_mo.real)) +"\n")
    for s in range(nstates):
        print("          charge on MO "+str(s)+" = "+str(2*m_mo[s][s].real))
   
    # Print real part of mean of number matrix in MO basis
    print("\n     (*********************************************************************************)")
    print("     (******************** Density matrix in AO basis from QMCPACK ********************)")
    print("     (*********************************************************************************)\n")
    matprint(2*m_ao.real)
    print( "     Total Charge of system: " + str(np.trace(2*m_ao.real)) +"\n")
    for a in range(nAtom):
        print("          S  charge on atom "+str(a)+" = "+str(2*m_ao[4*a][4*a].real))
        print("          Pz charge on atom "+str(a)+" = "+str(2*m_ao[4*a+1][4*a+1].real))
        print("          Px charge on atom "+str(a)+" = "+str(2*m_ao[4*a+2][4*a+2].real))
        print("          Py charge on atom "+str(a)+" = "+str(2*m_ao[4*a+3][4*a+3].real))

    # Print real part of mean of number matrix in MO basis
    print("\n     (*********************************************************************************)")
    print("     (******************* Density matrix in AO basis from PROJWFC.X *******************)")
    print("     (*********************************************************************************)\n")
    matprint(2*exct_nmqmcu.real)
    print( "     Total Charge of system: " + str(np.trace(2*exct_nmqmcu.real)) +"\n")
    for a in range(nAtom):
        print("          S  charge on atom "+str(a)+" = "+str(2*exct_nmqmcu[4*a][4*a].real))
        print("          Pz charge on atom "+str(a)+" = "+str(2*exct_nmqmcu[4*a+1][4*a+1].real))
        print("          Px charge on atom "+str(a)+" = "+str(2*exct_nmqmcu[4*a+2][4*a+2].real))
        print("          Py charge on atom "+str(a)+" = "+str(2*exct_nmqmcu[4*a+3][4*a+3].real))

    #  pop_s = []
    #  for b in range(nblocks):
    #      nmu_sample = nmu[b]
    #      # project onto AO basis using QE overlap matrix
    #      # get population for s state
    #      spop = 0.0
    #      pop_s.append(spop)
    #  #end for
    #  
    #  # get stats for population
    #  
    #  s_mean,s_var,s_err,s_kappa = simstats(pop_s)
    
