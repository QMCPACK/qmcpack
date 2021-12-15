#! /usr/bin/env python
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
            print("%0.5f" % element),
        print("\n")

if __name__ == '__main__':

    from developer import ci
    from qmcpack_analyzer import QmcpackAnalyzer
    from uncertainties import ufloat,unumpy

    # Exit if atomic_proj.xml and outdir locations not given
    if(len(sys.argv)<2):
        print("Usage: spillage.py <prefix> <outdir>")
        quit()

    # Collect values from atomic_proj.xml

    nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps = collectValuesFromAtomicProj(sys.argv[2]+"/"+sys.argv[1]+".save/atomic_proj.xml")

    # Collect values from <prefix>.xml

    nAtom,nElec,nOccUp,nOccDown = collectValuesFromXML(sys.argv[2]+"/"+sys.argv[1]+".xml")

    # Analyze QMC data for both spin manifolds

    qa = []
    nm = []
    for tn in range(nKpoints):

        qa_tmp = []

        qa_tmp_u = QmcpackAnalyzer('runs/vmc_1rdm_noJ/vmc_1rdm_noJ.g00'+str(tn)+'.twistnum_'+str(tn)+'.in.xml',verbose=False)
        qa_tmp_u.analyze()
        qa_tmp.append(qa_tmp_u)

        qa_tmp_d = QmcpackAnalyzer('runs/vmc_1rdm_noJ_dwn/vmc_1rdm_noJ_dwn.g00'+str(tn)+'.twistnum_'+str(tn)+'.in.xml',verbose=False)
        qa_tmp_d.analyze()
        qa_tmp.append(qa_tmp_d)

        qa.append(qa_tmp)

        # get the density matrix (called the number matrix here)

        nm_tmp = []

        nm_tmp.append(qa[tn][0].qmc[0].DensityMatrices.number_matrix.u.data)
        nm_tmp.append(qa[tn][1].qmc[0].DensityMatrices.number_matrix.d.data)

        nm.append(nm_tmp)

    nm = np.array(nm)
    print 'Shape of nm : (num_k,num_s,num_b,num_ks,num_ks) = ' + str(nm.shape)

    # shape of nmu[k][s] : 200 x 8 x 8
    #   number of basis (bloch states) = 8
    #   number of samples (blocks) from qmcpack = 200

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

    nmqmc = np.empty((nKpoints,nSpin,nblocks,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        for s in range(nSpin):
            for b in range(nblocks):
                nmqmc[k][s][b] = kWeights[k]*np.matmul(atomicProjections[k][s],np.matmul(nm[k][s][b],np.conj(atomicProjections[k][s].T)))
    m_ao,v_ao,e_ao,k_ao = simstats(nmqmc,dim=2)

    # Obtain exact number matrix corresponding to single determinant with no jastrow, projected
    # on AO basis.

    exct_nmqmc = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        for s in range(nSpin):
            exct_nmqmc[k][s] = kWeights[k]*np.matmul(atomicProjections[k][s][:,:nOccUp],np.conj(atomicProjections[k][s][:,:nOccUp].T))

    m_mo_avg = np.sum(unumpy.uarray(m_mo.real,e_mo.real),axis=0)
    m_ao_avg = np.sum(unumpy.uarray(m_ao.real,e_ao.real),axis=0)
    exavg = np.sum(exct_nmqmc,axis=0)

    #import matplotlib.pyplot as plt

    #fig, ax = plt.subplots()
    #index = np.arange(nAtomicWFC)
    #bar_width = 0.35

    #fig = plt.figure(figsize=(16,6))

    #xlabel = []
    #for i in range(nAtomicWFC):
    #    xlabel.append(str(i))

    s=0
    #aobars = plt.bar(index+bar_width/2,np.diagonal(unumpy.nominal_values(m_ao_avg[s])),bar_width,color='b',label='QMCPACK',yerr=np.diagonal(unumpy.std_devs(m_ao_avg[s])),capsize=4)
    #exbars = plt.bar(index+bar_width*3/2,np.diagonal(exavg[s]),bar_width,color='g',label='QE')

    #plt.xlabel('AO')
    #plt.ylabel('Up Population')
    #plt.xticks(index + bar_width, xlabel)
    #plt.legend(loc='lower left')
    #plt.title('K-grid = 1x2x2')

    #plt.tight_layout()
    #plt.show()
    #fig.savefig('Dwn_1x2x2.pdf')

    #plt.bar(x_pos,np.diagonal(exavg[0]),align='center',alpha=0.5)
    #plt.xticks(x_pos,xlabel)
    #plt.ylabel('Charge')
    #plt.title('Populations')
    #plt.show()

    # Print real part of mean of number matrix in MO basis

    for sp in range(nSpin):

    #    print '\n\n\nSPIN '+str(sp)

    #    print("\n     (*********************************************************************************)")
    #    print("     (******************** Density matrix in MO basis from QMCPACK ********************)")
    #    print("     (*********************************************************************************)\n")
    #    #print(m_mo_avg[sp])
    #    print( "     Total Charge of system: " + str(np.trace(m_mo_avg[sp])) +"\n")

    #    for s in range(nstates):
    #        print("          charge on MO "+str(s)+" = "+str(m_mo_avg[sp][s][s]))


    #    print("\n     (*********************************************************************************)")
    #    print("     (******************** Density matrix in AO basis from QMCPACK ********************)")
    #    print("     (*********************************************************************************)\n")
    #    #print(m_ao_avg[sp])
        print( "     Total Charge of system: " + str(np.trace(m_ao_avg[sp])) +"\n")
    #    for a in range(nAtomicWFC):
    #        print("          charge on AO "+str(a)+" = "+str(m_ao_avg[sp][a][a]))
    #        #print("          S  charge on atom "+str(a)+" = "+str(m_ao[4*a][4*a].real))
    #        #print("          Pz charge on atom "+str(a)+" = "+str(m_ao[4*a+1][4*a+1].real))
    #        #print("          Px charge on atom "+str(a)+" = "+str(m_ao[4*a+2][4*a+2].real))
    #        #print("          Py charge on atom "+str(a)+" = "+str(m_ao[4*a+3][4*a+3].real))
    #

    #    print("\n     (*********************************************************************************)")
    #    print("     (******************* Density matrix in AO basis from PROJWFC.X *******************)")
    #    print("     (*********************************************************************************)\n")
    #    #matprint(exavg[sp].real)
        print( "     Total Charge of system (QE): " + str(np.trace(exavg[sp].real)) +"\n")
    #    for a in range(nAtomicWFC):
    #        print("          charge on AO "+str(a)+" = "+str(exavg[sp][a][a].real))
    #        #print("          S  charge on atom "+str(a)+" = "+str(exct_nmqmcu[4*a][4*a].real))
    #        #print("          Pz charge on atom "+str(a)+" = "+str(exct_nmqmcu[4*a+1][4*a+1].real))
    #        #print("          Px charge on atom "+str(a)+" = "+str(exct_nmqmcu[4*a+2][4*a+2].real))
    #        #print("          Py charge on atom "+str(a)+" = "+str(exct_nmqmcu[4*a+3][4*a+3].real))

    #   pop_s = []
    #   for b in range(nblocks):
    #       nmu_sample = nmu[b]
    #       # project onto AO basis using QE overlap matrix
    #       # get population for s state
    #       spop = 0.0
    #       pop_s.append(spop)
    #   #end for
    #
    #   # get stats for population
    #
    #   s_mean,s_var,s_err,s_kappa = simstats(pop_s)

