import h5py
import numpy as np

#     <determinantset> \
#       <multideterminant optimize=\"yes\" spo_0=\"spo-up\" spo_1=\"spo-dn\" spo_2=\"spo-ps\"> \
#         <detlist size=\"2\" type=\"DETS\" nc0=\"0\" nc1=\"0\" nc2=\"0\" ne0=\"2\" ne1=\"2\" ne2=\"1\" nstates=\"5\" cutoff=\"1e-20\"> \
#           <ci id=\"CIcoeff_0\" coeff=\"0.7071\" qchem_coeff=\"0.7071\" occ0=\"11000\" occ1=\"11000\" occ2=\"10000\"/> \
#           <ci id=\"CIcoeff_1\" coeff=\"-0.7071\" qchem_coeff=\"-0.7071\" occ0=\"10100\" occ1=\"11000\" occ2=\"00100\" /> \
#         </detlist> \
#       </multideterminant> \
#     </determinantset> \

with h5py.File("LiH.Multidet.h5", "w") as hdf5_file:
    # Multidet group
    md = hdf5_file.create_group("MultiDet")
    # Determinants
    dets_0 = np.array([[3],[5]], dtype=np.uint)
    dets_1 = np.array([[3],[3]], dtype=np.uint)
    dets_2 = np.array([[1],[4]], dtype=np.uint)
    
    detset_0 = md.create_dataset('CI_0', data=dets_0)
    detset_1 = md.create_dataset('CI_1', data=dets_1)
    detset_2 = md.create_dataset('CI_2', data=dets_2)
    # Coeff
    coeff = np.array([0.7071, -0.7071])
    
    coeffset = md.create_dataset('Coeff', data=coeff)
    # misc bits
    numdets = md.create_dataset('NbDet', data=np.array([2]))
    numbits = md.create_dataset('Nbits', data=np.array([1]))
    nexited = md.create_dataset('mexcotedstate', data=np.array([1]))
    nstate = md.create_dataset('nstate', data=np.array([5]))
    
