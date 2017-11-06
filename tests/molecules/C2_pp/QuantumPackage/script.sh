#!/bin/bash

#On cooley, one needs to source quantum_package
source /soft/applications/quantum_package/quantum_package.rc

#Create the input file for Quantum Package(QP) using the BFD ECP and the corresponding Basis Set for CC-pvTz 
qp_create_ezfio_from_xyz C2.xyz -b vtz-bfd -p bfd

#Forcing Quantum Package to create every variable in the Input. (Avoid using VIM to edit input)
qp_edit -c C2.ezfio

#Sets QP reader to use the C2.ezfio directory as input
ezfio set_file C2.ezfio

# Save the AO orbitals (save time for restart)
sed -i s/"None"/"Write"/ C2.ezfio/integrals_bielec/disk_access_ao_integrals

# Save the MO orbitals (save time for restart)
sed -i s/"None"/"Write"/ C2.ezfio/integrals_bielec/disk_access_mo_integrals

# Set the max determinants to 6M instead of 1M
sed -i s/"             1000000"/"             6000000"/ C2.ezfio/determinants/n_det_max

###Run Hartree Fock in parallel with 10 nodes (1 Master + 9 slaves)
#Master needs one node only While the slave nodes will attach to this job
mpirun -n 1 qp_run SCF C2.ezfio &> C2.HF.out &
sleep 10
mpirun -n 9 qp_run -slave qp_ao_ints C2.ezfio &> C2.HF-slave.out 
wait

#run the 4 index transformation (MO transformation) Serial run
qp_run four_idx_transform C2.ezfio 


###Run CIPSI in parallel with 10 nodes (1 Master + 9 slaves)
#Master needs one node only While the slave nodes will attach to this job
mpirun -n 1 qp_run fci_zmq C2.ezfio &> C2.FCI.out &
sleep 10
mpirun -n 9 qp_run -slave selection_davidson_slave C2.ezfio &> C2.FCI-slave.out 
wait

#Truncate the determinants to keep only weights of 10^-3 and higher. 
ezfio set QMC ci_threshold 1e-3 
mpirun -n 1 qp_run truncate_wf_spin C2.ezfio &> C2.Trunc-1e-3.out


qp_run save_for_qmcpack C2.ezfio &> C2-Dump-1e-3.out

convert4qmc -QP C2-Dump-1e-3.out -add3BodyJ 

