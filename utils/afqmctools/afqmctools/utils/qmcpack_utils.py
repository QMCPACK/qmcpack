import h5py
from afqmctools.wavefunction.mol import read_qmcpack_wfn

def write_skeleton_input(qmc_in, hamil_file, wfn_file=None, series=0,
                         blocks=10000, nelec=None):
    walker_types = ['NONE', 'CLOSED', 'COLLINEAR', 'NONCOLLINEAR']
    if wfn_file is None:
        wfn_file = hamil_file
    with h5py.File(hamil_file, 'r') as  fh5:
        dims = fh5['Hamiltonian/dims'][:]
        nmo = dims[3]
        nalpha = dims[4]
        nbeta = dims[5]
        if nelec is not None:
            nalpha, nbeta = nelec
        try:
            dset = fh5['Hamiltonian/KPFactorized']
            hamil_type = 'KPFactorized'
        except KeyError:
            hamil_type = 'Factorized'
    with h5py.File(wfn_file, 'r') as  fh5:
        try:
            dims = fh5['Wavefunction/PHMSD/dims'][:]
            wfn_type = 'PHMSD'
        except KeyError:
            dims = fh5['Wavefunction/NOMSD/dims'][:]
            wfn_type = 'NOMSD'
        walker_type = walker_types[dims[3]]

    xml_string = """<?xml version="1.0"?>
<simulation method="afqmc">
    <project id="qmc" series="{:d}"/>

    <AFQMCInfo name="info0">
        <parameter name="NMO">{:d}</parameter>
        <parameter name="NAEA">{:d}</parameter>
        <parameter name="NAEB">{:d}</parameter>
    </AFQMCInfo>

    <Hamiltonian name="ham0" type="{:s}" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">{:s}</parameter>
    </Hamiltonian>

    <Wavefunction name="wfn0" type="{:s}" info="info0">
      <parameter name="filetype">hdf5</parameter>
      <parameter name="filename">{:s}</parameter>
      <parameter name="cutoff">1e-8</parameter>
    </Wavefunction>

    <WalkerSet name="wset0" type="shared">
      <parameter name="walker_type">{:s}</parameter>
    </WalkerSet>

    <Propagator name="prop0" info="info0">
      <parameter name="hybrid">yes</parameter>
    </Propagator>

    <execute wset="wset0" ham="ham0" wfn="wfn0" prop="prop0" info="info0">
      <parameter name="ncores">1</parameter>
      <parameter name="timestep">0.005</parameter>
      <parameter name="blocks">{:d}</parameter>
      <parameter name="steps">10</parameter>
      <parameter name="nWalkers">10</parameter>
      <Estimator name="energy">
          <parameter name="print_components">true</parameter>
      </Estimator>
   </execute>
</simulation>""".format(series, nmo, nalpha, nbeta, hamil_type,
                        hamil_file, wfn_type, wfn_file,
                        walker_type, blocks)
    with open(qmc_in, 'w') as f:
        f.write(xml_string)
