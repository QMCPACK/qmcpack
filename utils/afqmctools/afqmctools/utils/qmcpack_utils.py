import h5py
from afqmctools.wavefunction.mol import read_qmcpack_wfn

def write_skeleton_input(qmc_in, hamil_file, wfn_file=None, series=0,
                         blocks=10000):
    if wfn_file is None:
        wfn_file = 'wfn.dat'
    with h5py.File(hamil_file, 'r') as  fh5:
        dims = fh5['Hamiltonian/dims'][:]
        nmo = dims[3]
        nalpha = dims[4]
        nbeta = dims[5]
        try:
            dset = fh5['Hamiltonian/KPFactorized']
            hamil_type = 'KPFactorized'
        except KeyError:
            hamil_type = 'Factorized'
        try:
            walker_type = fh5['Wavefunction/walker_type'][()]
            wfn_type = fh5['Wavefunction/type'][()]
        except KeyError:
            if wfn_file is not None:
                with open(wfn_file) as f:
                    header = [next(f) for i in range(5)]
                    uhf = len([s for s in header if 'UHF = 1' in s]) == 1
                    phmsd = len([s for s in header if 'TYPE = occ' in s]) == 1
                    walker_type = 'COLLINEAR' if uhf else 'CLOSED'
                    wfn_type = 'PHMSD' if phmsd else 'NOMSD'
                    if wfn_type == 'PHMSD':
                        walker_type = 'COLLINEAR'
    # TODO: Fix if GHF every comes back.

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
      <parameter name="filetype">ascii</parameter>
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
