import h5py
import os
import unittest
import xml.etree.ElementTree as et
from afqmctools.utils.qmcpack_utils import write_xml_input


class TestXMLWrite(unittest.TestCase):

    def test_write(self):
        options = {
            "execute": {
                "nWalkers": 10,
                "blocks": 1000,
                "timestep": 0.01,
                "Estimator": {
                    "back_propagation": {
                        "ortho": 1,
                        "naverages": 4,
                        "obs": {
                            "OneRDM": {
                                    "filename": "this.h5"
                                },
                            "TwoRDM": {
                                    "filename": "that.h5"
                                }
                            },
                        "block_size": 2,
                        "nsteps": 200
                    }
                }
            }
        }

        with h5py.File("afqmc.h5", 'w') as fh5:
            fh5['Wavefunction/NOMSD/dims'] = [37, 3, 4, 1]

        write_xml_input("afqmc.xml", "afqmc.h5", "afqmc.h5",
                        id_name="afqmc", options=options, rng_seed=7)

        xml_file = et.parse("afqmc.xml")
        nmo = xml_file.find("./AFQMCInfo/parameter[@name='NMO']").text
        self.assertEqual(nmo,"37")
        nalpha = xml_file.find("./AFQMCInfo/parameter[@name='NAEA']").text
        self.assertEqual(nalpha,"3")
        obs = xml_file.find("./execute/Estimator[@name='back_propagation']")
        fname = obs.find("./TwoRDM/parameter[@name='filename']").text
        self.assertEqual(fname, "that.h5")

    def tearDown(self):
        cwd = os.getcwd()
        files = ['afqmc.h5', 'afqmc.xml']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass

if __name__ == '__main__':
    unittest.main()
