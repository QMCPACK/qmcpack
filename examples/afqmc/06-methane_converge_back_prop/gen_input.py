#!/usr/bin/env python3

from afqmctools.utils.qmcpack_utils import write_xml_input

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
                    "OneRDM": {}
                    },
                "block_size": 2,
                "nsteps": 200
            }
        }
    }
}

write_xml_input("afqmc.xml", "afqmc.h5", "afqmc.h5",
                options=options, rng_seed=7)
