#!/usr/bin/env python3

from afqmctools.utils.qmcpack_utils import write_xml_input

options = {
    "Wavefunction": {
        "rediag": True, # Recompute ci coefficients
        "ndet": 10 # Truncate determinant expansion further
    },
    "execute": {
        "nWalkers": 10,
        "blocks": 100,
        "timestep": 0.01
    }
}

write_xml_input("afqmc.xml", "afqmc.h5", "afqmc.h5",
                options=options, rng_seed=7)
