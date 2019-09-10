from afqmctools.utils.qmcpack_utils import write_xml_input

options = {
    "execute": {
        "nWalkers": 10,
        "blocks": 100,
        "timestep": 0.01
    }
}

write_xml_input("afqmc.xml", "afqmc.h5", "afqmc.h5",
                options=options, rng_seed=7)
