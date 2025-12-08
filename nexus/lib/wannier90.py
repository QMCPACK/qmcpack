 from wannier90_input import Wannier90Input


class Wannier90(Simulation):
    input_type = Wannier90Input
    # analyzer_type = Wannier90Analyzer
    generic_identifier = 'wannier90'
    infile_extension = '.win'
    application = 'wannier90'
    application_properties = set(['serial','omp','mpi'])