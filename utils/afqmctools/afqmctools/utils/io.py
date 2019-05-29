import numpy

def format_fixed_width_strings(strings):
    return ' '.join('{:>17}'.format(s) for s in strings)

def format_fixed_width_floats(floats):
    return ' '.join('{: .10e}'.format(f) for f in floats)

def to_qmcpack_complex(array):
    shape = array.shape
    return array.view(numpy.float64).reshape(shape+(2,))
