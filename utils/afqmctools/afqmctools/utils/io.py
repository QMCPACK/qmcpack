import numpy

def format_fixed_width_strings(strings):
    return ' '.join('{:>17}'.format(s) for s in strings)

def format_fixed_width_floats(floats):
    return ' '.join('{: .10e}'.format(f) for f in floats)

def to_qmcpack_complex(array):
    shape = array.shape
    return array.view(numpy.float64).reshape(shape+(2,))

def from_qmcpack_complex(data, shape=None):
    if shape is not None:
        return data.view(numpy.complex128).ravel().reshape(shape)
    else:
        return data.view(numpy.complex128).ravel()

def add_dataset(fh5, name, value):
    try:
        fh5[name] = value
    except:
        del fh5[name]
        fh5[name] = value

def add_group(fh5, name):
    try:
        group = fh5.create_group(name)
    except:
        del fh5[name]
        group = fh5.create_group(name)
    return group
