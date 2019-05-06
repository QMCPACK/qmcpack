def format_fixed_width_strings(strings):
    return ' '.join('{:>17}'.format(s) for s in strings)


def format_fixed_width_floats(floats):
    return ' '.join('{: .10e}'.format(f) for f in floats)
