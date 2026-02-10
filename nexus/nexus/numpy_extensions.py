import numpy as np


def reshape_inplace(arr: np.ndarray, *new_shape) -> None:
    if len(new_shape) == 1 and hasattr(new_shape[0], "__iter__"):
        new_shape = new_shape[0]
    if hasattr(arr, '_set_shape'):
        arr._set_shape(new_shape)
    else:
        arr.shape = new_shape
#end def reshape_inplace
