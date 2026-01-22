import numpy as np

def reshape_array(arr: np.ndarray, *new_shape) -> None:
    """Reshape a numpy array in-place without changing the size.

    Parameters
    ----------
        arr : NDArray
            Array to be reshaped in-place.
        new_shape : int or tuple of int
            Shape of the reshaped array.
    """

    if len(new_shape)==1 and isinstance(new_shape[0], (tuple, list)):
        new_shape = new_shape[0]

    if np.prod(new_shape) != arr.size:
        raise ValueError(
            "Reshaping an array must not change the size.\n"
            "Array size before reshape: {0}\n"
            "Array size after reshape : {1}\n"
            "Array shape before reshape: {2}\n"
            "Array shape after reshape : {3}".format(
                arr.size, np.prod(new_shape), arr.shape, new_shape
            )
        )
    #end if

    prev_size     = arr.size
    prev_shape    = arr.shape
    prev_mem_addr = arr.ctypes.data

    arr.resize(new_shape, refcheck=False)

    if arr.size != prev_size:
        raise RuntimeError(
            "Reshaping an array must not change the size.\n"
            "Array size before reshape: {0}\n"
            "Array size after reshape : {1}\n"
            "Array shape before reshape: {2}\n"
            "Array shape after reshape : {3}".format(
                prev_size, arr.size, prev_shape, arr.shape
            )
        )
    #end if

    if arr.ctypes.data != prev_mem_addr:
        raise RuntimeError(
            "Additional array data buffer created during attempted reshape.\n"
            "Reshape is not performed in-place as required."
        )
    #end if
#end def reshape_array
