
from numpy.lib.stride_tricks import as_strided as ast

def block_view(A, block= (3, 3)):
    ## Credit to http://stackoverflow.com/a/5078155/1828289
    """Provide a 2D block view to 2D array. No error checking made.
    Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A."""
    shape= (A.shape[0]// block[0], A.shape[1]// block[1])+ block
    strides= (block[0]* A.strides[0], block[1]* A.strides[1])+ A.strides
    return ast(A, shape= shape, strides= strides)


# conv functions

def convLST(arr):
    return (arr * 0.02) + (-273)

def convVIM(arr):
    arr[arr < 0] = 0
    return (arr * 0.0001)
