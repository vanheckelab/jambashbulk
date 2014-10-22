import os

import numpy as np

import ctypes
from ctypes import Structure, c_bool, c_int, c_double, c_longdouble, cdll, byref

if os.name == "nt":
    c_LDBL = c_double
    c_LDBL_p = ctypes.POINTER(c_double)
    np_LDBL = np.double
else:
    c_LDBL = c_longdouble
    c_LDBL_p = ctypes.POINTER(c_longdouble)
    np_LDBL = np.longdouble

sz_LDBL = np.dtype(np_LDBL).itemsize
    
class PACKINGPARAMS(Structure):
    _fields_ = [('P', c_LDBL),
                ('phi', c_LDBL),
                ('Z', c_LDBL),
                ('Ncorrected', c_int),
                ('sxx', c_LDBL),
                ('sxy', c_LDBL),
                ('syy', c_LDBL),
                ('U', c_LDBL),
                ('H', c_LDBL),
                ('gg', c_LDBL)]

PACKINGPARAMS_p = ctypes.POINTER(PACKINGPARAMS)

dll_path = os.path.join(
        os.path.split(os.path.abspath(__file__))[0],
        '..',
        'bin',
        'jamBashbulk.so')
dll = cdll.LoadLibrary(dll_path)

__version__ = ctypes.c_char_p.in_dll(dll, "FILE_HEADER_C").value

class DLL(object):
    def __init__(self):
        raise Exception("Class used as namespace, cannot initialize")

    pack_loaded = False
    
    # void relax_packing(bool alphaFree, bool deltaFree, bool LFree) {
    relax_packing = dll.relax_packing
    relax_packing.argtypes = [c_bool, c_bool, c_bool]

    # void import_packing(int _N, LDBL _P0,
    #                     LDBL *_x, LDBL *_y, LDBL *_r,
    #                     LDBL _alpha, LDBL _delta, LDBL _L)
    import_packing = dll.import_packing
    import_packing.argtypes = [c_int, c_LDBL,
                                 c_LDBL_p, c_LDBL_p, c_LDBL_p,
                                 c_LDBL, c_LDBL, c_LDBL]

    # void export_packing(int _N, LDBL * _P0,
    #                     LDBL *_x, LDBL *_y, LDBL *_r,
    #                     LDBL *_alpha, LDBL *_delta, LDBL *_L)
    export_packing = dll.export_packing
    export_packing.argtypes = [c_int, c_LDBL_p,
                               c_LDBL_p, c_LDBL_p, c_LDBL_p,
                               c_LDBL_p, c_LDBL_p, c_LDBL_p]

    # int get_packing_size()
    get_packing_size = dll.get_packing_size
    get_packing_size.argtypes = []
    get_packing_size.restype = c_int

    # void get_packing_data(packingparams *out)
    get_packing_data = dll.get_packing_data
    get_packing_data.argtypes = [PACKINGPARAMS_p]

def relax_packing(alphaFree=True, deltaFree=True, LFree=True):
    DLL.relax_packing(alphaFree, deltaFree, LFree)

def get_packing_size():
    return DLL.get_packing_size()
    
def convertLvectors(L1, L2):
    L = np.sqrt(L1[0] * L2[1])
    alpha = L2[0] / L
    delta = np.sqrt(L2[1] / L1[0]) - 1.0;

    return {'alpha': alpha, 'delta': delta, 'L': L}

def import_packing(pack):
    """ Imports a pack into the simulation code.
    dict format:
    { "L1": periodic box vector 1, (aligned with x)
      "L2": periodic box vector 2,
      "N": numer of particles,
      "P0": target pressure,
      "particles": array or dict with keys=["x", "y", "r"] and optional keys ["x_err", "y_err"]
                   for particle positions.
    }
    """
    N = pack["N"]
    P0 = pack["P0"]
    
    # calculate alpha, delta and L from L1 and L2
    res = convertLvectors(pack["L1"], pack["L2"])
    alpha = res["alpha"]
    delta = res["delta"]
    L = res["L"]
    
    x = np_LDBL(pack["particles"]["x"]) + np_LDBL(pack["particles"]["x_err"])
    y = np_LDBL(pack["particles"]["y"]) + np_LDBL(pack["particles"]["y_err"])
    r = np_LDBL(pack["particles"]["r"]) + 0 # force copy to make sure strides == 8

    assert(x.strides == (sz_LDBL,))
    assert(y.strides == (sz_LDBL,))
    assert(r.strides == (sz_LDBL,))
    
    x_ptr = x.ctypes.data_as(ctypes.POINTER(c_LDBL))
    y_ptr = y.ctypes.data_as(ctypes.POINTER(c_LDBL))
    r_ptr = r.ctypes.data_as(ctypes.POINTER(c_LDBL))
    
    DLL.import_packing(N, P0, x_ptr, y_ptr, r_ptr, alpha, delta, L)
    DLL.pack_loaded = True

def export_packing():
    """ Returns a packing dict from the loaded packing. Format:
    {'L': 26.21776889062624,
 'L1': array([ 25.77606417,   0.        ]),
 'L2': array([ -0.12999204,  26.66704277]),
 'N': 128,
 'P': 0.009995008664793113,
 'P0': 0.01000000000000001,
 'particles': array([(9.349917925369324, 0.0, 12.224353792603848, 0.0, 1.000000000000001),
       ...], 
      dtype=[('x', '<f8'), ('x_err', '<f8'), ('y', '<f8'), ('y_err', '<f8'), ('r', '<f8')])}
    """
    if not DLL.pack_loaded:
        raise Exception("No packing loaded in simulation!")

    # First, we get the packing size and create arrays to store particle positions
    N = get_packing_size()
    
    x = np.tile(np_LDBL(np.nan), N)
    y = np.tile(np_LDBL(np.nan), N)
    r = np.tile(np_LDBL(np.nan), N)
    
    assert(x.strides == (sz_LDBL,))
    assert(y.strides == (sz_LDBL,))
    assert(r.strides == (sz_LDBL,))  

    x_ptr = x.ctypes.data_as(ctypes.POINTER(c_LDBL))
    y_ptr = y.ctypes.data_as(ctypes.POINTER(c_LDBL))
    r_ptr = r.ctypes.data_as(ctypes.POINTER(c_LDBL))
    
    # Create storage for byref parameters
    P0 = c_LDBL()
    alpha = c_LDBL()
    delta = c_LDBL()
    L = c_LDBL()

    DLL.export_packing(N, byref(P0), x_ptr, y_ptr, r_ptr, byref(alpha), byref(delta), byref(L))
    
    P0 = P0.value
    alpha = alpha.value
    delta = delta.value
    L = L.value
    
    # now we just need to get current pressure P using get_packing_data
    P = get_packing_data()["P"]
    
    # reshape particle positions in x + x_err, y + y_err format
    x_major = np.float64(x); x_minor = np.float64(x-np_LDBL(x_major))
    y_major = np.float64(y); y_minor = np.float64(y-np_LDBL(y_major))
    r = np.float64(r)

    particles = np.array(zip(x_major, x_minor, y_major, y_minor, r), dtype=[('x', np.float64), ('x_err', np.float64),
                                                                            ('y', np.float64), ('y_err', np.float64),
                                                                            ('r', np.float64)])
                                                                            
    # and calculate L1 and L2 from alpha, delta and L
    
    lxx = L / (1.0 + delta);
    lxy = L * 0.0;
    lyx = L * alpha;
    lyy = L * (1.0 + delta);
    
    return {'P0': P0,
            'L': L,
            'N': N,
            'L1': np.array([lxx, lxy]),
            'L2': np.array([lyx, lyy]),
            'P': P,
            'particles': particles}
    
def get_packing_data():
    """ Retrieve packing data from an imported packing""" 
    if not DLL.pack_loaded:
        raise Exception("No packing loaded in simulation!")
    
    packingdata = PACKINGPARAMS()
    DLL.get_packing_data(byref(packingdata))
    return dict((field[0], getattr(packingdata, field[0])) for field in packingdata._fields_)

if __name__ == '__main__':
    import sys
    this_file_dir = os.path.split(__file__)[0]
    sys.path.append(os.path.join(this_file_dir, "..", "..", "phd-library"))
    from packing_tools.parser import parser
    packing = parser.read_packings(os.path.join(
        this_file_dir, "..", "..", "phd-library", "packing_tools", "parser", "N16~P1e-3~9000.txt"
    )).next()

    print "%s, using:" % os.path.split(__file__)[1]
    print __version__
    import_packing(packing)
    print export_packing()
    converge_packing
    print export_packing()
