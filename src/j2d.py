""" Wrapper for Jo's code, to allow easy access from python """

import ctypes
import numpy as np

class _j2d(object):
    def __init__(self):
        self.dll = ctypes.CDLL("./j2d_dll.so")
        
        self._screenOutput = ctypes.c_bool.in_dll(self.dll, "screenOutput")
        self._p_ptr = ctypes.c_void_p.in_dll(self.dll, "p_ptr") 
        
    @property
    def output(self):
        return self._screenOutput.value
    
    @output.setter
    def output(self, value):
        if value:
            self._screenOutput.value = True
        else:
            self._screenOutput.value = False
    
    @property
    def N(self):
        return ctypes.c_int.in_dll(self.dll, "N").value
           
    @property
    def p(self):
        p_array_ptr = ctypes.POINTER(ctypes.c_longdouble)()
        size = dll.vector_as_array(_p_ptr, ctypes.byref(p_array_ptr))
        return np.ctypeslib.as_array(p_array_ptr, shape=size)
    
    @property
    def x(self):
        return p[:self.N]
    
    @property
    def y(self):
        return p[self.N:2*self.N]
    
    @property
    def alpha(self):
        return p[2*self.N]
    
    @alpha.setter
    def alpha(self, value):
        p[2*self.N] = value
    
    @property
    def delta(self):
        return p[(2*self.N)+1]
    
    @delta.setter
    def delta(self, value):
        p[(2*self.N)+1] = value
        
    @property
    def L(self):
        return p[(2*self.N)+2]
    
    @L.setter
    def L(self, value):
        p[(2*self.N)+2] = value
        
    def __getattr__(self, attr):
        return getattr(self.dll, attr)
        
    
j2d = _j2d()
        