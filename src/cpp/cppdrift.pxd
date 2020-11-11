
from ..cython.constants cimport*
from .cppelement cimport CppElement
from libc.math cimport pi,atan
cdef cppclass CppDrift(CppElement):
    __init__(double* kwd)nogil:
        this.elem_type = DRIFT
        this.l = kwd[L]
        this.init_matrixM()
        this.update_matrixM()
        this.update_matrixT()
    
    inline void update_matrixM()nogil:
        l = this.l
        this.M[0][1] = l
        this.M[2][3] = l
        this.M[4][5] = l
        
    
    inline void get_phase(double* tws0, double* tws)nogil:
        l = this.l
        tws[NUX] = tws0[NUX] + (atan(tws0[GAMMAX]*l - tws0[ALPHAX]) + atan(tws0[ALPHAX]))/(2*pi)
        tws[NUY] = tws0[NUY] + (atan(tws0[GAMMAY]*l - tws0[ALPHAY]) + atan(tws0[ALPHAY]))/(2*pi)
    
    
    inline void update(double* kwd)nogil:
        this.l = kwd[L]
        this.update_matrixM()
        this.update_matrixT()
    