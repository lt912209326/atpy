
from .cppelement cimport CppElement
from libc.math cimport pi,atan
cdef cppclass CppDrift(CppElement):
    __init__(double* parms)nogil:
        this.elem_type = 1
        this.l = parms[0]
        this.init_matrixM()
        this.update_matrixM()
        this.update_matrixT()
    
    inline void update_matrixM()nogil:
        l = this.l
        this.M[0][1] = l
        this.M[2][3] = l
        this.M[4][5] = l
        
    
    inline void get_phase(double* parms0, double* parms)nogil:
        l = this.l
        parms[9] = parms0[9] + (atan(parms0[2]*l - parms0[1]) + atan(parms0[1]))/(2*pi)
        parms[10] = parms0[10] + (atan(parms0[5]*l - parms0[4]) + atan(parms0[4]))/(2*pi)
#         parms[8] = parms0[8] + (atan(parms0[2]*l - parms0[1]) + atan(parms0[1]))/(2*pi)
#         parms[9] = parms0[9] + (atan(parms0[5]*l - parms0[4]) + atan(parms0[4]))/(2*pi)
    
    
    inline void update(double* parms)nogil:
        this.l = parms[0]
        this.update_matrixM()
        this.update_matrixT()
    