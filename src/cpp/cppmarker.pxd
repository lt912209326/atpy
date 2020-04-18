from .cppelement cimport CppElement
cdef cppclass CppMarker(CppElement):
    __init__():
        this.elem_type = 0
        this.l = 0.
        this.init_matrixM()
        this.update_matrixT()
    
    
    inline void get_phase(double* parms0, double* parms)nogil:
        parms[9] = parms0[9]
        parms[10] = parms0[10]
