
from ..constants cimport*
from .cppelement cimport CppElement
cdef cppclass CppMarker(CppElement):
    __init__()nogil:
        this.elem_type = MARKER
        this.l = 0.
        this.init_matrixM()
        this.update_matrixT()
    
    
    inline void get_phase(double* tws0, double* tws)nogil:
        tws[NUX] = tws0[NUX]
        tws[NUY] = tws0[NUY]