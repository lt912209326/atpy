
from .cppelement cimport CppElement
from ..constants cimport*
from libc.math cimport sin,cos,sinh,cosh,tan,atan,exp,sqrt,pi
cdef cppclass CppSextupole(CppElement):
    __init__(double* kwd)nogil:
        this.elem_type = SEXTUPOLE
        this.l = kwd[L]
        this.k2 = kwd[K2]
#         this.rho = this.l/this.angle
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
        this.k2 = kwd[K2]
        this.update_matrixM()
        this.update_matrixT()