from ..cython.constants cimport*
from .cppelement cimport CppElement
from libc.math cimport pi,cos,sin,cosh,sinh,tan,atan,sqrt,exp
cdef cppclass CppQuadrupole(CppElement):
    __init__(double* kwd)nogil:
        this.elem_type = QUADRUPOLE
        this.l = kwd[L]
        this.k1 = kwd[K1]
        this.init_matrixM()
        this.update_matrixM()
        this.update_matrixT()
    
        
    inline void update_matrixM()nogil:
        cdef int i
        l = this.l
        if this.k1 > 0:
            k1 = sqrt(this.k1)
            i = 1
        elif this.k1 < 0:
            k1 = sqrt(-this.k1)
            i = -1
        else:
            i = 0
        
        if i == 0:
            this.init_matrixM()
            this.M[0][1] = l
            this.M[2][3] = l
            this.M[4][5] = l
        elif i==1 or i==-1:
            this.M[1-i][1-i] = cos(k1*l)
            this.M[1-i][2-i] = sin(k1*l)/k1
            this.M[2-i][1-i] = -k1*sin(k1*l)
            this.M[2-i][2-i] = cos(k1*l)
            
            this.M[1+i][1+i] = cosh(k1*l)
            this.M[1+i][2+i] = sinh(k1*l)/k1
            this.M[2+i][1+i] = k1*sinh(k1*l)
            this.M[2+i][2+i] = cosh(k1*l)
            
    

    inline void get_phase(double* tws0, double* tws)nogil:
        l = this.l

        if this.k1 > 0:
            k1 = sqrt(this.k1)
            tws[NUX] = tws0[NUX] + (atan(tws0[GAMMAX]*tan(k1*l)/k1 - tws0[ALPHAX]) + atan(tws0[ALPHAX]))/(2*pi)

            tws[NUY] = tws0[NUY] + (atan((k1*tws0[BETAY]*cosh(k1*l) + tws0[GAMMAY]/k1*sinh(k1*l) - tws0[ALPHAY]*exp(k1*l))*exp(k1*l))
                                    -atan(k1*tws0[BETAY] - tws0[ALPHAY]))/(2*pi)
        elif this.k1 < 0:
            k1 = sqrt(-this.k1)
            tws[NUX] = tws0[NUX] + (atan(( (k1*tws0[BETAX]-2*tws0[ALPHAX])*cosh(k1*l) + tws0[GAMMAX]/k1*sinh(k1*l))*exp(k1*l)+tws0[ALPHAX])
                                    - atan(k1*tws0[BETAX] - tws0[ALPHAX]))/(2*pi)
            tws[NUY] = tws0[NUY] + (atan(tws0[GAMMAY]*tan(k1*l)/k1 - tws0[ALPHAY]) + atan(tws0[ALPHAY]))/(2*pi)
        else:
            tws[NUX] = tws0[NUX] + (atan(tws0[GAMMAX]*l - tws0[ALPHAX]) + atan(tws0[ALPHAX]))/(2*pi)
            tws[NUY] = tws0[NUY] + (atan(tws0[GAMMAY]*l - tws0[ALPHAY]) + atan(tws0[ALPHAY]))/(2*pi)
    
    
    inline void update(double* kwd)nogil:
        this.l = kwd[L]
        this.k1= kwd[K1]
        this.update_matrixM()
        this.update_matrixT()