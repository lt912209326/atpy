
from .cppelement cimport CppElement
from libc.math cimport pi,cos,sin,cosh,sinh,tan,atan,sqrt,exp
cdef cppclass CppQuadrupole(CppElement):
    __init__(double* parms)nogil:
        this.elem_type = 4
        this.l = parms[0]
        this.k1 = parms[2]
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
            
        

    inline void get_phase(double* parms0, double* parms)nogil:
        l = this.l
#         if this.k1 > 0:
#             k1 = sqrt(this.k1)
#             parms[8] = parms0[8] + (atan(parms0[2]*tan(k1*l)/k1 - parms0[1]) + atan(parms0[1]))/(2*pi)
#             parms[9] = parms0[9] + (atan((k1*parms0[3]*cosh(k1*l) + parms0[5]/k1*sinh(k1*l) - parms0[4]*exp(k1*l))*exp(k1*l))
#                                     -atan(k1*parms0[3] - parms0[4]))/(2*pi)
#         elif this.k1 < 0:
#             k1 = sqrt(-this.k1)
#             parms[8] = parms0[8] + (atan(( (k1*parms0[0]-2*parms0[1])*cosh(k1*l) + parms0[2]/k1*sinh(k1*l))*exp(k1*l)+parms0[1])
#                                     - atan(k1*parms0[0] - parms0[1]))/(2*pi)
#             parms[9] = parms0[9] + (atan(parms0[5]*tan(k1*l)/k1 - parms0[4]) + atan(parms0[4]))/(2*pi)
#         else:
#             parms[8] = parms0[8] + (atan(parms0[2]*l - parms0[1]) + atan(parms0[1]))/(2*pi)
#             parms[9] = parms0[9] + (atan(parms0[5]*l - parms0[4]) + atan(parms0[4]))/(2*pi)

        if this.k1 > 0:
            k1 = sqrt(this.k1)
            parms[9] = parms0[9] + (atan(parms0[2]*tan(k1*l)/k1 - parms0[1]) + atan(parms0[1]))/(2*pi)

            parms[10] = parms0[10] + (atan((k1*parms0[3]*cosh(k1*l) + parms0[5]/k1*sinh(k1*l) - parms0[4]*exp(k1*l))*exp(k1*l))
                                    -atan(k1*parms0[3] - parms0[4]))/(2*pi)
        elif this.k1 < 0:
            k1 = sqrt(-this.k1)
            parms[9] = parms0[9] + (atan(( (k1*parms0[0]-2*parms0[1])*cosh(k1*l) + parms0[2]/k1*sinh(k1*l))*exp(k1*l)+parms0[1])
                                    - atan(k1*parms0[0] - parms0[1]))/(2*pi)
            parms[10] = parms0[10] + (atan(parms0[5]*tan(k1*l)/k1 - parms0[4]) + atan(parms0[4]))/(2*pi)
        else:
            parms[9] = parms0[9] + (atan(parms0[2]*l - parms0[1]) + atan(parms0[1]))/(2*pi)
            parms[10] = parms0[10] + (atan(parms0[5]*l - parms0[4]) + atan(parms0[4]))/(2*pi)
    
    
    inline void update(double* parms)nogil:
        this.l = parms[0]
        this.k1= parms[2]
        this.update_matrixM()
        this.update_matrixT()