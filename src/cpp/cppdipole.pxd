from .cppelement cimport CppElement
from libc.math cimport sin,cos,sinh,cosh,tan,atan,exp,sqrt,pi
cdef cppclass CppDipole(CppElement):
    __init__(double* parms):
        this.elem_type = 2
        this.l = parms[0]
        this.angle = parms[1]
#         this.rho = this.l/this.angle
        this.init_matrixM()
        this.update_matrixM()
        this.update_matrixT()
    
    inline void update_matrixM()nogil:
        l = this.l
        angle = this.angle
        rho =this.l/angle
        
        this.M[0][1] = rho*sin(angle)
        this.M[0][5] = rho*(1-cos(angle))
        
        this.M[1][5] = 2*tan(angle/2)
        
        this.M[2][2] = 1-angle*tan(angle/2)
        this.M[2][3] = 1
        
        this.M[3][2] = -2*tan(angle/2)*(1-angle/2*tan(angle/2))/rho
        this.M[3][3] = 1-angle*tan(angle/2)
        
    
    inline void get_phase(double* parms0, double* parms)nogil:
        l = this.l
        angle = this.angle
        parms[9] = parms0[9] + (atan(parms0[2]*tan(angle)*l/angle - parms0[1]) + atan(parms0[1]))/(2*pi)
        parms[10] = parms0[10] + (atan(parms0[5]*l - parms0[4]) + atan(parms0[4]))/(2*pi)
#         parms[8] = parms0[8] + (atan(parms0[2]*tan(angle)*l/angle - parms0[1]) + atan(parms0[1]))/(2*pi)
#         parms[9] = parms0[9] + (atan(parms0[5]*l - parms0[4]) + atan(parms0[4]))/(2*pi)
    
    
    inline void update(double* parms)nogil:
        this.l = parms[0]
        this.angle=parms[1]
        this.update_matrixM()
        this.update_matrixT()