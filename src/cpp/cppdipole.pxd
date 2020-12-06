
from .cppelement cimport CppElement
from ..cython.constants cimport*
from libc.math cimport sin,cos,sinh,cosh,tan,atan,exp,sqrt,pi
cdef cppclass CppDipole(CppElement):
    __init__(double* kwd)nogil:
        this.elem_type = DIPOLE
        this.l = kwd[L]
        this.angle = kwd[ANGLE]
        this.k1 = kwd[K1]
#         this.rho = this.l/this.angle
        this.init_matrixM()
        this.update_matrixM()
        this.update_matrixT()
    
    inline void update_matrixM()nogil:
        l = this.l
        angle = this.angle
        halfangle = angle/2.0
        rho =this.l/angle
        K = 1/rho**2-this.k1
        tanhalfangle = tan(halfangle)
        if K>0.0:
            k1 = sqrt(K)
            C=cos(k1*l)
            S=sin(k1*l)
            this.M[0][0]=C+S*tanhalfangle/(rho*k1)
            this.M[0][1]=S/k1
            this.M[1][0]=2.0*C*tanhalfangle/rho + tanhalfangle**2*S/(rho*rho*k1) - k1*S
            this.M[1][1]=C+S*tanhalfangle/(rho*k1)
        elif K<0.0:
            k1 = sqrt(-K)
            C=cosh(k1*l)
            S=sinh(k1*l)
            this.M[0][0]=C+S*tanhalfangle/(rho*k1)
            this.M[0][1]=S/k1
            this.M[1][0]=2.0*C*tanhalfangle/rho + tanhalfangle**2*S/(rho*rho*k1) + k1*S
            this.M[1][1]=C+S*tanhalfangle/(rho*k1)
        else:
            this.M[0][0]=1+angle*tanhalfangle
            this.M[0][1]=l
            this.M[1][0]=2*tan(halfangle)*(1+halfangle*tan(halfangle))/rho
            this.M[1][1]=1+angle*tanhalfangle
            
        
        
        # this.M[0][1] = rho*sin(angle)
        
        this.M[0][5] = rho*(1-cos(angle))
        
        this.M[1][5] = 2*tan(halfangle)
        
        this.M[2][2] = 1-angle*tan(halfangle)
        this.M[2][3] = this.l
        
        this.M[3][2] = -2*tan(halfangle)*(1-halfangle*tan(halfangle))/rho
        this.M[3][3] = 1-angle*tan(halfangle)
        #this.M[4][1] =-rho*(1-cos(angle))
        
    
    inline void get_phase(double* tws0, double* tws)nogil:
        l = this.l
        angle = this.angle
        tws[NUX] = tws0[NUX] + (atan(tws0[GAMMAX]*tan(angle)*l/angle - tws0[ALPHAX]) + atan(tws0[ALPHAX]))/(2*pi)
        tws[NUY] = tws0[NUY] + (atan(tws0[GAMMAY]*l - tws0[ALPHAY]) + atan(tws0[ALPHAY]))/(2*pi)
    
    
    inline void update(double* kwd)nogil:
        this.l = kwd[L]
        this.angle=kwd[ANGLE]
        this.k1 = kwd[K1]
        this.update_matrixM()
        this.update_matrixT()