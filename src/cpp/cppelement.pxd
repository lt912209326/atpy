
# from libc.stdlib cimport free
from .radiation_integral cimport calc_sync_int,calc_chrom
from .constants cimport eye
from libc.math cimport fabs
from libc.string cimport memcpy

cdef cppclass CppElement:
    int elem_type
    double l
    double angle
    double k1
    double k2
    double e1
    double e2
    double parameters[6]
    double M[6][6]
    double T[6][6]
    
    __init__()nogil:
        pass
    
    inline double* get_parameters_ptr()nogil:
        return &this.parameters[0]
    
    inline void init_matrixM()nogil:
        memcpy(&this.M[0][0], &eye[0][0],sizeof(eye))
                
    
    inline void update_matrixM()nogil:
        pass
    
    inline void update_matrixT()nogil:
        M = this.M
        cdef int i
        for i in range(2):
            this.T[3*i][3*i+0] =  M[2*i][2*i]**2
            this.T[3*i][3*i+1] =  -2*M[2*i][2*i]*M[2*i][2*i+1]
            this.T[3*i][3*i+2] =  M[2*i][2*i+1]**2
            
            this.T[3*i+1][3*i+0] =  -M[2*i][2*i]*M[2*i+1][2*i]
            this.T[3*i+1][3*i+1] =  M[2*i][2*i]*M[2*i+1][2*i+1] + M[2*i][2*i+1]*M[2*i+1][2*i]
            this.T[3*i+1][3*i+2] =  -M[2*i][2*i+1]*M[2*i+1][2*i+1]
            
            this.T[3*i+2][3*i+0] =  M[2*i+1][2*i]**2
            this.T[3*i+2][3*i+1] =  -2*M[2*i+1][2*i]*M[2*i+1][2*i+1]
            this.T[3*i+2][3*i+2] =  M[2*i+1][2*i+1]**2
            
            
            
#         this.T = [[M[0][0]**2, -2*M[0][0]*M[0][1], M[0][1]**2, 0., 0., 0.],
#                   [-M[0][0]*M[1][0], M[0][0]*M[1][1] + M[0][1]*M[1][0], -M[0][1]*M[1][1], 0., 0., 0.],
#                   [M[1][0]**2, -2*M[1][0]*M[1][1], M[1][1]**2, 0., 0., 0.],
#                   [0., 0., 0., M[2][2]**2, -2*M[2][2]*M[2][3], M[2][3]**2],
#                   [0., 0., 0., -M[2][2]*M[3][2], M[2][2]*M[3][3] + M[2][3]*M[3][2], -M[2][3]*M[3][3]],
#                   [0., 0., 0., M[3][2]**2, -2*M[3][2]*M[3][3], M[3][3]**2]]
    
    inline void get_twiss(double* parms0, double* parms)nogil:
        parms[0] = this.T[0][0]*parms0[0] + this.T[0][1]*parms0[1] + this.T[0][2]*parms0[2]
        parms[1] = this.T[1][0]*parms0[0] + this.T[1][1]*parms0[1] + this.T[1][2]*parms0[2]
        parms[2] = this.T[2][0]*parms0[0] + this.T[2][1]*parms0[1] + this.T[2][2]*parms0[2]
        
        parms[3] = this.T[3][3]*parms0[3] + this.T[3][4]*parms0[4] + this.T[3][5]*parms0[5]
        parms[4] = this.T[4][3]*parms0[3] + this.T[4][4]*parms0[4] + this.T[4][5]*parms0[5]
        parms[5] = this.T[5][3]*parms0[3] + this.T[5][4]*parms0[4] + this.T[5][5]*parms0[5]
    
    
    inline void get_phase(double* parms0, double* parms)nogil:
        pass
    
    inline void get_dispersion(double* parms0, double* parms)nogil:
        parms[12] = this.M[0][0]*parms0[12] + this.M[0][1]*parms0[13] + this.M[0][5]
        parms[13] = this.M[1][0]*parms0[12] + this.M[1][1]*parms0[13] + this.M[1][5]
    
    inline void update(double* parms)nogil:
        pass
    
    inline void get_rad_integral(double* parms0, double* I)nogil:
        if this.l>1.0e-6 and fabs(this.angle) >1.0e-6:
            calc_sync_int(this.angle/this.l, this.l, this.k1, this.e1, this.e2, parms0[0], parms0[1], parms0[12], parms0[13], I)
        #else:
        #    memcpy(I,0,8*sizeof(double) )
    
    inline void get_chrom(double* parms0, double* parms)nogil:
        if (this.elem_type==4 or this.elem_type==6) and this.l>1.0e-6:
        #fabs(this.k1 )>1.0e-6 or fabs(this.k2)>1.0e-6:
            parms[18] = calc_chrom(this.l, this.k1, this.k2, parms0[0], parms0[1], parms0[12], parms0[13])
            parms[19] = calc_chrom(this.l, -this.k1, -this.k2, parms0[3], parms0[4], parms0[12], parms0[13])
        else:
            parms[18] = 0.0
            parms[19] = 0.0