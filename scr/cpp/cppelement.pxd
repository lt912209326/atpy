
# from libc.stdlib cimport free

cdef cppclass CppElement:
    int elem_type
    double l
    double angle
    double k1
    double k2
    double e1
    double e2
    double M[6][6]
    double T[6][6]
    
    __init__():
        pass
    
    inline void init_matrixM():
        this.M = [[1., 0., 0., 0., 0., 0.],
                  [0., 1., 0., 0., 0., 0.],
                  [0., 0., 1., 0., 0., 0.],
                  [0., 0., 0., 1., 0., 0.],
                  [0., 0., 0., 0., 1., 0.],
                  [0., 0., 0., 0., 0., 1.]]
    
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
#         for i in range()
        parms[0] = this.T[0][0]*parms0[0] + this.T[0][1]*parms0[1] + this.T[0][2]*parms0[2]
        parms[1] = this.T[1][0]*parms0[0] + this.T[1][1]*parms0[1] + this.T[1][2]*parms0[2]
        parms[2] = this.T[2][0]*parms0[0] + this.T[2][1]*parms0[1] + this.T[2][2]*parms0[2]
        
        parms[3] = this.T[3][3]*parms0[3] + this.T[3][4]*parms0[4] + this.T[3][5]*parms0[5]
        parms[4] = this.T[4][3]*parms0[3] + this.T[4][4]*parms0[4] + this.T[4][5]*parms0[5]
        parms[5] = this.T[5][3]*parms0[3] + this.T[5][4]*parms0[4] + this.T[5][5]*parms0[5]
    
    
    inline void get_phase(double* parms0, double* parms)nogil:
        pass
    
    inline void get_dispersion(double* parms0, double* parms)nogil:
#         parms[6] = this.M[0][0]*parms0[6] + this.M[0][1]*parms0[7] + this.M[0][5]
#         parms[7] = this.M[1][0]*parms0[6] + this.M[1][1]*parms0[7] + this.M[1][5]
        parms[12] = this.M[0][0]*parms0[12] + this.M[0][1]*parms0[13] + this.M[0][5]
        parms[13] = this.M[1][0]*parms0[12] + this.M[1][1]*parms0[13] + this.M[1][5]
    
    inline void update(double* parms)nogil:
        pass
