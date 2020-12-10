
# from libc.stdlib cimport free
from .radiation_integral cimport calc_sync_int,calc_chrom
from ..constants cimport*
from libc.math cimport fabs,sin,cos
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
        memcpy(&this.M[0][0], &EYE66[0][0],sizeof(EYE66))
                
    
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
            
            
    
    
    inline void get_twiss(double* tws0, double* tws)nogil:
        # pass
        tws[BETAX]  = this.T[0][0]*tws0[BETAX] + this.T[0][1]*tws0[ALPHAX] + this.T[0][2]*tws0[GAMMAX]
        tws[ALPHAX] = this.T[1][0]*tws0[BETAX] + this.T[1][1]*tws0[ALPHAX] + this.T[1][2]*tws0[GAMMAX]
        tws[GAMMAX] = this.T[2][0]*tws0[BETAX] + this.T[2][1]*tws0[ALPHAX] + this.T[2][2]*tws0[GAMMAX]
        
        tws[BETAY]  = this.T[3][3]*tws0[BETAY] + this.T[3][4]*tws0[ALPHAY] + this.T[3][5]*tws0[GAMMAY]
        tws[ALPHAY] = this.T[4][3]*tws0[BETAY] + this.T[4][4]*tws0[ALPHAY] + this.T[4][5]*tws0[GAMMAY]
        tws[GAMMAY] = this.T[5][3]*tws0[BETAY] + this.T[5][4]*tws0[ALPHAY] + this.T[5][5]*tws0[GAMMAY]
        this.get_phase(tws0, tws)
        this.get_dispersion(tws0, tws)
        this.get_chrom(tws0, tws)
    
    inline void get_global()nogil:
        pass
    
    
    inline void get_local(double* loc0, double* loc)nogil:
        # pass
        cdef double theta = loc0[THETA]+0.5*this.angle
        loc[S] = loc0[S] + this.l
        loc[X] = loc0[X] + this.l*cos(theta)
        loc[Y] = loc0[Y] + this.l*sin(theta)
        loc[THETA]=theta + 0.5*this.angle
        
        
    inline void get_phase(double* tws0, double* tws)nogil:
        pass
    
    inline void get_dispersion(double* tws0, double* tws)nogil:
        # pass
        tws[ETAX]  = this.M[0][0]*tws0[ETAX] + this.M[0][1]*tws0[ETAPX] + this.M[0][5]
        tws[ETAPX] = this.M[1][0]*tws0[ETAX] + this.M[1][1]*tws0[ETAPX] + this.M[1][5]
    
    inline void update(double* tws)nogil:
        pass
    
    inline void get_rad_integral(double* tws0, double* I)nogil:
        if this.l>1.0e-6 and fabs(this.angle) >1.0e-6:
            calc_sync_int(this.angle/this.l, this.l, this.k1, this.e1, this.e2, tws0[BETAX], tws0[ALPHAX], tws0[ETAX], tws0[ETAPX], I)
        #else:
        #    memcpy(I,0,8*sizeof(double) )
    
    inline void get_chrom(double* tws0, double* tws)nogil:
        if (this.elem_type==QUADRUPOLE or this.elem_type==SEXTUPOLE) and this.l>1.0e-6:
        #fabs(this.k1 )>1.0e-6 or fabs(this.k2)>1.0e-6:
            tws[CHROMX] = calc_chrom(this.l, this.k1, this.k2, tws0[BETAX], tws0[ALPHAX], tws0[ETAX], tws0[ETAPX])
            tws[CHROMY] = calc_chrom(this.l, -this.k1, -this.k2, tws0[BETAY], tws0[ALPHAY], tws0[ETAX], tws0[ETAPX])
        elif this.elem_type==DIPOLE and fabs(this.k1)>1.0e-6:
            tws[CHROMX] = calc_chrom(this.l, this.k1, this.k2, tws0[BETAX], tws0[ALPHAX], tws0[ETAX], tws0[ETAPX])
            tws[CHROMY] = calc_chrom(this.l, -this.k1, -this.k2, tws0[BETAY], tws0[ALPHAY], tws0[ETAX], tws0[ETAPX])
    