# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 10:04:35 2020

@author: LT
"""

import numpy as np


cdef class BeamLine:
    
    def __cinit__(self,str name,*args, bint cell=True):
        cdef Element temp_elem
        self.mem = Pool()
        self.cell = cell
        self.name = name if name is not None else 'default'
        self.energy = 2500
        self.mass   = 0.511
        
        self.parameters2index = {'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5,
                            'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,'beta_z':6,'alpha_z':7,'gamma_z':8,
                            'nu_x':9,'nu_y':10,'nu_z':11,'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17,
                            'emit':0,'chromx':6,'chromy':7,'angle':1,'s':2}
        self.parameters = ['l','angle','k1','k2','e1','e2','beta_x','alpha_x','beta_y','alpha_y','beta_z','alpha_z',
                                    'eta_x','etap_x','eta_y','etap_y','eta_z','etap_z','nu_x','nu_y','nu_z','emit','chromx','chromy','angle','s']
        self.pyseq = []
        self.pyelems = []
        self.nseq = len(args)
        
        self.twiss =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.parms =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.elems =<CppElement**>self.mem.alloc(self.nseq,sizeof(CppElement*))
        self.loc_properties = <Matrix*>self.mem.alloc(self.nseq,sizeof(Matrix) )
        #cdef Element elem
        assert all([isinstance(arg,Element) for arg in args]), 'Args must be Element class or its subclasses'
        for index,elem in enumerate(args):
            temp_elem = elem if elem not in self.pyelems else elem.copy()
            self.pyseq.append(temp_elem)
            if elem not in self.pyelems:
                self.pyelems.append(elem)
                (<Element?>elem).eids=[index]
            else:
                (<Element?>elem).eids.append(index)
            self.twiss[index] = (<Element?>self.pyseq[index]).bind2twiss()
            self.elems[index] = (<Element?>self.pyseq[index]).bind2element()
            self.parms[index] = (<Element?>self.pyseq[index]).bind2parameter()
        

    
    cdef void _matmul(self,double* returned, double* lm, double* rm)nogil:
        cdef int i,j,k
        for i in range(6):
            for j in range(6):
                returned[6*i+j] = (lm[6*i]*rm[j] + lm[6*i+1]*rm[6+j] + lm[6*i+2]*rm[12+j] 
                                    + lm[6*i+3]*rm[18+j] + lm[6*i+4]*rm[24+j] + lm[6*i+5]*rm[30+j] )
    
    
    cdef void _init_twiss(self,double* transfermatrix)nogil:
        cdef double[6][6] r66
        cdef double cscmux,cscmuy,cosmux,cosmuy
        memcpy(&r66[0][0], transfermatrix, sizeof(r66))
        if self.cell==True and fabs(r66[0][0] + r66[1][1]) <2.0 and fabs(r66[2][2] + r66[3][3]) <2.0:
            cosmux = 0.5*(r66[0][0] + r66[1][1])
            cosmuy = 0.5*(r66[2][2] + r66[3][3])
            cscmux = r66[0][1]/( fabs(r66[0][1])*cysqrt(1-cosmux*cosmux) ) #1/sin(mux)
            cscmuy = r66[2][3]/( fabs(r66[2][3])*cysqrt(1-cosmuy*cosmuy) ) #1/sin(muy)
            self.twiss[0][0] = r66[0][1]*cscmux
            self.twiss[0][1] = 0.5*(r66[0][0]-r66[1][1])*cscmux
            self.twiss[0][2] = -r66[1][0]*cscmux
            self.twiss[0][3] = r66[2][3]*cscmuy
            self.twiss[0][4] = 0.5*(r66[2][2]-r66[3][3])*cscmuy
            self.twiss[0][5] = -r66[3][2]*cscmuy
            self.twiss[0][12]= 0.5*(r66[0][5]*(1.0-r66[1][1]) + r66[0][1]*r66[1][5])/(1.0-cosmux)
            self.twiss[0][13]= 0.5*(r66[1][5]*(1.0-r66[0][0]) + r66[1][0]*r66[0][5])/(1.0-cosmux)
    
    
    cdef void _evaluate_lattice(self)nogil:
        cdef int i,j,k,p, n_integral = 8
        self.circumference = 0.0
        cdef double[:,:] r66
        for i in range(8):
            self.rad_integrals[i] = 0
            self.global_parameters[i]=0
        if self.cell==True:
            memcpy( self.loc_properties[0].r, self.elems[0].M, sizeof(self.elems[0].M) )
            for i in range(1,self.nseq):
                self._matmul(&self.loc_properties[i].r[0][0], &self.elems[i].M[0][0], &self.loc_properties[i-1].r[0][0])
        
            self._init_twiss(&self.loc_properties[self.nseq-1].r[0][0])
        self.loc_properties.angle = 0.0
        self.loc_properties.s = 0.0
        for i in range(1,self.nseq):
            self.elems[i].get_twiss(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_dispersion(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_phase(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_rad_integral(self.twiss[i-1], &self.rad_integrals[0])
            self.elems[i].get_chrom(self.twiss[i-1],self.twiss[i])
            self.circumference += self.parms[i][0]
            self.loc_properties[i].s = self.loc_properties[i-1].s + self.elems[i].l
            self.loc_properties[i].angle = self.loc_properties[i-1].angle + self.elems[i].angle
            #total horizontal chromaticity
            self.global_parameters[6] += self.twiss[i][18]
            #total vertical chromaticity
            self.global_parameters[7] += self.twiss[i][19]
        # self._get_global_parameters()
    
    cdef void _get_global_parameters(self)nogil:
        cdef double* I = &self.rad_integrals[0]
        cdef:
            double Cq = 3.83E-13
            double E0 = self.energy
            double gamma0 = self.energy/self.mass
        #natural emittance
        self.global_parameters[0] = Cq*gamma0*gamma0*I[4]/(I[1] - I[3])
        #momentum compaction factor
        self.global_parameters[1] = I[0]/self.circumference
        #energy loss U0 [GeV]
        self.global_parameters[2] = 1.404*E0**4*I[1]
        #square of Energy dispersion
        self.global_parameters[3] = Cq*gamma0*gamma0*I[2]/(2*I[1] + I[3])
        #spin polarization I[6]:I3a in S.Y. Lee's book
        self.global_parameters[4] = -1.6/cysqrt(3)*I[6]/I[2]
        #damping factor D Jx=1-D,Jy=1,Jz=2+D
        self.global_parameters[5] = I[3]/I[1]
    
    
    
    def evaluate_lattice(self):
        self._evaluate_lattice()
    
    
    
    def __neg__(self):
        cdef list args = self.pyseq[::-1]
        cls = self.__class__
        return cls.__new__(cls,*args )
    
    
        
    def __add__(BeamLine self, BeamLine other not None):
        cdef list args = self.pyseq + other.pyseq
        cls = self.__class__
        return cls.__new__(cls, *args)
    
    
    def __mul__(int other, BeamLine self):
        cdef list args
        if other > 0:
            args = other*self.pyseq
        elif other < 0:
            args = -other*self.pyseq[::-1]
        else:
            raise ValueError('BeamLine element cannot multiply with 0!')
        cls = self.__class__
        return cls.__new__(cls, *args)
