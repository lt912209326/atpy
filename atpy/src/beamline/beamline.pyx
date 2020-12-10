# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 10:04:35 2020

@author: LT
"""

import numpy as np

__all__ = ['BeamLine']

cdef class BeamLine:
    
    def __cinit__(self,str name,*args, bint cell=True,**kargs):
        cdef Element temp_elem
        self.mem = Pool()
        self.cell = cell
        self.name = name if name is not None else 'default'
        
        self.pyseq = []
        self.pyelems = []
        self.elems_index = {}
        self.nseq = len(args)
        
        self.kwd_properties =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.elems =<CppElement**>self.mem.alloc(self.nseq,sizeof(CppElement*))
        self.tws_properties =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.loc_properties = <double**>self.mem.alloc(self.nseq,sizeof(double*) )
        
        self.glb_properties = <double*>self.mem.alloc(GLB_NUM,sizeof(double) )
        
        self.glb_properties[ENERGY] = kargs['energy'] if 'energy' in kargs.keys() else 2500
        self.glb_properties[MASS0]  = kargs['mass0'] if 'mass0' in kargs.keys() else 0.511
        self.glb_properties[GAMMA] = self.glb_properties[ENERGY]/self.glb_properties[MASS0]
        
        cdef Element elem
        assert all([isinstance(arg,Element) for arg in args]), 'Args must be Element class or its subclasses'
        
        for index,elem in enumerate(args):
            self.pyseq.append(elem)
            if elem not in self.pyelems:
                self.pyelems.append(elem)
                (<Element?>elem).eids=[index]
                self.elems_index[elem.name] = [index]
            else:
                (<Element?>elem).eids.append(index)
                self.elems_index[elem.name].append(index)
            self.loc_properties[index] = <double*>self.mem.alloc(LOC_NUM,sizeof(double))
            self.tws_properties[index] = <double*>self.mem.alloc(TWS_NUM,sizeof(double))
            self.elems[index] = (<Element?>self.pyseq[index]).bind2element()
            self.kwd_properties[index] = (<Element?>self.pyseq[index]).bind2parameter()
        for key,value in kargs.items():
            if key in TWS_INDEX.keys():
                self.tws_properties[0][TWS_INDEX[key] ]=value

    
    cdef void _matmul(self,double* returned, double* lm, double* rm)nogil:
        cdef int i,j,k
        for i in range(6):
            for j in range(6):
                returned[6*i+j] = (lm[6*i]*rm[j] + lm[6*i+1]*rm[6+j] + lm[6*i+2]*rm[12+j] 
                                    + lm[6*i+3]*rm[18+j] + lm[6*i+4]*rm[24+j] + lm[6*i+5]*rm[30+j] )
    
    
    cdef void _init_tws_properties(self,double* transfermatrix)nogil:
        cdef double[6][6] r66
        cdef double cscmux,cscmuy,cosmux,cosmuy
        memcpy(&r66[0][0], transfermatrix, sizeof(r66))
        if self.cell==True and fabs(r66[0][0] + r66[1][1]) <2.0 and fabs(r66[2][2] + r66[3][3]) <2.0:
            cosmux = 0.5*(r66[0][0] + r66[1][1])
            cosmuy = 0.5*(r66[2][2] + r66[3][3])
            cscmux = r66[0][1]/( fabs(r66[0][1])*cysqrt(1-cosmux*cosmux) ) #1/sin(mux)
            cscmuy = r66[2][3]/( fabs(r66[2][3])*cysqrt(1-cosmuy*cosmuy) ) #1/sin(muy)
            self.tws_properties[0][BETAX ] = r66[0][1]*cscmux
            self.tws_properties[0][ALPHAX] = 0.5*(r66[0][0]-r66[1][1])*cscmux
            self.tws_properties[0][GAMMAX] = -r66[1][0]*cscmux
            self.tws_properties[0][BETAY ] = r66[2][3]*cscmuy
            self.tws_properties[0][ALPHAY] = 0.5*(r66[2][2]-r66[3][3])*cscmuy
            self.tws_properties[0][GAMMAY] = -r66[3][2]*cscmuy
            self.tws_properties[0][ETAX  ] = 0.5*(r66[0][5]*(1.0-r66[1][1]) + r66[0][1]*r66[1][5])/(1.0-cosmux)
            self.tws_properties[0][ETAPX ] = 0.5*(r66[1][5]*(1.0-r66[0][0]) + r66[1][0]*r66[0][5])/(1.0-cosmux)
    
    
    cdef void _evaluate_lattice(self)nogil:
        cdef: 
            int i,j,k,p
            double k2
        
        for i in range(ENERGY+1, GLB_NUM):
            self.glb_properties[i] = 0.0
        for i in range(LOC_NUM):
            self.loc_properties[0][i] = 0.0
        
        self.tws_properties[0][GAMMAX] = (1+ self.tws_properties[0][ALPHAX]**2)/self.tws_properties[0][BETAX]
        self.tws_properties[0][GAMMAY] = (1+ self.tws_properties[0][ALPHAY]**2)/self.tws_properties[0][BETAY]
        if self.cell==True:
            memcpy( &self.loc_properties[0][R11], &self.elems[0].M[0][0], sizeof(self.elems[0].M) )
            for i in range(1,self.nseq):
                self._matmul(&self.loc_properties[i][R11], &self.elems[i].M[0][0], &self.loc_properties[i-1][R11])
            self._init_tws_properties(&self.loc_properties[self.nseq-1][R11])
        
        for i in range(1,self.nseq):
            self.elems[i].get_twiss(self.tws_properties[i-1], self.tws_properties[i])
            # self.elems[i].get_dispersion(self.tws_properties[i-1], self.tws_properties[i])
            # self.elems[i].get_phase(self.tws_properties[i-1], self.tws_properties[i])
            self.elems[i].get_rad_integral(self.tws_properties[i-1], &self.glb_properties[0])
            # self.elems[i].get_chrom(self.tws_properties[i-1], self.tws_properties[i])
            
            self.elems[i].get_local(&self.loc_properties[i-1][0], &self.loc_properties[i][0])
            
            if self.elems[i].elem_type ==QUADRUPOLE:
                self.glb_properties[NAT_CHROMX] += self.tws_properties[i][CHROMX]
                self.glb_properties[NAT_CHROMY] += self.tws_properties[i][CHROMY]
            elif self.elems[i].elem_type ==SEXTUPOLE:
                k2 = self.kwd_properties[i][K2]*self.kwd_properties[i][L]
                self.glb_properties[TOTAL_K2X] += k2 if k2>0 else 0.0
                self.glb_properties[TOTAL_K2Y] += k2 if k2<0 else 0.0
            self.glb_properties[TOTALCHROMX] += self.tws_properties[i][CHROMX]
            self.glb_properties[TOTALCHROMY] += self.tws_properties[i][CHROMY]
        
        self.glb_properties[CIRCUMFERENCE] = self.loc_properties[self.nseq-1 ][S]
        self._get_glb_properties()
    
    cdef void _get_glb_properties(self)nogil:
        cdef double* I = &self.glb_properties[0]
        cdef:
            double Cq = 3.83E-13
            double E0 = self.glb_properties[ENERGY]
            double gamma0 = self.glb_properties[GAMMA]
        #natural emittance
        self.glb_properties[EMITX] = Cq*gamma0*gamma0*I[RI5]/(I[RI2] - I[RI4])
        #momentum compaction factor
        self.glb_properties[MOMENTUM_COMPACT_FACTOR] = I[RI1]/self.glb_properties[CIRCUMFERENCE ]
        #energy loss U0 [GeV]
        self.glb_properties[U0] = 1.404*E0**4*I[RI2]
        #square of Energy dispersion
        self.glb_properties[ENERGY_DISP2] = Cq*gamma0*gamma0*I[RI3]/(2*I[RI2] + I[RI4])
        #spin polarization I[6]:I3a in S.Y. Lee's book
        self.glb_properties[SPIN] = -1.6/cysqrt(3)*I[RI3A]/I[RI3]
        #damping factor D Jx=1-D,Jy=1,Jz=2+D
        self.glb_properties[DAMP_FACTOR] = I[RI4]/I[RI2]
    
    
    
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
