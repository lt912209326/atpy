# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 09:59:11 2020

@author: LT
"""


from ..cython cimport Element,Marker,Drift, Dipole, Quadrupole

from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole

from libc.math cimport fmax,fabs,pi
from libc.math cimport sqrt as cysqrt
from libc.string cimport memcpy
from cymem.cymem cimport Pool
cimport cython


cdef struct Matrix:
    double[6][6] r
    double       s
    double       angle
    double[3]    loc


cdef class BeamLine:
    '''
    需要修改相应myproblem.py
    '''
    cdef:
        Pool mem
        CppElement** elems
        double**     parms
        double**     twiss
        Matrix*         loc_properties
        double[8]       rad_integrals
        #natural emittance, Synchrotron energy loss U0 [GeV], momentum compact factor alpha, 
        double[8]       global_parameters
        bint     cell
        int      nseq
        
        double          circumference
        double          energy
        double          mass
        
        dict            parameters2index
        list            parameters
        str             name
        list            pyseq, pyelems
            

    
    cdef void _matmul(self,double* returned, double* lm, double* rm)nogil
    
    
    cdef void _init_twiss(self,double* transfermatrix)nogil
    
    
    cdef void _evaluate_lattice(self)nogil
    
    
    cdef void _get_global_parameters(self)nogil
            
    
