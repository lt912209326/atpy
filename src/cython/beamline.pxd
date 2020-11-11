# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 09:59:11 2020

@author: LT
"""

from libc.math cimport fmax,fabs,pi
from libc.math cimport sqrt as cysqrt
from libc.string cimport memcpy
from cymem.cymem cimport Pool
cimport cython

from ..cython.constants cimport*
from ..cython cimport Element,Marker,Drift, Dipole, Quadrupole
from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole



cdef class BeamLine:
    '''
    需要修改相应myproblem.py
    '''
    cdef:
        Pool mem
        str             name
        CppElement**    elems
        double**        kwd_properties
        double**        tws_properties
        double**        loc_properties
        #natural emittance, Synchrotron energy loss U0 [GeV], momentum compact factor alpha, 
        double*         glb_properties
        bint            cell
        int             nseq
        
        dict            elems_index
        list            pyseq, pyelems

    
    cdef void _matmul(self,double* returned, double* lm, double* rm)nogil
    
    
    cdef void _init_tws_properties(self,double* transfermatrix)nogil
    
    
    cdef void _evaluate_lattice(self)nogil
    
    
    cdef void _get_glb_properties(self)nogil
            
    
