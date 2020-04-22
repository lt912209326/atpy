from libcpp.vector cimport vector
from libcpp.set cimport set as cppset

cdef class Variables:
    cdef:
        vector[vector[int]] index
        vector[vector[double]] bounds
        cppset[int]  indexset
        list var_parameters
        dict parameters2index
        
    
    cdef void get_variables(self, tuple args, dict kargs)
