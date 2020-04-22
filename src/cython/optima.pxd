from libcpp.vector cimport vector

cdef class Optima:
    cdef:
        vector[vector[int]] index
        list opt_parameters
        dict parameters2index
    
    
    cdef void get_optima(self, tuple args, dict kargs)