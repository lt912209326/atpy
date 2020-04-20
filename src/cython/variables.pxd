from libcpp.vector cimport vector

cdef class Variables:
    cdef:
        vector[vector[int]] index
        vector[vector[double]] bounds
        list var_parameters
        dict parameters2index
        
    
    cdef void get_variables(self, tuple args, dict kargs)
