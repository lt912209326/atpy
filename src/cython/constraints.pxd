from libcpp.vector cimport vector

cdef class Constraints:
    cdef:
        vector[vector[int]] index
        vector[vector[double]] bounds
        list cst_parameters
        dict twiss2index
        
    
    cdef void get_constraints(self, tuple args, dict kargs)