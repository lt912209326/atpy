from cymem.cymem cimport Pool
from libc.stdlib cimport free
from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole,eye

cdef class Element:
    cdef CppElement*  elem
    cdef Pool mem
    cdef str name
    cdef bint owner
    cdef str element_type
    cdef double* parms
    cdef double[20] twiss
    cdef int    ntwiss
    cdef list  eids



    cdef CppElement* bind2element(self)

    cdef double* bind2twiss(self)

    cdef double* bind2parameter(self)
    
    
#     cdef void bind2element(self, CppElement* elem)

#     cdef void bind2twiss(self, double* twiss)

#     cdef void bind2parameter(self, double* parms)




cdef class Marker(Element):
    pass


cdef class Drift(Element):
    pass
    

cdef class Dipole(Element):
    pass


cdef class Quadrupole(Element):
    pass


cdef class Sextupole(Element):
    pass

    
cdef class Octupole(Element):
    pass

