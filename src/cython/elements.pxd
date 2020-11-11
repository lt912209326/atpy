from cymem.cymem cimport Pool
from libc.stdlib cimport free
from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole,eye
from ..cython.constants cimport*

cdef class Element:
    cdef CppElement*  elem
    cdef Pool       mem
    cdef str        name
    cdef bint       owner
    cdef str        elem_kind
    cdef double*    parms
    cdef list       eids



    cdef CppElement* bind2element(self)

    cdef double* bind2parameter(self)




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

