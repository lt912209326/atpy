from .src.element.elements cimport Marker

cdef Marker MK0=Marker('MK0')
print(MK0)
print('test_cimport!')