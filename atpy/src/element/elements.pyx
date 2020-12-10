import numpy as np
cdef class Element:
    def __cinit__(self,str name, **kargs):
        self.mem = Pool()
        self.name = name if name is not None else 'default'
        self.eids = []
        self.parms = <double*>self.mem.alloc(KWD_NUM, sizeof(double))
        cdef int i
        for i in range(KWD_NUM):
            self.parms[i]=0.0
        self.elem_kind = self.__class__.__name__
        for key,value in kargs.items():
            if key in KWD_INDEX.keys():
                self.parms[ KWD_INDEX[key] ]= value
            else:
                raise ValueError('Error arg was input!')
            
        
    def __init__(self,*arg, **kargs):
        self.elem = new CppElement()
        self.owner = True
    
    def __dealloc__(self):
        if self.elem is not NULL and self.owner:
            del self.elem
        self.elem = NULL
    
    def copy(self):
        cls = self.__class__
        cdef Element cp = cls.__new__(cls,self.name)
        if cp.elem is not NULL:
            del cp.elem
        cp.eids = self.eids
        cp.elem = self.elem
        cp.parms= self.parms
        cp.owner = False
        return cp
    
    
    #def __eq__(self, Element other not None):
    #    return True if self.name == other.name else False
    
        
    def __getitem__(self,index not None):
        if index in KWD_INDEX.keys():
            return self.parms[ KWD_INDEX[index] ]
        elif index == 'name':
            return self.name
        elif index == 'elem_kind':
            return self.elem_kind
        elif index == 'keywords':
            return [self.elem_kind]+[self.parms[i] for i in KWD_INDEX.values()]
        elif type(index) is int:
            if index<len(self.eids):
                return self.eids[index]
            else:
                raise ValueError('Index is too big!')
        else:
            raise ValueError('Error arg was input!')

    cdef CppElement* bind2element(self):
        return self.elem

    cdef double* bind2parameter(self):
        return self.parms



cdef class Marker(Element):
        
    def __init__(self, *args, **kargs):
        self.elem = new CppMarker()
        self.owner = True

cdef class Drift(Element):
    def __init__(self, *args, **kargs):
        self.elem = new CppDrift(self.parms)
        self.owner = True
    

cdef class Dipole(Element):
    def __init__(self, *args, **kargs):
        self.elem = new CppDipole(self.parms)
        self.owner = True

cdef class Quadrupole(Element):
    def __init__(self, *args, **kargs):
        self.elem = new CppQuadrupole(self.parms)
        self.owner = True

cdef class Sextupole(Element):
    def __init__(self, *args, **kargs):
        self.elem = new CppSextupole(self.parms)
        self.owner = True
    
cdef class Octupole(Element):
    def __init__(self, *args, **kargs):
        pass
