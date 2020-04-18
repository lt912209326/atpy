import numpy as np
cdef class Element:
    def __cinit__(self,**kargs):
        cdef dict type2index = {'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,
                                'beta_z':6,'alpha_z':7,'gamma_z':8,'nu_x':9,'nu_y':10,'nu_z':11,
                                'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17}
        cdef dict parms2index ={'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5}
        
        cdef dict class2type ={Marker:'Marker', Drift:'Drift', Dipole:'Dipole', Quadrupole:'Quadrupole', Sextupole:'Sextupole', Octupole:'Octupole'}
        
        self.mem = Pool()
        self.parms = <double*>self.mem.alloc(6,sizeof(double))
        if self.__class__ in class2type.keys():
            self.element_type = class2type[self.__class__]
        for key,value in kargs.items():
            if key in type2index.keys():
                self.twiss[type2index[key]]= value
            elif key in parms2index.keys():
                self.parms[parms2index[key]]=value
            else:
                print('Error arg was input!')
            
        
    def __init__(self,**kargs):
        self.elem = new CppElement()
        self.owner = True
    
    def __dealloc__(self):
        if self.elem is not NULL and self.owner:
            del self.elem
        self.elem = NULL
    
    def copy(self):
        cls = self.__class__
        cdef Element cp = cls.__new__(cls)
        if cp.elem is not NULL:
            del cp.elem
        cp.elem = self.elem
        cp.parms= self.parms
        cp.owner = False
        return cp
    
    def __neg__(self):
        cdef int[6] index = [1,4,7,13,15,17]
        cdef int i
        cls = self.__class__
        cdef Element cp = cls.__new__(cls)
        if cp.elem is not NULL:
            del cp.elem
        cp.elem = self.elem
        cp.parms= self.parms
        for i in range(6):
            cp.twiss[index[i]] = -self.twiss[index[i]]
        cp.owner = False
        return cp

    cdef CppElement* bind2element(self):
        return self.elem

    cdef double* bind2twiss(self):
        return &self.twiss[0]

    cdef double* bind2parameter(self):
        return self.parms
    
    def get_value(self, str parms=None):
        assert parms is not None, 'Please input arg!'
        cdef dict type2index = {'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,
                                'beta_z':6,'alpha_z':7,'gamma_z':8,'nu_x':9,'nu_y':10,'nu_z':11,
                                'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17}
        cdef dict parms2index ={'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5}
        if parms in type2index.keys():
            return self.twiss[type2index[parms]]
        elif parms in parms2index.keys():
            return self.parms[parms2index[parms]]
        elif parms == 'twiss':
            return np.array([self.twiss[i] for i in range(18)])
        elif parms == 'parameters':
            return [self.element_type]+[self.parms[i] for i in range(6)]
        else:
            print('Error arg was input!')



cdef class Marker(Element):
        
    def __init__(self, **kargs):
        self.elem = new CppMarker()
        self.owner = True

cdef class Drift(Element):
    def __init__(self,**kargs):
        self.elem = new CppDrift(self.parms)
        self.owner = True
    

cdef class Dipole(Element):
    def __init__(self,**kargs):
        self.elem = new CppDipole(self.parms)
        self.owner = True

cdef class Quadrupole(Element):
    def __init__(self,**kargs):
        self.elem = new CppQuadrupole(self.parms)
        self.owner = True

cdef class Sextupole(Element):
    def __init__(self,**kargs):
        pass
    
cdef class Octupole(Element):
    def __init__(self,**kargs):
        pass
