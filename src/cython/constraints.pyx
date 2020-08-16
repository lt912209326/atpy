
cdef class Constraints:
    '''
    This class is to set constraints for a lattice
    '''
        
    def __init__(self,*args,tol = 1.0e-4):
        '''
        def __init__(self,*args):
        args:( {'place':(0,15),'constraints':{'alpha_x':(1,4) [,'beta_x':(0.5,1)]} },
               {},
               {})
        'place': n or (n1,n2)
        'constraints': {'alpha_x':1 or (1,4) [,'beta_x':1 or (0.5,1)]}
        '''
        self.tolerance = tol
        self.cst_parameters = ['beta_x','alpha_x','beta_y','alpha_y','beta_z','alpha_z','eta_x',
                                    'etap_x','eta_y','etap_y','eta_z','etap_z','nu_x','nu_y','nu_z','telescope_x','telescope_y']
        
        self.twiss2index = {'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,
                                'beta_z':6,'alpha_z':7,'gamma_z':8,'nu_x':9,'nu_y':10,'nu_z':11,
                                'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17,'telescope_x':0,'telescope_y':2}
        self.index.resize(4)
        self.bounds.resize(2)
        for arg in args:
            assert type(arg) is dict, 'Wrong args was input, please check the check the init doc!'
            place = arg.pop('place')
            constraints = arg.pop('constraints')
            assert type(constraints) is dict, 'Wrong args was input for constraints, please check the check the init doc!'
            
            assert all([key in self.cst_parameters  for key in constraints.keys()]), 'Wrong parameter as constraint!'
            if type(place) is tuple:
                assert 0<=place[0]<place[1], 'The second number is smaller than the first for the place tuple!'
                self.get_constraints(place, constraints)
            elif type(place) is int:
                assert place > 0, 'Place parameter must be integer > 0, you might need variables for this element!'
                self.get_constraints((place, place), constraints)
    
    
    
    cdef void get_constraints(self, tuple args, dict kargs):
        cdef int i, j, k=0, m,start = args[0], end = args[1], place1, place2, index4twiss
        cdef double lower, upper, tol=self.tolerance
            
        for key,value in kargs.items():
            if type(value) is not tuple:
                lower, upper = value - tol, value + tol
            else:
                lower,upper = value
            assert lower < upper, 'Constraint lower is not smaller than upper!'
            index4twiss = self.twiss2index[key]
            
            index = 0 if key in self.cst_parameters[12:15] or (start != end and type(value) is not tuple) else 1
            if key in self.cst_parameters[15:]:
                index = 2
                
            place1 = start
            if key in self.cst_parameters[12:15]:
                place1 = 0 if start == end else start
            
            place2 = end
            self.bounds[0].push_back(lower)
            self.bounds[1].push_back(upper)
            self.index[0].push_back(index) # 0/1 mean difference/value
            self.index[1].push_back(index4twiss) # specific twiss
            self.index[2].push_back(place1) # place1
            self.index[3].push_back(place2) # place2
                    
    
    
    @property
    def bounds(self):
        cdef int i, num_vars = self.bounds[0].size()
        cdef list bnds, lower=[], upper=[]
        for i in range(num_vars):
            lower.append(self.bounds[0][i])
            upper.append(self.bounds[1][i])
        bnds = [lower,upper]
        return bnds
    
    
    @property
    def index(self):
        cdef int i, num_vars = self.index[0].size()
        cdef list index, indx1=[], indx2=[], indx3=[], indx4=[]
        for i in range(num_vars):
            indx1.append(self.index[0][i])
            indx2.append(self.index[1][i])
            indx3.append(self.index[2][i])
            indx4.append(self.index[3][i])
        index = [indx1,indx2,indx3,indx4]
        return index
    
    
    def __dealloc__(self):
        cdef int i
        for i in range(2):
            vector[double]().swap(self.bounds[i])
        for i in range(4):
            vector[int]().swap(self.index[i])
        vector[vector[double]]().swap(self.bounds)
        vector[vector[int]]().swap(self.index)
        
        
        
