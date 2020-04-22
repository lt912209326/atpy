

cdef class Variables:
    def __init__(self,*args):
        '''
        def __init__(self,*args):
        args:( {'place':(0,15),'variables':{'l':(1,4) [,'beta_x':(0.5,1)]} },
               {},
               {})
        'place': n or (n1,n2)
        'variables': {'alpha_x':1 or (1,4) [,'beta_x':1 or (0.5,1)]}
        '''
        self.parameters2index = {'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5,
                            'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,'beta_z':6,'alpha_z':7,'gamma_z':8,
                            'nu_x':9,'nu_y':10,'nu_z':11,'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17}
        
        self.var_parameters = ['l','angle','k1','k2','e1','e2','beta_x','alpha_x','beta_y','alpha_y','beta_z','alpha_z',
                                    'eta_x','etap_x','eta_y','etap_y','eta_z','etap_z']
        
        self.index.resize(3)
        self.bounds.resize(2)
        self.indexset = cppset[int]()
        for arg in args:
            assert type(arg) is dict, 'Wrong args was input, please check the check the init doc!'
            place = arg.pop('place')
            variables = arg.pop('variables')
            assert type(variables) is dict, 'Wrong args was input for variables, please check the check the init doc!'
            
            assert all([key in self.var_parameters  for key in variables.keys()]), 'Wrong parameter as variable!'
            if type(place) is tuple:
                assert 0<=place[0]<place[1], 'The second number is smaller than the first for the place tuple!'
                self.get_variables(place, variables)
            elif type(place) is int:
                assert place >= 0, 'Place parameter must be integer > 0, you might need variables for this element!'
                if place !=0 : assert all([key in self.var_parameters[:6]  for key in variables.keys()]), 'Wrong parameter as variable at this place!'
                self.get_variables((place, place), variables)
        
    
    
    cdef void get_variables(self, tuple args, dict kargs):
    
        cdef int i, j,start = args[0], end = args[1], index4twiss
        cdef double lower, upper, tol=10e-8
            
        for key,value in kargs.items():
            lower,upper = value
            assert lower < upper, 'Variable lower is not smaller than upper!'
            index4twiss = self.parameters2index[key]
            index = 0 if key in self.var_parameters[:6] else 1
            for i in range(start,end+1):
                self.bounds[0].push_back(lower)
                self.bounds[1].push_back(upper)
                
                self.indexset.insert(i)
                self.index[0].push_back(index) # 0/1 mean parms/twiss
                self.index[1].push_back(index4twiss) # specific parms/twiss
                self.index[2].push_back(i) # place
        
    
    
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
        cdef int i, num_vars = self.bounds[0].size()
        cdef list index, indx1=[], indx2=[], indx3=[]
        for i in range(num_vars):
            indx1.append(self.index[0][i])
            indx2.append(self.index[1][i])
            indx3.append(self.index[2][i])
        index = [indx1,indx2,indx3]
        return index
    
    
    def __dealloc__(self):
        cdef int i
        for i in range(2):
            vector[double]().swap(self.bounds[i])
        for i in range(3):
            vector[int]().swap(self.index[i])
        vector[vector[double]]().swap(self.bounds)
        vector[vector[int]]().swap(self.index)
        
        
        
