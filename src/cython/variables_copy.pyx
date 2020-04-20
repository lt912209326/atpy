import re

cdef class Variables:
    def __init__(self,*args):
        '''
        args:( {'place':'1~2','variables':{'l':'1~4' [,'beta_x':0.5~1]} },
               {},
               {},)
               'place':n or 'n1~n2'
               'variables':{{'l': '1~4' [,'beta_x':0.5~1]}}
        '''
        cdef int place,end
        cdef list var_parameters = ['l','angle','k1','k2','e1','e2','beta_x','alpha_x','beta_y','alpha_y','beta_z','alpha_z',
                                    'eta_x','etap_x','eta_y','etap_y','eta_z','etap_z']
        self.variable_index.resize(3)
        self.variable_bounds.resize(2)
        for arg in args:
            assert all([key in var_parameters  for key in arg['variables'].keys()]), 'Wrong parameter as variable!'
            if type(arg['place']) is str:
                assert all([key in var_parameters[:6]  for key in arg['variables'].keys()]), 'Wrong place for some parameter as variable!'
                pattern = re.compile(r'(?P<place>\d+)\s*~s*(?P<end>\d+)')
                mat = pattern.match(arg['place']).groupdict()
                place = int(mat['place'])
                end = int(mat['end'])
                self.get_vars(arg['variables'],place, end)
            elif type(arg['place']) is int:
                assert arg['place'] is 0 and all([key in var_parameters  for key in arg['variables'].keys()]), 'Wrong place for some parameter as variable!'
                self.get_vars(arg['variables'], arg['place'],0)
    
    
    cdef void get_vars(self, dict kargs, int place, int end):
        cdef int i, j, k=0, m
        cdef float lower, upper
        cdef dict parms2index = {'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5}
        cdef dict twiss2index = {'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,
                                'beta_z':6,'alpha_z':7,'gamma_z':8,'nu_x':9,'nu_y':10,'nu_z':11,
                                'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17}
        
        pattern = re.compile(r'(?P<lower>\+?-?\d+\.?\d*e*\+?-?\d*)\s*~s*(?P<upper>\+?-?\d+\.?\d*e*\+?-?\d*)')
        if end != 0:
            k = 0
        elif place != 0:
            end = place
            k = 0
        for i in range(place,end+1):
            for key,value in kargs.items():
                mat = pattern.match(value).groupdict()
                lower = float(mat['lower'])
                upper = float(mat['upper'])
                self.variable_bounds[0].push_back(lower)
                self.variable_bounds[1].push_back(upper)
                if end != 0:
                    m = parms2index[key]
                if key in parms2index.keys() and end == 0:
                    k = 0
                    m = parms2index[key]
                elif key in twiss2index.keys() and end == 0:
                    k = 1
                    m = twiss2index[key]
                elif end == 0:
                    print('Error: Wrong variable with ')
                self.variable_index[0].push_back(k) # 0/1 mean parms/twiss
                self.variable_index[1].push_back(i) # place
                self.variable_index[2].push_back(m) # specific parms/twiss
    
    
    @property
    def bounds(self):
        cdef int i, num_vars = self.variable_bounds[0].size()
        cdef list bnds, lower=[], upper=[]
        for i in range(num_vars):
            lower.append(self.variable_bounds[0][i])
            upper.append(self.variable_bounds[1][i])
        bnds = [lower,upper]
        return bnds
    
    
    @property
    def variable_index(self):
        cdef int i, num_vars = self.variable_bounds[0].size()
        cdef list index, indx1=[], indx2=[], indx3=[]
        for i in range(num_vars):
            indx1.append(self.variable_index[0][i])
            indx2.append(self.variable_index[1][i])
            indx3.append(self.variable_index[2][i])
        index = [indx1,indx2,indx3]
        return index
    
    
    def __dealloc__(self):
        cdef int i
        for i in range(2):
            vector[double]().swap(self.variable_bounds[i])
        for i in range(3):
            vector[int]().swap(self.variable_index[i])
        vector[vector[double]]().swap(self.variable_bounds)
        vector[vector[int]]().swap(self.variable_index)
        
        
        
