
cdef class Optima:
    '''
    This class is to set optima for a lattice
    '''
        
    def __init__(self,*args):
        '''
        def __init__(self,*args):
        args:( {'place':(0,15),'optima':{'alpha_x':1 [,'beta_x':-1]} },
               {},
               {})
        'place': n or (n1,n2)
        'optima': {'alpha_x':1 or -1 [,'beta_x':1 or -1]} (1/-1:min/max)
        '''
        self.parameters2index = {'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5,
                            'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,'beta_z':6,'alpha_z':7,'gamma_z':8,
                            'nu_x':9,'nu_y':10,'nu_z':11,'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17}
        
        self.opt_parameters = ['l','angle','k1','k2','e1','e2','beta_x','alpha_x','beta_y','alpha_y','beta_z','alpha_z',
                                    'eta_x','etap_x','eta_y','etap_y','eta_z','etap_z','nu_x','nu_y','nu_z']
        self.index.resize(5)
        
        for arg in args:
            assert type(arg) is dict, 'Wrong args was input, please check the init doc!'
            place = arg.pop('place')
            optima = arg.pop('optima')
            assert type(optima) is dict, 'Wrong args was input for optima, please check the check the init doc!'
            
            assert all([key in self.opt_parameters  for key in optima.keys()]), 'Wrong parameter as optima!'
            if type(place) is tuple:
                assert 0<=place[0]<place[1], 'The second number is smaller than the first for the place tuple!'
                self.get_optima(place, optima)
            elif type(place) is int:
                assert place > 0, 'Place parameter must be integer > 0, you might need variables for this element!'
                self.get_optima((place, place), optima)
    
    
    
    cdef void get_optima(self, tuple args, dict kargs):
        cdef int i, j, start = args[0], end = args[1], place1, place2, index4parameters,minormax
        cdef double lower, upper, tol=10e-8
        
        for key,value in kargs.items():
            assert value==-1 or value== 1,'Minormax can only be value 1 or -1'
            minormax = value
            index4parameters = self.parameters2index[key]
            index = 0 if key in self.opt_parameters[18:] else 1
            place2= 0 if key in self.opt_parameters[:6] else 1
            for i in range(start,end+1):
                if key in self.opt_parameters[18:]:
                    place1,place2 = 0 if start == end else start, end
                elif key in self.opt_parameters[:18]:
                    place1 = i
                self.index[0].push_back(index) # 0/1 mean difference/value
                self.index[1].push_back(index4parameters) # specific parms/twiss
                self.index[2].push_back(place1) # place1
                self.index[3].push_back(place2) # place2 for nu or 0/1 for parms/twiss
                self.index[4].push_back(minormax) # place2 for nu or 0/1 for parms/twiss
                if key in self.opt_parameters[18:]:
                    break #nu_y2 - nu_y1 don't need loop
                    
    
    
    @property
    def index(self):
        cdef int i, num_vars = self.index[0].size()
        cdef list index, indx1=[], indx2=[], indx3=[], indx4=[], indx5=[]
        for i in range(num_vars):
            indx1.append(self.index[0][i])
            indx2.append(self.index[1][i])
            indx3.append(self.index[2][i])
            indx4.append(self.index[3][i])
            indx5.append(self.index[4][i])
        index = [indx1,indx2,indx3,indx4,indx5]
        return index
    
    
    def __dealloc__(self):
        cdef int i
        for i in range(5):
            vector[int]().swap(self.index[i])
        vector[vector[int]]().swap(self.index)
        
        
        
