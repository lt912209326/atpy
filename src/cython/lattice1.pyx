from ..cython cimport Element,Marker,Drift, Dipole, Quadrupole
from .variables cimport Variables
from .constraints cimport Constraints


from .optima cimport Optima

from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole

from libc.math cimport fmax,fabs,pi
from libc.math cimport sqrt as cysqrt
from libc.string cimport memcpy
from cymem.cymem cimport Pool
cimport cython
import numpy as np

cdef struct Matrix:
    double[6][6] r
    double       angle
    double[3]    loc


cdef class Lattice:
    '''
    需要修改相应myproblem.py
    '''
    cdef:
        Pool mem
        CppElement** elems
        double**     parms
        double**     twiss
        Matrix*         loc_properties
        double[8]       rad_integrals
        #natural emittance, Synchrotron energy loss U0 [GeV], momentum compact factor alpha, 
        double[8]       global_parameters
        Variables    variables
        Constraints  constraints
        Optima       optima
        bint     cell
        int      nseq
        
        double          circumference
        double          energy
        double          mass
        
        dict            parameters2index
        list            parameters
        #readonly 
        list         pyseq, pyelems
    
    
    def __cinit__(self,*args, bint cell=True):
        cdef Element temp_elem
        self.mem = Pool()
        self.cell = cell
        
        self.energy = 2500
        self.mass   = 0.511
        
        self.parameters2index = {'l':0, 'angle':1, 'k1':2, 'k2':3, 'e1':4, 'e2':5,
                            'beta_x':0,'alpha_x':1,'gamma_x':2,'beta_y':3,'alpha_y':4,'gamma_y':5,'beta_z':6,'alpha_z':7,'gamma_z':8,
                            'nu_x':9,'nu_y':10,'nu_z':11,'eta_x':12,'etap_x':13,'eta_y':14,'etap_y':15,'eta_z':16,'etap_z':17,
                            'emit':0,'chromx':6,'chromy':7}
        self.parameters = ['l','angle','k1','k2','e1','e2','beta_x','alpha_x','beta_y','alpha_y','beta_z','alpha_z',
                                    'eta_x','etap_x','eta_y','etap_y','eta_z','etap_z','nu_x','nu_y','nu_z','emit','chromx','chromy']
        self.pyseq = []
        self.pyelems = []
        self.nseq = len(args)
        
        self.twiss =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.parms =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.elems =<CppElement**>self.mem.alloc(self.nseq,sizeof(CppElement*))
        if self.cell==True:
            self.loc_properties = <Matrix*>self.mem.alloc(self.nseq,sizeof(Matrix) )
        
        assert all([isinstance(arg,Element) for arg in args]), 'Args must be Element class or its subclasses'
        for index,elem in enumerate(args):
            temp_elem = elem if elem not in self.pyelems else elem.copy()
            self.pyseq.append(temp_elem)
            if elem not in self.pyelems:
                self.pyelems.append(temp_elem)
            self.twiss[index] = (<Element?>self.pyseq[index]).bind2twiss()
            self.elems[index] = (<Element?>self.pyseq[index]).bind2element()
            self.parms[index] = (<Element?>self.pyseq[index]).bind2parameter()
            

    
    cdef void _matmul(self,double* returned, double* lm, double* rm)nogil:
        cdef int i,j,k
        for i in range(6):
            for j in range(6):
                returned[6*i+j] = (lm[6*i]*rm[j] + lm[6*i+1]*rm[6+j] + lm[6*i+2]*rm[12+j] 
                                    + lm[6*i+3]*rm[18+j] + lm[6*i+4]*rm[24+j] + lm[6*i+5]*rm[30+j] )
    
    
    cdef void _init_twiss(self,double* transfermatrix)nogil:
        cdef double[6][6] r66
        cdef double cscmux,cscmuy,cosmux,cosmuy
        memcpy(&r66[0][0], transfermatrix, sizeof(r66))
        if self.cell==True and fabs(r66[0][0] + r66[1][1]) <2.0 and fabs(r66[2][2] + r66[3][3]) <2.0:
            cosmux = 0.5*(r66[0][0] + r66[1][1])
            cosmuy = 0.5*(r66[2][2] + r66[3][3])
            cscmux = r66[0][1]/( fabs(r66[0][1])*cysqrt(1-cosmux*cosmux) ) #1/sin(mux)
            cscmuy = r66[2][3]/( fabs(r66[2][3])*cysqrt(1-cosmuy*cosmuy) ) #1/sin(muy)
            self.twiss[0][0] = r66[0][1]*cscmux
            self.twiss[0][1] = 0.5*(r66[0][0]-r66[1][1])*cscmux
            self.twiss[0][2] = -r66[1][0]*cscmux
            self.twiss[0][3] = r66[2][3]*cscmuy
            self.twiss[0][4] = 0.5*(r66[2][2]-r66[3][3])*cscmuy
            self.twiss[0][5] = -r66[3][2]*cscmuy
            self.twiss[0][12]= 0.5*(r66[0][5]*(1.0-r66[1][1]) + r66[0][1]*r66[1][5])/(1.0-cosmux)
            self.twiss[0][13]= 0.5*(r66[1][5]*(1.0-r66[0][0]) + r66[1][0]*r66[0][5])/(1.0-cosmux)
    
    
    cdef void _evaluate_lattice(self)nogil:
        cdef int i,j,k,p, n_integral = 8
        self.circumference = 0.0
        cdef double[:,:] r66
        for i in range(8):
            self.rad_integrals[i] = 0
            self.global_parameters[i]=0
        if self.cell==True:
            memcpy( self.loc_properties[0].r, self.elems[0].M, sizeof(self.elems[0].M) )
            for i in range(1,self.nseq):
                self._matmul(&self.loc_properties[i].r[0][0], &self.elems[i].M[0][0], &self.loc_properties[i-1].r[0][0])
        
            self._init_twiss(&self.loc_properties[self.nseq-1].r[0][0])
        
        for i in range(1,self.nseq):
            self.elems[i].get_twiss(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_dispersion(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_phase(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_rad_integral(self.twiss[i-1], &self.rad_integrals[0])
            self.elems[i].get_chrom(self.twiss[i-1],self.twiss[i])
            self.circumference += self.parms[i][0]
            self.loc_properties[i].angle = self.loc_properties[i-1].angle + self.elems[i].angle
            #total horizontal chromaticity
            self.global_parameters[6] += self.twiss[i][18]
            #total vertical chromaticity
            self.global_parameters[7] += self.twiss[i][19]
        self._get_global_parameters()
    
    cdef void _get_global_parameters(self)nogil:
        cdef double* I = &self.rad_integrals[0]
        cdef:
            double Cq = 3.83E-13
            double E0 = self.energy
            double gamma0 = self.energy/self.mass
        #natural emittance
        self.global_parameters[0] = Cq*gamma0*gamma0*I[4]/(I[1] - I[3])
        #momentum compaction factor
        self.global_parameters[1] = I[0]/self.circumference
        #energy loss U0 [GeV]
        self.global_parameters[2] = 1.404*E0**4*I[1]
        #square of Energy dispersion
        self.global_parameters[3] = Cq*gamma0*gamma0*I[2]/(2*I[1] + I[3])
        #spin polarization I[6]:I3a in S.Y. Lee's book
        self.global_parameters[4] = -1.6/cysqrt(3)*I[6]/I[2]
        #damping factor D Jx=1-D,Jy=1,Jz=2+D
        self.global_parameters[5] = I[3]/I[1]
        
        
    cdef void update_variables(self, double[:] variables, int num_variables)nogil:
        cdef int i, j
        for j in range(num_variables):
            if self.variables.index[0][j]==1:
                self.twiss[0][ self.variables.index[1][j] ] = variables[j]
            elif self.variables.index[0][j]==0:
                self.parms[ self.variables.index[2][j] ][ self.variables.index[1][j] ] = variables[j]
        self.twiss[0][2] = (1 + self.twiss[0][1]**2)/self.twiss[0][0]
        self.twiss[0][5] = (1 + self.twiss[0][4]**2)/self.twiss[0][3]
        for j in self.variables.indexset:
            self.elems[j].update( &self.parms[j][0] )
    
    
    cdef void update_constraints(self, double[:] CV, int num_constraints, double[:] scale)nogil:
        cdef int i, j, k, l, index, start,end, sizeofmatrix = sizeof(self.elems[0].M)
        cdef double lower, upper,temp_max
        for j in range(num_constraints):
            start = self.constraints.index[2][j] 
            end = self.constraints.index[3][j] 
            lower = self.constraints.bounds[0][j] 
            upper = self.constraints.bounds[1][j]
            index = self.constraints.index[1][j]
            if self.constraints.index[0][j]==0:
                #phase advance or twiss differnce between two position
                #CV[j] = fmax(0.0, fabs(self.twiss[ end ][index] - self.twiss[ start ][index] - (upper + lower)/2.0)-(upper - lower)/2.0 )
                CV[j] = scale[j]*fmax(0.0, fabs(self.twiss[ end ][index] - self.twiss[ start ][index] - 0.5*(upper + lower))-0.5*(upper - lower) )
            elif self.constraints.index[0][j]==1:
                temp_max = -1.0e8
                if end > start:
                    for i in range(start,end+1):
                        temp_max = fmax(temp_max, fabs(self.twiss[i][index] - (upper + lower)/2.0)-(upper - lower)/2.0)
                else:
                    temp_max = fmax(0.0, fabs(self.twiss[start][index] - (upper + lower)/2.0)-(upper - lower)/2.0 )
                CV[j] = scale[j]*fmax(0.0, temp_max)
            #elif self.constraints.index[0][j]==1:
                #twiss parameter at start position
            #    CV[j] = scale[j]*fmax(0.0, fabs(self.twiss[start][index] - 0.5*(upper + lower))-0.5*(upper - lower))
            elif self.constraints.index[0][j]==2:
                #global parameters such as emit, U0...
                CV[j] = scale[j]*fmax(0.0, fabs(self.global_parameters[index] - 0.5*(upper + lower))-0.5*(upper - lower) )
            elif self.constraints.index[0][j]==3:
                CV[j] = scale[j]*fmax(0.0, fabs(self.loc_properties[end].angle - self.loc_properties[start].angle - 0.5*(upper + lower))-0.5*(upper - lower) )
        if self.cell==True:
            CV[num_constraints] = (fmax(0.0,fabs( self.loc_properties[self.nseq].r[0][0] + self.loc_properties[self.nseq].r[1][1] ) - 2.0) 
                                   +fmax(0.0,fabs( self.loc_properties[self.nseq].r[2][2] + self.loc_properties[self.nseq].r[3][3] ) - 2.0) )
            #if CV[num_constraints]> 1.0e-20:
            #    for j in range(num_constraints):
            #        CV[j] = 1.0e8

    

    cdef void update_optima(self, double[:] objectives, int num_optima)nogil:
        cdef int i, j, start,end,index
        for j in range(num_optima):
            index = self.optima.index[1][j]
            start = self.optima.index[2][j]
            end = self.optima.index[3][j]
            
            if self.optima.index[0][j]==0:
                objectives[j]= (self.twiss[ end ][ index ] - self.twiss[ start ][ index ]) #only the twiss parameters nu
            elif self.optima.index[0][j]==1:
                if end == 0:
                    objectives[j]= self.parms[ start ][ index ]
                if end == 1:
                    objectives[j]= self.twiss[ start ][ index ]
            elif self.optima.index[0][j]==2:
                objectives[j] = fabs( self.global_parameters[index] )
    
    
    def moea_without_CV(self, double[:,:] variables, double[:,:] objectives, double[:] scale, double factor=1.0e8, bint addCV=False):
        '''
        moea_with_CV(double[:,:] variables, double[:,:] objectives, double[:] scale)
        '''
        cdef:
            int i, num_variables= variables.shape[1], num_optima = objectives.shape[1], num_pop = variables.shape[0]
            int NCV = scale.shape[0]
            double totalCV = 0.0
        cdef double[:] CV = np.zeros( NCV , dtype='float')
        
        assert num_variables == self.variables.index[0].size(), "Input numpy array doesn't match the variables!"
        if addCV is True:
            assert num_optima == self.optima.index[0].size(), "Input numpy array doesn't match the optima!"
        else:
            assert num_optima == self.optima.index[0].size()+1, "Input numpy array doesn't match the optima!"
        for i in range(num_pop):
            self.update_variables(variables[i,:], num_variables)
            self._evaluate_lattice()
            self.update_constraints(CV, NCV-1, scale)
            self.update_optima(objectives[i,:], self.optima.index[0].size())
            totalCV = 0.0
            for j in range(NCV):
                totalCV += CV[j]
            if addCV is True:
                for j in range(num_optima):
                    objectives[i,j]+=factor*totalCV
            elif addCV is False:
                objectives[i,num_optima-1]=totalCV
            
    
    
    def get_results(self,double[:,:] variables, double[:,:] objectives, double[:,:] CV, double[:] scale):
        '''
        get_results(variables, objectives, CV)
        '''
        cdef int i, num_variables= variables.shape[1], num_optima = objectives.shape[1], num_constraints=CV.shape[1], num_pop = variables.shape[0]
        
        assert num_variables == self.variables.index[0].size(), "Input numpy array doesn't match the variables!"
        assert num_optima == self.optima.index[0].size(), "Input numpy array doesn't match the optima!"
        if self.cell==True:
            assert num_constraints == self.constraints.index[0].size()+1 , "Input numpy array: {} doesn't match the constraints\
            :{}!".format(num_constraints,self.constraints.index[0].size())
        else:
            assert num_constraints == self.constraints.index[0].size() , "Input numpy array: {} doesn't match the constraints\
            :{}!".format(num_constraints,self.constraints.index[0].size())
        for i in range(num_pop):
            self.update_variables(variables[i,:], num_variables)
            self._evaluate_lattice()
            self.update_constraints(CV[i,:], self.constraints.index[0].size(), scale)
            self.update_optima(objectives[i,:], num_optima)
             
    
    def evaluate_lattice(self):
        self._evaluate_lattice()
    
    
    cdef double _getitem(self,name,index=0):
        cdef int kind = self.parameters2index[name]
        if name in self.parameters[:6]:
            return self.parms[index][kind]
        elif name in self.parameters[6:21]:
            return self.twiss[index][kind]
        elif name in self.parameters[21:]:
            return self.global_parameters[kind]
    
    def __getitem__(self,arg):
        cdef int column, row
        if arg == 'cell':
            return self.cell
        elif type(arg) is dict:
            assert len(arg.keys())==1,'ArgsError:argument only receive dict with one key-value!'
            for key in arg.keys():
                position = key
            parameters=arg[key]
            assert type( position ) is tuple or type( position ) is int, \
            'key for the dict is tuple or int which are both index of element in the beamline!'
            assert type( parameters ) is list or type( parameters ) is str, \
            'value for the dict is list or str which are both name of the local parameters!'
            (row,position) = (len(position),position) if type(position) is tuple else (1,[position])
            (column,parameters) = (len(parameters),parameters) if type(parameters) is list else (1,[parameters])
            
            return_array = np.zeros((row,column))
            for i,value in enumerate(position):
                for j, name in enumerate(parameters):
                    
                    return_array[i,j]=self._getitem(name,value)
            return return_array
            
        elif type(arg) is list or type(arg) is tuple:
            return_array = np.zeros(len(arg))
            for index,name in enumerate(arg):
                assert type(name) is str,'arg must be the same types:like (str,str,str)'
                return_array[index] = self._getitem(name,self.nseq-1)
            return return_array
        elif type(arg) is str:
            if arg == 'variables':
                return self.variables
            elif arg == 'constraints':
                return self.constraints
            elif arg == 'optima':
                return self.optima
            elif arg == 'elements':
                return self.pyelems
            elif arg == 'beamline':
                return self.pyseq
            elif arg == 'rad_integrals':
                return self.rad_integrals
            elif arg in self.parameters[21:]:
                return self._getitem(arg)
            else:
                raise ValueError('Error str delivered to __getitem__()!')
    
    
    def __setitem__(self, index, value):
        if index == 'variables':
            assert isinstance(value,Variables),'Wrong type is given for variables!'
            self.variables = value
        elif index == 'constraints':
            assert isinstance(value,Constraints),'Wrong type is given for constraints!'
            self.constraints = value
        elif index == 'optima':
            assert isinstance(value,Optima),'Wrong type is given for optima!'
            self.optima = value
        
        
    def __neg__(self):
        cdef list args = self.pyseq[::-1]
        cls = self.__class__
        return cls.__new__(cls,*args )
    
    
    def add(self, Lattice other not None):
        cdef list args = self.pyseq + other.pyseq
        cls = self.__class__
        return cls.__new__(cls, *args)
        
    def __add__(self, Lattice other not None):
        return self.add(other)
    
    
    def write2file(self, str filename, str filetype='SAD'):
        cdef dict code2kind = {0:'MARK ', 1:'DRIFT', 2:'BEND ', 4:'QUAD ', 6:'SEXT ', 8:'OCTU '}
        cdef dict namehead = {0:'M', 1:'D', 2:'B', 4:'Q', 6:'S', 8:'O'}
        cdef dict index2parms={0:'L', 1:'ANGLE', 2:'K1', 3:'K2', 4:'E1', 5:'E2'}
        cdef double tol = 1.0e-8,value
        cdef int i, index, code, nparms=6
        cdef str elemstr
        cdef Element elem
        print('File name is : ', filename)
        with open(filename,'w+') as fn:
            for index,elem in enumerate(self.pyelems):
                code = elem.elem.elem_type
                elemstr= code2kind[ code ] + '    {}{:<4} = ('.format(namehead[code], index)
                for i in range(nparms):
                    if fabs(elem.parms[i]) > tol:
                        value = elem.parms[i]*elem.parms[0] if i==2 or i==3 else elem.parms[i]
                        elemstr = elemstr + '{:<5}={:.6e}  '.format(index2parms[i], value)
                elemstr = elemstr + ');\n'
                fn.write(elemstr)
            # elemstr = 'LINE    LINE1 = ('
            # for elem in self.pyseq:
            #     elemstr = elemstr + elem.name + ' '
            # elemstr = elemstr + ');\n'
            # fn.write(elemstr)
        
