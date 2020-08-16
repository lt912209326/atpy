from ..cython cimport Element,Marker,Drift, Dipole, Quadrupole
from .variables cimport Variables
from .constraints cimport Constraints


from .optima cimport Optima

from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole

from libc.math cimport fmax,fabs,pi
from libc.string cimport memcpy
from cymem.cymem cimport Pool
cimport cython
import numpy as np

cdef class Lattice:
    cdef:
        Pool mem
        CppElement** elems
        double**     parms
        double**     twiss
        double[6][6]   matrix,temp_matrix
        Variables    variable
        Constraints  constraint
        Optima       optimum
        bint     cell
        int      nseq
        
        readonly Variables    pyvariable
        readonly Constraints  pyconstraint
        readonly Optima       pyoptimum
        readonly list         pyseq, pyelems
    
    
    def __cinit__(self,*args, bint cell=True):
        cdef Element temp_elem
        self.mem = Pool()
        self.cell = cell
        self.pyseq = []
        self.pyelems = []
        self.nseq = len(args)
        
        self.twiss =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.parms =<double**>self.mem.alloc(self.nseq,sizeof(double*))
        self.elems =<CppElement**>self.mem.alloc(self.nseq,sizeof(CppElement*))
        
        
        assert all([isinstance(arg,Element) for arg in args]), 'Args must be Element class or its subclasses'
        for index,elem in enumerate(args):
            temp_elem = elem if elem not in self.pyelems else elem.copy()
            self.pyseq.append(temp_elem)
            if elem not in self.pyelems:
                self.pyelems.append(temp_elem)
            self.twiss[index] = (<Element?>self.pyseq[index]).bind2twiss()
            self.elems[index] = (<Element?>self.pyseq[index]).bind2element()
            self.parms[index] = (<Element?>self.pyseq[index]).bind2parameter()
            

    def set_variables(self, Variables var):
        self.variable = var
        self.pyvariable = var

        
    def set_constraints(self, Constraints constraint):
        self.constraint = constraint
        self.pyconstraint = constraint

        
    def set_optima(self, Optima optimum):
        self.optimum = optimum
        self.pyoptimum = optimum
        
    
    cdef void _matrix_multiply(self,double* returned, double* leftmatrix, double* temp)nogil:
        cdef int i,j,k
        for i in range(6):
            for j in range(6):
                returned[6*i + j] = 0
                for k in range(6):
                    returned[6*i + j] += leftmatrix[6*i + k]*temp[6*k + j]
        
    
    cdef void _evaluate_lattice(self)nogil:
        cdef int i
        for i in range(1,self.nseq):
            self.elems[i].get_twiss(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_dispersion(self.twiss[i-1],self.twiss[i])
            self.elems[i].get_phase(self.twiss[i-1],self.twiss[i])
            #self.elems[i].get_emit(self.twiss[i-1],self.twiss[i])
    
    
    cdef void update_variables(self, double[:] variables, int num_variables)nogil:
        cdef int i, j
        for j in range(num_variables):
            if self.variable.index[0][j]==1:
                self.twiss[0][ self.variable.index[1][j] ] = variables[j]
            elif self.variable.index[0][j]==0:
                self.parms[ self.variable.index[2][j] ][ self.variable.index[1][j] ] = variables[j]
        self.twiss[0][2] = (1 + self.twiss[0][1]**2)/self.twiss[0][0]
        self.twiss[0][5] = (1 + self.twiss[0][4]**2)/self.twiss[0][3]
        for j in self.variable.indexset:
            self.elems[j].update( &self.parms[j][0] )
    
    
    cdef void update_constraints(self, double[:] CV, int num_constraints)nogil:
        cdef int i, j, k, l, index, start,end, sizeofmatrix = sizeof(self.elems[0].M)
        cdef double lower, upper,temp_max
        for j in range(num_constraints):
            start = self.constraint.index[2][j] 
            end = self.constraint.index[3][j] 
            lower = self.constraint.bounds[0][j] 
            upper = self.constraint.bounds[1][j]
            index = self.constraint.index[1][j]
            if self.constraint.index[0][j]==0:
                CV[j] = fmax(0.0, fabs(self.twiss[ end ][index] - self.twiss[ start ][index] - (upper + lower)/2.0)-(upper - lower)/2.0 )
            if self.constraint.index[0][j]==1:
                temp_max = -1.0e8
                if end > start:
                    for i in range(start,end+1):
                        temp_max = fmax(temp_max, fabs(self.twiss[i][index] - (upper + lower)/2.0)-(upper - lower)/2.0)
                else:
                    temp_max = fmax(0.0, fabs(self.twiss[start][index] - (upper + lower)/2.0)-(upper - lower)/2.0 )
                CV[j] = temp_max
            if index==9 or index==10:
                CV[j]*=200.0*pi
            elif index==1 or index==4:
                CV[j]*=100.0
#             if self.constraint.index[0][j]==0:
#                 CV[j] = fmax(0.0,fabs(self.twiss[ end ][index] - self.twiss[ start ][index] - (upper + lower)/2.0) - (upper - lower)/2.0)
#             if self.constraint.index[0][j]==1:
#                 temp_max = -1.0e8
#                 if end > start:
#                     for i in range(start,end+1):
#                         temp_max = fmax(temp_max, fabs(self.twiss[i][index] - (upper + lower)/2.0) - (upper - lower)/2.0)
#                 else:
#                     temp_max = fabs(self.twiss[start][index] - (upper + lower)/2.0) - (upper - lower)/2.0
#                 CV[j] = fmax(0.0, temp_max)

                
    def print_matrix(self):
        print(self.matrix)

    cdef void update_optima(self, double[:] objectives, int num_optima)nogil:
        cdef int i, j, start,end
        for j in range(num_optima):
            if self.optimum.index[0][j]==0:
                objectives[j]= (self.twiss[ self.optimum.index[3][j] ][ self.optimum.index[1][j] ]
                                  -self.twiss[ self.optimum.index[2][j] ][ self.optimum.index[1][j] ]) #only the twiss parameters nu
            elif self.optimum.index[0][j]==1:
                if self.optimum.index[3][j] == 0:
                    objectives[j]= self.parms[ self.optimum.index[2][j] ][ self.optimum.index[1][j] ]
                if self.optimum.index[3][j] == 1:
                    objectives[j]= self.twiss[ self.optimum.index[2][j] ][ self.optimum.index[1][j] ]

    
    def match(self,double[:,:] variables, double[:,:] ObjV):
        '''
        get_results(variables, objectives, CV)
        '''
        cdef int i,j, num_variables= variables.shape[1],num_pop = variables.shape[0], num_constraints = self.constraint.index[0].size()
        cdef double* CV_ptr
        if CV_ptr is NULL: CV_ptr = <double*>self.mem.alloc(num_constraints,sizeof(double))
        cdef double[:] CV =<double[:num_constraints]>CV_ptr
        
        assert num_variables == self.variable.index[0].size(), "Input numpy array doesn't match the variables!"
        ObjV[:,0]=0.0
        for i in range(num_pop):
            self.update_variables(variables[i,:], num_variables)
            self._evaluate_lattice()
            self.update_constraints(CV[:], num_constraints)
            for j in range(num_constraints):
                ObjV[i,0]+=CV[j]
    
    def _match(self,double[:] variables ):
        cdef int num_variables=variables.shape[0],nobj=2
        cdef double[:] CV=np.ones(self.variable.index[0].size(),dtype='float')
#         cdef double[:] optimize=np.ones(nobj,dtype='float')
        self.update_variables(variables, num_variables)
        self._evaluate_lattice()
        
        self.update_constraints(CV, self.variable.index[0].size())
        return np.sum(CV)
#         self.update_optima(optimize,nobj)
#         return np.sum(optimize)
    
    def _constr(self,double[:] variables):
        cdef int num_variables=variables.shape[0]
        cdef double[:] CV=np.ones(self.variable.index[0].size(),dtype='float')
        self.update_variables(variables, num_variables)
        self._evaluate_lattice()
        self.update_constraints(CV, self.variable.index[0].size())
        return np.array(CV)
#         pass
    
    def get_results(self,double[:,:] variables, double[:,:] objectives, double[:,:] CV):
        '''
        get_results(variables, objectives, CV)
        '''
        cdef int i, num_variables= variables.shape[1], num_optima = objectives.shape[1], num_constraints=CV.shape[1], num_pop = variables.shape[0]
        
        assert num_variables == self.variable.index[0].size(), "Input numpy array doesn't match the variables!"
        assert num_optima == self.optimum.index[0].size(), "Input numpy array doesn't match the optima!"
        assert num_constraints == self.constraint.index[0].size(), "Input numpy array: {} doesn't match the constraints\
        :{}!".format(num_constraints,self.constraint.index[0].size())
        for i in range(num_pop):
            self.update_variables(variables[i,:], num_variables)
            self._evaluate_lattice()
            self.update_constraints(CV[i,:], num_constraints)
            self.update_optima(objectives[i,:], num_optima)
             
    
    def evaluate_lattice(self):
        self._evaluate_lattice()
    

