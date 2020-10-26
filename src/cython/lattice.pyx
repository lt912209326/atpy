# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 10:11:19 2020

@author: LT
"""

from ..cython cimport Element,Marker,Drift, Dipole, Quadrupole
from .variables cimport Variables
from .constraints cimport Constraints


from .optima cimport Optima

from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole
from .beamline cimport BeamLine

from libc.math cimport fmax,fabs,pi
from libc.math cimport sqrt as cysqrt
from libc.string cimport memcpy

cimport cython
import numpy as np


cdef class Lattice(BeamLine):
    cdef:
        Variables    variables
        Constraints  constraints
        Optima       optima
    def __init__(self,*args,**kargs):
        if 'cell' in kargs.keys():
            self.cell=<bint?>kargs['cell']
    
    
    
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
    
    
    cdef double update_constraints(self, double[:] CV, int num_constraints, double[:] scale,  bint isarray=True)nogil:
        cdef int i, j, k, l, index, start,end, sizeofmatrix = sizeof(self.elems[0].M)
        cdef double lower, upper,temp_max,constr_value=0.0, fitness = 0.0
        for j in range(num_constraints):
            start = self.constraints.index[2][j] 
            end = self.constraints.index[3][j] 
            lower = self.constraints.bounds[0][j] 
            upper = self.constraints.bounds[1][j]
            index = self.constraints.index[1][j]
            if self.constraints.index[0][j]==0:
                #phase advance or twiss differnce between two position
                #CV[j] = fmax(0.0, fabs(self.twiss[ end ][index] - self.twiss[ start ][index] - (upper + lower)/2.0)-(upper - lower)/2.0 )
                constr_value = scale[j]*fmax(0.0, fabs(self.twiss[ end ][index] - self.twiss[ start ][index] - 0.5*(upper + lower))-0.5*(upper - lower) )
            elif self.constraints.index[0][j]==1:
                temp_max = -1.0e8
                if end > start:
                    for i in range(start,end+1):
                        temp_max = fmax(temp_max, fabs(self.twiss[i][index] - (upper + lower)/2.0)-(upper - lower)/2.0)
                else:
                    temp_max = fmax(0.0, fabs(self.twiss[start][index] - (upper + lower)/2.0)-(upper - lower)/2.0 )
                constr_value = scale[j]*fmax(0.0, temp_max)
            #elif self.constraints.index[0][j]==1:
                #twiss parameter at start position
            #    CV[j] = scale[j]*fmax(0.0, fabs(self.twiss[start][index] - 0.5*(upper + lower))-0.5*(upper - lower))
            elif self.constraints.index[0][j]==2:
                #global parameters such as emit, U0...
                constr_value = scale[j]*fmax(0.0, fabs(self.global_parameters[index] - 0.5*(upper + lower))-0.5*(upper - lower) )
            elif self.constraints.index[0][j]==3:
                if index == 1:
                    constr_value = scale[j]*fmax(0.0, fabs(self.loc_properties[end].angle - self.loc_properties[start].angle - 0.5*(upper + lower))-0.5*(upper - lower) )
                if index == 2:
                    constr_value = scale[j]*fmax(0.0, fabs(self.loc_properties[end].s - self.loc_properties[start].s - 0.5*(upper + lower))-0.5*(upper - lower) )
            if isarray==True:
                CV[j] = constr_value
            else:
                fitness+=constr_value
        if self.cell==True:
            constr_value = scale[num_constraints]*(fmax(0.0,fabs( self.loc_properties[self.nseq].r[0][0] + self.loc_properties[self.nseq].r[1][1] ) - 2.0) 
                                   +fmax(0.0,fabs( self.loc_properties[self.nseq].r[2][2] + self.loc_properties[self.nseq].r[3][3] ) - 2.0) )
            if isarray==True:
                CV[num_constraints] = constr_value
                return 0.0
            else:
                fitness+=constr_value
                return fitness
        else:
            return fitness

    

    cdef void update_optima(self, double[:] objectives, int num_optima, double initoptima=0.0)nogil:
        cdef int i, j, start,end,index
        for j in range(num_optima):
            index = self.optima.index[1][j]
            start = self.optima.index[2][j]
            end = self.optima.index[3][j]
            
            if self.optima.index[0][j]==0:
                objectives[j]= initoptima + (self.twiss[ end ][ index ] - self.twiss[ start ][ index ]) #only the twiss parameters nu
            elif self.optima.index[0][j]==1:
                if end == 0:
                    objectives[j]= initoptima + self.parms[ start ][ index ]
                if end == 1:
                    objectives[j]= initoptima + self.twiss[ start ][ index ]
            elif self.optima.index[0][j]==2:
                objectives[j] = initoptima + fabs( self.global_parameters[index] )
    
    
    def moea_without_CV(self, double[:,:] variables, double[:,:] objectives, double[:] scale, double factor=1.0e8, bint addCV=False):
        '''
        moea_with_CV(double[:,:] variables, double[:,:] objectives, double[:] scale)
        '''
        cdef:
            int i, num_variables= variables.shape[1], num_optima = objectives.shape[1], num_pop = variables.shape[0]
            int NCV = scale.shape[0]
            double fitness = 0.0, initoptima=0.0
        cdef double[:] CV = np.zeros( NCV , dtype='float')
        
        if self.cell==True:
            assert NCV == self.constraints.index[0].size()+1 , "Input numpy array: {} doesn't match the constraints\
            :{}!".format(NCV,self.constraints.index[0].size())
        else:
            assert NCV == self.constraints.index[0].size() , "Input numpy array: {} doesn't match the constraints\
            :{}!".format(NCV,self.constraints.index[0].size())
            
        assert num_variables == self.variables.index[0].size(), "Input numpy array doesn't match the variables!"
        if addCV is True:
            assert num_optima == self.optima.index[0].size(), "Input numpy array doesn't match the optima!"
        else:
            assert num_optima == self.optima.index[0].size()+1, "Input numpy array doesn't match the optima!"
            
        if self.cell==True:
            NCV = NCV-1
        for i in range(num_pop):
            self.update_variables(variables[i,:], num_variables)
            self._evaluate_lattice()
            fitness = self.update_constraints(CV, NCV, scale, isarray=False)
            if addCV is True:
                initoptima = factor*fitness
            self.update_optima(objectives[i,:], self.optima.index[0].size(),initoptima)
            if addCV is False:
                objectives[i,num_optima-1]=fitness
            
    
    
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
        elif name in self.parameters[21:24]:
            return self.global_parameters[kind]
        elif name in self.parameters[24:25]:
            return self.loc_properties[index].angle
        elif name in self.parameters[25:26]:
            return self.loc_properties[index].s
    
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
    
    
    
    # def write2file(self, str filename, str filetype='SAD'):
        # cdef dict code2kind = {0:'MARK ', 1:'DRIFT', 2:'BEND ', 4:'QUAD ', 6:'SEXT ', 8:'OCTU '}
        # cdef dict namehead = {0:'M', 1:'L', 2:'B', 4:'Q', 6:'S', 8:'O'}
        # cdef dict index2parms={0:'L', 1:'ANGLE', 2:'K1', 3:'K2', 4:'E1', 5:'E2'}
        # cdef double tol = 1.0e-8,value
        # cdef int i, index, code, nparms=6
        # cdef str elemstr
        # cdef Element elem
        # print('File name is : ', filename)
        # with open(filename,'w+') as fn:
            # for index,elem in enumerate(self.pyelems):
                # code = elem.elem.elem_type
                # elemstr= code2kind[ code ] + '    {}{:<4} = ('.format(namehead[code], index)
                # for i in range(nparms):
                    # if fabs(elem.parms[i]) > tol:
                        # value = elem.parms[i]*elem.parms[0] if i==2 or i==3 else elem.parms[i]
                        # elemstr = elemstr + '{:<5}={:.6e}  '.format(index2parms[i], value)
                # elemstr = elemstr + ');\n'
                # fn.write(elemstr)
            # # elemstr = 'LINE    LINE1 = ('
            # # for elem in self.pyseq:
            # #     elemstr = elemstr + elem.name + ' '
            # # elemstr = elemstr + ');\n'
            # # fn.write(elemstr)
    
    
    def write2file(self, str filename, str filetype='SAD'):
        cdef dict code2kind = {0:'MARK ', 1:'DRIFT', 2:'BEND ', 4:'QUAD ', 6:'SEXT ', 8:'OCTU '}
        cdef dict namehead = {0:'M', 1:'L', 2:'B', 4:'Q', 6:'S', 8:'O'}
        cdef dict index2parms={0:'L', 1:'ANGLE', 2:'K1', 3:'K2', 4:'E1', 5:'E2'}
        cdef double tol = 1.0e-8,value
        cdef int i, index, code, nparms=6
        cdef str elemstr
        cdef Element elem
        print('File name is : ', filename)
        with open(filename,'w+') as fn:
            for index,elem in enumerate(self.pyelems):
                code = elem.elem.elem_type
                elemstr= code2kind[ code ] + '    {:<8} = ('.format(elem.name)
                for i in range(nparms):
                    if fabs(elem.parms[i]) > tol:
                        value = elem.parms[i]*elem.parms[0] if i==2 or i==3 else elem.parms[i]
                        elemstr = elemstr + '{:<5}={:.6e}  '.format(index2parms[i], value)
                elemstr = elemstr + ');\n'
                fn.write(elemstr)
            
            elemstr= 'LINE    LINE0 = ('
            for index,elem in enumerate(self.pyseq):
                code = elem.elem.elem_type
                elemstr = elemstr+ '{:<8}'.format(elem.name)
            elemstr = elemstr+ ');'
            fn.write(elemstr)
    
    def write2py(self, str filename):
        cdef dict code2kind = {0:'Marker', 1:'Drift', 2:'Dipole', 4:'Quadruple', 6:'Sextupole', 8:'Octupole'}
        cdef dict namehead = {0:'M', 1:'L', 2:'B', 4:'Q', 6:'S', 8:'O'}
        cdef dict index2parms={0:'l', 1:'angle', 2:'k1', 3:'k2', 4:'e1', 5:'e2'}
        cdef double tol = 1.0e-8,value
        cdef int i, index, code, nparms=6
        cdef str elemstr
        cdef Element elem
        print('File name is : ', filename)
        with open(filename,'w+') as fn:
            fn.write('from atpy import *\n')
            for index,elem in enumerate(self.pyelems):
                elemstr= '{0:<10} = {1:>10}({2:<10},'.format(elem.name, elem.__class__.__name__, f"'{elem.name}'")
                for i in range(nparms):
                    if fabs(elem.parms[i]) > tol:
                        value =  elem.parms[i]
                        elemstr = elemstr + '{:<5}={:.6e}  ,'.format(index2parms[i], value)
                elemstr = elemstr[:-1] + ')\n' if elemstr[-1]==',' else elemstr + ')\n'
                fn.write(elemstr)
            clsname = self.__class__.__name__
            elemstr= '{0:<10} = {1:>10}({2:<10},'.format(self.name, clsname, f"'{self.name}'")
            for index,elem in enumerate(self.pyseq):
                elemstr = elemstr+ '{:<10},'.format(elem.name)
            elemstr = elemstr[:-1] + ')\n' if elemstr[-1]==',' else elemstr + ')\n'
            fn.write(elemstr)
        

