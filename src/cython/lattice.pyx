# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 10:11:19 2020

@author: LT
"""

from ..cython cimport Element,Marker,Drift, Dipole, Quadrupole
# from .variables cimport Variables
# from .constraints cimport Constraints

from ..cython.parser cimport*
from ..cython.structures cimport*
from ..cython.constants cimport*

# from .optima cimport Optima

from ..cpp cimport CppElement, CppMarker, CppDrift, CppDipole, CppQuadrupole
from .beamline cimport BeamLine

from libc.math cimport fmax,fabs,pi
from libc.string cimport memcpy
cimport cython
import numpy as np
import re


cdef class Lattice(BeamLine):
    cdef:
        int             constraints_num, vary_num,covary_num,optima_num
        bint            init_constraints, init_variables, init_optima, begin_match
        Parser          parser
        Variable*        variables
        Optimize*       optima
        Constraint*     constraints
    def __init__(self,*args,**kargs):
        if 'cell' in kargs.keys():
            self.cell=<bint?>kargs['cell']
            self.parser = Parser(self.elems_index)
            self.parser.set_database(self.kwd_properties, self.tws_properties, self.loc_properties, self.glb_properties)
            self.init_constraints =False
            self.init_variables =False
            self.init_optima =False
            self.begin_match =False
    
    
    cdef void set_constraints(self, dict constraints):
        cdef: 
            int i = 0
       
        self.init_constraints =True
        self.constraints_num = len(constraints)+1  if self.cell else len(constraints)
            
        self.constraints = <Constraint*>self.mem.alloc(self.constraints_num,sizeof(Constraint))
        for key,value in constraints.items():
            self.constraints[i] = CONSTRAINT0
            self.constraints[i].expr = self.parser.parse(key)
            if '<' in value.keys() and '>' in value.keys() and len(value)==2:
                if value['<'] < value['>']:
                    raise ValueError(f'Upper is less than upper in {i+1}-th constraint!')
                self.constraints[i].type = BOTH
                self.constraints[i].lower = value['>']
                self.constraints[i].upper = value['<']
            elif '>' in value.keys() and len(value)==1:
                self.constraints[i].type = LOWER
                self.constraints[i].lower = value['>']
            elif '<' in value.keys() and len(value)==1:
                self.constraints[i].type = UPPER
                self.constraints[i].upper = value['<']
            elif '=' in value.keys() and len(value)==1:
                self.constraints[i].type = EQUIV
                self.constraints[i].equiv = value['=']
            i += 1
        if self.cell:
            name = (<Element?>self.pyseq[-1]).name
            i_th = len( (<Element?>self.pyseq[-1]).eids )-1
            cell_token = 'DIM(ABS({0}[{1}].R11 + {0}[{1}].R22),2)+DIM(ABS({0}[{1}].R33 + {0}[{1}].R44),2)'.format(name,i_th)
            self.constraints[self.constraints_num-1] = CONSTRAINT0
            self.constraints[self.constraints_num-1].expr = self.parser.parse(cell_token)
            self.constraints[self.constraints_num-1].type = UPPER
            self.constraints[self.constraints_num-1].upper= 0
    
    
    cdef void set_variables(self, dict variables):
        cdef:
            int i, position, eid, index, data_kind=-1
            str name,index_name
            double** database
            Variable*   pvar
            Variable    tmp_var
            Variable*   phead
            Variable*   ptail
        self.init_variables =True
        self.vary_num = 0
        self.covary_num = 0
        self.variables = <Variable*>self.mem.alloc(1,sizeof(Variable))
        memcpy(self.variables,&VARIABLE0,sizeof(Variable))
        pvar = self.variables
        loc_token = re.compile(r'(?P<name>.+)\[(?P<position>\d+)\]\.(?P<index>\w+)')
        glb_token = re.compile(r'(?P<index>\w+)')
        
        for key,value in variables.items():
            tmp_var = VARIABLE0
            ##==============================================================================================##
            if '<' in value.keys() and '>' in value.keys() and len(value)==2:
                assert type(value['<']) is int or type(value['<']) is float,f'the value of {value} is not number'
                assert type(value['>']) is int or type(value['>']) is float,f'the value of {value} is not number'
                if value['<'] < value['>']:
                    raise ValueError(f'Upper is less than upper in {i+1}-th constraint!')
                tmp_var.type = BOTH
                tmp_var.lower = value['>']
                tmp_var.upper = value['<']
            elif '>' in value.keys() and len(value)==1:
                assert type(value['>']) is int or type(value['>']) is float,f'the value of {value} is not number'
                tmp_var.type = LOWER
                tmp_var.lower = value['>']
            elif '<' in value.keys() and len(value)==1:
                assert type(value['<']) is int or type(value['<']) is float,f'the value of {value} is not number'
                tmp_var.type = UPPER
                tmp_var.upper = value['<']
            elif '=' in value.keys() and len(value)==1:
                assert type(value['=']) is str,f'the value of {value} is not an expression!'
                tmp_var.type = EQUIV
                tmp_var.expr = self.parser.parse(value['='])
            ##==============================================================================================##
            if loc_token.search(key) is not None:
                elem_name,position_str,index_name = loc_token.search(key).group('name','position','index')     
                position = int(position_str)
            elif glb_token.search(key):
                elem_name = None
                index_name = glb_token.search(key).group('index')
            else:
                raise ValueError(f'Variable str: "{key}" is not correct!')
            ##==============================================================================================##
            if elem_name is None:
                if index_name not in GLB_INDEX.keys():
                    raise ValueError(f'{index_name} is not a valid global property')
                if pvar.vary is not NULL:
                    pvar.next=<Variable*>self.mem.alloc(1,sizeof(Variable))
                    pvar.next.prev=pvar
                    pvar = pvar.next
                index = GLB_INDEX[index_name ]
                tmp_var.data_kind = GLB
                tmp_var.index = index
                tmp_var.prev = pvar.prev
                memcpy(pvar,&tmp_var,sizeof(Variable))
                pvar.vary = &self.glb_properties[index]
                if tmp_var.type== EQUIV:
                    pvar.expr = tmp_var.expr
                    self.covary_num+=1
                else:
                    self.vary_num+=1
            else:
                if index_name in LOC_INDEX.keys():
                    tmp_var.data_kind = LOC
                    index = LOC_INDEX[index_name ]
                    database = self.loc_properties
                elif index_name in TWS_INDEX.keys():
                    tmp_var.data_kind = TWS
                    index = TWS_INDEX[index_name ]
                    database = self.tws_properties
                elif index_name in KWD_INDEX.keys():
                    tmp_var.data_kind = KWD
                    index = KWD_INDEX[index_name ]
                    database = self.kwd_properties
                else:
                    raise ValueError(f'{index_name} is not a valid keyword!')
                    
                name_token = re.compile(elem_name+'$')
                for elem_name in self.elems_index.keys():
                    name_match = name_token.search(elem_name)
                    if name_match is None:
                        continue
                    name = name_match.group()
                    eid = self.elems_index[name][position]
                    if pvar.vary is not NULL:
                        pvar.next=<Variable*>self.mem.alloc(1,sizeof(Variable))
                        pvar.next.prev=pvar
                        pvar = pvar.next
                    tmp_var.index = index
                    tmp_var.position = eid
                    tmp_var.prev = pvar.prev
                    memcpy(pvar,&tmp_var,sizeof(Variable))
                    pvar.vary = &database[eid][index]
                    if tmp_var.type== EQUIV:
                        pvar.expr = tmp_var.expr
                        self.covary_num+=1
                    else:
                        self.vary_num+=1
        ptail = pvar
        pvar = self.variables
        phead= self.variables
        cdef Variable* ptmp
        for i in range(self.vary_num+self.covary_num):
            if pvar is NULL:
                raise RuntimeError('linklist is outof memory!')
            # print( (<Element>self.pyseq[pvar.position]).name )
            ptmp = pvar
            pvar = pvar.next
            if ptmp.type==EQUIV:
                if ptmp.next is NULL:
                    continue
                if ptmp.prev is NULL:
                    phead=ptmp.next
                    ptmp.next.prev = NULL
                else:
                    ptmp.next.prev = ptmp.prev
                    ptmp.prev.next = ptmp.next
                        
                ptmp.prev=ptail
                ptail.next=ptmp
                ptail = ptmp
                ptmp.next=NULL
        self.variables = phead
        # pvar=phead 
        # while pvar is not NULL:
            # print(f'name: {(<Element?>self.pyseq[pvar.position]).name} {pvar.prev is not NULL}' )
            # pvar=pvar.next
    
    
    cdef void set_optima(self, dict optima):
        cdef: 
            int i = 0
        self.init_optima =True
        self.optima_num = len(optima)
        self.optima = <Optimize*>self.mem.alloc(self.optima_num, sizeof(Optimize))
        for key,value in optima.items():
            self.optima[i] = OPTIMIZE0
            self.optima[i].expr = self.parser.parse(key)
            if value =='min':
                self.optima[i].minormax = 1.0
            elif value=='max':
                self.optima[i].minormax = -1.0
            else:
                raise ValueError(f'Wrong minor max was input in {key}:{value} !')
            i+=1
        

    
    
    cdef inline void update_variables(self, double* variables):
        cdef: 
            int i
            Variable* pvar=self.variables
        for i in range(self.vary_num):
            pvar.vary[0] = variables[i]
            pvar=pvar.next
        for i in range(self.covary_num):
            pvar.vary[0] = pvar.expr.calc()
            pvar=pvar.next
        for i in range(self.nseq):
            self.elems[i].update(&self.kwd_properties[i][0])
        # print('update_variable',variables[0])
    
    
    cdef double update_constraints(self, double[:] CV, int num_constraints)nogil:
        return collect(self.constraints, self.constraints_num, &CV[0])


    cdef void update_optima(self, double* optima):
        cdef: 
            int i
            double  value
        for i in range(self.optima_num):
            value = self.optima[i].minormax*self.optima[i].expr.calc()
            optima[i] = value
    
    
    
    def get_results(self,double[:,:] variables, double[:,:] objectives, double[:,:] CV):
        '''
        get_results(variables, objectives, CV)
        '''
        cdef:
            int i, num_variables= variables.shape[1], num_optima = objectives.shape[1], num_constraints=CV.shape[1], num_pop = variables.shape[0]
            double** parameters
        if self.begin_match is False:
            if self.init_constraints and self.init_variables and self.init_optima:
                self.begin_match = True
            else:
                raise RuntimeError('Match conndition was not OK !')
        if num_variables != self.vary_num:
            raise ValueError("Input numpy array doesn't match the variables!")
        if num_optima != self.optima_num:
            raise ValueError("Input numpy array doesn't match the optima!")
        if num_constraints != self.constraints_num :
            raise ValueError("Input numpy array: {num_constraints} doesn't match the constraints:{self.constraints_num}!")
        for i in range(num_pop):
            self.update_variables(&variables[i,0])
            self._evaluate_lattice()
            collect(self.constraints, self.constraints_num, &CV[i,0])
            # parameters = self.constraints[0].expr.database()
            # print(' max OR min token:', self.constraints[0].expr.token)
            # print('range number: ',self.constraints[0].expr.item(0))
            # print(' range.token :',self.constraints[0].expr.item(1))
            # print(parameters[0][0],parameters[1][0],parameters[2][0],parameters[3][0])
            # print('QA1E.k1:',self.kwd_properties[2][K1],' Phen[0]:',variables[i,0])
            # print('optima[0]:',self.optima[0].expr.calc(),'optima[1]:',self.optima[1].expr.calc())
            # print('constraints[0]:',self.constraints[3].expr.calc())
            self.update_optima(&objectives[i,0])
             
    
    def evaluate_lattice(self):
        self._evaluate_lattice()
        
    
    cdef list get_variables(self, bint ifprint=False):
        cdef:
            int i
            list variables 
            dict index_kwd,index_tws,index_loc
            Variable* pvar
        variables = []
        index_bound={LOWER:'lower', UPPER:'upper', BOTH:'both', EQUIV:'covary'}
        index_kwd = {value:key for key,value in KWD_INDEX.items()}
        index_tws = {value:key for key,value in TWS_INDEX.items()}
        index_loc = {value:key for key,value in LOC_INDEX.items()}
        pvar = self.variables
        if ifprint:
            print(f'{70*"="}')
            print(f'{"Total Variables":<20}:{self.vary_num+self.covary_num:<10}, {"variable":<10}: {self.vary_num:<10}, {"covary":<10}:{self.covary_num:<10}')
            print(f'{"Position":<10} {"keyword":<8} {"bound":<6} {"lower":>13} {"upper":>13}')
        for i in range(self.vary_num+self.covary_num):
            position = (<Element?>self.pyseq[ pvar.position ]).name
            lower = pvar.lower
            upper=pvar.upper
            data_kind= pvar.data_kind
            if data_kind==KWD:
                keyword = index_kwd[pvar.index]
            elif data_kind==TWS:
                keyword = index_tws[pvar.index]
            elif data_kind==LOC:
                keyword = index_loc[pvar.index]
            bound= index_bound[pvar.type]
            if ifprint:
                print(f'{position:<10} {keyword:<8} {bound:<6} {lower:>13.2e} {upper:>13.2e}')
            if pvar.type != EQUIV:
                variables.append( {position:[ keyword, bound, lower, upper] } )
            pvar=pvar.next
        return variables
    
    
    
    cdef list get_constraints(self, bint ifprint=False):
        cdef:
            int i
            list constraints 
        constraints = []
        index_bound={LOWER:'lower', UPPER:'upper', BOTH:'both', EQUIV:'equiv'}
        if ifprint:
            print(f'{70*"="}')
            print(f'{"Total Constraints":<20}:{self.constraints_num:<10}, cell: {self.cell:<10}')
            print(f'{"Index":<10} {"bound":<6} {"lower":>13} {"upper":>13}')
        for i in range(self.constraints_num):
            lower = self.constraints[i].lower
            upper=self.constraints[i].upper
            bound = index_bound[ self.constraints[i].type ]
            if ifprint:
                print(f'{i:<10} {bound:<6} {lower:>13.2e} {upper:>13.2e}')
            constraints.append( {i:[ bound, lower, upper] } )
        return constraints
    
    
    
    cdef list get_optima(self, bint ifprint=False):
        cdef:
            int i
            list optima
        optima = []
        if ifprint:
            print(f'{70*"="}')
            print(f'{"Total Optima":<20}:{self.optima_num:<10}')
            print(f'{"Index":<10} {"min|max":<10}')
        for i in range(self.optima_num):
            minormax = 'min' if self.optima[i].minormax==1.0 else 'max'
            if ifprint:
                print(f'{i:<10} {minormax:<10}')
            optima.append( {i:int(self.optima[i].minormax) } )
        return optima
    
    
    cdef double getitem(self,name,position=0):
        cdef int kind
        if name in KWD_INDEX.keys():
            kind = KWD_INDEX[name]
            return self.kwd_properties[position][kind]
        elif name in TWS_INDEX.keys():
            kind =TWS_INDEX[name]
            return self.tws_properties[position][kind]
        elif name in LOC_INDEX.keys():
            kind = LOC_INDEX[name]
            return self.loc_properties[position][kind]
        elif name in GLB_INDEX.keys():
            kind = GLB_INDEX[name]
            return self.glb_properties[kind]
    
    
    def __getitem__(self,arg):
        cdef int column, row
        
        if type(arg) is dict:
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
                    
                    return_array[i,j]=self.getitem(name,value)
            return return_array
            
        elif type(arg) is list or type(arg) is tuple:
            return_array = np.zeros(len(arg))
            for index,name in enumerate(arg):
                assert type(name) is str,'arg must be the same types:like (str,str,str)'
                return_array[index] = self.getitem(name,self.nseq-1)
            return return_array
        elif type(arg) is str:
            
            if arg == 'cell':
                return self.cell
            elif arg =='kwd':
                return KWD_INDEX
            elif arg =='tws':
                return TWS_INDEX
            elif arg =='loc':
                return LOC_INDEX
            elif arg =='glb':
                return GLB_INDEX
            
            elif arg == 'variables_num':
                return self.vary_num
            elif arg == 'constraints_num':
                return self.constraints_num
            elif arg == 'optima_num':
                return self.optima_num
            
            elif arg == 'variables':
                return self.get_variables()
            elif arg == 'constraints':
                return self.get_constraints()
            elif arg == 'optima':
                return self.get_optima()
            
            elif arg == 'elements':
                return self.pyelems
            elif arg == 'beamline':
                return self.pyseq
            elif arg in TOTAL_INDEX.keys():
                return self.getitem(arg)
            else:
                raise ValueError('Error str delivered to __getitem__()!')
    
    
    def __setitem__(self, index, value):
        if index == 'variables':
            assert isinstance(value,dict),'Wrong type is given for variables!'
            self.set_variables(value)
            self.get_variables(ifprint=True)
        elif index == 'constraints':
            assert isinstance(value,dict),'Wrong type is given for constraints!'
            self.set_constraints(value)
            self.get_constraints(ifprint=True)
        elif index == 'optima':
            assert isinstance(value,dict),'Wrong type is given for optima!'
            self.set_optima(value)
            self.get_optima(ifprint=True)
    
    
    
    def write2file(self, str filename, str filetype='SAD'):
        cdef dict code2kind = {0:'MARK ', 1:'DRIFT', 2:'BEND ', 4:'QUAD ', 6:'SEXT ', 8:'OCTU '}
        cdef dict namehead = {0:'M', 1:'L', 2:'B', 4:'Q', 6:'S', 8:'O'}
        cdef dict index2parms={value:key.upper() for key,value in KWD_INDEX.items()}
        cdef double tol = 1.0e-8,value
        cdef int i, index, code, nparms=6
        cdef str elemstr
        cdef Element elem
        print('File name is : ', filename)
        with open(filename,'w+') as fn:
            for index,elem in enumerate(self.pyelems):
                code = elem.elem.elem_type
                elemstr= code2kind[ code ] + '    {:<8} = ('.format(elem.name)
                for i in range(KWD_NUM):
                    if fabs(elem.parms[i]) > tol:
                        value = elem.parms[i]*elem.parms[0] if i==K1 or i==K2 else elem.parms[i]
                        elemstr = elemstr + f'{index2parms[i]:<5}={value:.6e}  '
                elemstr = elemstr + ');\n'
                fn.write(elemstr)
            
            elemstr= 'LINE    LINE0 = ('
            for index,elem in enumerate(self.pyseq):
                code = elem.elem.elem_type
                elemstr = elemstr+ '{:<8}'.format(elem.name)
            elemstr = elemstr+ ');'
            fn.write(elemstr)
    
    def write2py(self, str filename):
        cdef dict namehead = {0:'M', 1:'L', 2:'B', 4:'Q', 6:'S', 8:'O'}
        cdef dict index2parms={value:key for key,value in KWD_INDEX.items()}
        cdef double tol = 1.0e-8,value
        cdef int i, index, code
        cdef str elemstr
        cdef Element elem
        print('File name is : ', filename)
        with open(filename,'w+') as fn:
            fn.write('from atpy import *\n')
            for index,elem in enumerate(self.pyelems):
                elemstr= '{0:<10} = {1:>10}({2:<10},'.format(elem.name, elem.elem_kind, f"'{elem.name}'")
                for i in range(KWD_NUM):
                    if fabs(elem.parms[i]) > tol:
                        value =  elem.parms[i]
                        elemstr = elemstr + f'{index2parms[i]:<5}={value:.6e}  ,'
                elemstr = elemstr[:-1] + ')\n' if elemstr[-1]==',' else elemstr + ')\n'
                fn.write(elemstr)
            clsname = self.__class__.__name__
            elemstr= '{0:<10} = {1:>10}({2:<10},'.format(self.name, clsname, f"'{self.name}'")
            for index,elem in enumerate(self.pyseq):
                elemstr = elemstr+ '{:<10},'.format(elem.name)
            elemstr = elemstr[:-1] + ')\n' if elemstr[-1]==',' else elemstr + ')\n'
            fn.write(elemstr)
        

