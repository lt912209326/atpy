
from ..cython.constants cimport*

from libc.math cimport fabs,fmod,floor,pow,fdim,sqrt,fmax, fmin
from libc.stdlib cimport calloc,free
import re

cdef enum:
    ADD
    SUB
    MUL
    DIV
    POW
    MOD
    FLOOR
    DOT
    SLICE
    COMMA
    NUMBER
    PROPERTY
    ABS
    SQRT
    DIM
    MAX
    MIN
    LPAREN
    RPAREN
    LBRA
    RBRA
    POSITION
    INDEX
    END
    

cdef cppclass AST:
    int token
    inline void __init__()nogil:
        pass
    inline double calc()nogil:
        pass



cdef cppclass Node(AST):
    AST* left
    AST* right
    inline void __init__(AST* left, int token, AST* right)nogil:
        this.left = left
        this.token = token
        this.right = right
    
    inline double calc()nogil:
        if this.token == ADD:
            return this.left.calc() + this.right.calc()
        elif this.token == SUB:
            return this.left.calc() - this.right.calc()
        elif this.token == MUL:
            return this.left.calc()*this.right.calc()
        elif this.token == DIV:
            return this.left.calc()/this.right.calc()
        elif this.token == POW:
            return this.left.calc()**this.right.calc()
        elif this.token == MOD:
            return fmod(this.left.calc(), this.right.calc() )
        elif this.token == FLOOR:
            return floor( this.left.calc()/this.right.calc() )



cdef cppclass Number(AST):
    double value
    inline void __init__(int token, double number)nogil:
        this.token = token
        this.value = number
    
    inline double calc()nogil:
        return this.value


cdef cppclass Property(AST):
    int  token, index, position
    double* value
    inline void __init__(int token, int position, int index, double* database)nogil:
        this.token = token
        this.position = position
        this.index = index
        this.value = database
    
    inline double calc()nogil:
        return this.value[0]
    
        


cdef cppclass MonoFunction(AST):
    AST* arg
    AST* right
    inline void __init__(int token, AST* arg)nogil:
        this.token = token
        this.arg = arg
            
    inline double calc()nogil:
        if this.token == ABS:
            return fabs(this.arg.calc() )
        elif this.token==SQRT:
            return sqrt(this.arg.calc())


cdef cppclass BiFunction(AST):
    AST* arg1
    AST* arg2
    __init__(AST* arg1, int token, AST* arg2)nogil:
        this.arg1 = arg1
        this.token = token
        this.arg2 = arg2
        
    inline double calc()nogil:
        if this.token==MIN:
            return fmin(this.arg1.calc(), this.arg2.calc() )
        elif this.token==MAX:
            return fmax(this.arg1.calc(), this.arg2.calc() )
        elif this.token==DIM:
            return fdim(this.arg1.calc(), this.arg2.calc() )
    
    
    
cdef cppclass Slice(AST):
    int     start, end, index
    double**    database
    __init__(int token, int start, int end, int index, double** database)nogil:
        cdef:
            int i
        this.start = start
        this.token = token
        this.end = end
        this.index=index
        this.database = <double**>calloc(this.end-this.start+1,sizeof(double*))
        for i in range(start,end+1):
            this.database[i-start] = &database[i][index]
    
    
    inline double calc()nogil:
        cdef:
            int i, total_properties= this.end-this.start+1
            double value
        value = this.database[0][0]
        if this.token==MIN:
            for i in range(1,total_properties):
                value = fmin(value,this.database[i][0])
        elif this.token==MAX:
            for i in range(1,total_properties):
                value = fmax(value,this.database[i][0])
        return value
    
    __dealloc__():
        free(this.database)

