
from ..cython.constants cimport*

from libc.math cimport fabs,fmod,floor,pow,fdim,sqrt,fmax, fmin
from libc.stdlib cimport calloc
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
    NAME
    INDEX
    END
    

cdef cppclass AST:
    int token
    inline void __init__()nogil:
        pass
    inline double calc()nogil:
        pass
    int item(int code)nogil:
        pass
    double** database()nogil:
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

cdef cppclass Parameter(AST):
    int  index,parm,data_kind
    double* value
    double** properties
    inline void __init__(int index, int parm, int data_kind, double** kwd_properties, double** tws_properties, double** loc_properties, double* glb_properties,)nogil:
        this.token = PROPERTY
        this.index = index
        this.data_kind = data_kind
        if data_kind==GLB:
            this.value = &glb_properties[parm]
            this.parm  = parm
        elif data_kind==LOC:
            this.value = &loc_properties[index][parm]
            this.properties = loc_properties
            this.parm  = parm
        elif data_kind==TWS:
            this.value = &tws_properties[index][parm]
            this.properties = tws_properties
            this.parm  = parm
        elif data_kind==KWD:
            this.value = &kwd_properties[index][parm]
            this.properties = kwd_properties
            this.parm  = parm
    
    inline double calc()nogil:
        return this.value[0]
    
    inline int item(int code)nogil:
        if code ==0:
            return this.index
        elif code==1:
            return this.parm
        elif code==2:
            return this.data_kind
        
    inline double** database()nogil:
        return this.properties
        


cdef cppclass Function(AST):
    int parms_num
    double** parameters
    AST* left
    AST* right
    inline void __init__(int token, AST* left)nogil:
        this.left = left
        this.token = token
        this.parms_num = 0
        if left.token== SLICE:
            this.parms_num = left.item(0)
            this.parameters = left.database()
            
    inline double calc()nogil:
        cdef: 
            int i
            double minval, maxval
        if this.token == ABS:
            return fabs(this.left.calc() )
        elif this.token==MIN:
            minval = this.parameters[0][0]
            for i in range(1,this.parms_num):
                minval = fmin(minval,this.parameters[i][0])
            return minval
        elif this.token==MAX:
            maxval = this.parameters[0][0]
            for i in range(1,this.parms_num):
                maxval = fmax(maxval,this.parameters[i][0])
            return maxval
        elif this.token==SQRT:
            return sqrt(this.left.calc())
    
    inline double** database()nogil:
        return this.parameters
    
    inline int item(int code)nogil:
        if this.token==MAX or this.token==MIN:
            if code==0:
                return this.parms_num
            if code==1:
                return this.left.item(1)
        
cdef cppclass Dim(AST):
    AST* left
    AST* right
    __init__(AST* left, int token, AST* right)nogil:
        this.left = left
        this.token = token
        this.right = right
    inline double calc()nogil:
        return fdim( this.left.calc(), this.right.calc() )
        
        
    
cdef cppclass Slice(AST):
    AST* left
    AST* right
    inline void __init__(AST* left, int token, AST* right):
        this.left = left
        this.token = token
        this.right = right
    
    inline int item(int code)nogil:
        if code == 0:
            return this.right.item(0) - this.left.item(0) + 1
        elif code==1:
            return this.token
    
    inline double** database()nogil:
        cdef: 
            int i, start, end, parm
            double** tmp_ptr = left.database()
            double** parameters
        
        start = this.left.item(0)
        end = this.right.item(0)
        parm = this.left.item(1)
        parameters = <double**>calloc(end-start+1,sizeof(double*))
        for i in range(start,end+1):
            parameters[i-start] = &tmp_ptr[i][parm]
        return parameters
