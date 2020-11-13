from ..cython.constants cimport*
from ..cython.ast cimport*
from libc.math cimport fabs,fmod,floor,pow,fdim,sqrt,fmax, fmin
from libc.stdlib cimport calloc

cdef dict token_enum
cdef dict enum_token

cdef class Token:
    cdef:
        int kind, data_kind
        str value
        int lineno, column

cdef class Lexer:
    cdef:
        int count
        int token_num
        int line_num
        int column
        list tokens
        dict elems_index
        
    
    
    cdef Token get_current_token(self)
    
    cdef get_next_token(self)
    
    cdef int check_next_token(self)
    
    cdef void tokenize(self, str code)

    
cdef class Parser:  # 定义语法分析器的类
    cdef:
        Lexer       lexer
        Token       current_token
        bint        isdatabase
        dict        elems_index
        double**    kwd_properties
        double**    tws_properties
        double**    loc_properties
        double*     glb_properties
        
    
    cdef void set_database(self, double** kwd_properties, double** tws_properties, double** loc_properties, double* glb_properties)nogil
        
    cdef void error(self,int input, int current)

    cdef void eat(self, int kind)

    cdef AST* property(self)
    
    cdef AST* function(self, int func)
    
    cdef AST* slice(self,int func, int start, int end, int data_kind, int index)
        
    cdef AST* factor(self)
    
    cdef AST* term(self)
    
    cdef AST* expr(self)

    cdef AST* parse(self, str code)

    
    
    
