from libc.math cimport fabs
from ..cython.constants cimport*
from ..cython.ast cimport*

ctypedef unsigned int   uint



cdef struct Constraint:
    AST*        expr
    double      scale
    double      lower
    double      upper
    double      equiv
    int         type        #1:lower 2:upper 3: lower&upper 4:equiv

cdef Constraint CONSTRAINT0


cdef double get_constraint_value(Constraint* pconst)nogil

cdef double collect(Constraint* constr, int constr_num, double* fvect)nogil





cdef struct Variable:
    double*     vary
    double      lower
    double      upper
    AST*        expr
    int         type        #LOWER,UPPER,BOTH, EQUIV
    int         data_kind
    int         position
    int         index
    Variable*   prev
    Variable*   next

cdef Variable   VARIABLE0


cdef struct Optimize:
    AST*        expr
    double      minormax

cdef Optimize   OPTIMIZE0