
__all__=[]

CONSTRAINT0 = Constraint(NULL, 1.0, 0.0, 0.0, 0.0, 0)


cdef double get_constraint_value(Constraint* pconst)nogil:
    cdef double value
    value = pconst.expr.calc()
    
    if pconst.type == LOWER:
        return pconst.scale*fdim(pconst.lower, value)
    elif pconst.type == UPPER:
        return pconst.scale*fdim(value, pconst.upper)
    elif pconst.type == BOTH:
        return pconst.scale*fdim(value, pconst.upper) + pconst.scale*fdim(pconst.lower, value)
    elif pconst.type == EQUIV:
        return pconst.scale*fabs(value - pconst.equiv)

cdef double collect(Constraint* constr, int constr_num, double* fvect)nogil:
    cdef:
        int i
        double          value
        double          fsum = 0.0
    for i in range(constr_num):
        value = get_constraint_value( &constr[i] )
        fsum += value*value
        fvect[i] = value
        # print(value)
    return fsum




VARIABLE0  =Variable(NULL, -1.0e10, 1.0e10, NULL, -1, -1, -1, -1, NULL, NULL)


OPTIMIZE0=Optimize(NULL, 1)
