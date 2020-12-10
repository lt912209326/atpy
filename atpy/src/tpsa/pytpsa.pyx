
#Initialize the DA domain.
def init_base(int da_order, int da_dim, int n_vec):
    print("Init the base with order,number or vars and number of DAVectors.")
    da_init(da_order, da_dim, n_vec)
    
def test_print(int index):
    da[index].print()
def test_add(int case):
    if case>0:
        pass
    else:
        return
    print("Fundamental calculations of DA vectors.")
    cdef DAVector x = 1 + da[0] + 2*da[1] + 5*da[2]
    x.print()
    
    if case>1:
        pass
    else:
        return
    cdef DAVector y = exp(x)
    y+=12.0
    y.print()
    if case>2:
        pass
    else:
        return
    
#     #Substitute a number for a base.
    print("Substitute a number for a base.")
    cdef DAVector z
    da_substitute_const(y, 0, 1, z)
    z.print()
    if case>3:
        pass
    else:
        return

    #Substitute a DA vector for a base.
    print("Substitute a DA vector for a base.")
    da_substitute(y, 0, x, z)
    z.print()
    if case>4:
        pass
    else:
        return

    #Substitute multiple DA vectors for bases at once.
    print("Substitute multiple DA vectors for bases at once.")
    cdef vector[DAVector] lv=vector[DAVector](2)
    lv[0] = sin(x)
    lv[1] = cos(x)
    
    if case>5:
        pass
    else:
        return

    cdef vector[unsigned int] idx=[0,1]

    da_substitute(y, idx, lv, z)
    z.print()

    if case>6:
        pass
    else:
        return
    #Bunch processing for substitutions.
    print("Bunch processing for substitutions.")
    cdef vector[DAVector] lx=vector[DAVector](3)
    cdef vector[DAVector] ly=vector[DAVector](3)
    lx[0] = x
    lx[1] = y
    lx[2] = sinh(x)

    if case>7:
        pass
    else:
        return
    da_substitute(lx, idx, lv, ly)
    ly[0].print()
    ly[1].print()
    ly[2].print()

    
    if case>8:
        pass
    else:
        return
    #Composition of DA vectors with numbers.
    print("Composition of DA vectors with numbers.")
    cdef vector[double] lm=[0.1, 2, 1]
    cdef vector[double] ln=vector[double](3)
    da_composition(lx, lm, ln)

#     # for(auto& x: ln) print(x<<' ')
#     # print(std::endl<<std::endl

    if case>9:
        pass
    else:
        return

    #Composition of DA vectors with DA vectors.
    print("Composition of DA vectors with DA vectors.")
    cdef vector[DAVector] lu=vector[DAVector](3)
    lu[0] = sin(x)
    lu[1] = cos(x)
    lu[2] = tan(x)
    
    
    if case>10:
        pass
    else:
        return
    da_composition(lx, lu, ly)
    ly[0].print()
    ly[1].print()
    ly[2].print()
