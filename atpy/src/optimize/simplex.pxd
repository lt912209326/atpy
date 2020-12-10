from libcpp.vector cimport vector
cdef extern from "./include/simplex.h" namespace "BT":
    vector[D] Simplex[D,OP](OP f, vector[D] init) except *
    vector[D] Simplex[D,OP](OP f, vector[D] init, D tol) except *
    vector[D] Simplex[D,OP](OP f, vector[D] init, D tol,int iterations) except *
    vector[D] Simplex[D,OP](OP f, vector[D] init, D tol,int iterations, vector[vector[D] ] x) except *