/*
 * Copyright(C) 2008 by Lingyun Yang
 * lingyun(.dot.]yang@gmail.com
 * http://www.lingyunyang.com
 *
 * Please get permission from Lingyun Yang before you redistribute this file.
 *
 */
//! \brief Automatic Differentiation
//! \file tpsa.h
//! \version $Id: tpsa.h,v 1.4 2009-04-17 17:32:23 frs Exp $
//! \author Lingyun Yang, http://www.lingyunyang.com/



#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#ifndef AD_HH
#define AD_HH
//! Type of order and number of variables.
//typedef unsigned char TNVND;
typedef unsigned int TNVND;
typedef unsigned int TVEC;
#ifdef __cplusplus
extern "C" {
#endif

 #ifdef MSVC_DLL
    _declspec(dllexport) void _stdcall ad_reserve(const TVEC* n);
    _declspec(dllexport) void _stdcall ad_init(const TNVND* nv, const TNVND* nd);
    _declspec(dllexport) void _stdcall ad_resetvars(const TNVND* nv);
    _declspec(dllexport) void _stdcall ad_alloc(TVEC* i);
    _declspec(dllexport) void _stdcall ad_free(const TVEC* i);
    _declspec(dllexport) void _stdcall ad_poolsize(size_t* n);

    _declspec(dllexport) void _stdcall ad_count(TVEC* n);
    _declspec(dllexport) void _stdcall ad_nvar(TVEC* n);
    _declspec(dllexport) void _stdcall ad_length(const TVEC* iv, unsigned int* n);
    _declspec(dllexport) void _stdcall ad_copy(const TVEC* i, const TVEC* j);
    _declspec(dllexport) void _stdcall ad_elem(const TVEC* ivec, unsigned int* idx, unsigned int* c, double* x);
    _declspec(dllexport) void _stdcall ad_pek(const TVEC* ivec, int* c, size_t* n, double* x);
    _declspec(dllexport) void _stdcall ad_pok(const TVEC* ivec, int* c, size_t* n, double* x);
    _declspec(dllexport) void _stdcall ad_var(const TVEC* ii, const double* x, unsigned int* iv);
    _declspec(dllexport) void _stdcall ad_abs(const TVEC* iv, double* r);
    _declspec(dllexport) void _stdcall ad_truncate(const TVEC* iv, const TNVND* d);

    _declspec(dllexport) void _stdcall ad_clean(const TVEC* iv, const double* eps);
    _declspec(dllexport) void _stdcall ad_reset(const TVEC* iv);
    _declspec(dllexport) void _stdcall ad_const(const TVEC* ii, const double* r);
    _declspec(dllexport) void _stdcall ad_fill_ran(const TVEC* iv, const double* ratio, const double* xm);

    _declspec(dllexport) void _stdcall ad_add(const TVEC* i, const TVEC* j);
    _declspec(dllexport) void _stdcall ad_add_const(const TVEC* i, double *r);

    _declspec(dllexport) void _stdcall ad_sub(const TVEC* i, const TVEC* j);
    //! internal multiplication, dst should be different from lhs and rhs.
    _declspec(dllexport) void _stdcall ad_mult(const TVEC* ivlhs, const TVEC* ivrhs, TVEC* ivdst);
    _declspec(dllexport) void _stdcall ad_mult_const(const TVEC* iv, double* c);

    _declspec(dllexport) void _stdcall ad_div_c(const TVEC* iv, const double* c);
    _declspec(dllexport) void _stdcall ad_c_div(const TVEC* iv, const double* c, TVEC* ivret);
    _declspec(dllexport) void _stdcall ad_div(const TVEC* ilhs, const TVEC* irhs, TVEC* idst);

    _declspec(dllexport) void _stdcall ad_sqrt(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_exp(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_log(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_sin(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_cos(const TVEC* iv, const TVEC* iret);

    _declspec(dllexport) void _stdcall ad_derivative(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_tra(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_shift(const TVEC* iv, unsigned int* ishift, const TVEC* iret, const double* eps);

    _declspec(dllexport) void _stdcall ad_read_block(const TVEC* iv, double* v, TNVND* J, const unsigned int* N);
    _declspec(dllexport) void _stdcall ad_save_block(const TVEC* iv, const double* v, const TNVND* J, const unsigned int* N);

    _declspec(dllexport) void _stdcall ad_rev(const TVEC* iv);

    _declspec(dllexport) void _stdcall ad_subst(const TVEC* iv, const TVEC* ibv, const TNVND* nbv, const TVEC* iret);
    //void ad_inverse(const TVEC* iv, const TNVND* nbv, const TVEC* iret, const TNVND* nret);
    _declspec(dllexport) void _stdcall ad_inverse(const TVEC* iva, const TNVND* nva, const TVEC* iret, const TNVND* nret);

    _declspec(dllexport) void _stdcall print_index(std::ostream& os);
    _declspec(dllexport) void _stdcall ad_print(const TVEC* iv);
    _declspec(dllexport) void _stdcall ad_print_array(const TVEC* iv, const TVEC* nv);
 #else
    void ad_reserve(const TVEC* n);
    void ad_init(const TNVND* nv, const TNVND* nd);
    void ad_resetvars(const TNVND* nv);
    void ad_alloc(TVEC* i);
    void ad_free(const TVEC* i);
    void ad_poolsize(size_t* n);

    void ad_count(TVEC* n);
    void ad_nvar(TVEC* n);
    void ad_length(const TVEC* iv, unsigned int* n);
    void ad_copy(const TVEC* i, const TVEC* j);
    void ad_elem(const TVEC* ivec, unsigned int* idx, unsigned int* c, double* x);
    void ad_pek(const TVEC* ivec, int* c, size_t* n, double* x);
    void ad_pok(const TVEC* ivec, int* c, size_t* n, double* x);
    void ad_var(const TVEC* ii, const double* x, unsigned int* iv);
    void ad_abs(const TVEC* iv, double* r);
    void ad_truncate(const TVEC* iv, const TNVND* d);

    void ad_clean(const TVEC* iv, const double* eps);
    void ad_reset(const TVEC* iv);
    void ad_const(const TVEC* ii, const double* r);
    void ad_fill_ran(const TVEC* iv, const double* ratio, const double* xm);

    void ad_add(const TVEC* i, const TVEC* j);
    void ad_add_const(const TVEC* i, double *r);

    void ad_sub(const TVEC* i, const TVEC* j);
    //! internal multiplication, dst should be different from lhs and rhs.
    void ad_mult(const TVEC* ivlhs, const TVEC* ivrhs, TVEC* ivdst);
    void ad_mult_const(const TVEC* iv, double* c);

    void ad_div_c(const TVEC* iv, const double* c);
    void ad_c_div(const TVEC* iv, const double* c, TVEC* ivret);
    void ad_div(const TVEC* ilhs, const TVEC* irhs, TVEC* idst);

    void ad_sqrt(const TVEC* iv, const TVEC* iret);
    void ad_exp(const TVEC* iv, const TVEC* iret);
    void ad_log(const TVEC* iv, const TVEC* iret);
    void ad_sin(const TVEC* iv, const TVEC* iret);
    void ad_cos(const TVEC* iv, const TVEC* iret);

    void ad_derivative(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    void ad_tra(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    void ad_shift(const TVEC* iv, unsigned int* ishift, const TVEC* iret, const double* eps);

    void ad_read_block(const TVEC* iv, double* v, TNVND* J, const unsigned int* N);
    void ad_save_block(const TVEC* iv, const double* v, const TNVND* J, const unsigned int* N);

    void ad_rev(const TVEC* iv);

    void ad_subst(const TVEC* iv, const TVEC* ibv, const TNVND* nbv, const TVEC* iret);
    void ad_inverse(const TVEC* iva, const TNVND* nva, const TVEC* iret, const TNVND* nret);

    void print_index(std::ostream& os);
    void ad_print(const TVEC* iv);
    void ad_print_array(const TVEC* iv, const TVEC* nv);
 #endif

#ifdef __cplusplus
}
#endif

#endif
