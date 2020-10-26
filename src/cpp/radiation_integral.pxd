from libc.math cimport sqrt, cos, sin,tan, fabs,sinh,cosh,pi
cdef inline void calc_sync_int(double rhoinv,double blen,double k1,double e1,double e2,double betxi,double alfxi,double dxi,double dpxi,double* I)nogil:
    '''
        !----------------------------------------------------------------------*
        !     Purpose:                                                         *
        !     Calculate synchrotron radiation integrals contribution of        *
        !     single element with parameters passed as input.                  *
        !                                                                      *
        !     Input:                                                           *
        !     rhoinv (double) inverse radius of curvature                      *
        !     blen (double) length of element                                  *
        !     k1 (double) gradient of element                                  *
        !     e1, e2 (double) pole face rotations at entrance and exit         *
        !     betxi, alfxi, dxi, dpxi (double) twiss parameters in x plane     *
        !     Output:                                                          *
        !     I[8] synchrotron radiation integral                              *
        !                                                                      *
        !     Author: Ghislain Roy - June 2014                                 *
        !----------------------------------------------------------------------*
    '''

    # local variables
    cdef double  dx2, gamx, dispaverage, curlyhaverage, lq
    cdef double  betx, alfx, dx, dpx, u0x, u1x, u2x
    cdef double  gammai, betxaverage, k1n
    cdef double k2, k, kl
    cdef double zero, one, two

    betx = betxi
    dx = dxi
    k1 = 0
    zero, one, two = 0.0, 1.0, 2.0

    # effect of poleface rotation
    alfx = alfxi - betxi*rhoinv*tan(e1)
    dpx = dpxi + dxi*rhoinv*tan(e1)

    gamx = (1+alfx**2)/betx

    # global gradient combining weak focusing and dipole gradient
    # k2 can be positive or negative and k can be real or imaginary
    k2 = rhoinv*rhoinv + 2*k1
    k = sqrt(k2)
    kl = k*blen

    # propagation of dispersion at exit
    dx2 = dx*cos(kl) + dpx*sin(kl)/k + rhoinv*(1-cos(kl))/(k*k)

    dispaverage = dx * sin(kl)/kl \
             + dpx * (1 - cos(kl))/(k*kl) \
             + rhoinv * (kl - sin(kl))/(k2*kl)

    curlyhaverage =  gamx*dx*dx + 2*alfx*dx*dpx + betx*dpx*dpx \
                + 2*rhoinv*blen*( -(gamx*dx + alfx*dpx)*(kl-sin(kl))/(kl*kl*k) \
                                 + (alfx*dx + betx*dpx)*(1-cos(kl))/(kl*kl)) \
                + blen*blen*rhoinv*rhoinv*( \
                       gamx*(3*kl - 4*sin(kl) + sin(kl)*cos(kl))/(2*k2*kl**3) \
                     - alfx*(1-cos(kl))**2/(k*kl**3) \
                     + betx*(kl-cos(kl)*sin(kl))/(2*kl**3))

    if (rhoinv != 0.0):
        I[1-1] += dispaverage * rhoinv * blen
        I[2-1] += rhoinv*rhoinv * blen
        I[3-1] += fabs(rhoinv)**3 * blen
        I[7-1] += rhoinv**3*blen
        I[4-1] += dispaverage*rhoinv*(rhoinv**2 + 2*k1) * blen \
               - rhoinv*rhoinv*(dx*tan(e1) + dx2*tan(e2))
        I[5-1] += curlyhaverage * fabs(rhoinv)**3 * blen

    if(k1 != 0.0):
        lq = blen
        gammai = (one+alfxi*alfxi)/betxi;
        dx2 = dx*cos(kl) + dpx*sin(kl)/k
        dispaverage = dx * sin(kl)/kl \
                 + dpxi * (1 - cos(kl))/(k*kl)
        if(k1 > zero):
            k1n = k1
            u0x = (one + sin(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/two
            u1x = sin(sqrt(k1n)*lq)**two/(k1n*lq)
            u2x = (one - sin(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/(two*k1n)
            dx2 = cos(sqrt(k1n)*lq)*dxi + (one/sqrt(k1n))*sin(sqrt(k1n)*lq)*dpxi
            dispaverage = (dxi+dx2)/two
        else:
            k1n = -k1
            u0x = (one + sinh(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/two
            u1x = sinh(sqrt(k1n)*lq)**two/(k1n*lq)
            u2x = -(one - sinh(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/(two*k1n)
            dx2 = cosh(sqrt(k1n)*lq)*dxi + (one/sqrt(k1n))*sinh(sqrt(k1n)*lq)*dpxi
            dispaverage = (dxi+dx2)/two

        betxaverage = betxi*u0x - alfxi*u1x  + gammai*u2x

        I[6-1] += (k1n**2)*betxaverage*lq
        I[8-1] += (k1n**2)*dispaverage**2*lq

cdef inline double calc_chrom(double l, double K1, double K2, double beta0, double alpha0, double eta0, double etap0)nogil:
    cdef double chrom = 0.0
    cdef double int_beta = 0
    cdef double gamma0 = (1+alpha0*alpha0)/beta0
    cdef double k1
    if K1>0:
        k1 = sqrt(K1)
        chrom = -0.0625*( (beta0*k1 -gamma0/k1 )*sin(2*k1*l) + 2*alpha0*cos(2*k1*l) + 2*(beta0*K1 + gamma0)*l - 2*alpha0 )/pi
    elif K1<0:
        k1 = sqrt(-K1)
        chrom = -0.0625*( (-beta0*k1 -gamma0/k1 )*sinh(2*k1*l) + 2*alpha0*cosh(2*k1*l) + 2*(beta0*K1 + gamma0)*l - 2*alpha0 )/pi
    if K2 !=0.0:
        chrom += 0.25*K2*(etap0*gamma0*l**3 + (eta0*gamma0 -2*etap0*alpha0 )*l*l + (etap0*beta0 -2*eta0*alpha0 )*l + eta0*beta0 )/pi
    return chrom
        
        
        