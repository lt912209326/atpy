# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 19:26:35 2020

@author: LT

This file defines all the constants used in this program
###############################################################################
Don't change it unless you know what you are doing!!!
###############################################################################
"""

cdef dict KWD_INDEX

cdef dict TWS_INDEX

cdef dict LOC_INDEX

cdef dict GLB_INDEX


cdef dict TOTAL_INDEX




#ELEMENTS CODE
cdef enum:
    MARKER = 0
    DRIFT = MARKER + 1
    DIPOLE = DRIFT + 1
    QUADRUPOLE = DIPOLE + 2
    SEXTUPOLE = QUADRUPOLE + 2
    OCTUOPLE = SEXTUPOLE + 2
    DECAPOLE = OCTUOPLE + 2

#KEYWORD CODE
cdef enum:
    L
#     K0 = L + 1
    ANGLE
    K1
    K2
    K3
    K4
    E1
    E2
    TILT
    KWD_NUM


cdef enum:
    KWD
    TWS
    LOC
    GLB


#TWISS CODE
cdef enum:
    BETAX
    ALPHAX
    GAMMAX
    NUX
    BETAY
    ALPHAY
    GAMMAY
    NUY
    ETAX
    ETAPX
    CHROMX
    CHROMY
    TWS_NUM

#GEOMETRY CODE
cdef enum:
    R11
    R12
    R13
    R14
    R15
    R16
    R21
    R22
    R23
    R24
    R25
    R26
    R31
    R32
    R33
    R34
    R35
    R36
    R41
    R42
    R43
    R44
    R45
    R46
    R51
    R52
    R53
    R54
    R55
    R56
    R61
    R62
    R63
    R64
    R65
    R66
    S  
    X
    Y
    Z
    THETA
    PHI
    PSI
    LOC_NUM
    
    


                    
#GLOBAL PARAMETER CODE
cdef enum:
    MASS0
    GAMMA
    ENERGY
    EMITX
    NAT_CHROMX
    NAT_CHROMY
    TOTAL_K2X
    TOTAL_K2Y
    TOTALCHROMX
    TOTALCHROMY
    CIRCUMFERENCE
    U0
    MOMENTUM_COMPACT_FACTOR
    DAMP_FACTOR
    ENERGY_DISP2
    SPIN
    RI1
    RI2
    RI3
    RI3A
    RI4
    RI5
    RI6
    RI7
    GLB_NUM
    
    

# =============================================================================
# CONSTRAINT CONST
# =============================================================================
cdef enum:
    ZERO
    ONE
    TWO
    THREE
    FOUR
    FIVE
    SIX
    SEVEN
    EIGHT
    NINE
    TEN

cdef enum:
    LOWER
    UPPER
    BOTH
    EQUIV