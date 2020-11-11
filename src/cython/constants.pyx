
KWD_INDEX = {'l': L, 'angle': ANGLE, 'k1': K1, 'k2': K2, 'k3': K3, 'k4': K4, 'e1': E1, 'e2': E2, 'tilt': TILT}

    #'marker': MARKER, 'drift': DRIFT, 'bend': BEND, 'quadrupole': QUADRUPOLE, 'sextupole': SEXTUPOLE, 'octuople': OCTUOPLE, 'decapole': DECAPOLE, 
    
    
    
TWS_INDEX = {'betax': BETAX, 'alphax': ALPHAX, 'gammax': GAMMAX, 'nux': NUX, 'betay': BETAY, 'alphay': ALPHAY, 'gammay': GAMMAY, 
             'nuy': NUY, 'etax': ETAX, 'etapx': ETAPX, 'chromx': CHROMX, 'chromy': CHROMY }

LOC_INDEX = {'R11': R11, 'R12': R12, 'R13': R13, 'R14': R14, 'R15': R15, 'R16': R16, 
            'R21': R21, 'R22': R22, 'R23': R23, 'R24': R24, 'R25': R25, 'R26': R26, 
            'R31': R31, 'R32': R32, 'R33': R33, 'R34': R34, 'R35': R35, 'R36': R36, 
            'R41': R41, 'R42': R42, 'R43': R43, 'R44': R44, 'R45': R45, 'R46': R46, 
            'R51': R51, 'R52': R52, 'R53': R53, 'R54': R54, 'R55': R55, 'R56': R56, 
            'R61': R61, 'R62': R62, 'R63': R63, 'R64': R64, 'R65': R65, 'R66': R66, 
            's': S, 'x': X, 'y': Y, 'z': Z, 'theta': THETA, 'phi': PHI, 'psi': PSI }

GLB_INDEX = {'mass0':MASS0, 'energy': ENERGY, 'gamma':GAMMA, 'emitx': EMITX, 'NAT_CHROMX': NAT_CHROMX, 'NAT_CHROMY': NAT_CHROMY, 'TOTAL_K2X': TOTAL_K2X, 
             'TOTAL_K2Y': TOTAL_K2Y, 'totalchromx': TOTALCHROMX, 'totalchromy': TOTALCHROMY, 'circumference': CIRCUMFERENCE, 
             'u0': U0, 'momentum_compact_factor': MOMENTUM_COMPACT_FACTOR, 'damp_factor': DAMP_FACTOR, 'energy_disp2': ENERGY_DISP2, 'spin': SPIN,
             'RI1': RI1, 'RI2': RI2, 'RI3': RI3, 'RI3A': RI3A, 'RI4': RI4, 'RI5': RI5, 'RI6': RI6, 'RI7': RI7 }


TOTAL_INDEX = {**KWD_INDEX,**TWS_INDEX,**LOC_INDEX,**GLB_INDEX}