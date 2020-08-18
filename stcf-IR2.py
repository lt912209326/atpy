from atpy import *
# ARCFOFO lattice definition
#---ELEMENTS DIFINE---#
MK0 = Marker(beta_x=0.04,alpha_x=0,gamma_x=25,beta_y=0.0008,alpha_y=0,gamma_y=1250)
L1 = Drift(l=0.935)
L2 = Drift(l=0.2)
QD1 = Quadrupole(l=0.3,k1=-5.734671860546432)
QF2 = Quadrupole(l=0.4,k1= 2.8060832378782052)
L3 = Drift(l=6)
L31 = Drift(l=4.5)
L32 = Drift(l=0.5)
QD3 = Quadrupole(l=0.4,k1=-1.6923916161556376)
L4 = Drift(l=0.22)
QF4 = Quadrupole(l=0.4,k1= 1.705134294804317)
L51 = Drift(l=0.5)
L52 = Drift(l=1.5-0.0)

QD5 = Quadrupole(l=0.4,k1=-2.1692589345478774)
L6 = Drift(l=0.2)
L61 = Drift(l=1.5)
L62 = Drift(l=1.0)
QF6 = Quadrupole(l=0.4,k1= -1.663370351531717+0.5)
L7 = Drift(l=1.0)
L71 = Drift(l=1.0)
L72 = Drift(l=1.0)


QD7 = Quadrupole(l=0.4,k1=-1.3)
L8 = Drift(l=0.8)
L81 = Drift(l=0.5+0.5)
L82 = Drift(l=0.5+0.5)
QF8 = Quadrupole(l=0.4,k1= 1.6)
L9 = Drift(l=1.0)

QD9 = Quadrupole(l=0.4,k1= -0.6)
L10= Drift(l=1.0)

B1 = Dipole(l=1,angle=0.05)
B2 = Dipole(l=1,angle=0.05)

#---CELL DEFINE---#
ARCFODO=Lattice(MK0,L1,QD1,L2,QF2,L31,B1,L32,QD3,L4,QF4,L51)#,L52)#,QD5,L6,QF6,L7,QD7,L81,L82,QF8,L9,QD9,L10,L10,QD9,L9,QF8,L82,L81)

#---OPTIMIZATION---#
var= (#{'place':0,'variables':{'beta_x':(0.03,0.045)}},#,'beta_y':(0.00075,0.00085)
#       {'place':1,'variables':{'l':(0.9,0.99)}},
      {'place':2,'variables':{'k1':(-9.2, -3.1)}},#
#       {'place':3,'variables':{'l':(0.1,1.0/2)}},
      {'place':4,'variables':{'k1':(2.01, 8)}},#
      {'place':5,'variables':{'l':(2.5,7.5)}},
      {'place':6+2,'variables':{'k1':(-4, 0)}},#
#       {'place':7+2,'variables':{'l':(0.1,1)}},
      {'place':8+2,'variables':{'k1':(0.01, 4)}},
      {'place':7+4,'variables':{'l':(0.1,2)}},
     )

CV= (
#     {'place':2 ,'constraints':{'beta_y':(200,2000)}},
#      {'place':7 ,'constraints':{'alpha_x':(1,600),'alpha_y':(1,600)} },
     {'place':11 ,'constraints':{'alpha_x':0,'alpha_y':0,"nu_x":0.5,"nu_y":0.5}},#'beta_x':(0,10),'beta_y':(0,10),
#      {'place':9+2 ,'constraints':{'alpha_x':0,'alpha_y':0} },
#      {'place':(0,9+2) ,'constraints':{'telescope_x':-15,'telescope_y':-15}},
#      {'place':9+2 ,'constraints':{'beta_x':(0,4),'beta_y':(0,0.8)}},
     
#      {'place':18 ,'constraints':{'alpha_x':0,'alpha_y':0,'beta_y':(0,0.01)} },
#      {'place':23 ,'constraints':{'alpha_x':0,'alpha_y':0,'beta_x':(0,0.01)} },
#      {'place':(5,11) ,'constraints':{'beta_x':(0,400),'beta_y':(0,450)}},
#      {'place':16 ,'constraints':{'beta_x':(0,30),'beta_y':(0,20)}},
#      ,'beta_x':(20,150),'alpha_x':0,'alpha_y':0,'eta_x':0,'
    )

opti = ({'place':4,'optima':{'beta_x':1}},
        {'place':4,'optima':{'beta_y':1}},
#         {'place':8,'optima':{'alpha_y':-1}},
#         {'place':8,'optima':{'alpha_x':-1}},
       )

variables = Variables(*var)
constraints = Constraints(*CV)
optima = Optima(*opti)

ARCFODO.set_variables(variables)
ARCFODO.set_constraints(constraints)
ARCFODO.set_optima(optima)









from myproblem import MyProblem as Problem
from main import main as main
from cProfile import run

problem = Problem(ARCFODO)
NDSet,pop_trace = main(problem,NIND=300,MAXGEN=200,drawing=1,Parallel=True)

# run('main(problem,NIND=300,MAXGEN=200,drawing=1,Parallel=True)',)





























