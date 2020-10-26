# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示了一个带约束连续决策变量的最小化目标的双目标优化问题。
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self, line):
        self.line = line
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        Dim = len(self.line['variables'].index[0])#1 # 初始化Dim（决策变量维数）
        
        self.M = len(self.line['optima'].index[0])
        maxormins = self.line['optima'].index[4]  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        if self.line['cell']==True:
            self.NCV=len(self.line['constraints'].index[0])+1
        else:
            self.NCV=len(self.line['constraints'].index[0])
        
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = self.line['variables'].bounds[0]  # 决策变量下界
        ub = self.line['variables'].bounds[1]  # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, self.M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        self.scale = np.array([1]*self.NCV,dtype=float)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        nind = Vars.shape[0]
        if pop.ObjV is None:
            pop.ObjV = np.zeros((nind,self.M),dtype=float)
            pop.CV = np.zeros((nind, self.NCV),dtype=float)
        self.line.get_results(Vars, pop.ObjV,pop.CV,self.scale)
