# -*- coding: utf-8 -*-

import geatpy as ea
from NSGA3_DE import NSGA3_DE
import time

def main(problem,population=None, NIND=200,MAXGEN=200, drawing=1, CVdrawing=None, F=None, CR=0.8, Parallel=False):
    '''
    problem                     自定义问题类
    NIND = 200                  # 种群规模
    MAXGEN=200                  最大进化代数
    drawing=0                   不作帕累托前沿图
            1                   作目标空间结果图
            2                   目标空间动态图
            3                   决策空间图
    
    
    '''
    """==================================种群设置==============================="""
    """================================实例化问题对象==========================="""
    if population is None:
        Encoding = 'RI'           # 编码方式
        Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
        population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """=================================算法参数设置============================"""
    myAlgorithm = NSGA3_DE(problem, population) # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = MAXGEN   # 最大进化代数
    myAlgorithm.drawing= drawing
    myAlgorithm.CVdrawing= CVdrawing
    myAlgorithm.mutOper.F= F
    myAlgorithm.mutOper.Parallel= Parallel
    myAlgorithm.recOper.XOVR=CR
    myAlgorithm.recOper.Parallel= Parallel
    
    """============================调用算法模板进行种群进化======================"""
    NDSet=myAlgorithm.run()
    NDSet.save()              # 把结果保存到文件中
    print('非支配个体数：%s 个'%(NDSet.sizes))
    print('单位时间找到帕累托前沿点个数：%s 个'%(int(NDSet.sizes // myAlgorithm.passTime)))
    print('用时：%s 秒'%(myAlgorithm.passTime))
    return NDSet,myAlgorithm.pop_trace