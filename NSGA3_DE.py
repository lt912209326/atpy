import geatpy as ea
from tqdm import tqdm
import numpy as np
import time
import warnings

class NSGA3_DE(ea.MoeaAlgorithm):
        
    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population) # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'NSGA3-DE'
        if self.problem.M < 10 and self.problem.NCV < 10:
            self.ndSort = ea.ndsortESS # 采用ENS_SS进行非支配排序
        else:
            self.ndSort = ea.ndsortTNS # 高维目标采用T_ENS进行非支配排序，速度一般会比ENS_SS要快
        self.selFunc = 'tour' # 基向量选择方式，采用锦标赛选择
        if population.Encoding == 'RI':
            self.mutOper = ea.Mutde(F = None) # 生成差分变异算子对象
            self.recOper = ea.Xovbd(XOVR = 0.8, Half = True) # 生成二项式分布交叉算子对象，这里的XOVR即为DE中的Cr
        else:
            raise RuntimeError('编码方式必须为''RI''.')
        self.F = 0.5 # 差分变异缩放因子（可以设置为一个数也可以设置为一个列数与种群规模数目相等的列向量）
        self.pc = 0.2 # 交叉概率
        self.CVdrawing = None
    
    
    
    
    def stat(self, population): # 分析记录，更新进化记录器，population为传入的种群对象
        feasible = np.where(np.all(population.CV <= 0, 1))[0] if population.CV is not None else np.arange(population.sizes) # 找到可行解个体的下标
        if len(feasible) > 0:
            tempPop = population[feasible] # 获取可行解个体
            self.pop_trace.append(tempPop) # 添加记录（只添加可行解个体到种群记录器中）
            self.forgetCount = 0 # “遗忘策略”计数器清零
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2:
                # 绘制目标空间动态图
                self.ax = ea.moeaplot(tempPop.ObjV, 'objective values', False, self.ax, self.currentGen, gridFlag = True)
            elif self.drawing == 3:
                # 绘制决策空间动态图
                self.ax = ea.varplot(tempPop.Phen, 'decision variables', False, self.ax, self.currentGen, gridFlag = False)
            self.timeSlot = time.time() # 更新时间戳
        else:
            self.currentGen -= 1 # 忽略这一代
            self.forgetCount += 1 # “遗忘策略”计数器加1
            if self.CVdrawing == 2:
                # 绘制目标空间动态图
                self.ax = ea.moeaplot(population.CV, 'CV values', False, self.ax, self.forgetCount, gridFlag = True)
            elif self.CVdrawing == 3:
                # 绘制决策空间动态图
                self.ax = ea.varplot(population.Phen, 'decision variables', False, self.ax, self.forgetCount, gridFlag = False)
        
    
    
    def reinsertion(self, population, offspring, NUM, uniformPoint):
        
        """
        描述:
            重插入个体产生新一代种群（采用父子合并选择的策略）。
            NUM为所需要保留到下一代的个体数目。
            
        """
        
        # 父子两代合并
        population = population + offspring
        # 选择个体保留到下一代
        [levels, criLevel] = self.ndSort(population.ObjV, NUM, None, population.CV, self.problem.maxormins) # 对NUM个个体进行非支配分层
        chooseFlag = ea.refselect(population.ObjV, levels, criLevel, NUM, uniformPoint, self.problem.maxormins) # 根据参考点的“入龛”个体筛选
        
        selected_flag = np.all([levels<float('inf'), chooseFlag],axis=0) # chooseFlag中True并且levels中被NUM个进行非支配分层则为True，确保选择策略在非支配占优和被转中个体
        selected_id = np.where(selected_flag)[0] #将进行selecting的部分在种群中的索引
        FitN = 1/levels[selected_id].reshape(len(selected_id),1)#**2 # 虚拟适应度 如：1/1^2 1/2^2, 1/3^2……
        NewChrId = ea.selecting('rws',FitN,NUM) #返回被选择个体在FitN的索引 0-len(FitN)
        
        rbest = selected_id[NewChrId] # best solution在种群中的索引
        
        return population[chooseFlag],population[rbest]
    
    def run(self, prophetPop = None): # prophetPop为先知种群（即包含先验知识的种群）
        #==========================初始化配置===========================
        population = self.population
        self.initialization() # 初始化算法模板的一些动态参数
        #===========================准备进化============================
        uniformPoint, NIND = ea.crtup(self.problem.M, population.sizes) # 生成在单位目标维度上均匀分布的参考点集
        population.initChrom(NIND) # 初始化种群染色体矩阵，此时种群规模将调整为uniformPoint点集的大小，initChrom函数会把种群规模给重置
        self.call_aimFunc(population) # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND] # 插入先知种群
        #===========================开始进化============================
        
        r0_index = np.arange(NIND)
        rbest_index=np.random.randint(NIND,size=NIND)
        rbest = population[rbest_index]
        while self.terminated(population) == False:
            # 进行差分进化操作
            offspring = population.copy() # 存储子代种群
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field, [offspring.Chrom,None,None,rbest.Chrom,offspring.Chrom]) # 变异
            tempPop = population + offspring # 当代种群个体与变异个体进行合并（为的是后面用于重组）
            offspring.Chrom = self.recOper.do(tempPop.Chrom) # 重组
            self.call_aimFunc(offspring) # 计算目标函数值
            # 重插入生成新一代种群
            population, rbest = self.reinsertion(population, offspring, NIND, uniformPoint)
        return self.finishing(population) # 调用finishing完成后续工作并返回结果