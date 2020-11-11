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
    
    
        
    def initialization(self):
        
        """
        描述: 该函数用于在进化前对算法模板的参数进行初始化操作。
        该函数需要在执行算法模板的run()方法的一开始被调用，同时开始计时，
        以确保所有这些参数能够被正确初始化。
        
        """
        self.ax = None # 重置ax
        self.passTime = 0 # 初始化计时器
        self.forgetCount = 0 # 初始化“遗忘策略”计数器
        self.maxForgetCount = 3000 # 初始化“遗忘策略”计数器最大上限值，当超过这个上限时将终止进化
        self.pop_trace = [] # 初始化种群记录器
        self.currentGen = 0 # 设置初始为第0代
        self.evalsNum = 0 # 设置评价次数为0
        self.timeSlot = time.time() # 开始计时
    
    
    def call_aimFunc(self, pop):
        
        """
        使用注意:
        本函数调用的目标函数形如：aimFunc(pop), (在自定义问题类中实现)。
        其中pop为种群类的对象，代表一个种群，
        pop对象的Phen属性（即种群染色体的表现型）等价于种群所有个体的决策变量组成的矩阵，
        该函数根据该Phen计算得到种群所有个体的目标函数值组成的矩阵，并将其赋值给pop对象的ObjV属性。
        若有约束条件，则在计算违反约束程度矩阵CV后赋值给pop对象的CV属性（详见Geatpy数据结构）。
        该函数不返回任何的返回值，求得的目标函数值保存在种群对象的ObjV属性中，
                              违反约束程度矩阵保存在种群对象的CV属性中。
        例如：population为一个种群对象，则调用call_aimFunc(population)即可完成目标函数值的计算。
             之后可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。
        若不符合上述规范，则请修改算法模板或自定义新算法模板。
        
        """
        
        pop.Phen = pop.decoding() # 染色体解码
        if self.problem is None:
            raise RuntimeError('error: problem has not been initialized. (算法模板中的问题对象未被初始化。)')
        self.problem.aimFunc(pop,self.currentGen,self.MAXGEN) # 调用问题类的aimFunc()
        self.evalsNum = self.evalsNum + pop.sizes if self.evalsNum is not None else pop.sizes # 更新评价次数
        if type(pop.ObjV) != np.ndarray or pop.ObjV.ndim != 2 or pop.ObjV.shape[0] != pop.sizes or pop.ObjV.shape[1] != self.problem.M:
            raise RuntimeError('error: ObjV is illegal. (目标函数值矩阵ObjV的数据格式不合法，请检查目标函数的计算。)')
        if pop.CV is not None:
            if type(pop.CV) != np.ndarray or pop.CV.ndim != 2 or pop.CV.shape[0] != pop.sizes:
                raise RuntimeError('error: CV is illegal. (违反约束程度矩阵CV的数据格式不合法，请检查CV的计算。)')
    
    
    def stat(self, population): # 分析记录，更新进化记录器，population为传入的种群对象
        feasible = np.where(np.all(population.CV <= 1.0e-20, 1))[0] if population.CV is not None else np.arange(population.sizes) # 找到可行解个体的下标
        if len(feasible) > 0:
            tempPop = population[feasible] # 获取可行解个体
            self.pop_trace.append(tempPop) # 添加记录（只添加可行解个体到种群记录器中）
            self.forgetCount = 0 # “遗忘策略”计数器清零
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2 and self.currentGen%5==0:
                # 绘制目标空间动态图
                self.ax = ea.moeaplot(tempPop.ObjV, 'objective values', False, self.ax, self.currentGen, gridFlag = True)
            elif self.drawing == 3:
                # 绘制决策空间动态图
                self.ax = ea.varplot(tempPop.Phen, 'decision variables', False, self.ax, self.currentGen, gridFlag = False)
            self.timeSlot = time.time() # 更新时间戳
        else:
            self.currentGen -= 1 # 忽略这一代
            self.forgetCount += 1 # “遗忘策略”计数器加1
            if self.CVdrawing == 2 and self.forgetCount%5==0:
                # 绘制目标空间动态图
                self.ax = ea.moeaplot(population.CV, 'CV values', False, self.ax, self.forgetCount, gridFlag = True)
            elif self.CVdrawing == 3 and self.forgetCount%5==0:
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
#         CV=population.CV.T
#         prior=np.array([17,18,19,20],dtype=int)
#         print(0,0)
#         rank = self.ndSort(CV,NUM,prior, CV.shape[0], CV.shape[1],len(prior) )
#         levels = rank[:,0]
#         print( len( np.where( rank[:,0]>0 )[0] ) )
#         Inf = float('inf')
#         if len( np.where( rank[:,0]>0 )[0] )<NUM:
#             criLevel=np.max( rank[:,0] )+1
#             levels=np.array([i if i>0 else float('inf') for i in levels])#float('inf')
#             less = NUM-len( np.where( levels<Inf )[0] )
#             levels[ np.where(levels == Inf)[0][:less] ]=criLevel
#             print(1,criLevel,len( np.where( levels<Inf )[0] ))
#         else:
#             criLevel=np.max( rank[:,0] )
#             levels=np.array([i if i>0 else float('inf') for i in levels])#float('inf')
#             more = len( np.where( levels<Inf )[0] )-NUM
#             levels[ np.where(levels == criLevel)[0][:more] ]=Inf
#             print(2,criLevel,len( np.where( levels<Inf )[0] ))
#         levels = np.array(rank[:,0])
        [levels, criLevel] = self.ndSort(population.ObjV, NUM, None, population.CV, self.problem.maxormins) # 对NUM个个体进行非支配分层
#         print(np.max(levels[np.where( levels<float('inf') )[0] ] )  )
#         print( len( np.where( levels <float('inf') )[0] ) )
#         print(3,3)
        chooseFlag = ea.refselect(population.ObjV, levels, criLevel, NUM, uniformPoint, self.problem.maxormins) # 根据参考点的“入龛”个体筛选
#         #
#         print(4,4)
#         print(np.max(levels[np.where( levels<float('inf') )[0] ] )  )
        selected_flag = np.all([levels< (criLevel+1), chooseFlag],axis=0) # chooseFlag中True并且levels中被NUM个进行非支配分层则为True，确保选择策略在非支配占优和被转中个体
        selected_id = np.where(selected_flag)[0] #将进行selecting的部分在种群中的索引
        FitN = 1/np.sqrt(levels[selected_id].reshape(len(selected_id),1))  # 虚拟适应度 如：1/1^2 1/2^2, 1/3^2……
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
        
        #lbounds=self.problem.ranges[0]
        #ubounds=self.problem.ranges[1]
        #kNIND = int(NIND/20)
        #r0 = np.arange(NIND)
        while self.terminated(population) == False:
            # 进行差分进化操作
            #rrand = np.random.uniform(lbounds,ubounds,size=(kNIND,self.problem.Dim) )
            #np.random.shuffle(r0)
            #r0 = r0[:NIND-kNIND]
            
            offspring = population.copy() # 存储子代种群
            
            #rrand = np.vstack(( offspring.Chrom[r0], rrand ))
            
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field, [offspring.Chrom,None,None,rbest.Chrom,offspring.Chrom]) # 变异
            tempPop = population + offspring # 当代种群个体与变异个体进行合并（为的是后面用于重组）
            offspring.Chrom = self.recOper.do(tempPop.Chrom) # 重组
            self.call_aimFunc(offspring) # 计算目标函数值
            # 重插入生成新一代种群
            population, rbest = self.reinsertion(population, offspring, NIND, uniformPoint)
#             if self.currentGen%100==0 and self.currentGen>=200: 
#                 try:
#                     print('Please press Ctrl + V to stop in 4s and get Optimization!')
#                     for i in range(20):
#                         time.sleep(0.2)
#                 except KeyboardInterrupt:
#                     break
# #                 ifstop = raw_input("是否停止进化?:[True/False]")
# #                 if ifstop:
# #                     break
        try:
            return self.finishing(population) # 调用finishing完成后续工作并返回结果
        except RuntimeError:
            print('CV are not less than zero!')
            [levels, criLevel] = ea.ndsortDED(population.ObjV, None, 1, population.CV, self.problem.maxormins) # 非支配分层
            NDSet = population[np.where(levels == 1)[0]] # 只保留种群中的非支配个体，形成一个非支配种群
            return NDSet
    
    
    def finishing(self, population):
        
        """
        进化完成后调用的函数。
        
        """
        
        # 得到非支配种群
        [levels, criLevel] = ea.ndsortDED(population.ObjV, None, 1, population.CV, self.problem.maxormins) # 非支配分层
        NDSet = population[np.where(levels == 1)[0]] # 只保留种群中的非支配个体，形成一个非支配种群
        if NDSet.CV is not None: # CV不为None说明有设置约束条件
            NDSet = NDSet[np.where(np.all(NDSet.CV <= 1.0e-20, 1))[0]] # 最后要彻底排除非可行解
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        if NDSet.sizes == 0:
            raise RuntimeError('error: No feasible solution. (没找到可行解。)')
        # 绘图
        if self.drawing != 0:
            if NDSet.ObjV.shape[1] == 2 or NDSet.ObjV.shape[1] == 3:
                ea.moeaplot(NDSet.ObjV, 'Pareto Front', saveFlag = True, gridFlag = True)
            else:
                ea.moeaplot(NDSet.ObjV, 'Value Path', saveFlag = True, gridFlag = False)
        # 返回帕累托最优集
        return NDSet