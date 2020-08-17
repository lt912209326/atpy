import geatpy as ea
from tqdm import tqdm

class NSGA3_DE(ea.moea_NSGA3_DE_templet):
    
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
        pbar = tqdm(total=self.MAXGEN)
        while self.terminated(population) == False:
            # 进行差分进化操作
            offspring = population.copy() # 存储子代种群
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field) # 变异
            tempPop = population + offspring # 当代种群个体与变异个体进行合并（为的是后面用于重组）
            offspring.Chrom = self.recOper.do(tempPop.Chrom) # 重组
            self.call_aimFunc(offspring) # 计算目标函数值
            # 重插入生成新一代种群
            population = self.reinsertion(population, offspring, NIND, uniformPoint)
            pbar.update(1)
        pbar.close()
        print()
            
        return self.finishing(population) # 调用finishing完成后续工作并返回结果
