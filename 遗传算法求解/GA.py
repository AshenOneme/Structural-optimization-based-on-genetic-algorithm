# -*- coding: UTF-8 -*-
import os
# 创建文件夹
def Create_folder(filename):
    filename = filename.strip()
    filename = filename.rstrip("\\")
    isExists = os.path.exists(filename)
    if not isExists:
        os.makedirs(filename)
        return  True
    else:
        return False

filename = "./Phen"
Create_folder(filename)
filename = "./ObjV"
Create_folder(filename)

import numpy as np
import geatpy as ea  # import geatpy

from bridge import MyProblem  # 导入自定义问题接口

# 实例化问题对象
problem = MyProblem()
# 构建算法

algorithm = ea.soea_EGA_templet(problem,
                                ea.Population(Encoding='BG', NIND=25),
                                MAXGEN=20,  # 最大进化代数
                                logTras=1)  # 表示每隔多少代记录一次日志信息

# 求解
ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True, dirname='result')
np.savetxt(f"./ObjV/out_ObjV{problem.rank}.txt",problem.out_ObjV)
np.savetxt(f"./Phen/out_Phen{problem.rank}.txt",problem.out_Phen)