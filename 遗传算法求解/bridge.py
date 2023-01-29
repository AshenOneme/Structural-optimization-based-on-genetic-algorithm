# -*- coding: UTF-8 -*-
import openseespy.opensees as ops
import vfo.vfo as vfo
import matplotlib.pyplot as plt
import math
import numpy as np
import geatpy as ea
from mpi4py import MPI

class MyProblem(ea.Problem):  # 继承Problem父类
    def __init__(self):
        name = 'MyProblem'  # 初始化name（函数名称，可以随意设置）
        M = 1  # 初始化M（目标维数）
        maxormins = [1]*M  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 3  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim  # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb =[0.01,0.01,0.01]    # 决策变量下界
        ub =[100,500,1]  # 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界
        ubin = [1] * Dim  # 决策变量上边界
        self.out_ObjV = np.array([])
        self.out_Phen = np.empty([0, Dim])
        self.comm = MPI.COMM_WORLD
        self.size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def aimFunc(self, pop):  # 目标函数
        x1 = pop.Phen[:, [0]] # 获取表现型矩阵的第一列，得到所有个体的x1的值
        x2 = pop.Phen[:, [1]]
        x3 = pop.Phen[:, [2]]
        accel_temp = np.array([])
        for s in range(len(x1)):
            ops.wipe()
            ops.model('basic', '-ndm', 3, '-ndf', 6)

            print("****************************************", "Material", "****************************************")
            IDSteel = 1
            Fy_Steel = 360
            E0_Steel = 210000
            bs_Steel = 0.005
            R0 = 12.5
            cR1 = 0.925
            cR2 = 0.15
            ops.uniaxialMaterial('Steel02', IDSteel, Fy_Steel, E0_Steel, bs_Steel, R0, cR1, cR2)

            IDPT = 2
            Fy_PT = 1860
            E0_PT = 210000
            bs_PT = 0.001
            R0 = 18
            cR1 = 0.925
            cR2 = 0.15
            a1 = 0
            a2 = 1
            a3 = 0
            a4 = 1
            sigInit = 188.71
            ops.uniaxialMaterial('Steel02', IDPT, Fy_PT, E0_PT, bs_PT, R0, cR1, cR2, a1, a2, a3, a4, sigInit)

            IDConcrete = 3
            fpc = -40
            epsc0 = -0.0024
            fpcu = -18
            epsU = -0.006
            ops.uniaxialMaterial('Concrete01', IDConcrete, fpc, epsc0, fpcu, epsU)

            IDPlate = 4
            ops.nDMaterial('ElasticIsotropic', IDPlate, 32500, 0.2, 0)

            IDFiber_Plate = 5
            ops.nDMaterial('PlateFiber', IDFiber_Plate, IDPlate)

            IDBearingZ = 6
            ops.uniaxialMaterial('Elastic', IDBearingZ, 1900000)

            IDBearingXY = 7
            ops.uniaxialMaterial('Elastic', IDBearingXY, 13640)

            IDXX = 8
            ops.uniaxialMaterial('ElasticPP', IDXX, 1000, 100000, -100000, 0)

            IDZZ = 9
            ops.uniaxialMaterial('ElasticPP', IDZZ, 1000e10, 2000000, -200000, 0)
            print("****************************************", "Section", "****************************************")
            Fiber_X_Section = 1
            ops.section('Fiber', Fiber_X_Section, '-GJ', 1e15)
            ops.patch('rect', IDSteel, 1, 4, -400, -150, -370, 150)
            ops.patch('rect', IDSteel, 4, 1, -370, -150, 370, -125)
            ops.patch('rect', IDSteel, 4, 1, -370, 125, 370, 150)
            ops.patch('rect', IDSteel, 1, 4, 370, -150, 400, 150)

            Fiber_Y_Section = 2
            HI = 368
            WI = 250
            cx = 16
            cy = 12
            A = [-HI / 2 - cx, -WI / 2]
            B = [-HI / 2, WI / 2]
            C = [-HI / 2, -cy / 2]
            D = [HI / 2, cy / 2]
            E = [HI / 2, -WI / 2]
            F = [HI / 2 + cx, WI / 2]
            ops.section('Fiber', Fiber_Y_Section, '-GJ', 1e15)
            ops.patch('rect', IDSteel, 1, 4, A[0], A[1], B[0], B[1])
            ops.patch('rect', IDSteel, 4, 1, C[0], C[1], D[0], D[1])
            ops.patch('rect', IDSteel, 1, 4, E[0], E[1], F[0], F[1])

            R90 = 3
            nc1, nr1 = 8, 1
            ri1, re1 = 90, 80
            a_beg1, a_end1 = 0., 360.
            ops.section('Fiber', R90, '-GJ', 1e15)
            ops.patch('circ', IDSteel, nc1, nr1, 0., 0., ri1, re1, a_beg1, a_end1)

            PlateFiber = 4
            ops.section('PlateFiber', PlateFiber, IDFiber_Plate, 120)

            print("****************************************", "Node", "****************************************")
            Node_Messages = np.loadtxt("Node_Messages.txt")
            Len_Node = Node_Messages.shape[0]
            for i in range(Len_Node):
                ops.node(i + 1, Node_Messages[i, 0], Node_Messages[i, 1], Node_Messages[i, 2])

            Node_Mass = np.loadtxt("Mass.txt")
            Len_Mass = Node_Mass.shape[0]
            for i in range(Len_Mass):
                ops.mass(i + 1, Node_Mass[i], Node_Mass[i], Node_Mass[i], 0, 0, 0)

            ops.node(101, 0, 0, -400)
            ops.node(125, 1152.1, -1600, -400)
            ops.node(143, 2304.21, -3200, -400)
            ops.node(166, 3456.32, -4800, -400)

            ops.fix(101, 1, 1, 1, 1, 1, 1)
            ops.fix(125, 1, 1, 1, 1, 1, 1)
            ops.fix(143, 1, 1, 1, 1, 1, 1)
            ops.fix(166, 1, 1, 1, 1, 1, 1)

            # =========================================================================================================
            ops.node(118, 40480.15, 0, 0)
            ops.node(142, 40951.1, -1600, 0)
            ops.node(159, 41422.05, -3200, 0)
            ops.node(182, 41892.99, -4800, 0)

            ops.fix(18, 0, 1, 0, 0, 0, 0)
            ops.fix(42, 0, 1, 0, 0, 0, 0)
            ops.fix(59, 0, 1, 0, 0, 0, 0)
            ops.fix(82, 0, 1, 0, 0, 0, 0)

            ops.fix(118, 1, 1, 1, 1, 1, 1)
            ops.fix(142, 1, 1, 1, 1, 1, 1)
            ops.fix(159, 1, 1, 1, 1, 1, 1)
            ops.fix(182, 1, 1, 1, 1, 1, 1)

            print("****************************************", "coordTransf", "****************************************")
            coordTransf = "PDelta"  # Linear, PDelta, Corotational
            IDXTransf = 1
            ops.geomTransf(coordTransf, IDXTransf, 0, 0, 1)
            IDYTransf = 2
            ops.geomTransf(coordTransf, IDYTransf, 0, 0, 1)
            IDRTransf = 3
            ops.geomTransf(coordTransf, IDRTransf, 1, 0, 0)

            print("****************************************", "beamIntegration",
                  "****************************************")
            IDXIntegration = 1
            ops.beamIntegration('Trapezoidal', IDXIntegration, Fiber_X_Section, 4)
            IDYIntegration = 2
            ops.beamIntegration('Trapezoidal', IDYIntegration, Fiber_Y_Section, 4)
            IDRIntegration = 3
            ops.beamIntegration('Trapezoidal', IDRIntegration, R90, 4)

            print("****************************************", "X", "****************************************")
            for i in range(17):
                ops.element('forceBeamColumn', i + 1, i + 1, i + 2, IDXTransf, IDXIntegration)

            ops.element('forceBeamColumn', 18, 19, 20, IDXTransf, IDXIntegration)
            ops.element('forceBeamColumn', 19, 20, 21, IDXTransf, IDXIntegration)
            ops.element('forceBeamColumn', 20, 22, 23, IDXTransf, IDXIntegration)
            ops.element('forceBeamColumn', 21, 23, 24, IDXTransf, IDXIntegration)

            for i in range(17):
                ops.element('forceBeamColumn', i + 22, i + 25, i + 26, IDXTransf, IDXIntegration)

            for i in range(16):
                ops.element('forceBeamColumn', i + 39, i + 43, i + 44, IDXTransf, IDXIntegration)

            ops.element('forceBeamColumn', 55, 60, 61, IDXTransf, IDXIntegration)
            ops.element('forceBeamColumn', 56, 61, 62, IDXTransf, IDXIntegration)
            ops.element('forceBeamColumn', 57, 63, 64, IDXTransf, IDXIntegration)
            ops.element('forceBeamColumn', 58, 64, 65, IDXTransf, IDXIntegration)

            for i in range(16):
                ops.element('forceBeamColumn', i + 59, i + 66, i + 67, IDXTransf, IDXIntegration)

            print("****************************************", "Y", "****************************************")
            ops.element('forceBeamColumn', 101, 1, 25, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 102, 25, 43, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 103, 43, 60, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 104, 60, 66, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 105, 2, 19, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 106, 19, 26, IDYTransf, IDYIntegration)

            ops.element('forceBeamColumn', 107, 3, 20, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 108, 20, 27, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 109, 27, 44, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 110, 44, 61, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 111, 61, 67, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 112, 4, 21, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 113, 21, 28, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 114, 28, 45, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 115, 45, 62, IDYTransf, IDYIntegration)
            ops.element('forceBeamColumn', 116, 62, 68, IDYTransf, IDYIntegration)

            for i in range(11):
                ops.element('forceBeamColumn', 117 + i * 3, 5 + i, 29 + i, IDYTransf, IDYIntegration)
                ops.element('forceBeamColumn', 118 + i * 3, 29 + i, 46 + i, IDYTransf, IDYIntegration)
                ops.element('forceBeamColumn', 119 + i * 3, 46 + i, 69 + i, IDYTransf, IDYIntegration)

            for i in range(3):
                ops.element('forceBeamColumn', 150 + i * 5, 16 + i, 22 + i, IDYTransf, IDYIntegration)
                ops.element('forceBeamColumn', 151 + i * 5, 22 + i, 40 + i, IDYTransf, IDYIntegration)
                ops.element('forceBeamColumn', 152 + i * 5, 40 + i, 57 + i, IDYTransf, IDYIntegration)
                ops.element('forceBeamColumn', 153 + i * 5, 57 + i, 63 + i, IDYTransf, IDYIntegration)
                ops.element('forceBeamColumn', 154 + i * 5, 63 + i, 80 + i, IDYTransf, IDYIntegration)

            print("****************************************", "Rod Piece", "****************************************")
            # for i in range(7):
            #     ops.element('forceBeamColumn', 201+4*i,83+2*i,4+2*i, IDRTransf, IDRIntegration)
            #     ops.element('forceBeamColumn', 202+4*i,83+2*i,28+2*i, IDRTransf, IDRIntegration)
            #     ops.element('forceBeamColumn', 203+4*i,84+2*i,45+2*i, IDRTransf, IDRIntegration)
            #     ops.element('forceBeamColumn', 204+4*i,84+2*i,68+2*i, IDRTransf, IDRIntegration)
            for i in range(7):
                ops.element('elasticBeamColumn', 201 + 4 * i, 83 + 2 * i, 4 + 2 * i, 5341, 206000e5, 0.79e5, 38720129,
                            19360065, 19360065, IDRTransf)
                ops.element('elasticBeamColumn', 202 + 4 * i, 83 + 2 * i, 28 + 2 * i, 5341, 206000e5, 0.79e5, 38720129,
                            19360065, 19360065, IDRTransf)
                ops.element('elasticBeamColumn', 203 + 4 * i, 84 + 2 * i, 45 + 2 * i, 5341, 206000e5, 0.79e5, 38720129,
                            19360065, 19360065, IDRTransf)
                ops.element('elasticBeamColumn', 204 + 4 * i, 84 + 2 * i, 68 + 2 * i, 5341, 206000e5, 0.79e5, 38720129,
                            19360065, 19360065, IDRTransf)

            print("****************************************", "PT", "****************************************")
            ops.element('truss', 229, 20, 83, 6359, IDPT)
            ops.element('truss', 230, 61, 84, 6359, IDPT)
            ops.element('truss', 231, 23, 95, 6359, IDPT)
            ops.element('truss', 232, 64, 96, 6359, IDPT)

            for i in range(6):
                ops.element('truss', 233 + 4 * i, 83 + 2 * i, 85 + 2 * i, 6358, IDPT)
                ops.element('truss', 234 + 4 * i, 84 + 2 * i, 86 + 2 * i, 6358, IDPT)

            print("****************************************", "ShellMITC4", "****************************************")
            ops.element('ShellMITC4', 301, 1, 2, 26, 25, PlateFiber)
            ops.element('ShellMITC4', 302, 2, 3, 20, 19, PlateFiber)
            ops.element('ShellMITC4', 303, 19, 20, 27, 26, PlateFiber)
            ops.element('ShellMITC4', 304, 25, 27, 44, 43, PlateFiber)
            ops.element('ShellMITC4', 305, 43, 44, 61, 60, PlateFiber)
            ops.element('ShellMITC4', 306, 60, 61, 67, 66, PlateFiber)
            ops.element('ShellMITC4', 307, 3, 4, 28, 27, PlateFiber)
            ops.element('ShellMITC4', 308, 27, 28, 45, 44, PlateFiber)
            ops.element('ShellMITC4', 309, 44, 45, 68, 67, PlateFiber)

            for i in range(12):
                ops.element('ShellMITC4', 310 + 3 * i, 4 + i, 5 + i, 29 + i, 28 + i, PlateFiber)
                ops.element('ShellMITC4', 311 + 3 * i, 28 + i, 29 + i, 46 + i, 45 + i, PlateFiber)
                ops.element('ShellMITC4', 312 + 3 * i, 45 + i, 46 + i, 69 + i, 68 + i, PlateFiber)

            for i in range(2):
                ops.element('ShellMITC4', 346 + 3 * i, 16 + i, 17 + i, 41 + i, 40 + i, PlateFiber)
                ops.element('ShellMITC4', 347 + 3 * i, 40 + i, 41 + i, 58 + i, 57 + i, PlateFiber)
                ops.element('ShellMITC4', 348 + 3 * i, 57 + i, 58 + i, 81 + i, 80 + i, PlateFiber)

            print("****************************************", "Bears", "****************************************")
            ops.element('twoNodeLink', 401, 101, 1, '-mat', IDBearingZ, IDBearingXY, IDBearingXY, '-dir', 1, 2, 3,
                        '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 402, 125, 25, '-mat', IDBearingZ, IDBearingXY, IDBearingXY, '-dir', 1, 2, 3,
                        '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 403, 143, 43, '-mat', IDBearingZ, IDBearingXY, IDBearingXY, '-dir', 1, 2, 3,
                        '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 404, 166, 66, '-mat', IDBearingZ, IDBearingXY, IDBearingXY, '-dir', 1, 2, 3,
                        '-orient', 0, 0, 1, 1, 0, 0)

            ops.element('twoNodeLink', 405, 18, 118, '-mat', IDXX, IDZZ, '-dir', 1, 3)
            ops.element('twoNodeLink', 406, 42, 142, '-mat', IDXX, IDZZ, '-dir', 1, 3)
            ops.element('twoNodeLink', 407, 59, 159, '-mat', IDXX, IDZZ, '-dir', 1, 3)
            ops.element('twoNodeLink', 408, 82, 182, '-mat', IDXX, IDZZ, '-dir', 1, 3)

            print("****************************************", "TMD", "****************************************")
            IDViscous = 10
            C_damper = x1[s,:][0]
            ops.uniaxialMaterial('Viscous', IDViscous, C_damper, 1)

            IDElastic = 11
            K_damper = x2[s,:][0]
            ops.uniaxialMaterial('Elastic', IDElastic, K_damper, 0.0, K_damper)

            M_damper = x3[s,:][0]

            ops.node(133, 19566.64, -1600, -500)
            ops.node(150, 19566.64, -3200, -500)
            ops.node(134, 22066.64, -1600, -500)
            ops.node(151, 22066.64, -3200, -500)
            ops.node(135, 24566.64, -1600, -500)
            ops.node(152, 24566.64, -3200, -500)

            ops.mass(133, 0, 0, M_damper, 0, 0, 0)
            ops.mass(150, 0, 0, M_damper, 0, 0, 0)
            ops.mass(134, 0, 0, M_damper, 0, 0, 0)
            ops.mass(151, 0, 0, M_damper, 0, 0, 0)
            ops.mass(135, 0, 0, M_damper, 0, 0, 0)
            ops.mass(152, 0, 0, M_damper, 0, 0, 0)

            ops.fix(133, 1, 1, 0, 1, 1, 1)
            ops.fix(150, 1, 1, 0, 1, 1, 1)
            ops.fix(134, 1, 1, 0, 1, 1, 1)
            ops.fix(151, 1, 1, 0, 1, 1, 1)
            ops.fix(135, 1, 1, 0, 1, 1, 1)
            ops.fix(152, 1, 1, 0, 1, 1, 1)

            ops.element('twoNodeLink', 501, 33, 133, '-mat', IDViscous, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 502, 33, 133, '-mat', IDElastic, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 503, 50, 150, '-mat', IDViscous, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 504, 50, 150, '-mat', IDElastic, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)

            ops.element('twoNodeLink', 505, 34, 134, '-mat', IDViscous, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 506, 34, 134, '-mat', IDElastic, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 507, 51, 151, '-mat', IDViscous, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 508, 51, 151, '-mat', IDElastic, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)

            ops.element('twoNodeLink', 509, 35, 135, '-mat', IDViscous, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 510, 35, 135, '-mat', IDElastic, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 511, 52, 152, '-mat', IDViscous, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)
            ops.element('twoNodeLink', 512, 52, 152, '-mat', IDElastic, '-dir', 1, '-orient', 0, 0, 1, 1, 0, 0)

            print("****************************************", "Model_View", "****************************************")
            # vfo.plot_model(show_nodes='yes', show_nodetags='yes', show_eletags='no', line_width=1)
            #
            # from ModeAnalyze import *
            # modal(10)
            print("****************************************", "Rayleigh", "****************************************")
            xDamper = 0.02
            nEigenI = 1
            nEigenJ = 2
            lambdaN = ops.eigen('-fullGenLapack', nEigenJ)
            lambdaI = lambdaN[nEigenJ - 2]
            lambdaJ = lambdaN[nEigenJ - 1]
            omegaI = math.pow(lambdaI, 0.5)
            omegaJ = math.pow(lambdaJ, 0.5)
            alphaM = xDamper * (2 * omegaI * omegaJ) / (omegaI + omegaJ)
            betaKcomm = 2.0 * xDamper / (omegaI + omegaJ)
            ops.rayleigh(alphaM, 0, 0, betaKcomm)

            print("****************************************", "Load", "****************************************")
            ops.timeSeries('Linear', 11)
            ops.pattern('Plain', 100, 11)
            Node_Mass = np.loadtxt("Mass.txt")
            Len_Mass = Node_Mass.shape[0]
            for i in range(Len_Mass):
                ops.load(i + 1, 0, 0, Node_Mass[i] * 10000, 0, 0, 0)
            ops.constraints("Penalty", 1e5, 1e5)
            ops.numberer("RCM")
            # ops.system('BandGen')
            ops.system('UmfPack')
            ops.test('NormDispIncr', 1.e-5, 2000)
            # ops.algorithm("KrylovNewton")
            ops.algorithm('NewtonLineSearch')
            ops.integrator("LoadControl", 0.1)
            ops.analysis("Static")
            ops.analyze(10)
            ops.loadConst("-time", 0.0)

            ops.recorder('Node', '-file', "Disp134.txt", "-time", '-node', 134, '-dof', 3, 'disp')
            ops.recorder('Node', '-file', "Disp34.txt", "-time", '-node', 34, '-dof', 3, 'disp')
            ops.recorder('Node', '-file', "ACC34.txt", "-time", '-node', 34, '-dof', 3, 'accel')
            ops.recorder('Node', '-file', "ACC51.txt", "-time", '-node', 51, '-dof', 3, 'accel')
            ops.wipeAnalysis()
            ops.timeSeries('Path', 22, '-dt', 0.01, '-filePath', 'Load.txt', '-factor',
                           2.5)
            ops.pattern("Plain", 200, 22)
            for i in range(82):
                ops.load(i + 1, 0, 0, 1, 0, 0, 0)
                ops.load(i + 1, 0, 0, 1, 0, 0, 0)
            ops.constraints("Penalty", 1e5, 1e5)
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test("NormDispIncr", 1.e-5, 1000)
            # ops.algorithm('KrylovNewton')
            # ops.algorithm('ExpressNewton', 2, 1.0, '-currentTangent', '-factorOnce')
            ops.algorithm('NewtonLineSearch',0.75)
            ops.integrator('Newmark', 0.5, 0.25)
            ops.analysis('Transient')
            ops.analyze(1000, 0.01)
            ops.wipe()
            accel = np.loadtxt("ACC34.txt", usecols=1)
            a = np.mean(np.abs(accel))
            accel_temp = np.append(accel_temp, a)
        pop.ObjV =accel_temp.reshape(-1, 1)
        out_temp = np.min(pop.ObjV)
        NO = np.argmin(pop.ObjV)
        self.out_ObjV = np.append(self.out_ObjV, out_temp)
        out_temp = pop.Phen[NO, :].reshape(1, 3)
        self.out_Phen = np.append(self.out_Phen, out_temp, axis=0)




