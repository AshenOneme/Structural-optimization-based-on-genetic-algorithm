# Structural-optimization-based-on-genetic-algorithm
# 基于遗传算法的人致振动优化 #

<div align=center>
  <img height="250px" src="https://user-images.githubusercontent.com/98397090/215322917-744ecc35-fb79-4c79-a06d-11c333d31674.png" />
</div>

* ### 最优化调谐质量阻尼器(Tuned Mass Damper,TMD)参数设计
```diff
#结构实际总质量 222.624t
#一阶振型参与质量 182.552t
#一阶质量占比 82%
import math
#一阶模态质量
mass_moda=182.552
#结构阻尼比
damping_ratio_structure=0.02
#结构一阶频率
frequency_structure = 1.758
ω=2*math.pi*frequency_structure

def tmd(ratio):
    mass_ratio_TMD=ratio/100
    mass_total_TMD=mass_ratio_TMD*mass_moda
    δ=(1-damping_ratio_structure*(math.sqrt(mass_ratio_TMD/(1+mass_ratio_TMD))))/(1+mass_ratio_TMD)
    damping_ratio_TMD=damping_ratio_structure/(1+mass_ratio_TMD)+math.sqrt(mass_ratio_TMD/(1+mass_ratio_TMD))
    ω_TMD=δ*ω
    K_TMD=mass_total_TMD*math.pow(ω_TMD,2)
    C_TMD=2*mass_total_TMD* ω_TMD*damping_ratio_TMD
    return K_TMD,C_TMD,mass_total_TMD
---> K=109.95 N/mm, C=1.81 N/(mm/s), M=0.91 ton
```

_结构模态信息_
<div align=center>
    
|_MODE NO_|__FREQUENCY__|__X-COMPONENT__|__Y-COMPONENT__|__Z-COMPONENT__|
|:---:|:---:|:---:|:---:|:---:|
|1|1.758|0.215%|0.017%|82.005%|
|MODAL MASS| |X:0.479|Y:0.038|Z:182.552|
|2|2.564|99.513%|0.000%|0.172%|
|MODAL MASS| |X:221.527|Y:0.000|Z:0.382|
|3|3.679|0.010%|82.141%|0.024%|
|MODAL MASS| |X:0.022|Y:182.856|Z:0.054|
    
</div>

* ### 人行荷载计算
```
#单人移动荷载
fp=2
G=750
Alpha1=0.5
Alpha2=0.2
Alpha3=0.1
f1=1.69
f2=2.41
f3=2.99
t=np.arange(0,3.01,0.01)
Fp=G+G*(Alpha1*np.sin(2*np.pi*1*fp*t)+Alpha2*np.sin(2*np.pi*2*fp*t-np.pi/2)+Alpha3*np.sin(2*np.pi*3*fp*t-np.pi/2))
```

<div align=center>
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215421119-97556e0e-394a-45e9-b519-b3bdfbc83521.png" />
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215420032-ddf91613-c362-4a32-ab27-794237389084.png" />
</div>



<div align=center>
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215321260-7b035c60-78ed-47d2-8f28-d39df2b7decb.svg" />
</div>

* ### GA算法优化结果展示 ###
<div align=center>

<img height="180px" src="https://user-images.githubusercontent.com/98397090/215423884-1d9e3a94-7ab4-4a6f-82e1-af1736f2a63a.png" />
<img height="180px" src="https://user-images.githubusercontent.com/98397090/215424145-294a33a9-3688-4ad5-9f25-f192ba299aff.png" />

</div>



<div align=center>
(a) 迭代进化20轮
</div>

