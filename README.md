# Structural-optimization-based-on-genetic-algorithm
# 基于遗传算法的人致振动优化 #
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

```

<div align=center>
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215421119-97556e0e-394a-45e9-b519-b3bdfbc83521.png" />
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215420032-ddf91613-c362-4a32-ab27-794237389084.png" />
</div>



<div align=center>
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215322917-744ecc35-fb79-4c79-a06d-11c333d31674.png" />
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215321260-7b035c60-78ed-47d2-8f28-d39df2b7decb.svg" />
</div>

* ### GA算法优化结果展示 ###
<div align=center>

<img height="180px" src="https://user-images.githubusercontent.com/98397090/215321299-c6573987-1c64-4cdd-bfd5-95a27623a227.png" />
<img height="180px" src="https://user-images.githubusercontent.com/98397090/215321295-a7a48243-a3d6-42aa-b668-b911d0bfb14f.png" />

</div>
<div align=center>
(a) 迭代进化20轮
</div>

