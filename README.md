# Structural-optimization-based-on-genetic-algorithm
# 基于遗传算法的人致振动优化 #
* ### 最优化调谐质量阻尼器(Tuned Mass Damper,TMD)参数设计
```
#结构实际总质量 233.17t
#一阶振型参与质量 186.536t
#一阶质量占比 80%
import math
#TMD质量比
mass_moda=1.6
#结构阻尼比
damping_ratio_structure=0.02
#结构一阶频率
frequency_structure = 2.41
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
```

<div align=center>
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215322917-744ecc35-fb79-4c79-a06d-11c333d31674.png" />
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215321260-7b035c60-78ed-47d2-8f28-d39df2b7decb.svg" />
</div>

* ### GA算法优化结果展示 ###
<div align=center>
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215321299-c6573987-1c64-4cdd-bfd5-95a27623a227.png" />
  <img height="180px" src="https://user-images.githubusercontent.com/98397090/215321295-a7a48243-a3d6-42aa-b668-b911d0bfb14f.png" />
</div>

