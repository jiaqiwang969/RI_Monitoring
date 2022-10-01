# RI_Monitoring
目的：建立以小波为特征的失稳监测模型

## 进程
- [x] 频谱横坐标等间距问题，通过fk换算回scale解决
- [x] 提取频带的能量作为特征：可以发现直接用归一化RI频带的PI特征，效果最佳
- [x] 不同频带：RI频带、1 BPF频带、联合指标
- [x] 不同位置：证明RI频带的与之前PI有一致的效果，但归一化的PI无效
- [x] 不同转速
- [x] 不同采样率

## 文档整理
- [ ] 分析内容整理



## 主代码介绍：
Main.mlx



## 数据处理的规律


### 任务规划：

- [ ] 1. 局部图：参考
https://blog.csdn.net/tengweitw/article/details/103555475

- [ ] 2. PI幅值的坐标修改（物理意义和单位明确）

- [ ] 3. 讨论16000、17000rpm时RI频带的能量幅值问题

- [ ] 4. 将PI随开度的变化曲线，加上一条失速和失稳的区分线,能量突变

- [ ] 5. 关于频带选择的分析和对比，可能需要做出3张图，包含重点传感器、不同转速的频带分布规律图

- [ ] 6. R1和R2主要的最佳传感器

- [ ] 7. 检查代码的统一处理问题

- [ ] 8. PI曲线画的太粗了

- [ ] 9. 对于不同实验的情况，统一阀门开度

- [ ] 10. 表述PI用能量的物理意义明确

- [ ] 11. 把所有的曲线都列出来，分析不同未知的RI频带曲线，最终将其分析集中在R1和R2。R1和R2的对失稳的灵敏度最高。位置和方法一起分析。额外，给个三个典型转速：9000rpm、12000rpm、16000rpm。

- [ ] 12. 采样率（可能要分析小波的coi）、转速，之前的分析

- [ ] 13. 统一基准(不要乱改)、研究选择问题

- [ ] 14. 模态的问题（12000rpm）专门有一个小节介绍这个问题，然后再下一章节对应到这个问题。发现了并总结这个现象

- [ ] 15. 在44%存在一个台阶的问题。有可能时旋转失稳和临近失速的状态变化，最好找到一个物理解释。研究该位置的能量的变化规律，总结






# 章节安排
## 1. PI指标的凝练（频带和传感器位置的选择）

### 如何提炼PI指标
在实验中，我们是通过减小压气机出口阀门开度，进而减小通道流量、增加背压的方式，来模拟制造压气机的失稳状态。在此过程中，将采集压气机信号，通过前述小波变换，提取出小波特征，累加失稳敏感频带得到能量指标，最终建立评估压气机运行状态的性能指数PI。

重点交代完整的实验步骤和过程：

*图1*
画一个PI提取的流程图：
<img src="https://cdn.mathpix.com/snip/images/6FR8nHTM5rxypqVrYjtsfk0-thj3NBXuBNqyELDacgg.original.fullsize.png" width="600px">
信号--在线-->2圈作为一个截断-小波变换->提取频带特征-->统一做归一化处理得到PI



### 频带的分析（不同频带的图）

- 选择合适的RI频带

如何选择合适的频带范围？

<img src="https://cdn.mathpix.com/snip/images/tSbnGFDIhBMvsvTcM5dGjaiKLB3u86Yo4GcrMVJTjJI.original.fullsize.png"  width="600px">

<img src="https://cdn.mathpix.com/snip/images/cWl5OODbKVSJOvWGyVSlWUzCmdJ7YEjMQcu5Un3I4kc.original.fullsize.png"  width="600px">

<img src="https://cdn.mathpix.com/snip/images/q4j1Yby_f5CDC8sXRU6Yu-dxY-AhWNyhIasSS3110KU.original.fullsize.png"  width="600px">


*图2-2*

画出9000rpm、12000rpm、16000rpm放在一块的PI对比图
说明能量的变化趋势、解释16000rpm的特点（同样符合频带选择的要求）



关于这个RI频带，早在2002年文献[1]，具有该特征的频谱的流动现象被称为“旋转失稳（Rotating Instability, RI）”。值得注意的是，虽然很多文献都提到该现象与压气机失稳密切关联[2-3]，但却没有很好地挖掘这一特征在失稳监测的应用。通过上述分析，我们发现RI频带的小波幅值，表现出连续且稳定（形态不变）的失稳跟随特性。抓住这一点，进一步引入HMM概率统计框架，可凝练成用于失稳监测的PI指标。



- 分析不同频带指标随阀门开度的变化过程（1BPF和2BPF、可以再将1阶次、2阶次、等频带一起对比）

	如图3-24给出了PI随阀门的变化曲线。压气机从阀门开度为64%左右进入失稳区，随着阀门开度的减小，其不稳定性递增，直到失速先兆SI（阀门开度为29%）增至最大。可以清晰看到，唯有基于RI频带小波特征的PI曲线，从L.1到L.3的幅值明显递增，整个曲线变化与压气机失稳演进的过程完全吻合，同时也证明了PI的灵敏性和有效性。相比而言，其余频带的PI曲线虽有所下降，但无法满足与失稳的一致性规律，这与图3-23的解释一致，同时也说明PI指标的表达更加简单有效。

### 传感器位置的分析（选择R1和R2的位置）

*图3*

<img src="https://cdn.mathpix.com/snip/images/iqCYBHuJxrj-tYrK1lhdUFHiuq71YcIBNVGDu7DBteo.original.fullsize.png" width="200px">

<img src="https://cdn.mathpix.com/snip/images/rGc_Zd11TFtyG5kwbrFhxb3JiGusEjJ4LHbe8S3pfiY.original.fullsize.png" width="200px">

<img src="https://cdn.mathpix.com/snip/images/DW47H9GOg9LmkLbmdxU6T6kd_MxLXxi7yE0DBhmPVF4.original.fullsize.png" width="200px">

- 不同位置的PI随阀门开度的变化（选择3个典型的转速）
图的布置：固定转速（3个典型转速）、不同位置的10个传感器


图3–25 给出了压气机从正常运行状态(阀门全开)到失速(阀门开度28%以下)的整个过程中，B1、R1-R8、C1等10个传感器的PI曲线。从图中可以清晰地看到，伴随着阀门开度的逐渐减小，所有传感器的PI值均有不同程度的变化，呈波动下降趋势。其中，R2传感器最为敏感，当阀门开度减小到64-65%左右时，其PI值已经先于所有其它传感器，出现具有一定显度的下降。对比R2压力传感器所处的位置和第3.2.2节中实验分析所得出叶顶压力泄漏涡流最先出现的位置(参照图3-7-(b))，可以发现两者完全吻合，表明，在进入失速失稳的临界状态时，R2最先感受到了压力扰动的变化。同时，也可以看到，R1传感器的PI曲线也有仅次于R2的较为明显趋势性变化。

随着阀门开度的进一步减小(压气机的失稳倾向进一步加剧)，当阀门开度减小到55%附近时，R1的PI曲线与R2形成交点。之后，伴随阀门开度减小(失稳加剧)，R1传感器的PI曲线呈现出迅速下降的特点。对照第3.2.2节实验分析、第3.2.3节中仿真分析的结论，我们已经知道，随着失稳的(增强)发展，与失稳密切关联的泄漏涡流，在叶顶前部(R2传感器位置)形成后，逐步朝向叶顶的前端(R1位置)移动、并不断扩大范围，其中心逐步移至R1位置。因此，可以认为，R1位置的压力传感器，很好地跟踪和反映了压气机失稳演进的过程，R1的PI值(曲线)非常突出和有效地指示了压气机的失稳变化过程。换句话说，叶顶泄漏涡流的压力云图(图3-7)浓缩成了一个的简单指标PI，可以替代斜排阵列测试的压力云图用作失稳监测。
另一方面，我们也看到，伴随着阀门开度的减小(失稳加强)，处于叶顶位置的其余R3-R8等传感器、以及位于叶顶后缘外面的C1传感器，其PI曲线均表现出一定的下降趋势性变化，但相较R1和R2，其灵敏度和突出性就逊色很多，不足以判断压气机失稳。而且，还可以看到，传感器的位置距离叶顶前缘越远，所监测到的压力扰动强度就越弱。事实上，这正好吻合了前面理论、仿真分析的泄漏涡流(即失稳)向叶顶前端发展的结果(参照图3-12)。

再来看看位于叶顶前缘外部的B1位置的传感器的情况。在失稳发展的过程中，B1传感器的PI曲线几乎没有任何明显变化，直到阀门开度降至大约30%(严重失稳状态)时，其PI曲线开始出现陡降。这再一次同实验分析的结果(见第3.2.1节的图3-3、第3.2.2节的图3-7-(e))相一致，即当压气机完全失稳之后，代表压力扰动的泄漏涡流会从叶顶前缘移出，正好到达B1传感器所在的区域。
依据上述分析结果，我们可以对传感器的安装位置和布局数量给出如下建议。(1)压力传感器应当安装在叶顶内靠近前缘的位置；(2)在条件允许的情况下，可以考虑安装两个传感器，例如，放在R1、R2位置。这样可以更早期地反映压力扰动，进而早期给出失稳预警；(3)若仅允许安装一个传感器，则可选择叶顶内靠近前缘的位置，例如R1传感器位置附近。

 

## 2. 采样率的影响
### 采样率的阶次分析
图3–26 给出了压气机从正常运行状态(阀门全开)到失速(阀门开度28%以下)的整个过程中，采样率从600点/圈到100点/圈(以50为间隔)的11条PI曲线。可以清晰地看到，伴随着采样率的降低， PI曲线所表现的下降趋势的梯度递减，直到100点/圈时，甚至完全不足以判断压气机失稳。

	根据香浓采样定理，低于两倍特征频带的采样率会导致信号混叠失真，这会直接影响PI结果的可靠性。当前所关心的RI频带范围为8-25阶次(f/f_"rot " ，参考图3-23)，这相当于采样25点/圈，所以理论上等价采样率至少要大于50点/圈。上述分析的等价采样率均满足该要求，但伴随着采样率的降低，PI表征失稳的能力会不断下降。这说明，信号所表现的非稳态随机特性，实际上无法仅通过采样定理的这一基础标准来进行衡量。从分析来看，尤其是针对等价采样率为100点/圈时，此时对应到每个叶道只有3-4个采样点（29个叶道）。这已经完全导致传感器无法捕捉到信号的峰值点，也无法反映转子的压力面和吸力面之间的压力梯度变化，造成关键的时频失稳特征严重失真，进而无法准确反映压气机的失稳过程。
	
	另一方面，我们也看到，随着采样率增加PI曲线的下降梯度增大，到500和600点/圈时，两者基本保持一致。可以发现，500点/圈已经对失稳监测具有明显的效果，已完全满足分析要求。倘若采样率进一步提高，不仅不会继续提高PI曲线针对失稳监测的灵敏度，还会增加分析计算的时间成本。
依据上述分析结果，我们可以对采样率的选择给出如下建议：(1)采样率设置在250-350点/圈已满足失速监测的要求；(2)在条件允许的情况下，考虑继续提高采样率到500点/圈，可以进一步提高监测的灵敏度。


## 3. 转速的影响




## 4. 12000rpm的特殊情况

## 5. 台阶的解释





