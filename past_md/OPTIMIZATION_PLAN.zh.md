# ECG1D_CPP 优化方案与问题诊断

## 1. 文档目的

这份文档总结当前 `ecg1d_cpp` 代码在两体基态 benchmark 和后续 many-body dynamical localization (MBDL) 研究上的主要问题，并给出一套分阶段、可执行的优化路线。

目标不是只把某个 benchmark 从 `1e-3` 调到 `1e-6`，而是回答两个更根本的问题：

1. 当前误差主要来自哪里？
2. 这套程序要如何演化成真正能够研究论文
   `Observation of many-body dynamical localization` 中实验现象内在机制的工具？

---

## 2. 当前代码的真实状态

### 2.1 现在擅长什么

当前程序已经具备以下基础模块：

- ECG/Gaussian 变分基底
- 重叠、动能、谐振子势、delta 接触相互作用、Gaussian 相互作用、`cos(2k_L z)` 势能项的矩阵元
- 对应的一阶/二阶导数与 metric tensor `C`
- 基于 `C dz = -g` 的 imaginary-time 风格 TDVP 更新
- 一个两体专用的 SVM 基底搜索器

这说明它已经是一个“静态变分优化器”的原型。

### 2.2 现在不擅长什么

当前程序**还不是**论文目标所需要的“周期驱动多体实时间传播器”：

- `kicking` 目前只是静态项 `cos(2k_L z)` 的期望值
- TDVP 采用“每一步要求能量下降”的线搜索
- 这等价于虚时间松弛/能量最小化，不是实时间演化
- 还没有真正实现周期脉冲 `V_kick(t) ~ \sum_n delta_{T_p}(t-nT)`
- 还没有实现 Floquet map
- 还没有实现与论文核心现象直接对应的观测量：
  - `n(k,t)`
  - 吸收能量随 kick 数的演化
  - entropy / von Neumann entropy
  - 一阶关联函数 `G^(1)(z,t)`

所以从方法论上讲，当前程序离论文目标还有一个明显的架构层差距。

---

## 3. 已观察到的数值现象

我在本地直接运行了现有可执行程序，得到以下结果。

### 3.1 一体谐振子

命令：

```bash
./build/ecg1d --tdvp-only
```

结果：

- 最终能量：`0.5000021995`
- 相对精确值 `0.5` 的误差：`2.1995e-06`

说明：

- 即使在一体谐振子这种最容易的问题上，当前 TDVP 也没有快速、干净地收敛到精确值
- 它最后进入了一个很慢的数值平台
- 这表明优化器本身就存在效率问题

### 3.2 两体 delta 相互作用

命令：

```bash
./build/ecg1d --delta
```

结果摘要：

- SVM 阶段：`E = 1.3098796772`
- TDVP 后：`E = 1.309750842`
- 目标值：`1.3067455`
- 最终误差：`3.0053e-03`

说明：

- `delta` benchmark 误差停在 `1e-3`
- TDVP 确实还能降能量，但只能降一点
- 到后期 `|dz|` 已很小，但能量平台仍然明显高于精确值
- 这不像单纯是步长问题，更像是“变分空间 + 优化策略”一起卡住了

### 3.3 两体 Gaussian 相互作用

命令：

```bash
./build/ecg1d --gaussian
```

结果摘要：

- SVM 阶段：`E = 1.5267309252`
- TDVP 后：`E = 1.526707895`
- 目标值：`1.526699831`
- 最终误差：`8.0640e-06`

说明：

- 对平滑相互作用，当前基底和优化器表现明显好于 `delta`
- 这非常像“普通 Gaussian 基底对平滑问题够用，但对 cusp/contact 问题表达力不足”

### 3.4 Python reference 的状态

本地运行：

```bash
python3 ../ecg1d/test_correctness.py
```

结果：

- `210 PASSED, 0 FAILED`

说明：

- 从“公式层面的导数、自洽性、有限差分检验”看，Python reference 当前是自洽的
- 因此当前 C++ 的主要问题不太像“明显公式写错”
- 更像是：
  - 变分流形太窄
  - 优化路径太弱
  - 目标物理与当前程序的架构不匹配

---

## 4. 主要问题归因

## 4.1 问题一：当前 `kicking` 不是论文里的真正 kick

论文中的哈密顿量核心是：

\[
H(t)=\sum_i \left[\frac{p_i^2}{2m} + \hbar \kappa \cos(2k_L z_i)\sum_n \delta_{T_p}(t-nT) + V_{ext}(z_i)\right]
+ g_{1D}\sum_{i<j}\delta(z_i-z_j)
\]

而当前代码里的 `kicking` 是静态的：

- 只在 Hamiltonian functional 中加入 `cos(2k_L z)` 项
- 没有时间依赖
- 没有 pulse duration
- 没有 period
- 没有 free evolution / kick evolution 的分裂

这意味着：

- 当前 `--kicking` 不能研究 dynamical localization
- 更不能研究 many-body dynamical localization
- 现在最多只是研究“含有余弦外势的静态基态”

这和论文的物理问题不是一回事。

## 4.2 问题二：TDVP 目前本质上是虚时间能量最小化

当前 TDVP 的关键逻辑是：

- 解一个最小二乘方向
- 试探更新
- 如果新能量不下降，就不断减小步长
- 直到能量下降为止

这适合做：

- 基态搜索
- 变分参数松弛

但不适合做：

- 实时间量子动力学
- 周期驱动
- 长时间 Floquet 演化

因为对实时间演化来说：

- 能量不一定单调下降
- kick 期间和 free evolution 期间的相位积累至关重要
- 需要保持相位、干涉、规范、归一化与长期稳定性

如果继续沿用现在这套“必须降能量”的更新规则，就无法回答 MBDL 的核心问题。

## 4.3 问题三：`delta` 问题对当前 Gaussian 基底不友好

从 benchmark 的结果可以看到：

- Gaussian interaction：误差约 `8e-6`
- delta interaction：误差约 `3e-3`

这非常说明问题。

原因是：

- `delta` 接触相互作用对应波函数在相对坐标上的 cusp / 非解析结构
- 普通 Gaussian 基底本质上是解析且光滑的
- 用少量普通 Gaussian 去逼近 cusp，通常收敛很慢

所以在两体 `delta` 问题上，如果继续使用当前这套普通 Gaussian 参数化：

- 想靠 TDVP 多走几步把 `1e-3` 压到 `1e-6`
- 往往是不现实的

这不是因为 TDVP 完全没用，而是因为表示能力先卡住了。

## 4.4 问题四：当前优化流程把 `u` 基本冻结了

当前流程是：

1. SVM 构建基底
2. 解一次广义本征值，设置 `u`
3. TDVP 只更新 `A/B/R`
4. 最后再用新的 basis 重新算一次能量

问题在于：

- 一旦非线性参数 `A/B/R` 变化，最优线性系数 `u` 也应该同步变化
- 当前做法相当于在一个“过时的线性组合”上优化非线性参数

结果就是：

- 真正应该走的下降方向被扭曲
- TDVP 可能在一个次优流形上收敛
- 尤其在多基函数情况下，这会很明显

更合理的做法是：

- 每次非线性参数变动后，都重新解一次广义本征值问题
- 把最低本征值 `E0(A,B,R)` 作为真正的目标函数

## 4.5 问题五：C++ 版优化器比 Python reference 更保守

当前 C++ 版有几个会明显影响收敛的平台因素：

- stochastic refinement 在 driver 里被关闭了
- `C` 矩阵没有加小正则
- SVD 截断固定为 `rcond = 1e-4`
- 步长增长机制比 Python reference 弱

这些设置会导致：

- 降能方向被过早截断
- 小奇异值方向完全不走
- 后期很容易卡在平台

因此即使公式都对，优化流程本身也会把精度限制在 `1e-6` 或 `1e-3` 平台附近。

## 4.6 问题六：当前 SVM 是“两体专用”的，不可直接扩展到 `N~10`

当前实现里最关键的限制是：

- driver 直接把 `N=2` 硬编码
- 随机基底生成器就是 `random_basis_2particle`
- 参数化默认依赖两体 CM/relative 宽度结构

所以现在的问题不是“怎么优化到 10 粒子”

而是：

- 主程序路径目前连一般 `N` 的 variational ansatz 都没有准备好

如果不先重构这层，后面谈 `N=10` 只会把误差和耗时同时放大。

## 4.7 问题七：显式 `N!` 置换求和会把大粒子数彻底卡死

当前很多模块都显式遍历全部 permutation：

- overlap
- Hamiltonian functional
- 各类导数
- metric tensor `C`

这意味着复杂度里有非常危险的 `N!` 因子。

对小体系还可以，
但对 `N=10`：

- `10! = 3628800`

这还只是 permutation 数，还没乘上：

- basis pair 数量
- 参数个数
- `C` 矩阵装配
- 每一步的线性代数代价

所以如果不改变对称化策略，`N~10` 在当前架构下几乎不可行。

## 4.8 问题八：当前程序没有实现论文关键观测量

论文真正关心的是：

- 动量分布 `n(k)` 是否冻结
- 能量是否饱和
- entropy 是否饱和
- 一阶关联函数 `G^(1)(z)` 的衰减形式是否发生变化

但当前代码还没有：

- momentum-space observable
- 熵
- Jensen-Shannon divergence
- `G^(1)(z)`

因此就算静态基态做得再准，也还不能直接支撑论文目标。

---

## 5. 为什么当前误差不能简单理解为“代码 bug”

当前更合理的判断是：

### 5.1 不像单点显性 bug

原因：

- Python reference 通过了完整 correctness test
- C++ 结果的误差模式有明显物理结构：
  - 平滑相互作用更准
  - 接触相互作用更差
  - 一体系统还能收敛到很小误差

这更像“方法限制”，而不是“某一行公式错了”

### 5.2 但仍然可能存在隐藏数值/实现问题

以下地方仍值得重点检查：

- `u` 不同步更新导致的伪平台
- `C` 的病态性与截断阈值
- complex Hermitian 结构被强行投成实矩阵的步骤
- basis 函数过于相似时的线性相关问题
- 大步长和线搜索失败后导致的早停

所以更准确的说法是：

- 当前主因不像“明显 bug”
- 但实现层仍然存在若干会放大误差的平台因素

---

## 6. 优化路线总览

建议分成四个阶段推进，而不是一次性冲击论文全部目标。

### 阶段 A：把“静态 benchmark 平台”做扎实

目标：

- 先把两体静态问题彻底吃透
- 明确误差究竟来自表示问题、优化问题还是实现问题

建议任务：

1. 把 C++ driver 改回完整优化流程
2. 每次 `A/B/R` 更新后重新求最优 `u`
3. 给 `C` 矩阵加小正则项
4. 把 SVD 阈值改成可配置
5. 恢复 stochastic refinement
6. 对 delta benchmark 做系统扫描：
   - basis size
   - regularization
   - `rcond`
   - 初始步长
   - 是否每步重解 `u`

阶段 A 的验收标准：

- 一体谐振子误差稳定到机器精度附近
- 两体 Gaussian 误差稳定压到 `1e-7 ~ 1e-8`
- 两体 delta 明确知道：是优化卡住，还是基底表达不够

### 阶段 B：专门解决 `delta` 的表示能力问题

目标：

- 不再用普通 Gaussian 硬逼 cusp

建议任务：

1. 把两体问题改写到 CM + relative coordinate
2. CM 方向直接固定在谐振子精确基态
3. 只对 relative coordinate 做变分
4. 为 relative basis 增加 cusp-aware 参数化

可选路线：

- 路线 B1：保留 Gaussian，但使用更强的 relative-coordinate 多尺度采样
- 路线 B2：在 relative coordinate 上加入 cusp-corrected basis
- 路线 B3：先只把 `delta` 当作静态 benchmark 工具，不要求它直接作为后续 many-body 主 ansatz

我更推荐：

- 先做 B1
- 如果 B1 仍明显卡在 `1e-4` 或 `1e-5`
- 再做 B2

阶段 B 的验收标准：

- 两体 delta benchmark 至少稳定进入 `1e-6` 量级
- 并且误差随 basis 扩展有清晰、可解释的收敛趋势

### 阶段 C：把程序从“静态松弛器”升级成“真正的 kicked dynamics 工具”

目标：

- 对齐论文里的物理问题

需要新增的核心能力：

1. 实时间 TDVP，而不是只做 imaginary-time
2. 显式时间依赖哈密顿量
3. pulse 序列：
   - free evolution
   - kick evolution
   - 或分裂算符/Floquet one-period propagator
4. 长时间稳定传播
5. 复相位结构保真

更具体地说，应把演化拆成：

- trap + interaction + kinetic 的 free evolution
- 短脉冲 kick 的驱动阶段

这一步完成后，程序才能真正回答：

- 为什么 `n(k)` 不继续扩散？
- 为什么能量和 entropy 会饱和？
- 为什么强相互作用与非相互作用的 `G^(1)` 衰减不同？

### 阶段 D：从小体系物理验证走向论文问题

目标：

- 做出能解释 MBDL 机制的最小可信结果

建议优先顺序：

1. `N=2`，有 trap + kick + interaction
2. `N=3`
3. `N=4`
4. 再考虑更大粒子数

重点不是盲目追求 `N=10`，而是先建立以下现象是否稳定出现：

- `n(k,t)` 冻结
- 吸收能量饱和
- entropy 饱和
- `G^(1)(z)` 衰减型变化

如果这些在 `N=2~4` 就已经有清晰趋势，那已经是非常有价值的理论结果。

---

## 7. 对当前代码最重要的具体修改建议

下面是我认为最值得优先做的改动。

## 7.1 第一优先级：重构静态优化器

### 建议 1：把目标函数改成“非线性参数上的最低本征值”

即：

- 变 `A/B/R`
- 每次都重新构造 `H,S`
- 重新解广义本征值
- 用最低本征值作为优化目标

优点：

- 避免 `u` 冻住带来的伪平台
- 物理上也更自然

### 建议 2：恢复 stochastic refinement

当前它在 driver 里虽然留了接口，但实际上没用起来。

建议：

- 在 SVM 之后先做若干轮 refine
- 再进入 TDVP

这通常对多峰/多尺度问题很重要。

### 建议 3：给 `C` 矩阵加正则，并把阈值参数化

建议最少引入：

- `lambda_C = 1e-8, 1e-7, 1e-6, 1e-5` 扫描
- `rcond` 扫描

需要实际记录：

- smallest singular value
- effective rank
- line search failure 次数
- 每一步能量下降量

### 建议 4：每一步记录优化诊断

至少输出：

- `E`
- `dE`
- `|dz|`
- `cond(C)`
- effective rank
- overlap matrix 最小特征值

否则很难区分“收敛”与“病态早停”。

## 7.2 第二优先级：做两体 delta 的专门基底

### 建议 5：切到 CM/relative coordinate

两体问题本身就允许：

- `R_cm = (x1+x2)/2`
- `r = x1-x2`

对谐振子 + 两体相互作用，这样做非常自然。

优点：

- 直接把问题拆成“容易部分 + 难部分”
- 相互作用只作用在相对坐标上
- basis 设计会清晰很多

### 建议 6：relative coordinate 用多尺度宽度采样

建议 relative Gaussian width 做 log-uniform 多尺度分布：

- 很窄的 Gaussian：负责 cusp 和短程结构
- 中等宽度：负责中程
- 宽的 Gaussian：负责尾部

当前两体随机基底虽然已经有一点这个味道，但还不够针对 `delta`。

## 7.3 第三优先级：重构为真正的 real-time kick 程序

### 建议 7：把 kick 变成显式时间依赖

至少要引入：

- kick period `T`
- pulse duration `Tp`
- kick number `N_p`
- time-dependent evolution loop

可以从最小实现开始：

- free evolution 一段时间
- 施加一次短脉冲
- 重复

### 建议 8：先做 `gamma = 0` 与 `gamma -> infinity` 两个极限

因为这两个极限最容易 benchmark：

- `gamma = 0`：非相互作用极限
- `gamma -> infinity`：Tonks-Girardeau 极限

如果这两个极限都做不稳，就不应该急着推进中间相互作用。

## 7.4 第四优先级：重新考虑大粒子数策略

### 建议 9：不要在现架构下直接冲 `N=10`

当前最大瓶颈之一就是显式 permutation 求和。

如果后面确实要走到更大粒子数，必须至少考虑以下之一：

- 用更适合 Bose 对称化的坐标与基底
- 减少显式 permutation 的代价
- 转向别的 many-body 表示方法

否则计算量和误差会同时爆炸。

---

## 8. 针对论文目标的研究建议

你提到最终目标是理解论文观察到的 many-body dynamical localization 的内在原因。

我建议不要把“找到精确基态”当成最终目标，而要把它看成：

- 准备可信初态
- 校准数值方法

真正重要的是后续这些问题：

### 8.1 现象层

- 为什么动量分布在多次 kick 后停止扩散？
- 为什么吸收能量会饱和？
- 为什么 entropy 不继续增长？
- 为什么强相互作用极限与非相互作用极限的 `G^(1)` 衰减形态不同？

### 8.2 机制层

值得重点测试的机制假设包括：

- 量子干涉在 momentum space 上维持了某种动态局域化
- 强相互作用并没有简单破坏局域化，而是把它转化成另一种 many-body 相干冻结
- trap、interaction、kick 三者共同决定了一个有限能量吸收窗口
- 局域化并不是单纯“没有热化”，而是 Floquet 结构、相干与约束共同作用的结果

### 8.3 数值策略层

最值得先做的不是“大体系高精度基态”，而是：

1. 先做最小体系的真实 kicked dynamics
2. 用可观测量直接判断是否出现 MBDL 类行为
3. 再分析相互作用强度、trap 深度、kick strength、pulse period 对该行为的影响

这条路线比一开始就冲大粒子数更稳。

---

## 9. 推荐的执行顺序

我建议按下面顺序推进。

### Step 1

重构静态优化器：

- 启用 refine
- `C` 正则化
- 每步重解 `u`
- 诊断输出完整化

### Step 2

专门解决两体 delta：

- CM/relative coordinate
- 多尺度 relative basis
- 确认误差能否稳定到 `1e-6`

### Step 3

实现真正的 kicked real-time TDVP / Floquet 演化

### Step 4

加入 observables：

- `n(k,t)`
- 能量
- entropy
- `G^(1)(z,t)`

### Step 5

先做 `N=2` 的物理图像，再尝试 `N=3,4`

### Step 6

只有在以上都稳定后，再讨论 `N~10`

---

## 10. 我对当前工作的总体判断

我的总体判断是：

1. 你现在遇到的 `1e-3 ~ 1e-6` 平台，主因不像显性 bug
2. `delta` 问题的主要矛盾是基底表达力不足
3. 当前优化器还有明显可改进空间，尤其是 `u` 的处理、`C` 的正则化和 refine 流程
4. 更大的问题是，当前程序还没有真正对齐论文里的驱动多体动力学
5. 如果继续在当前静态框架里只追求更准的“kicking benchmark”，很可能会投入很多算力，但并不能直接回答 MBDL 的内在原因

更值得做的是：

- 先把静态部分变成一个可靠的初态制备器
- 再把程序升级成真实的 kicked real-time many-body 演化工具

这才是通向论文目标的主路线。

---

## 11. 如果接下来要我继续做，最推荐的两个任务

如果下一步继续由我动手，我最推荐以下两个方向二选一：

### 方向 A：先把静态求解器做扎实

我可以直接修改代码，做：

- 每步重解 `u`
- 重新启用 refine
- 给 `C` 加正则
- 把 SVD 阈值与步长参数化
- 加上诊断输出

这条路线适合先解决你现在最直观的不满：静态 benchmark 精度和收敛平台。

### 方向 B：直接搭 kicked real-time 骨架

我可以直接开始重构程序，做：

- 显式时间依赖 kick
- one-period 演化框架
- 基本 observables 接口

这条路线更直接对齐 MBDL 论文目标。

---

## 12. 一句话总结

当前最核心的问题不是“代码对不对”，而是：

**你现在手上的程序，本质上还是一个两体静态 Gaussian 变分优化器；而你真正需要的是一个能处理 `trap + interaction + periodic kick` 的多体实时间/Floquet 传播器。**

