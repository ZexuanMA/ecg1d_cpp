# 组会报告：ECG 方法模拟 Kicked Rotor 的 Dynamical Localization

## 1. 研究背景与目标

### 1.1 物理问题

参考论文：Guo et al., "Observation of many-body dynamical localization", *Science* 389, 716 (2025)

该论文在一维玻色气体（Cs 原子）中实验观测到了多体 dynamical localization (MBDL)，但理论建模（量子 Monte Carlo + Floquet 方法）只能做到**定性**符合。（这个我暂时不太清楚，没仔细读过）

**我们的目标**：用 ECG (Effective Core Gaussian) + TDVP 方法提供**定量**模拟。

### 1.2 系统哈密顿量

$$
H(t) = \sum_i \left[\frac{p_i^2}{2m} + \frac{1}{2}m\omega^2 z_i^2 + \hbar\kappa\cos(2k_L z_i)\sum_n\delta(t-nT)\right] + g_{1D}\sum_{i<j}\delta(z_i-z_j)
$$

四个核心组成部分：

| 项 | 物理含义 | 代码状态 |
|----|---------|---------|
| $p^2/2m$ | 动能 | 已实现 |
| $\frac{1}{2}m\omega^2 z^2$ | 谐振子势阱 | 已实现 |
| $\kappa\cos(2k_L z)\sum_n\delta(t-nT)$ | 周期性 kick | 已实现 |
| $g_{1D}\delta(z_i-z_j)$ | 接触相互作用 | 已实现 |

### 1.3 Kick 项对应的幺正算子推导

哈密顿量中 kick 项含 $\delta(t-nT)$，意味着 kick 势只在 $t=nT$ 瞬间作用。考虑第 $n$ 次 kick 附近一个无穷小时间窗 $[nT-\epsilon,\, nT+\epsilon]$，薛定谔方程为：

$$i\hbar\frac{\partial\psi}{\partial t} = \left[H_{\text{free}} + \hbar\kappa\cos(2k_L z)\,\delta(t-nT)\right]\psi$$

在 $\epsilon\to 0$ 的极限下，$H_{\text{free}}$ 的贡献 $\sim O(\epsilon)\to 0$（有限值乘以无穷小时间），而 $\delta$ 函数项的贡献是有限的。对两边在 $[nT-\epsilon,\, nT+\epsilon]$ 上积分：

$$i\hbar\left[\psi(nT^+) - \psi(nT^-)\right] = \hbar\kappa\cos(2k_L z)\,\psi(nT)$$

即：

$$\psi(nT^+) = e^{-i\kappa\cos(2k_L z)}\,\psi(nT^-)$$

这就是 kick 算子：

$$\hat{U}_{\text{kick}} = e^{-i\kappa\cos(2k_L z)}$$

它是幺正的（$\hat{U}^\dagger\hat{U} = I$），因为指数的指数是纯虚数乘以实函数。物理上，kick 等价于在位置空间给波函数逐点乘一个相位——不改变概率密度 $|\psi|^2$，只改变相位结构（即动量分布）。

因此完整的一个 Floquet 周期的时间演化可以分解为：

$$\psi(t_{n+1}) = e^{-iH_{\text{free}}T/\hbar}\;\hat{U}_{\text{kick}}\;\psi(t_n)$$

先瞬时 kick（改变动量分布），再自由演化一个周期 $T$。

### 1.4 策略：先解决无相互作用情况

设 $g_{1D}=0$（$\gamma=0$），此时 N 体问题分解为 N 个独立的单粒子问题：

$$\Psi(z_1, z_2, t) = \psi(z_1, t)\cdot\psi(z_2, t)$$

**三大优势**：
- **有精确解**：单粒子可用有限差分精确求解，$E_2(n) = 2\times E_1(n)$
- **物理已知**：量子 kicked rotor 是量子混沌经典问题（短时线性增长 → 长时 dynamical localization）
- **隔离误差**：排除相互作用后，任何误差只来自 ECG 方法本身

---

## 2. Phase A 成果：虚时间演化获得高精度基态

### 2.1 方法流程

```
SVM 随机搜索增加基底 → Stochastic Refine 逐个微扰 → TDVP 虚时间演化精调
```

1. **SVM 阶段**：随机生成 Gaussian 基函数，只接受使基态能量下降的候选，直至能量下降阈值饱和
2. **Stochastic Refine**：逐个微扰已有基底参数，只接受能量继续下降的扰动
3. **TDVP 虚时间演化**：用 $C\dot{z} = -g$ 做连续参数优化，进一步提高精度

### 2.2 各 benchmark 精度

```
误差 (log scale)
  |
1e-2 ─ ■                                          ← Delta (2.0e-3)
  |
1e-3 ─
  |
1e-4 ─
  |
1e-5 ─           ■                    ■           ← Kicking (3.8e-6), Gaussian (2.3e-6)
  |
1e-6 ─
  |
  └────────────────────────────────────────
       Delta     Kicking    Gaussian
       (K=12)    (K=33)     (K=20)
```

| Benchmark | 精确能量 | SVM+Refine+TDVP 能量 | 误差 | 基底数 K | 耗时 |
|-----------|---------|---------------------|------|---------|------|
| 1-particle harmonic | 0.5 | 0.4999999999982 | **1.8e-12** | - | - |
| Gaussian interaction (N=2) | 1.5266998310 | 1.5267021 | **2.3e-6** | 20 | 8s |
| Kicking cos potential (N=2) | 2.4452547216 | 2.4452586 | **3.8e-6** | 33 | 23s |
| Delta contact (N=2) | 1.3067455 | 1.3088 | **2.0e-3** | 12 | 2.4s |

- 光滑势（Gaussian, Kicking）：ECG 指数收敛，K=20~33 即达 $10^{-6}$
- 奇异势（Delta cusp）：代数收敛 $\sim K^{-1/2}$，精度受限

---

## 3. 进入实时间演化：Kick 算子的处理

### 3.1 Kick 算子的 Bessel 展开

kick 算子是一个幺正算子：

$$\hat{U}_{\text{kick}} = \exp\left[-i\kappa\cos(2k_L x)\right]$$

利用 Jacobi-Anger 展开：

$$e^{-i\kappa\cos\theta} = \sum_{n=-\infty}^{\infty}(-i)^n J_n(\kappa)\,e^{in\theta}$$

取 $\theta = 2k_L x$，kick 算子变为：

$$\hat{U}_{\text{kick}} = \sum_{n=-\infty}^{\infty}(-i)^n J_n(\kappa)\,e^{in\cdot 2k_L x}$$

每一项 $e^{in\cdot 2k_L x}$ 是一个平移动量的算子。设 kick 前的波函数在动量表象下为 $\tilde\psi(p)$（即 $\psi(x)$ 的傅里叶变换），考虑其中一项的作用：

$$\psi'_n(x) = e^{in\cdot 2k_L x}\,\psi(x)$$

对 $\psi'_n$ 做傅里叶变换：

$$\tilde\psi'_n(p) = \int e^{-ipx}\,e^{in\cdot 2k_L x}\,\psi(x)\,dx = \int e^{-i(p - 2nk_L)x}\,\psi(x)\,dx = \tilde\psi(p - 2nk_L)$$

即 $e^{in\cdot 2k_L x}$ 在动量空间中将分布整体平移了 $2nk_L$。因此 kick 后的总动量分布为：

$$\tilde\psi'(p) = \sum_n (-i)^n J_n(\kappa)\,\tilde\psi(p - 2nk_L)$$

**物理含义**：kick 将原始动量分布复制到 $p = p_0 + 2nk_L$ 的一系列位置上，每个副本的权重由 Bessel 函数 $J_n(\kappa)$ 决定。对 $\kappa=1$：$J_0 = 0.765$, $J_{\pm1} = 0.440$, $J_{\pm2} = 0.115$, $J_{\pm3} = 0.020$, ...，高阶项快速衰减。

### 3.2 Kick 矩阵元的解析计算

对 Gaussian 基函数，平面波矩阵元有解析闭合形式：

$$\frac{\langle\phi_i|e^{in\cdot 2k_L x_a}|\phi_j\rangle}{M_G} = \exp\left(-n^2 k_L^2 K^{-1}(a,a) + in\cdot 2k_L\mu(a)\right)$$

其中 $K^{-1}$ 和 $\mu$ 来自 PairCache 已有的缓存量。高斯衰减因子 $e^{-n^2 k_L^2 K^{-1}}$ 保证高阶 Bessel 项快速衰减，截断到 $n_{\max}=20$ 绰绰有余。

### 3.3 投影更新线性系数

Kick 前的态是基底的线性组合：

$$|\Psi\rangle = \sum_j u_j |\phi_j\rangle$$

精确 kick 后的态为：

$$|\Psi'_{\text{ex}}\rangle = \hat{U}_{\text{kick}}|\Psi\rangle$$

一般来说，这个态**不再严格落在**当前有限基底子空间内（kick 注入了基底表达不了的高动量成分）。我们要找一个子空间内的最佳近似：

$$|\Psi'_{\text{pr}}\rangle = \sum_j u'_j |\phi_j\rangle$$

采用 Galerkin 条件——要求近似态在每个基函数方向上的投影与精确态一致：

$$\langle \phi_i | \Psi'_{\text{pr}} \rangle = \langle \phi_i | \Psi'_{\text{ex}} \rangle, \qquad \forall\, i$$

展开左边：

$$\langle \phi_i | \Psi'_{\text{pr}} \rangle = \sum_j \langle \phi_i | \phi_j \rangle\, u'_j = \sum_j S_{ij}\, u'_j$$

展开右边：

$$\langle \phi_i | \Psi'_{\text{ex}} \rangle = \sum_j \langle \phi_i | \hat{U}_{\text{kick}} | \phi_j \rangle\, u_j = \sum_j K_{ij}\, u_j$$

其中 $K_{ij} = \langle \phi_i | \hat{U}_{\text{kick}} | \phi_j \rangle$ 就是 kick 矩阵。两边对应即得：

$$S\mathbf{u}' = K\mathbf{u} \quad\Longrightarrow\quad \mathbf{u}' = S^{-1}K\mathbf{u}$$

**注意**：$K$ 矩阵不是 Hermitian（$\hat{U}_{\text{kick}}$ 是酉算符，$K_{ji} \neq \overline{K_{ij}}$），必须计算全部 $K\times K$ 个矩阵元。

---

## 4. 基底的表达力问题

### 4.1 核心矛盾

基态基底针对位置空间优化（B=0, R=0），不携带动量信息。每次 kick 会注入高动量成分，但这些成分不在基底子空间内——**投影时被永久丢弃**。

```
|ψ'⟩ = e^{-iκV} |ψ⟩          ← 无穷维，包含所有动量分量
    ↓ 投影到 K 维基底
|ψ'_proj⟩                      ← 丢失了高动量分量
    ↓ 自由演化（精确，零损失）
|ψ(T)⟩                         ← K 维内精确，但信息已经少了
    ↓ 下一次 kick 投影
|ψ'_proj⟩                      ← 又丢了一些 ...
```

### 4.2 验证指标

1. **Norm 守恒**：kick 是幺正算子，精确 kick 后 $\langle\psi|\psi\rangle$ 不变。投影后 norm 的偏差 $= 1 - \text{fidelity}$ 直接反映丢失量
2. **$\langle x^2\rangle$ 守恒**：$[U_{\text{kick}}, x^2]=0$，精确 kick 后 $\langle x^2\rangle$ 不变。投影后的偏差说明位置密度被扭曲
3. **能量 $E(n)$**：与有限差分精确解逐 kick 对比

### 4.3 当前最佳结果：动量增广基底

**策略**：在基态基底上添加携带动量的 Gaussian：

$$\phi_{\text{mom}}(x) \sim \text{Gaussian}(x)\cdot e^{i\cdot 2mk_L x}$$

通过 $R_a = i\cdot m\cdot k_L / b$ 编码动量。

**最优配置**：K=5（位置基）+ 12（动量基）= K=17

| 配置 | K | kick 1 误差 | kick 50 | 评价 |
|------|---|------------|---------|------|
| 纯位置 K=10 | 10 | 6% | 823% | 单调发散 |
| 纯位置 K=30 | 30 | 7% | 928% | 推迟但仍发散 |
| **K=5 + 动量(n=2, b=0.3)** | **17** | **5.5%** | **饱和 ~138%** | **最佳** |

前 3 个 kick 误差 <10%，但之后能量漂移。根本原因：每次 kick 投影丢失高动量成分，误差逐 kick 累积。

### 4.4 关于"饱和"的澄清：有限基底的谱上界

在固定的 K 维子空间里，态永远是这 K 个基函数的线性组合，其能量被子空间中哈密顿量的最大本征值封顶：

$$E = \frac{\mathbf{u}^\dagger H\mathbf{u}}{\mathbf{u}^\dagger S\mathbf{u}} \leq E_{\max}(H, S)$$

不管 kick 多少次、投影丢失多少信息，态的能量不可能超过 $E_{\max}$。因此上表中动量增广基底（K=17）在 kick 50 时能量趋于平稳，并不是物理上的 dynamical localization，而是**被基底的谱范围截断**了。纯位置基底（K=10, K=30）之所以看起来"不饱和"一直涨，只是因为它们的 $E_{\max}$ 更高，还没到顶。

### 4.5 应对策略：逐 kick 增扩动量基底

上述分析表明，固定基底方案有不可逾越的天花板。我们的策略是**每做一次 kick 就增扩一次动量基底**：

```
kick n → 投影 → 增扩新的动量基函数 → 自由演化 → kick n+1 → ...
```

每次 kick 后，波函数获得了新的高动量成分（$p = \pm 2k_L, \pm 4k_L, \ldots$）。如果基底中没有对应的基函数，这些成分在投影时就被丢弃。逐 kick 增扩确保子空间始终跟得上波函数在动量空间的扩展，从而：

1. **提高每次 kick 投影的 fidelity**——新增的动量基函数能承接 kick 注入的高动量成分
2. **打破谱上界限制**——基底维度 K 随 kick 数增长，$E_{\max}$ 随之提高，不再被截断

这也是为什么我们要在第一步（虚时间基态准备）就把精度做到尽可能高的原因：**每次 kick 的投影误差会逐 kick 累积，初始基态的误差是这条误差链的起点**。如果基态本身就有 $10^{-3}$ 的误差，后续再怎么增扩动量基底也无法补回；反之，如果基态精度达到 $10^{-6}$，我们就有更大的误差预算留给 kick 投影和长时间演化。

---

## 5. 下一步方向

回顾 H+K（含 kicking 项的虚时间演化）如何达到 3.8e-6 精度——它的 K=33 基底里有什么结构？

- 检查 kicking benchmark 基底的 A 矩阵分布：是否已经隐含了动量成分？
- 如果是，可以直接借用这些基底作为 kicked dynamics 的起点，也可以作为后边扩增基底的参考

---

## 8. 总结

| 里程碑 | 状态 |
|--------|------|
| 静态基态（各 benchmark）| 精度 $10^{-6}$ ~ $10^{-3}$ |
| 解析 kick 算子 | 已完成（Bessel 展开 + 解析矩阵元）|
| 有限差分精确解 | 已完成（单粒子 kicked rotor）|
| 固定基底精确传播 | 已完成（自由演化零误差）|
| 动量增广基底 | 已完成（前 3 kicks 误差 <10%）|

**核心结论**：ECG 基底的表达力是足够的（静态精度证明了这一点），当前瓶颈在于如何在 kick 投影时保留更多高动量信息。
