# Q&A: 为什么 Kicked Harmonic Oscillator 的 TDVP 实时间演化误差这么大？

## 问题

2体无相互作用、1D 谐振子势里的 kicked harmonic oscillator，kick 用 cos 函数。纯谐振子的实时间演化明明不难，加了瞬时 kick 之后误差就很大，为什么？

## 核心答案：cos kick 破坏了 Gaussian 基的"天然优势"

### 为什么纯谐振子很简单

谐振子的时间演化本质上是相空间里的**刚性旋转**。一个 Gaussian 波包在谐振子势中演化，**永远是 Gaussian**：

$$\psi(x,0) = e^{-ax^2} \quad \xrightarrow{e^{-iHt}} \quad \psi(x,t) = e^{-a(t)x^2 + b(t)x + c(t)}$$

所以 Gaussian 变分基函数 `phi(x) = exp(-x^T(A+B)x + 2R^T B x - ...)` 对谐振子来说是**精确表示**——不是近似，是精确的。TDVP 在 Gaussian manifold 上做投影的误差为零。

### 瞬时 kick 做了什么

kick 算子是：

```
U_kick = exp(-i * kappa * cos(2 * k_L * x))
```

用 Bessel 函数展开（`kick_operator.cpp` 里的实现）：

```
exp(-i*kappa*cos(theta)) = sum_{n=-inf}^{inf} (-i)^n * J_n(kappa) * exp(i*n*theta)
```

每个 `exp(i * n * 2k_L * x)` 项乘以一个 Gaussian，得到的是一个**带平面波因子的 Gaussian**——在动量空间里，这等于把动量平移 `n * 2k_L`。

所以一次 kick 之后，原来集中在 p ~ 0 的波函数变成了在 `p = 0, +/-2k_L, +/-4k_L, ...` 处都有分量的叠加态。

## 三层困难

### 1. 基的表达能力不足（最根本的问题）

静态优化基（SVM + TDVP 虚时间）是为谐振子基态量身定做的——所有基函数都集中在位置空间原点附近、动量 ~ 0。kick 之后波函数有高动量分量，但基里**没有能表达这些分量的函数**。

这就是 `augment_basis_with_momentum()` 存在的原因——它添加带动量的基函数。但这个方法有上限：加太多基函数会导致 overlap 矩阵 S 条件数爆炸。

### 2. 投影误差在每次 kick 累积

即使用精确 kick 算子（`apply_analytic_kick`），操作是：

```
u_new = S^{-1} * K * u     // K = kick matrix, S = overlap
```

这是一个**有损投影**：把精确的 kicked 态投影回有限基。误差量级取决于基对 kicked 态的完备程度。每次 kick 丢失信息，n 次 kick 之后误差累积。

### 3. TDVP 方法更糟糕

`kicked_evolution.cpp` 里的 TDVP 方法（`realtime_tdvp_step`），把 kick 当成 Hamiltonian 的一项来做无穷小时间步进。但 kick 是**有限变换**，不是无穷小的！

TDVP 方程 `C * dz/dt = -i * g` 假设状态沿切空间做无穷小运动，而 kick 算子一下把状态踢到切空间外面去了。这就像试图用泰勒展开来算 `exp(i*pi)`——你需要无穷多项才能收敛到 -1。

## 定量估算

用 `k_L = 0.5, kappa = 1.0` 来估算各阶 Bessel 函数权重：

| n | J_n(1.0) | 动量平移 n*2k_L |
|---|----------|-----------------|
| 0 | 0.765    | 0               |
| 1 | 0.440    | 1.0             |
| 2 | 0.115    | 2.0             |
| 3 | 0.020    | 3.0             |

n=1 项的权重 0.44 **不可忽略**。谐振子基态的动量展宽 `sigma_p ~ sqrt(m*omega/2) = 1/sqrt(2) ~ 0.7`，而 kick 平移动量 1.0，已经远超基态宽度。这意味着原始基函数（宽度 ~1 的 Gaussian）几乎完全"看不到"被踢走的那部分波函数。

## 根本出路

1. **精确 kick + 固定基精确演化**（`--kicked-gamma0` 模式）：基要**预先包含动量分量**，kick 用解析 Bessel 展开来精确投影，free evolution 在固定基里精确对角化。误差只来自于基的有限完备性。

2. **动态扩展基**：每次 kick 后，检查 kicked 态在当前基里的 fidelity，如果太低，添加新的动量基函数。但这导致 K 越来越大。

3. **更大的 kappa 或更多 kick 会更难**：kappa 越大，高阶 Bessel 函数权重越大，需要更多动量谐波。多次 kick 后动量分布会弥散到越来越宽的范围。

## 总结

谐振子简单是因为 Gaussian 基精确覆盖其 Hilbert 空间；cos kick 把状态踢出这个子空间。这本质上是**基的表达能力**问题，不是时间步长或算法的问题。

---

# Q&A: 初始态、实验对应、以及简化方案

## 问题

参考论文 Guo et al., Science 389, 716 (2025) 是冷原子实验，初始态是不是基态？实时间演化能不能简化？

## 1. 初始态确实是基态

论文明确说：

> "The experiment starts by loading a Bose-Einstein condensate (BEC) of Cs atoms ... the 1D systems are in equilibrium with a temperature of ~2 nK"

BEC 就是多体基态。实验流程：
1. 制备 BEC（基态）→ 装入 1D 管（~10^4 根管，每根 ~18 个原子）
2. 周期性地施加 cos 光晶格 kick
3. 测量动量分布 n(k)、动能、熵

所以模拟的初始态应该是 `H_0 = 动能 + 阱势 + delta相互作用` 的多体基态。这一步用虚时间 TDVP 已经能做。

## 2. 实验中阱势是 flat-bottom，不是纯谐振子

论文用蓝失谐 anti-trap 激光部分补偿了谐振约束，形成 flat-bottom potential。但弱谐振子 (nu_z = 14.7 Hz) 作为近似是合理的第一步。

## 3. gamma=0（无相互作用）不需要 ECG

这是最关键的简化。**无相互作用的 BEC** 中所有粒子处于同一个单粒子态：

```
Psi(z_1, ..., z_N) = phi(z_1) * phi(z_2) * ... * phi(z_N)
```

所有可观测量是单粒子的 N 倍：

```
E_N = N * E_1
n(k)_N = N * n(k)_1
S_N = N * S_1
```

所以 gamma=0 直接用 `kicked_exact_1particle`（有限差分网格，已实现为 `--kicked-exact`）就是精确解，根本不需要 ECG。

## 4. gamma>0 的简化策略

ECG 是必要的，但**不应该用 TDVP 来处理 kick**。正确策略是把 kick 和 free evolution 分离：

### Free evolution（两次 kick 之间）

Hamiltonian 是 `H_0 = 动能 + 阱势 + delta相互作用`，所有项在 ECG 框架内都有精确矩阵元。在**固定基**里做精确对角化传播：

```
u(t) = S^{-1/2} U exp(-i Lambda t) U† S^{1/2} u(0)
```

这就是 `free_evolve_fixed_basis`。误差只取决于基的完备性，**不随时间步数积累**（和 TDVP 逐步积累误差完全不同）。

### Kick

用解析 Bessel 展开（`apply_analytic_kick`）精确计算 kick 算子在基里的矩阵元，然后做投影。不要用 TDVP 来处理 kick。

### 完整流程

```
1. 虚时间 TDVP 找基态 → 得到基函数 {phi_k} 和系数 u
2. augment_basis_with_momentum() → 添加动量基函数
3. 循环 N_p 次：
   a. apply_analytic_kick(basis, kappa, k_L)   // Bessel 展开，精确投影
   b. free_evolve_fixed_basis(basis, H, S, T)  // 对角化精确传播
   c. 记录可观测量
```

这就是 `--kicked-gamma0` 模式的策略。

## 5. 好消息：MBDL 物理限制了基的需求

核心困难是基的完备性——每次 kick 把波函数推到更高动量，需要更多基函数。

但 **MBDL 的物理本身就是说动量分布会冻结**（localization）！论文 Fig. 2 清楚地显示：

- gamma=0：n(k) 在 ~300-500 kicks 后冻结（Anderson localization in momentum space）
- gamma=11：n(k) 在第一次 kick 后就几乎冻结了

这意味着动量不会无限扩散，**基的需求有上限**。具体来说：

- 实验参数 K = 3.3, k_L 对应的最大动量大约在 +/-10 * k_L 范围内
- 只需要预先准备覆盖这个动量范围的基函数就够了

## 6. 实际参数对应

论文实验参数（自然单位转换后供参考）：

| 参数 | 实验值 | 说明 |
|------|--------|------|
| 原子 | Cs-133 | m = 133 u |
| k_L | pi/a, a=532.25 nm | 光晶格波矢 |
| K (kick strength) | 3.3 或 4.4 | K = V_L T_p / (2 hbar) |
| T (kick period) | 60 或 80 μs | |
| T_p (pulse duration) | 10 μs | 有限脉冲宽度 |
| V_L (lattice depth) | 20 E_r | E_r = pi^2 hbar^2 / (2ma^2) |
| gamma | 0 到 11 | gamma = m g_{1D} / (n_{1D} hbar^2) |
| N (每根管) | ~18 | 加权平均 |
| N_p (kick 数) | 1 到 1101 | |

注意 **kick 不是理想 delta 函数**，脉冲宽度 T_p = 10 μs 是有限的。但 T_p << T，所以近似为瞬时 kick 是合理的第一步。

## 总结

1. **初始态 = 基态**，用虚时间 TDVP 准备，这已经能做
2. **gamma=0 用网格精确解**，不需要 ECG
3. **gamma>0 用 ECG**，但 kick 用 Bessel 精确投影，free evolution 用对角化精确传播，**不用 TDVP 做实时间演化**
4. MBDL 物理本身保证动量不会无限扩散，基的需求有天然上限

---

# Q&A: 没有精确解时如何验证？

## 问题

gamma=0（无相互作用）ECG 误差已经很大（~100-300%），对于没有精确解的 gamma>0，怎么判断做得对不对？用能量守恒还是动量守恒？

## 首先：能量和动量在 kicked 系统中都不守恒

- **能量不守恒**：kick 本身会注入能量，这是 kicked rotor 的核心物理（能量吸收 → 饱和 = dynamical localization）
- **动量不守恒**：cos(2k_L z) 势破坏了连续平移对称性（只有离散平移对称性 z → z + pi/k_L）
- 所以这两个**都不能直接作为验证判据**

## 当前 gamma=0 的误差诊断

`--kicked-gamma0` 输出显示：

```
[kick] norm before=6.28319, norm after=8.80389, ratio=1.40118, S cond=8.643869e+05
```

| 指标 | 当前值 | 理想值 | 含义 |
|------|--------|--------|------|
| kick norm ratio | 1.40 | 1.00 | 每次 kick 有 ~40% 的信息丢失/扭曲 |
| S condition number | 8.6e5 | < 1e4 | overlap 矩阵病态 |
| Energy error (kick 10) | 130% | < 1% | 能量完全错误 |

**根本原因：基只有 5 个基态函数 + 少量动量函数，完全不够表达 kicked 态。**

## 正确的验证策略：三层保险

### 第一层：Kick Fidelity（不需要精确解）

`norm_after / norm_before` 是最直接的基质量指标。kick 算子是幺正的，对精确表示应有 ratio = 1.0。

- ratio = 1.00 ± 0.01：基对 kick 表达很好
- ratio > 1.05：基不够，结果不可信
- ratio = 1.40（当前值）：完全不行

**这个指标对 gamma>0 同样有效**，因为它只检查基能不能表达 kicked 态，和相互作用无关。

### 第二层：基收敛性（不需要精确解）

对同一物理问题，改变基的大小 K，观察结果是否收敛：

```
K = 20: E(kick=10) = ?
K = 40: E(kick=10) = ?
K = 60: E(kick=10) = ?
```

如果 K=40 和 K=60 的差异 < 1%，说明基已收敛。如果还在剧烈变化，说明 K 不够。

### 第三层：gamma→0 极限（最强验证）

gamma>0 代码取 gamma → 0 时，结果必须趋近精确的网格解。这同时验证了：
- 基的完备性
- kick 投影的正确性
- free evolution 的正确性
- 整个代码链路

## 正确的工作顺序

```
Step 0: 有限差分网格精确解（gamma=0）—— 已完成
        这是 ground truth，不可能错

Step 1: ECG gamma=0 对齐精确解
        目标：kick fidelity > 0.99，能量误差 < 1%
        手段：增大 K，优化动量基函数的选取
        这一步不过关，不要进入 Step 2

Step 2: ECG gamma → 0 极限验证
        在有相互作用的代码中取 g_contact → 0
        结果必须和 Step 0 一致
        这验证了整个代码链路

Step 3: ECG gamma > 0 物理结果
        此时有三层信心保障：
        - kick fidelity 指示基的质量
        - 基收敛性指示结果的可靠性
        - gamma→0 极限已经通过验证
```

## 关于两次 kick 之间的能量守恒

虽然总能量不守恒，但有一个弱版本可以用：**两次 kick 之间**，Hamiltonian 是 H_0（时间无关），能量应该精确守恒。

如果用 `free_evolve_fixed_basis`（对角化精确传播），这个自动满足（by construction）。

如果用 TDVP 做 free evolution，可以检查 `E(t_after_kick) vs E(t_before_next_kick)` 是否一致。但推荐用固定基对角化，这样这个问题就不存在了。

## Kick Fidelity 的物理原理

### 幺正算子保持 norm

kick 算子 `U = exp(-i*kappa*cos(2k_L*x))` 是幺正的（指数上是 Hermitian 算子乘以 -i），幺正意味着 U†U = I，所以：

```
<ψ'|ψ'> = <ψ|U†U|ψ> = <ψ|ψ>
```

这是量子力学的基本要求：幺正演化保持概率守恒。norm = 总概率 = 1。

### 有限基中的投影问题

波函数在 K 个基函数张成的子空间 V 中表示：`|ψ> = Σ u_k |φ_k>`。kick 后精确的态 `U|ψ>` 不一定还在 V 里。代码做的是求解 `S u' = K u`（即 `u' = S^{-1} K u`），这等于把 kicked 态投影回子空间。

把精确的 kicked 态分解：

```
U|ψ> = P_V U|ψ> + (I - P_V) U|ψ>
        ↑                ↑
    基能表达的部分    丢失的部分（基外分量）
```

两部分正交，所以：

```
||U|ψ>||² = ||P_V U|ψ>||² + ||(I-P_V) U|ψ>||²
```

如果基完备，丢失部分 = 0，投影后 norm 不变。

### 为什么 ratio > 1（norm 增大了）

正交投影只会让 norm 减小或不变，不可能增大。但 `S^{-1} K` 不是正交投影——当 S 病态时，S^{-1} 会放大小奇异值方向的分量，导致 norm 膨胀。

ratio ≠ 1 同时反映两件事：
1. **基不完备**：kicked 态有分量在子空间之外
2. **S 病态**：S^{-1} 放大了数值误差

ratio = 1.40 意味着每次 kick 有 ~40% 的信息被扭曲，结果完全不可信。目标：ratio = 1.00 ± 0.01。

## 当前最紧迫的任务

**不是去做 gamma>0，而是让 gamma=0 的 ECG 结果和精确解对齐。**

---

# Q&A: 什么是动量基函数

## ECG 基函数的一般形式

$$\varphi(\mathbf{x}) = \exp\!\bigl[-\mathbf{x}^T (A+B)\,\mathbf{x} + 2\,\mathbf{R}^T B\,\mathbf{x} - \mathbf{R}^T B\,\mathbf{R}\bigr]$$

其中 $A, B$ 是 $N \times N$ 矩阵，$\mathbf{R}$ 是 $N$ 维向量，$\mathbf{x} = (x_1, \dots, x_N)^T$ 是各粒子坐标。

## 静态基函数（$B=0, R=0$）

SVM 和虚时间 TDVP 优化出来的基函数，取 $B=0$, $\mathbf{R}=0$，退化为纯 Gaussian：

$$\varphi(\mathbf{x}) = \exp\!\bigl[-\mathbf{x}^T A\,\mathbf{x}\bigr]$$

这是中心在原点、**零动量**的 Gaussian 波包。只能描述基态附近、动量 $\approx 0$ 的状态。

## 动量基函数（$B \neq 0$, $R$ 有虚部）

设 $B_{aa} = b$（实数），$R_a = i\,p/b$（纯虚），代入基函数后，对粒子 $a$ 的相关项展开：

$$\exp\!\bigl[\dots + 2 \cdot \tfrac{ip}{b} \cdot b \cdot x_a + \dots\bigr] = \exp\!\bigl[\dots + 2ip\,x_a + \dots\bigr]$$

这个 $e^{2ip\,x_a}$ 就是动量为 $2p$ 的**平面波因子**。整个基函数变成：

$$\varphi(\mathbf{x}) = \underbrace{\exp\!\bigl[-\mathbf{x}^T (A+B)\,\mathbf{x} + \dots\bigr]}_{\text{Gaussian 包络}} \;\times\; e^{2ip\,x_a}$$

**物理意义**：一个动量中心在 $2p$ 处的 Gaussian 波包。

## 为什么 kicked 演化需要动量基函数

Kick 算子用 Jacobi-Anger / Bessel 展开：

$$e^{-i\kappa\cos(2k_L x)} = \sum_{n=-\infty}^{\infty} (-i)^n J_n(\kappa)\, e^{in \cdot 2k_L x}$$

每个 $e^{in \cdot 2k_L x}$ 项把粒子动量平移 $n \cdot 2k_L$。kick 后的波函数在动量空间散布到 $p = 0, \pm 2k_L, \pm 4k_L, \dots$。

零动量的 Gaussian 基函数根本"看不到"这些高动量分量（内积指数衰减）。所以需要在 $\mathbf{R}$ 里编码不同的动量值，让基函数能覆盖 kicked 态的动量分布。

## 代码中的对应

```cpp
// 粒子 a 携带动量 m × 2k_L
B_new(a, a) = b_val;                           // B 对角元
R_new(a) = Cd(0.0, m * k_L / b_val);           // R 纯虚部
// → 基函数包含因子 exp(i·m·2k_L·x_a)
```

设 $p = m \cdot k_L$，则 $R_a = ip/b$，基函数包含 $e^{2ip\,x_a} = e^{i \cdot m \cdot 2k_L \cdot x_a}$，正好对应动量平移 $m \cdot 2k_L$ 的 Bessel 项。

## 一句话总结

**动量基函数 = Gaussian 包络 × 平面波**，是在 ECG 变分框架内表达非零动量态的方式。

---

# Q&A: 什么是 $(p,p)$ 型基函数，什么叫"所有粒子共享同一个动量"

## Kick 算子对多粒子的作用

Kick 对每个粒子**独立**作用（$\cos$ 势是单体算子的和）：

$$U_{\text{kick}} = \prod_{a=1}^{N} e^{-i\kappa\cos(2k_L x_a)}$$

对 $N=2$，用 Bessel 展开每个因子：

$$U_{\text{kick}} = \biggl[\sum_m (-i)^m J_m(\kappa)\,e^{im \cdot 2k_L x_1}\biggr] \times \biggl[\sum_n (-i)^n J_n(\kappa)\,e^{in \cdot 2k_L x_2}\biggr]$$

所以 kicked 态是所有 $(m, n)$ **组合**的叠加，权重为 $J_m(\kappa) \cdot J_n(\kappa)$：

$$U_{\text{kick}} |\psi\rangle = \sum_{m,n} (-i)^{m+n} J_m J_n \; e^{im \cdot 2k_L x_1}\,e^{in \cdot 2k_L x_2} \;|\psi\rangle$$

## $(p,p)$ 型 = 两个粒子携带相同动量

一个动量基函数的 $\mathbf{R}$ 向量决定了每个粒子的动量。$N=2$ 时 $\mathbf{R} = (R_1, R_2)$。

**旧代码**对动量阶数 $n$ 设置：

```cpp
for (int a = 0; a < N; a++) {
    R_new(a) = i * n * k_L / b_val;   // 所有粒子的 R 相同！
}
```

即 $R_1 = R_2 = i \cdot n \cdot k_L / b$，两个粒子都携带**同一个**动量 $n \cdot 2k_L$。这就是 $(p, p)$ 型——上面展开式中 $m = n$ 的那些项。

旧代码遍历 $n = -2, -1, +1, +2$，生成的基函数对应：

| 基函数 | 粒子1动量 | 粒子2动量 | 对应 Bessel 项 |
|--------|-----------|-----------|----------------|
| $n=1$  | $+2k_L$   | $+2k_L$   | $(m,n)=(1,1)$, 权重 $J_1^2 = 0.194$ |
| $n=2$  | $+4k_L$   | $+4k_L$   | $(m,n)=(2,2)$, 权重 $J_2^2 = 0.013$ |
| $n=-1$ | $-2k_L$   | $-2k_L$   | $(m,n)=(-1,-1)$, 权重 $J_1^2 = 0.194$ |
| ...    | ...       | ...       | ... |

## 缺失的混合动量项

kicked 态中 $m \neq n$ 的项——比如粒子1不动（$m=0$）、粒子2被踢（$n=1$）——旧代码**完全不生成**：

| 缺失的基函数 | 粒子1动量 | 粒子2动量 | 权重 |
|-------------|-----------|-----------|------|
| $(0, +1)$   | $0$       | $+2k_L$   | $J_0 \cdot J_1 = 0.765 \times 0.440 = \mathbf{0.337}$ |
| $(0, -1)$   | $0$       | $-2k_L$   | $J_0 \cdot J_1 = \mathbf{0.337}$ |
| $(0, +2)$   | $0$       | $+4k_L$   | $J_0 \cdot J_2 = \mathbf{0.088}$ |
| $(+1, -1)$  | $+2k_L$   | $-2k_L$   | $J_1^2 = \mathbf{0.194}$ |
| $(+1, +2)$  | $+2k_L$   | $+4k_L$   | $J_1 \cdot J_2 = \mathbf{0.051}$ |

**混合项 $(0, \pm1)$ 的权重 0.337 比同动量项 $(\pm1, \pm1)$ 的 0.194 还大！** 这是 kicked 态中最重要的分量，但旧基完全没有能表达它的函数。

## 新代码的修正

新代码生成所有**无序对** $(m_1, m_2)$，$m_1 \leq m_2$：

```cpp
for (int m1 = -n_mom; m1 <= n_mom; m1++) {
    for (int m2 = m1; m2 <= n_mom; m2++) {
        if (m1 == 0 && m2 == 0) continue;  // 跳过基态
        R_new(0) = i * m1 * k_L / b_val;   // 粒子1独立动量
        R_new(1) = i * m2 * k_L / b_val;   // 粒子2独立动量
    }
}
```

这样就能覆盖所有 $(m, n)$ 组合，包括混合动量项。

## 为什么取无序对（$m_1 \leq m_2$）而不是所有 $(m_1, m_2)$

因为波函数要玻色对称化。$(m_1, m_2) = (0, 1)$ 和 $(1, 0)$ 经过对称化后是同一个态——代码中 `PermutationSet` 会自动处理粒子交换。所以只需要 $m_1 \leq m_2$ 避免冗余。

---

# Q&A: 什么是动量格点，为什么截断到 $n_\text{mom}=4$，为什么需要 $(2n_\text{mom}+1)^2$ 个

## 从 Bessel 展开说起

Kick 算子的精确展开是**无穷级数**：

$$e^{-i\kappa\cos(2k_L x)} = \sum_{n=-\infty}^{+\infty} (-i)^n J_n(\kappa)\, e^{in \cdot 2k_L x}$$

每一项 $e^{in \cdot 2k_L x}$ 把粒子动量平移 $n \cdot 2k_L$（下面证明），前面的系数 $J_n(\kappa)$ 是第一类 Bessel 函数，决定了这一项的**权重**（振幅大小）。

### 为什么 $e^{iqx}$ 是动量平移算子——Fourier 变换的平移定理

一个粒子的波函数 $\psi(x)$ 在动量空间的表示是 Fourier 变换：

$$\tilde{\psi}(p) = \int_{-\infty}^{+\infty} e^{-ipx}\,\psi(x)\,dx$$

（取 $\hbar = 1$）。现在考虑在位置空间**乘以平面波** $e^{iqx}$，得到新波函数 $\phi(x) = e^{iqx}\psi(x)$。它的动量表示是：

$$\tilde{\phi}(p) = \int e^{-ipx}\,e^{iqx}\,\psi(x)\,dx = \int e^{-i(p-q)x}\,\psi(x)\,dx = \tilde{\psi}(p - q)$$

这就是 Fourier 变换的**平移定理**：

$$\boxed{\text{位置空间乘以 } e^{iqx} \;\Longleftrightarrow\; \text{动量空间平移 } p \to p + q}$$

直觉理解：$e^{iqx}$ 本身就是动量为 $q$ 的平面波（$\hat{p}\,e^{iqx} = -i\partial_x e^{iqx} = q\,e^{iqx}$）。乘上它相当于给粒子"叠加"了一个额外的动量 $q$。

### 应用到 kick 展开

取 $q = n \cdot 2k_L$，Bessel 展开中的第 $n$ 项：

$$e^{in \cdot 2k_L \cdot x}\,\psi(x) \;\longleftrightarrow\; \tilde{\psi}(p - n \cdot 2k_L)$$

即把整个动量分布**刚性平移** $n \cdot 2k_L$。

举例：基态波函数是 Gaussian，动量分布集中在 $p = 0$ 附近。$n = 1$ 项乘上 $e^{i \cdot 2k_L \cdot x}$ 后，动量分布变为集中在 $p = 2k_L$ 附近的 Gaussian：

```
原始 ψ(p):     n=1 项:         n=-1 项:        n=2 项:

   ╱╲             ╱╲              ╱╲                ╱╲
  ╱  ╲           ╱  ╲            ╱  ╲              ╱  ╲
─╱────╲──    ──╱────╲──    ──╱────╲──    ────╱────╲──
  p=0           p=2k_L         p=-2k_L          p=4k_L
```

kick 后的总波函数是这些平移分量的**加权叠加**（权重 $J_n(\kappa)$），在动量空间形成一系列等间距的峰。

## Bessel 函数随 $n$ 快速衰减——这是截断的物理依据

$J_n(\kappa)$ 有一个关键性质：**当 $|n|$ 明显大于 $\kappa$ 时，$J_n$ 指数级衰减趋近于零**。直觉上，$\kappa$ 控制了 kick 的"强度"，kick 越强（$\kappa$ 越大），能把粒子踢到越高的动量，需要保留越多的项。

以当前代码使用的 $\kappa = 1.0$ 为例，各阶 Bessel 函数的数值：

| $n$ | $J_n(1.0)$ | $J_n^2$（概率权重） | 累积概率 | 说明 |
|-----|------------|---------------------|----------|------|
| 0   | 0.7652     | 0.5856              | 58.6%    | 不动的分量，最大 |
| 1   | 0.4401     | 0.1937              | 97.3%    | $n=0,\pm1$ 已覆盖 97% |
| 2   | 0.1149     | 0.0132              | 99.9%    | |
| 3   | 0.0196     | 0.0004              | ~100%    | |
| 4   | 0.0025     | 0.000006            | ~100%    | 几乎为零 |
| 5   | 0.0002     | —                   | —        | 完全可忽略 |

（这里"累积概率"是 $J_0^2 + 2\sum_{k=1}^{n} J_k^2$，因为 $\pm n$ 对称）

**$n=0$ 和 $n=\pm 1$ 就占了 97% 的权重。** 到 $n=\pm 2$ 已经 99.9%。$n \geq 3$ 的项加起来不到 0.1%。

## 为什么选 $n_\text{mom} = 4$

选 $n_\text{mom} = 4$ 意味着保留 $|n| \leq 4$ 的所有项，丢弃 $|n| \geq 5$ 的项。对 $\kappa = 1.0$：

- 保留的项覆盖了 >99.999% 的权重
- 丢弃的最大项 $J_5(1.0) = 0.0002$，平方后 $\sim 10^{-8}$，完全可忽略

这个截断引入的误差远小于其他误差来源（比如基不完备导致的 fidelity 问题）。

**经验法则**：$n_\text{mom} \approx \kappa + 几$ 就够了。$\kappa = 1$ 时 $n_\text{mom} = 2$ 其实就够，取 4 是留余量。

如果将来用论文中的实验参数 $K = 3.3$（对应更大的 $\kappa$），Bessel 函数衰减更慢，就需要更大的 $n_\text{mom}$（大概 $n_\text{mom} \sim 6\text{--}8$）。

## 截断之后：动量空间变成离散格点

截断到 $|n| \leq n_\text{mom}$ 后，kick 只能把粒子踢到有限个离散动量值上：

$$p \in \{-n_\text{mom} \cdot 2k_L, \;\dots,\; -2k_L, \; 0, \; +2k_L, \;\dots,\; +n_\text{mom} \cdot 2k_L\}$$

共 $2n_\text{mom} + 1$ 个值。这些值在动量空间等间距排列，就像一维网格上的格点，所以叫**动量格点**。

以 $n_\text{mom} = 4$ 为例，单粒子有 $2 \times 4 + 1 = 9$ 个格点：

```
p / (2k_L):   -4   -3   -2   -1    0   +1   +2   +3   +4
               ·    ·    ·    ·    ★    ·    ·    ·    ·
              极小  极小  小   重要  最大  重要  小   极小  极小
```

注意这些格点的**权重差别巨大**：中心（$n=0$）最大，越远越小。但只要权重不为零，对应的基函数就不能完全缺失，否则 kick 投影会把那部分波函数丢掉。

## 多粒子：动量格点的直积

$N$ 个粒子的 kick 是**独立**作用的（见上节）。粒子1可以在任何格点上，粒子2也可以在任何格点上，两者独立。所以 $N=2$ 的总状态空间是两个粒子动量格点的**直积**：

$$\text{总格点数} = \underbrace{(2n_\text{mom}+1)}_{\text{粒子1}} \;\times\; \underbrace{(2n_\text{mom}+1)}_{\text{粒子2}} = (2n_\text{mom}+1)^2$$

$n_\text{mom} = 4$ 时：$9 \times 9 = 81$ 个格点。

这 81 个格点可以画成二维网格，每个格点 $(m_1, m_2)$ 代表"粒子1动量 $m_1 \cdot 2k_L$，粒子2动量 $m_2 \cdot 2k_L$"：

```
m2  +4 | ·  ·  ·  ·  ·  ·  ·  ·  ·
    +3 | ·  ·  ·  ·  ·  ·  ·  ·  ·
    +2 | ·  ·  ·  ·  ·  ·  ·  ·  ·
    +1 | ·  ·  ·  ·  ·  ·  ·  ·  ·
     0 | ·  ·  ·  ·  ★  ·  ·  ·  ·    ← 基态在 (0,0)
    -1 | ·  ·  ·  ·  ·  ·  ·  ·  ·
    -2 | ·  ·  ·  ·  ·  ·  ·  ·  ·
    -3 | ·  ·  ·  ·  ·  ·  ·  ·  ·
    -4 | ·  ·  ·  ·  ·  ·  ·  ·  ·
       +-----------------------------
  m1    -4  -3  -2  -1   0  +1  +2  +3  +4
```

旧代码只覆盖**对角线**上的格点 $(m, m)$——即上图的对角线，共 8 个（除去原点）。新代码覆盖上三角（$m_1 \leq m_2$，玻色对称化），共 $\frac{9 \times 10}{2} - 1 = 44$ 个（减去基态 $(0,0)$）。

## 为什么实际需要的基函数远多于 44 个

每个动量格点 $(m_1, m_2)$ 只确定了基函数的**动量**（由 $\mathbf{R}$ 编码），但还没确定**空间宽度**（由 $A$ 矩阵对角元决定）。

不同宽度的 Gaussian 描述不同的物理：
- 窄 Gaussian（$A$ 大）→ 位置空间局域，动量空间弥散
- 宽 Gaussian（$A$ 小）→ 位置空间弥散，动量空间局域

为了充分描述 kicked 态在每个动量格点上的空间结构，每个格点通常需要**多个不同宽度**的基函数。如果用 $n_w$ 种宽度，总共需要：

$$K_{\text{total}} = n_w \times 44 \;\text{（动量格点数）}$$

用 3 种宽度就是 $3 \times 44 = 132$ 个基函数。这就是为什么估计 $K > 100$。

## 当前的瓶颈

新代码生成了 $15 \times 44 = 660$ 个候选，但 S 条件数约束（$< 10^6$）只允许 12 个通过。22 个基函数（10 基态 + 12 动量）要覆盖 44 个动量格点 × 多种宽度，严重不足。这就是 fidelity 达到 1.79、能量误差 >400% 的根本原因。

## 一般情形：$N$ 粒子

$N$ 个粒子的动量格点数是 $(2n_\text{mom}+1)^N$，指数增长。这就是多体问题的困难所在：

| $N$ | $n_\text{mom}=4$ 格点数 |
|-----|------------------------|
| 1   | 9                      |
| 2   | 81                     |
| 3   | 729                    |
| 4   | 6561                   |

但玻色对称化后独立格点数大约是 $\binom{2n_\text{mom}+N}{N}$（组合数），比直积小很多。$N=2$ 时 $\binom{9+1}{2} = 45$（含原点），$N=3$ 时 $\binom{9+2}{3} = 165$。

---

# Q&A: gamma=0 kicked 基质量改进尝试（2026-04-01）

## 做了什么

### 问题诊断

旧代码 `augment_basis_with_momentum()` 的核心缺陷：只生成 $(p,p)$ 型基函数（所有粒子共享同一动量），缺失权重最大的混合动量项。详见上一节。

kick 后各分量的 Bessel 权重：

| 分量 | 权重 |
|------|------|
| (0, ±1) 混合动量 | J_0 × J_1 = 0.765 × 0.440 = **0.337** |
| (±1, ±1) 同动量 | J_1² = **0.194** |
| (0, ±2) 混合动量 | J_0 × J_2 = 0.765 × 0.115 = **0.088** |

**混合动量项权重最大**，但旧代码完全不生成这些基函数。

### 代码修改

#### 1. `kick_operator.cpp` — `augment_basis_with_momentum()` 重写

- **旧**：只生成 R(0)=R(1)=同一动量值（所有粒子共享动量 n×k_L）
- **新**：生成所有**无序动量对** (m1, m2)，m1 ≤ m2，跳过 (0,0)
  ```
  R(0) = i * m1 * k_L / b_val
  R(1) = i * m2 * k_L / b_val
  ```
- 新增：从 ground state basis 提取 A 对角值加入 width 候选列表（使动量函数的空间宽度匹配基态）
- `max_cond` 参数化（原来硬编码 1e6）

#### 2. `kick_operator.hpp` — 签名更新

- `augment_basis_with_momentum()` 新增 `max_cond` 参数
- `apply_analytic_kick()` 返回值 `void` → `double`（返回 fidelity = norm_after/norm_before）

#### 3. `main.cpp` — `run_kicked_gamma0_test()` 参数更新

| 参数 | 旧 | 新 |
|------|-----|-----|
| K_max (基态基函数) | 5 | 10 |
| n_mom (动量阶数) | 2 | 4 |
| b_val (B 值) | 0.3 | 0.5 |
| SVM n_trials | 5000 | 8000 |

新增 kick fidelity 汇总（min/max/mean），使用 `apply_analytic_kick` 返回的真实 fidelity 值而非重归一化后的 norm ratio。

#### 4. `kick_operator.cpp` — `apply_analytic_kick()` 诊断改进

简化输出，每次 kick 打印 fidelity + PASS/FAIL 判定（阈值 ±5%）。

### 结果

```
Candidate pool: 660 (15 widths × 44 momentum pairs)
Augmented basis: 10 original + 12 momentum = 22 total   ← 660 个候选只有 12 个通过
S conditioning: cond = 9.9e5 (接近 1e6 上限)

Kick Fidelity: min=1.008, max=1.792, mean=1.610         ← 远超 1.05 阈值
Kick 50 energy: ECG=8.97, exact=1.74, rel_error=416%    ← 远超 10% 阈值
```

### 分析：为什么 660 个候选只通过了 12 个

候选池 = 15 widths × 44 momentum pairs = 660。贪心选择（按候选顺序逐个检查 S 条件数 < 1e6）只接受了 12 个。

原因是同一 width 下的多个动量对之间有很高的 overlap（特别是宽的 Gaussian 函数之间），导致 S 迅速变病态。每加一个新基函数都推高 S 条件数，很快就碰到 1e6 上限。

### 第二次尝试：放宽 max_cond + 减少 widths

修改：
- widths：15 个 → 5 个（{0.3, 0.8, 2.0, 6.0} + ground state A 对角值）
- max_cond：1e6 → 1e10

结果：

```
Candidate pool: 528 (12 widths × 44 momentum pairs)
Augmented basis: 10 original + 14 momentum = 24 total   ← 只多了 2 个
S conditioning: cond = 1e10 (打满上限)

Kick 1-17: fidelity PASS (< 1.05)     ← 短期有改善！
Kick 18+:  fidelity 急速恶化到 1.55
Kick 50 energy: ECG=16.83, exact=1.74, rel_error=869%   ← 比之前还差
```

**分析**：

1. **短期改善**：前 17 次 kick 的 fidelity 都 < 1.05，说明混合动量基函数确实有帮助——kicked 态的主要分量能被表达了
2. **长期恶化**：每次 kick 把波函数推到更高动量，14 个动量函数只覆盖了有限的格点。到 kick 18，波函数已经扩散到基覆盖不了的高动量区域
3. **条件数代价**：S cond = 1e10 意味着 double 精度（16 位）只剩 6 位有效数字。S^{-1} 放大了数值噪声，反而让长期误差比 max_cond=1e6 时更大
4. **只多了 2 个基函数**：即使放宽到 1e10，贪心选择还是很快碰壁

### 结论：固定 Gaussian 基 + 贪心 S conditioning 这条路走不通

根本矛盾：
- **需要**：~44 个动量格点 × ≥1 种宽度 = 至少 44 个动量基函数
- **能加**：贪心选择最多加 ~14 个（因为 Gaussian 之间 overlap 太高，S 条件数快速爆炸）
- **差距**：3 倍以上

Gaussian 基函数之间天然有很高的 overlap（特别是宽度接近的函数）。这不是调参数能解决的——是 Gaussian 型基的内禀限制。

### 可能的出路

1. **正交化基底**：不用原始 Gaussian，而是先对候选池做 S^{-1/2} 正交化，在正交基里工作。S 变成单位矩阵，条件数问题消失。代价：失去 Gaussian 形式的解析矩阵元优势
2. **动态基扩展**：每次 kick 后，不用固定基，而是针对当前 kicked 态新生成/优化基函数。每次 kick 的基可以不同
3. **平面波 + Gaussian 混合基**：高动量分量用平面波（自然正交），低动量/空间结构用 Gaussian
4. **放弃 ECG 做 gamma=0**：正如 qa.md 之前所述，gamma=0 直接用网格精确解。ECG 的价值只在 gamma>0（有相互作用时才需要 Gaussian 基来处理 delta 势的矩阵元）

---

# Q&A: 当前代码状态复核与建议（2026-04-01）

## 实际复跑结果

直接运行：

```bash
./build/ecg1d --kicked-gamma0
```

得到的关键结果与上节结论一致：

```text
Candidate pool: 528
Augmented basis: 10 original + 14 momentum = 24 total
S conditioning: cond = 9.98e9

Kick 1-17: fidelity 基本都在 1.05 以内
Kick 18+:  fidelity 开始持续失控
Kick 50 energy: ECG = 16.83, exact = 1.737, rel_error = 8.69
```

这说明：

1. **混合动量项修复已经落地**：`augment_basis_with_momentum()` 现在确实生成所有无序 `(m1,m2)` 动量对，而不再只是旧的 `(p,p)` 型
2. **当前主瓶颈已经不是“漏掉了主要动量分量”**：前十几次 kick fidelity 明显改善，说明主要低阶分量已经被覆盖
3. **真正卡住的是固定原始 Gaussian 基 + 贪心 cond 筛选**：基函数数目加不上去，且 S 已经病态到 `~1e10`

## 新判断

到这一步，继续微调下面这些参数，预计不会从根本上解决问题：

- `K_max`
- `n_mom`
- `b_val`
- `max_cond`
- width 候选列表

原因不是“参数还没调到最佳”，而是当前方法的结构性瓶颈：

- 原始 Gaussian 动量基之间 overlap 太高
- 贪心逐个加入候选时，很快就触到 S 条件数上限
- 放宽 cond 上限虽然能多留少量基函数，但 `S^{-1}` 会更强地放大数值误差

所以当前路线更像是在一个已经饱和的方法上继续拧参数，而不是解决核心表示问题。

## 还存在的工程阻塞

### 1. 物理参数被硬编码

`g_contact`, `g_gauss`, `kappa`, `k_L` 等都还是 `physical_constants.hpp` 里的编译期常量。

后果：

- 不能方便地做 `gamma -> 0` 极限验证
- 不能系统扫描参数
- 每次改参数都要重新编译，不利于做收敛性实验

### 2. 主程序入口太“嘈杂”

`main()` 现在不管跑什么模式，都会先跑：

- functionals
- derivatives
- Hamiltonian gradients

这些输出对 `--kicked-gamma0` 这种实验模式不是必须的，会干扰批量测试和结果记录。

## 建议的推进顺序

### Step A：先做工程整理，不先做新物理

先把下面这些改成运行时参数，而不是编译期常量：

- `g_contact`
- `g_gauss`
- `kappa`
- `k_L`
- `T_period`
- `n_kicks`
- `K_max`
- `n_mom`
- `b_val`
- overlap truncation / conditioning 参数

同时把 `--kicked-gamma0` 做成一个干净入口，只输出：

- kick fidelity
- `E(n)`
- `cond(S)`
- basis size

这样才能真正做系统实验，而不是靠手工改代码试参数。

### Step B：先测试“贪心选基”是不是根本瓶颈

优先做一个**正交化子空间原型**：

1. 先生成完整 momentum candidate pool
2. 在候选池上构造 overlap 矩阵 `S`
3. 用 `S` 的本征分解或 SVD 做全局正交化 / 截断
4. 在这个正交化子空间里做 kick + free evolution
5. 再与 exact gamma=0 结果比较

这个实验最关键，因为它能区分两种情况：

- **情况 A**：问题主要是贪心选基太差，换成全局正交化后 gamma=0 会明显改善
- **情况 B**：即使做了正交化，gamma=0 仍然很差，那说明 Gaussian 子空间本身就不适合长期 kicked dynamics

这一步的信息价值远高于继续手调 `n_mom/max_cond`。

### Step C：通过 gamma=0 后，再谈 gamma>0

只有当下面三条成立时，才值得认真推进相互作用情形：

1. gamma=0 的 kick fidelity 能长期稳定接近 1
2. gamma=0 的结果对 basis size / truncation 参数收敛
3. 有相互作用代码在 `g_contact -> 0` 时能回到 gamma=0 exact 参考解

在这之前，不建议直接进入 gamma>0 的物理解读，因为链路验证还不够闭合。

## 当前最推荐的方向

如果只选一个最值得投入的下一步，我的建议是：

**先做“运行时参数化 + 正交化子空间 gamma=0 原型”**。

理由：

1. 这一步实现成本比 Rothe / 动态基扩展低很多
2. 它能最快回答“ECG 路线还有没有继续做下去的价值”
3. 如果连这一步都不能显著改善 gamma=0，就应该停止在固定 Gaussian 基上继续投入

## 如果正交化原型仍然失败

那后续策略应调整为：

1. **gamma=0 永远使用网格 exact 作为基准线**
2. **ECG 只服务于 gamma>0**
3. **先做 `g_contact -> 0` 极限验证**
4. 再决定是否值得上更重的方法：
   - Rothe 方法
   - 动态基扩展
   - 平面波 + Gaussian 混合基

一句话总结：

**当前最该做的不是继续调固定 Gaussian 基的参数，而是先验证“正交化后的 Gaussian 子空间”是否还值得继续。**

---

# Q&A: Step A 缩减版与 Step B 的真正意义（2026-04-01）

## 关于 Step A：已经按“最小改动”收缩

后续讨论中，认为没必要在这一步就把所有物理常数都改成运行时参数。因为当前最紧迫的不是扫物理参数，而是先把 gamma=0 的实验链路变得可重复、可比较。

所以 Step A 已经收缩为：

1. **保留物理常数硬编码**
   - `g_contact`
   - `g_gauss`
   - `kappa`
   - `k_L`

2. **只开放最直接影响 gamma=0 验证的实验参数**
   - `--n-kicks`
   - `--n-mom`
   - `--max-cond`

3. **保留默认值不变**

也就是说，下面这条命令的默认行为不变：

```bash
./build/ecg1d --kicked-gamma0
```

只是现在也可以方便地做最小实验：

```bash
./build/ecg1d --kicked-gamma0 --n-kicks 20 --n-mom 6 --max-cond 1e8
```

4. **把 `--kicked-gamma0` 的输出改成实验模式**

保留：

- 参数摘要
- 基态能量摘要
- candidate pool / basis size
- `cond(S)`
- 每个 kick 的 `fidelity`
- 每个 kick 的 `E_ECG / E_exact / rel_error`

去掉：

- functionals
- derivatives
- Hamiltonian gradients

因为这些对 kicked gamma=0 验证不是关键指标，只会淹没真正要看的信号。

## Step B 的意义，不只是“再试一个新方法”

Step B 的核心目的不是为了把代码写得更复杂，而是为了回答一个**决定路线是否继续投入**的问题：

> 当前失败，到底是因为“贪心选基策略太差”，还是因为“Gaussian 子空间本身就不适合长期 kicked dynamics”？

这两个原因表面看起来很像，但工程含义完全不同。

### 可能性 1：贪心选基是主要问题

当前 `augment_basis_with_momentum()` 的做法是：

1. 生成很多候选动量 Gaussian
2. 按顺序一个个尝试加入
3. 只要一加入就让 `cond(S)` 超阈值，就拒掉

这个策略很局部、很保守。它的问题是：

- 候选的质量不是按全局最优排序的
- 早期留下的一些基函数，可能并不是对 kicked dynamics 最重要的
- 一旦前面选得不理想，后面真正有用的候选会因为 `cond(S)` 被卡掉

如果问题主要在这里，那么：

- 同一批 Gaussian 候选
- 不用贪心逐个筛
- 改成“先全部收进来，再做全局正交化/截断”

就可能明显改善 gamma=0 结果。

这意味着：

**ECG 路线本身还有价值，当前只是选基算法不够好。**

### 可能性 2：Gaussian 子空间本身就不合适

另一种可能是：即使把候选池整体正交化，结果还是很差。

这说明问题不再是“你怎么选这些 Gaussian”，而是：

- cos kick 会不断把态推向更宽的动量分布
- 原始 Gaussian 包络乘平面波的表达方式，对长期 kicked 态仍然不够自然
- 即使你数值上把它们正交化，子空间的物理表达能力还是不够

如果验证到这一步，那么结论就很明确：

**不应该继续在固定 Gaussian 基上投入大量时间调参。**

此时更合理的路线会变成：

- gamma=0 永远用 exact 网格做参考
- ECG 只保留给 gamma>0
- 或者转向更重的方法：Rothe、动态基扩展、混合基

## 为什么 Step B 值得先做

因为它的信息价值极高，但实现成本还没有高到失控。

和继续手调 `n_mom`、`b_val`、`max_cond` 相比，Step B 能更快回答“这条路还有没有必要继续”。

也就是说，Step B 不是为了追求一个更 fancy 的算法，而是为了避免我们在错误方向上继续耗时间。

## Step B 实际上在做什么

所谓“正交化子空间原型”，并不是要改物理模型，而是把同一批候选基换一种数值组织方式：

1. 先生成完整 momentum candidate pool
2. 对这整个候选池构造 overlap 矩阵 `S`
3. 对 `S` 做本征分解或 SVD
4. 丢掉特别小的奇异值方向
5. 在剩余的正交子空间中做 kick 和 free evolution

这样做之后：

- `S` 不再病态
- 不需要靠贪心一个个拒基
- 可以更公平地测试“这批 Gaussian 候选整体上到底有没有表达 kicked 态的能力”

## 一句话理解 Step B

Step B 本质上是在做一次“路线体检”：

- **如果它成功**：说明应该继续优化 ECG 的子空间构造
- **如果它失败**：说明应尽快停止在固定 Gaussian 基上继续深挖

所以 Step B 的最大价值不是直接给出最终物理解，而是帮助判断**这条方法论值不值得继续投入**。

## 对话记录策略

从这一步开始，和项目主线直接相关的讨论结论会继续追加到 `qa.md` 里，作为工作日志和决策记录。

---

# Q&A: 当前整个流程到底是什么（2026-04-01）

## 短答案

是的，**当前主线流程就是**：

1. 先用**静态 Gaussian 基**做 SVM + 虚时间演化，把**基态**找出来
2. 再在这个基态基的基础上，**补充动量基函数**
3. 然后在这个**增广后的固定基**里做 kicked real-time evolution

但要注意一个关键点：

**实时间阶段并不是再用 TDVP 去处理 kick。**

当前更推荐、也已经实现的做法是：

- kick：`apply_analytic_kick()`
- 两次 kick 之间的自由演化：`free_evolve_fixed_basis()`

也就是说：

**虚时间 TDVP 用来准备基态；实时间 kicked dynamics 用固定基精确传播。**

## 现在代码里的实际顺序

以 `--kicked-gamma0` 为例，当前代码顺序是：

### Step 1：用静态基找基态

在 `B=0, R=0` 的静态 Gaussian 模式下：

1. `svm_build_basis()` 做初始基构造
2. `stochastic_refine()` 做随机优化
3. `evolution(..., terms_free, ...)` 做虚时间 polish

这一步的目标不是直接做 kicked dynamics，而是得到一套适合**基态/低动量态**的 Gaussian 基。

### Step 2：在基态基上加动量基

然后调用：

```cpp
augment_basis_with_momentum(refined.basis, k_L, n_mom, b_val, max_cond)
```

把原来的静态基，扩展成：

- 原始基态 Gaussian
- 带不同动量对 `(m1,m2)` 的 Gaussian × plane-wave 基

也就是说，**动量基函数是在基态基已经准备好之后再加的**，不是一开始就和基态优化混在一起做。

## Step 3：在增广后的固定基里做实时间 kicked evolution

接下来：

1. 先用增广基重新构造 `H_aug` 和 `S_aug`
2. 用 `set_u_from_eigenvector()` 把初态重新放回增广基的基态
3. 进入 kicked real-time loop

循环里每一步做的是：

```text
kick:          apply_analytic_kick(...)
free evolve:   free_evolve_fixed_basis(...)
measure:       E, fidelity, rel_error
```

这一步里：

- **不会**再用实时 TDVP 去近似 kick
- **不会**更新 `A/B/R`
- 只是在固定基里更新线性系数 `u`

## 所以整条思路可以写成

```text
静态 Gaussian 基
    ↓
SVM 找初始基
    ↓
stochastic refine
    ↓
虚时间 TDVP polish 基态
    ↓
在基态基上增广动量基函数
    ↓
重建增广基的 H, S
    ↓
把初态放到增广基的基态上
    ↓
[ kick + fixed-basis free evolution ] × n_kicks
```

## 为什么不是“一开始就把动量基也放进去一起找基态”

因为两类基函数服务的物理目标不同：

1. **静态 Gaussian 基**
   - 主要负责把 trapped ground state 做准
   - 最擅长低动量、局域态

2. **动量基函数**
   - 主要负责表示 kick 之后出现的非零动量分量
   - 是为 real-time kicked dynamics 补上的

如果一开始就把大量动量基也混进基态优化：

- overlap 矩阵更容易病态
- SVM 搜索更难
- 会把“基态优化”和“kicked 表示能力”两个问题缠在一起

所以当前流程把它们拆开，是合理的。

## 当前流程的物理含义

可以把它理解成：

- 前半段：先找到“系统在不踢的时候”的最佳静态表示
- 后半段：再给这套表示补上“被 kick 踢出零动量子空间之后”所需的动量方向

这比“直接拿基态基去做 kicked dynamics”更合理，也比“用 TDVP 硬扛 kick”更稳定。

## 当前真正的问题不在流程顺序，而在增广后子空间的质量

也就是说，**“先找基态，再加动量基”这个总流程本身并没有问题。**

当前卡住的是后半段：

- 增广后的 raw Gaussian 动量基之间 overlap 太高
- 贪心筛选后能留下来的动量基太少
- 导致 kick 几次之后 fidelity 和能量误差迅速恶化

所以接下来的 Step B，不是推翻这个总体流程，而是改进其中的这一步：

> “基态基之后，如何构造一个足够好的 kicked 子空间？”

一句话总结：

**当前流程是对的：先基态，再补动量，再做固定基实时间传播；问题出在‘补出来的 kicked 子空间质量不够’，而不是这个先后顺序本身。**

---

# Q&A: 在基态基础上加动量基，凭什么认为“加得好”（2026-04-01）

## 短答案

严格说，**现在并不能先验认定“加得好”**。

当前做法只能说是：

- **物理上合理**
- **数值上可实现**
- **但本质上仍然是启发式构造**

真正能判断“加得好不好”的，不是“看起来像合理的动量基”，而是**事后验证指标**。

## 当前代码到底凭什么去“加”

现在的动量基构造，依据的是三层直觉：

### 1. 动量位置来自 kick 的 Bessel 展开

因为

$$
e^{-i\kappa\cos(2k_L x)}=\sum_n (-i)^n J_n(\kappa)e^{i n 2k_L x}
$$

所以 kick 后自然会出现：

- 0
- ±2k_L
- ±4k_L
- ...

这些离散动量谐波。

因此当前代码会去生成对应 `(m_1,m_2)` 的动量对。这一步的物理依据是明确的，不是瞎猜。

### 2. 空间宽度借用基态基的信息

当前实现不会随便乱造宽度，而是：

- 先给一组稀疏的默认 width
- 再把基态优化出来的 `A` 对角元吸收进 width 候选

意思是：希望新加的动量基，在空间包络上不要离基态基太远。

这也是一个**合理启发式**，但仍然不是最优性证明。

### 3. 用 `cond(S)` 做数值可解性约束

如果新基函数一加进去就让 overlap 矩阵太病态，就拒掉。

这一步保证的是：

- 不是“物理最好”
- 而是“至少数值上还能算”

所以 `max_cond` 本质上是**数值稳定性阈值**，不是“物理正确性阈值”。

## 所以：当前加基的逻辑是“合理候选生成”，不是“已经证明最佳”

这点必须说清楚。

当前方法能回答的问题只有：

> 这些候选动量基，是否值得拿来试？

它回答不了的问题是：

> 这些候选是不是已经足够好、足够完备、足够接近 kicked 态的最佳子空间？

后者必须靠验证。

## 那么“加得好”到底该怎么判断

### 判据 1：Kick fidelity

这是最直接的标准。

如果基真的能表达 kicked 态，那么：

```text
fidelity = norm_after / norm_before ≈ 1
```

也就是说，kick 之后投影回基空间时，不应该丢掉太多信息，也不应该因为病态矩阵把 norm 扭曲掉。

所以：

- fidelity 长时间接近 1：说明加得比较好
- fidelity 很快偏离 1：说明加得不够

这是当前最重要的判据。

### 判据 2：gamma=0 exact benchmark

因为 gamma=0 有精确网格解，所以我们不能只说“fidelity 看起来还行”，还要看：

- `E_ECG` 和 `E_exact` 是否一致
- 随 kick 数增长后误差是否可控

如果一套动量基真的加得好，gamma=0 不应该在十几个 kick 后就彻底漂走。

### 判据 3：基收敛性

不是看某一组参数跑得像不像，而是看结果是否稳定：

- `n_mom = 4`
- `n_mom = 6`
- `n_mom = 8`

或者不同的 truncation / conditioning 设定下，结果是否收敛。

如果一改参数结果就乱跳，说明“加得好”这个说法根本站不住。

### 判据 4：一阶增量价值

更细一点的看法是：

一个新加的候选基函数，应该能显著改善下列量中的至少一个：

- 单次 kick fidelity
- 与 exact 的一步误差
- 长期误差增长速度

如果加进去以后这些都没改善，只是把 `cond(S)` 变差，那它就不是一个“好”的新增方向。

## 当前代码实际上只能做“事后判定”

也就是说，现在我们不是：

> 我知道这些动量基一定好，所以我加它们

而是：

> 这些动量基在物理上有理由出现，所以我先把它们作为候选；它们到底好不好，要靠 fidelity / exact / 收敛性来验

这正是为什么 Step B 很重要。

## Step B 和这个问题的关系

你问“拿什么认定加得好”，其实正好打到 Step B 的核心。

因为当前方法里，“加得不好”可能有两层原因：

1. **候选本身没问题，但贪心选法太差**
2. **候选子空间本身就不够好**

Step B 做正交化子空间原型，就是为了把这两层分开。

如果 Step B 成功，说明：

- 这些候选动量基本身未必错
- 错的是当前“怎么挑、怎么组织”它们

如果 Step B 失败，说明：

- 不是挑选顺序的问题
- 而是这类 Gaussian 动量基整体上就不够适合

## 一句话总结

**当前并不是“我确信这些动量基加得好”，而是“我知道它们为什么应该被纳入候选；至于是不是真的加得好，必须靠 kick fidelity、exact benchmark 和收敛性来证明”。**

---

# Q&A: 当前方法的完整数学过程（2026-04-01）

这部分专门回答下面几个问题：

1. **我们一开始的波函数 ansatz 到底是什么**
2. **基态是怎么找出来的**
3. **动量基函数是怎么加的**
4. **候选库到底有哪些**
5. **加完之后到底做了什么**
6. **`u' = S^{-1}Ku` 这个投影公式是怎么来的**

---

## 1. 总波函数 ansatz

当前代码里，整个多体波函数并不是直接在网格上存，而是写成有限个基函数的线性组合：

$$
\Psi(\mathbf{x}) = \sum_{j=1}^{K} u_j \,\phi_j(\mathbf{x})
$$

其中：

- $\mathbf{x} = (x_1,\dots,x_N)^T$
- $u_j \in \mathbb{C}$ 是线性系数
- $\phi_j$ 是第 $j$ 个 ECG 型基函数

代码里对应：

- `BasisParams::u` = $u_j$
- `BasisParams::A, B, R` 决定 $\phi_j$

## 2. 单个基函数的形式

当前使用的基函数形式是

$$
\phi(\mathbf{x})
=
\exp\!\Big[
-\mathbf{x}^T(A+B)\mathbf{x}
+ 2\,\mathbf{R}^T B \mathbf{x}
- \mathbf{R}^T B \mathbf{R}
\Big]
$$

这里：

- $A$ 是对称矩阵，控制 Gaussian 宽度和粒子间相关
- $B$ 在当前实现里通常取对角形式
- $R$ 是复向量，可以编码中心位移或动量相位

### 为什么 $R$ 的虚部会编码动量

看第 $a$ 个粒子的相关部分：

$$
\exp\!\bigl[2 R_a B_{aa} x_a \bigr]
$$

如果取

$$
B_{aa}=b,\qquad R_a = i\,\frac{m_a k_L}{b}
$$

那么

$$
2R_a B_{aa} x_a
=
2\left(i\frac{m_a k_L}{b}\right) b x_a
=
i\, 2 m_a k_L x_a
$$

于是基函数里就出现了平面波因子

$$
e^{i\,2m_a k_L x_a}
$$

这就是“Gaussian 包络 × 动量平面波”。

所以当前代码中的动量基，其实就是：

$$
\text{动量基} = \text{Gaussian 包络} \times e^{i\,2m_1k_L x_1} \times e^{i\,2m_2k_L x_2}\times \cdots
$$

---

## 3. 在固定一组基函数下，基态怎么求

给定一组基函数 $\{\phi_j\}_{j=1}^K$ 后，代码会构造：

$$
S_{ij} = \langle \phi_i | \phi_j \rangle
$$

和

$$
H_{ij} = \langle \phi_i | \hat H | \phi_j \rangle
$$

然后解广义本征值问题：

$$
H \mathbf{u} = E\, S \mathbf{u}
$$

最低本征值对应基态能量，最低本征向量给出线性系数 $\mathbf{u}$。

代码里对应：

- `build_HS(...)`
- `lowest_energy(...)`
- `set_u_from_eigenvector(...)`

### 为什么是广义本征值问题

因为基函数不是正交的。

如果基是正交归一的，那么就直接是普通本征值问题：

$$
H \mathbf{u} = E \mathbf{u}
$$

但现在

$$
\langle \phi_i | \phi_j \rangle \neq \delta_{ij}
$$

所以必须带上 overlap 矩阵 $S$。

### 数值上怎么解

代码实际做的是先把 $S$ 正交化：

$$
S = V \Lambda V^\dagger
$$

定义

$$
S^{-1/2} = V \Lambda^{-1/2} V^\dagger
$$

然后把广义本征值问题化成普通 Hermitian 本征值问题：

$$
\widetilde H = S^{-1/2} H S^{-1/2}
$$

解

$$
\widetilde H \mathbf{c} = E \mathbf{c}
$$

再回到原系数：

$$
\mathbf{u} = S^{-1/2}\mathbf{c}
$$

---

## 4. 基态基是怎么“找出来”的，而不是直接给定的

当前 `--kicked-gamma0` 的前半段，不是拿一组固定基直接算，而是三步：

### Step 1: SVM 构造初始静态基

先随机生成很多静态 Gaussian 候选（此时基本是 $B=0, R=0$ 或近似静态形式），然后逐步挑出能降低最低能量的基函数。

目标是：

$$
\min_{\text{basis set}} E_0(\text{basis})
$$

这里优化的是**基函数集合本身**。

### Step 2: stochastic refine

对已有基函数做随机微扰，如果新基函数替换旧基函数后能让最低能量更低，就接受。

目标还是把这组基变得更适合表示基态。

### Step 3: 虚时间 TDVP polish

在已有基上，对非线性参数（主要是 $A,B,R$ 的某些分量）做虚时间方向的小步更新，使能量继续下降。

形式上是求解

$$
C\,\delta z = -g
$$

这里：

- $z$ 表示所有可变参数
- $g$ 是对能量的梯度
- $C$ 是 TDVP 的几何矩阵

这一步的目的不是做 real-time kicked dynamics，而是把“**表示基态的这组静态基**”再打磨一下。

---

## 5. 然后怎么加动量基

这一部分才是你关心的核心。

### 5.1 出发点

先有一组已经优化好的基态基：

$$
\mathcal{B}_{\text{gs}} = \{\phi_1,\dots,\phi_{K_0}\}
$$

这组基善于表示：

- 基态
- 低动量
- 靠近 trap 中心的态

但 kick 之后会出现离散高动量分量，所以要再补一个动量子空间：

$$
\mathcal{B}_{\text{mom}}
$$

最后工作空间是：

$$
\mathcal{B}_{\text{aug}} = \mathcal{B}_{\text{gs}} \cup \mathcal{B}_{\text{mom}}
$$

### 5.2 备选基底库有哪些

当前代码不是“想加一个就加一个”，而是先生成一个**候选库**：

$$
\mathcal{C} = \{\chi_\alpha\}
$$

候选库由两部分笛卡尔积构成：

1. **width 候选**
2. **动量对候选**

#### (a) width 候选

当前 width 集合大致是：

$$
W = \{0.3,\ 0.8,\ 2.0,\ 6.0\} \cup \{\text{ground-state basis 中 } A_{aa}\text{ 的对角值}\}
$$

后半部分会去重，只保留和已有 width 差异足够大的那些。

意思是：

- 一部分 width 是人工给的稀疏采样
- 一部分 width 直接借用基态基已经学到的空间尺度

#### (b) 动量对候选

对 $N=2$，当前代码会枚举

$$
(m_1,m_2),\qquad m_1 \le m_2,\qquad m_1,m_2 \in \{-n_{\text{mom}},\dots,n_{\text{mom}}\}
$$

并跳过

$$
(0,0)
$$

因为 $(0,0)$ 已经被静态基态基覆盖了。

所以当前一个候选基函数，实际上由

$$
(w,m_1,m_2)
$$

三元组决定。

### 5.3 候选基函数的具体形式

对每个 $(w,m_1,m_2)$，当前代码构造：

$$
A_{\text{new}}
=
\begin{pmatrix}
w & 0.1w \\
0.1w & w
\end{pmatrix},
\qquad
B_{\text{new}}
=
\begin{pmatrix}
b & 0 \\
0 & b
\end{pmatrix}
$$

并取

$$
R_{\text{new}}
=
\begin{pmatrix}
i\,m_1 k_L / b \\
i\,m_2 k_L / b
\end{pmatrix}
$$

于是这个候选基函数带有：

$$
e^{i\,2m_1k_L x_1}\,e^{i\,2m_2k_L x_2}
$$

的动量结构。

也就是说，当前代码的候选库是：

$$
\mathcal{C}

=
\left\{
\chi_{w,m_1,m_2}
\;\middle|\;
w\in W,\;
m_1\le m_2,\;
(m_1,m_2)\neq(0,0)
\right\}
$$

---

## 6. 候选库生成后，怎么决定“哪些真加进去”

当前代码不是把整个候选库都收下，而是做一个**贪心筛选**。

设当前已经接受的基集是

$$
\mathcal{B}^{(n)}
$$

对每个候选 $\chi$：

1. 形成试探基集

$$
\mathcal{B}_{\text{trial}} = \mathcal{B}^{(n)} \cup \{\chi\}
$$

2. 计算它的 overlap 矩阵

$$
S_{\text{trial}}
$$

3. 如果

$$
\kappa(S_{\text{trial}})
=
\frac{\lambda_{\max}(S_{\text{trial}})}{\lambda_{\min}(S_{\text{trial}})}
<
\texttt{max\_cond}
$$

就接受；否则拒绝。

所以当前并不是“按贡献最大去选”，而是：

> 按候选遍历顺序，谁不会让 $S$ 太病态，就先留下谁

这就是为什么我一直强调：

**当前方法的瓶颈，很可能不只是候选不够，而是这个选法本身太弱。**

---

## 7. 加完之后到底做了什么

假设最终留下的增广基是

$$
\mathcal{B}_{\text{aug}} = \{\phi_1,\dots,\phi_K\}
$$

然后做三件事：

### Step 1: 重建增广基的 $H$ 和 $S$

$$
S_{ij} = \langle \phi_i | \phi_j \rangle,
\qquad
H_{ij} = \langle \phi_i | \hat H_0 | \phi_j \rangle
$$

这里 $\hat H_0$ 是 free Hamiltonian，不含 kick。

### Step 2: 把初态重新放回增广基的基态

也就是在这套新基上再解一次

$$
H \mathbf{u}_0 = E_0 S \mathbf{u}_0
$$

得到增广基里的初始系数 $\mathbf{u}_0$。

注意，这一步不是说“动量基也参与了基态物理”，而是：

> 同一个物理基态，现在改用更大的工作空间来表示

### Step 3: 开始 kicked real-time loop

之后每一步循环做：

1. kick
2. free evolution
3. 测量

---

## 8. Kick 的数学过程：为什么会出现 $u' = S^{-1} K u$

这是最关键的一步。

### 8.1 真实想做的事情

假设当前态是

$$
|\Psi\rangle = \sum_j u_j |\phi_j\rangle
$$

真正的 kick 后态是

$$
|\Psi_{\text{exact}}'\rangle = \hat U_{\text{kick}} |\Psi\rangle
$$

其中

$$
\hat U_{\text{kick}} = e^{-i\kappa \sum_a \cos(2k_L x_a)}
$$

### 8.2 问题：kick 后态一般不在当前有限子空间里

也就是说，通常不存在一个系数向量 $\mathbf{u}'$，使得

$$
|\Psi_{\text{exact}}'\rangle
=
\sum_j u'_j |\phi_j\rangle
$$

精确成立。

所以我们必须做“投影回当前子空间”。

### 8.3 当前代码采用的投影条件

代码采用的是 Galerkin 条件：

$$
\langle \phi_i | \Psi'_{\text{proj}} \rangle
=
\langle \phi_i | \Psi'_{\text{exact}} \rangle,
\qquad \forall i
$$

其中

$$
|\Psi'_{\text{proj}}\rangle
=
\sum_j u'_j |\phi_j\rangle
$$

代入左边：

$$
\langle \phi_i | \Psi'_{\text{proj}} \rangle
=
\sum_j \langle \phi_i | \phi_j \rangle u'_j
=
\sum_j S_{ij} u'_j
$$

代入右边：

$$
\langle \phi_i | \Psi'_{\text{exact}} \rangle
=
\left\langle \phi_i \middle| \hat U_{\text{kick}} \middle| \Psi \right\rangle
=
\sum_j \langle \phi_i | \hat U_{\text{kick}} | \phi_j \rangle u_j
=
\sum_j K_{ij} u_j
$$

这里定义 kick 矩阵

$$
K_{ij} = \langle \phi_i | \hat U_{\text{kick}} | \phi_j \rangle
$$

于是就得到线性方程组：

$$
S \mathbf{u}' = K \mathbf{u}
$$

如果 $S$ 可逆，就有

$$
\mathbf{u}' = S^{-1} K \mathbf{u}
$$

这就是代码里 `apply_analytic_kick()` 的核心。

### 8.4 这不是“正交投影”，而是非正交基下的 Galerkin 投影

这一点很重要。

如果基是正交归一的，那么 $S=I$，就直接得到：

$$
\mathbf{u}' = K \mathbf{u}
$$

但现在基不正交，所以必须先经过 $S^{-1}$。

这也是病态问题的根源之一：

- 如果 $S$ 很病态
- 那么 $S^{-1}$ 会强烈放大小奇异值方向

于是投影误差和数值噪声都会被放大。

---

## 9. Kick 矩阵 $K_{ij}$ 是怎么解析算出来的

当前代码没有数值积分 `exp(-i κ cos)`，而是用 Bessel 展开：

$$
e^{-i\kappa\cos\theta}
=
\sum_{n=-\infty}^{\infty}
(-i)^n J_n(\kappa)e^{in\theta}
$$

取

$$
\theta = 2k_L x_a
$$

就得到单粒子 kick 的谐波展开。

对于 Gaussian 基函数之间的矩阵元，代码使用了解析结果：

$$
\frac{
\langle \phi_i | e^{i n 2k_L x_a} | \phi_j \rangle
}{
\langle \phi_i | \phi_j \rangle
}
=
\exp\!\left(
-n^2 k_L^2 K^{-1}_{aa}
+ i n\,2k_L\,\mu_a
\right)
$$

这里：

- $K$ 来自左右两个 Gaussian 指数相加
- $\mu = \tfrac12 K^{-1} b$

于是每个粒子的 kick kernel 可以解析求和，最后总 kick kernel 对不同粒子做乘积，就得到整体系的 $K_{ij}$。

---

## 10. Free evolution 又是怎么做的

在两次 kick 之间，Hamiltonian 不含时间依赖，且基固定不变。

所以在有限子空间里，系数满足：

$$
i S \frac{d\mathbf{u}}{dt} = H \mathbf{u}
$$

这也是一个广义线性演化问题。

代码同样先把 $S$ 正交化，解广义本征值问题：

$$
H \mathbf{v}_n = E_n S \mathbf{v}_n
$$

然后把初态展开到这些广义本征模上：

$$
\mathbf{u}(0) = \sum_n c_n \mathbf{v}_n
$$

自由演化就是

$$
\mathbf{u}(t)
=
\sum_n c_n e^{-iE_n t}\mathbf{v}_n
$$

所以当前 free evolution 是**固定基中的精确传播**，不是 TDVP 逐步近似。

---

## 11. 整个流程压缩成一套公式

### 基态准备阶段

1. 选一组静态 Gaussian 基 $\mathcal{B}_{\text{gs}}$
2. 构造

$$
H^{(\text{gs})},\ S^{(\text{gs})}
$$

3. 解

$$
H^{(\text{gs})}\mathbf{u}_0 = E_0 S^{(\text{gs})}\mathbf{u}_0
$$

得到基态

### 动量增广阶段

4. 生成候选库

$$
\mathcal{C} = \{\chi_{w,m_1,m_2}\}
$$

5. 通过 `cond(S)` 贪心筛选，得到

$$
\mathcal{B}_{\text{aug}} = \mathcal{B}_{\text{gs}} \cup \mathcal{B}_{\text{mom}}
$$

6. 在增广基上重新解基态：

$$
H_{\text{aug}} \mathbf{u}_0^{\text{aug}}
=
E_0^{\text{aug}} S_{\text{aug}} \mathbf{u}_0^{\text{aug}}
$$

### Kicked evolution 阶段

每一步：

7. kick 投影

$$
S_{\text{aug}} \mathbf{u}^{+} = K_{\text{aug}} \mathbf{u}^{-}
$$

8. free evolution

$$
i S_{\text{aug}} \frac{d\mathbf{u}}{dt}
=
H_{\text{aug}} \mathbf{u}
$$

并精确推进一周期

---

## 12. 当前最大的不确定性到底在哪

通过上面这套数学过程，现在最不确定的，不是公式本身，而是这一步：

$$
\mathcal{B}_{\text{gs}}
\longrightarrow
\mathcal{B}_{\text{aug}}
$$

也就是：

> 你构造出来的增广 kicked 子空间，到底是不是一个“足够好”的有限维工作空间？

当前代码给出的答案是启发式的，不是严格最优的。

所以后面的验证重点才是：

- fidelity
- exact benchmark
- 收敛性
- Step B 的正交化子空间测试

一句话总结：

**当前整条链路的数学是清楚的；真正不清楚的，是“我们给 kick 准备的那个有限维子空间，究竟是不是对的”。**

---

# Q&A: 针对当前方法的 5 个基础问题（2026-04-01）

## 问题 1：单个基函数里加了动量因子，有什么意义？总波函数最后又不是相乘，而是有限个基函数的线性组合

这是一个非常关键的问题。

### 先把两层结构分开

当前波函数有两层：

1. **单个基函数内部的结构**
   - 例如 `exp(-a x^2) * exp(i q x)`
   - 这里“乘上一个平面波”是在**单个基函数内部**发生的

2. **不同基函数之间的组合**
   - 总波函数是
     `Psi = sum_j u_j phi_j`
   - 这里是**线性组合**

这两件事并不矛盾。

### 为什么“线性组合”仍然能表示动量

因为我们并不是想让整个总波函数都“乘一个统一的动量相位”，而是想让它的可表示空间里，**出现若干个带不同动量中心的方向**。

你可以把它想成普通二维向量：

- 原来只有一个方向 `e_x`
- 你后来加了一个新方向 `e_y`

最终任何向量都还是

`v = a * e_x + b * e_y`

它当然不是“把整个向量乘一个 y”，但它已经**拥有了朝 y 方向延伸的能力**。

在这里也是一样：

- 原来的静态 Gaussian 基，主要张成“零动量附近”的子空间
- 加入 `exp(i q x)` 的 Gaussian 后，基空间里就多了一个“中心动量在 q 附近”的方向

所以“系统知道你加了个动量”这句话，更准确地说应该是：

> 你的有限维工作子空间里，多了一个能表示该动量分量的基向量方向

### 一个最简单的 1D 例子

设原来只有一个基函数：

`phi_0(x) = exp(-a x^2)`

它的动量分布中心在 `p = 0`。

现在加一个新基函数：

`phi_q(x) = exp(-a x^2) * exp(i q x)`

它在动量空间里的分布，不是还在 `p = 0`，而是整体平移到 `p = q` 附近。

于是如果总波函数写成

`Psi(x) = c_0 * phi_0(x) + c_1 * phi_q(x)`

那么这个 `Psi` 的动量分布就可以同时含有：

- `p ~ 0` 的成分
- `p ~ q` 的成分

这正是 kick 之后会发生的事情。

### 关键点

所以加动量基的意义不是：

> “让整个波函数统一获得一个动量”

而是：

> “让有限维基空间能够容纳 kicked 态中出现的不同动量分量”

这和 Fourier 展开、平面波展开本质上是一样的思路：

- 没有某个模式，就表示不了那个模式
- 有了那个模式，线性组合才有机会表示它

### 还要再加一句诚实的话

这并不自动保证“加了以后就一定表示得很好”。

它只说明：

- **如果不加这些带动量的方向，你几乎肯定表示不好**
- **加了以后，才至少有可能表示好**

所以“加动量基”是**必要方向**，但不是“自动成功”的充分条件。

---

## 问题 2：第五节里说“已经优化好的基态基”，这个到底是哪来的？是不是之前虚时间演化得到的？

是的，这里的“已经优化好的基态基”，指的就是 kicked real-time 开始之前那一整段静态准备流程的输出。

当前代码里，顺序是：

1. `svm_build_basis()`
   - 先构造一组适合基态的静态 Gaussian 基

2. `stochastic_refine()`
   - 在这组基上继续做随机优化

3. 虚时间 `evolution(...)`
   - 再对这组基做 imaginary-time polish

这三步结束以后，你手里得到的不是一个单独的波函数，而是两样东西：

1. **一组基函数**
   - 这就叫“基态基”

2. **在这组基函数上的最低能量系数 `u`**
   - 这组系数对应基态波函数

所以“基态基”这个词容易让人误会。

它不是说：

> 这些基函数本身就是基态

而是说：

> 这些基函数是经过静态优化后，比较适合表示基态的那一组 basis functions

然后才在这组 basis 的基础上去加动量基函数。

---

## 问题 3：width 候选是什么意思？为什么要这样做？我还是不理解你构造出来的动量基为什么有效

这个问题其实分成两层：

1. 为什么动量基除了“动量”之外还要管 width
2. 为什么这种 “Gaussian 包络 × 平面波” 的形式会有效

### 3.1 为什么不能只管动量，不管 width

因为一个波包不只有“中心动量”一个属性，它还有“包络的形状和宽度”。

例如

`exp(-a x^2) * exp(i q x)`

里面：

- `q` 决定动量中心在哪里
- `a` 决定包络有多宽

在动量空间里：

- `q` 决定“分布中心”
- `a` 决定“分布有多宽”

所以两个基函数即使有相同的动量中心 `q`：

- 如果 width 对了，它和真实 kicked 成分的 overlap 会比较大
- 如果 width 错得很厉害，它虽然“方向上像”，但还是抓不住真实波包

### 3.2 为什么 width 候选会借用基态基里的 `A_aa`

因为第一次 kick 最直接做的事情，其实是：

`psi(x) -> exp(-i kappa cos(2 k_L x)) * psi(x)`

这里是“乘一个相位因子”，并没有先把原来波函数的空间包络彻底改掉。

这意味着：

- kick 首先主要改变的是动量结构
- 原来基态波包的空间尺度，通常仍然是重要信息

所以把基态基里已经学到的 `A` 对角元拿来当 width 候选，是在表达一个想法：

> 如果基态基已经学到了系统在空间上的自然尺度，那么 kicked 后的新分量，很可能也还是围绕这些尺度展开

这不是严格证明，但它是很合理的启发式。

### 3.3 为什么还要再加一组人工给定的稀疏 width

因为只用基态里现成的 width 太死了。

kick 后实际需要的空间包络，未必完全等于基态那一套。

所以当前做法是两头都要：

- 一组人工给的稀疏 width：给一点灵活度
- 一组从基态学到的 width：保留系统原本的自然尺度

### 3.4 为什么 “Gaussian 包络 × 平面波” 会有效

因为它正好对应 kicked 后最典型的局部形状：

- 空间上还是一个局域波包
- 但动量中心已经被平移到 `0, ±2k_L, ±4k_L, ...`

换句话说，kick 后最自然出现的对象，并不是“纯平面波”，而更像是：

`局域包络 * 动量相位`

所以这个 ansatz 的有效性来自于：

> 它和 kicked 后局部波包的形状是相匹配的

但同样要强调：

这只能说明“形式上有道理”，不等于已经证明“数量上够了”。

---

## 问题 4：这个候选库会不会太少了

这个问题要分开回答。

### 从“生成数量”看

当前候选库在数量上其实不算特别少。

比如：

- `n_mom = 2` 时，当前一轮测试看到的是 `168` 个候选
- `n_mom = 4` 时，当前一轮测试看到的是 `528` 个候选

所以问题不完全是“库一开始就太小”。

### 从“最终有效留下来的基”看

真正的问题是：

> 候选库不少，但能通过 `cond(S)` 筛选真正留下来的基太少

例如之前的结果里，经常是：

- 候选几百个
- 真正留下的动量基只有十几个

于是对长期 kicked dynamics 来说，**有效库**是偏少的。

所以更准确的话应该是：

- **备选库本身未必太少**
- **但通过筛选后留下的有效子空间太小**

### 还要再补一句

如果 `n_mom` 取得很小，那么物理上候选库也可能确实不够。

因为这等于你只允许：

- `0`
- `±2k_L`
- `±4k_L`

之类有限几个动量谐波。

如果真实态已经扩散到更高动量，那就算数值上一个都不拒，你的候选库物理上还是截断得太狠。

所以这里有两层限制：

1. **物理截断限制**：`n_mom` 太小
2. **数值筛选限制**：`cond(S)` 让很多候选进不来

当前更突出的瓶颈是第 2 个。

---

## 问题 5：选好基底后，你又解了一次广义本征方程。你能保证它还在基态吗？还是说这回这个态已经是被踢后的态了？

这一步非常容易混。

答案是：

### 5.1 在“重新解广义本征方程”这一刻，它还不是 kicked 后的态

因为 kick 还没有施加。

这一刻你做的是：

1. 先把工作空间从
   - 原来的基态基 `B_gs`
   扩大到
   - 增广基 `B_aug`

2. 然后在新的更大空间里，再求一次 `H_0` 的最低本征态

所以这里算出来的是：

> **free Hamiltonian `H_0` 在增广基中的最低能量态**

按设计目标，它应该还是“同一个物理基态”，只是现在用更大的基空间来重新表示它。

### 5.2 为什么理论上它应该还是基态

因为你解的还是同一个 `H_0`：

- kinetic
- harmonic
- （gamma>0 时还有 interaction）

你并没有把 kick 加进这个本征值问题里。

所以求出来的最低本征向量，目标仍然是系统的基态，而不是 kicked 态。

### 5.3 能不能“绝对保证”它就是对的

严格说，**不能只靠一句话绝对保证**。

原因有两个：

1. 基是有限的、非正交的
2. 如果 `S` 病态，可能出现 ghost state 或数值污染

所以真正的判断方式只能是验证：

- 增广前后基态能量是否几乎不变
- 是否仍然接近已知基态能量
- 是否出现明显更低的非物理解

### 5.4 当前代码里为什么我们暂时认为这一步没出大问题

因为在 `gamma=0` 的当前测试里，增广后重新求出来的最低能量仍然接近基态值：

- `E_svm ~ 1`
- `E_refined ~ 1`
- `E_polished ~ 1`
- `E_aug ~ 1`

所以当前我们把它理解为：

> 这一步仍然是在同一个物理基态附近，只是换了一套更大的 basis 来表示

### 5.4.1 但这只是“弱验证”，不是强验证

这一点必须单独强调。

到目前为止，**我们并没有做足够强的数值验证，来严格证明“增广后重新解出来的确实还是同一个基态”**。

目前真正做过、并能算作验证的，只有比较弱的几项：

1. **能量检查**
   - `E_svm`
   - `E_refined`
   - `E_polished`
   - `E_aug`
   都接近已知的 gamma=0 基态能量 `1`

2. **Hamiltonian 没变**
   - 在重新解广义本征方程那一步，解的还是同一个 free Hamiltonian `H_0`
   - kick 还没有施加

这两点只能说明：

> “它看起来像是同一个基态的重新表示”

但还不能说明：

> “它已经被充分数值确认就是同一个物理态”

### 5.4.2 目前还没做、但应该做的强验证

如果要把这个问题说得更硬，我们至少还应该补下面几项数值检查：

#### (a) 增广前后态的 overlap / fidelity

要检查：

`| <psi_old | psi_aug> |^2`

是否接近 `1`。

这比单纯比较能量更强，因为两个不同态也可能碰巧有接近的能量。

#### (b) 基态可观测量是否一致

例如比较增广前后：

- `<x^2>`
- `<p^2>`
- 动能
- 势能

如果这些都几乎不变，才更能说明“它还是同一个物理基态”。

#### (c) 能量方差

要看增广后最低态的能量方差是不是仍然很小：

`sigma^2 = <H^2> - <H>^2`

如果方差很小，说明它确实接近某个真正的本征态。

#### (d) 不加 kick 时的静止性测试

一个非常直接的测试是：

1. 在增广基里重新求得基态
2. **不施加 kick**
3. 只做 `free_evolve_fixed_basis()`

理论上基态在自由演化下只会得到一个整体相位，不应改变任何物理可观测量。

所以如果这一步都不稳，那就说明“增广后还是基态”的说法站不住。

### 5.4.3 当前结论应该怎么表述才准确

所以关于这个问题，最准确的说法不是：

> 我们已经验证它就是基态

而应该是：

> 目前只有弱证据支持“它还是基态附近的同一物理态”，但还缺少更强的数值验证

这也是后续应该补上的一项基础检查。

### 5.5 kicked 态是什么时候才真正出现

真正 kicked 态是从这一步之后才出现的：

`u^+ = S^{-1} K u^-`

也就是调用 `apply_analytic_kick()` 的那一刻。

所以时间顺序要分清：

1. 增广基
2. 在增广基中重新求 `H_0` 基态
3. 然后才施加 kick

### 5.6 一个很重要的细节

就算增广后基态展开系数里，某些“动量基函数”的系数不是严格 0，也**不代表这个态已经被踢了**。

因为基是非正交的。

在非正交基里：

- “某个动量基的系数非零”
- 不等于“物理上就一定存在对应的独立动量成分”

它只说明：

> 在这组非正交 basis 的坐标表示里，这个最低本征态需要用这些 basis direction 一起拼出来

真正的 kicked 与否，还是要看：

- 它是不是 `H_0` 的最低本征态
- kick 是否已经施加
- 可观测量是否改变

## 一句话总总结

这 5 个问题可以压成一句话：

**当前“加动量基”的意义，不是说我们已经知道 kicked 子空间一定对；而是说我们先把物理上最可能需要的方向放进候选库，再通过重新表示基态、施加 kick、看 fidelity 和 exact benchmark，去检验这套子空间到底有没有用。**

---

# Q&A: 第 5 个问题的数值验证补充（2026-04-01）

前面关于“增广后重新解广义本征方程，是否仍然是基态”的回答，最初只给了**弱验证**和理论解释。

现在已经把更强的数值检查补进 `--kicked-gamma0` 的流程，并实际跑了一次。

## 补做了哪些检查

对“增广前后的基态是否还是同一个物理态”，目前补了 4 个检查：

1. **增广前后态的 cross-basis fidelity**
2. **增广前后关键可观测量比较**
   - `E`
   - `<x^2>`
   - `<p^2>`
3. **增广前后最低态的能量方差**
4. **不加 kick 的自由演化静止性测试**
   - 在增广基里只做 `free_evolve_fixed_basis()`
   - 看态是否只积累整体相位，而不改变物理量

## 当前默认参数下的结果

使用当前默认的增广参数：

- `n_mom = 4`
- `max_cond = 1e10`

运行得到：

```text
--- Ground-State Consistency Checks ---
Cross-basis fidelity (refined vs augmented ground state) = 1
Variance: refined = -2.66e-15, augmented = -2.11e-12
Refined moments:   E = 1, <x^2> = 1, <p^2> = 1
Augmented moments: E = 1, <x^2> = 1, <p^2> = 1
No-kick free test (5 periods): fidelity = 1,
    dE = -1.11e-15,
    d<x^2> = -9.25e-13,
    d<p^2> =  9.22e-13
```

## 这些结果说明什么

### 1. Cross-basis fidelity = 1

说明：

- 增广前重新表示出来的基态
- 增广后重新求得的最低态

在当前数值精度下几乎是同一个物理态。

这比“只看能量差不多”强很多。

### 2. `<x^2>` 和 `<p^2>` 完全一致

说明增广前后：

- 空间尺度没变
- 动量展宽没变

所以当前这一步并没有在“还没 kick 之前”就把态搞歪。

### 3. 能量方差接近 0

说明无论在 refined basis 还是 augmented basis 中，这个态都非常接近真正的本征态。

负号只是浮点舍入误差量级，不是物理负方差。

### 4. No-kick free test 基本不变

只做自由演化、完全不施加 kick 时：

- fidelity = 1
- `dE`, `d<x^2>`, `d<p^2>` 都在 `1e-12 ~ 1e-15` 量级

这正是“基态只积累整体相位”的表现。

## 现在对第 5 个问题应该怎么回答

现在可以比之前更硬一点地回答：

> 对于当前默认参数下的 gamma=0 增广基，重新解出来的最低态已经做过较强的数值验证，结果支持它仍然是原来同一个物理基态，而不是已经被 kick 后的态。

更准确地说：

- **kick 之前重新解本征方程得到的仍是 `H_0` 的基态**
- **真正的 kicked 态，是调用 `apply_analytic_kick()` 之后才出现**

## 还剩下什么没有验证

这并不代表整个 kicked dynamics 链路已经正确。

当前这些检查只能说明：

> “增广基这一步没有先把初态弄坏”

但后面的主要误差仍然在：

- kick 投影误差
- 长时间 kicked evolution 中子空间不够完备
- `S^{-1}` 对病态方向的放大

所以：

- 第 5 个问题现在基本补上了
- 但整个方法是否最终可靠，仍然要看 kick fidelity 和 exact benchmark

---

# Q&A: 为什么增广后还是基态？是不是新加的动量基系数都变成 0 了（2026-04-01）

不是，**不能把“新加基函数的系数是不是 0”当成判断依据**。

## 短答案

增广后仍然是基态，不是因为：

> “所有新加的动量基系数都变成 0 了”

而是因为：

> 在更大的基空间里重新解 `H_0` 的最低本征态后，得到的整个波函数和原来的基态是同一个物理态（至少在当前数值精度下）

这个判断依赖的是：

- cross-basis fidelity
- 基态可观测量
- 能量方差
- no-kick free test

而**不是**看某几个 basis coefficient。

## 为什么“系数是不是 0”不可靠

因为当前基是**非正交基**。

在非正交基里，总波函数写成：

`Psi = sum_j u_j phi_j`

这里的 `u_j` 不是“某个物理模式的独立振幅”，原因是：

- `phi_i` 和 `phi_j` 彼此有 overlap
- 一个基函数的那部分内容，可能已经部分包含在别的基函数张成的方向里

所以：

- 某个新基函数的系数不为 0
  不代表“物理上已经有一个独立的新动量成分”
- 某个系数恰好接近 0
  也不代表“这个方向就完全没参与”

## 一个简单类比

想象二维平面里有两组基：

1. 正交基：`e_x, e_y`
2. 非正交基：`v_1 = e_x`, `v_2 = e_x + 0.001 e_y`

同一个向量在这两组基下的坐标会差很多。

特别是在非正交基里，某个系数的大小很容易被“基向量之间互相重叠”影响，不能直接当作物理含义。

在当前 ECG 基里也是一样。

## 那为什么增广后还是基态

更本质的原因是：

### 1. 你解的是同一个 Hamiltonian

增广前后，重新解的都是 free Hamiltonian `H_0`：

- kinetic
- harmonic
- gamma>0 时再加 interaction

这里**没有 kick 项**。

所以从理论上讲，最低本征态的物理身份不应该变。

### 2. 增广只是“换了一个更大的表示空间”

你可以把它理解成：

- 原来用 10 个 basis function 表示基态
- 现在允许用 24 个 basis function 表示同一个态

如果新增方向对这个基态没有必要，它们可以通过非正交组合被吸收进整体表示里，而不改变最终物理态。

### 3. 当前已经做了强于“看系数”的数值验证

当前默认参数下，验证结果是：

- cross-basis fidelity = 1
- `E`, `<x^2>`, `<p^2>` 完全一致
- 能量方差仍接近 0
- no-kick free test 只积累整体相位

这些才是真正说明“它还是基态”的证据。

## 那新加的动量基系数到底会不会非零

**很可能会非零。**

这并不奇怪，也不危险。

原因有两个：

### 原因 1：非正交展开下，系数没有唯一物理解释

同一个物理态，在不同的非正交 basis 中，坐标可以变化很大。

所以“某个新加 basis 的系数非零”完全可能只是表示方式变化，而不是物理内容变化。

### 原因 2：增广基的最低态是在整个更大子空间里重新最优化得到的

最低本征向量会自动在这组更大的 basis 中寻找最优表示。

这个最优表示可能自然会用到一些新基函数，即使最终物理态仍然和旧基态几乎一样。

所以真正该问的不是：

> 新加基函数的系数是不是 0？

而是：

> 用这些系数拼出来的整个波函数，和原来是不是同一个物理态？

## 真正会让态离开基态的，不是“换基表示”，而是施加 kick

这点再强调一次：

1. 增广基
2. 在增广基中重新求 `H_0` 的最低态
3. 这时仍然是基态
4. 真正让它离开基态的是

`u^+ = S^{-1} K u^-`

也就是 `apply_analytic_kick()` 那一步。

所以：

- **增广基**：只是换一个更大的表示空间
- **kick 投影**：才是真正把态推离基态的动力学步骤

## 一句话总结

**增广后仍然是基态，不是因为“新加基的系数都变成 0 了”，而是因为虽然坐标表示变了，但整个波函数作为一个物理态并没有变；这件事要靠态的 fidelity 和可观测量去判断，而不是看单个非正交基函数的系数。**

---

# Q&A: 第一次 kick 的投影误差有多大？能不能在第一次 kick 后立刻优化（2026-04-01）

## 短答案

### 1. 当前已经量到的“第一次 kick 投影误差”

在当前默认参数下，代码直接量到的最直接指标是：

```text
first-kick fidelity = 0.998956
```

也就是说，按当前 `fidelity = norm_after / norm_before` 的定义，第一次 kick 的**norm 级别误差大约是**

```text
|1 - fidelity| ≈ 1.0e-3
```

约 `0.1%` 量级。

### 2. 但这不是完整的“态误差”

这个 `fidelity` 只告诉我们：

> 把 kicked 态投影回当前 basis 后，norm 扭曲有多大

它**不能完整告诉我们**：

> kicked 后的整个波函数形状到底错了多少

因为一个态可以：

- norm 几乎不变
- 但方向已经偏了不少

所以这个量是一个**弱指标**，不是最终答案。

### 3. 第一次 kick 后立刻优化：可以做，而且值得做

可以。

但它更准确的意义是：

> 看“开头这一步的 basis mismatch 能不能先压下去”

而不是：

> 一次优化就把整个长期 kicked dynamics 都解决

因为当前代码的主要崩坏，并不是只发生在第一次 kick，而是随着 kick 次数增加，态不断扩散到更高动量后，子空间越来越不够用。

所以：

- **第一次 kick 后优化** 可能明显改善前几步
- 但不一定能单独解决长期问题

## 1. 当前第一次 kick 的误差，已经知道多少

在当前默认参数下，`--kicked-gamma0 --n-kicks 1` 的结果是：

```text
Kick 1: fidelity = 0.998956
Kick 1 after one period: rel_error ≈ 7.1%
```

这里要分清两个量：

### (a) `fidelity = 0.998956`

这是**kick 投影本身**的弱误差指标。

表示：

- 如果只看 norm，第一次 kick 并没有炸
- 当前 basis 对第一次 kick 的主要分量，其实已经抓到不少

### (b) `rel_error ≈ 7.1%`

这是第一次“完整一步”之后的误差指标。

它包含了：

- kick 投影误差
- 后续自由演化中的表示误差
- 以及当前 finite basis 对 kicked 态整体表示的不足

所以它已经不是“单纯第一次 kick 投影误差”。

## 2. 为什么第一次 kick 看起来还行，但后面还是会崩

因为现在的问题更像是：

### 第一次 kick

主要只需要表示：

- 基态附近
- 再加上一些低阶动量谐波

这时候当前动量基库还能凑合覆盖。

### 多次 kick 之后

态会逐渐扩散到更高动量方向：

- ±2k_L
- ±4k_L
- ±6k_L
- ...

而当前有效留下来的动量基数量太少。

所以即使第一次 kick 的投影误差不大，后面还是会越来越差。

这就是为什么：

- **第一次 kick 优化是有意义的**
- 但**它不太可能单独解决整个长期问题**

## 3. 第一次 kick 后优化，数学上到底是在干什么

如果要做这件事，本质上是在做：

> 不再只让基态 basis 去“硬接”第一次 kick，而是让 basis 针对第一次 kicked state 再适配一次

也就是把目标从

```text
“把基态表示好”
```

变成

```text
“把第一次 kick 后的态也表示好”
```

### 当前流程

当前是：

1. 用静态流程优化基态基
2. 加动量候选
3. 重新解基态
4. 直接开始 kick

### 可以尝试的新流程

可以改成：

1. 用静态流程优化基态基
2. 加动量候选
3. 重新解基态
4. **先做一次 kick**
5. **针对这次 kick 后的态，再优化 basis**
6. 然后再开始正式 kicked evolution

这个思路可以叫：

> “one-kick adapted basis”  
> 或者  
> “first-kick bootstrapping”

## 4. 怎么优化，才不是拍脑袋

这一步最好不要只说“再调调参数”，而要给一个明确目标函数。

### 方案 A：gamma=0 下，用 exact 一步 kicked state 做目标

这是最强、也最干净的方案。

因为 gamma=0 已经有 exact reference，所以可以：

1. 从 exact solver 得到第一次 kick 后的精确态 `psi_exact^(1)`
2. 在 ECG basis 中构造一个 trial state `psi_basis(z)`
3. 调整新增动量基的参数，使

```text
F = |<psi_exact^(1) | psi_basis(z)>|^2
```

最大

或者等价地最小化

```text
1 - F
```

这就把“第一次 kick 后优化”变成了一个非常明确的拟合问题。

### 方案 B：没有 exact 时，用 kick 投影质量做目标

如果没有 exact，可以定义一个弱一点但仍有意义的目标：

- 第一次 kick 的 norm fidelity 尽量接近 1
- 同时惩罚过大的 `cond(S)`

例如目标函数可以是

```text
J = |1 - fidelity_1| + lambda * penalty(cond(S))
```

这样做虽然比 exact overlap 弱，但至少比“手调参数”更有原则。

### 方案 C：第一次 kick 后做一次 basis augmentation / refine

更工程化的做法是：

1. 先做第一次 kick
2. 看哪些动量方向最缺
3. 只为这一步新加一些 basis
4. 然后再重新求投影

这相当于做一次“针对第一次 kicked state 的局部 basis repair”

## 5. 我对这件事的判断

我的判断是：

### 值得做

因为它可以直接回答一个很实际的问题：

> 当前的问题是不是在一开始就已经埋下了？

如果第一次 kick 后做一次定向优化，前几步精度显著提高，那说明：

- 初始 kicked 子空间的准备确实不够好
- 从开头修一下是有价值的

### 但不要指望它单独解决最终问题

因为当前长期误差的根源更深：

- 态会随着 kick 次数继续扩散
- 有效 basis 仍然会逐渐不够

所以最可能的结果是：

- **第一次 kick 优化** 能改善 early-time behavior
- **长期是否还行** 仍然要靠更强的子空间构造办法

## 6. 现在最推荐怎么做

如果真的要试，我建议按这个顺序：

### Step 1：先把“第一次 kick 后的精确态误差”量出来

也就是不要只看 `norm fidelity`，而是补一个更强的 gamma=0 指标：

- 第一次 kick 后 ECG 态和 exact 态的 overlap

先知道“第一次 kick 的真实态误差到底有多大”。

### Step 2：做一次 one-kick adapted basis 原型

先只针对第一次 kick 后的态做一次 basis 适配，看：

- 第一次 kick overlap 有没有显著提升
- 接下来前几步误差有没有明显下降

### Step 3：再决定要不要推广成更一般的动态适配

如果第一步和第二步都显示效果明显，再考虑：

- 每隔几步做一次 re-adaptation
- 或者直接进入 Step B 的正交化子空间方案

## 一句话总结

**第一次 kick 后做一次定向优化是值得试的，它能帮我们判断“问题是不是从第一步就开始了”；但它更像是 early-time repair，而不是长期 kicked dynamics 的最终解法。**

---

# Q&A: 现在的问题是不是因为一开始动量基截断太小？每次 kick 前会重新加动量基吗（2026-04-01）

这个问题里，其实混了 3 个不同的量，必须先分开：

## 1. 这里有哪几个不同的 “N”

### (a) 粒子数 `N = 2`

这是系统本身的粒子数。

当前 `--kicked-gamma0` 测试里，确实固定是：

- `N = 2` 粒子

这不是“动量截断”，只是物理系统大小。

### (b) 动量截断 `n_mom`

这才是“动量基候选库截到几阶”的参数。

例如：

- `n_mom = 2` 表示只考虑 `m = -2,-1,0,1,2`
- 对 `N=2` 粒子来说，就是所有 `(m1,m2)` 的组合

这会决定：

- 候选库里有哪些动量对
- 也就是你允许 basis 去表示到多高的离散动量谐波

### (c) Bessel 截断 `n_bessel`

这又是另一回事。

在解析 kick 矩阵里，`exp(-i kappa cos)` 用 Bessel 展开时，也有一个和式截断。

当前 `apply_analytic_kick()` 默认是：

- `n_bessel = 20`

这和 `n_mom` 不是一回事。

## 2. 如果你说的是“旧代码里动量截断只取到 2”

那这个怀疑是有道理的，但只能算**部分原因**，不是完整原因。

### 为什么它有道理

如果 `n_mom = 2`，那你一开始允许的动量谐波就只有：

- `0`
- `±2k_L`
- `±4k_L`

对于多次 kick 之后逐渐扩散的态，这个库物理上可能确实太窄。

### 为什么它不是完整原因

因为当前我们已经把默认值从旧的 `n_mom = 2` 提到 `n_mom = 4` 了。

在这个情况下：

- 候选库已经有 `528` 个候选
- 但最后真正保留下来的动量基只有 `14` 个左右

所以当前更尖锐的瓶颈不是：

> “一开始候选库完全没生成出来”

而是：

> “候选库生成出来了，但因为 overlap / cond(S) 问题，真正能留下来的有效 basis 太少”

换句话说：

- `n_mom` 太小，确实可能带来物理截断误差
- 但当前更主要的问题，是**贪心筛选之后有效子空间太小**

## 3. 当前代码会不会在每次 kick 前重新加一次动量基

**不会。当前代码只加一次。**

这是现在流程里非常关键的一点。

### 当前真实流程

现在的顺序是：

1. 准备基态基
2. 调一次 `augment_basis_with_momentum(...)`
3. 得到一个固定的增广基
4. 后面所有 kick 都在这同一个固定基里做

也就是说，后面的循环里：

- 会更新 `u`
- 不会重新生成新的动量基
- 也不会在每次 kick 前再调用一次动量增广

### 后面的循环里到底在变什么

后面的每一步只做：

1. `apply_analytic_kick(...)`
   - 更新线性系数 `u`

2. `free_evolve_fixed_basis(...)`
   - 还是只更新线性系数 `u`

而下面这些都保持不变：

- basis 的个数
- 每个基函数的 `A`
- 每个基函数的 `B`
- 每个基函数的 `R`

所以当前方法本质上是：

> **一次性准备一个固定的 kicked 子空间，然后所有 kick 都在这个固定子空间里传播**

## 4. 这会带来什么后果

这正是当前长期误差变大的重要原因之一。

### 第一次、前几次 kick

当前固定基还能大致覆盖主要分量，所以：

- first-kick fidelity 往往还不错
- 前几步误差也未必立刻爆炸

### 后面的 kick

随着态不断扩散到更高动量方向：

- 当前固定子空间里没有新的 basis direction 被加入
- 所以它只能继续用原来那几十个 basis 硬撑

结果就是：

- kicked 态越来越多的部分落到当前子空间外面
- 投影误差逐步累积
- fidelity 和能量误差开始恶化

所以你这个直觉本身是对的：

> 问题确实和“初始动量子空间准备得不够”有关

但更准确的话应该是：

> 不是简单地说“因为当初只取了 n_mom=2”，而是“当前方法只在一开始做一次有限动量增广，之后不再自适应更新子空间，所以随着 kick 次数增加，固定子空间越来越不够用”

## 5. 如果每次 kick 前都重新加动量基，会不会更好

从原理上说，**很可能会更好**。

因为这就变成了：

> basis 跟着当前态一起扩展

而不是：

> 用第一次准备好的 basis 去硬撑全部后续动力学

这类思路可以叫：

- 动态基扩展
- adaptive basis update
- kick-by-kick basis repair

### 但代价也很明显

如果每次 kick 前都重新加基：

1. 基的大小会越来越大
2. `S` 的病态问题会更严重
3. 你还要处理“不同 kick 之间 basis 改了，态怎么转移”这个问题

所以它不是白来的改进。

## 6. 当前最准确的结论

可以把这个问题总结成三句话：

1. **当前代码不会在每次 kick 前重新加动量基，只在开始时加一次**
2. **固定子空间不再更新，确实是长期误差变大的重要原因之一**
3. **旧代码里 `n_mom=2` 可能太小，但当前更深层的问题是：即使候选库更大，最后有效保留下来的 basis 仍然太少，而且之后不再自适应更新**

## 一句话总结

**当前方法不是“每次 kick 前都重新准备动量基”，而是“一开始准备一次固定动量子空间，然后一路用到底”；这正是为什么它在前几步还能工作、但后面越来越吃力。**

---

# Q&A: 后面的每一步是不是“直接乘算子，然后投影”（2026-04-01）

可以这么理解，但要分成两半说：

## 短答案

### 对 kick 这一步

**是的。**

当前固定基方法里，每次 kick 的核心就是：

1. 先让 kick 算子作用在当前态上
2. 再把结果投影回当前固定增广基

数学上就是：

`S u_after = K u_before`

也就是：

`u_after = S^(-1) K u_before`

其中：

- `S_ij = <phi_i | phi_j>`
- `K_ij = <phi_i | U_kick | phi_j>`

所以就 kick 而言，你说的“直接乘算子，然后投影”是对的。

### 但对整个一个周期

**不完全是。**

因为每个周期里还有第二步：

- 两次 kick 之间的自由演化 `free_evolve_fixed_basis()`

这一步不是“再乘一次 kick 算子然后再投影”，而是：

- 在同一个固定基里
- 对 `i S du/dt = H u`
- 做广义本征模下的精确传播

所以一个完整周期更准确地写是：

1. **kick + 投影**
2. **fixed-basis free evolution**

## 当前代码里后面每一步到底有没有重新加 basis

没有。

当前后续循环里：

- 只更新 `u`
- 不更新 basis 的 `A/B/R`
- 不重新调用 `augment_basis_with_momentum()`

所以它不是：

> 每次 kick 前重新准备一个新 basis

而是：

> 用同一个固定增广基，重复执行“kick 投影 + 自由传播”

## 为什么这件事重要

因为这正是当前方法误差的来源之一。

如果某一次 kick 后，真实态已经有很大一部分跑到了当前子空间外面，那么：

- `K u_before` 代表“kick 后的理想结果”
- 但我们没有无限维空间去承接它
- 所以只能解 `S u_after = K u_before`，把它压回当前有限子空间

这一压回去，就会有投影误差。

而且当前 basis 不会自适应更新，所以这个误差会随着 kick 次数继续累积。

## 一句话总结

**对，当前后面的 kick 步骤本质上就是“算子作用 + 投影回固定基”；但完整一个周期还包含一次固定基中的精确自由演化，而且不会在每一步重新加新的动量基。**

---

# Q&A: 目前可考虑的解决办法与取舍（2026-04-01）

下面这些不是互斥关系，而是几条不同层级的路线。

## 方案 1：先做 one-kick adapted basis

### 核心想法

不要直接拿“基态基 + 一次性动量增广”去承接全部后续动力学，而是在**第一次 kick 之后**，专门针对第一次 kicked state 再做一次 basis 适配。

### 适合解决什么

- 第一脚就已经存在的 basis mismatch
- 早期 few-kick 误差过大

### 优点

- 实现成本相对低
- 物理直觉清楚
- 很适合 gamma=0，因为有 exact 一步结果可以做目标

### 缺点

- 更像 early-time repair
- 不一定能解决长期 kicked dynamics 的问题

## 方案 2：做正交化子空间（Step B）

### 核心想法

不是贪心地一个个加 raw Gaussian，而是：

1. 先生成完整候选库
2. 构造整个候选库的 overlap 矩阵
3. 做 SVD / 本征分解
4. 丢掉太小奇异值方向
5. 在正交化后的 reduced subspace 里做 kick 和 free evolution

### 适合解决什么

- 贪心选基太差
- `S^{-1}` 病态放大

### 优点

- 信息价值最高
- 最能区分“选基算法问题”还是“Gaussian 子空间本身不行”
- 不需要一开始就改成全动态方法

### 缺点

- 实现复杂度中等
- 得重新组织当前 `raw basis -> overlap -> propagation` 的数值流程

## 方案 3：动态 basis expansion

### 核心想法

不是只在开始时增广一次 basis，而是：

- 每次 kick 后
- 或每隔几次 kick
- 根据当前态的误差指标再加一些新基

### 适合解决什么

- 长时间动力学里动量支撑持续外扩

### 优点

- 最贴合当前误差的来源
- 原理上最有希望改善长期问题

### 缺点

- basis 会越来越大
- `S` 会越来越病态
- 需要设计增长、截断、甚至 pruning 机制

## 方案 4：残差驱动的加基，而不是盲目枚举

### 核心想法

当前是“先枚举一大堆候选，再看 cond(S)”。更聪明的做法是直接围绕 kick residual 来选：

`r = U_kick psi - P_basis(U_kick psi)`

优先加入那些能最大幅度降低残差的候选基函数。

### 优点

- 更有针对性
- 比纯几何上的候选池筛选更接近真正需要的方向

### 缺点

- 需要定义和近似计算 residual
- 实现比当前贪心法复杂

## 方案 5：Rothe / 每步重新优化非线性参数

### 核心想法

不是固定 `A/B/R`，而是在动力学中不断重新优化非线性参数，让 basis 跟着态走。

### 优点

- 理论上最强
- 更接近“变分 manifold 跟踪真实态”

### 缺点

- 实现和调试成本最高
- 很容易掉进新的稳定性问题
- 现在直接上这条，风险偏大

## 方案 6：对 gamma=0 停止使用 ECG，当作基准线

### 核心想法

承认 gamma=0 本来就有 exact 网格解，因此：

- gamma=0 永远用 exact
- ECG 只服务于 gamma>0

### 优点

- 最稳妥
- 工程上最干净

### 缺点

- 解决的是“项目路线”，不是“ECG kicked dynamics 本身的问题”
- 不能回答“ECG 到底能不能处理 kicked dynamics”

## 当前我最推荐的顺序

如果要按性价比排序，我会建议：

1. **先补 first-kick stronger diagnostic**
   - 不只看 norm fidelity，还要看第一脚后的更强态误差

2. **做方案 1：one-kick adapted basis**
   - 看 early-time 是否明显改善

3. **做方案 2：正交化子空间**
   - 判断问题到底出在贪心选基，还是 Gaussian 子空间本身

4. **如果长期还是不行，再考虑方案 3 或方案 5**

## 一句话总结

最务实的路线不是一上来就做最重的方法，而是：

**先做 first-kick 诊断，再做 one-kick repair，再做正交化子空间判别。**

---

# Q&A: 第一次 kick 后那一瞬间的能量是多少？第二次 kick 的动量基该怎么继续加（2026-04-01）

## 1. 第一次 kick 后“立刻”的能量是多少

在当前默认参数下：

- `kappa = 1`
- `k_L = 0.5`
- `N = 2`
- 初态是 `H_0` 的基态

当前可以把 `kick 1` 的能量理解为“第一次 kick 后那一瞬间的能量”。

### 为什么可以这么理解

因为：

1. 对 exact 参考解来说，第一次 kick 之前那一段 free evolution 只是基态积累整体相位，不改变 `H_0` 能量
2. 对当前 ECG 固定基传播来说，kick 之后再做 `free_evolve_fixed_basis()`，`H_0` 能量也是守恒的

所以：

> 第一次 kick 后立刻的 `H_0` 能量  
> = 代码里 `kick 1` 打印出来的能量

### 当前数值结果

在默认参数下：

- **Exact:** `E_exact^(1) = 1.315334`
- **ECG:**  `E_ECG^(1)   = 1.4088286`

所以第一次 kick 后的误差大约是：

- 绝对误差 `~ 0.0935`
- 相对误差 `~ 7.1%`

### 一个解析上的参考值

对连续单粒子谐振子基态，可以直接算出：

单粒子第一次 kick 后的能量增量大约是：

`Delta E_1 = kappa^2 * k_L^2 * (1 - exp(-4*k_L^2))`

代入 `kappa=1, k_L=0.5`：

`Delta E_1 ≈ 0.15803`

单粒子：

`E_1^(after kick) ≈ 0.5 + 0.15803 = 0.65803`

两粒子非相互作用总能量：

`E_2^(after kick) ≈ 1.31606`

这和 exact 网格算出来的 `1.315334` 是一致的。

## 2. 为什么第二次 kick 后，动量表示不能还只停在原来那些模式

这个直觉是完全对的。

当前 kick 算子按 Bessel 展开是：

`U_kick = sum_n c_n exp(i * 2 * n * k_L * x)`

其中：

`c_n = (-i)^n J_n(kappa)`

### 如果当前一个基函数带有动量标签 `m`

例如它内部有因子：

`exp(i * 2 * m * k_L * x)`

那么再施加一次 kick 后，会耦合到：

`exp(i * 2 * (m + n) * k_L * x)`

也就是说：

> kick 会把原来的动量标签 `m` 推到 `m+n`

### 对两粒子就是向量相加

如果当前模式标签是：

`m = (m1, m2)`

再施加一次 kick 后，就会耦合到：

`m' = (m1 + n1, m2 + n2)`

其中：

- `n1`
- `n2`

来自单粒子 kick 的显著 Bessel 谐波。

所以第二次 kick 后的 basis 支撑，绝对不能还只看第一次的那些模式。

## 3. 那“应该怎么加”才对

这里最自然的规则其实很清楚：

### 规则：按 kick 的卷积闭包去扩张

如果当前 basis 中已经有一个动量标签：

`m = (m1, m2)`

而 kick 的显著谐波集合是：

`N_sig = { n : |J_n(kappa)| > eps }`

那么下一轮应该考虑加入：

`m' = (m1 + n1, m2 + n2)`

其中：

- `n1 in N_sig`
- `n2 in N_sig`

这就是最直接、最物理的“第二次 kick 应该怎么补 basis”规则。

## 4. 具体可以有哪几种加法方案

### 方案 A：静态闭包法

一开始就把未来若干次 kick 可能需要的模式都准备进去。

例如：

- 如果认为显著谐波大约到 `|n| <= n_sig`
- 想支持 `P` 次 kick

那么单粒子动量范围至少准备到：

`|m| <= P * n_sig`

两粒子就准备所有：

`(m1, m2),  |m1|, |m2| <= P * n_sig`

#### 优点

- 最简单
- 不需要中途改 basis

#### 缺点

- basis 会非常大
- overlap 病态会更严重

### 方案 B：frontier expansion（我更推荐）

每次 kick 之后，只从当前“已经显著占据”的模式出发，向外扩一圈。

也就是：

1. 找当前系数里重要的动量标签 `m`
2. 对每个 `m`，加上所有 `n in N_sig`
3. 生成它们的“孩子模式”
4. 只把最重要的一部分新模式加入 basis

#### 优点

- 比一开始全塞进去节省很多
- 更贴近当前态真正需要的方向

#### 缺点

- basis 会随时间变化
- 需要处理不同 kick 之间的 basis 转移

### 方案 C：残差驱动加法

不直接按动量标签扩一圈，而是先计算：

`r = U_kick psi - P_basis(U_kick psi)`

然后只加那些最能降低这个残差的 basis 方向。

#### 优点

- 最有针对性

#### 缺点

- 实现更复杂

## 5. 当前代码属于哪一种

当前代码其实属于：

> “方案 A 的一个非常弱的、只做一次的版本”

也就是：

- 一开始用 `n_mom` 生成一个有限候选库
- 筛选一次
- 后面再也不扩

所以它的问题不是难理解，而是它确实太静态了。

## 6. 我对下一步实现的建议

如果要往前推，我建议别一上来就做最复杂的动态自适应。

更合理的顺序是：

### Step 1：先把“显著谐波集合”定出来

例如对 `kappa=1`，按 Bessel 权重取：

- `n = 0, ±1, ±2` 为主
- `±3` 视阈值决定要不要纳入

先明确一个 `N_sig`。

### Step 2：做一个最小版 frontier expansion

只做：

- 第一次 kick 后
- 根据当前重要模式
- 扩一圈 `m -> m + n`

也就是只做“一次动态补基”，先看前几步能不能显著改善。

### Step 3：再决定要不要推广成每几步扩一次

如果这一步明显有效，再走向更一般的 dynamic basis expansion。

## 一句话总结

**第一次 kick 后的能量，当前默认参数下 exact 是 `1.315334`、ECG 是 `1.4088286`；而第二次 kick 之后的动量基，原则上应该按 `m -> m+n` 的卷积闭包规则继续往外加，而不是一直死守第一次准备好的那一小圈模式。**

### 第一次 kick 后“瞬时”误差的实际数值（当前默认参数）

现在已经把“第一次 kick 后、还没做自由演化”的瞬时诊断补进代码，当前默认参数下结果是：

```text
Projection fidelity = 0.9989555239
|1 - fidelity|      = 1.0444760547e-03

E after kick:
  ECG   = 1.408828597
  exact = 1.315334023
  dE    = 9.3494573850e-02

<x^2> after kick:
  ECG   = 1.027879977
  exact = 1
  dx2   = 2.7879977035e-02

<p^2> after kick:
  ECG   = 1.789777216
  exact = 1.630668046
  dp2   = 1.5910917067e-01
```

这组结果说明：

1. **norm 级别的投影误差不大**
   - `|1-fidelity| ~ 1e-3`

2. **但物理量误差已经不算小**
   - `dE ~ 9.35e-2`
   - `dx2 ~ 2.79e-2`
   - `dp2 ~ 1.59e-1`

3. 特别是 `x^2` 的误差很说明问题
   - 对 exact kick 而言，kick 只是位置依赖相位，`|psi(x)|^2` 不变
   - 所以 exact 的 `<x^2>` 在 kick 瞬间应保持不变
   - 但 ECG 投影后出现了 `dx2 ~ 2.8e-2`
   - 这意味着投影误差不只是“norm 稍微不准”，而是已经改变了位置分布相关的物理量

也就是说：

> **第一次 kick 的问题并不是单纯一个很小的数值误差，而是投影回有限基后，态本身的物理形状已经出现了可见偏差。**

## 那这个误差到底从哪来

这个问题的答案其实比“basis 不够大”更细一点。

当前第一次 kick 的瞬时误差，主要来自 3 层：

### 来源 1：kicked 后的真实态不在当前有限子空间里

第一次 kick 后的真实态是：

`psi_exact^(1) = U_kick psi_0`

而我们当前只能在固定增广基张成的有限维子空间 `V_aug` 里表示态。

所以实际上发生的是：

`psi_exact^(1) = psi_parallel + psi_perp`

其中：

- `psi_parallel` 在当前子空间里
- `psi_perp` 在当前子空间外

投影步骤做的，只是把 `psi_exact^(1)` 压回 `psi_parallel`。

所以只要 `psi_perp != 0`，就一定有投影误差。

### 来源 2：当前动量基并不是“第一次 kick 后基态的精确闭包”

这一点非常关键，而且往往容易被忽略。

第一次 kick 本质上是：

`psi_0(x) -> exp(-i kappa cos(2 k_L x)) * psi_0(x)`

也就是：

- 原来的空间包络还在
- 只是乘上了一个位置依赖相位

如果我们真的想让“第一次 kick 后的态”尽量被 basis 精确承接，那么最自然的 basis 应该接近：

`exp(i 2 n k_L x) * psi_0(x)`

或者更细一点说，是：

`exp(i 2 n k_L x) * phi_j(x)`

这种“对原有基函数逐个做动量平移”的闭包。

但当前代码不是这样构造的。

当前 `augment_basis_with_momentum()` 做的是：

- 从一组 width 候选里取一个 `w`
- 设定一个固定风格的 `A_new`
- 设定一个固定 `B_new`
- 再把动量标签塞进 `R_new`

这是一种**启发式动量包络基**，不是“原基函数在 kick 下的精确 image”。

所以即使只看第一次 kick，也会有 mismatch：

- 真实 kicked 态想要的是“原始空间包络 + 新相位”
- 当前 basis 给的是“一组近似相位化的 Gaussian 候选”

这就是为什么第一次 kick 后，连 `<x^2>` 都已经被投影误差改掉了。

### 来源 3：overlap 矩阵病态，`S^{-1}` 会放大误差

即使物理子空间已经差不多，如果 `S` 很病态：

- 求 `u_after = S^{-1} K u_before`
- 这一步也会放大小方向上的数值误差

当前默认参数下：

- `cond(S) ~ 1e10`

这已经非常大了。

所以第一次 kick 的误差并不是纯粹的“物理截断误差”，还混着一部分数值放大误差。

## 为什么 `<x^2>` 的误差特别能说明问题

因为 exact kick 是：

`psi(x) -> exp(-i phase(x)) psi(x)`

它只改相位，不改模方：

`|psi_after(x)|^2 = |psi_before(x)|^2`

所以在 exact 情况下，kick 那一瞬间：

- 所有只依赖 `|psi(x)|^2` 的位置空间观测量都不该变

例如：

- `<x^2>`
- 更一般地任何 `f(x)` 的期望值

但我们现在看到：

- exact: `<x^2> = 1`
- ECG projection 后: `<x^2> = 1.02788`

这说明：

> 当前投影不是只把“相位结构”近似了一下，而是已经把位置空间的振幅分布也扭曲了

这正是“basis 不能很好承接第一次 kicked state”的直接证据。

## 所以一句话回答“误差从哪来”

**第一次 kick 的误差，来自“真实 kicked 态不在当前有限子空间里”这一物理截断误差；而这个误差之所以在第一脚就已经明显，是因为当前动量基并不是第一次 kicked image 的精确闭包，再叠加病态 `S^{-1}` 的数值放大，于是连 `<x^2>` 这样的相位不变量都被投影扭曲了。**

---

# Q&A: 既然 `U_kick^dagger x^2 U_kick = x^2`，为什么现在 `<x^2>` 还是变了？这个误差怎么克服、怎么优化（2026-04-01）

## 1. 你这个判断是对的：在 exact 情况下，`<x^2>` 的确不该变

对单粒子 kick：

`U_kick = exp[-i kappa cos(2 k_L x)]`

这是一个只依赖位置算符 `x` 的函数。

所以：

- `U_kick` 和 `x` 对易
- `U_kick` 和任何 `f(x)` 都对易

特别地：

`[U_kick, x^2] = 0`

因此：

`U_kick^dagger x^2 U_kick = x^2`

所以对于 exact kicked state：

`<x^2>_after = <psi| U_kick^dagger x^2 U_kick |psi> = <psi|x^2|psi>`

这点完全正确。

## 2. 那为什么当前 ECG 投影后 `<x^2>` 还是变了

因为我们真正算的并不是 exact 态：

`psi_exact_after = U_kick psi`

而是**投影回有限 basis 后的近似态**：

`psi_proj_after = P_basis( U_kick psi )`

或者更准确地说，是非正交基下的 Galerkin 投影对应的那个有限维表示。

### 关键点

即使 exact 的 `U_kick` 与 `x^2` 对易，

**投影算符** `P_basis` 一般并不与 `x^2` 对易，也不与 `U_kick` 对易。

所以：

`P_basis U_kick psi`

和

`U_kick psi`

不是同一个态。

于是就不再能写成：

`<psi_proj_after| x^2 |psi_proj_after> = <psi|x^2|psi>`

### 更直白一点的写法

把 exact kicked state 拆成：

`U_kick psi = psi_parallel + psi_perp`

其中：

- `psi_parallel` 在当前 basis 子空间里
- `psi_perp` 在当前 basis 子空间外

如果我们真的保留的是 exact 态，那么 `<x^2>` 不变。

但现在实际保留的是：

- `psi_parallel`
- 甚至还可能带有非正交投影和再归一化带来的额外扭曲

那当然会改掉：

- norm
- 干涉项
- 位置分布

于是 `<x^2>` 就会变。

所以：

> 问题不在于 `U_kick` 不是幺正  
> 问题在于我们没有真的保留 `U_kick psi`，而是把它压回了一个不够好的有限维子空间

## 3. 这也解释了为什么只看 norm fidelity 不够

因为 `fidelity = norm_after / norm_before` 只能说明：

- 投影后的态 norm 扭曲了多少

但它不能说明：

- 位置分布是不是还对
- 动量分布是不是还对
- 相位结构和干涉是不是已经错了

所以现在第一次 kick 的结果才会出现：

- `|1 - fidelity| ~ 1e-3`
- 但 `<x^2>` 已经错了 `~ 2.8e-2`

这说明：

> 即使 norm 看起来还不错，态的物理形状也可能已经被投影扭曲了

## 4. 那这个误差要怎么克服

核心思路不是“让 `U_kick` 更幺正”，因为它本来就幺正。

真正要做的是：

> **让 basis 更擅长承接 `U_kick psi` 这个 kicked state**

也就是说，要优化的不是 kick 算子本身，而是：

- basis 的构造
- basis 的选择
- 投影的目标

## 5. 最应该优化什么

从这个角度看，当前最值得优化的目标函数其实很明确：

### 目标 A：直接最大化第一次 kick 后的态 overlap

如果是 gamma=0，有 exact 参考态：

`psi_exact^(1) = U_kick psi_0`

那么最合理的优化目标就是：

`F_1 = | <psi_exact^(1) | psi_basis^(1)> |^2`

最大化 `F_1`。

这比只看：

- norm fidelity
- 能量误差

都更直接，因为它盯住的是整个 kicked state 本身。

### 目标 B：至少要把“应当不变的量”守住

如果暂时不想直接做 overlap 优化，那么至少应该要求：

- 第一次 kick 后 `<x^2>` 尽量不变
- 其它纯位置分布量尽量不变

因为这是 exact kick 的硬约束。

所以可以把目标函数写成类似：

`J = w1 * |1 - F_1| + w2 * |<x^2>_basis - <x^2>_exact| + w3 * penalty(cond(S))`

这里：

- 第一项管态的整体 fidelity
- 第二项强制位置分布别被投影扭曲
- 第三项防止病态 basis

## 6. 具体应该优化 basis 的什么部分

我觉得可以分三层，从低成本到高成本：

### 方案 1：先优化“候选库的生成方式”

当前动量基不是原基函数的 exact kicked image，只是启发式的 `A/B/R` 候选。

更自然的候选应该是：

`phi_j^(n) = exp(i 2 n k_L x) * phi_j`

或者两粒子版本：

`phi_j^(n1,n2) = exp(i 2 n1 k_L x1) exp(i 2 n2 k_L x2) * phi_j`

这类候选有一个非常重要的优点：

- 它们和第一次 kick 后的真实态结构更贴近
- 因为第一次 kick 后的态本来就是这类对象的叠加

也就是说：

> 不再是“重新发明一组动量 Gaussian”
> 而是“直接取原 basis 在 kick 下最自然出现的 image”

这是我觉得最应该优先尝试的优化。

### 方案 2：做 one-kick adapted basis

流程改成：

1. 准备基态基
2. 生成 first-kick image 候选库
3. 用第一次 kicked state 作为目标
4. 选出最能承接第一次 kick 的 basis

也就是让 basis 不是只为基态优化，而是也为第一次 kicked state 优化。

### 方案 3：在这些候选上做正交化 / 截断

因为即使候选本身更对了，如果 `S` 还是很病态，`S^{-1}` 还是会把误差放大。

所以很可能需要：

1. 先生成更物理的 kicked-image 候选
2. 再对这些候选做正交化 / SVD 截断
3. 在这个 reduced subspace 里传播

## 7. 我现在最推荐的优化路线

如果按“最有希望直接改善第一次 kick 的 `<x^2>` 误差”的顺序，我会这样排：

### 第一步

先不要再用现在这种纯启发式动量基作为唯一候选。

改成加入：

`exp(i 2 n k_L x) * phi_j`

这类 **first-kick image basis**。

### 第二步

对这批候选，直接用第一次 kick 后的 exact gamma=0 态做目标，比较：

- overlap
- `<x^2>`
- `<p^2>`
- `E`

### 第三步

如果候选更多后 `S` 太病态，再引入正交化子空间。

## 一句话总结

**`U_kick` 本身没有问题；问题在于我们把 `U_kick psi` 投影回了一个不够合适的有限 basis，所以连本来应当不变的 `<x^2>` 都被扭曲了。要克服这个误差，最该优化的不是幺正性，而是“第一次 kicked state 的 basis 承接能力”，最自然的方向就是引入 first-kick image basis 并围绕第一次 kick 态本身做优化。**

---

# Q&A: kick 之后到底是“基底在变”，还是“同一组基底里的系数在变”（2026-04-01）

这个问题非常重要，因为它决定我们到底在讨论：

- **固定基方法**
还是
- **动态基方法**

## 短答案

### 你说对的部分

你说得对的一点是：

> kick 之后，系统的“主导动量内容”确实变了

也就是说，波函数本身的物理状态在变。

原来主要是：

- 零动量附近

kick 之后会逐步变成：

- 0
- ±2k_L
- ±4k_L
- ...

这些分量的叠加。

### 你需要修正的部分

但在**当前代码实现**里，这不意味着：

> basis 本身在每一步自动变化了

当前实现里，**basis 是固定的**。

不变的是：

- basis 的个数
- 每个基函数的 `A`
- 每个基函数的 `B`
- 每个基函数的 `R`

变的是：

- 线性系数 `u`

所以更准确的话是：

> **物理态在变化，但承载这个态的 basis 本身没有变化；变化的是这个态在固定 basis 中的坐标表示。**

## 为什么看起来像“基底在变”

因为当前 basis 不是普通的坐标轴，而是一组自带物理结构的函数：

- 有的 basis 函数偏零动量
- 有的 basis 函数偏 `+2k_L`
- 有的 basis 函数偏 `-2k_L`
- 有的 basis 函数偏更高动量

所以当系数 `u_j` 改变时，你看到的不是“抽象数字变了”，而是：

> 某些带特定动量标签的 basis direction 被更强地占据了

这会让人感觉：

> “好像 basis 自己在动”

但严格来说不是 basis 在动，而是：

> **固定 basis 上的占据权重在重新分配**

## 一个最直接的类比

把当前固定 basis 想成一排固定好的乐器：

- 钢琴
- 小提琴
- 长笛
- 大提琴

这些乐器没有变。

变的是当前这首曲子里：

- 哪个乐器声音更大
- 哪个乐器几乎不响

所以音乐在变，但乐器没变。

当前 ECG 固定基方法里：

- basis functions = 那排固定乐器
- coefficients `u_j` = 每个乐器此刻的音量

kick 后，变的是“这首曲子怎么分配到各个乐器上”，不是乐器本身换了。

## 在数学上到底是什么没变，什么变了

### 没变的

固定基方法里，始终是同一个有限维子空间：

`V = span{phi_1, ..., phi_K}`

这里的 `phi_j` 是一开始准备好的增广基函数。

这个 `V` 在整个 kicked evolution 里不变。

### 变的

变的是状态向量在这个子空间里的坐标：

`Psi(t) = sum_j u_j(t) phi_j`

所以随着时间变化：

- `Psi(t)` 在变
- `u(t)` 在变

但：

- `phi_j` 不变
- `span{phi_j}` 不变

## 那“动量变了”体现在哪里

体现在：

1. 系数 `u_j(t)` 重新分配到不同动量标签的 basis function 上
2. 不同 basis function 之间的干涉项也在变化

所以即使 basis 固定，波函数的动量分布也可以显著变化。

这就像 Fourier 展开里：

`f(x,t) = sum_n c_n(t) e^{i n x}`

这里平面波 `e^{i n x}` 没变，但只要 `c_n(t)` 变了，整个函数的谱分布就变了。

当前 ECG 的情况更复杂一些，因为 basis 不是正交平面波，而是带 Gaussian 包络的非正交函数，但逻辑是一样的。

## 什么时候才叫“基底在变”

只有当下面这些东西真的被更新时，才能说 basis 在变：

- 新增基函数
- 删除基函数
- 修改某个基函数的 `A`
- 修改某个基函数的 `B`
- 修改某个基函数的 `R`

这类方法才是：

- dynamic basis expansion
- adaptive basis update
- time-dependent basis

### 当前代码不是这种方法

当前代码在 kicked evolution 阶段：

- 不新增 basis
- 不修改 `A/B/R`
- 只更新 `u`

所以它仍然是一个**固定基方法**。

## 那你刚才的直觉哪里特别有价值

你真正抓住的重点是：

> 虽然当前代码是固定基方法，但系统真正需要的“有效动量支撑”在随 kick 次数持续变化

这正是为什么固定 basis 迟早会吃力。

也就是说，你的直觉虽然在术语上需要修正：

- 不是“basis 已经变了”
- 而是“basis 应该跟着态变，但当前代码没有变”

这其实正好指出了当前方法的深层局限。

## 一句话总结

**当前代码里，kick 之后变的是波函数和系数，不是 basis 本身；但你说对了一件更重要的事：真实物理态需要的动量支撑是在持续变化的，而当前固定 basis 并没有跟上这一点。**

---

# Q&A: 一次失败的尝试：直接引入 first-kick image candidates（2026-04-01）

已经实际尝试过一版最直接的实现：

- 对每个原始基函数 `phi_j`
- 直接生成它在离散动量标签下的 image：
  - 单粒子：`exp(i 2 n k_L x) * phi_j`
  - 两粒子：`exp(i 2 n1 k_L x1) exp(i 2 n2 k_L x2) * phi_j`

从物理上，这个方向很自然，因为第一次 kick 后的真实态本来就是这类对象的叠加。

## 但在“当前贪心筛选 + 非正交 raw basis”框架下，结果更差

### 情况 A：`max_cond = 1e10`

```text
Candidate pool: 440
Accepted momentum basis: 51
K_total = 61
First-kick projection fidelity = 12.82
E after kick: ECG = 66.23, exact = 1.315
```

这一版已经完全数值失控。

### 情况 B：`max_cond = 1e6`

```text
Candidate pool: 440
Accepted momentum basis: 17
K_total = 27
First-kick projection fidelity = 1.00322
E after kick:    ECG = 1.44427, exact = 1.31533, dE  = 0.12894
<x^2> after kick: ECG = 1.03279, exact = 1,        dx2 = 0.03279
<p^2> after kick: ECG = 1.85576, exact = 1.63067,  dp2 = 0.22509
```

这比当前原始启发式候选的基线结果还差：

```text
原始启发式候选（默认设置）：
E after kick:    dE  = 0.09349
<x^2> after kick: dx2 = 0.02788
<p^2> after kick: dp2 = 0.15911
```

## 这说明什么

说明至少在“当前候选筛选框架”下：

- **物理上更对的候选**
- 不等于
- **数值上更好用的候选**

当前这版失败，主要不是因为物理想法完全错，而是因为：

1. first-kick image candidates 之间高度相关
2. 进入当前非正交贪心筛选后，`S` 更容易接近病态
3. `u_after = S^{-1} K u_before` 的数值放大比原来更严重

## 当前阶段应得出的结论

这一轮实验的结论不是：

> “first-kick image basis 这个方向错了”

而是：

> “first-kick image basis 不能直接塞进当前这套 raw nonorthogonal greedy pipeline；如果要走这条路，几乎肯定要和正交化 / 更强筛选一起做”

这反而进一步支持了 Step B 的必要性。

---

# Q&A: `x^2` 和 `p^2` 的误差到底从哪来，有没有数学推导（2026-04-01）

有，而且这个推导能把问题说得很清楚。

## 1. 先定义 3 个态

设初态是：

`psi_0`

第一次 kick 后的 **exact** 态是：

`psi_ex = U_kick psi_0`

而代码真正保留下来的，是投影回当前有限子空间后的近似态：

`psi_pr = P_V (U_kick psi_0)`

如果再考虑代码里后面的重归一化，那么更严格地说应写成归一化后的 `psi_pr`，这里记号就不再细分。

定义误差态：

`delta psi = psi_pr - psi_ex`

那么所有物理量的误差，本质上都来自这个 `delta psi`。

## 2. 任意可观测量的误差公式

对任何算符 `A`，投影后和 exact 的期望值差是：

`Delta<A> = <psi_pr|A|psi_pr> - <psi_ex|A|psi_ex>`

把 `psi_pr = psi_ex + delta psi` 代进去，得到：

`Delta<A> = <delta psi|A|psi_ex> + <psi_ex|A|delta psi> + <delta psi|A|delta psi>`

也就是：

`Delta<A> = 2 Re <delta psi|A|psi_ex> + <delta psi|A|delta psi>`

这条公式很重要。

它说明：

- 如果 `delta psi = 0`，所有物理量都完全正确
- 只要投影造成 `delta psi != 0`，任何物理量都可能出错

所以从最一般的层面说：

> **是的，误差的根源就是投影后的态和 exact kicked state 不一样**

但不同算符对这个误差的敏感程度不同。

## 3. 为什么 `x^2` 的误差特别干净、特别能说明问题

对 exact kick：

`U_kick = exp[-i kappa cos(2 k_L x)]`

它只依赖 `x`，所以与任何位置函数都对易：

`[U_kick, f(x)] = 0`

特别地：

`U_kick^dagger x^2 U_kick = x^2`

因此 exact 情况下：

`<x^2>_ex = <psi_0|x^2|psi_0>`

也就是说：

> **第一次 kick 在 exact 情况下不会改变 `<x^2>`**

所以如果我们现在看到：

- exact: `<x^2> = 1`
- ECG projection 后: `<x^2> = 1.02788`

那这个误差就没有任何“物理本来就会变”的借口。

它只能来自：

`delta psi != 0`

也就是：

> **投影后的态已经把本来不该变的空间概率分布也扭曲了**

这也是为什么我说 `<x^2>` 的误差是一个非常强的诊断指标。

### 用密度的角度看更直观

exact kick 后：

`psi_ex(x) = e^{-i phi(x)} psi_0(x)`

所以：

`|psi_ex(x)|^2 = |psi_0(x)|^2`

位置密度完全不变。

但投影后我们保留的是 `psi_pr`，一般不再能写成：

`psi_pr(x) = e^{-i phi(x)} psi_0(x)`

于是：

`|psi_pr(x)|^2 != |psi_0(x)|^2`

这就直接导致：

- `<x^2>`
- 乃至所有只依赖位置密度的观测量

都开始出错。

## 4. 为什么 `p^2` 的误差通常更大

因为 `p^2` 本来就会在 kick 后改变，而且它对相位梯度非常敏感。

### 4.1 exact 情况下 `p` 怎么变

由

`U_kick = exp[-i kappa cos(2 k_L x)]`

可以算出：

`U_kick^dagger p U_kick = p + 2 kappa k_L sin(2 k_L x)`

于是：

`U_kick^dagger p^2 U_kick = (p + 2 kappa k_L sin(2 k_L x))^2`

展开为：

`p^2 + 2 kappa k_L {p, sin(2 k_L x)} + 4 kappa^2 k_L^2 sin^2(2 k_L x)`

所以 exact 情况下，`p^2` 会真实地发生变化。

### 4.2 为什么投影误差会把 `p^2` 搞得更厉害

因为 `p^2` 不只是看振幅，还看导数：

`p^2 ~ - d^2/dx^2`

而 kick 的本质正是在波函数里注入位置依赖相位：

`psi_ex(x) = e^{-i phi(x)} psi_0(x)`

这个相位的导数直接进入动量分布。

因此如果投影没有把这个相位结构表示好，那么：

- 振幅会错
- 相位会错
- 相位梯度会更错

于是 `p^2` 的误差通常会比 `x^2` 更敏感、更大。

这和当前数据完全一致：

```text
dx2  ~ 2.79e-2
dp2  ~ 1.59e-1
```

## 5. 这是不是说明“只是投影效果太差了”

严格说，应该分成两层：

### 第一层：物理截断误差

也就是：

- exact kicked state 不在当前有限子空间里
- 所以即使做“完美投影”，仍然会有 `delta psi != 0`

这是最本质的误差来源。

### 第二层：数值放大误差

当前我们又是在非正交 basis 里解：

`S u_after = K u_before`

当 `S` 很病态时，`S^{-1}` 会把某些方向上的误差放大。

所以现在看到的误差是：

`总误差 = 子空间不完备带来的误差 + 病态数值放大`

所以不能只说一句“投影效果太差了”就结束，应该更准确地说：

> **根本原因是当前 basis 子空间承接不了 exact kicked state，而病态 `S^{-1}` 又把这个截断误差进一步放大。**

## 6. 那应该优化什么

既然误差来源已经拆清楚了，优化目标也就更清楚了。

### 对 `x^2` 来说

因为 exact 第一脚后 `<x^2>` 应该严格不变，所以可以把

`| <x^2>_pr - <x^2>_ex |`

直接当成一个硬指标。

这几乎就是在要求：

> 投影后不要把位置密度扭曲掉

### 对 `p^2` 来说

则是在要求：

> 投影后不仅振幅要对，局部相位梯度也要对

这比 `x^2` 更难。

### 一个更完整的 first-kick 目标函数

所以第一次 kick 后如果要优化 basis，一个合理目标可以是：

`J = w1 * (1 - F_state) + w2 * |dx2| + w3 * |dp2| + w4 * penalty(cond(S))`

其中：

- `F_state` 是 first-kick 态的 overlap / fidelity
- `dx2` 约束位置分布别被扭曲
- `dp2` 约束动量结构别偏得太厉害
- `penalty(cond(S))` 防止数值病态

## 一句话总结

**`x^2` 和 `p^2` 的误差都来自投影后的误差态 `delta psi`；其中 `x^2` 的误差尤其说明投影已经把本来不该变的空间概率分布扭曲了，而 `p^2` 因为还对相位梯度敏感，所以通常会错得更厉害。真正要优化的不是 `U_kick` 的幺正性，而是让 basis 更好地承接 `U_kick psi`。**

---

# Q&A: 混合选基思路行不行？能不能在原有基底方法上再叠加新的选基方式（2026-04-01）

这个思路是**可行的，而且我觉得比“完全替换原方法”更稳妥**。

也就是说，可以考虑：

> **保留当前原始启发式动量基作为主干**
> +
> **再额外加入一层更贴近 first-kick 的候选基**

然后不要求这些新候选一定大量占据，而是让广义本征解和后续投影自己决定它们要不要被用上。

## 1. 这个思路为什么有吸引力

它的优点在于：

### 优点 1：保留当前方法的稳定主干

原始启发式动量基虽然不够好，但至少：

- 数值上还能跑
- 第一脚不会立刻炸
- 基态一致性检查也能过

所以它可以继续作为“稳定底座”。

### 优点 2：把新的 first-kick information 当成“修补层”

新加入的 first-kick-oriented 候选，不是用来完全取代原始 basis，而是用来补那些原有 basis 最接不住的方向。

这就很像：

- 老 basis = 主骨架
- 新候选 = repair layer

### 优点 3：如果某些新候选并不重要，理论上系数确实可以自动变小

这一点你的直觉是对的。

如果新增 basis 真的不重要，那么在最后的最优表示里：

- 它们的系数可以很小
- 甚至接近 0

所以在“物理表示”层面，混合候选库是说得通的。

## 2. 但有一个关键风险：不能指望“小系数”自动解决数值问题

这一点一定要说清楚。

在当前非正交 basis 里：

> **即使某个新基函数最后的系数很小，它仍然可能先把 overlap 矩阵 `S` 搞得很病态。**

而一旦 `S` 病态：

- `u_after = S^{-1} K u_before`
- 这里的 `S^{-1}` 就会放大误差

所以：

### 物理层面

“不重要的 basis 系数会自动变小”

这个说法是合理的。

### 数值层面

“既然系数会变小，那就不会有害”

这个说法**不成立**。

这就是为什么我们不能简单地说：

> “把所有新候选都丢进去，反正系统自己会把它们系数压小”

在非正交广义本征值问题里，这样很容易先把 `S` 弄坏。

## 3. 所以最合理的不是“全塞进去”，而是“受控混合”

我觉得可行方案有三档。

## 方案 A：保守型混合候选库（我最推荐）

### 结构

候选库分成两层：

1. **主干层**
   - 当前原始启发式动量 Gaussian

2. **修补层**
   - 少量 first-kick-oriented 候选

### 修补层怎么选

不要把全部 first-kick image 都扔进去，而只取：

- 最重要的少数原始基函数
  - 例如按 `|u_j|` 排名前 `M` 个
- 最显著的低阶 kick 谐波
  - 例如 `n = ±1, ±2`

也就是说，只生成：

`exp(i 2 n k_L x) * phi_j`

其中：

- `j` 只取重要的几个
- `n` 只取最显著的几个

### 优点

- 风险最小
- 先看看第一脚能不能改善
- 不会像之前那次实验一样一下子把候选做得太大

### 缺点

- 改善幅度可能有限

## 方案 B：first-kick 目标函数驱动的贪心选基

这个更系统一点。

### 核心想法

对每个候选基 `chi`，不是只看 `cond(S)`，而是看加入它后第一脚误差有没有下降。

### 目标函数

我认为可以直接用一个 first-kick objective：

`J = w1 * (1 - F_state) + w2 * |dx2| + w3 * |dp2| + w4 * |dE| + w5 * penalty(cond(S)) + w6 * gs_penalty`

其中：

- `F_state`：第一脚后态和 exact 态的 overlap（如果暂时没有，就先用弱一点的替代量）
- `dx2`：第一脚后 `<x^2>` 误差
- `dp2`：第一脚后 `<p^2>` 误差
- `dE`：第一脚后能量误差
- `penalty(cond(S))`：惩罚病态
- `gs_penalty`：惩罚增广后把初始基态搞坏

### 具体做法

对每个候选：

1. 先加入当前 basis
2. 重新解增广后基态
3. 做第一脚 kick
4. 计算 `J`
5. 只接受那些让 `J` 真正下降的候选

### 优点

- 完全围绕你真正关心的误差在选基
- 比单纯几何上的 `cond(S)` 筛选更有物理针对性

### 缺点

- 计算量会明显变大
- 需要先把第一脚的 stronger diagnostics 都稳定接好

## 方案 C：混合候选 + 正交化修补层

这其实是我觉得从长远看最靠谱的。

### 结构

1. 保留当前原始启发式 basis 作为主干
2. 额外生成一小批 first-kick image repair candidates
3. **只对 repair layer 做局部正交化 / SVD 截断**
4. 再把这批“已经去掉近重复方向”的修补基并入主干

### 这比直接塞 raw image candidates 好在哪里

因为上一次实验已经说明：

- raw first-kick image candidates 物理上有道理
- 但数值上高度相关

所以真正危险的不是“新候选方向本身”，而是：

> 它们在非正交表示里太接近，直接把 `S` 弄坏

局部正交化正好对付这个问题。

### 优点

- 保留你这个“混合方案”的核心想法
- 又避免 raw image candidates 直接炸掉

### 缺点

- 实现复杂度高于方案 A
- 但低于完全重做整个 Step B

## 4. 我觉得最现实的落地顺序

如果要按风险和收益平衡，我建议：

### 第一步：先做方案 A

也就是一个**小规模受控混合库**：

- 原始启发式 basis 保留
- 只加少量最重要的 first-kick repair candidates

这一步主要看：

- `dx2` 能不能降
- `dp2` 能不能降
- `dE` 能不能降
- `cond(S)` 会不会立刻恶化

### 第二步：如果 A 有改善，再升级成方案 B

把“是否接受新 basis”的规则，从现在的纯 `cond(S)` 筛选，升级成：

- `first-kick objective + cond penalty`

### 第三步：如果 raw 候选还是太病态，再上方案 C

也就是给 repair layer 做局部正交化。

## 5. 所以我对你这个思路的判断

### 我认为这个方向是对的

因为它避免了两个极端：

- 极端 1：完全信任原始 basis，不补任何新的 first-kick 信息
- 极端 2：完全替换成一套新的 image basis，结果数值上直接炸掉

你的想法正好是中间路线：

> 让原始 basis 继续当底座，再用少量更物理的候选去修补第一脚误差

### 但实现上必须带“约束”

最重要的不是“把新候选加进去”，而是：

- 控制数量
- 控制相关性
- 用 first-kick objective 来筛
- 同时保住基态一致性和 `cond(S)`

## 一句话总结

**这个混合思路是有前途的，但不能简单理解成“多加点 basis，反正不重要的系数会自动变小”；在非正交广义本征值问题里，必须把数值病态一起纳入目标函数。最现实的做法，是从“小规模受控混合候选库 + first-kick 目标函数筛选”开始。**

## 已做的一次实验：小规模受控混合候选库 + first-kick 目标函数筛选

已经按这个思路做了一个最小原型：

### 原型设置

1. **保留原始启发式增广基** 作为稳定底座
2. 只生成一小批 repair candidates：
   - 取最重要的 `3` 个原始静态基函数
   - 只取少量低阶动量标签
   - 总 repair candidate 数量 = `18`
3. 对每个 repair candidate：
   - 加到当前 basis 中
   - 重新求增广后的基态
   - 计算 first-kick 诊断
   - 用下面这个目标函数判断是否接受：

`J = rel(dE) + rel(dx2) + rel(dp2) + 0.25*|1-fidelity| + 5*(1-gs_fid) + 0.05*(cond/max_cond)`

这里：

- `dE`：第一脚瞬时能量误差
- `dx2`：第一脚瞬时 `<x^2>` 误差
- `dp2`：第一脚瞬时 `<p^2>` 误差
- `fidelity`：第一脚投影 fidelity
- `gs_fid`：增广后基态与原基态的一致性
- `cond/max_cond`：病态惩罚

### 实际结果

运行结果是：

```text
Baseline first-kick objective = 0.24671917

--- First-Kick Repair Search ---
candidate count = 18
no improving repair candidate found after 0 accepted candidates

Post-repair basis: K = 24
first-kick objective = 0.24671917
```

也就是说：

- 在当前这套“小规模 repair candidate + first-kick 目标函数”原型下
- **没有任何一个新增候选能同时在 first-kick 物理误差和数值稳定性之间带来净改善**

### 这说明什么

说明当前问题不是“随便补几个更物理的方向就行”。

更准确地说：

1. **当前启发式基底已经把可接受的低成本方向差不多用掉了**
2. 新增的少量 repair candidate，要么改善太小，要么带来的 `cond(S)` / ground-state 扰动抵消了收益
3. 因此只靠“小修小补”很可能不够

### 这轮实验给出的判断

这轮实验并没有证明“混合候选库方向错了”，但它说明：

> **如果不同时改进 basis 的数值组织方式，仅靠少量 raw repair candidates，很难把 first-kick 误差再明显压下去。**

这进一步支持两个判断：

1. first-kick 误差不是一个特别局部的小问题
2. 真正的突破口更可能在：
   - 正交化子空间
   - 或更系统的动态 basis update

---

# Q&A: 这次 repair 原型里的候选基底是不是不够多（2026-04-01）

## 短答案

**是，从“表达能力”角度看，这次 repair 原型的候选库大概率是偏少的。**

但同时：

**不能简单地靠“多塞一些 raw candidate”来解决，因为当前主要瓶颈还包括病态 `S`。**

所以更准确的话应该是：

> repair candidate 的数量和多样性确实不够，但问题也不只是“数量不够”，而是“数量一上去以后，当前非正交贪心框架会先变得数值不稳定”。

## 1. 为什么说这次 repair candidate 确实偏少

这次原型里，我们是故意做得非常保守的：

- 只取最重要的 `3` 个原始静态基函数
- 只取 `6` 个低阶动量标签
- 所以总 repair candidate 只有：

`3 * 6 = 18`

这和原始启发式候选库比起来其实很小：

- 原始启发式候选库：`528`
- repair 原型候选库：`18`

所以从纯组合规模上讲，18 个候选当然是非常保守的。

## 2. 为什么它不只是“数量少”，还包括“种类少”

更关键的是，这 18 个 repair candidate 只在两件事上做了变化：

1. 选了少量原始 basis function
2. 给它们乘了少量低阶动量标签

但它**没有**去探索下面这些自由度：

- 不同 width
- 不同 `A` 的相关结构
- 不同 `B` 的分配方式
- 更高阶动量标签
- 不同粒子之间更丰富的动量组合

所以这次 repair 原型不仅“数量少”，而且“形状种类也少”。

## 3. 但为什么我没有一开始就把 repair candidate 做很多

因为前面的实验已经说明：

### 情况 A：raw image candidate 做很多

结果是：

- 候选更物理
- 但 basis 之间高度相关
- `S` 迅速接近病态
- 第一脚直接数值爆炸

### 情况 B：少量 repair candidate

结果是：

- 数值还稳
- 但改善不明显

所以我们现在卡在一个典型张力里：

- **候选太少**：表达能力不够
- **候选太多**：当前 nonorthogonal greedy pipeline 顶不住

这就是为什么我一直说，真正的问题已经不只是“多不多”，而是：

> **在当前框架里，candidate richness 和 numerical stability 是互相打架的**

## 4. 所以这个问题该怎么更准确地理解

最准确的说法不是：

> “这次没改善，因为候选库太少”

而是：

> “这次候选库既偏少，又被限制成了非常保守的低阶修补层；它不足以跨过当前误差瓶颈，但如果在当前框架里粗暴加大，又很可能先被病态 `S` 拖垮。”

## 5. 下一步如果还想沿这个方向推进，应该怎么加大

如果还想继续探索混合 repair 路线，我认为应该按下面顺序扩展，而不是直接无脑增大数量：

### Step 1：先增加 source basis 数量

例如：

- 从前 `3` 个重要原始 basis
- 提高到前 `5` 个或前 `6` 个

这样增加的是“原始基态几何信息”的覆盖度。

### Step 2：再增加动量标签的层数

例如：

- 现在只看很低阶标签
- 下一步扩到更多 `(m1,m2)` 组合

这样增加的是“动量方向”的覆盖度。

### Step 3：最后才增加包络形状自由度

也就是再考虑：

- width 变化
- `A/B` 结构变化

这是最容易把 `S` 搞坏的一层，应该最后动。

## 6. 但这条路的前提是什么

前提是：

> 每次扩展都要同时盯住 `cond(S)` 和 first-kick objective

不然你会很容易看到一个假象：

- 候选更多了
- 理论上更完备了
- 实际结果却更差

这并不一定是物理候选错了，而是数值组织方式先崩了。

## 一句话总结

**这次 repair 原型的候选库确实偏少，但当前问题不能靠简单“堆更多 raw candidate”来解决；如果要继续走这条路，必须按“先扩 source basis，再扩动量标签，最后扩包络形状”的顺序，和 `cond(S)` 约束一起做。**

## 已做的一次扩展实验：把 repair candidate 放大一档

按照上面的思路，已经把 repair 池从非常保守的版本扩大到一档中等规模：

- source basis: `3 -> 6`
- 动量标签数: `6 -> 12`

所以 repair candidate 总数从：

`18 -> 72`

## 实际结果

运行后，repair search 的结果是：

```text
--- First-Kick Repair Search ---
candidate count = 72
accept base=0, m=(0,3)   -> objective = 0.23603745
accept base=3, m=(-2,1)  -> objective = 0.23395217
no improving repair candidate found after 2 accepted candidates

Post-repair basis: K = 26
```

也就是说：

- 这次不再是“一个都选不进来”
- 而是确实找到并接受了 `2` 个 repair candidates

## 第一脚瞬时诊断的变化

### 原基线（不加 repair）

```text
Projection fidelity = 0.9989555239
dE   = 9.3494573850e-02
dx2  = 2.7879977035e-02
dp2  = 1.5910917067e-01
```

### 扩大 repair 池后

```text
Projection fidelity = 0.9991319039
dE   = 8.6927463067e-02
dx2  = 2.8591686588e-02
dp2  = 1.4526323955e-01
```

## 这说明什么

### 改善的地方

- `projection fidelity` 略好
- `dE` 变小
- `dp2` 变小

所以这次扩大 repair 池后，**第一脚的能量和动量误差确实有小幅改善**。

### 没改善、甚至略差的地方

- `dx2` 没有改善，反而略微变大

这说明：

> 当前这组 repair candidates 更像是在修补 kick 后的相位 / 动量结构，
> 但对“保持位置密度不变”这一约束，帮助还不够。

## 当前阶段的判断

这次结果很重要，因为它介于前两次实验之间：

1. **太小的 repair 池**
   - 一个都选不进来
   - 几乎没效果

2. **太激进的 raw image 替换**
   - 直接数值爆炸

3. **中等规模、受控 repair 池**
   - 已经能带来一些改善
   - 但改善还不够，尤其是 `dx2`

所以现在可以更有把握地说：

> **“混合基底 + first-kick objective” 这条路不是完全没用，确实能带来改善；但当前 repair 候选的方向还不够对，尤其还没有真正抓住 `<x^2>` 这个位置密度约束。**

## 下一步最自然的改进方向

如果继续走这条线，最自然的下一步不是盲目继续扩数量，而是：

1. **在目标函数里提高 `dx2` 的权重**
   - 因为 exact 第一脚后 `<x^2>` 本来就必须守住

2. **专门设计更“保位置密度”的 repair candidates**
   - 当前这轮接受的两个 candidate 更偏动量修补
   - 还缺少对位置密度不变这一约束更敏感的方向

一句话总结：

**把 repair 池从 18 扩到 72 以后，第一脚误差已经出现了“小幅但真实”的改善；这说明混合 repair 路线是有信号的，但下一步应该重点围绕 `dx2` 去继续优化，而不是只一味加大 candidate 数量。**

## 已做的又一轮实验：repair 池再放大一档，并把 `dx2` 权重加大

接着又试了一轮更激进但仍受控的版本：

- source basis: `6 -> 8`
- 动量标签上限：进一步放大
- repair candidate 总数：`72 -> 160`
- 同时把 first-kick objective 里 `dx2` 的权重显著提高

## 实际结果

repair search 的结果是：

```text
--- First-Kick Repair Search ---
candidate count = 160
accept base=0, m=(0,3)
accept base=3, m=(-2,1)
no improving repair candidate found after 2 accepted candidates
```

也就是说：

- repair 候选库确实变大了
- `dx2` 权重也变大了
- 但最终被选中的仍然是和上一轮同一组 repair candidate

### 第一脚瞬时诊断

结果也几乎没有进一步变化：

```text
Projection fidelity = 0.9991319039
dE   = 8.6927463067e-02
dx2  = 2.8591686588e-02
dp2  = 1.4526323955e-01
```

和上一轮 `72` 个 repair candidate 的结果基本一样。

## 这说明什么

这说明当前 repair 路线已经出现了一个很清晰的信号：

> **问题不是简单“候选不够多”**

更准确地说：

1. repair 池从 `18 -> 72` 时，确实出现了改善
2. 但从 `72 -> 160`，并且还加大了 `dx2` 权重之后，结果几乎不再变化
3. 这说明当前 first-kick objective 在现有 raw repair candidate 家族里，已经很快选到了它认为最有价值的那几个方向

也就是说，当前瓶颈更像是：

- **候选家族本身的表达能力有限**
- 而不是“只要继续堆更多同类 candidate 就会持续改善”

## 当前应得出的判断

到这一步，可以更明确地说：

> **在当前 nonorthogonal greedy + raw repair candidate 的框架下，继续单纯扩大 repair 候选库，收益已经明显趋于饱和。**

所以后续要继续突破，优先级更高的不是再继续堆量，而是：

1. 改 repair candidate 的“类型”
2. 或改这些 candidate 的数值组织方式（例如正交化 / 局部子空间化）

一句话总结：

**候选库从 72 扩到 160、同时提高 `dx2` 权重之后，结果几乎没有继续改善，这说明当前问题已经不主要是“数量不够”，而更像是“同一类 raw candidate 的可用信息已经被吃得差不多了”。**

---

# Q&A: 当前主线重新聚焦到什么问题（2026-04-01）

需要明确纠正一件事：

> 当前 kicked 主线里，真正卡住的问题不是“TDVP 导数是不是对”，而是“怎么把第一脚的 `dx2` 压下来，以及怎么在多次 kick 中动态更新 basis 同时保持高投影精度”。

原因是：

- 现在 real-time kicked 主线用的是
  - `apply_analytic_kick()`
  - `free_evolve_fixed_basis()`
- 并没有在实时间演化里继续依赖 TDVP 导数去推进 kick

所以后续讨论应聚焦在两件事：

1. **怎么降低 first-kick 的 `dx2`**
2. **怎么在 repeated kicks 下更新 basis，同时维持高投影精度**

---

# Q&A: “先对 repair layer 做局部正交化”到底是在干什么（2026-04-01）

## 短答案

它的意思不是“换一套全新的基底”，而是：

> **先保留当前已经稳定的主干 basis 不动，只把新增的 repair candidates 单独拿出来，去掉其中彼此几乎重复、几乎平行的方向，再把剩下那些真正独立的修补方向接回主干。**

所以它解决的不是“候选方向对不对”，而是：

> **这些候选方向虽然有用，但 raw nonorthogonal 形式太互相重叠了，直接拿来算会把 `S` 搞病态。**

## 1. 为什么要这么做

当前 repair candidates 的问题不是完全没有信息，而是：

- 很多 candidate 彼此非常像
- 它们在非正交表示里高度相关
- 一旦直接并进主干 basis，`S` 就更容易病态

结果就是：

- 物理上也许加了“对的方向”
- 数值上却先被 `S^{-1}` 放大误差搞坏

所以局部正交化的目的，是把 repair 候选里的“重复信息”先压缩掉。

## 2. 为什么叫“局部”

因为我们不是把**整个** basis 全部推倒重来做正交化，而只是对：

- 原有稳定主干 basis 之外的
- 那一小层 repair candidates

单独处理。

也就是说，分层看是：

1. **主干层**
   - 当前启发式 basis
   - 先不动

2. **repair 层**
   - 新加的修补候选
   - 只对这一层做正交化/压缩

所以它叫“局部正交化”。

## 3. 直觉版理解

你可以把 raw repair candidates 想成一堆很像的箭头：

- 方向大致都指向“第一脚 kick 后需要补的区域”
- 但很多箭头彼此几乎重合

如果你把 20 支几乎重合的箭头全都拿进来：

- 你并没有真的多出 20 个独立信息方向
- 你只是制造了一个非常病态的坐标系

局部正交化做的事情就是：

> 把这 20 支几乎重合的箭头，压缩成 2 到 5 个真正独立的方向

这样：

- 信息尽量保留
- 冗余尽量去掉
- 数值稳定性好很多

## 4. 数学上在做什么

设当前主干 basis 是：

`B_main = {phi_1, ..., phi_K}`

repair candidates 是：

`R_raw = {chi_1, ..., chi_M}`

我们先构造这些 block 矩阵：

- `S_mm`：主干和主干的 overlap
- `S_mr`：主干和 repair 的 overlap
- `S_rr`：repair 和 repair 的 overlap

如果只看 repair 内部，最简单的局部正交化就是对 `S_rr` 做本征分解或 SVD：

`S_rr = U Lambda U^dagger`

然后丢掉很小的特征值方向，只保留：

`lambda_i > cutoff`

再把保留下来的方向重新组合成新的 repair basis：

`eta_alpha = sum_a U_{a alpha} / sqrt(lambda_alpha) * chi_a`

这样新 repair basis 在 repair 子空间内部就近似正交了。

## 5. 更严格一点：还要去掉与主干重复的部分

只正交化 `S_rr` 还不够，因为 repair candidates 还可能和主干 basis 高度重叠。

所以更合理的是先把 repair 层投影到“主干的补空间”里，再做正交化。

形式上可以写成：

`chi_tilde = (I - P_main) chi`

这里 `P_main` 是主干 basis 张成子空间上的投影。

在矩阵语言里，这对应于一个 Schur complement / 有效 overlap：

`S_eff = S_rr - S_rm * S_mm^{-1} * S_mr`

然后对这个 `S_eff` 再做本征分解或 SVD。

这一步的含义是：

> **先去掉 repair candidates 里那些主干已经会表示的部分，再对剩下真正“新的信息”做压缩。**

## 6. 做完之后得到了什么

做完以后，不再直接把 raw candidates `{chi_a}` 并进来，而是得到一小组新的 repair directions：

`R_orth = {eta_1, ..., eta_r}`

其中：

- `r` 通常远小于原始候选数 `M`
- `eta_alpha` 是 raw repair candidates 的线性组合
- 这些方向彼此更独立
- 也更不容易把 `S` 弄病态

然后最终工作空间变成：

`B_total = B_main + R_orth`

## 7. 它的好处是什么

### 好处 1：保留候选信息

不是简单粗暴删掉很多 candidate，而是尽量保留它们张成的子空间。

### 好处 2：去掉近重复方向

这正是当前 raw repair candidates 最大的问题。

### 好处 3：更适合后续做目标函数筛选

因为你筛选的是“真正独立的修补方向”，而不是一堆互相高度重复的 raw candidate。

### 好处 4：比“整个 basis 全局正交化”更保守

它不推翻当前主干方法，只是给 repair 层做数值整理。

## 8. 它的代价是什么

代价也很明确：

### 代价 1：新的 repair direction 不再是简单单个 Gaussian

它会变成：

`eta_alpha = linear combination of raw chi_a`

所以它的解释性会差一点。

### 代价 2：实现比现在复杂

因为你不能再只把它当作“多了几个 `BasisParams`”来处理，而要开始显式地处理：

- block overlap
- block Hamiltonian
- 子空间变换矩阵

### 代价 3：它更像是“子空间算法”，而不是“再加几个 basis function”

这意味着代码组织上会更接近 Step B。

## 9. 一句话理解

所以“先对 repair layer 做局部正交化”的真正意思是：

> **不要把一堆彼此高度相关的 raw repair candidates 直接塞进当前 nonorthogonal basis；先把这批修补候选压缩成少数几个真正独立的方向，再看看这些独立修补方向能不能把 first-kick 误差往下压。**

## 一句话总结

**局部正交化 repair layer，不是在“换物理模型”，而是在“先把 repair 候选里的重复信息清掉，再拿真正独立的修补方向去帮主干 basis 接第一脚 kick”。**

---

# 实现评估（修订版）：First-Kick 目标函数驱动的贪心选基（2026-04-01，Claude 评估）

## 0. 修订说明

初版评估只看了 qa.md 末尾的目标函数提案和代码库结构，没有读完前面两轮实验和 first-kick image candidates 失败尝试的完整上下文。现在基于全文重新评估。

## 1. 核心判断：first-kick 目标函数 **能做、但不能解决根本问题**

先说结论：

- **技术上完全可以实现**（代码积木齐全，计算量不是瓶颈）
- **但它解决的是”在 ~14 个可接受候选里选得更好”，而不是”让可接受候选变多”**
- **当前真正的结构性瓶颈，不在选基准则，而在非正交 Gaussian 基之间 overlap 太高导致只有 ~14 个候选能通过任何合理的 cond(S) 门槛**

两轮实验已经非常清楚地表明了这一点：

| 实验 | 候选池 | 最终接受 | 结论 |
|------|--------|----------|------|
| 第一次（n_mom=4, max_cond=1e6） | 660 | 12 | 贪心 + cond 硬截断，有效基太少 |
| 第二次（max_cond=1e10, 减少 widths） | 528 | 14 | 放宽到 1e10 也只多了 2 个；但 S^{-1} 放大让长期更差 |
| first-kick image candidates（max_cond=1e6）| 440 | 17 | 物理上更对的候选反而数值上更差 |

这三组数据说明：

> **问题不在”选了哪 14 个”，而在”只能选 14 个”。**

用 first-kick objective 替代纯 cond(S) 作为选基准则，即使能选出”更优的 14 个”，也不太可能从根本上改变 kicked dynamics 的精度——因为 14 个动量基对 44 个动量格点来说差距太大。

## 2. 那它还有什么价值？

有，但比初版评估里说的要弱。具体来说：

### 2.1 early-time（前几步）可能有改善

如果在那 ~14 个名额里，选出对**第一脚**最有价值的候选（而不是碰巧先遍历到的），那前几步的 dE、dx2、dp2 有可能改善一些。

但根据 first-kick image candidates 失败的经验：物理上更对的候选，在当前非正交框架下反而更容易搞坏 S。所以 first-kick objective 的改善幅度很可能有限。

### 2.2 它可以作为 Step B（正交化子空间）的配套工具

如果先做 Step B（对整个候选池做全局正交化 / SVD 截断），然后在正交化后的子空间里再用 first-kick objective 来做 basis truncation，那它的价值就大很多。因为：

- 正交化解决了”只能接受 14 个”的结构性瓶颈
- first-kick objective 解决了”在更大可选集里怎么取舍”的质量问题

**这两件事配合才有意义；单独做 first-kick objective 意义有限。**

## 3. 为什么 first-kick image candidates 会失败——这直接约束了 first-kick objective 的上限

qa.md 里已经记录了这次失败尝试。核心教训是：

first-kick image candidates = `exp(i 2 n k_L x) * phi_j`，是物理上最自然、最正确的候选。但它们之间的 overlap 极高——因为同一个 phi_j 乘上不同的 n 之后，Gaussian 包络几乎一样，只是相位不同。在非正交表示里，这些相位差异表现为 S 矩阵接近奇异。

first-kick objective 不能改变这一点。**候选之间的 overlap 结构是由候选本身决定的，不是由选基准则决定的。**

这是最关键的约束：即使你有一个完美的目标函数 J，它也不能让那些因为 S 病态而不能共存的候选突然变得可以共存。

## 4. 如果还是想做，实现上的调整建议

考虑到上述限制，如果仍然要做 first-kick objective（例如作为 Step B 的前期铺垫，或者想确认它到底能带来多少改善），建议做以下调整：

### 4.1 用 forward stepwise selection 代替顺序贪心

当前贪心是”按候选遍历顺序，第一个不超阈值的就加”。这对 cond(S) 来说还凑合，但如果用 J 作为准则，顺序敏感性会更强。

改成每轮扫描所有剩余候选，选 J 下降最大的那个。代价 O(candidates² × eval_cost)，对 ~500 候选 + 单次评估 ~1ms，总计 ~几秒，完全可接受。

### 4.2 保留 cond(S) 硬截断作为 pre-filter

不要去掉 cond(S) 门槛。先用 cond < max_cond 预筛，通过预筛的再用 J 排序。理由：

- cond(S) check 很便宜（只算 S 本征值）
- J 的完整评估要做一次 kick + 观测量计算，贵 ~100 倍
- 用 cond 预筛可以快速砍掉数值上不可能的候选

### 4.3 目标函数简化

鉴于当前能接受的候选只有 ~14 个，目标函数不需要太复杂。建议只保留最强的 3 项：

```
J = w_dx2 * |dx2| / x2_exact + w_dE * |dE| / |E_exact| + w_cond * log10(cond(S)) / 10
```

理由：
- `dx2` 是最强的物理硬约束（exact kick 下应当不变）
- `dE` 是最直接的总误差指标
- `cond(S)` 防止病态
- `F_state`（norm fidelity）在已有 cond 预筛后信息增量不大
- `dp2` 在 exact kick 下会真实变化，其误差和 dE 高度相关，不需要单独加
- `gs_penalty` 在 14 个候选的量级下，增广后基态不太容易被搞坏（前面实验也验证了这一点：cross-basis fidelity = 1）

### 4.4 代码积木确认

以下组件已经存在且可直接复用：

| 组件 | 位置 | 状态 |
|------|------|------|
| `apply_analytic_kick()` 返回 fidelity | `kick_operator.cpp:128` | 已有 |
| `summarize_state()` 算 E, x², p² | `main.cpp` 静态函数 | 已有，需移到 .hpp |
| `kicked_exact_1particle()` 给 exact 参考 | `kicked_exact.cpp` | 已有 |
| 贪心选基框架 | `augment_basis_with_momentum()` | 已有 |
| cond(S) 计算 | 贪心循环内 | 已有 |

预估工作量：~150 行新代码 + ~30 行改动。

## 5. 真正的优先级建议

基于整个 qa.md 两轮实验 + first-kick image 失败尝试的上下文，我的优先级建议是：

### 最优先：Step B（正交化子空间原型）

这是信息价值最高的下一步。它能回答”Gaussian 子空间本身行不行”这个路线级问题。具体来说：

1. 把整个候选池（~500 个）的 overlap 矩阵构造出来
2. 做 SVD，保留奇异值 > threshold 的方向
3. 在这个正交化后的子空间里做 kick + free evolution
4. 和 exact gamma=0 比

如果正交化后 gamma=0 精度大幅改善 → 说明当前 Gaussian 候选本身没问题，只是贪心选法太差
如果正交化后仍然很差 → 说明 Gaussian 子空间本身不够，应该转向其他方法

**如果 Step B 也失败，那所有在贪心框架上的优化（包括 first-kick objective）都是白做。**

### 第二优先：First-kick objective 作为 Step B 的配套

在 Step B 的正交化子空间里，如果需要做 truncation（比如保留多少个 SVD 方向），first-kick objective 作为 truncation 准则是有价值的。此时可接受的基函数数量不再被 cond(S) 卡死在 ~14 个，J 的筛选能力才真正发挥出来。

### 第三优先：单独在当前贪心框架里做 first-kick objective

可以做，但预期改善幅度有限。更多是作为”把诊断代码重构成可复用函数”的工程价值，而不是物理精度上的突破。

## 6. 一句话总结

**First-kick 目标函数在技术上可行、代码积木齐全，但两轮实验已经清楚表明当前结构性瓶颈在”只能接受 ~14 个候选”而不在”这 14 个选得不够好”。把选基准则从 cond(S) 换成 J，最多改善 early-time 几步的精度，不能从根本上改变长期 kicked dynamics 的质量。它的真正价值在于作为 Step B（正交化子空间）的配套工具——在正交化消除了 cond(S) 瓶颈之后，才能在更大的可选集上发挥筛选作用。**

---

# Q&A: 方案 A 是不是只是把截断 `n` 取大一些（2026-04-01）

不是，**不只是把截断 `n` 取大一些**。

这两件事表面上有点像，但数学上不是一回事。

## 1. “把截断 `n` 取大一些”是什么意思

这指的是：

- 还是沿用现在的“单个谐波 candidate”思路
- 只不过允许更多动量标签

例如原来只考虑：

`n = 0, ±1, ±2`

现在改成：

`n = 0, ±1, ±2, ±3, ±4, ...`

这种做法本质上还是：

> 一个 candidate 对应一个单独的谐波模式

它只是把候选库范围扩大。

## 2. 方案 A 真正想做的是什么

方案 A 的核心不是“多几个单谐波模式”，而是：

> **把一个 candidate 做成一个完整的 kick packet**

也就是一个 candidate 本身就是：

`chi_j^(L) = sum_{|n|<=L} c_n * exp(i 2 n k_L x) * phi_j`

这里：

- `phi_j` 是某个原始 basis function
- `c_n = (-i)^n J_n(kappa)` 是 kick 的 Bessel 权重
- `L` 是 packet 内部的截断阶数

所以一个 candidate 里同时已经包含：

- `n = 0`
- `n = ±1`
- `n = ±2`
- ...

这些谐波，而且相对权重不是手调的，而是直接由 kick 算子决定。

## 3. 为什么这和“单纯把 `n` 取大”不一样

### 方法 1：只把 `n` 取大

你得到的是很多分开的 candidate：

- `exp(i 2 k_L x) * phi_j`
- `exp(i 4 k_L x) * phi_j`
- `exp(i 6 k_L x) * phi_j`
- ...

然后再让算法自己线性组合它们。

### 方法 2：kick packet

你直接给一个 candidate：

`J_0 phi_j + (-i) J_1 e^{i 2 k_L x} phi_j + (-1) J_2 e^{i 4 k_L x} phi_j + ...`

这个 candidate 本身就已经很像：

`U_kick phi_j`

也就是说：

- 方法 1 是“给零散零件，让系统自己拼”
- 方法 2 是“直接给一个按 kick 结构拼好的小模块”

## 4. 为什么我会觉得方案 A 可能更适合降 `dx2`

因为 exact 第一脚 kick 做的事情是：

`psi(x) -> exp[-i kappa cos(2 k_L x)] * psi(x)`

这相当于：

- 原来的空间包络不变
- 只是乘上一个完整的相位因子

如果你只用单个谐波去逼近，它更容易：

- 抓住某些动量变化
- 但把位置密度扭歪

而如果你给的是一个完整的 kick packet，它从一开始就更接近：

`U_kick phi_j`

因此更有希望守住：

- 位置密度
- 也就是更有希望把 `dx2` 压下来

## 5. 所以方案 A 和“把 `n` 取大”之间的关系

关系是：

- 方案 A **也需要一个截断 `|n| <= L`**
- 但这个截断只是 packet 内部的级数截断
- 它的重点不是“截断取多大”
- 而是“candidate 的组织方式已经从单个谐波，变成了一个整包谐波”

所以更准确地说：

> **方案 A 是“有截断的 kick packet candidate”，而不是“单纯把单个谐波的截断范围变大”。**

## 一句话总结

**不是简单把 `n` 截断取大，而是把一个 candidate 从“单个谐波模式”升级成“按 Bessel 权重打包好的多谐波 kick image”。**

## 已做的一次实验：方案 A 的最小 packet 原型

已经实际做了一个最小版的 Scheme A 原型：

- 不再用“单个 repair candidate”
- 而是用 **packet block**
- 每个 packet 对应：
  - 一个重要的原始 basis function
  - 加上一组按 kick 结构组织的低阶 image functions

### 这次原型的设置

- 主干 basis 保持不动
- 只测试少量 packet block
- packet 总数：`32`
- 每个 packet 是一小组 first-kick image functions 的组合块

### 实际结果

运行结果是：

```text
--- First-Kick Packet Search ---
packet count = 32
no improving repair packet found after 0 accepted packets
```

对应的第一脚瞬时诊断完全没有变化：

```text
Projection fidelity = 0.9989555239
dE   = 9.3494573850e-02
dx2  = 2.7879977035e-02
dp2  = 1.5910917067e-01
```

## 这说明什么

这说明至少在“当前这个最小 packet 构造方式”下：

- packet 的想法本身并没有立刻带来改善
- 当前 objective 甚至没有愿意接受任何一个 packet block

所以当前应得出的结论是：

> **Scheme A 的直觉并不荒谬，但这个最小版 packet 原型没有比现有启发式 basis 更好。**

更进一步说，它提示我们：

- 要么当前 packet 的构造方式还太粗
- 要么当前 fixed-basis + greedy 框架并不擅长利用这种 packet block

无论哪一种，都说明：

> **至少“最小版 Scheme A”目前没有证据支持它能降低 first-kick 的 `dx2`。**

---

# Q&A: 这个误差到底怎么改善？是不是数学表示本身有问题（2026-04-01）

## 当前最核心的判断

到目前为止，所有实验给出的信号都越来越一致：

> **当前问题更像是“表示空间 + 投影结构”的问题，而不是某个简单公式写错了。**

更准确地说：

- `U_kick` 本身没有问题
- `U_kick^dagger x^2 U_kick = x^2` 也没有问题
- 问题出在：

> **我们选择的有限维 basis/manifold 不是 kick 作用下近似封闭的，因此把 exact kicked state 压回这个空间时，会不可避免地扭曲态。**

## 1. 为什么说更像“表示问题”

因为现在已经试过很多种“小修小补”：

1. 原始启发式动量基
2. raw first-kick image candidates
3. 混合 repair layer
4. 扩大 repair pool
5. 提高 `dx2` 权重
6. 最小 packet 原型

这些尝试的结果都说明：

- 误差能小幅动一动
- 但很快就遇到瓶颈

这类行为更像：

> “你在同一个不太合适的表示空间里，反复重排基向量”

而不是：

> “只要把某个公式修对就会明显变好”

## 2. 这个“表示问题”具体指什么

当前 fixed-basis 方法假设：

`U_kick psi`

可以被当前增广基张成的有限维子空间很好地承接。

但从 first-kick 的结果看：

- exact 第一脚后 `<x^2>` 应该严格不变
- 而 projection 后 `<x^2>` 已经变了

这说明当前子空间对 `U_kick psi` 的表示不是“只差一点”，而是：

> **它连最基本的“纯相位 kick 不应改变位置密度”这个结构都没有守住。**

所以问题本质上是：

> **当前 basis family 对 kicked state 的几何结构不对。**

## 3. 所以误差真正要怎么改善

从现在的证据看，后续的有效改善大致只有三类。

### 路线 A：换“表示方式”

也就是承认：

> 当前这套 raw ECG fixed-basis family 并不适合承接 repeated kicks

那就改用更自然的表示：

- gamma=0：exact grid / harmonic oscillator eigenbasis / plane-wave-like basis
- gamma>0：再考虑 ECG 或混合基

这是最根本的改善。

### 路线 B：保留候选 family，但改“数值组织方式”

这就是正交化子空间那条线。

它不是换掉物理候选，而是：

- 先把冗余方向压掉
- 再在真正独立的子空间里传播

如果这条线成功，说明：

- 问题主要是 raw nonorthogonal 表示太差

如果这条线也不成功，说明：

- 问题更接近 candidate family 本身不够

### 路线 C：允许 basis 真正随 kick 动态更新

也就是不再坚持：

- 一次增广
- 全程固定 basis

而是承认：

> kicked state 的有效支撑在不断变化，basis 也必须跟着变

这条线理论上很合理，但工程代价最高。

## 4. 哪些办法大概率只是局部缓解

从现在的实验看，这些方法更像“局部缓解”，而不是根本解决：

- 单纯扩大 raw candidate 数量
- 单纯调目标函数权重
- 单纯加几个 first-kick repair basis
- 单纯做更花的贪心筛选

它们可以改善一点点，但很难跨过当前的结构性瓶颈。

## 5. 所以现在最准确的一句话

不是说：

> “数学推导错了”

而是说：

> **当前我们用来承接 kicked state 的有限维数学表示不够合适，因此投影本身就会系统性地产生物理量误差。**

## 一句话总结

**这个误差的根源更像“表示空间不对”，不是“公式错了”；要真正把误差压下去，最终要么换更合适的表示，要么至少把当前候选池放进更合适的数值组织方式（例如正交化子空间），再不然就必须让 basis 随 kick 动态更新。**

---

# Q&A: 当前最关键的判断：问题不在演化，而在 kick 后瞬时表示失败（2026-04-01）

这个判断现在已经可以说得很明确：

> **当前主要问题不是后续 free evolution 累积误差，而是在第一次 kick 后那一瞬间，当前 basis 就已经不能很好地表示 kicked state。**

## 证据

第一次 kick 后、还没做自由演化时，已经看到：

- `Projection fidelity = 0.998955...`
- `dE  = 9.35e-2`
- `dx2 = 2.79e-2`
- `dp2 = 1.59e-1`

其中最关键的是：

- exact kick 后 `<x^2>` 本来应严格不变
- 但当前投影后 `<x^2>` 已经变了

这说明：

> **问题在 kick 后瞬时投影这一步就已经出现，而不是要等演化很多步之后才暴露出来。**

## 这意味着什么

这意味着当前 basis 的问题不是“长期动力学支撑不够”这么简单，而是更前面的一层：

> **当前 basis family 连“第一次 kicked image”都没有很好地承接住。**

所以当前更准确的任务不是：

- “怎么改善后续演化”

而是：

- “怎么构造一个能正确承接 `U_kick psi` 的表示空间”

## 这会怎么改变后续思路

既然问题在 kick 瞬间就已经存在，那么：

1. **先修 first-kick representation** 的优先级最高
2. 如果第一脚都接不住，后面的 repeated-kick basis update 讨论都建立在不稳的基础上
3. 也就是说，后续路线应先回答：

> “什么样的 basis / subspace 能让第一脚 kick 后的态不扭曲位置密度？”

## 一句话总结

**对，当前最该承认的事实是：问题首先不在演化，而在于当前 basis 没有很好地表示 kick 之后的情况；如果第一脚都接不住，后面的所有长期动力学问题都会被放大。**

---

# Q&A: “最优先 Step B（正交化子空间）”这段判断到底对不对（2026-04-01）

## 结论

### 大方向：对

把 **Step B（正交化子空间原型）** 放在高优先级，方向上是对的。

因为它确实能回答一个非常关键的问题：

> **在当前这套候选库里，如果把非正交冗余和病态 `S` 的问题先去掉，误差有没有可能明显降下来？**

换句话说，它是在测试：

- 问题主要是“当前贪心选基 + 非正交病态”造成的
- 还是“当前候选家族本身就不够表达 kicked state”

这点判断是成立的。

### 但表述里有一句话太满

原文这句：

> **如果 Step B 也失败，那所有在贪心框架上的优化（包括 first-kick objective）都是白做。**

这个说法**方向上接近，但严格说太绝对了**。

更准确的版本应该是：

> **如果 Step B 也失败，那么继续在“当前候选家族 + 当前固定基框架”里做贪心小修小补，回报大概率会很有限。**

而不是说“所有后续工作都白做”。

## 为什么 Step B 很有价值

因为它给了我们一个很有用的“上界测试”：

### 当前贪心框架在做什么

现在我们是：

- 先生成一堆 raw candidate
- 由于 `S` 病态
- 最后只能保留很少一部分

所以当前结果其实受两种因素混在一起影响：

1. 候选是否足够好
2. 这些候选能不能在数值上安全地一起使用

### Step B 在做什么

Step B 的核心就是：

- 先不急着贪心删掉候选
- 先对整个候选池做正交化 / SVD 截断
- 再在这个“去冗余后的子空间”里做 kick 和传播

这样就把问题改写成：

> **如果候选家族里真的有足够的信息，那么在正交化后，它们应该能更充分地发挥出来。**

所以 Step B 本质上是在测：

> **当前候选池在“最佳数值组织方式”下，最多能把误差降到什么程度。**

这就是为什么它是一个非常有价值的判别实验。

## 它到底能回答什么

Step B 真正能回答的是：

### 如果成功

如果正交化子空间后，gamma=0 的 first-kick 和后续 kicked dynamics 明显变好，那么说明：

- 当前候选家族本身不是主要问题
- 主要问题在于：
  - raw nonorthogonal 表示
  - 贪心筛选
  - 病态 `S`

这时可以说：

> **理论上，在当前候选家族内部，误差确实有明显下压空间。**

### 如果失败

如果正交化之后仍然明显不好，那么说明：

- 至少在**当前候选家族 + 固定基表示**这个框架里
- 问题不是单纯数值病态
- 而是候选方向本身就不够承接 kicked state

这时我们就会更有理由转向：

- 新的 candidate 类型（例如 kick packet）
- 动态 basis update
- 甚至更根本的表示方式

## 它不能回答什么

Step B 不能严格回答下面这些更强的问题：

### 不能证明“所有 Gaussian 路线理论上都不行”

因为 Step B 用的仍然是：

- 当前这批候选的 span
- 当前固定基 / 截断框架

如果它失败，最多说明：

> **当前候选池的 span 不够，或者当前固定基思路不够。**

它不能证明：

- 所有可能的 Gaussian candidate family 都不行
- 所有动态 Gaussian 方法都不行
- 所有 packetized / adaptive Gaussian 方法都不行

### 不能直接给出“最终最优答案”

它更像一个路线判别器，而不是最终物理解法。

## 所以最准确的说法是什么

我认为最准确的表述应该是：

> **Step B 可以用来判断：在当前候选家族里，如果先消除非正交冗余和病态问题，误差理论上还能不能明显往下压。**

以及：

> **如果 Step B 成功，说明当前候选家族有潜力；如果 Step B 失败，说明问题不只是贪心筛选，而更可能出在候选家族本身或固定基框架。**

## 对那段“真正的优先级建议”的最终评价

### 我认为正确的部分

- 把 Step B 放高优先级：对
- 把它看成路线判别实验：对
- 认为它的信息价值高于继续堆小修小补：对

### 我认为需要收一收的部分

- “如果 Step B 失败，那所有贪心优化都是白做”

更好的说法是：

- **如果 Step B 失败，那么继续在当前 raw candidate 家族里做贪心小修小补，大概率收益有限。**

## 一句话总结

**可以，正交化子空间确实能用来判断“在当前候选家族里，理论上还有没有把误差明显压下去的空间”；但它给的是“当前候选池的上界测试”，不是“所有 Gaussian 路线的终极判决”。**
