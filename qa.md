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
