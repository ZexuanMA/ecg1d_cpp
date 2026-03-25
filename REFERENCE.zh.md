# ECG1D 项目参考手册

汇总自项目开发过程中的全部技术文档。原始文件保留在 `past_md/` 中。

---

## 目录

1. [方法论基础：Varga/Suzuki SVM](#1-方法论基础vargasuzuki-svm)
2. [我们的基函数形式](#2-我们的基函数形式)
3. [广义本征值问题与变分原理](#3-广义本征值问题与变分原理)
4. [能量方差作为质量指标](#4-能量方差作为质量指标)
5. [Ghost State 事件：调查与修复](#5-ghost-state-事件调查与修复)
6. [Benchmark 结果](#6-benchmark-结果)
7. [Delta 接触势的特殊困难](#7-delta-接触势的特殊困难)
8. [N 粒子扩展性](#8-n-粒子扩展性)
9. [项目路线图](#9-项目路线图)
10. [关键参考文献](#10-关键参考文献)

---

## 1. 方法论基础：Varga/Suzuki SVM

参考书：Suzuki & Varga, *Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems* (Springer, 1998)

### 1.1 核心思想

- 用**显式关联高斯函数 (ECG)** 作为变分基底
- 用**随机试探 + 选择**优化非线性参数（SVM = Stochastic Variational Method）
- 线性系数通过广义本征值问题自动确定
- 适用于 few-body (N=2~7) 的高精度 bound-state 计算

### 1.2 SVM 三阶段流程（书 p.51-56）

```
1. SVM Growth (competitive selection): K=1 → K=K_target
   - 逐个添加基函数，每个从 n_trials 个随机候选中选最优

2. SVM Refinement: 多轮，直到能量不再改善
   - 逐个替换每个基函数，保留能量更低的替代

3. Fine tuning (TDVP/Powell): 连续微调非线性参数
```

**Refine 对精度的贡献是决定性的**（书 Table 4.6: K=100 refine 10 轮将误差从 1.2e-5 降到 6.0e-7）。

### 1.3 线性相关问题（书 p.48-50）

书中明确警告："the overlap of the random basis functions should be smaller than a prescribed limit"。我们实现了 `s_well_conditioned()` 检查 S 矩阵条件数（全局度量，优于 pair-wise overlap 检查）。

### 1.4 ECG 基函数类型

| | Varga Type I | Varga Type II | 我们的参数化 |
|---|---|---|---|
| 坐标 | Jacobi | 单粒子 | 单粒子 |
| 矩阵大小 | (N-1)×(N-1) | N×N | N×N |
| 参数 | 实对称正定 A | 实 A+B+R | **复** A+B+R |
| 用途 | 静态 | 静态 | 静态 + 动力学 |

**我们的参数化是 Type II 的复数推广，用于支持 TDVP 实时间演化。**

### 1.5 对当前项目的关键启发

1. **SVM 基底构建是主战场**——TDVP 是后期抛光，不能替代前期的 basis search
2. **静态优化应只用实参数（B=0, R=0）**——复数参数浪费搜索空间
3. **Refine 必须做**——SVM growth 后的基函数质量不够，refine 是关键补充
4. **物理启发式约束正常且必要**——对 bosonic ground state 做对称采样、对 trap 做中心在原点的采样

---

## 2. 我们的基函数形式

### 2.1 参数（`basis_params.hpp`）

| 参数 | 类型 | 含义 |
|------|------|------|
| u | 复标量 | 线性组合系数 |
| A | N×N 复矩阵 | 二次型主要部分（宽度/关联） |
| B | N×N 对角复矩阵 | 额外二次型 + 位移耦合 |
| R | N 维复向量 | 位移/中心参数 |

### 2.2 Ket 侧基函数

$$\phi_k(\mathbf{x}) = \exp\!\Big(-\mathbf{x}^T (A_k + B_k)\,\mathbf{x} + 2\,\mathbf{R}_k^T B_k\,\mathbf{x} - \mathbf{R}_k^T B_k\,\mathbf{R}_k\Big)$$

Bra 侧取复共轭。

### 2.3 Overlap 矩阵元（`pair_cache.cpp`）

两个高斯指数合并后做 N 维高斯积分：

$$K = \overline{A_i + B_i} + P^T(A_j + B_j)P, \quad M_G = \pi^{N/2} \det(K)^{-1/2} \exp(C + b^T K^{-1} b / 4)$$

其中 P 是置换矩阵，b 和 C 来自线性项和常数项的合并。

### 2.4 B 和 R 的角色

| 阶段 | B | R | 说明 |
|------|---|---|------|
| 静态基态 | **= 0** | **= 0** | 基态是实函数、零动量、原点对称 |
| TDVP 动力学 | 非零 | 非零 | 编码动量和相位，由方程演化 |

**关键：静态阶段的 SVM/Refine 应该冻结 B=0, R=0**，只搜索 A。这减少搜索维度（N=2: 从 7 个参数降到 1 个），大幅提升 basis 质量。

---

## 3. 广义本征值问题与变分原理

### 3.1 来源

```
H|ψ⟩ = E|ψ⟩     ← 薛定谔方程（抽象 Hilbert 空间）
    ↓ 选基底展开
Hu = ESu          ← 广义本征值问题（非正交基底）
    ↓ S = I 时
Hu = Eu           ← 普通本征值问题（特例）
```

路线 A（投影方程）和路线 B（Rayleigh 商 $R(u) = u^\dagger Hu / u^\dagger Su$ 求极值）都给出 $Hu = ESu$。

### 3.2 本征态是 Rayleigh 商的驻点

$R(u)$ 在方向空间上的地形：全局最小值 = 基态 $E_0$，全局最大值 = 最高本征值，中间是鞍点（激发态）。

变分原理 = 加权平均不小于最小值：$R(u) = \sum_n |c_n|^2 E_n / \sum_n |c_n|^2 \geq E_0$。

### 3.3 标准解法（S^{-1/2} 变换）

1. S 谱分解 → 截断小特征值 → 构造 $S^{-1/2}$
2. $\tilde{H} = S^{-1/2} H S^{-1/2}$ → 标准对角化
3. 回代 $u = S^{-1/2} \tilde{u}$

S 近奇异时，$S^{-1/2}$ 在差态方向上放大噪声（放大因子 $1/\sqrt{w_\min}$），但我们的实验证明**当前 w_min (~1e-5) 只导致 ~1e-12 的噪声，远不是问题**。

---

## 4. 能量方差作为质量指标

### 4.1 定义

$$\sigma^2 = \langle H^2 \rangle - E^2 = \sum_n |c_n|^2 (E_n - E)^2$$

$\sigma^2 = 0 \iff |\psi\rangle$ 是精确本征态。

### 4.2 实际计算

在 $S^{-1/2}$ 变换空间中：$\sigma^2 = \|\tilde{H}\tilde{c}\|^2 - E_0^2$（一次矩阵-向量乘法）。

### 4.3 Weinstein 定理

至少有一个精确本征值在 $[E-\sigma, E+\sigma]$ 内。

### 4.4 Temple 下界

$$E_0 \geq E - \frac{\sigma^2}{\varepsilon - E}$$

其中 $\varepsilon \leq E_1$（第一激发态能量）。给出 $O(\sigma^2)$ 精度的下界，比 Weinstein 紧得多。

---

## 5. Ghost State 事件：调查与修复

### 5.1 现象

Stochastic refine 导致能量从 1.527 跌到 1.000（低于真实基态 1.527），违反变分原理。

### 5.2 排除过程

| 假设 | 结果 | 证据 |
|------|------|------|
| A: S 近奇异 → 噪声放大 | **排除** | w_min=6.7e-5，噪声仅 2e-12 |
| B: `.real()` 丢虚部 | **排除** | 改复 Hermitian 后行为完全相同 |
| C: incremental H/S 更新 bug | **已修复** | 改为全量 `build_HS` 重建 |
| D: H 矩阵元不 Hermitian | **确认为症状** | A 不对称 → K 不对称 → kernel 不对称 |
| **E: Eigen aliasing in perturb_basis** | **根本原因** | `.eval()` 一行修复 |

### 5.3 根本原因

`perturb_basis()` 中：

```cpp
A_new = 0.5 * (A_new + A_new.transpose());  // BUG: Eigen aliasing
// 修复：
A_new = (0.5 * (A_new + A_new.transpose())).eval();
```

Eigen 的惰性求值导致 `A_new` 在被写入时 `.transpose()` 同时被读取，结果不对称。

### 5.4 因果链

A 不对称 → K 矩阵不对称 → K_inv 不对称 → H kernel 不满足 Hermitian 对称性 → build_HS 用 conj 填下三角掩盖了不对称 → H 矩阵不代表真实 Hamiltonian → 变分原理被打破 → refine 找到 E < E_true

### 5.5 附带改进

| 改动 | 说明 |
|------|------|
| 复 Hermitian 求解 | `lowest_energy_full` 用 `SelfAdjointEigenSolver<MatrixXcd>`，去掉 `.real()` |
| 全量重建 | `stochastic_refine` 用 `build_HS` 替代增量行/列更新 |
| S 条件数检查 | `s_well_conditioned()` 替代 pair-wise overlap 检查 |

### 5.6 教训

**Eigen 中 `A = f(A)` 形式的自赋值必须加 `.eval()`。** 尤其是涉及 `.transpose()`、`.conjugate()` 等惰性操作时。

---

## 6. Benchmark 结果

### 6.1 1-particle harmonic oscillator (E_exact = 0.5)

| 版本 | TDVP 误差 |
|------|-----------|
| Phase A 前 | 2.2e-6 |
| Phase A 后 | **1.8e-12** |

### 6.2 2-particle Gaussian interaction (E_exact = 1.5266998310)

| 阶段 | K | 误差 |
|------|---|------|
| SVM | 20 | 3.0e-6 |
| + Refine 10 轮 | 20 | ~2.8e-6 |

精度瓶颈在 `random_basis_2particle` 的采样分布，不在 K 大小。增加 K（至 30）、增加 n_trials（至 20000）、设 B=0/R=0 均未突破 2.5e-6。

### 6.3 2-particle delta contact (E_exact = 1.3067455)

| 阶段 | K | 误差 |
|------|---|------|
| SVM | 15 | 2.26e-3 |
| + Refine 10 轮 | 15 | 2.26e-3 |
| TDVP | — | 发散（C 矩阵条件数 3e+11） |

误差卡在 2.3e-3 是**代数收敛（K^{-1/2}）的预期行为**，不是 bug。Gaussian 基底逼近 delta cusp 效率极低。

### 6.4 2-particle kicking (E_exact = 2.4452547216)

精确值由有限差分法计算：$E = 2 \times E_\text{single}$，其中 $H_\text{single} = -\frac{1}{2}\frac{d^2}{dx^2} + \frac{1}{2}x^2 + \cos(x)$（无粒子间相互作用）。

| 阶段 | K | 误差 |
|------|---|------|
| SVM | 20 | 3.7e-5 |
| + Refine 10 轮 | 20 | 1.8e-5 |
| + TDVP 5 步 | 20 | 1.5e-5 |

Kicking 比 Gaussian interaction 难约一个量级（3.7e-5 vs 3.0e-6），因为 cos(x) 势需要更多基函数来表达非高斯形状。Refine 效果显著（误差砍半），TDVP 每步改善约 5e-7。

### 6.5 精度总结

| Benchmark | E_exact | SVM 误差 | + Refine | 收敛特性 |
|-----------|---------|---------|----------|---------|
| 1-particle harmonic | 0.5 | — | — | 1.8e-12（TDVP） |
| Gaussian interaction | 1.5267 | 3.0e-6 | 2.8e-6 | 平滑势，指数收敛 |
| Kicking | 2.4453 | 3.7e-5 | 1.8e-5 → 1.5e-5 (TDVP) | cos 势，收敛较慢 |
| Delta contact | 1.3067 | 2.3e-3 | 2.3e-3 | cusp，代数收敛 K^{-1/2} |

**对 MBDL 动力学，Gaussian 和 Kicking 的精度（1e-5 ~ 1e-6）足够做初态准备。**

---

## 7. Delta 接触势的特殊困难

### 7.1 收敛率

| 势函数 | 收敛率 | K=13 精度 |
|--------|--------|-----------|
| Gaussian（平滑）| 指数 $e^{-cK^\alpha}$ | 2.5e-6 |
| Delta（cusp）| 代数 $K^{-1/2}$ | 2.3e-3 |

### 7.2 文献方案

| 方案 | 收敛提升 | 实现难度 |
|------|---------|---------|
| 增加 K（暴力） | 无（仍是 K^{-1/2}）| 低 |
| r_ij 前因子 ECG | 大 | 中（新矩阵元公式）|
| Transcorrelated | 极大（K^{-3}）| 高（非厄米 H）|
| **有限程近似 + 外推** | 恢复指数 | **低（已有实现）** |

### 7.3 对 MBDL 项目的实用策略

Delta 精度不是当前瓶颈。实验中的冷原子相互作用本身就是有限程的。用 Gaussian 势做动力学是物理上合理的选择。

---

## 8. N 粒子扩展性

### 8.1 计算量

| N | N! | 单矩阵元代价 | K 需求 | 典型精度 |
|---|-----|-------------|--------|---------|
| 2 | 2 | 基准 | 10-20 | 1e-6 ~ 1e-10 |
| 3 | 6 | 3× | 50-200 | 1e-4 ~ 1e-6 |
| 4 | 24 | 12× | 200-500 | 1e-3 ~ 1e-4 |
| 5-6 | 120-720 | 60-360× | 500-5000 | 1e-2 ~ 1e-3 |
| 10 | 3.6M | 不可行 | — | — |

### 8.2 瓶颈

显式 N! 置换求和是硬性瓶颈。N=10 需要完全不同的方法（利用 Bose 对称性减少置换数、cluster 分解等）。

### 8.3 实际策略

先做 N=2 的完整动力学，再 N=3, 4。如果 N=2~4 已经看到 MBDL 特征信号，就是有价值的理论贡献。

---

## 9. 项目路线图

### 当前状态

- [x] Phase A: 静态优化器（SVM + Refine + TDVP）
- [x] Ghost state 修复（Eigen aliasing）
- [ ] **SVM 用 B=0, R=0 提升基底质量** ← 当前任务
- [ ] Phase C: 实时间 kicked dynamics
- [ ] Phase D: N=3 扩展

### 优先级

```
[立即] B=0, R=0 优化：SVM 基底质量提升（Gaussian 1e-8 目标）
[高]   Phase C: 实时间 TDVP + 周期 kick + observables
[中]   N=3: 通用 random_basis_Nparticle
[低]   Delta 精度提升（有限程外推或 cusp-aware 基底）
```

### Phase C 核心需求

论文目标（MBDL 模拟）需要：
1. 实时间 TDVP（不是 imaginary-time 能量最小化）
2. 显式时间依赖 Hamiltonian（kick 脉冲序列）
3. 观测量：n(k,t)、能量吸收、entropy、G^(1)(z,t)

---

## 10. 关键参考文献

### 书

1. Suzuki & Varga, *Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems*, Springer (1998)

### 原始论文

2. Varga & Suzuki, Phys. Rev. C **52**, 2885 (1995) — SVM 原始论文
3. Varga & Suzuki, Phys. Rev. A **53**, 1907 (1996) — PRA 版本

### 综述

4. Mitroy et al., Rev. Mod. Phys. **85**, 693 (2013) — ECG 综述

### 1D ECG 矩阵元

5. Zaklama et al., Few-Body Syst. **61**, 6 (2020) — 1D ECG 矩阵元公式

### Delta/Cusp 处理

6. Jeszenszki et al., Phys. Rev. A **98**, 053627 (2018) — Transcorrelated 方法
7. Jeszenszki et al., Phys. Rev. Research **2**, 043270 (2020) — 消除奇异性
8. Koscik, Phys. Lett. A (2018) — 接触势优化基底
9. Cencek & Kutzelnigg, Chem. Phys. Lett. (2004) — cusp 条件高斯基底

### 梯度方法

10. arXiv:2408.08522 (2024) — SVM + GVM 混合最优
11. arXiv:2512.07323 (2025) — 扩展 SVM（muonic 系统）

### 目标论文

12. Guo et al., Science **389**, 716 (2025) — MBDL 实验观测
