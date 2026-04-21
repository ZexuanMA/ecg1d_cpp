# N=2 无相互作用 Kicked Dynamics：第一个里程碑

## 1. 为什么是这个问题

### 1.1 目标：验证实时间演化的完整管线

我们最终要模拟的是 Science 论文中的 many-body dynamical localization (MBDL)：

$$
H(t) = \sum_i \left[\frac{p_i^2}{2m} + V_\text{trap}(z_i) + \hbar\kappa\cos(2k_L z_i)\sum_n\delta(t-nT)\right] + g_{1D}\sum_{i<j}\delta(z_i-z_j)
$$

这里面有四个核心模块需要验证：

| 模块 | 说明 |
|------|------|
| 实时间 TDVP | $C\dot{z} = -ig$，不是虚时间 |
| 周期性 kick | 瞬时施加 $\kappa\cos(2k_L z)$ |
| 自由演化 | kick 之间的 $H_\text{free}$ 演化 |
| 观测量 | E(n), n(k,n) 等随 kick 数的变化 |

**如果我们直接上完整问题（有相互作用 + 多粒子），出了错就不知道是哪个模块的问题。**

### 1.2 为什么 γ=0（无相互作用）是最好的起点

设 $g_{1D} = 0$，两个粒子之间没有相互作用。此时：

$$
H(t) = \underbrace{\left[\frac{p_1^2}{2m} + V(z_1) + \kappa\cos(2k_L z_1)\sum_n\delta(t-nT)\right]}_{H_1(t)} + \underbrace{\left[\frac{p_2^2}{2m} + V(z_2) + \kappa\cos(2k_L z_2)\sum_n\delta(t-nT)\right]}_{H_2(t)}
$$

**两个粒子完全独立！** $N$-体波函数分解为：

$$
\Psi(z_1, z_2, t) = \psi(z_1, t) \cdot \psi(z_2, t)
$$

其中 $\psi(z, t)$ 满足单粒子含阱 kicked rotor (KR) 方程：

$$
i\hbar\frac{\partial\psi}{\partial t} = \left[\frac{p^2}{2m} + \frac{1}{2}m\omega^2 z^2 + \hbar\kappa\cos(2k_L z)\sum_n\delta(t-nT)\right]\psi
$$

### 1.3 三个关键优势

**优势 1：精确解可算。** 单粒子问题可以用有限差分/分裂算符法在位置空间精确求解。N=2 的精确结果就是单粒子结果乘以 2：

$$
E_2(n) = 2 \times E_1(n), \quad n_2(k, n) = n_1(k, n) * n_1(k, n)
$$

这给了我们一个**逐 kick 逐时刻**的精确参考，不只是终态能量。

**优势 2：物理已知。** 单粒子 quantum kicked rotor (QKR) 的行为是量子混沌领域最经典的结果之一：
- 短时间：能量线性增长（经典扩散）
- 长时间：能量饱和（dynamical localization，Anderson 局域化在动量空间的类比）
- 动量分布 n(k) 冻结为指数衰减形状

我们知道该看到什么，如果看不到就说明代码有问题。

**优势 3：隔离误差来源。** γ=0 排除了相互作用项，任何误差都只能来自：
- 实时间 TDVP 本身的数值精度
- Kick 施加方式
- 基底表达力
- 长时间相位累积稳定性

这正是我们需要验证的模块。

---

## 2. 物理图像：Quantum Kicked Rotor

### 2.1 经典 kicked rotor

考虑一个被周期性踢的粒子（先忽略阱）：

$$
H = \frac{p^2}{2} + K\cos(x)\sum_n\delta(t - nT)
$$

经典动力学是标准映射（Chirikov standard map）：

$$
p_{n+1} = p_n + K\sin(x_n), \quad x_{n+1} = x_n + p_{n+1}T
$$

当 $K$ 足够大时，动量扩散：$\langle p^2 \rangle \propto n$（线性增长），类似随机行走。

### 2.2 量子 kicked rotor：dynamical localization

量子力学中，由于量子干涉效应，动量扩散**在一定 kick 数后停止**：

$$
\langle p^2(n) \rangle \approx \begin{cases}
D \cdot n & n < n^* \quad \text{（经典扩散阶段）} \\
\text{const} & n > n^* \quad \text{（局域化阶段）}
\end{cases}
$$

其中 $n^* \sim D / \hbar^2$ 是局域化时间，$D$ 是经典扩散系数。

**动量分布的特征：**

$$
n(k) \propto \exp\left(-\frac{|k|}{k_\text{loc}}\right)
$$

指数衰减，特征宽度 $k_\text{loc}$ 在局域化后不再增长。

### 2.3 加谐振子阱的影响

我们的系统还有一个谐振子阱 $V = \frac{1}{2}\omega^2 z^2$。阱的效果是：

- 阱束缚了粒子，位置空间不能无限扩展
- 阱的频率 $\omega$ 与 kick 周期 $T$ 的关系决定了准能量谱的结构
- 当 $\omega T / (2\pi)$ 是无理数时，行为最接近"标准" QKR

对我们来说，阱的存在是好事——它保证波函数始终局域，Gaussian 基底有限个就能表示。

---

## 3. 我们现有的代码基础

### 3.1 已经实现的

| 模块 | 文件 | 状态 |
|------|------|------|
| 实时间 TDVP step | `kicked_evolution.cpp:realtime_tdvp_step()` | 已有 |
| 自由演化 | `kicked_evolution.cpp:free_evolve()` | 已有 |
| 施加 kick | `kicked_evolution.cpp:apply_kick()` | 已有 |
| kick 循环 | `kicked_evolution.cpp:kicked_evolution()` | 已有 |
| 观测量记录 | `kicked_evolution.cpp:Observable` | 已有（E, T, norm）|
| 测试入口 | `main.cpp:run_kicked_evolution_test()` | 已有 |

### 3.2 实时间 TDVP 与虚时间的区别

虚时间（能量最小化）：

```cpp
// tdvp_solver.cpp
VectorXcd rhs = -g_bar_update;    // C dz = -g
```

实时间（薛定谔演化）：

```cpp
// kicked_evolution.cpp
VectorXcd rhs = Cd(0, -1) * g_bar_update;    // C dz = -ig
```

就是多了一个 $-i$ 因子。但**没有线搜索**——实时间演化中能量可以涨可以跌，不能用"能量必须下降"来验证步长。

### 3.3 当前 `run_kicked_evolution_test` 的流程

```
1. SVM + Refine：准备 N=2 + delta 的基态（K=15）
2. 虚时间 TDVP：100 步 polish
3. 切换到实时间：10 个 kick，T=1.0，ideal delta kick
4. 输出：(time, energy, kinetic_energy, overlap_norm)
```

### 3.4 需要调整的

| 项目 | 现状 | 目标 |
|------|------|------|
| 相互作用 | delta（g=1.0） | **γ=0（无相互作用）** |
| 基态准备 | K=15 SVM + refine | K=20 SVM (B=0,R=0) + refine |
| TDVP 参数 | optimize_B=false, optimize_R=false | **dynamics 模式：optimize_B=true, optimize_R=true** |
| Kick 数 | 10 | **50-100**（需要看到饱和） |
| dt | 0.01 | 需要收敛性测试 |
| 验证 | 无精确解对比 | **有限差分精确解对比** |
| 观测量 | E, T, norm | **+ n(k), ΔE** |

---

## 4. 精确解的计算方法

### 4.1 单粒子有限差分法

在位置空间离散化：z_j = z_min + j * dz, j = 0, 1, ..., N_grid - 1。

H_free = -(1/2) d^2/dz^2 + (1/2) omega^2 z^2

用中心差分近似二阶导数，得到 N_grid x N_grid 的三对角矩阵。

### 4.2 Kick 的施加

#### 4.2.1 位置空间（有限差分精确解）

理想 delta kick 等价于瞬时乘以相位：

$$\psi(z, t^+) = e^{-i\kappa\cos(2k_L z)} \cdot \psi(z, t^-)$$

这在位置空间是**逐点乘法**，不涉及矩阵运算。

#### 4.2.2 ECG 基底中的解析 kick（Jacobi-Anger 展开）

对 ECG 基底，kick 不是逐点乘法——需要算 kick 矩阵元 $K_{ij} = \langle\phi_i|U_{\text{kick}}|\phi_j\rangle$。

**第 1 步：Jacobi-Anger 展开。** 把指数中的 cos 展开为平面波求和：

$$e^{-i\kappa\cos\theta} = \sum_{n=-\infty}^{\infty} (-i)^n J_n(\kappa) \, e^{in\theta}$$

其中 $J_n$ 是第一类 Bessel 函数。$\kappa=1$ 时 $J_0=0.765$, $J_1=0.440$, $J_2=0.115$, $J_3=0.020$, $J_4=0.002$... 截断到 n_max=20 绰绰有余。

**第 2 步：平面波的高斯矩阵元。** 需要 $\langle\phi_i|e^{in \cdot 2k_L x_a}|\phi_j\rangle$（高斯夹平面波）。

两个高斯的 overlap 为：

$$M_G = \frac{\pi^{N/2}}{\sqrt{\det K}} \exp\left(C + \frac{1}{4}b^T K^{-1} b\right)$$

其中 K, b, C 是 bra 和 ket 参数合并的结果（PairCache 中已算好）。

**关键技巧**：插入平面波 $e^{in \cdot 2k_L x_a}$ 等价于给线性项 b 的第 a 分量加一个虚部偏移：

$$b \to b' = b + 2ink_L \, \mathbf{e}_a$$

修改后的 overlap $M_G'$ 和原始 $M_G$ 的比值：

$$\frac{\langle\phi_i|e^{in \cdot 2k_L x_a}|\phi_j\rangle}{M_G} = \frac{M_G'}{M_G} = \exp\left(\frac{1}{4}\left[b'^T K^{-1} b' - b^T K^{-1} b\right]\right)$$

**第 3 步：化简。** 设 $\delta b = 2ink_L \mathbf{e}_a$，利用 $\mu = \frac{1}{2}K^{-1}b$（PairCache 中已缓存）：

$$\delta b^T K^{-1} b = 2ink_L \cdot 2\mu_a = 4ink_L\mu_a$$

$$\delta b^T K^{-1} \delta b = -4n^2k_L^2 K^{-1}(a,a)$$

合并得：

$$\frac{\langle\phi_i|e^{in \cdot 2k_L x_a}|\phi_j\rangle}{M_G} = \exp\left(-n^2 k_L^2 K^{-1}(a,a) + in \cdot 2k_L \mu(a)\right)$$

**这是一个解析闭合表达式**，只用到 PairCache 中已有的 $K^{-1}(a,a)$ 和 $\mu(a)$，无需额外积分。

**第 4 步：单粒子 kick kernel。** 把 Bessel 展开和高斯矩阵元合起来：

$$\text{kick}_a = \sum_{n=-n_{\max}}^{n_{\max}} (-i)^n J_n(\kappa) \cdot \exp\left(-n^2 k_L^2 K^{-1}(a,a) + in \cdot 2k_L \mu(a)\right)$$

对应代码中的 `single_particle_kick_kernel()`。高斯的 $e^{-n^2 k_L^2 K^{-1}}$ 因子保证高阶项快速衰减——n 越大，指数衰减越快，Bessel 函数也越小。双重衰减使截断非常有效。

**第 5 步：多粒子 + 置换求和。** N 个粒子的 kick 因子相乘（粒子间独立），加上 Bose 对称化的置换求和：

$$K_{ij} = \sum_p \text{sign}(p) \cdot M_G^{(p)} \cdot \prod_{a=1}^{N} \text{kick}_a^{(p)}$$

其中 p 遍历所有 N! 个置换。

**第 6 步：更新线性系数。** 求 $u'$ 使 $\sum_k u'_k |\phi_k\rangle$ 最优逼近 $U_{\text{kick}}\sum_k u_k |\phi_k\rangle$：

$$u' = S^{-1} K u$$

**注意**：kick 矩阵 K 不是 Hermitian！$U_{\text{kick}}$ 是酉算符，$K_{ji} \neq \overline{K_{ij}}$。必须计算全部 K×K 个矩阵元。

**整个 kick 过程不涉及任何 ODE 求解、C 矩阵、时间积分器**。只是矩阵元计算 + 线性方程组求解。

### 4.3 自由演化

$$
\psi(z, t+T) = e^{-iH_\text{free}T} \cdot \psi(z, t)
$$

可以预先对角化 $H_\text{free}$：

$$
H_\text{free} = U \Lambda U^\dagger
$$

则：

$$
e^{-iH_\text{free}T} = U \cdot e^{-i\Lambda T} \cdot U^\dagger
$$

一次对角化，之后每个 kick 周期只需两次矩阵-向量乘法。

### 4.4 完整 Floquet map

```
每个周期：
    1. ψ → exp(-iκ cos(2k_L z)) · ψ        （kick，逐点乘法）
    2. ψ → exp(-iH_free T) · ψ              （自由演化，矩阵乘法）
    3. 计算 E(n) = <ψ|H_free|ψ>
    4. 计算 n(k) = |FT(ψ)|²
```

### 4.5 N=2 的精确结果

无相互作用时：

$$
E_2(n) = 2 \times E_1(n)
$$

$$
\langle p^2 \rangle_2(n) = 2 \times \langle p^2 \rangle_1(n)
$$

---

## 5. 验证方案

### 5.1 逐 kick 能量对比

最直接的验证：画 $E(n)$ vs $n$（kick 数），ECG+TDVP 的曲线应该和有限差分的曲线重合。

```
E(n)
 ^
 |         ___________  ← 饱和（dynamical localization）
 |        /
 |       /  ← 线性增长（经典扩散）
 |      /
 |     /
 |____/
 |
 └─────────────────→ n (kick 数)
```

**判据：** $|E_\text{ECG}(n) - E_\text{exact}(n)| / E_\text{exact}(n) < 10^{-3}$ 对所有 $n$。

### 5.2 归一化守恒

实时间演化应该保持归一化 $\langle\psi|\psi\rangle = 1$。但 TDVP 是在变分流形上投影，不精确保持归一化。

**判据：** $|\langle\psi|\psi\rangle - 1| < 10^{-4}$ 在整个演化过程中。

如果归一化漂移严重，说明 dt 太大或基底不够。

### 5.3 动量分布 n(k) 冻结

如果 dynamical localization 出现：
- 前几个 kick：n(k) 快速展宽
- 后续 kick：n(k) 不再变化（冻结）

**判据：** $\int |n(k, n) - n(k, n-10)|^2 dk \to 0$ 当 $n$ 足够大。

### 5.4 时间步长收敛性

用 dt = 0.01, 0.005, 0.002 跑相同的物理，如果结果不随 dt 变化，说明时间步长足够小。

**判据：** $|E(n; dt=0.01) - E(n; dt=0.002)| < 10^{-4}$。

### 5.5 基底收敛性

用 K=15, 20, 25 跑相同物理。如果结果不随 K 变化，说明基底足够大。

---

## 6. 预期结果

### 6.1 短时间行为（n < 5-10 kicks）

- 能量近似线性增长
- ECG+TDVP 应该能很好地跟踪（波函数变化不大，基底够用）

### 6.2 中期行为（n ~ 10-30 kicks）

- 能量增长减缓，开始偏离线性
- ECG+TDVP 可能开始出现偏差（波函数在动量空间展宽，需要 B 和 R 参数来编码）
- 这是最关键的测试区间

### 6.3 长时间行为（n > 30 kicks）

- 能量饱和（如果 dynamical localization 出现）
- ECG+TDVP 的挑战：长时间相位累积误差
- 可能需要增加 K 或改进时间积分方案

### 6.4 可能出现的问题

| 问题 | 原因 | 解决方案 |
|------|------|---------|
| 能量爆炸 | dt 太大 | 减小 dt |
| 归一化漂移 | TDVP 投影误差 | 减小 dt，增加 K |
| n(k) 不冻结 | 基底不够大 | 增加 K |
| 与精确解偏差大 | B/R 不够灵活 | 确认 dynamics 模式开启 |
| TDVP 步进太慢 | C 矩阵条件差 | 调整 lambda_C, rcond |

---

## 7. 实现计划

### 7.1 Step 1：编写单粒子有限差分精确解（Python）

```python
# 有限差分 + 分裂算符
# 输入：omega, k_L, kappa, T_period, n_kicks
# 输出：E(n), <p^2>(n), n(k,n) 对每个 kick n
```

这是验证的基准，必须先有。

### 7.2 Step 2：修改 `run_kicked_evolution_test`

- 去掉 delta 相互作用（γ=0）
- 用 B=0,R=0 的 SVM + refine 准备基态
- kick 后切换到 dynamics 模式（B, R 放开）
- 增加 kick 数到 50-100
- 输出格式与 Python 精确解对齐

### 7.3 Step 3：对比验证

- 画 E(n) 曲线：ECG vs exact
- 画 n(k, n) 快照
- 检查归一化守恒
- 做 dt 收敛性测试

### 7.4 Step 4：确认 dynamical localization 信号

- 能量饱和？
- n(k) 冻结？
- 如果 N=2, γ=0 都看不到 → 基底或 TDVP 有问题
- 如果看到了 → 进入下一步：加入相互作用

---

## 8. 物理参数选择

### 8.1 推荐初始参数

| 参数 | 值 | 依据 |
|------|-----|------|
| $\omega$ | 1.0 | 自然单位 |
| $k_L$ | 0.5 | 现有代码默认值 |
| $\kappa$ | 1.0 | 中等 kick 强度，能看到局域化 |
| $T$ | 1.0 | kick 周期 |
| $n_\text{kicks}$ | 50 | 足够看到饱和 |
| dt | 0.01 | 自由演化时间步长（需验证） |
| K | 20 | SVM 基函数数量 |

### 8.2 为什么这些参数

- $\kappa = 1.0$：足够强以看到扩散，但不至于太剧烈导致基底不够
- $T = 1.0$：与阱频率 $\omega = 1$ 不共振（$\omega T / 2\pi \approx 0.159$ 是无理数）
- K=20：静态 kicking benchmark 已验证 K=20 给出 ~3.7e-5 精度

### 8.3 论文中的典型参数

论文 Fig. 2 使用：
- 约 40 次 kick 看到能量饱和
- 不同 $\gamma$ (0 到 11) 对比
- 我们先做 $\gamma = 0$ 的参考线

---

## 9. 成功标准

| 标准 | 阈值 | 说明 |
|------|------|------|
| E(n) 误差 | < 1% per kick | 与有限差分对比 |
| 归一化守恒 | $|1 - \langle\psi|\psi\rangle| < 10^{-3}$ | 整个演化过程 |
| 能量饱和 | 可观察到 | 50 kicks 内 E(n) 趋于常数 |
| dt 收敛 | 结果不随 dt 变化 | dt=0.01 vs dt=0.002 |

**达到这些标准后，我们就可以确信实时间 TDVP + kick 管线是正确的，可以进入下一步：加入相互作用。**

---

## 10. 从 γ=0 到 MBDL 的路线

```
γ=0, N=2  ←  我们在这里
    │
    │  验证通过
    ↓
γ>0, N=2  ←  加入 Gaussian 相互作用（近似 δ 接触）
    │
    │  观察：相互作用如何改变局域化？
    ↓
γ>0, N=3  ←  第一个真正的"多体"
    │
    │  观察：N=3 的 MBDL 信号是否与 N=2 不同？
    ↓
γ>0, N=4  ←  趋势确认
    │
    │  定量分析：MBDL 信号 vs γ, N
    ↓
论文结果  ←  与实验对比
```

---

## 11. 实验记录

### 11.1 实验 A：TDVP kick（失败）

**方法**：用 TDVP 做一步实时间演化来模拟 kick（原始 `apply_kick` 实现）。

**参数**：K=10, dt=0.05, 20 kicks

**结果**：

| Kick | ECG E | 精确 E |
|------|-------|--------|
| 0 | 1.000 | 1.000 |
| 1 | 1.001 | **1.315** |
| 2 | 1.001 | **1.705** |
| 5 | 1.001 | **1.073** |

**分析**：ECG 能量几乎不变（0.1% 变化），而精确解第一次 kick 后涨了 30%。

**原因**：TDVP kick 在基底可表达的范围内做"最优近似"——但 $e^{-i\kappa\cos(x)}|\psi\rangle$ 不在 Gaussian 基底内，TDVP 无法表达 kick 注入的高动量成分，被困在初态附近。

**结论**：δ-kick 是瞬时的，不能用时间演化来模拟。

### 11.2 实验 B：解析 kick（进行中）

**方法**：Jacobi-Anger 展开计算 kick 矩阵元 $K_{ij} = \langle\phi_i|e^{-i\kappa V}|\phi_j\rangle$，通过 $u' = S^{-1}Ku$ 精确更新线性系数。kick 之间用 TDVP 做自由演化。

**数学**：

$$e^{-i\kappa\cos\theta} = \sum_{n=-\infty}^{\infty} (-i)^n J_n(\kappa) e^{in\theta}$$

对高斯基函数：

$$\frac{\langle\phi_i|e^{in \cdot 2k_L x_a}|\phi_j\rangle}{M_G} = e^{-n^2 k_L^2 K^{-1}(a,a) + in \cdot 2k_L \mu(a)}$$

N 个粒子的 kick kernel 相乘（独立坐标）。截断到 n_max=20 项（$J_n(\kappa=1)$ 在 n>5 时已极小）。

**关键实现细节**：

- kick 矩阵不是 Hermitian（$U = e^{-i\kappa V}$ 是酉算符，不是 Hermitian），必须计算全部 K×K 矩阵元
- 首次实现误用了 `K(j,i) = conj(K(i,j))`（Hermitian 假设），导致 E → -427，已修复

**参数**：K=10, dt=0.05, 20 kicks

**结果（kick-only，无自由演化）**：

```
Kick  0: E=1.000,  S=6.283  (初态)
Kick  1: E=1.399,  S=6.184  (精确值 1.315，误差 6%)
Kick  2: E=2.569,  S=5.871  (精确值 1.705，已严重偏离)
Kick  5: E=8.721,  S=4.225  (精确值 1.073)
Kick 10: E=13.81,  S=3.322
Kick 20: E=11.54,  S=44.26  (norm 发散)
```

第一次 kick 方向正确（能量涨了 40%），但没有自由演化时误差纯累积。Norm 从 6.28 发散到 44.26。

**关键诊断输出解读**：

```
[kick] norm before=6.28, norm after=6.18, ratio=0.984, |u|=1 -> 2.59, S cond=2577
[post-kick] E=1.399, S=6.28
```

| 参数 | 含义 | 当前值 | 评价 |
|------|------|--------|------|
| norm before | kick 前 u†Su（波函数模方） | 6.283 | 正常（= 2π，归一化约定） |
| norm after | kick 后 u'†Su'（投影后模方） | 6.184 | 掉了 1.6%，因为 kicked 态不完全在基底空间内 |
| ratio | norm_after/norm_before | 0.984 | 每次 kick 损失 ~2% 的范数 |
| \|u\| | u 向量的欧几里得范数 | 1→2.59 | kick 注入了大的复数分量 |
| S cond | overlap 矩阵 S 的条件数 | 2577 | 适中，不是问题 |
| post-kick E | kick 后立刻测量的能量 | 1.399 | 精确值 1.315，误差 6%（K=10 基底不够） |
| post-kick S | 重归一化后的 norm | 6.283 | 已恢复（加了 renormalize 后） |

**结论**：
- kick 矩阵本身基本正确（第一次 kick 能量变化方向和量级对）
- 6% 误差来自 K=10 基底的投影损失（kicked 态有高动量成分，Gaussian 基底表达不了）
- 需要更大的 K 和/或自由演化来让基底"适应" kick 后的波函数

### 11.3 实验 C：解析 kick + TDVP 自由演化（进行中）

**改进**：
1. Kick 后重归一化 u（防止 TDVP 梯度爆炸）
2. 减小 dt 从 0.05 到 0.01（100 步/周期）

**参数**：K=10, dt=0.01, 5 kicks, static_mode (params=40, 只有 A)

**结果**：

```
Kick  0: t=0, E=1.000, S=6.283
  [kick] norm=6.283→6.184, |u|=1→2.59
  [post-kick] E=1.399, S=6.283 (renormalized)
Kick  1: t=1, E=1.847, S=27.326  ← norm 从 6.28 飙到 27.3（自由演化导致）
  [kick] norm=27.3→NaN  ← 第二次 kick 就 NaN
```

**分析**：free_evolve 不保持 norm。TDVP 更新 A 参数后 S 矩阵变了，但 u 没有同步调整，导致 u†Su 从 6.28 飙到 27.3（4 倍）。到第二次 kick 时 S 矩阵已经病态，LDLT 分解失败。

**修复**：在 `realtime_tdvp_step` 末尾加归一化——每步更新 A 后重算 S，调整 u 使 u†Su 保持不变。代价是每步多两次 `build_HS`（速度减半），但保证了 norm 守恒。

### 11.4 实验 D：解析 kick + 归一化 TDVP（进行中）

**改进**：
1. kick 后 renormalize u
2. TDVP 每步后 renormalize u（保持 u†Su 不变）
3. static_mode (params=40)

**参数**：K=10, dt=0.01, 5 kicks

**结果**：运行中...

### 11.5 TDVP 崩溃的根本原因：前向 Euler 对振荡问题不稳定

#### 实时间 vs 虚时间

TDVP 方程：
- 虚时间（能量最小化）：C dz/dt = -g → 对应 dx/dt = -ωx（衰减）
- 实时间（薛定谔演化）：C dz/dt = -ig → 对应 dx/dt = -iωx（振荡）

前向 Euler 更新：z(t+dt) = z(t) + dz * dt

对 dx/dt = -iωx，Euler 给 x(t+dt) = x(t)(1 - iω*dt)，幅度为：

    |1 - iω*dt|² = 1 + ω²dt² > 1  （每步都在增长！）

对 dx/dt = -ωx，Euler 给 x(t+dt) = x(t)(1 - ω*dt)，幅度为：

    |1 - ω*dt| < 1  （当 dt < 2/ω，衰减，稳定）

**虚时间是衰减 → Euler 稳定。实时间是振荡 → Euler 永远不稳定。**

#### 具体数值估算

Kick 后波函数在基底中有效频率 ω ~ E_max（基底中最高本征值）。
对 K=10，E_max ~ 50。dt = 0.01 时：

    每步增长 = √(1 + 50² × 0.01²) = √1.25 ≈ 1.118
    100 步（一个周期）后 = 1.118^100 ≈ 10^5

**Norm 增长 10 万倍。** 这就是为什么一个周期后 norm 从 6.28 飙到 27 然后 NaN。

减小 dt 没用——Euler 对振荡问题**无条件不稳定**。dt 再小也只是减缓增长速度，不能消除。

#### 为什么虚时间 TDVP 一直正常

虚时间 TDVP (C dz = -g) 对应衰减型 ODE。前向 Euler 对衰减问题条件稳定（dt < 2/ω_max）。我们用的 dt ~ 0.01-1.0 足够小，所以虚时间从来没崩过。

#### 修复方案：辛积分器

用 **RK2 midpoint 方法**替代 Euler：

    1. 半步: z_half = z(t) + dt/2 * f(z(t))
    2. 在 z_half 处重算 C 和 g
    3. 全步: z(t+dt) = z(t) + dt * f(z_half)

Midpoint 法对振荡问题的幅度变化：

    |1 - iωdt + (-iωdt)²/2| ≈ 1 + O(ω⁴dt⁴)

dt=0.01, ω=50 时：每步增长 ~ 1 + 3e-4。100 步后增长 ~3%，而非 Euler 的 10^5 倍。

### 11.6 实验 E：固定基底 + 精确传播（K=10 和 K=30）

放弃 TDVP 自由演化，改用固定基底精确传播：u(t) = exp(-i S^{-1}H t) u(0)

基底 A/B/R 不变，只有 u 随时间演化。通过特征分解一次预计算，之后每个周期只需矩阵乘向量。

#### K=10 结果

| Kick | ECG E | 精确 E | 误差 |
|------|-------|--------|------|
| 1 | 1.399 | 1.315 | +6% |
| 3 | 1.333 | 1.636 | -19% |
| 5 | 2.049 | 1.073 | +91% |
| 10 | 6.055 | 1.094 | +453% |
| 20 | 11.71 | 1.268 | +823% |

能量单调增长而精确解振荡。Norm 完美守恒。

#### K=30 结果

| Kick | ECG E | 精确 E | 误差 |
|------|-------|--------|------|
| 1 | 1.404 | 1.315 | +7% |
| 3 | **1.651** | **1.636** | **+0.9%** |
| 5 | 1.439 | 1.073 | +34% |
| 10 | 2.305 | 1.094 | +111% |
| 20 | 13.04 | 1.268 | +928% |

前 3-4 个 kick 明显改善（K=30 误差 <10%），但之后仍然发散。

#### 结论

固定基底方案的根本限制：基函数为基态优化（B=0, R=0），不携带动量信息。Kick 注入的高动量成分每次都被投影损失，误差逐 kick 累积。增大 K 只是推迟发散，不能消除。

**必须让基函数能演化 B 和 R（携带动量），即回到 TDVP 路线——但用辛积分器替代前向 Euler。**

### 11.6.1 误差来源分析：为什么自由演化误差为零，问题全在 kick

#### 自由演化是精确的

固定基底下，薛定谔方程变成常系数线性 ODE：

    i S du/dt = H u

S 和 H 都是常数矩阵（基底不动），精确解是：

    u(t) = Σ_n c_n exp(-i E_n t) v_n

其中 (E_n, v_n) 是广义本征对 Hv = ESv，c_n 是初态的展开系数。

**每个本征分量只乘了一个相位 exp(-iE_n t)**。这是解析精确的，不涉及任何数值近似（不是 Euler 也不是 RK4，是精确解）。

能量严格守恒：

    E(t) = Σ_n |c_n|² E_n / Σ_n |c_n|²

因为 |c_n exp(-iE_n t)|² = |c_n|²，相位不改变模。分子分母都不随时间变化。

**不是"误差很小"，是数学上精确为零。**

唯一的近似在于我们用了 K 个基函数而不是无穷个——但在这 K 维子空间内部，演化没有任何误差。

类比：真实波函数是一首完整的交响乐，我们只保留了 K 个音轨。在这 K 个音轨内演奏是完美的，但一开始就缺了一些音轨——这个缺失不是演化造成的。

#### Kick 是全部误差的来源

Kick 算符 exp(-iκcos(2k_L x)) 在位置空间逐点乘相位。Bessel 展开包含所有频率：

    e^{-iκcos(θ)} = Σ_{n=-∞}^{∞} (-i)^n J_n(κ) e^{inθ}

乘完之后，波函数有了基底中**不存在**的高动量分量（n × 2k_L）。

投影 u' = S^{-1}Ku 是把无穷维的 kicked 态压缩到 K 维子空间——基底外的分量被永久丢弃：

```
|ψ'⟩ = e^{-iκV} |ψ⟩          ← 无穷维，包含所有动量分量
    ↓ 投影到 K 维基底
|ψ'_proj⟩ = Σ u'_k |φ_k⟩     ← 丢失了基底外的高动量分量
    ↓ 自由演化（精确，零损失）
|ψ(T)⟩                         ← K 维内精确，但信息已经少了
    ↓ 下一次 kick 投影
|ψ'_proj⟩                      ← 又丢了一些
    ...
```

#### 误差累积是单向的

```
kick 1: 投影丢失 ~2% → 剩余 98% 的信息
   ↓ 自由演化（精确，零损失，98% 不变）
kick 2: 再投影丢失 ~2% → 剩余 96%
   ↓ 自由演化（精确，零损失）
kick 3: 再丢失 ~2% → 剩余 94%
   ...
kick 50: 剩余 ~36%
```

每次 kick 都在缩小波函数在基底中的"有效信息量"，自由演化不损失也不恢复。这是**不可逆的投影损失**。

#### 为什么精确解（有限差分）没有这个问题

有限差分用 N_grid = 512 个网格点，相当于 K=512 的"基底"（位置空间 delta 函数）。这个基底包含所有动量成分直到 Nyquist 频率 π/dz。Kick 的 cos(2k_L x) 注入的最高动量 ~ 2κ × k_L，远低于 Nyquist 频率，投影损失接近零。

我们的 K=10~30 个 Gaussian 基函数远不够覆盖 kick 需要的动量空间。

#### 总结

> **自由演化 = K 维子空间内的精确旋转（零误差）。**
>
> **Kick = 从无穷维到 K 维的投影（不可避免的信息损失）。**
>
> **改善精度的唯一途径：让基底能更好地表达 kicked 态（更大的 K，或包含动量分量的基函数）。自由演化不需要改。**

### 11.7 实验 F：RK2 midpoint + dynamics mode（崩溃）

**改进**：用 RK2 midpoint 替代前向 Euler

**参数**：K=10, dt=0.02, 5 kicks, dynamics mode (params=80)

**结果**：

```
Kick 0: E=1.000, S=6.283
[kick] → E=1.399
Kick 1: E=5.738, S=51.08  ← norm 从 6.28 飙到 51
[kick] → NaN
```

RK2 没有解决问题。Norm 仍然在一个周期内增长 8 倍。

**深层原因**：问题不在时间积分器，而在 TDVP 的参数更新逻辑。

当前实现只更新 z = (A, B, R)，**不更新 u**（线性系数）。当 A/B/R 变了：
- 基函数 φ_k(z) 变了
- 但 u 还是旧的
- ψ = Σ u_k φ_k(z_new) 已经不是之前那个物理态
- u†S(z_new)u ≠ u†S(z_old)u → norm 漂移

虚时间没这个问题，因为每步后 `resolve_u=true` 重新解本征值。实时间不能 resolve_u（态不是本征态），所以 u 和 z 不匹配。

### 11.8 实验 G：RK2 + u 参与 TDVP 演化（进行中）

**修复**：将 `updata_constant` 设为 0，让 u 和 z 一起参与 TDVP 演化。

原理：完整的 McLachlan TDVP 应该同时优化所有参数（u 和 z）。之前排除 u 是因为虚时间用 resolve_u 处理 u。实时间应该让 u 也作为变分参数参与演化。

C 矩阵变成 (K+d_z) × (K+d_z) 维，包含 u-u、u-z、z-z 块。dz 向量同时包含 du 和 d(A,B,R)。

**参数**：K=10, dt=0.02, 5 kicks, dynamics mode (params=80), RK2 midpoint

**结果**：

```
Kick 0: E=1.000, S=6.283
[kick] → E=1.399
Kick 1: E=73.01, S=-0.551  ← norm 变负，崩溃
```

仍然不稳定。u 参与 TDVP 后问题变成：u（良态 ODE）和 z（病态 ODE）混在同一个 C 矩阵中，良态被病态污染。

---

## 12. 文献调研：为什么 TDVP 不行，以及该怎么办

### 12.1 核心发现：学界也没有让 TDVP + ECG 实时间演化成功

> "It is well-known that the TDVP approach suffers from severe numerical instabilities when applied to linear combinations of Gaussians with adjustable parameters."
> — Schrader, Kristiansen, Kvaal, Pedersen, JCP 162, 024109 (2025)

**Varga 书（1998）完全不涉及实时间动力学。** 全书只讨论束缚态（静态问题）。

**唯一成功的 ECG 实时间动力学**来自 Adamowicz/Kvaal 组（2024-2025），他们用的是 **Rothe 方法**——完全不用 TDVP。

### 12.2 TDVP 崩溃的根本原因

| 层次 | 问题 | 影响 |
|------|------|------|
| 数学 | C 矩阵（度量张量）对近线性相关的 Gaussian 严重病态 | C^{-1} 放大数值噪声 |
| 物理 | Kick 后的态不在变分流形上 | 梯度 g 巨大，TDVP 线性化假设失效 |
| 数值 | 前向 Euler/RK2 对振荡 ODE 不保辛 | 幅度每步增长 |
| 参数 | A 更新后可能丢失正定性 | 基函数变为不可归一化 |
| 耦合 | u（良态）和 z（病态）混合求解 | 良态被病态污染 |

### 12.3 文献中的解决方案

| 方案 | 来源 | 核心思路 | 难度 |
|------|------|---------|------|
| **固定基底 + 精确传播** | 本项目 | 不动 A/B/R，只演化 u | 已实现 |
| **动量增广基底** | 标准量子化学 | 预建包含高动量成分的基底 | 低 |
| **分步法（u 精确 + z 虚时间）** | Lubich (2015), MCTDH 社区 | 分离 u 和 z 的演化 | 中 |
| **基底增广** | AIMS (Ben-Nun & Martinez) | kick 后动态添加新基函数 | 中 |
| **Rothe 方法** | Schrader et al. (2025) | 每步重新优化所有参数 | 高 |
| **投影分裂积分器** | Kloss, Lubich (2017) | 避免 C^{-1}，保辛 | 高 |

### 12.4 推荐实施路线

**阶段 1（立即）**：动量增广基底 + 固定基底精确传播

在 SVM 基底中添加携带动量的 Gaussian（非零 B 和 R），使基底能表达 kick 后的高动量成分。然后用已有的 `free_evolve_fixed_basis` + `apply_analytic_kick`。

Kick 注入动量 p = n × 2k_L。基函数需要 Im(2R_a B_aa) = p 来编码动量。具体地，对每个基态基函数，生成 n = ±1, ±2 的动量副本：
- B_aa = b（某个正值）
- R_a = i × n × k_L（纯虚，编码动量）

**阶段 2（如果阶段 1 不够）**：分步法

```
每个时间子步 dt:
    1. u 步：固定 A/B/R，用特征分解精确传播 u（已有代码）
    2. z 步：固定 u，做一步虚时间 TDVP 调整 A/B/R（已有代码，稳定）
    3. 重投影 u：在新 A/B/R 下重解本征值，投影到最近的态
```

u 的实时间演化由精确求解器完成（稳定），z 的优化由虚时间 TDVP 完成（也稳定）。两者互不干扰。

**阶段 3（N>2 扩展时）**：考虑 Rothe 方法

### 12.6 动量增广基底实验结果

| 配置 | K_total | kick 1 误差 | kick 10 | kick 50 | 评价 |
|------|---------|------------|---------|---------|------|
| 纯位置 K=10 | 10 | 6% | 453% | 823% | 单调发散 |
| 纯位置 K=30 | 30 | 7% | 111% | 928% | 推迟但仍发散 |
| K=5+mom(n=2,b=0.3) | 17 | **5.5%** | 130% | **饱和 ~4** | **最佳** |
| K=8+mom(n=4,b=0.3) | 14 | 7.5% | 131% | 饱和 ~4.2 | 条件数阈值限制 |
| K=8+mom(n=4,b=1.0) | 19 | 7.5% | 378% | 9.1 | b_val 太大 |
| K=15+mom(n=3,b=0.3) | 29 | 7.1% | 702% | 8.8 | 冗余太多 |

**结论**：
1. 动量增广显著改善了前 3-5 kicks（误差从 >100% 降到 <30%）
2. 但长时间仍然漂移（50 kicks 后 ECG ~4-5 vs 精确 ~1.4）
3. 最优配置是**少量高质量基函数**（K=17），不是堆砌数量
4. S 条件数是核心瓶颈——动量基函数和位置基函数之间天然高度重叠
5. b_val=0.3 最优，太大（1.0）反而更差

**下一步**：固定基底方案已接近天花板。要进一步改善需要：
- 分步法（阶段 2）：允许 A/B/R 适应性调整
- 或 Rothe 方法（阶段 3）：彻底重新优化每个时间步的基底

### 12.5 Rothe 方法详解（Schrader et al., JCP 2025）

#### 核心思路：时间演化 → 优化问题

TDVP 把时间演化变成 ODE（$C\dot{z} = -ig$），然后用数值积分器解这个 ODE。问题是 C 矩阵病态，积分器不稳定。

Rothe 方法完全不同：**每个时间步都是一个独立的优化问题**，不解 ODE，不需要 C 矩阵。

#### 第 1 步：Crank-Nicolson 传播子

精确的时间演化 $\Psi(t+\Delta t) = e^{-iH\Delta t}\Psi(t)$ 用 Crank-Nicolson 近似：

$$\Psi(t_{i+1}) = \hat{A}_i^{-1} \hat{A}_i^\dagger \Psi(t_i)$$

其中

$$\hat{A}_i = \hat{I} + i\frac{\Delta t}{2}\hat{H}\left(t_i + \frac{\Delta t}{2}\right)$$

Crank-Nicolson 是二阶精度、保辛、无条件稳定的隐式格式。

#### 第 2 步：定义 Rothe 误差

不直接算 $\hat{A}_i^{-1}$（在无穷维空间中做不到），而是两边乘以 $\hat{A}_i$：

$$\hat{A}_i \Psi(t_{i+1}) = \hat{A}_i^\dagger \Psi(t_i)$$

定义 Rothe 误差为：

$$r^{i+1}(\alpha, c) = \left\| \sum_m c_m \hat{A}_i \phi_m(\alpha) - \hat{A}_i^\dagger \Psi(t_i) \right\|^2$$

其中 $\alpha$ 是非线性参数（A/B/R），$c$ 是线性系数。

**关键**：右端 $\hat{A}_i^\dagger\Psi(t_i)$ 是已知的（上一步的波函数乘以 $\hat{A}_i^\dagger$），左端是待优化的下一步波函数乘以 $\hat{A}_i$。最小化这个误差就是在找最接近精确 Crank-Nicolson 传播的参数化波函数。

#### 第 3 步：VarPro 分离线性和非线性参数

对给定的非线性参数 $\alpha$，最优线性系数 $c$ 有闭合解：

$$c^{i+1}(\alpha) = [S^{i+1}(\alpha)]^{-1} \rho^{i+1}(\alpha)$$

其中：

$$S^{i+1}_{mn}(\alpha) = \langle\hat{A}_i\phi_m(\alpha)|\hat{A}_i\phi_n(\alpha)\rangle$$

$$\rho^{i+1}_m(\alpha) = \langle\hat{A}_i\phi_m(\alpha)|\hat{A}_i^\dagger\Psi(t_i)\rangle$$

注意 $S^{i+1}$ 不是普通的 overlap 矩阵！它是 **"修正 overlap"**，包含了 $\hat{A}_i = I + i\frac{\Delta t}{2}H$ 的效应。

展开后：

$$\langle\hat{A}_i\phi_m|\hat{A}_i\phi_n\rangle = \langle\phi_m|\phi_n\rangle + i\frac{\Delta t}{2}[\langle\phi_m|H|\phi_n\rangle - \langle H\phi_m|\phi_n\rangle] + \frac{\Delta t^2}{4}\langle\phi_m|H^2|\phi_n\rangle$$

**需要 $\langle\phi_m|H^2|\phi_n\rangle$ 矩阵元！** 这是 Rothe 方法比 TDVP 多需要的东西。

#### 第 4 步：优化非线性参数

把 $c(\alpha)$ 代入后，Rothe 误差只是 $\alpha$ 的函数 $r^{i+1}(\alpha)$。用 BFGS 算法优化 $\alpha$。

初始猜测：$\alpha^{i+1}_{init} = \alpha^i_{opt} + \delta(\alpha^i_{opt} - \alpha^{i-1}_{opt})$（线性外推，$\delta$ 通过线搜索确定）。

#### 为什么 Rothe 比 TDVP 稳定

| | TDVP | Rothe |
|---|---|---|
| 每步做什么 | 解 ODE：$C\dot{z} = -ig$ | 解优化：$\min_\alpha r(\alpha)$ |
| 需要 $C^{-1}$ 吗 | **需要**（病态！C^{-1} 放大噪声） | **不需要** |
| 需要 $H^2$ 矩阵元吗 | 不需要 | **需要** |
| 时间积分器 | Euler/RK（外部积分器） | **Crank-Nicolson（内嵌在目标函数中）** |
| 稳定性 | 差（C 病态 + 积分器不保辛） | **好（优化问题天然稳定）** |
| 基底自适应 | 能但不稳定 | **能且稳定**（自然随优化调整） |
| norm/能量守恒 | 不保证 | 可通过约束优化强制保证 |

核心区别：TDVP 需要求 $C^{-1}$（度量张量的逆），这对近线性相关的 Gaussian 严重病态。Rothe 方法把 $C^{-1}$ 的角色转移给了 $S^{-1}$（修正 overlap 的逆），而 $S$ 可以加正则化（$S + \lambda I$）而不影响物理。

#### 他们的结果

| 系统 | 维度 | Gaussian 数 | 传播时间 | 精度 |
|------|------|-----------|---------|------|
| Henon-Heiles | 2D | **20** | t=100 | 网格级 |
| Henon-Heiles | 3D | **30** | t=100 | 定量吻合 |
| Henon-Heiles | 4D | **40** | t=100 | 近定量 |

dt=0.01，10000 步。代码用 Python + mpi4py 并行化矩阵元计算。

**20 个 Gaussian 在 2D 就达到网格级精度**——说明 ECG 基底的表达力是够的，问题纯粹在演化方法（TDVP 不行，Rothe 行）。

#### 实现 Rothe 方法需要什么

**已有（可复用）**：

| 组件 | 文件 | 状态 |
|------|------|------|
| $\langle\phi_i\|\phi_j\rangle$ overlap | `pair_cache.cpp` | 已有 |
| $\langle\phi_i\|H\|\phi_j\rangle$ Hamiltonian | `interaction_kernels.cpp` | 已有 |
| 置换求和框架 | `permutation.cpp` | 已有 |
| 基函数参数化 | `basis_params.cpp` | 已有 |
| 解析 kick | `kick_operator.cpp` | 已有 |

**需要新增**：

| 组件 | 工作量 | 说明 |
|------|--------|------|
| $\langle\phi_i\|H^2\|\phi_j\rangle$ 矩阵元 | **大** | 需要推导 $T^2$, $TV$, $V^2$ 等 kernel |
| BFGS 优化器 | 中 | 可用 Eigen unsupported 或手写 L-BFGS |
| VarPro 框架 | 中 | 分离 $c$ 和 $\alpha$ 的优化 |
| Rothe 误差函数及梯度 | 中 | $\partial r / \partial\alpha$ 需要矩阵元对参数的导数（已有部分） |
| norm/能量约束 | 小 | 约束优化或事后投影 |

最大工作量在 $H^2$ 矩阵元。但论文指出：对多项式势，$H^2$ 的矩阵元可以用 Isserli 定理（Wick 定理）从**高阶矩**计算。我们的 Hamiltonian 是动能 + 谐振子 + cos(x)，其中 cos(x) 可以用 Bessel/指数展开处理。

### 12.6 关键参考文献

| 主题 | 文献 |
|------|------|
| **Rothe 方法 + ECG（核心参考）** | Schrader et al., JCP 162, 024109 (2025) |
| Rothe 方法理论基础 | Kvaal, JCP 2023; Schrader et al., JCP 2024 |
| 投影分裂积分器 | Lubich, AMRX 2015; Kloss et al., JCP 2017 |
| vMCG 正则化 | Richings et al., IRPC 2015 |
| 保辛几何积分器 | Begusic & Vanicek, arXiv:2310.05633 (2023) |
| 约束变分法 | Shalashilin, JPCA 2013 |
| DF vs McLachlan VP | Raab, CPL 2000 |
