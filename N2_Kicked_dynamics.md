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

在位置空间离散化：$z_j = z_\min + j \cdot dz$, $j = 0, 1, ..., N_\text{grid}-1$。

$$
H_\text{free} = -\frac{1}{2}\frac{d^2}{dz^2} + \frac{1}{2}\omega^2 z^2
$$

用中心差分近似 $d^2/dz^2$，得到 $N_\text{grid} \times N_\text{grid}$ 的三对角矩阵。

### 4.2 Kick 的施加

理想 delta kick 等价于瞬时乘以相位：

$$
\psi(z, t^+) = e^{-i\kappa\cos(2k_L z)} \cdot \psi(z, t^-)
$$

这在位置空间是**逐点乘法**，不涉及矩阵运算。

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
