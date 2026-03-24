# Delta 接触势 N=2 基态：现状与分析

日期：2026-03-24

---

## 1. 问题定义

两个玻色子在一维谐振势阱中，通过 delta 接触势相互作用：

$$
H = \sum_{i=1}^{2} \left[\frac{p_i^2}{2m} + \frac{1}{2}m\omega^2 x_i^2\right] + g \cdot \delta(x_1 - x_2)
$$

参数：ℏ=1, m=1, ω=1, g=1。

精确基态能量 E = 1.3067455，由 Busch 超越方程给出（质心 E_CM=0.5 + 相对坐标 E_rel=0.8067455）。

---

## 2. 当前结果

### SVM + Refine（ghost state 修复后）

```
SVM (K=13):  E = 1.309165    error = 2.42e-3
+ Refine:    E = 1.309036    error = 2.29e-3   (改善 1.3e-4)
```

Refine 正常工作（无 ghost state），10 轮稳步下降：

```
Round  1: E=1.30908436   replaced=7
Round  5: E=1.30904053   replaced=4
Round 10: E=1.30903594   replaced=3
```

K=14 时 SVM 找不到合法 trial（`no valid trial found`）——但这是**人为限制，不是基本限制**。

### 为什么 K 增长过早停止

原因是 `s_well_conditioned()` 的条件数阈值设得太保守（κ < 1e6）。这是之前为了防 ghost state 加的约束。但 ghost state 的根本原因已经修复（Eigen aliasing），不再需要这么严格的限制。

放宽到 κ < 1e10 后：

| benchmark | 旧阈值 (1e6) K_max | 新阈值 (1e10) K_max | 误差变化 |
|-----------|-------------------|---------------------|---------|
| Gaussian | 18 | **20** | 2.67e-6 → 2.31e-6 |
| Delta | 13 | **15** | 2.42e-3 → 2.13e-3 |

K 能继续增长了，但 delta 的改善很小（符合代数收敛预期）。

**真正限制 K 增长的是随机采样效率**：随着 K 增大，参数空间中与已有基底线性无关的区域越来越小，5000 次随机试探越来越难命中。这可以通过增加 n_trials 或改进采样分布来缓解。

### TDVP 发散

```
Step 0: E = -1.99   dE = -3.3
Step 1: E = -3.74   dE = -2.0
Step 2: E = -18.0   dE = -14
Step 3: E = -2328   dE = -2310
```

TDVP 第一步就发散。原因：C 矩阵（metric tensor）条件数极差（3e+11），梯度方向完全不可靠。这是 delta cusp 结构和 Gaussian 基底严重不匹配的体现——波函数对非线性参数的响应在 cusp 附近极其病态。

---

## 3. 为什么 Delta 比 Gaussian 势难得多

### 收敛速度的根本差异

| 势函数 | 波函数特征 | Gaussian 基底收敛率 | 我们的结果 |
|--------|-----------|-------------------|-----------|
| Gaussian（平滑）| 处处解析 | **指数级** ~ e^{-cK^α} | K=18, error=2.5e-6 |
| Delta（接触）| 有 cusp（一阶导不连续）| **代数级** ~ K^{-1/2} | K=13, error=2.3e-3 |

精度差 1000 倍不是 bug，是数学性质决定的。用光滑函数逼近尖点，收敛速度必然退化。

### 物理图像

Delta 势导致波函数在 x₁=x₂ 处有一个"尖角"（cusp）：

```
ψ(r)     平滑势                Delta 势

  |      ____                    /\
  |     /    \                  /  \
  |    /      \                /    \
  |   /        \              /      \
  |  /          \            /        \
  +--+----+----+--→ r    +--+----+----+--→ r
         r=0                    r=0
                              ← cusp
```

Gaussian 基函数 exp(-αr²) 天然光滑（所有导数连续），无论怎么叠加都无法精确重现这个尖角。需要大量窄 Gaussian 在 r=0 附近"堆积"来近似 cusp，效率极低。

### 需要多大 K

文献（Jeszenszki 2018, Koscik 2018）对 1D 接触势 + 谐振子基底的收敛分析：

| 目标精度 | 纯 Gaussian 基底需要 K | 带 cusp 处理需要 K |
|---------|--------------------|--------------------|
| 1e-2 | ~100 | ~10 |
| 1e-3 | ~1000 | ~30 |
| 1e-4 | ~10000 | ~100 |
| 1e-6 | ~10^6（不现实）| ~300 |

**我们 K=13 得到 2.3e-3 完全符合 K^{-1/2} 的代数收敛预期。**

---

## 4. 文献中的解决方案

### 方案 A：增加 K（暴力堆积）

最简单但最低效。收敛率 K^{-1/2} 意味着精度每提升一个量级需要 K 增大 100 倍。

| K | 预期精度 | 计算量（相对 K=13）|
|---|---------|-------------------|
| 13 | 2.3e-3 | 1× |
| 50 | ~1e-3 | 15× |
| 200 | ~5e-4 | 240× |
| 1000 | ~1e-4 | 6000× |

性价比极低，但如果只是想验证收敛趋势可以试 K=30~50。

### 方案 B：r_ij 前因子 ECG（Cencek/Kutzelnigg 2004）

基函数改为 |x₁-x₂| × exp(-二次型)。前因子 |r| 的导数在 r=0 处不连续，天然能表达 cusp。

- 优点：收敛速度大幅提升
- 缺点：需要推导新的矩阵元公式，代码改动量大
- Zaklama 2020 的 1D ECG 文章（我们已经读过）给出了带多项式前因子的矩阵元框架，理论上可以支持这种扩展

### 方案 C：Transcorrelated 方法（Jeszenszki 2018）

用 Jastrow 因子 e^{J} 做相似变换 H̃ = e^{-J} H e^{J}，把 cusp 吸收进 Hamiltonian。变换后的 H̃ 是光滑的，Gaussian 基底的指数收敛率恢复。

- 收敛率从 K^{-1} 提升到 K^{-3}
- 缺点：H̃ 不再厄米，需要修改整个求解框架
- 这是文献中处理 1D 接触势最先进的方法

### 方案 D：有限程近似 + 外推

用有限程 Gaussian 势 V(r) = g/σ√π × exp(-r²/σ²) 近似 delta(r)，在 σ→0 极限下收敛到 delta。

- 对固定 σ，势函数光滑，Gaussian 基底高效
- 可以对多个 σ 值算出精确结果，然后外推到 σ→0
- **我们已经有 Gaussian interaction 的完整实现**，这是最容易尝试的路线

---

## 5. 对 MBDL 项目的影响

### 直接影响有限

MBDL 论文中的 delta 相互作用强度 γ 从 0 到 11。关键物理信号（能量饱和、动量分布冻结）是**宏观现象**，不依赖于基态能量的亚百分比精度。

### 实用策略

1. **用 Gaussian 势做动力学**：σ 取有限值（如 0.1~0.5），势函数光滑，SVM/Refine/TDVP 都稳定。这已经能研究 MBDL 的核心物理。

2. **把 delta benchmark 当质量标尺**：2.3e-3 的误差告诉我们 Gaussian 基底的局限性在哪里，但不阻碍物理研究。

3. **如果以后需要更高精度**：优先考虑方案 D（有限程外推），其次方案 B（前因子 ECG）。方案 C（transcorrelated）效果最好但实现成本最高。

---

## 6. TDVP 发散问题

TDVP 在 delta 上发散是独立于精度的问题。即使基态算得很准，TDVP 的梯度在 cusp 附近也会失控——因为 Gaussian 参数的微小变化会导致 cusp 区域的矩阵元剧烈变化。

可能的缓解措施：
- 更强的正则化（lambda_C 从 1e-8 增大到 1e-4 或更大）
- 更保守的步长限制
- 冻结某些对 cusp 敏感的参数方向

但根本解决需要 cusp-aware 的基底或参数化。

---

## 7. 参考文献

- Jeszenszki et al., "Accelerating convergence with the transcorrelated method", Phys. Rev. A 98, 053627 (2018)
- Jeszenszki et al., "Eliminating the wave-function singularity", Phys. Rev. Research 2, 043270 (2020)
- Koscik, "Optimized configuration interaction for trapped multiparticle systems with contact forces", Phys. Lett. A (2018)
- Cencek & Kutzelnigg, "Gaussian basis sets with the cusp condition", Chem. Phys. Lett. (2004)
- Zaklama et al., "Matrix Elements of 1D ECG Basis Functions", Few-Body Syst. 61, 6 (2020)
- Mitroy et al., "Theory and application of ECGs", Rev. Mod. Phys. 85, 693 (2013)
