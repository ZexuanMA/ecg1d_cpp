# Varga-Suzuki SVM 方法学习笔记

参考书：Suzuki & Varga, *Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems* (Springer, Lecture Notes in Physics Monographs Vol. 54, 1998)

综述文章：Mitroy, Bubin, Horiuchi, Suzuki, Adamowicz, Cencek, Szalewicz, Komasa, Blume, Varga, "Theory and application of explicitly correlated Gaussians", *Rev. Mod. Phys.* **85**, 693 (2013)

---

## 1. 核心思想

SVM (Stochastic Variational Method) 的核心是：

- 用 **显式关联高斯函数 (ECG, Explicitly Correlated Gaussians)** 作为变分基底
- 用 **随机试探 + 选择** 的方式优化非线性参数
- 线性系数通过求解广义本征值问题自动确定

这套方法的优势在于：

1. 矩阵元全部解析可算
2. 基底紧凑（K 很小就能高精度）
3. 不需要事先知道波函数结构
4. 对不同类型的相互作用（Coulomb、Gaussian、delta等）都适用

---

## 2. ECG 基函数形式

### 2.1 Type I: 关联高斯（Correlated Gaussians）— 书中 Eq. (6.1)

用相对坐标 $\mathbf{x} = (x_1, ..., x_{N-1})$（Jacobi 坐标），基函数为：

$$
\text{Type I:} \quad \Psi = \exp\left(-\frac{1}{2}\tilde{x}Ax\right) \theta(x)
$$

其中 $A$ 是 $(N-1)\times(N-1)$ 的正定对称矩阵，$A = \tilde{G}DG$（正交 $G$，正对角 $D$）。$\theta(x)$ 是角函数部分。

### 2.2 Type II: 关联高斯型 geminals — 书中 Eq. (6.2)

用单粒子坐标 $\mathbf{r} = (r_1, ..., r_N)$，基函数为：

$$
\text{Type II:} \quad \Psi = \exp\left\{-\frac{1}{2}\tilde{r}Ar - \frac{1}{2}(\widehat{r-R})B(r-R)\right\} \theta(r-R)
$$

其中 $A_{ij} = -\alpha_{ij}$ ($i < j$)，$A_{ii} = \sum_{j \neq i} \alpha_{ji}$（注意 $\det A = 0$，$A$ 不正定），$B$ 是正对角矩阵（$B_{ii} = \beta_i$），$R$ 是位移参数。$C = A + B$ 必须正定。

**关键性质** (p.78)：Eq. (6.2) 可以简化为 $\exp\{-\frac{1}{2}(r-S)C(r-S)\}$（常数因子除外），其中 $C = A+B$，$S = (A+B)^{-1}BR$。因此 $A+B$ 必须正定。

### 2.3 我们代码中的形式

我们的 `BasisParams` 使用的是 $N \times N$ 矩阵 $A$、$B$（而非 $(N-1)\times(N-1)$ Jacobi 坐标矩阵）。我们的高斯函数形式是：

$$
\phi_k(\mathbf{x}) = \exp\left(-\mathbf{x}^T (A_k + i B_k) \mathbf{x} + i \mathbf{R}_k^T \mathbf{x}\right)
$$

这比标准 ECG 多了虚部 $B$ 和线性项 $R$，是为了支持实时间 TDVP 演化中的相位和动量。

### 2.4 关键区别

| 特征 | Varga Type I | Varga Type II | 我们的参数化 |
|------|------------|--------------|------------|
| 坐标 | Jacobi 坐标 | 单粒子坐标 | 单粒子坐标 |
| 矩阵大小 | $(N-1)\times(N-1)$ | $N\times N$ | $N\times N$ |
| A 矩阵约束 | 正定 | det A=0, A+B 正定 | 实部正定 |
| 虚部 | 无（实 A） | 无（实 A） | 有（复数 A+iB） |
| 线性项 | 无 | 位移 R | 有（iR） |
| 对称化 | 显式 P 求和 | 显式 P 求和 | 显式 P 求和 |

**我们的代码本质上是 Type II 的复数推广**，用于支持实时动力学。静态优化时 $B \to 0$, $R \to 0$，退化为标准 Type II。

---

## 3. 变分方法理论基础（书第 3 章）

### 3.1 能量方差作为质量指标 — Eq. (3.26)

能量方差定义为：

$$
\sigma^2 = \frac{\langle\Psi|H^2|\Psi\rangle}{\langle\Psi|\Psi\rangle} - E^2
$$

**Weinstein 定理 (Theorem 3.6)**：至少有一个本征值位于区间 $[E-\sigma, E+\sigma]$。

**Temple 下界 (Eq. 3.30)**：

$$
E - \frac{\sigma^2}{\varepsilon - E} \leq E_1
$$

其中 $\varepsilon < E_2$（第一激发态能量）。这给出了基态能量的严格下界。

**对我们的启示**：我们目前只用变分能量 $E$ 来评估基底质量。应该同时计算能量方差 $\sigma^2$ 和 Temple 下界，这样可以：
1. 量化结果的可靠性
2. 判断基底是否足够完备
3. 区分"真正的低能量"和"数值伪影"

### 3.2 Virial 定理 — Theorem 3.9

对谐振子势（$s=2$ 的齐次势），virial 定理给出：

$$
\langle T \rangle = \langle V \rangle = \frac{1}{2}E
$$

virial 比 $\eta = |\langle V_A \rangle / (2\langle T \rangle) - 1|$ 对精确解为零。**可以用来检查我们变分解的质量**。

### 3.3 Hellmann-Feynman 定理 — Theorem 3.7

$$
\frac{dE(\lambda)}{d\lambda} = \langle\Phi(\lambda)|\frac{\partial H(\lambda)}{\partial\lambda}|\Phi(\lambda)\rangle
$$

可以用来检查变分解的一致性：例如改变相互作用强度 $\lambda$，检查能量变化是否符合预期。

---

## 4. SVM 算法详解（书第 4 章）

**这是 Varga 方法最核心的要素。**

### 4.1 参数规模分析 (p.41)

对于最简单的关联高斯 Eq. (4.1)，每个基函数有 $N(N-1)/2$ 个非线性参数（$\alpha_{ij}$）。$K$ 个基函数的线性组合有 $K \times N(N-1)/2$ 个非线性参数。

例如 $N=4$, $K=200$：有 1200 个参数。**直接优化这么多参数几乎不可能**，局部极小值无处不在。

**这正是 SVM 随机搜索的核心优势**：避免被局部极小值困住。

### 4.2 SVM Competitive Selection（竞争选择）— 书 p.51

这是 SVM 的核心算法。逐个增加基函数（growth 阶段）：

```
s1. 随机生成 n 个候选参数集 (A_k^1, ..., A_k^n)
s2. 解 n 个 k 维广义本征值问题，得到能量 (E_k^1, ..., E_k^n)
s3. 选择能量最低的参数集 A_k^m 作为第 k 个基函数
s4. k → k+1
```

**关键效率优化 (Theorem 3.5, p.44)**：添加第 $k+1$ 个基函数时，不需要重新计算完整的 $H$ 和 $S$ 矩阵，只需要计算新行/列的矩阵元（$O(K)$ 而非 $O(K^2)$），也不需要重新对角化——可以直接在 $(K+1)$ 维空间中解本征值问题。

**我们的 `svm_build_basis` 已经正确实现了这一步。**

### 4.3 Utility Selection（效用选择）— 书 p.55

竞争选择的替代方案，增加了能量增益阈值：

```
s1. 随机生成参数集 A_k
s2. 计算 k 维本征值问题得到 ε₁(k)
s3. 只有当能量增益超过阈值 δε 时才接受：
    ε₁(k) < ε₁(k-1) - δε
s4. 如果连续 n 次都未通过测试，将 δε 减半，重新尝试
```

**优势**：避免接受与已有基底高度重叠的冗余基函数。**竞争选择**总是接受最佳候选，即使改善微乎其微，这可能导致线性相关。

**对我们的启示**：我们的 `svm_build_basis` 使用竞争选择。对大 K，应考虑切换到效用选择或添加重叠检查。

### 4.4 Refinement Cycle（精炼循环）— 书 p.53

**这是 SVM 精度的关键第二步。**

```
r1. 随机生成 n 个替代参数集 (A_i^1, ..., A_i^n)
r2. 将第 i 个基函数替换为每个候选，计算能量 E_K^1, ..., E_K^n
r3. 如果最佳候选的能量低于当前能量 → 接受替换
r4. 对 i = 1, 2, ..., K 循环执行 r1-r3
```

重复多轮直到一整轮都没有改进。

**书中的定量结果 (Table 4.6, p.54)**：

| Refinement 轮次 | Ps⁻ 能量 (a.u.) |
|---------|------------|
| 初始值（SVM growth） | -0.261992767 |
| 第 1 轮 | -0.262002742 |
| 第 2 轮 | -0.262003677 |
| 第 3 轮 | -0.262004063 |
| 第 10 轮 | -0.262004398 |
| 精确值 | -0.262005 |

K=100 基底，每轮 5 个 trial。10 轮 refinement 将误差从 1.2e-5 降到 6.0e-7。

### 4.5 不同优化策略的定量比较 (Tables 4.4, 4.5)

#### K=10 (Table 4.4, p.52)

| 方法 | 能量 (a.u.) | 矩阵元计算次数 | 时间 |
|------|----------|-----------|------|
| Powell 全优化 | -0.261251 | 214500 | 60 |
| SVM (n=1) | -0.244364 | 55 | 0.8 |
| SVM (n=10) | -0.251083 | 550 | 1.0 |
| Refining (n=10, 10轮) | **-0.261180** | 88055 | 11 |

#### K=100 (Table 4.5, p.53)

| 方法 | 能量 (a.u.) | 矩阵元计算次数 | 时间 |
|------|----------|-----------|------|
| Powell 全优化 | -0.26200016 | 35466150 | 7200 |
| SVM (n=1) | -0.26176191 | 11615 | 27 |
| SVM (n=10) | -0.26199427 | 86150 | 43 |
| Refining (n=10, 1轮) | **-0.26200231** | 805050 | 195 |

**关键结论**：
1. **SVM + Refining 以远少的计算量超越了 Powell 全优化**（K=100: 195 秒 vs 7200 秒，精度更好）
2. 增加 trial 数 n 从 1 到 10 能显著改善 SVM 质量
3. **Refining 对精度的贡献是决定性的**

### 4.6 Fine Tuning（微调）— 书 p.56

在 SVM 选中一个基函数后，可以用**确定性局部优化**（如 Powell 方法）在其邻域内进一步微调。当矩阵元计算不太昂贵时推荐使用。

**对我们的启示**：我们的 TDVP 实质上就在做这种 fine tuning（梯度下降优化非线性参数）。最优流程应该是：SVM growth → Refinement → TDVP fine tuning。

### 4.7 Excited States（激发态）— 书 Sect. 4.3, p.56

Mini-Max 定理保证：在 K 维基底中对角化 H，不仅给出基态的上界，也给出激发态的上界。可以用同样的 SVM 方法优化激发态，但需要在选择函数时考虑多个本征值。

### 4.8 线性相关问题 — 书 p.48, 50

**Random tempering** 虽然效果好，但"often leads to (almost) linearly dependent bases"，进一步改善困难。

**解决方案**：
1. 完全随机基底较少出现线性相关（p.48）
2. 重叠过大的 trial 应该被拒绝（p.50）："the overlap of the random basis functions should be smaller than a prescribed limit"
3. 使用效用选择而非竞争选择

**对我们的启示**：这正是我们 `stochastic_refine` 不稳定的根源。应该添加重叠检查：新 trial 与已有基函数的重叠 $|\langle\phi_{\text{new}}|\phi_i\rangle| / \sqrt{\langle\phi_{\text{new}}|\phi_{\text{new}}\rangle \langle\phi_i|\phi_i\rangle}$ 不应超过阈值（如 0.99）。

---

## 5. 坐标选择与变换（书第 2 章）

### 5.1 相对坐标的定义 — Eq. (2.4)

$$
x_i = \sum_j U_{ij} r_j, \quad r_i = \sum_j (U^{-1})_{ij} x_j
$$

### 5.2 两种标准坐标集

**Jacobi 坐标** $U_J$ (Eq. 2.5)：递归定义，$x_1 = r_1 - r_2$，$x_2 = (m_1 r_1 + m_2 r_2)/(m_1+m_2) - r_3$，...

**Heavy-particle center 坐标** $U_C$ (Eq. 2.7)：$x_i = r_i - r_N$（所有粒子相对于第 N 个粒子）。

### 5.3 粒子间距离表示 — Eq. (2.13)

$$
r_i - r_j = \widehat{w^{(ij)}} x
$$

其中 $\widehat{w^{(ij)}}$ 是 $1\times(N-1)$ 行矩阵。两体势矩阵元的计算需要用到这个表达式。

### 5.4 对称化 — Sect. 2.3

置换算符 $P$ 在相对坐标空间中表示为线性变换 $x' = T_P x$，其中 $T_P = U P U^{-1}$（截断到 $(N-1)$ 维）。对称化后的基函数：

$$
\Psi^{(\pm)} = \sum_P (\pm 1)^{[P]} \phi(T_P x)
$$

**我们的代码**（`PermutationSet`）已正确实现了这一步。

---

## 6. 矩阵元公式（书附录 A）

### 6.1 Overlap — Eq. (A.6)-(A.7)

对 $K = K' = 0$（纯高斯，无多项式前因子）：

$$
\langle f_0 | f_0 \rangle = \frac{(2L+1)!!}{4\pi} \left(\frac{(2\pi)^{N-1}}{\det B}\right)^{3/2} \rho^L
$$

其中 $B = A + A'$，$\rho = \tilde{u}'B^{-1}u$。

### 6.2 Kinetic Energy — Eq. (A.10)-(A.11)

对 $K = K' = 0$：

$$
\langle f_0 | T | f_0 \rangle = \frac{\hbar^2}{2}(R + LQ\rho^{-1}) \frac{(2L+1)!!}{4\pi} \left(\frac{(2\pi)^{N-1}}{\det B}\right)^{3/2} \rho^L
$$

其中 $R = 3\text{Tr}(AB^{-1}A'\Lambda)$，$P = -\tilde{u}B^{-1}A'\Lambda A'B^{-1}u$，$Q = -2\tilde{u}'B^{-1}A\Lambda AB^{-1}u$。

### 6.3 两体相互作用 — Eq. (A.21)

中心势 $V(|r_i - r_j|)$ 的矩阵元可以用 $w^{(ij)}$ 向量和积分 $I(n, c)$ 表示。

**与我们 1D 代码的对应**：我们的 `hamiltonian.cpp` 中的矩阵元公式是这些 3D 公式的 1D 特化。1D 情况下很多公式大幅简化（没有角动量耦合，$L=0$），所以我们的实现是正确的。

---

## 7. 与我们 TDVP 的关系

### 7.1 TDVP ≈ SVM 中的 gradient optimization

我们的 TDVP 步骤本质上是：
1. 计算 C 矩阵（metric tensor）和 g 向量（能量梯度）
2. 解 C dz = -g
3. 沿 dz 方向做线搜索

这与 GVM (Gradient Variational Method) 的思路非常相似：
- 都是在固定 K 的情况下优化非线性参数
- 都利用了能量关于参数的梯度信息

### 7.2 TDVP 不能替代 SVM

但 TDVP 不能做 SVM 的 growth 和 refinement：
- TDVP 不能增加基函数数量
- TDVP 不能做"整个替换"某个基函数的离散跳跃
- TDVP 容易卡在局部极小值

### 7.3 书中的最优流程 (p.55-56)

```
1. SVM Growth (competitive/utility selection): K=1 → K=K_target
2. SVM Refinement: 多轮，直到能量不再改善
3. Fine tuning (deterministic local opt / TDVP): 连续微调非线性参数
4. 如果还不够：增大 K_target，重复 1-3
```

我们目前只做了 1 和 3（而且 3 的效果受限于 1 的质量）。

---

## 8. 我们代码的现状与差距分析

### 8.1 已实现

| 功能 | 对应书中内容 | 状态 |
|------|----------|------|
| Competitive selection growth | Sect. 4.2.5 | `svm_build_basis` |
| Refinement cycle | Sect. 4.2.6 | `stochastic_refine`（已实现但有 bug） |
| Gradient optimization | Powell/GVM 的类似物 | TDVP solver |
| 对称化 | Sect. 2.3 | `PermutationSet` |
| 矩阵元 | Appendix A (1D 特化) | `hamiltonian.cpp` |
| S^{-1/2} 本征值求解 | Theorem 3.5 | `lowest_energy` |

### 8.2 缺失或需改进

| 功能 | 对应书中内容 | 紧迫性 |
|------|----------|--------|
| **Overlap 检查** | p.50: "overlap should be < prescribed limit" | **P0 — 根本原因** |
| **特征值截断** | — | **P0 — 改善 lowest_energy()** |
| **Utility selection** | Sect. 4.2.7 | P1 — 大 K 时需要 |
| **能量方差 σ²** | Sect. 3.2 | P1 — 质量诊断 |
| **Virial 检查** | Sect. 3.3 | P2 — 可选诊断 |
| **Temple 下界** | Eq. (3.30) | P2 — 可选诊断 |
| **多安排通道** | Sect. 4.2.1 | P3 — 三体以上 |

---

## 9. 具体代码修复建议

### 9.1 P0: 修复 `lowest_energy()` 的数值稳定性

当前 bug：refinement 时接受虚假负能量。根据书中的经验，有三个互补措施：

#### 措施 A：添加重叠检查（防止线性相关）

根据书 p.50 的建议，在 `stochastic_refine` 和 `svm_build_basis` 中添加：

```cpp
// 检查 trial 与已有基函数的归一化重叠
for (int i = 0; i < n; i++) {
    if (i == k) continue;  // 跳过被替换的
    double overlap_norm = std::abs(S_new(i, k))
                        / std::sqrt(std::abs(S_new(i,i)) * std::abs(S_new(k,k)));
    if (overlap_norm > 0.99) {
        // 拒绝此 trial
        break;
    }
}
```

#### 措施 B：S 矩阵特征值截断

当前 `lowest_energy()` 只做全/无的条件数检查。应改为**截断小特征值后在子空间求解**：

```cpp
// 截断阈值：w_i < w_max * rcond
double rcond = 1e-10;
int n_keep = 0;
for (int i = n-1; i >= 0; i--) {
    if (w(i) > w_max * rcond) n_keep++;
    else break;
}
// 只在 n_keep 维子空间中求解
```

这比当前的"全部保留或全部拒绝"更细致。

#### 措施 C：Rayleigh 商交叉验证（已实现）

当前代码已有此检查。但阈值 `1e-6 * |E0| + 1e-12` 可能对小能量差太宽松。

### 9.2 P1: 增大 K 并启用 Refinement

根据书中的定量数据（Table 4.4-4.5），我们应该：

```
当前设置:  K=5,  refine_rounds=0
目标设置:  K=20, refine_rounds=10, refine_trials=500
```

对应 main.cpp 中 Gaussian benchmark 的修改：

```cpp
run_svm_tdvp("2-particle Gaussian interaction", terms, 1.5266998310,
             /*K_max=*/20, /*svm_trials=*/5000,
             /*refine_trials=*/500, /*refine_rounds=*/10,
             /*tdvp_steps=*/1000,
             /*E_lower_bound=*/1.0);
```

**预期精度估计**：
- K=10 + no refine: ~1e-6
- K=10 + refine: ~1e-7
- K=20 + refine: ~1e-8
- K=20 + refine + TDVP: ~1e-9 或更好

### 9.3 P2: 添加能量方差诊断

实现 $\sigma^2 = \langle H^2 \rangle / \langle\Psi|\Psi\rangle - E^2$。需要计算 $\langle H^2 \rangle$ 矩阵元，这比 $\langle H \rangle$ 更复杂（包含 $T^2$, $TV$, $V^2$ 项）。

替代方案：**数值方差**。可以用 S^{-1/2} 变换后的 H_tilde 矩阵，$\sigma^2 = ||H_{\text{tilde}} c||^2 - E_0^2$，其中 $c$ 是归一化本征向量。这在代码中很容易实现。

### 9.4 P3: N-粒子推广的准备

当前 `random_basis_2particle` 是硬编码的 N=2。需要推广为 `random_basis_Nparticle`：

```cpp
BasisParams random_basis_Nparticle(int N, std::mt19937_64& rng, int name) {
    // 1. 生成 N(N-1)/2 个关联参数 α_ij（log-uniform）
    // 2. 构建 A 矩阵：A_ij = -α_ij (i<j), A_ii = Σ_j α_ji
    // 3. 加 diagonal B（trap 参数）
    // 4. 随机 R（位移）
    // 5. 确保 A+B 正定
}
```

---

## 10. 文献中的典型精度与基底大小

来自 Varga & Suzuki (1998), Tables 4.3-4.6 和 Mitroy et al. (2013) 的经验：

| 系统 | K（基函数数） | 精度 | 备注 |
|------|-------------|------|------|
| Ps⁻ (3体, 3D) K=100 SVM | 100 | 1.4e-5 | Table 4.3 |
| Ps⁻ (3体, 3D) K=100 SVM+refine | 100 | 6.0e-7 | Table 4.6, 10轮 |
| Ps⁻ (3体, 3D) K=600 SVM+refine | 600 | ~1e-10 | 书 p.54 |
| 2 体（平滑势, 1D） | 10~20 | 1e-8 ~ 1e-10 | 收敛快 |
| 2 体（Coulomb, 3D） | 20~50 | 1e-6 ~ 1e-8 | 有 cusp |
| 3 体 (3D) | 50~200 | 1e-4 ~ 1e-6 | N! = 6 |
| 4 体 (3D) | 200~500 | 1e-3 ~ 1e-5 | N! = 24 |
| 5-7 体 (3D) | 500~5000 | 1e-2 ~ 1e-4 | 需要大量基函数 |

**注意**：这些数据大多是 3D Coulomb 系统。对 1D 系统（如我们的问题），所需基函数通常更少。

关键观察：
- 增加 trial 数 n 从 1 到 10 是性价比最高的改善（Table 4.4: 能量从 -0.244 到 -0.251）
- Refinement 是精度的决定性因素（Table 4.4: 从 -0.251 到 -0.261）
- cusp 条件（如 delta 相互作用的 Kato cusp）会减慢收敛
- **1D 系统的 delta 相互作用会引入 cusp**，收敛速度类似于 3D Coulomb

---

## 11. 1D ECG 特殊公式

2019 年 Zhen et al. 发表了专门针对一维 ECG 的矩阵元公式：

> "Matrix Elements of One Dimensional Explicitly Correlated Gaussian Basis Functions"
> *Few-Body Systems* (2019)

该文给出了：
- 1D ECG 乘多项式前因子的 overlap、kinetic、potential 矩阵元
- 可用于直接计算 1D cold atom 系统的能量
- 也可用 tensor product 推广到 2D/3D 非球对称势

这对我们的 1D 代码直接适用。我们的矩阵元公式（`hamiltonian.cpp`）与此文献的结果应该一致，可以用来验证。

---

## 12. 关键参考文献

1. **书**：Y. Suzuki, K. Varga, *Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems*, Springer (1998)
   - [Springer](https://link.springer.com/book/10.1007/3-540-49541-X)

2. **原始论文**：K. Varga, Y. Suzuki, Phys. Rev. C **52**, 2885 (1995)
   - "Precise solution of few-body problems with the stochastic variational method on a correlated Gaussian basis"
   - [APS](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.52.2885)

3. **PRA 版本**：K. Varga, Y. Suzuki, Phys. Rev. A **53**, 1907 (1996)
   - [APS](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.53.1907)

4. **综述**：J. Mitroy et al., Rev. Mod. Phys. **85**, 693 (2013)
   - "Theory and application of explicitly correlated Gaussians"
   - [APS](https://link.aps.org/doi/10.1103/RevModPhys.85.693)

5. **1D 矩阵元**：Zhen et al., Few-Body Syst. (2019)
   - "Matrix Elements of One Dimensional Explicitly Correlated Gaussian Basis Functions"
   - [Springer](https://link.springer.com/article/10.1007/s00601-019-1539-3)

6. **梯度方法**：arXiv:2408.08522 (2024)
   - "Gradient Variational Methods for the Quantum Few-body Problem"
   - SVM + GVM 混合方法最优

7. **扩展 SVM**：arXiv:2512.07323 (2025)
   - "Extended Stochastic Variational Approach" for muonic systems
   - 详细描述了 growth + refinement 算法

8. **学士论文**（含 SVM 教程级别说明）：
   - [Study of small atoms using the Stochastic Variational Method](https://phys.au.dk/~fedorov/subatom/bachelor/david-svm.pdf)

---

## 13. 总结与行动计划

### 当前差距

我们的 Gaussian interaction 误差卡在 4.7e-6 的根本原因是 **K 太小 + refinement 未启用**。书中的定量数据表明 SVM+Refining 以极少的额外计算代价就能带来巨大精度提升。

### 行动优先级

```
[紧急] P0: 修复线性相关问题 → 添加 overlap 检查 + S 矩阵特征值截断
        → 这是启用 refinement 的前提

[高]   P1: 增大 K (5→20), 启用 refinement (10 轮)
        → 预期误差从 4.7e-6 降到 ~1e-8

[中]   P2: 添加能量方差诊断
        → 用 σ² 和 Temple 下界量化结果可靠性

[低]   P3: N-粒子推广
        → random_basis_Nparticle, 多安排通道
```

### 预期效果

| 配置 | 预期误差 | 计算量增加 |
|------|---------|---------|
| 当前 (K=5, no refine) | 4.7e-6 | baseline |
| K=15, no refine | ~1e-6 | 3x |
| K=15, refine 10轮 | ~1e-7 | 10x |
| K=20, refine 10轮 | ~1e-8 | 15x |
| K=20, refine + TDVP | ~1e-9 | 20x |
