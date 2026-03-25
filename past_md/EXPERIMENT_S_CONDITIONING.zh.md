# 实验记录：S 矩阵条件数控制能否防止 Ghost State

日期：2026-03-24
分支：SVM_avoid_linear

---

## 假设 A：S 近奇异是 Ghost State 的根本原因

Ghost state 的假设机制：S 近奇异 → S^{-1/2} 放大数值误差 → H̃ 本征值有噪声 → 贪心接受的棘轮效应累积噪声 → 能量系统性下漂。

如果在 SVM growth 和 stochastic refine 中，**拒绝所有会让 S 条件数超标的 trial**，就能从源头防止 ghost state。

### 实现

#### 新增函数 `s_well_conditioned()`

```cpp
bool s_well_conditioned(const MatrixXcd& S, double max_cond = 1e6, double* w_min_out = nullptr);
```

- 对 S 做 Hermitian 对称化后特征分解
- 检查 κ(S) = w_max / w_min < max_cond
- max_cond 默认 1e6

#### 修改 `svm_build_basis()` 和 `stochastic_refine()`

```diff
- if (has_excessive_overlap(S_ext, n, 0.99)) continue;
+ if (!s_well_conditioned(S_ext)) continue;
```

同时移除了 `max_round_drop` / `E_floor` 机制，启用 refine_rounds=10。

### 测试结果

```
--- Phase 1: SVM ---
K=18: E = 1.526702503  w_min=6.7237234823e-05

--- Phase 1.5: Stochastic Refinement ---
Refine: K=18, initial E=1.5267025028, w_min=6.7237234823e-05
Round 1: E = 1.000009925  w_min=3.0484201574e-04  replaced=5
Round 2: E = 1.000005519  w_min=2.9348107844e-04  replaced=1
Round 3: No improvement. Stopping.
```

### 结果：假设 A 被证伪

1. **S 的条件数始终很好。** w_min 从 6.7e-5 升到 3.0e-4（refine 后条件数反而更好了）
2. **Ghost state 仍然出现。** 能量从 1.527 降到 1.000，偏差 0.527
3. **sigma 没有报警。** σ = 3.7e-7，看起来很小

#### 误差量级分析

w_min = 6.7e-5 → S^{-1/2} 最大放大 = 1/√(6.7e-5) ≈ 122。

| 误差来源 | δ_H 量级 | H̃ 噪声 = δ_H / w_min | 能解释 0.5 漂移？ |
|---------|---------|---------------------|----------------|
| 机器精度 | 1e-16 | 2e-12 | 不能 |
| 置换求和累积 | 1e-14 | 2e-10 | 不能 |

**S 病态导致的噪声只有 ~2e-12，不可能产生 0.5 的漂移。S 条件数不是 ghost 的原因。**

---

## 假设 B：取实部 `.real()` 丢弃了虚部信息

`lowest_energy_full` 中原来做了 `(0.5*(S+S†)).real()`，丢弃 S 和 H 的虚部。如果虚部不为零，这会引入误差。

### 实现

将三个函数（`lowest_energy_full`、`set_u_from_eigenvector`、`s_well_conditioned`）从：

```cpp
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S.real());
```

改为：

```cpp
Eigen::SelfAdjointEigenSolver<MatrixXcd> es(S);  // 复 Hermitian，特征值仍是实数
```

所有中间变量（V_keep、S_inv_half、H_tilde、c）改为复数类型。S^{-1/2} = V diag(1/√w) V†（adjoint，不是 transpose）。

### 测试结果

```
--- Phase 1: SVM ---
K=18: E = 1.526702503  w_min=6.7237234821e-05  (与假设A几乎相同)

--- Phase 1.5: Stochastic Refinement ---
Round 1: E = 1.000009925  sigma=1.8821527846e-07  w_min=3.0484201574e-04  replaced=5
Round 2: E = 1.000005519  sigma=5.5716284647e-07  w_min=2.9348107844e-04  replaced=1
Round 3: No improvement. Stopping.
```

### 结果：假设 B 也被证伪

行为与取实部版本**完全相同**——相同的能量、相同的 w_min、相同的 ghost。

这证实了：**SVM/refine 阶段的基函数参数确实是实的**（B 虚部 = 0，R 虚部 = 0），所以 `.real()` 本来就没丢信息。去掉 `.real()` 是代码质量改进（更正确、更通用），但不影响 ghost state。

---

## 排除清单

| 假设 | 结果 | 证据 |
|------|------|------|
| S 近奇异 → 噪声放大 | **排除** | w_min=6.7e-5，噪声仅 2e-12 |
| `.real()` 丢虚部 | **排除** | 去掉后行为完全相同 |
| Rayleigh 商交叉验证 | 无效 | 自我验证，抓不住 ghost |
| sigma 方差检查 | 无效 | σ=1.9e-7，没有报警 |

---

## 假设 C：incremental H/S 更新有 bug

### 发现

在 refine 后添加诊断：从 basis 重建 H/S，与 incremental 存储的 H/S 对比。

```
Refine result (stored H/S): E = 1.000005519
Refine result (rebuilt H/S): E = -5.724651375
*** MISMATCH: differ by 6.72 ***
```

### 修复

将 incremental 更新（只重算第 k 行/列）替换为全量重建 `build_HS()`。修复后 stored = rebuilt。

### 结果：bug 修复，但 ghost 仍在

修复后 E 仍然 = 1.000。矩阵一致了，说明**这个 E=1.000 是当前 H/S 矩阵的真实最低本征值**——但它违反变分原理（E_true = 1.527）。

---

## 假设 D：H 矩阵元计算不满足 Hermitian 对称性

### 发现 1：能量分解异常

```
<T>     =  2.673    (正，合理)
<V_har> =  1.244    (正，合理)
<V_int> = -2.917    (负！物理上 Gaussian 排斥势的期望值必须 > 0)
<T+V>   =  3.917
<Total> =  1.000
```

**Gaussian 相互作用能量为负**，物理上不可能。H_g 矩阵有负特征值（-0.069, -0.003）。

### 发现 2：逐 kernel Hermitian 对称性测试

```
         kinetic: h(a,b)=2.3305  h(b,a)=2.3369  |diff|=6.4e-3  ✗
        harmonic: h(a,b)=0.8132  h(b,a)=0.8132  |diff|=2.2e-16 ✓
        gaussian: h(a,b)=1.9623  h(b,a)=1.9508  |diff|=1.1e-2  ✗
```

kinetic 和 gaussian 不对称，harmonic 完美对称。

### 发现 3：逐置换分析

```
P=0 (identity): M_G 对称 ✓, kernel 对称 ✓
P=1 (swap):     M_G 对称 ✓, kernel 不对称 ✗ (diff=0.008)
```

只有 swap 置换的 kernel 不对称。

### 发现 4：K 矩阵不对称

```
K(0,1) = -0.4129
K(1,0) = -0.0564    ← 差 0.357！
```

K = conj(A_i+B_i) + P^T(A_j+B_j)P。如果 A 对称，K 应该对称。

### 发现 5（根本原因）：A 矩阵不对称

```
basis[0].A = [[0.911, -0.088],
              [0.134,  0.528]]   ← A(0,1)=-0.088, A(1,0)=0.134, 差 0.22!

basis[1].A = [[0.869, -0.325],
              [-0.190, 0.431]]   ← A(0,1)=-0.325, A(1,0)=-0.190, 差 0.135!
```

`random_basis_2particle` 和 `perturb_basis` 都显式对称化了 A。所以 A 的不对称一定是在其他地方引入的。

### 完整因果链

```
A 矩阵不对称（源头待查）
    → K = conj(A_i+B_i) + P^T(A_j+B_j)P 不对称
    → K_inv 不对称
    → kernel(K_inv 元素, mu 元素) 对交换 bra/ket 不对称
    → h(i,j) ≠ conj(h(j,i))
    → build_HS 用 conj 填下三角，掩盖了真实的不对称
    → H 矩阵不代表真实 Hamiltonian
    → H_g 出现负特征值
    → 变分原理被打破
    → refine 找到 E < E_true 的"解"
```

---

## 排除清单（更新）

| 假设 | 结果 | 证据 |
|------|------|------|
| A: S 近奇异 → 噪声放大 | **排除** | w_min=6.7e-5，噪声仅 2e-12 |
| B: `.real()` 丢虚部 | **排除** | 去掉后行为完全相同 |
| C: incremental H/S 更新 bug | **已修复** | 改为全量重建 |
| D: H 矩阵元不 Hermitian | **确认** | A 矩阵不对称是根本原因 |

---

## 假设 E（最终确认）：Eigen aliasing bug 导致 A 矩阵不对称

### 定位过程

1. SVM 结束后 A 完全对称（||A-A^T|| = 0）
2. Refine 接受的 trial A 不对称（||A-A^T|| ~ 0.2）
3. `random_basis_2particle` 无 BUG 报告
4. **`perturb_basis` 大量产出不对称 A**

### 根本原因

`perturb_basis` 第 320 行：

```cpp
A_new = 0.5 * (A_new + A_new.transpose());  // 意图对称化
```

**Eigen 的惰性求值（lazy evaluation）导致 aliasing**：`A_new` 在表达式右侧通过 `.transpose()` 被读取的同时，左侧的 `A_new =` 正在写入。Eigen 不会自动检测这种自赋值 aliasing，导致结果不对称。

### 修复

```cpp
A_new = (0.5 * (A_new + A_new.transpose())).eval();  // .eval() 强制先算完右侧
```

### 验证结果

```
修复前（ghost state）:
  SVM:     E = 1.526702503
  Round 1: E = 1.000007499  ← 跳到无相互作用极限

修复后（正常工作）:
  SVM:     E = 1.526702503
  Round 1: E = 1.526702494  ← 稳步下降 9e-9
  Round 4: E = 1.526702348  ← 累计改善 1.5e-7
```

**变分原理恢复，stochastic refine 正常工作。**

---

## 最终排除清单

| 假设 | 结果 | 证据 |
|------|------|------|
| A: S 近奇异 → 噪声放大 | **排除** | w_min=6.7e-5，噪声仅 2e-12 |
| B: `.real()` 丢虚部 | **排除** | 去掉后行为完全相同 |
| C: incremental H/S 更新 bug | **已修复** | 改为全量重建（独立于根本原因） |
| D: H 矩阵元不 Hermitian | **确认为症状** | A 不对称 → K 不对称 → kernel 不对称 |
| **E: Eigen aliasing in perturb_basis** | **根本原因** | `.eval()` 一行修复一切 |
