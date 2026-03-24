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

## 仍然待查

Ghost 在 S 条件很好、没有取实部误差的情况下仍然出现（Round 1 就替换 5 个基函数后能量从 1.527 降到 1.000）。**问题不是数值精度，而是广义本征值问题本身给出了一个"合法但非物理"的解。**

下一步方向：理解什么是"非物理态"——为什么替换 5 个基函数后，数学上正确的最低本征值对应的是一个不同于目标基态的态。
