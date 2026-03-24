# 实验记录：S 矩阵条件数控制能否防止 Ghost State

日期：2026-03-24
分支：SVM_avoid_linear

---

## 假设

Ghost state 的机制是：S 近奇异 → S^{-1/2} 放大数值误差 → H̃ 本征值有噪声 → 贪心接受的棘轮效应累积噪声 → 能量系统性下漂。

如果在 SVM growth 和 stochastic refine 中，**拒绝所有会让 S 条件数超标的 trial**，就能从源头防止 ghost state。

## 实现

### 新增函数 `s_well_conditioned()`

```cpp
bool s_well_conditioned(const MatrixXcd& S, double max_cond = 1e6, double* w_min_out = nullptr);
```

- 对 S 做实部对称化后特征分解
- 检查 κ(S) = w_max / w_min < max_cond
- max_cond 默认 1e6

### 修改 `svm_build_basis()`

```diff
- if (has_excessive_overlap(S_ext, n, 0.99)) continue;
+ if (!s_well_conditioned(S_ext)) continue;
```

### 修改 `stochastic_refine()`

```diff
- if (has_excessive_overlap(S_new, k, 0.95)) continue;
- double E_trial = lowest_energy(H_new, S_new, 1e8, E_floor, 1e-8);
+ if (!s_well_conditioned(S_new)) continue;
+ double E_trial = lowest_energy(H_new, S_new, 1e8, E_lower_bound);
```

同时移除了 `max_round_drop` / `E_floor` 机制。

### 启用 refine

Gaussian 和 Delta benchmark 的 `refine_rounds` 从 0 改为 10。

## 测试结果（Gaussian benchmark）

```
--- Phase 1: SVM ---
K=18: E = 1.526702503  w_min=6.7237234823e-05

--- Phase 1.5: Stochastic Refinement ---
Refine: K=18, initial E=1.5267025028, w_min=6.7237234823e-05
Round 1: E = 1.000009925  w_min=3.0484201574e-04  replaced=5
Round 2: E = 1.000005519  w_min=2.9348107844e-04  replaced=1
Round 3: No improvement. Stopping.
```

## 关键观察

1. **S 的条件数始终很好。** w_min 从 6.7e-5 升到 3.0e-4（refine 后条件数反而更好了），κ 从 ~1.5e4 降到 ~3.3e3。
2. **Ghost state 仍然出现了。** 能量从 1.527 降到 1.000，偏差 0.527，远超精确值。
3. **sigma（能量方差）显示正常。** σ = 3.7e-7，看起来很小，没有报警。

## 误差量级分析

w_min = 6.7e-5 时，S^{-1/2} 的最大放大因子 = 1/√(6.7e-5) ≈ 122。

| 误差来源 | δ_H 量级 | H̃ 噪声 = δ_H / w_min | 能解释 0.5 漂移？ |
|---------|---------|---------------------|----------------|
| 机器精度 | 1e-16 | 2e-12 | 不能 |
| 置换求和累积 | 1e-14 | 2e-10 | 不能 |
| 取实部（若 Im(S) ~ 1e-3）| 1e-3 | 15 | **远超** |

**结论：w_min = 6.7e-5 对应的 S 条件数完全不构成问题。纯 S 病态导致的噪声量级只有 1e-12，不可能产生 0.5 的漂移。**

## 假设被证伪

原始假设（"S 近奇异是 ghost state 的根本原因"）被实验否定。S 条件数好的情况下 ghost 仍然出现。

## 新的怀疑方向

1. **取实部操作。** `lowest_energy_full` 中 `(0.5*(S+S†)).real()` 丢弃了 S 和 H 的虚部。如果基函数有非零的 B 或 R 参数，虚部不为零，取实部是一个 O(|Im|) 的误差。这个误差经过 S^{-1/2} 放大后可能很大——但这与 S 条件数无关，它是矩阵元本身的误差。

2. **基函数替换破坏了波函数结构。** Refine 替换了 18 个基函数中的 5 个。新的广义本征值问题可能给出一个完全不同的基态——不是原来波函数的改进版，而是一个全新的（可能非物理的）态。

3. **能量方差 σ² 没有报警的原因。** σ² 是在变换后的空间中算的，如果变换本身有问题（取实部引入的误差），σ² 也会跟着错。这和 Rayleigh 商交叉验证失效的原因相同：都是用同一套有偏的数据自我验证。

## 下一步

- [ ] 检查 refine 过程中 S 和 H 矩阵的虚部大小
- [ ] 尝试在 refine 阶段冻结 B=0, R=0（回到 Varga 的实参数设定）
- [ ] 对比 refine 前后波函数的结构变化
