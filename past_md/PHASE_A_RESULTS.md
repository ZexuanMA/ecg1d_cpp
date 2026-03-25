# Phase A: Robust Static Optimizer — 实施结果

## 改动内容

### 1. `src/tdvp_solver.hpp` — SolverConfig 增强
- 新增 `energy_tol = 1e-12`：基于能量变化的收敛判据
- 新增 `stagnation_window = 50`：连续无改进步数触发恢复机制
- 新增 `adaptive_lambda = false`：自适应正则化开关
- 新增 `lambda_max = 1e-4`：正则化上界
- `dtao_max` 默认改为 `100 * dtao`（原为 `dtao`）

### 2. `src/tdvp_solver.cpp` — TDVP 求解器改进

**`tdvp_step()` 改进：**
- 当 `resolve_u=true` 时，每次试探步都重新求解 u（原来只在能量改善时才解）
- Line search 失败后增加 exploratory tiny step（dtao=1e-6），避免过早放弃

**`evolution()` 改进：**
- `dtao_max` 默认 = 100 * 初始步长，允许步长增长
- 能量收敛检测：连续 5 步 |dE| < energy_tol 时停止
- 停滞检测与恢复：连续 stagnation_window 步无改进时，增大 lambda_C 和 dtao
- Adaptive lambda：line search 失败时自动增大正则化
- 结束时打印诊断摘要（初始/最终能量、总 dE、line search 失败次数、停滞事件数、最小条件数）

### 3. `main.cpp` — 驱动参数
- TDVP 步数从 300 增加到 1000
- 启用 `resolve_u = true`、`adaptive_lambda = true`
- 设置 `dtao_max = 10.0`
- 注意：Stochastic refine 保持关闭（`refine_rounds = 0`），因为存在数值稳定性问题

### 4. 未修改 `src/svm.cpp`
- `lowest_energy()` 保持原始阈值（`w_min < 1e-10`, `cond > 1e12`）
- 尝试过放松阈值，但会导致 stochastic refine 接受虚假负能量（E ~ -1e+200）
- 根本原因：shift/truncate overlap eigenvalues 后用原始 H 做广义本征值问题，数学上不自洽

## Benchmark 结果对比

### 1-particle harmonic oscillator (精确值 E = 0.5)

| 版本 | TDVP 误差 |
|------|-----------|
| 改前 | 2.2e-6 |
| **改后** | **1.8e-12** |

提升约 **100 万倍**，接近机器精度。能量收敛检测在第 123 步正确触发。

### 2-particle Gaussian interaction (精确值 E = 1.5266998310)

| 阶段 | 改前 error | 改后 error |
|------|-----------|-----------|
| SVM (K=5) | 3.1e-5 | 3.1e-5 |
| TDVP | 8.1e-6 | **4.7e-6** |

TDVP 误差改善约 40%。1000 步结束时 dE 仍在 ~1e-10 量级，缓慢下降中。

### 2-particle delta contact (精确值 E = 1.3067455)

| 阶段 | 改前 error | 改后 error |
|------|-----------|-----------|
| SVM (K=5) | 3.1e-3 | 3.1e-3 |
| TDVP | 3.0e-3 | **3.0e-3** |

基本无改善。TDVP 1000 步只降了 1.4e-4，最终 dE ~1e-9 但距目标还有 3e-3。

## 核心结论

1. **TDVP 优化器改进有效**：1-particle 从 2.2e-6 到 1.8e-12 证明优化器本身大幅改善
2. **Gaussian interaction 受益有限**：K=5 基函数数量限制了进一步改善空间
3. **Delta 瓶颈确认是基底表达力**：优化器已收敛（dE~1e-9），但能量距目标 3e-3——Gaussian 基底无法有效表达 cusp
4. **Stochastic refine 有稳定性 bug**：放松 `lowest_energy()` 阈值会导致接受虚假负能量，需要单独修复

## Ghost State Bug 修复（2026-03-24）

### 根本原因

`perturb_basis()` 中的 A 矩阵对称化语句：

```cpp
A_new = 0.5 * (A_new + A_new.transpose());  // BUG: Eigen aliasing
```

Eigen 的惰性求值导致 `A_new.transpose()` 在 `A_new` 被写入时同时被读取，结果不对称。

**修复**：加 `.eval()` 强制先求值：

```cpp
A_new = (0.5 * (A_new + A_new.transpose())).eval();
```

### 因果链

A 不对称 → K 矩阵不对称 → K_inv 不对称 → H kernel 不满足 Hermitian 对称性 → build_HS 用 conj 填下三角掩盖了真实不对称 → H 矩阵不代表真实 Hamiltonian → H_g 出现负特征值 → 变分原理被打破 → refine 找到 E < E_true 的"解"

详见 `EXPERIMENT_S_CONDITIONING.zh.md`。

### 附带改进

| 改动 | 说明 |
|------|------|
| 复 Hermitian 求解 | `lowest_energy_full` 去掉 `.real()`，用 `SelfAdjointEigenSolver<MatrixXcd>` |
| 全量重建 | `stochastic_refine` 用 `build_HS` 替代增量行/列更新 |
| S 条件数检查 | `s_well_conditioned()` 替代 pair-wise overlap 检查 |

## 修复后 Benchmark 结果（2026-03-24）

### 2-particle Gaussian interaction (精确值 E = 1.5266998310)

| 阶段 | K | 误差 | 耗时 |
|------|---|------|------|
| SVM | 20 | 2.31e-6 | ~8s |
| + Refine (10 轮) | 20 | ~2.1e-6 | ~70s |
| + TDVP (收敛中) | 20 | ~2.0e-6 | — |

Refine 每轮稳步改善约 2e-8，10 轮累计 ~2e-7。TDVP 每步约 5e-10。
变分原理始终满足（E 单调递减且 > E_exact）。

### K 增大的收益分析

测试 K_max=30（条件数阈值 1e10）：

```
K=20: E=1.526702145  error=2.31e-6  w_min=8.1e-7
K=21: E=1.526702143  error=2.31e-6  w_min=8.1e-7  (改善 2e-9)
K=22: E=1.526702143  error=2.31e-6  w_min=8.1e-7  (无改善)
K=23: E=1.526702143  error=2.31e-6  w_min=8.1e-7  (无改善)
K=24: no valid trial found
```

**K=20 之后收益为零**——误差卡在 2.31e-6。原因不是 K 不够大，而是 `random_basis_2particle` 的采样分布已经覆盖不到能继续降低能量的参数区域。K=21~23 加的基函数对能量几乎没有边际贡献。

突破 2.3e-6 的可能方向：
- 更多轮 refine（在已有 K=20 基础上精炼参数，比盲目加基函数更高效）
- 改进采样分布（pair-coordinate-aware sampling）
- 增大 n_trials（当前 5000，可以试 20000）

### 收敛曲线（Refine, K=18 版本）

```
Refine初始: E = 1.5267025028
Round  1:   E = 1.526702494   replaced=7
Round  2:   E = 1.52670242    replaced=9
Round  5:   E = 1.52670233    replaced=5
Round 10:   E = 1.526702308   replaced=5
```

## N 粒子扩展性评估

| N | N! | K 需求 | 典型精度 | 可行性 |
|---|-----|--------|---------|--------|
| 2 | 2 | 10-20 | 1e-6 ~ 1e-10 | **已实现** |
| 3 | 6 | 50-200 | 1e-4 ~ 1e-6 | 可行（分钟级） |
| 4 | 24 | 200-500 | 1e-3 ~ 1e-4 | 有挑战（小时级） |
| 5-6 | 120-720 | 500-5000 | 1e-2 ~ 1e-3 | 极限（天级） |
| 10 | 3.6M | — | — | 不可行（需换方法） |

当前 N=2 精度 2.5e-6 对静态基态和动力学初态准备都够用。
瓶颈不在精度，在 N! 置换求和的计算量。

## 下一步方向

1. **跑 Delta 和 Kicking benchmark**：验证 refine 在其他 Hamiltonian 下也正常
2. **推进实时间动力学**：kicked TDVP 演化（Phase C）
3. **N=3 扩展**：通用 `random_basis_Nparticle`，验证 N=3 可行性
