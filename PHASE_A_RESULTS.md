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

## 已知问题

### Stochastic refine 数值不稳定
- 启用 refine 后，即使在原始阈值下，某些 trial basis 会产生 E ~ -1e+3 到 -1e+200 的虚假能量
- 原因：两个几乎线性相关的 basis function 使 overlap 矩阵接近奇异
- `lowest_energy()` 的阈值检查不够：通过了 `w_min > 1e-10` 和 `cond < 1e12`，但实际广义本征值问题仍然病态
- 修复方向：在 `stochastic_refine()` 中增加额外的物理合理性检查（如拒绝能量低于物理下界的 trial）

## 下一步方向

1. **修复 stochastic refine**：增加物理能量下界检查，安全增加 K 到 8-10
2. **改进 SVM 基底采样**：对 delta 问题用更多窄 Gaussian 覆盖 cusp 区域
3. **推进 Phase C**：实时间 kicked dynamics（不依赖 delta 精度的静态提升）
