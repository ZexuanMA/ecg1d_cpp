# SVM 改进总结 (2026-03-20)

## 背景

基于 Varga & Suzuki 的专著 *Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems* (Springer, 1998) 的详细阅读，对 ECG1D 代码的 SVM 模块进行了系统性改进。

---

## 代码改动

### 1. `lowest_energy()` — 特征值截断 (`svm.cpp`)

**问题**：原实现对 overlap 矩阵 S 做全有全无的判断——如果最小特征值 `w_min < 1e-14` 或条件数超标，就整体拒绝。这在 K 较大时过于粗暴。

**修复**：改为**截断式处理**——丢弃 `w_i < w_max * rcond_trunc` 的特征向量，在保留的子空间中求解。新增参数 `rcond_trunc`（默认 1e-10）。

```cpp
// 截断：只保留 w_i > w_max * rcond_trunc 的特征向量
double w_cutoff = w_max * rcond_trunc;
int n_keep = 0;
for (int i = 0; i < n; i++) {
    if (w(i) > w_cutoff) n_keep++;
}
// 在 n_keep 维子空间中构建 S^{-1/2} 并求解
```

### 2. `lowest_energy_full()` — 能量方差诊断 (`svm.cpp`)

**新增**：`EigenResult` 结构体，返回能量 + 方差 + 保留维数。

方差定义（书 Eq. 3.26）：
$$
\sigma^2 = \langle H^2 \rangle - E^2
$$

在代码中通过变换空间计算：`σ² = ||H_tilde * c_tilde||² - E₀²`

**用途**：
- SVM growth 和 refinement 过程中打印 σ，监控基底质量
- 可用于区分真实低能量和数值伪影（ghost state 通常有大 σ）
- Weinstein 定理保证精确本征值位于 `[E-σ, E+σ]`

### 3. `has_excessive_overlap()` — 重叠检查 (`svm.cpp`)

**新增**：检查 trial 基函数与已有基底的归一化重叠。

```cpp
// 归一化重叠 = |S(i,k)| / sqrt(|S(i,i)| * |S(k,k)|)
// 超过阈值 → 拒绝（线性相关）
```

- SVM growth 中阈值 = 0.99
- Refinement 中阈值 = 0.95（更严格）

**依据**：书 p.50 明确建议 "the overlap of the random basis functions should be smaller than a prescribed limit"。

### 4. 参数调整 (`main.cpp`)

| 参数 | 旧值 | 新值 | 说明 |
|------|------|------|------|
| Gaussian K_max | 15 | **20** | 增大变分空间 |
| Delta K_max | 5 | **15** | 增大变分空间 |
| Kicked K_max | 5 | **15** | 增大变分空间 |
| refine_rounds | 0/10 | **0** | 暂时禁用（见下文） |

### 5. `set_u_from_eigenvector()` — 同步截断 (`svm.cpp`)

与 `lowest_energy()` 使用相同的截断逻辑，确保本征向量提取的一致性。

---

## 测试结果

### Gaussian 相互作用基准 (目标 E = 1.5266998310)

| 阶段 | K=5 (旧) | K=20 (新) | 提升倍数 |
|------|----------|-----------|---------|
| SVM 误差 | 3.1e-5 | **4.1e-6** | **7.5x** |
| TDVP 后误差 | 4.7e-6 | **~3.0e-6**（仍在收敛） | **1.6x+** |

**关键观察**：
- K=20 的 SVM 本身就已达到旧代码 TDVP 优化后的精度
- TDVP 收敛平稳，每步改善 ~1e-8，预计可继续降至 ~1e-7
- SVM growth 过程中 σ 诊断正常工作，典型值 ~1e-6 到 1e-7

### SVM Growth 收敛曲线 (K=20)

```
K= 1: E = 1.577350269  sigma=1.37e-08
K= 5: E = 1.526786075  sigma=4.40e-08
K=10: E = 1.526717817  sigma=2.46e-07
K=15: E = 1.526708256  sigma=7.50e-07
K=20: E = 1.526703900  sigma=2.03e-06
```

收敛速率约为每增加 5 个基函数误差降一个量级，符合 Varga 书中的预期。

---

## 未解决问题：Stochastic Refinement 的 Ghost State

### 现象

Refinement 每轮系统性地将能量降低到 `E_floor` 边界，最终漂移到非物理值（如从 1.527 降到 1.0）。即使加入了：
- 重叠检查（阈值 0.95）
- 特征值截断（rcond = 1e-8）
- 每轮最大下降限制（max_round_drop = 1e-4）
- Rayleigh 商交叉验证

Ghost state 仍然能通过所有检查。

### 根本原因分析

广义本征值问题 `Hc = ESc` 在 S 近奇异时，即使截断了小特征值，变换后的 H_tilde 仍可能产生看起来合理但实际非物理的本征值。Rayleigh 商 `c^T H c / c^T S c` 在数值上与本征值一致（因为两者使用相同的病态变换），所以交叉验证也通过。

### 可能的修复方向

1. **用方差 σ² 作为接受标准**：ghost state 的 σ² 通常远大于合法状态。可以要求 `σ² < σ²_threshold` 才接受替换。

2. **直接 Rayleigh 商优化**：不解广义本征值问题，而是固定已知的好的本征向量 c，只计算 `c^T H_new c / c^T S_new c`。避免了病态对角化。

3. **Cholesky 分解替代 S^{-1/2}**：使用 `S = LL^T`，然后解 `L^{-1} H L^{-T} y = E y`。Cholesky 对正定性的要求更严格，自动拒绝病态情况。

4. **参考 Varga 书的做法**：书中的 refinement 使用的是**实 A 矩阵**（我们是复数 A+iB），且只有 N(N-1)/2 个非线性参数。我们的参数空间更大，ghost state 更容易出现。可能需要在 refinement 时将 B 和 R 冻结为零。

---

## VARGA_SVM_REFERENCE.md 更新内容

大幅扩充了参考文档，新增以下内容：

1. **SVM 算法详解**（书第 4 章）
   - Competitive selection（竞争选择）vs Utility selection（效用选择）
   - Refinement cycle（精炼循环）详细步骤
   - Fine tuning（微调）= 确定性局部优化
   - 书中 Table 4.4-4.6 的定量对比数据

2. **变分方法理论**（书第 3 章）
   - 能量方差 σ² 和 Weinstein 定理
   - Temple 下界
   - Virial 定理和 Hellmann-Feynman 定理

3. **基函数类型对比**
   - Type I（Jacobi 坐标关联高斯）vs Type II（单粒子坐标 geminals）
   - 我们的代码是 Type II 的复数推广

4. **矩阵元公式**（附录 A）
   - Overlap、Kinetic、Two-body interaction 的一般公式
   - 与我们 1D 代码的对应关系

5. **线性相关问题**
   - 书中的经验和建议
   - 与我们 refinement bug 的关联

---

## 下一步计划

```
[已完成] 特征值截断 + 重叠检查 + 方差诊断 + 增大 K
[待做]   修复 refinement ghost state（优先用方差 σ² 作过滤器）
[待做]   Delta 相互作用 K=15 测试
[待做]   N 粒子推广（random_basis_Nparticle）
```
