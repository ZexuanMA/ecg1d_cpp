# 能量方差作为变分计算质量指标：理论推导

参考：Suzuki & Varga, *Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems* (1998), Chapter 3

---

## 1. 符号与基本展开

设 Hamilton 量 $H$ 的精确本征态为 $\{|n\rangle\}$，对应本征值 $E_0 \leq E_1 \leq E_2 \leq \cdots$：

$$H|n\rangle = E_n|n\rangle$$

我们的试探波函数 $|\Psi\rangle$（由变分基底展开得到）可以用精确本征态展开（完备性）：

$$|\Psi\rangle = \sum_n c_n |n\rangle, \quad \sum_n |c_n|^2 = 1$$

**注意**：$c_n$ 是未知的——我们不知道精确本征态，自然也不知道展开系数。但 $c_n$ 的存在性由完备性保证，它们在理论推导中非常有用。

---

## 2. 变分能量

$$E \equiv \langle\Psi|H|\Psi\rangle = \sum_n |c_n|^2 E_n$$

这是各本征能量 $E_n$ 以概率 $|c_n|^2$ 加权的**统计平均值**。

**变分原理**（上界性质）：

$$E = \sum_n |c_n|^2 E_n \geq E_0 \sum_n |c_n|^2 = E_0$$

所以 $E \geq E_0$，变分能量永远是基态能量的**上界**。但仅知道上界，不知道差多少。

---

## 3. 能量方差的定义与物理含义

### 3.1 定义

$$\sigma^2 \equiv \langle\Psi|H^2|\Psi\rangle - E^2 = \langle\Psi|(H - E)^2|\Psi\rangle$$

用本征态展开：

$$\langle\Psi|H^2|\Psi\rangle = \sum_n |c_n|^2 E_n^2$$

所以：

$$\sigma^2 = \sum_n |c_n|^2 E_n^2 - \left(\sum_n |c_n|^2 E_n\right)^2$$

利用方差恒等式 $\mathrm{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}[X])^2 = \mathbb{E}[(X - \mathbb{E}[X])^2]$：

$$\boxed{\sigma^2 = \sum_n |c_n|^2 (E_n - E)^2}$$

### 3.2 物理含义

在量子力学中，如果对处于 $|\Psi\rangle$ 态的系统测量能量 $H$：
- 得到结果 $E_n$ 的概率是 $|c_n|^2$
- 测量结果的期望值是 $E = \sum_n |c_n|^2 E_n$
- 测量结果的方差就是 $\sigma^2 = \sum_n |c_n|^2 (E_n - E)^2$

**$\sigma^2$ 衡量的是：对 $|\Psi\rangle$ 做能量测量时，结果的分散程度。**

---

## 4. 定理：$\sigma^2 = 0$ 当且仅当 $|\Psi\rangle$ 是本征态

### 证明

$(\Rightarrow)$ 若 $|\Psi\rangle = |m\rangle$（某个本征态），则 $c_n = \delta_{nm}$，$E = E_m$，所以：

$$\sigma^2 = 1 \cdot (E_m - E_m)^2 = 0 \quad \checkmark$$

$(\Leftarrow)$ 若 $\sigma^2 = 0$。由 $\sigma^2 = \sum_n |c_n|^2(E_n - E)^2$，每一项 $\geq 0$，求和为零，故每项为零：

$$|c_n|^2(E_n - E)^2 = 0, \quad \forall n$$

对任意 $n$：要么 $c_n = 0$，要么 $E_n = E$。

设非零系数集合为 $S = \{n : c_n \neq 0\}$。对 $n \in S$，$E_n = E$。

若 $S$ 中的所有 $E_n$ 相同（$= E$），则 $|\Psi\rangle$ 是能量 $E$ 的本征子空间中的态，满足 $H|\Psi\rangle = E|\Psi\rangle$。

在非简并情况下，$|S| = 1$，即 $|\Psi\rangle$ 恰好是一个本征态。$\square$

### 推论

$$\sigma^2 \text{ 越小} \iff |\Psi\rangle \text{ 越接近精确本征态} \iff \text{变分能量越可信}$$

---

## 5. Weinstein 定理：$\sigma$ 给出误差上限

### 定理陈述

至少存在一个精确本征值 $E_m$ 落在区间 $[E - \sigma,\; E + \sigma]$ 内。

### 证明（反证法）

假设所有本征值都不在 $[E - \sigma, E + \sigma]$ 内，即对所有 $n$：

$$|E_n - E| > \sigma$$

代入方差表达式：

$$\sigma^2 = \sum_n |c_n|^2 (E_n - E)^2 > \sigma^2 \sum_n |c_n|^2 = \sigma^2$$

得到 $\sigma^2 > \sigma^2$，矛盾。$\square$

### 实际意义

如果 $|\Psi\rangle$ 主要是基态分量（$|c_0|^2 \approx 1$，变分法中通常成立），那么：

$$|E - E_0| \leq \sigma$$

**$\sigma$ 直接给出变分能量与最近本征值之差的上限。**

### 局限性

Weinstein 定理只保证"某个本征值"在区间内，不一定是基态。而且 $\sigma$ 给出的界往往偏宽松。更紧的界需要额外信息，引出下面的 Temple 下界。

---

## 6. Temple 下界：把 $E_0$ 夹在区间里

### 动机

变分原理给了上界 $E_0 \leq E$。如果能再找到一个**下界**，就可以把 $E_0$ 夹在一个已知区间内。

### 定理陈述

设 $\varepsilon$ 满足 $\varepsilon \leq E_1$（不超过第一激发态能量），且 $E < \varepsilon$，则：

$$\boxed{E_0 \geq E - \frac{\sigma^2}{\varepsilon - E}}$$

### 证明

考虑算符 $(H - E_0)(H - \varepsilon)$ 在 $|\Psi\rangle$ 上的期望值：

$$\langle\Psi|(H - E_0)(H - \varepsilon)|\Psi\rangle = \sum_n |c_n|^2 (E_n - E_0)(E_n - \varepsilon)$$

逐项分析：

| $n$ | $E_n - E_0$ | $E_n - \varepsilon$ | 乘积 |
|-----|-------------|---------------------|------|
| $n = 0$ | $= 0$ | $\leq 0$ | $= 0$ |
| $n \geq 1$ | $\geq 0$（因为 $E_n \geq E_1 \geq E_0$） | $\geq 0$（因为 $E_n \geq E_1 \geq \varepsilon$） | $\geq 0$ |

所以每一项 $\geq 0$，求和也 $\geq 0$：

$$\sum_n |c_n|^2 (E_n - E_0)(E_n - \varepsilon) \geq 0$$

展开左边：

$$\langle H^2 \rangle - (E_0 + \varepsilon)\langle H \rangle + E_0 \varepsilon \geq 0$$

代入 $\langle H \rangle = E$，$\langle H^2 \rangle = \sigma^2 + E^2$：

$$(\sigma^2 + E^2) - (E_0 + \varepsilon) E + E_0 \varepsilon \geq 0$$

整理各项：

$$\sigma^2 + E^2 - E_0 E - \varepsilon E + E_0 \varepsilon \geq 0$$

$$\sigma^2 + E(E - \varepsilon) - E_0(E - \varepsilon) \geq 0$$

$$\sigma^2 + (E - E_0)(E - \varepsilon) \geq 0$$

由条件 $E < \varepsilon$，所以 $(E - \varepsilon) < 0$。令 $\delta = E - E_0 \geq 0$（变分原理），移项：

$$\delta \cdot (\varepsilon - E) \leq \sigma^2$$

$$E - E_0 \leq \frac{\sigma^2}{\varepsilon - E}$$

$$\boxed{E_0 \geq E - \frac{\sigma^2}{\varepsilon - E}} \quad \square$$

### 最终结论：$E_0$ 的严格区间

结合变分上界，基态能量被夹在：

$$E - \frac{\sigma^2}{\varepsilon - E} \;\leq\; E_0 \;\leq\; E$$

区间宽度为 $\sigma^2 / (\varepsilon - E)$。

**区间变窄的两个途径**：
1. 减小 $\sigma^2$（增加基函数数量、做 refinement）
2. 增大 $\varepsilon - E$（能量间隙大的系统天然有利）

### $\varepsilon$ 的实际选取

$\varepsilon$ 需要满足 $\varepsilon \leq E_1$。实际中通常用**同一个变分计算的第二本征值**作为 $\varepsilon$。这是 $E_1$ 的上界（变分原理对激发态也成立），所以是一个保守但可用的选择。严格来说需要独立验证 $\varepsilon \leq E_1$，但在基底足够大时第二本征值通常是 $E_1$ 的好近似。

---

## 7. 数值举例

### 谐振子基态

精确值 $E_0 = 0.5$，$E_1 = 1.5$。用 $K = 5$ 个高斯基底做变分：

| 量 | 值 |
|---|---|
| 变分能量 $E$ | 0.500047 |
| $\sigma^2$ | $2.1 \times 10^{-6}$ |
| $\sigma$ | $1.45 \times 10^{-3}$ |
| 第二本征值 $\varepsilon$ | 1.52 |

**Weinstein**：

$$E_0 \in [0.500047 - 0.00145,\; 0.500047 + 0.00145] = [0.4986,\; 0.5015]$$

区间宽度 $\approx 2.9 \times 10^{-3}$。

**Temple 下界**：

$$E_0 \geq 0.500047 - \frac{2.1 \times 10^{-6}}{1.52 - 0.500047} = 0.500047 - 0.0000021 = 0.500045$$

$$E_0 \in [0.500045,\; 0.500047]$$

区间宽度 $\approx 2 \times 10^{-6}$，比 Weinstein 紧了三个数量级。

### 为什么 Temple 比 Weinstein 紧得多

- Weinstein 只用了 $\sigma$（方差的平方根），给出 $O(\sigma)$ 的界
- Temple 用了 $\sigma^2$（方差本身）除以能量间隙 $\varepsilon - E$，给出 $O(\sigma^2)$ 的界
- 当 $\sigma \ll 1$ 且间隙不太小时，$\sigma^2 / (\varepsilon - E) \ll \sigma$

---

## 8. 能量方差作为数值伪影检测器

### 问题背景

在 SVM 的 `stochastic_refine` 中，可能出现数值上近乎线性相关的基底。此时广义本征值问题 $Hc = ESc$ 病态，求出的"能量"可能是假的（虚假负能量）。

### 为什么 $\sigma^2$ 能区分真假

**真正的好结果**：$|\Psi\rangle$ 接近本征态 → $\sigma^2$ 小。

**数值伪影**：$S$ 矩阵病态 → 求得的"本征向量" $c$ 包含巨大的近似线性相关方向上的分量 → $\langle H^2 \rangle$ 爆炸 → $\sigma^2$ 巨大。

直觉上：假的低能量来自于线性相关分量的相消，这种精细的相消在计算 $H$ 时恰好给出一个假的低值，但计算 $H^2$ 时相消被平方项破坏，暴露出真实的数值不稳定性。

### 判据

$$\text{如果 } E \text{ 很低但 } \sigma^2 \gg |E|^2 \text{，则几乎确定是数值伪影。}$$

一个实用的阈值：$\sigma / |E| < 10^{-3}$（相对方差小于千分之一）表示结果可信。

---

## 9. 实际计算：为什么不需要知道精确本征态

### 9.1 一个自然的疑问

看到 $\sigma^2 = \sum_n |c_n|^2(E_n - E)^2$，第一反应可能是：要算 $\sigma^2$，岂不是得先知道精确本征值 $E_n$ 和展开系数 $c_n$？如果我已经知道了这些，还要变分法干什么？

**答案是：完全不需要。** 上面那个展开式只是理论推导用的（用来证明 Weinstein、Temple 等定理）。实际计算用的是**算符形式**，里面不含任何精确本征态的信息。

### 9.2 关键等式：两种等价表达

$$\underbrace{\sum_n |c_n|^2(E_n - E)^2}_{\text{本征态展开形式（理论推导用）}} \;=\; \underbrace{\langle\Psi|H^2|\Psi\rangle - E^2}_{\text{算符形式（实际计算用）}}$$

左边需要 $c_n$ 和 $E_n$（未知的），右边只需要算 $\langle H \rangle$ 和 $\langle H^2 \rangle$（在你自己的试探波函数上的期望值）。

**为什么相等？** 推导如下：

$$\langle\Psi|H^2|\Psi\rangle = \left(\sum_m c_m^* \langle m|\right) H^2 \left(\sum_n c_n |n\rangle\right)$$

因为 $H^2|n\rangle = E_n^2|n\rangle$，$\langle m|n\rangle = \delta_{mn}$：

$$= \sum_n |c_n|^2 E_n^2$$

所以：

$$\langle H^2\rangle - E^2 = \sum_n |c_n|^2 E_n^2 - \left(\sum_n |c_n|^2 E_n\right)^2 = \sum_n |c_n|^2(E_n - E)^2$$

推导过程插入了完备基 $\sum_n |n\rangle\langle n| = I$，但**计算右边时完全不需要这一步**——你只需要算算符 $H^2$ 在 $|\Psi\rangle$ 上的期望值。

**类比**：就像经典统计里算方差。你可以用定义式 $\mathrm{Var}(X) = \sum_i p_i(x_i - \mu)^2$（需要知道每个值 $x_i$ 和均值 $\mu$ 的偏差），也可以用 $\mathrm{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}[X])^2$（只需要一阶矩和二阶矩）。后者不需要拆解每个样本——你只要能算出总体的一阶矩和二阶矩就行。

### 9.3 $\langle H^2 \rangle$ 的具体计算

我们的变分波函数是 $|\Psi\rangle = \sum_k u_k |\phi_k\rangle$（$K$ 个基函数的线性组合），则：

$$\langle H \rangle = \frac{\mathbf{u}^\dagger H_{\text{mat}} \mathbf{u}}{\mathbf{u}^\dagger S_{\text{mat}} \mathbf{u}}, \qquad
\langle H^2 \rangle = \frac{\mathbf{u}^\dagger H_{\text{mat}} S_{\text{mat}}^{-1} H_{\text{mat}} \mathbf{u}}{\mathbf{u}^\dagger S_{\text{mat}} \mathbf{u}}$$

其中 $H_{\text{mat}}$ 和 $S_{\text{mat}}$ 是已有的 Hamiltonian 矩阵和 overlap 矩阵。

> **推导 $\langle H^2 \rangle$ 的表达式：**
>
> $$\langle\Psi|H^2|\Psi\rangle = \sum_{k,l} u_k^* u_l \langle\phi_k|H^2|\phi_l\rangle$$
>
> 但我们不想直接计算 $\langle\phi_k|H^2|\phi_l\rangle$（那需要新的矩阵元公式，涉及 $T^2$, $TV$, $VT$, $V^2$，非常复杂）。
>
> 技巧是插入单位算符（恒等分解）。在**我们自己的基底** $\{|\phi_k\rangle\}$ 中，恒等算符近似为 $I \approx \sum_{k,l} |\phi_k\rangle (S^{-1})_{kl} \langle\phi_l|$（这在基底完备时精确成立，在有限基底下是投影）。
>
> 因此：
>
> $$\langle\Psi|H \cdot I \cdot H|\Psi\rangle \approx \sum_{k,l,m,n} u_k^* u_n \langle\phi_k|H|\phi_l\rangle (S^{-1})_{lm} \langle\phi_m|H|\phi_n\rangle = \mathbf{u}^\dagger H_{\text{mat}} S_{\text{mat}}^{-1} H_{\text{mat}} \mathbf{u}$$
>
> 归一化后除以 $\mathbf{u}^\dagger S_{\text{mat}} \mathbf{u}$。

### 9.4 用 $S^{-1/2}$ 变换简化（代码中的实现方式）

代码中 `lowest_energy()` 已经在做 $S^{-1/2}$ 变换。设：

$$\tilde{H} = S^{-1/2} H_{\text{mat}} S^{-1/2}$$

对角化得到本征值 $E_0$ 和归一化本征向量 $\tilde{c}$（$\tilde{c}^\dagger \tilde{c} = 1$），则：

$$\langle H \rangle = E_0 = \tilde{c}^\dagger \tilde{H} \tilde{c}$$

$$\langle H^2 \rangle = \tilde{c}^\dagger \tilde{H}^2 \tilde{c} = \|\tilde{H}\tilde{c}\|^2$$

> **为什么 $\langle H^2 \rangle = \|\tilde{H}\tilde{c}\|^2$？**
>
> 原始波函数 $|\Psi\rangle = \sum_k u_k |\phi_k\rangle$，变换关系为 $\mathbf{u} = S^{-1/2}\tilde{c}$。
>
> $$\langle H^2 \rangle = \frac{\mathbf{u}^\dagger H S^{-1} H \mathbf{u}}{\mathbf{u}^\dagger S \mathbf{u}}$$
>
> 代入 $\mathbf{u} = S^{-1/2}\tilde{c}$：
> - 分子：$\tilde{c}^\dagger S^{-1/2} H S^{-1} H S^{-1/2} \tilde{c} = \tilde{c}^\dagger (S^{-1/2} H S^{-1/2})(S^{-1/2} H S^{-1/2}) \tilde{c} = \tilde{c}^\dagger \tilde{H}^2 \tilde{c}$
> - 分母：$\tilde{c}^\dagger S^{-1/2} S S^{-1/2} \tilde{c} = \tilde{c}^\dagger \tilde{c} = 1$
>
> 所以 $\langle H^2 \rangle = \tilde{c}^\dagger \tilde{H}^2 \tilde{c} = (\tilde{H}\tilde{c})^\dagger(\tilde{H}\tilde{c}) = \|\tilde{H}\tilde{c}\|^2$。

因此，能量方差只需一次矩阵-向量乘法：

$$\sigma^2 = \|\tilde{H}\tilde{c}\|^2 - E_0^2$$

```cpp
// 已有：H_tilde (变换后的 Hamiltonian), c_tilde (归一化本征向量), E0 (本征值)
VectorXcd Hc = H_tilde * c_tilde;          // 一次矩阵-向量乘法
double H2_expect = Hc.squaredNorm();        // = ||H̃c̃||² = ⟨H²⟩
double sigma2 = H2_expect - E0 * E0;        // 能量方差
double sigma = std::sqrt(std::max(sigma2, 0.0));  // 防止数值负值
```

### 9.5 总结：计算 $\sigma^2$ 的两条路

| 方法 | 需要什么 | 复杂度 | 精度 |
|------|---------|--------|------|
| 直接计算 $\langle\phi_k\|H^2\|\phi_l\rangle$ 矩阵元 | 新的矩阵元公式（$T^2$, $TV$, $V^2$） | 实现复杂 | 精确 |
| 用 $\tilde{H}^2$（即 $H S^{-1} H$ 的变换） | 已有的 $\tilde{H}$ 矩阵 | **一次矩阵乘向量** | 在基底完备性范围内精确 |

第二种方法几乎是免费的——你已经有了 $\tilde{H}$ 和 $\tilde{c}$，只需加两行代码。

---

## 10. 总结

| 工具 | 告诉你什么 | 需要什么 | 精度 |
|------|-----------|---------|------|
| 变分能量 $E$ | $E_0$ 的上界 | $\langle H \rangle$ | — |
| 能量方差 $\sigma^2$ | 波函数离本征态有多远 | $\langle H^2 \rangle$ | — |
| Weinstein 定理 | 某个 $E_n$ 在 $[E-\sigma, E+\sigma]$ 内 | $\sigma$ | $O(\sigma)$ |
| Temple 下界 | $E_0 \geq E - \sigma^2/(\varepsilon - E)$ | $\sigma^2$ + 能量间隙 | $O(\sigma^2)$ |

核心逻辑链：

$$\sigma^2 \text{ 小} \;\Rightarrow\; |\Psi\rangle \text{ 接近本征态} \;\Rightarrow\; E \approx E_0 \;\Rightarrow\; \text{结果可信}$$

变分能量告诉你"这是个上界"，能量方差告诉你"这个上界有多紧"。
