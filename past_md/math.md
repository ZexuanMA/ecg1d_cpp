# Kicked ECG 主链路数学说明

这份文档只讲当前主链路，不讲旁支实验。

目标是回答 5 个问题：

1. 虚时间之后我们手里到底有什么
2. 怎么从这个基态表示出发扩增基底
3. 扩增后怎么选基
4. kick 后的态具体怎么投影回当前基底
5. 后续固定基自由演化怎么做

---

## 0. 总体框架

当前总波函数始终写成有限个基函数的线性组合：

$$
\Psi(\mathbf{x}) = \sum_{j=1}^{K} u_j \,\phi_j(\mathbf{x})
$$

这里：

- $\phi_j$ 是第 $j$ 个基函数
- $u_j$ 是线性系数
- 当前代码里一个基函数由 `BasisParams = (u, A, B, R, name)` 描述

所以整个问题始终有两层：

1. 非线性层：basis functions 自身长什么样
2. 线性层：这些 basis functions 之间怎么线性组合

---

## 1. 单个基函数长什么样

当前基函数形式是：

$$
\phi(\mathbf{x})
=
\exp\!\Big[
-\mathbf{x}^T (A + B)\mathbf{x}
+ 2\,\mathbf{R}^T B \mathbf{x}
- \mathbf{R}^T B \mathbf{R}
\Big]
$$

其中：

- $A$ 是对称矩阵，控制 Gaussian 宽度和相关
- $B$ 当前实现里通常是对角矩阵
- $R$ 是复向量

### 1.1 为什么 $R$ 的虚部编码动量

看第 $a$ 个粒子的线性项：

$$
\exp\!\bigl[ 2 R_a B_{aa} x_a \bigr]
$$

如果取

$$
R_a = i\,\frac{m_a k_L}{b}, \qquad B_{aa}=b
$$

那么

$$
2 R_a B_{aa} x_a
=
2\left(i\frac{m_a k_L}{b}\right)b x_a
=
i\,2m_a k_L x_a
$$

于是基函数里出现

$$
e^{i\,2m_a k_L x_a}
$$

这就是平面波因子。

所以“动量基函数”的本质就是：

$$
\text{Gaussian 包络} \times \text{平面波}
$$

也就是

$$
\phi(x) \sim \text{Gaussian}(x)\,e^{iqx}.
$$

---

## 2. 虚时间之后我们手里到底有什么

在 kicked 主线开始之前，代码先做三步静态准备：

1. `svm_build_basis()`
2. `stochastic_refine()`
3. 虚时间 `evolution(...)`

这三步结束后，我们手里有两样东西。

### 2.1 一组适合表示基态的 basis functions

记为

$$
\mathcal{B}_{\text{gs}} = \{\phi_1,\dots,\phi_{K_0}\}
$$

这里 “gs” 的意思不是“某一个基函数就是基态”，而是：

> 这一组 basis functions 经过静态优化，比较适合表示 free Hamiltonian 的基态。

### 2.2 这组 basis 上表示基态的系数向量

记为

$$
\mathbf{u}^{(0)} = (u_1,\dots,u_{K_0})^T
$$

这不是手调的，而是通过广义本征值问题求出来的。

---

## 3. 在固定一组 basis 上，基态怎么求

给定 basis 集合 $\{\phi_j\}$，先构造：

$$
S_{ij} = \langle \phi_i | \phi_j \rangle
$$

和

$$
H_{ij} = \langle \phi_i | \hat H_0 | \phi_j \rangle
$$

其中 $\hat H_0$ 是 free Hamiltonian。

因为 basis 非正交，所以基态不是解普通本征值问题，而是解广义本征值问题：

$$
H \mathbf{u} = E\,S\mathbf{u}
$$

最低本征值 $E_0$ 对应基态能量，最低本征向量 $\mathbf{u}_0$ 对应基态系数。

### 3.1 数值上怎么解

代码不是直接硬解 $H u = E S u$，而是先正交化 $S$：

$$
S = V \Lambda V^\dagger
$$

构造

$$
S^{-1/2} = V \Lambda^{-1/2} V^\dagger
$$

定义

$$
\widetilde H = S^{-1/2} H S^{-1/2}
$$

于是变成普通 Hermitian 本征值问题：

$$
\widetilde H \mathbf{c} = E \mathbf{c}
$$

最后回到原系数：

$$
\mathbf{u} = S^{-1/2}\mathbf{c}
$$

当前代码里对应：

- `build_HS(...)`
- `lowest_energy_full(...)`
- `set_u_from_eigenvector(...)`

---

## 4. 从基态 basis 出发，怎么扩增基底

当前 kicked 主线采用：

$$
\mathcal{B}_{\text{gs}}
\longrightarrow
\mathcal{B}_{\text{aug}}
=
\mathcal{B}_{\text{gs}} \cup \mathcal{B}_{\text{mom}}
$$

也就是：

- 先保留原来的基态 basis
- 再加上一批带动量的 basis functions

### 4.1 为什么要加动量基

kick 算子是

$$
\hat U_{\text{kick}}
=
\exp\!\left[-i\kappa \sum_{a=1}^N \cos(2k_L x_a)\right]
$$

对单粒子那部分，用 Bessel 展开：

$$
e^{-i\kappa\cos\theta}
=
\sum_n (-i)^n J_n(\kappa)e^{in\theta}
$$

取 $\theta = 2k_L x_a$，就看到 kick 会产生这些谐波：

$$
0,\ \pm 2k_L,\ \pm 4k_L,\ \dots
$$

所以如果 basis 里没有这些非零动量方向，kick 后的态就承接不住。

---

## 5. 当前代码里，候选动量基是怎么生成的

当前主干代码走的是“启发式动量基生成”，不是 exact kick image。

### 5.1 width 候选

先准备一组 width 候选：

$$
W = \{0.3,\ 0.8,\ 2.0,\ 6.0\}
\cup
\{\text{原基态 basis 中 } A_{aa}\text{ 的对角值}\}
$$

这里第二部分的意思是：

> 从已经优化好的基态 basis 里，借用系统自己“学出来”的空间尺度。

### 5.2 动量标签候选

对 $N=1$：

$$
m \in \{-n_{\text{mom}},\dots,n_{\text{mom}}\}\setminus\{0\}
$$

对 $N=2$：

$$
(m_1,m_2),\qquad
m_1 \le m_2,\qquad
m_1,m_2 \in \{-n_{\text{mom}},\dots,n_{\text{mom}}\},
$$

并跳过

$$
(0,0)
$$

### 5.3 每个候选基函数具体长什么样

对每组 $(w,m_1,m_2)$，当前代码构造

$$
A_{\text{new}}
=
\begin{pmatrix}
w & 0.1w \\
0.1w & w
\end{pmatrix},
\qquad
B_{\text{new}}
=
\begin{pmatrix}
b & 0 \\
0 & b
\end{pmatrix}
$$

并取

$$
R_{\text{new}}
=
\begin{pmatrix}
i\,m_1 k_L / b \\
i\,m_2 k_L / b
\end{pmatrix}
$$

于是这个候选带有动量相位

$$
e^{i\,2m_1k_Lx_1} e^{i\,2m_2k_Lx_2}.
$$

当前代码里对应：

`augment_basis_with_momentum(...)`

---

## 6. 扩增后怎么选“合适的基底”

当前主干代码不是把所有候选都收下，而是做一个简单贪心。

设当前已经接受的 basis 是

$$
\mathcal{B}^{(n)}
$$

对某个候选 $\chi$：

1. 形成试探 basis

$$
\mathcal{B}_{\text{trial}} = \mathcal{B}^{(n)} \cup \{\chi\}
$$

2. 构造它的 overlap 矩阵

$$
S_{\text{trial}}
$$

3. 如果

$$
\lambda_{\min}(S_{\text{trial}}) > 0
$$

并且

$$
\mathrm{cond}(S_{\text{trial}})
=
\frac{\lambda_{\max}(S_{\text{trial}})}{\lambda_{\min}(S_{\text{trial}})}
<
\texttt{max\_cond}
$$

就接受；否则拒绝。

### 6.1 这一步真正优化的是什么

当前这个筛选优化的不是：

- first-kick 的 `dx2`
- first-kick 的 `dE`
- first-kick 的 `dp2`

它优化的只是：

> 不要把 overlap 矩阵搞得太病态，以至于后面的广义本征值问题和投影方程没法算。

所以这是一个“数值稳定性筛选准则”，不是“物理最优选基准则”。

---

## 7. 扩增后为什么还要再解一次基态

扩增完成后，代码在增广基上重新构造

$$
H_{\text{aug}},\quad S_{\text{aug}}
$$

然后再解一次广义本征值问题：

$$
H_{\text{aug}} \mathbf{u}^{(0)}_{\text{aug}}
=
E^{(0)}_{\text{aug}} S_{\text{aug}} \mathbf{u}^{(0)}_{\text{aug}}
$$

这一步的物理意义是：

> 现在工作空间变大了，所以要在这个更大的空间里，把同一个 free Hamiltonian 的基态重新表示一遍。

注意：

- 这一步还不是 kicked state
- 它仍然是 free Hamiltonian $H_0$ 的基态

真正的 kicked state 只有在下一步施加 $\hat U_{\text{kick}}$ 之后才出现。

---

## 8. kick 后的态，具体怎么投影回当前 basis

这是当前主链路最核心的一步。

### 8.1 exact kicked state

设 kick 前当前态是

$$
|\Psi\rangle = \sum_j u_j |\phi_j\rangle
$$

exact kick 后的态是

$$
|\Psi'_{\text{ex}}\rangle = \hat U_{\text{kick}} |\Psi\rangle
$$

一般来说，这个态不再严格落在当前有限 basis 子空间里。

### 8.2 我们实际保留的态

我们要找一个落在当前 basis 子空间中的近似态

$$
|\Psi'_{\text{pr}}\rangle = \sum_j u'_j |\phi_j\rangle
$$

去近似 $|\Psi'_{\text{ex}}\rangle$。

当前代码采用 Galerkin 条件：

$$
\langle \phi_i | \Psi'_{\text{pr}} \rangle
=
\langle \phi_i | \Psi'_{\text{ex}} \rangle,
\qquad \forall i
$$

把左边展开：

$$
\langle \phi_i | \Psi'_{\text{pr}} \rangle
=
\sum_j \langle \phi_i | \phi_j \rangle u'_j
=
\sum_j S_{ij} u'_j
$$

把右边展开：

$$
\langle \phi_i | \Psi'_{\text{ex}} \rangle
=
\sum_j \langle \phi_i | \hat U_{\text{kick}} | \phi_j \rangle u_j
=
\sum_j K_{ij} u_j
$$

这里定义 kick 矩阵：

$$
K_{ij} = \langle \phi_i | \hat U_{\text{kick}} | \phi_j \rangle
$$

于是得到线性方程组：

$$
S \mathbf{u}' = K \mathbf{u}
$$

如果 $S$ 可逆，形式上就是

$$
\mathbf{u}' = S^{-1} K \mathbf{u}
$$

当前代码里，对应 `apply_analytic_kick(...)`。

### 8.3 `apply_analytic_kick` 实际做的步骤

按顺序就是：

1. 构造 kick 矩阵 $K$
2. 构造 overlap 矩阵 $S$
3. 取当前系数向量 $\mathbf{u}$
4. 解

$$
S \mathbf{u}_{\text{new}} = K \mathbf{u}
$$

5. 再按 norm 做一次重归一化

---

## 9. kick 矩阵 $K$ 是怎么解析算出来的

代码没有直接数值积分 $\exp[-i\kappa\cos(2k_Lx)]$，而是先做 Bessel 展开。

对单粒子第 $a$ 个坐标：

$$
e^{-i\kappa \cos(2k_Lx_a)}
=
\sum_n (-i)^n J_n(\kappa) e^{i\,2n k_L x_a}
$$

所以问题转化成：如何计算

$$
\langle \phi_i | e^{i\,2n k_L x_a} | \phi_j \rangle
$$

对 Gaussian 基，这个矩阵元可以解析写成：

$$
\langle \phi_i | e^{i\,2n k_L x_a} | \phi_j \rangle
=
M_G \,
\exp\!\Big[
- n^2 k_L^2 (K^{-1})_{aa}
+ i\,2n k_L \mu_a
\Big]
$$

这里：

- $M_G$ 是普通 overlap 的 Gaussian prefactor
- $K$ 和 $\mu$ 来自 pair cache

因为总 kick 是各粒子单体 kick 的乘积，所以总 kernel 是各粒子 kernel 的乘积。

然后对 Bessel 级数做截断求和（当前 `n_bessel = 20`），得到 $K_{ij}$。

---

## 10. 为什么投影后会出误差

定义

$$
|\delta\Psi\rangle = |\Psi'_{\text{pr}}\rangle - |\Psi'_{\text{ex}}\rangle
$$

对任意可观测量 $A$，误差是：

$$
\Delta\langle A\rangle
=
\langle \Psi'_{\text{pr}}|A|\Psi'_{\text{pr}}\rangle
-
\langle \Psi'_{\text{ex}}|A|\Psi'_{\text{ex}}\rangle
$$

展开得到：

$$
\Delta\langle A\rangle
=
2\,\mathrm{Re}\,\langle \delta\Psi|A|\Psi'_{\text{ex}}\rangle
+
\langle \delta\Psi|A|\delta\Psi\rangle
$$

所以只要 exact kicked state 不在当前 basis 子空间里，就会有 $\delta\Psi \neq 0$，所有物理量都可能出错。

### 10.1 为什么 $\langle x^2\rangle$ 的误差特别关键

因为 exact kick 满足

$$
[U_{\text{kick}}, x^2] = 0
$$

所以 exact 情况下第一次 kick 后：

$$
\langle x^2\rangle_{\text{after}}
=
\langle x^2\rangle_{\text{before}}
$$

如果投影后 $\langle x^2\rangle$ 变了，说明：

> 投影不只是把相位结构近似了一下，而是连位置密度都扭曲了。

这正是当前 first-kick 误差里最关键的信号。

---

## 11. kick 之后的自由演化怎么做

在两次 kick 之间，当前代码不再用 TDVP 做实时间推进，而是在固定增广基里做精确传播。

系数满足：

$$
i S \frac{d\mathbf{u}}{dt} = H \mathbf{u}
$$

同样先正交化 $S$：

$$
S = V \Lambda V^\dagger,
\qquad
S^{-1/2} = V \Lambda^{-1/2} V^\dagger
$$

构造

$$
\widetilde H = S^{-1/2} H S^{-1/2}
$$

解普通本征值问题：

$$
\widetilde H \mathbf{c}_n = E_n \mathbf{c}_n
$$

于是

$$
\mathbf{u}(t) = \sum_n a_n e^{-iE_n t}\mathbf{v}_n
$$

就给出固定基里的精确自由演化。

当前代码里对应：

`free_evolve_fixed_basis(...)`

---

## 12. 把整条主链路压成一句话

从虚时间基态开始，当前 kicked 主链路就是：

$$
(\mathcal{B}_{\text{gs}}, \mathbf{u}^{(0)})
\quad \text{来自 imaginary-time / generalized eigenproblem}
$$

$$
\Downarrow
$$

$$
\text{生成动量候选} + \text{按 cond}(S)\text{ 贪心筛选}
\quad \Longrightarrow \quad \mathcal{B}_{\text{aug}}
$$

$$
\Downarrow
$$

$$
H_{\text{aug}} \mathbf{u}^{(0)}_{\text{aug}}
=
E^{(0)}_{\text{aug}} S_{\text{aug}} \mathbf{u}^{(0)}_{\text{aug}}
$$

$$
\Downarrow
$$

$$
\text{kick step: } S \mathbf{u}' = K \mathbf{u}
$$

$$
\Downarrow
$$

$$
\text{free evolution: } iS \dot{\mathbf{u}} = H \mathbf{u}
$$

---

## 13. 当前真正卡在哪里

当前证据最强的结论是：

> 问题首先不在后续 free evolution，而在第一次 kick 后那一瞬间，当前 basis 就已经不能很好地承接 exact kicked state。

也就是说：

- 不是“演化很多步以后才坏掉”
- 而是“kick 后投影这一步就已经有系统误差”

因此当前最关键的问题是：

$$
\text{如何构造一个对 } U_{\text{kick}}\Psi_0 \text{ 更近似封闭的表示空间}
$$

而不是先去优化长期传播器。
