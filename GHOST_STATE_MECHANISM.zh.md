# Ghost State 机制：从广义本征值问题到虚假能量

## 0. 广义本征值问题的来源

### 出发点：薛定谔方程

一切的起点是抽象 Hilbert 空间中的算符方程（薛定谔方程）：

$$
H|\psi\rangle = E|\psi\rangle
$$

这是量子力学的公设，不涉及任何基底选择，不涉及矩阵。**这是最先有的。**

能量的定义也是在这个层面：对任意态 $|\psi\rangle$，能量期望值是

$$
E = \frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle}
$$

当 $|\psi\rangle$ 恰好是本征态时，$H|\psi\rangle = E|\psi\rangle$，方差为零，期望值就是本征值。当 $|\psi\rangle$ 不是本征态时，期望值是各本征能量的加权平均，一定 $\geq E_0$（变分原理）。

### 选基底展开：两条路到同一个方程

我们无法在无穷维 Hilbert 空间里直接解方程，所以选一组有限基底 $\{|\phi_1\rangle, \dots, |\phi_K\rangle\}$，把波函数展开为 $|\psi\rangle = \sum_k u_k |\phi_k\rangle$。

然后有两条等价的路到达同一个方程：

**路线 A：投影薛定谔方程。** 把 $H|\psi\rangle = E|\psi\rangle$ 从左乘 $\langle\phi_i|$：

$$
\langle\phi_i|H|\psi\rangle = E\langle\phi_i|\psi\rangle
$$

$$
\sum_k \langle\phi_i|H|\phi_k\rangle \, u_k = E \sum_k \langle\phi_i|\phi_k\rangle \, u_k
$$

$$
\sum_k H_{ik} \, u_k = E \sum_k S_{ik} \, u_k
$$

即 $Hu = ESu$。

**路线 B：对 Rayleigh 商求极值。** 把展开代入能量定义 $E(u) = u^\dagger H u \,/\, u^\dagger S u$，令 $\partial E / \partial u^\dagger = 0$：

$$
\frac{Hu \cdot (u^\dagger S u) - Su \cdot (u^\dagger H u)}{(u^\dagger S u)^2} = 0 \quad\Rightarrow\quad Hu = E \cdot Su
$$

同一个方程。

### 三者的逻辑关系

```
H|ψ⟩ = E|ψ⟩          ← 抽象 Hilbert 空间中的薛定谔方程（公设，最先有的）
    │
    │  选一组基底展开 |ψ⟩ = Σ u_k |φ_k⟩
    ↓
┌──────────────────────────────┐
│  路线A: 投影方程              │  路线B: 能量期望值求极值
│  ⟨φ_i|H|ψ⟩ = E⟨φ_i|ψ⟩      │  ∂/∂u† (u†Hu / u†Su) = 0
└──────────────┬───────────────┘
               │  两条路到同一个终点
               ↓
           Hu = ESu              ← 广义本征值问题
               │
               │  如果基底恰好正交 (S = I)
               ↓
           Hu = Eu               ← 普通本征值问题（特例）
```

所以：

- **$Hu = Eu$ 是 $Hu = ESu$ 的特例**（当基底正交，$S = I$ 时）
- **$Hu = ESu$ 不是额外发明的**，它就是薛定谔方程在非正交基底下的自然表现
- **$E = u^\dagger Hu / u^\dagger Su$ 也不是额外定义的**，它就是 $\langle\psi|H|\psi\rangle / \langle\psi|\psi\rangle$ 在基底表示下的写法

### 本征态是 Rayleigh 商的驻点

路线 B 告诉我们 $Hu = ESu$ 是 Rayleigh 商 $R(u) = u^\dagger Hu / u^\dagger Su$ 的驻点条件。这个关系值得展开说清楚。

**正方向：本征态 → 驻点。** 如果 $u$ 是本征态（$Hu = E_n Su$），代入 Rayleigh 商：

$$
R(u) = \frac{u^\dagger H u}{u^\dagger S u} = \frac{u^\dagger (E_n S u)}{u^\dagger S u} = E_n
$$

而 $Hu = E_n Su$ 正是 $\partial R / \partial u^\dagger = 0$ 的条件，所以本征态确实是驻点，驻点值就是本征值。

**反方向：驻点 → 本征态。** 如果 $u$ 是 $R(u)$ 的驻点，则 $Hu = R(u) \cdot Su$，即 $u$ 满足广义本征值方程，$R(u)$ 就是对应的本征值。

**所以本征态和 Rayleigh 商的驻点一一对应。**

#### 几何图像

$R(u)$ 对 $u$ 的整体缩放不敏感（分子分母同阶），所以它实际上是定义在"方向空间"上的函数。可以想象它是高维球面上的一个地形：

```
R(u)
 ^
 |
E_K ──── ·                         ← 最高本征值 = 全局最大值
 |       / \
 |      /   \
E_1 ── ·     ·                     ← 第一激发态 = 鞍点
 |    / \   / \
 |   /   \ /   \
E_0 ·     ·     ·                  ← 基态 = 全局最小值
 |
 └──────────────────→ "方向"
```

- **全局最小值** = 基态能量 $E_0$，在基态方向取到
- **全局最大值** = 最高本征值 $E_K$
- **中间的鞍点** = 各激发态能量

#### 变分原理的本质

变分原理就是这个几何图像的直接推论：Rayleigh 商的全局最小值是基态能量，所以无论你选什么方向（什么线性组合），$R(u) \geq E_0$。

用本征态展开 $u = \sum_n c_n u_n$（$u_n$ 是第 $n$ 个本征态）可以直接验证：

$$
R(u) = \frac{\sum_n |c_n|^2 E_n}{\sum_n |c_n|^2} \geq \frac{E_0 \sum_n |c_n|^2}{\sum_n |c_n|^2} = E_0
$$

就是加权平均不小于最小值。等号当且仅当只有 $c_0 \neq 0$，即 $u$ 恰好是基态方向。

---

### 为什么我们的基底不正交

正交基底更方便（$S = I$，没有广义本征值问题的麻烦）。但高斯基函数 $\phi_k(x) = \exp(-x^T A_k x)$ 之间天然有非零 overlap——两个高斯函数的乘积仍是高斯函数，积分不为零。

我们**可以**先做正交化，得到正交基底后解普通 $Hu = Eu$。实际上 $S^{-1/2}$ 变换就是在做这件事——把非正交基底变换成正交基底。

所以 ghost state 问题也可以这样理解：**正交化过程（$S^{-1/2}$）本身在数值上不稳定。** 如果基函数几乎线性相关，"正交化"就要把几乎重叠的方向强行掰开——这个"掰开"的过程放大了所有数值误差。

---

## 1. 广义本征值问题的标准解法

标准解法是通过 $S^{-1/2}$ 变换消去右边的 $S$：

1. 对 $S$ 做谱分解：$S = V W V^T$（$W$ 是特征值对角矩阵，$V$ 是特征向量矩阵）
2. 构造 $S^{-1/2} = V \cdot \text{diag}(1/\sqrt{w_i}) \cdot V^T$
3. 令 $\tilde{u} = S^{1/2} u$，代入 $Hu = ESu$ 得到标准本征值问题：

$$
\tilde{H} \tilde{u} = E \tilde{u}, \quad \text{其中 } \tilde{H} = S^{-1/2} H S^{-1/2}
$$

**变分原理保证**：在精确算术下，$\tilde{H}$ 的最小本征值 $E_0 \geq E_\text{exact}$，变分能量永远是基态的上界。

下面用一个最小的 2×2 例子，展示 S 近奇异时这个保证如何被数值误差打破。

---

## 2. 设定：两个几乎相同的基函数

设两个高斯基函数 $\phi_1$ 和 $\phi_2$ 几乎相同，宽度参数仅差 $\varepsilon = 10^{-9}$。

### Overlap 矩阵

$$
S = \begin{pmatrix} 1 & 1-\varepsilon \\ 1-\varepsilon & 1 \end{pmatrix}
$$

对角元 = 1（归一化），非对角元 = $1 - \varepsilon \approx 0.999999999$，两个基函数几乎完全重叠。

### Hamiltonian 矩阵

因为 $\phi_1 \approx \phi_2$，H 矩阵元也几乎相同：

$$
H = \begin{pmatrix} h + \delta_1 & h + \delta_3 \\ h + \delta_3 & h + \delta_2 \end{pmatrix}
$$

其中 $h = 1.527$（接近我们的 Gaussian benchmark），$\delta_1, \delta_2, \delta_3$ 是 $O(\varepsilon)$ 量级的小修正。

取具体数值：$\delta_1 = 3\varepsilon$，$\delta_2 = \varepsilon$，$\delta_3 = 2\varepsilon$。

---

## 3. S 的谱分解

S 的特征值：

$$
w_1 = 1 - (1 - \varepsilon) = \varepsilon = 10^{-9}
$$

$$
w_2 = 1 + (1 - \varepsilon) = 2 - \varepsilon \approx 2
$$

对应的特征向量：

$$
v_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix} \quad \text{（差态，范数极小）}
$$

$$
v_2 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix} \quad \text{（和态，范数正常）}
$$

### 物理含义

- **和态** $\phi_1 + \phi_2$：两个几乎相同的函数相加，得到一个"加倍"的函数，范数正常（$\sim \sqrt{w_2} = \sqrt{2}$）
- **差态** $\phi_1 - \phi_2$：两个几乎相同的函数相减，几乎抵消，范数极小（$\sim \sqrt{w_1} = \sqrt{\varepsilon}$）

差态在物理上几乎不存在——它是两个基函数之间微小差异的放大。

---

## 4. 构造 $S^{-1/2}$

$$
S^{-1/2} = V \cdot \text{diag}\left(\frac{1}{\sqrt{w_1}}, \frac{1}{\sqrt{w_2}}\right) \cdot V^T
$$

代入数值：

$$
\frac{1}{\sqrt{w_1}} = \frac{1}{\sqrt{10^{-9}}} \approx 31623
$$

$$
\frac{1}{\sqrt{w_2}} = \frac{1}{\sqrt{2}} \approx 0.707
$$

**差态方向被放大了 31623 倍，和态方向基本不变。** 这是所有问题的根源。

---

## 5. 在精确算术下：一切正常

先在 S 的特征基下看 H。定义 $H_\text{rot} = V^T H V$：

$$
(H_\text{rot})_{00} = v_1^T H v_1 = \frac{1}{2}(H_{11} + H_{22} - 2H_{12})
= \frac{1}{2}(\delta_1 + \delta_2 - 2\delta_3)
= \frac{1}{2}(3\varepsilon + \varepsilon - 4\varepsilon) = 0
$$

$$
(H_\text{rot})_{11} = v_2^T H v_2 = \frac{1}{2}(H_{11} + H_{22} + 2H_{12})
= \frac{1}{2}(2h + \delta_1 + \delta_2 + 2\delta_3)
= h + 3\varepsilon \approx h
$$

$$
(H_\text{rot})_{01} = v_1^T H v_2 = \frac{1}{2}(H_{11} - H_{22})
= \frac{1}{2}(\delta_1 - \delta_2) = \varepsilon
$$

然后 $\tilde{H} = W^{-1/2} H_\text{rot} W^{-1/2}$：

$$
\tilde{H}_{00} = \frac{(H_\text{rot})_{00}}{w_1} = \frac{0}{\varepsilon} = 0
$$

$$
\tilde{H}_{11} = \frac{(H_\text{rot})_{11}}{w_2} = \frac{h + 3\varepsilon}{2 - \varepsilon} \approx \frac{h}{2}
$$

$$
\tilde{H}_{01} = \frac{(H_\text{rot})_{01}}{\sqrt{w_1 w_2}} = \frac{\varepsilon}{\sqrt{2\varepsilon}} = \sqrt{\varepsilon/2} \approx 2.2 \times 10^{-5}
$$

$\tilde{H}$ 的本征值可以精确算出来，约为 $0$ 和 $h/2 \approx 0.764$。回代后对应的能量分别约为 $h = 1.527$ 和某个激发态能量。

**精确算术下，变分能量确实 $\geq E_\text{exact}$，没有问题。**

---

## 6. 在浮点算术下：噪声被放大

现在考虑数值计算的现实。H 矩阵元的计算涉及积分、置换求和、复数运算，每个元素有 $O(\epsilon_\text{mach})$ 的舍入误差。设：

$$
H^\text{计算} = H^\text{精确} + \Delta H
$$

其中 $|\Delta H_{ij}| \sim \epsilon_\text{mach} \cdot |h| \approx 10^{-16} \times 1.5 \approx 1.5 \times 10^{-16}$。

### 关键步骤：$(H_\text{rot})_{00}$ 的误差

$$(H_\text{rot})_{00} = \frac{1}{2}(H_{11} + H_{22} - 2H_{12})$$

精确值为 $0$（在我们的例子中）。但数值计算时：

$$(H_\text{rot})_{00}^\text{计算} = 0 + \frac{1}{2}(\Delta H_{11} + \Delta H_{22} - 2\Delta H_{12})$$

误差量级：$\sim \frac{1}{2} \times 4 \times 1.5 \times 10^{-16} \approx 3 \times 10^{-16}$

**这是一个"大数减大数"的结果**——三个 $\approx 1.527$ 的数相加减，得到一个 $10^{-16}$ 量级的数。这就是 catastrophic cancellation。

### 致命一步：除以 $w_1$

$$
\tilde{H}_{00}^\text{计算} = \frac{(H_\text{rot})_{00}^\text{计算}}{w_1} = \frac{3 \times 10^{-16}}{10^{-9}} = 3 \times 10^{-7}
$$

精确值是 $0$，但数值结果是 $3 \times 10^{-7}$。

**这个 $3 \times 10^{-7}$ 的噪声直接出现在 $\tilde{H}$ 的对角元上。** 它会让 $\tilde{H}$ 的最小本征值偏移约 $\pm 3 \times 10^{-7}$。

如果偏移方向恰好是负的，回代后的能量就变成了 $1.527 - \text{something}$，低于真实基态。

---

## 7. 更真实的场景：误差来源不止机器精度

上面用的是最理想的情况：H 矩阵元精度达到 $10^{-16}$。实际中还有更大的误差来源：

### 7.1 取实部操作

我们的代码在 `lowest_energy_full` 中做：

```cpp
Eigen::MatrixXd Ss = (0.5 * (S + S.adjoint())).real();
```

如果 S 有 $O(|B|)$ 或 $O(|R|)$ 量级的虚部（因为基函数有复数参数），取实部操作引入了 $O(|B|)$ 的误差。当 $|B| \sim 10^{-3}$ 时：

$$
\tilde{H}_{00}^\text{误差} \sim \frac{10^{-3}}{10^{-9}} = 10^{6}
$$

**放大了一百万倍。** 这才是真正致命的。

### 7.2 多个近奇异方向的叠加

K=20 个基函数时，S 可能有多个较小的特征值。即使每个方向的误差不大，多个方向的误差叠加后效果更强。

### 7.3 矩阵元计算本身的精度

置换求和涉及多项复数加减，每一步都有舍入。最终矩阵元的有效精度可能远低于 $10^{-16}$。

---

## 8. 更一般的量级估计

设 S 的最小保留特征值为 $w_\min$，H 矩阵元的有效误差为 $\delta_H$。那么：

$$
\tilde{H} \text{ 中的噪声量级} \sim \frac{\delta_H}{w_\min}
$$

$$
\text{本征值的噪声量级} \sim \frac{\delta_H}{w_\min}
$$

| $w_\min$ | $\delta_H$ | 噪声量级 | 能否产生 ghost |
|-----------|------------|----------|--------------|
| $10^{-4}$ | $10^{-16}$ | $10^{-12}$ | 不能 |
| $10^{-8}$ | $10^{-16}$ | $10^{-8}$ | 累积后可能 |
| $10^{-8}$ | $10^{-12}$ | $10^{-4}$ | **很容易** |
| $10^{-8}$ | $10^{-6}$ | $10^{2}$ | **灾难性的** |
| $10^{-4}$ | $10^{-6}$ | $10^{-2}$ | 可能 |

我们的情况大约是：$w_\min \sim 10^{-8}$（截断阈值），$\delta_H \sim 10^{-12}$（复数取实部 + 置换求和的累积误差），所以噪声 $\sim 10^{-4}$。

这完全足以在 30 轮 refine 中累积出 $O(10^{-2})$ 到 $O(10^{-1})$ 的能量漂移。

---

## 9. 为什么漂移总是向下（棘轮效应）

上面的噪声是双向的——$\tilde{H}$ 的噪声可以让本征值偏高或偏低。但 stochastic refine 的接受准则是：

```
if (E_trial < E_current) → 接受
else                      → 拒绝
```

这是一个不对称滤波器：

- 噪声让 E 偏高 → 不满足 $E_\text{trial} < E_\text{current}$ → 被拒绝
- 噪声让 E 偏低 → 满足 $E_\text{trial} < E_\text{current}$ → 被接受

**系统性地选择了向下的噪声，忽略了向上的噪声。**

每次替换引入的 bias 量级约为噪声的一半（因为只保留了负偏差那一半分布）。30 轮 × K 次替换 × 每次 $\sim 10^{-5}$ 到 $10^{-4}$ 的 bias：

$$
\text{总漂移} \sim 30 \times 20 \times 5 \times 10^{-5} \approx 3 \times 10^{-2}
$$

这与观察到的"从 1.527 漂移到约 1.5"完全一致。

---

## 10. 为什么现有的五层保护都不够

### 10.1 $E_\text{lower\_bound}$

物理下界通常设得很保守（比如 0 或 1.0），不会拒绝"看起来合理但偏低 0.03"的能量。

### 10.2 $\text{max\_round\_drop} = 10^{-4}$

限制每轮的总降幅。但每轮内 K 个替换各降一点点（$< 10^{-4}/K$），累积到 $10^{-4}$ 仍然被允许。30 轮就是 $30 \times 10^{-4} = 3 \times 10^{-3}$。

而且这个阈值设大了不起作用，设小了又会阻止合法的能量改善。

### 10.3 Overlap 检查（阈值 0.95）

归一化 overlap = $|S_{ij}| / \sqrt{S_{ii} S_{jj}}$。

问题是：两个基函数可以"pair-wise overlap 都 < 0.95"，但 K 个基函数**合在一起**仍然导致 S 近奇异。这是因为线性相关可以是**多体效应**——三个基函数中任意两个都不太像，但三个合在一起张成一个几乎二维的空间。

Pair-wise 检查抓不住这种多体线性相关。

### 10.4 特征值截断（$\text{rcond} = 10^{-8}$）

截断确实丢掉了最病态的方向。但刚好在截断边界以上的方向仍然被保留，而它们的放大因子仍然很大（$1/\sqrt{10^{-8}} \approx 10^4$）。

而且截断本身改变了变分空间的维度。不同的 trial 可能导致不同的截断维度，使得能量之间的比较不在同一个空间中进行——这本身就可以引入 bias。

### 10.5 Rayleigh 商交叉验证

验证的是：$u^T H u / u^T S u \approx E_0$？

但 $u = S^{-1/2} \tilde{u}$，而 $S^{-1/2}$ 本身就是有噪声的。$u$ 在近零特征值方向上有大的 spurious 分量。计算 $u^T H u / u^T S u$ 时，分子和分母中这些大分量**精确配合**（因为数学上就是同一个问题的重写），给出和 $E_0$ 一致的结果。

**这是用同一把有偏差的尺子量了两次，两次结果当然一致——但不代表尺子是准的。**

---

## 11. Varga 为什么没被困住

Varga/Suzuki 书中使用的设定与我们有三个关键区别：

| | Varga | 我们 |
|---|---|---|
| 基函数参数 | 实对称正定 A | 复数 A + iB + iR |
| 非线性参数数（N=2） | 1 个（$\alpha_{12}$） | ~8 个 |
| S 矩阵元 | 实正数 | 复数（取实部处理） |

**结果**：

1. 他们的参数空间小 → 随机采样不容易撞上近线性相关的配置
2. 实正定 A 的约束很强 → S 天然条件更好
3. S 矩阵元都是正实数 → 没有"取实部"引入的额外误差
4. 简单的 overlap 阈值检查在他们那里就足够了

他们在书 p.48-50 确实提到了线性相关问题，说明他们也遇到过。但在实参数 + 小参数空间的设定下，这是"轻症"，简单处方就能治。

我们的复数参数化把这个问题变成了"重症"。

---

## 12. 小结

Ghost state 的完整因果链：

```
φ_i ≈ φ_j （基函数近线性相关）
    ↓
S 有一个很小的特征值 w ≈ ε
    ↓
S^{-1/2} 在差态方向上放大 1/√ε
    ↓
H 矩阵元在差态方向上是"大数减大数"（catastrophic cancellation）
    ↓
放大后的噪声 ~ δ_H / ε 出现在 H̃ 中
    ↓
H̃ 的本征值有 O(δ_H / ε) 的随机误差
    ↓
误差可正可负，但贪心接受准则只选负的（棘轮）
    ↓
多轮累积 → 能量系统性地低于真实基态
```
