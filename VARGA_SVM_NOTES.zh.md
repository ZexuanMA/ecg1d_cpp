# Varga/Suzuki SVM 参考笔记

## 目的

这份笔记整理了我基于以下参考资料对 `Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems`（Suzuki, Varga, 1998）所做的快速学习，重点回答两个问题：

1. 这本书的方法论对我们当前 `2-particle gaussian` benchmark 有什么直接启发？
2. 这些启发哪些可以纳入主线实现，哪些只能作为 `N=2` 的特化技巧，不能直接外推到以后 `N=10` 的代码框架？

## 说明

- 这里不是完整读书笔记。
- 目前我能直接访问到的内容主要是出版社/Google Books 的目录与简介、Varga/Suzuki 的原始论文摘要、以及后续综述的公开摘要或检索片段。
- 因此下面分成两类：
  - `文献直接支持`：能从可访问资料中直接确认的结论。
  - `结合本项目的推断`：基于文献和我们当前 benchmark 现象做的工程判断，不等同于书中的原话。

## 文献直接支持

### 1. 这本书的核心方法就是 correlated Gaussian + stochastic sampling

Springer 和 Google Books 对本书的简介都写得很明确：

- 方法主体是用 correlated Gaussian 展开波函数。
- 非线性参数通过 stochastic sampling 优化。
- 目标是 few-body bound-state problem 的高精度变分解。

这和我们当前程序的整体思路是一致的，不是旁门路线。

### 2. 这套方法本来就是围绕 few-body 而不是 generic many-body 搭建的

书的定位是 `few-body problems`，目录里也能看到它把大量篇幅放在：

- `Stochastic variational method`
- `Variational trial functions`
- `Matrix elements for spherical Gaussians`

原始论文和 1998 年综述也都把适用对象描述为 `N=2–7` 或 `N=2–8 and possibly more-particle problems`。

这说明：

- 它不是“对任意多体都自然扩展”的黑盒方案。
- 它对 `N=10` 不是不能想，但也绝对不能按“把 `N=2` 代码线性放大五倍”来理解。

### 3. correlated Gaussian 的优势在于矩阵元可解析、相关性表达自然

Varga/Suzuki 的论文摘要和后续综述都强调两点：

- 相关 Gaussian 使矩阵元计算能保持解析或半解析，数值效率高。
- 它对强相关 few-body 体系尤其自然。

这点对我们很重要，因为它说明：

- `gaussian interaction` 这类平滑势，本来就是 correlated Gaussian 应该擅长的对象。
- 如果 `gaussian` benchmark 仍然不够准，第一怀疑对象应当是优化流程、采样分布、广义本征值稳定性，而不是 ansatz 本身完全不对路。

### 4. SVM 的典型流程是逐个加入 basis function，并按能量准则筛选

公开综述和后续 review 对参数选择的描述非常一致：

- 新 basis function 通常是一个一个加入。
- 非线性参数的选择原则非常直接：能量越低越好。
- 优化是 stochastic trial-and-error，而不是一开始就做全局高维连续优化。

这和我们现在的 `svm_build_basis()` 大方向一致，但书/综述暗示得更强的一点是：

- **前期 basis construction 本身就是主战场**。
- TDVP 或其他连续优化更像后期抛光，而不是替代前期的 basis search。

### 5. 书和相关文献明确重视坐标选择、CM 分离、以及 trial function 设计

目录里单独有 `Variational trial functions` 和 `Matrix elements for spherical Gaussians`，后续相关 ECG 文献也一再强调：

- Jacobi/relative coordinate 的组织方式很关键。
- 对总角动量、宇称、平移不变性等约束，trial function 的写法会直接影响收敛速度和稳定性。

这说明“物理启发式地限制 trial space”在这个流派里是正常操作，不是投机取巧。

### 6. Zaklama 等人的 1D ECG 文章给了一个非常直接的 1D 解析矩阵元框架

你刚上传的文章
`Matrix Elements of One Dimensional Explicitly Correlated Gaussian Basis Functions`
比前面的书更贴近我们当前项目，因为它就是专门写 1D ECG 矩阵元的。

文中直接给出的几个关键点：

- 1D ECG 基函数写成
  `Phi_A^m(x) = (prod_k x_k^{m_k}) exp(-1/2 x^T A x)`。
- `A` 必须是对称正定矩阵，这样基函数才可归一。
- overlap、kinetic、one-body、two-body matrix elements 都可以通过生成函数系统推出来。
- 对 two-body potential `V(x_i - x_j)`，作者把它统一写成 `delta(w^(ij) x - x)` 的形式处理，其中 `w^(ij)` 就是选择粒子差坐标的权向量。
- 最终 two-body potential matrix element 会化成一个一维积分
  `∫ V(x) x^(...) exp(-c x^2 / 2) dx`，
  对很多势函数都能解析或高精度数值处理。

这和我们代码当前做的事情在结构上是对得上的：我们也在处理 1D、多粒子、pair potential，并且势函数也是 `x_a-x_b` 的函数。

### 7. 这篇文章说明“1D ECG 本身并不要求改坐标体系”，关键是如何组织 pair coordinate

这篇文章的推导是在实验室坐标 `x = (x1, ..., xN)` 下完成的，不是先整体改写成 Jacobi 坐标再做。

这点很关键，因为它直接支持一个判断：

- 我们**不需要为了 1D pair potential 的解析矩阵元而把主代码改成 CM/relative 专用坐标体系**。
- 在原始坐标中，通过 `w^(ij)` 这种“选取粒子差”的线性形式，同样可以把两体势矩阵元系统地做出来。

对我们当前项目来说，这进一步支持：

- 主框架可以继续保留原始坐标参数化。
- 真正需要强化的，是 basis family 和数值稳定性，而不是先重写坐标体系。

### 8. 文章里的 polynomial prefactor 值得认真考虑

这篇文章和我们当前实现最大的结构差异之一，是它的基函数允许乘上多项式前因子 `prod_k x_k^{m_k}`。

这件事的重要性在于：

- 作者明确说 polynomial prefactor 不仅能提高精度，还能支持更多算符矩阵元。
- 这意味着在 1D ECG 里，**单纯的“纯 Gaussian 核”并不是唯一自然的基函数选择**。

对我们项目的启发是：

- 如果将来 `gaussian` benchmark 仍然卡精度，除了继续加强 SVM 之外，一个正统方向是考虑更广的 ECG family，而不是只靠当前这套纯指数核去硬压误差。
- 对 `delta` 这类有 cusp 的问题，polynomial prefactor 虽然未必直接解决 cusp，但至少说明“扩充 basis family”是这个流派里完全正常的做法。

### 9. 这篇文章的基函数形式和我们当前实现不完全一样

这点必须单独写清楚。

文中 1D ECG 的主形式是：

- 对称正定 `A`
- 无显式位移中心
- 用多项式前因子 `m`

而我们当前代码里实际在用的是：

- `A`
- 额外的 `B`
- 位移/中心参数 `R`
- 没有多项式前因子

所以：

- 这篇文章**不是**我们现有公式的逐行对照版。
- 但它给出的核心思想仍然高度相关：
  - 1D ECG 完全可以在原始坐标中系统推导；
  - pair potential 的自然变量是 `x_i-x_j`；
  - 基函数家族可以比“纯 Gaussian without prefactor”更丰富；
  - 正定性约束是基础，不是可选项。

## 对当前项目新增的具体启发

### 1. `gaussian` benchmark 不一定该优先加“位移”，更应优先想 pair-coordinate-aware basis

从 Zaklama 这篇文章的表达方式看，1D 两体势的自然结构是 `w^(ij) x = x_i - x_j`。

这给出的启发是：

- 对 `gaussian interaction`，更重要的是把基函数在“相对坐标方向”的表达能力做强。
- 相比盲目增加一般性的 `R` 位移，自觉地增强与 pair coordinate 相关的宽度/相关矩阵采样，物理上更像正路。

这不等于必须改坐标求解，只是说明 sampling 应该更有方向感。

### 2. 现在的 `A/B/R` 参数化可能不是最标准的 1D ECG 形式

文章里的 1D ECG 形式比我们现在的写法更“教科书”一些。

这不说明我们现在错了，但说明：

- 我们当前的 `A/B/R` 分解更像是某种项目内的特定参数化。
- 以后如果要和 ECG 文献系统对齐，应该明确回答：
  - `A/B/R` 与文献中标准 ECG 参数之间的映射是什么？
  - 我们现在的 basis family 覆盖了文献里的哪一类函数？
  - 哪些文献中常见的自由度，我们现在实际上没有？

这个问题对 long term 很重要，因为它决定我们将来到底是在“优化一个标准 ECG 实现”，还是在“优化一个项目自定义的 ECG-like ansatz”。

### 3. polynomial-prefactor ECG 可能是后续比“坐标特化”更值得研究的扩展方向

如果以后要进一步提升静态精度，我现在更倾向于优先考虑：

1. 修数值稳定性；
2. 改善 SVM 采样；
3. 评估是否需要引入更丰富的 ECG family，例如 polynomial prefactor；

而不是优先考虑：

1. 为 `N=2 gaussian` 单独改成 CM/relative 主求解器。

## 结合本项目的推断

### 1. 不应该把“两体 gaussian 的 CM/relative 特化”当成主线架构

这是我现在最明确的判断。

原因：

- 对 `N=2` 且相互作用只依赖 `x1-x2` 的情形，CM/relative 分离当然成立。
- 但书和论文的主线方法不是“把所有问题都先拆成一堆独立坐标再解”，而是用 correlated Gaussian 在一般的 few-body 相对坐标空间里处理。
- 到 `N>2` 时，最多只是全局 CM 能分离，内部自由度仍是 `(N-1)` 维耦合问题，不会退化成一串独立的一维 relative 方程。

所以：

- `N=2 gaussian` 的 CM 特化可以作为 benchmark 辅助理解。
- 但不能把它升级成我们以后 `N=10` 主实现的基础设计。

### 2. 当前 `gaussian` 的主要问题更像是优化/接受机制，而不是表示能力

文献层面，correlated Gaussian 正是为这类平滑 few-body 势服务的。

结合我们当前现象：

- `delta` 停在 `1e-3`，很像 cusp/contact 结构对 plain Gaussian 的先天不友好。
- `gaussian` 只差到 `1e-6 ~ 1e-5`，更像 basis search 还不够强，或者 generalized eigenproblem / refine 的数值稳定性有缺陷。

我的工程判断是：

- `gaussian` 不该优先怀疑“坐标体系错了”。
- 应该优先怀疑：
  - basis 候选分布不够物理
  - overlap 近线性相关处理不稳
  - stochastic refine 接受了 ghost state
  - TDVP 在一个质量一般的 basis 上做了过多本不该承担的工作

### 3. 书里最值得迁移进我们代码主线的，不是“两体特化”，而是更标准的 SVM discipline

结合这本书和 Zaklama 这篇 1D ECG 文章，我认为最值得吸收的是下面几条：

1. **把 basis construction 当主战场。**
   `gaussian` 想做准，优先应该提高 Phase 1/1.5 的质量，而不是指望 TDVP 从较差 basis 上“救回来”。

2. **严格处理 overcompleteness / near-linear dependence。**
   这在 non-orthogonal Gaussian basis 里是基本问题，不是边角问题。我们现在 `refine` 出 ghost energy，方向上和这个问题高度一致。

3. **对 trial space 做物理启发式约束，但保持主框架通用。**
   例如：
   - 对 bosonic ground state，更偏向交换对称的候选；
   - 对 trap + smooth pair potential，更偏向加强 `x_i-x_j` 相关方向的候选；
   - 对 trap ground state，更偏向中心在原点、较少无意义位移的候选；
   - 但仍然保留原始坐标/一般 `A` 矩阵参数化，不把它改成 `N=2` 专用求解器。

4. **把局部 refinement 当作 basis polish，而不是容错垃圾桶。**
   文献主线是先把 basis 挑对，再做连续参数优化。我们当前的 `TDVP` 更像背了太多本该由 SVM 完成的工作。

### 4. 对 `N=10` 要有更保守的预期

文献摘要里说到 `N=2–8 and possibly more`，这个表述很微妙：

- 它表示方法不是只会做两三个粒子。
- 但它也没有承诺“完全 unrestricted 的 `N=10` 仍然轻松可做”。

对我们项目的含义是：

- `N=10` 不是不能想。
- 但很可能需要明显更强的结构化假设，例如：
  - 更强的对称性限制
  - 更精细的 basis family 分层
  - cluster / subsystem 思路
  - 更严格的 basis pruning 和数值稳定控制

不能把它理解成“把现有 `N=2` SVM/TDVP 直接推广到十粒子，只是算得更慢”。

## 对当前代码的具体启发

### 短期应优先做的

1. 修 `stochastic_refine()` 的 ghost-state 接受问题。
2. 把 generalized eigenproblem 的稳健性单独抽出来做诊断。
3. 把 `gaussian` benchmark 的 basis sampling 做得更物理，但仍留在当前通用参数化下。
4. 用更强的 Phase 1/1.5 去验证：`gaussian` 能否在不改坐标体系的前提下显著优于当前 `4e-6`。

### 不建议现在做成主线的

1. 把求解器整体改写成两体 CM/relative 坐标专用版。
2. 为了 `gaussian` benchmark 精度，牺牲未来 `N>2` 的统一参数化。
3. 先不修数值稳定性，就直接增加 `K`、开大 refine、或者把 TDVP 跑得更猛。

## 当前工作结论

一句话总结：

- **Varga/Suzuki 这条路本身是对的。**
- **但它支持的是“更规范的 correlated Gaussian + SVM few-body 框架”，不是“把两体 gaussian 特化成主算法”。**
- **对我们现在最直接的帮助，是修 basis search / overcompleteness / generalized eigenproblem 稳定性，而不是换整个坐标体系。**

## 参考链接

### 书

- Springer 书页：
  https://link.springer.com/book/10.1007/3-540-49541-X
- Google Books：
  https://books.google.com/books?id=SxnyCAAAQBAJ

### 原始论文与综述

- Varga, Suzuki, `Precise solution of few-body problems with the stochastic variational method on a correlated Gaussian basis`, Phys. Rev. C 52, 2885 (1995)：
  https://doi.org/10.1103/PhysRevC.52.2885
- Varga, Suzuki, `Stochastic variational method with a correlated Gaussian basis`, Phys. Rev. A 53, 1907 (1996)：
  https://doi.org/10.1103/PhysRevA.53.1907
- Suzuki, Varga, Usukura, `Stochastic variational methods in few-body systems`, Nucl. Phys. A 631, 91–110 (1998)：
  https://doi.org/10.1016/S0375-9474(98)00017-7
- Zaklama, Zhang, Rowan, Schatzki, Suzuki, Varga, `Matrix Elements of One Dimensional Explicitly Correlated Gaussian Basis Functions`, Few-Body Syst. 61, 6 (2020)：
  https://doi.org/10.1007/s00601-019-1539-3

### 后续综述/相关资料

- Pre-Born–Oppenheimer molecular structure theory（公开综述片段里明确提到 SVM 的 one-by-one basis selection 思路）：
  https://doi.org/10.1080/00268976.2018.1530461
- Few-body problem in terms of correlated Gaussians（说明 Suzuki/Varga 的书对 SVM+ECG 已经讲得很系统）：
  https://doi.org/10.1103/PhysRevE.76.046702
