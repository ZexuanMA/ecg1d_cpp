// ============================================================================
// main_simple.cpp — ECG1D 最小学习版
// ----------------------------------------------------------------------------
// 这个文件只做两件事（给 C++ 新手读的）：
//
//   【任务 A】用 ECG 变分方法 + 虚时间演化求谐振子基态能量 (N=1)
//             期待: E → 0.5 = ℏω/2
//
//   【任务 B】用网格有限差分法精确演化被踢谐振子 (N=1)
//             初始态: H₀ = T + V_ho 的数值基态（从网格对角化拿到）
//             每周期: 自由演化 exp(-i H₀ T) + 瞬时 kick exp(-iκcos(2kL z))
//
// 编译:
//   cd build && cmake .. && make ecg1d_simple
// 运行:
//   ./build/ecg1d_simple              # 两个任务都跑
//   ./build/ecg1d_simple --task1      # 只跑 A
//   ./build/ecg1d_simple --task2      # 只跑 B
//   ./build/ecg1d_simple --task2 --n-kicks 20   # 改 kick 次数
// ============================================================================

#pragma region Includes & 类型别名
// ---------- #include ----------
// #include 就是把另一个文件的声明拷到这里。"双引号" 找本项目的文件，<尖括号> 找系统/标准库。
#include "basis_params.hpp"       // ECG 基函数的数据结构 BasisParams
#include "permutation.hpp"        // PermutationSet: 粒子置换（玻色/费米对称化用）
#include "tdvp_solver.hpp"        // TDVP 求解器: evolution(), HamiltonianTerms, AlphaIndex
#include "svm.hpp"                // build_HS, lowest_energy, set_u_from_eigenvector
#include "kick_operator.hpp"      // apply_analytic_kick: 解析 kick 投影
#include "hamiltonian.hpp"        // overlap, kinetic_energy_functional, Harmonic_functional
#include "physical_constants.hpp" // ℏ, m, ω, κ, k_L 等物理常数

#include <Eigen/Dense>            // Eigen 线性代数库: MatrixXd, VectorXcd, 对角化...
#include <iostream>               // std::cout 输出
#include <iomanip>                // std::setprecision 格式化输出
#include <sstream>                // std::ostringstream 字符串流（格式化复数用）
#include <complex>                // std::complex 复数
#include <vector>                 // std::vector 可变长数组
#include <string>                 // std::string 字符串
#include <cmath>                  // std::cos, std::exp, std::sqrt...

// using namespace ecg1d; 让我们可以直接写 BasisParams 而不用写 ecg1d::BasisParams
using namespace ecg1d;

// 短名: Cd 代表 std::complex<double> (双精度复数)
//       在 types.hpp 里已经定义，这里重申一下方便阅读
using Cd = std::complex<double>;
#pragma endregion


#pragma region 辅助函数 (观测量 & 基函数打印)
// ============================================================================
// 辅助函数：打印当前态的能量和 ⟨x²⟩、⟨p²⟩
// ============================================================================
// 物理: 对 H = T + V_ho = p²/(2m) + (1/2)mω²x²
//   ⟨T⟩     = kinetic_energy_functional / overlap
//   ⟨V_ho⟩  = Harmonic_functional       / overlap
//   ⟨p²⟩    = 2m · ⟨T⟩
//   ⟨x²⟩    = 2·⟨V_ho⟩ / (m·ω²)
// 在 m=ω=1 单位下:  ⟨p²⟩ = 2⟨T⟩,  ⟨x²⟩ = 2⟨V_ho⟩
// ============================================================================
struct StateObs {
    double norm;   // ⟨ψ|ψ⟩
    double E;      // ⟨T⟩ + ⟨V_ho⟩
    double x2;     // ⟨x²⟩
    double p2;     // ⟨p²⟩
};

static StateObs compute_observables(const std::vector<BasisParams>& basis) {
    Cd S = overlap(basis);
    Cd T = kinetic_energy_functional(basis);
    Cd V = Harmonic_functional(basis);

    StateObs out;
    out.norm = S.real();
    double kinetic  = (T / S).real();
    double harmonic = (V / S).real();
    out.E  = kinetic + harmonic;
    out.p2 = 2.0 * mass * kinetic;
    out.x2 = 2.0 * harmonic / (mass * omega * omega);
    return out;
}

static void print_observables(const std::string& label, const StateObs& s) {
    std::cout << "  " << label << ":\n";
    std::cout << "    ⟨ψ|ψ⟩ = " << std::setprecision(10) << s.norm << "\n";
    std::cout << "    E     = " << std::setprecision(10) << s.E    << "\n";
    std::cout << "    ⟨x²⟩  = " << std::setprecision(10) << s.x2   << "\n";
    std::cout << "    ⟨p²⟩  = " << std::setprecision(10) << s.p2   << "\n";
}


// ============================================================================
// 辅助函数：打印 ECG 基函数参数
// ============================================================================
// 基函数形式: φ_i(x) = u_i · exp[-x^T(A_i+B_i)x + 2R_i^T·B_i·x - R_i^T·B_i·R_i]
//
// 对 N=1 单粒子：每个 BasisParams 只有 1 个 u, 1 个 A, 1 个 B, 1 个 R (都是复数)
// 对 N≥2: A 是 N×N, B 对角 N×N, R 是 N 向量
//
// 打印两种信息:
//   (a) 裸参数 u, A, B, R
//   (b) 物理量: 有效高斯宽度 = real(A+B), 中心位置 = real(R), 动量相位 = 2·B·imag(R)
//       对 N=1 很直观；对 N≥2 则矩阵形式
// ============================================================================
// 把复数格式化成固定宽度 14 字符的字符串, 如 " 1.000+0.000i" 或 "-0.123-4.567i"
// 规则: 实部 7 字符（宽度含符号）, 接 "+" 或 "-", 虚部 5 字符绝对值, 接 "i"
static std::string fmt_cd(Cd z, int prec = 3) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(prec);
    oss << std::setw(7) << z.real();               // " 1.000" or "-0.123"
    oss << (z.imag() >= 0 ? "+" : "-");            // 符号
    oss << std::setw(5) << std::abs(z.imag());     // "0.000"
    oss << "i";
    return oss.str();   // 总长 7 + 1 + 5 + 1 = 14
}

// 把实数格式化成固定宽度字符串
static std::string fmt_d(double x, int width = 8, int prec = 3) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(prec) << std::setw(width) << x;
    return oss.str();
}

static void print_basis(const std::string& label,
                        const std::vector<BasisParams>& basis) {
    int N = basis[0].N();
    int K = static_cast<int>(basis.size());

    std::cout << "  " << label << "  (K=" << K << ", N=" << N << ")\n";

    if (N == 1) {
        // Unicode 框线表格: 每行一个基函数, 列为 i, u, A, B, R, Re(A+B)
        // 每个复数列宽 14 (fmt_cd 输出), 左右各 1 空格 padding, 共 16
        // 实数列（i, Re(A+B)）另外处理
        const std::string sep_top =
            "  ┌─────┬────────────────┬────────────────┬────────────────┬────────────────┬──────────┐";
        const std::string sep_mid =
            "  ├─────┼────────────────┼────────────────┼────────────────┼────────────────┼──────────┤";
        const std::string sep_bot =
            "  └─────┴────────────────┴────────────────┴────────────────┴────────────────┴──────────┘";

        std::cout << sep_top << "\n";
        std::cout << "  │  i  │"
                  << "       u        │"
                  << "       A        │"
                  << "       B        │"
                  << "       R        │"
                  << " Re(A+B)  │\n";
        std::cout << sep_mid << "\n";

        for (int i = 0; i < K; i++) {
            Cd u = basis[i].u;
            Cd A = basis[i].A(0, 0);
            Cd B = basis[i].B(0, 0);
            Cd R = basis[i].R(0);
            double width = (A + B).real();

            std::cout << "  │ " << std::setw(3) << i << " │ "
                      << fmt_cd(u) << " │ "
                      << fmt_cd(A) << " │ "
                      << fmt_cd(B) << " │ "
                      << fmt_cd(R) << " │"
                      << fmt_d(width, 9, 4) << " │\n";
        }
        std::cout << sep_bot << "\n";
    } else {
        // 多粒子: 每个基函数按块打印，矩阵形式
        for (int i = 0; i < K; i++) {
            std::cout << "    basis[" << i << "]:  u = " << fmt_cd(basis[i].u) << "\n";
            std::cout << "      A =\n" << basis[i].A << "\n";
            std::cout << "      B =\n" << basis[i].B << "\n";
            std::cout << "      R = " << basis[i].R.transpose() << "\n";
        }
    }
}


#pragma endregion


#pragma region 任务 A: 虚时间演化求谐振子基态 (ECG)
// ============================================================================
// 【任务 A】虚时间演化求谐振子基态
// ============================================================================
// 物理:
//   哈密顿量 H = Σ_a [ -1/2 ∂_{x_a}² + 1/2 x_a² ]     (单位 ℏ=m=ω=1)
//   无相互作用玻色 N 粒子基态: 所有粒子都在 n=0 能级
//   精确基态能量: E₀ = N × 0.5  (每个粒子 ℏω/2)
//     N=1: 0.5
//     N=2: 1.0
//     N=3: 1.5   ...
//
// ECG 变分基函数:
//   φ(x) = u · exp[-x^T(A+B)x + 2R^T·B·x - R^T·B·R]
//   参数: u (复)，A (N×N 复对称)，B (N×N 对角复)，R (N 维复向量)
//   A 的 off-diagonal 编码 "粒子间关联"，对 N≥2 很关键
//   R 的实部是高斯中心，虚部是动量相位
//
// 多基函数 (K 个):
//   ψ(x) = Σ_{i=1}^K u_i · φ_i(x)
//   K 越多，变分流形越大，能逼近更复杂的波函数
//
// 玻色对称化:
//   对 N≥2，要对所有 N! 个粒子置换求和(+1 号玻色)
//   PermutationSet::generate(N) 给出 N! 个置换矩阵 + 符号
//
// 虚时间演化 (imaginary-time TDVP):
//   τ = it 替换后 ∂_τψ = -Hψ，高能态指数衰减 → 剩基态
//   TDVP 把这过程投影到变分流形上，梯度下降式地调 u, A, B, R
// ============================================================================
static void task1_harmonic_ground_state(int N, int K) {
    std::cout << "\n=====================================================\n";
    std::cout << " 任务 A: 纯谐振子虚时间演化求基态 (N=" << N << ", K=" << K << ")\n";
    std::cout << "=====================================================\n";

    const double E_exact = 0.5 * N;  // 玻色基态: 所有粒子都在 n=0

    // ---------------- 1. 用 SVM 方法建立一组好的初始 ECG 基 ----------------
    // SVM (Stochastic Variational Method) = 随机变分法:
    //   逐个添加基函数，每次从 n_trials 个随机候选里挑"使当前能量最低"的那一个
    //   (贪心选择)。相比纯随机初始化，SVM 给出的 basis 已经能量很低，
    //   TDVP 只需要在此基础上做微调，数值稳定性好得多。
    //
    // 这里用 K_max = K, n_trials = 2000, E_lower_bound = 0.0 (谐振子基态 ≥ 0)
    std::cout << "\n--- Phase 1: SVM 贪心建基 ---\n";
    SvmResult svm = svm_build_basis(N, /*K_max=*/K, /*n_trials=*/2000,
                                     /*terms=*/HamiltonianTerms::kinetic_harmonic(),
                                     /*seed=*/42, /*E_lower_bound=*/0.0);
    std::vector<BasisParams> basis = svm.basis;
    std::cout << "SVM 完成: basis 大小 = " << basis.size() << "\n";

    // ---------------- 2. 告诉 TDVP 要优化哪些参数 ----------------
    // AlphaIndex{a1, a2, a3, a4} 含义:
    //   a1=1 → u     参数  (a2=基函数编号, a3/a4=0)
    //   a1=2 → B 对角参数  (a2=基编号, a3=粒子编号)
    //   a1=3 → R 向量参数  (a2=基编号, a3=粒子编号)
    //   a1=4 → A 矩阵元    (a2=基编号, a3≤a4 矩阵索引)
    //
    // 四种参数全部打开: u, B, R, A 都参与 TDVP 优化
    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K; i++)                                    // u
        alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K; i++)                                    // B
        for (int j = 0; j < N; j++)
            alpha_z_list.push_back({2, i, j, 0});
    for (int i = 0; i < K; i++)                                    // R
        for (int j = 0; j < N; j++)
            alpha_z_list.push_back({3, i, j, 0});
    for (int i = 0; i < K; i++)                                    // A (对称，只存 j≤k)
        for (int j = 0; j < N; j++)
            for (int k = j; k < N; k++)
                alpha_z_list.push_back({4, i, j, k});

    std::cout << "变分参数总数: " << alpha_z_list.size()
              << "  (u:" << K
              << ", B:" << K*N
              << ", R:" << K*N
              << ", A:" << K*N*(N+1)/2 << ")\n";

    // ---------------- 3. 哈密顿量: T + V_ho ----------------
    auto terms = HamiltonianTerms::kinetic_harmonic();

    // ---------------- 4. 置换集合（玻色对称化用）----------------
    // N=1 时只有 identity 一个置换；N=2 有 2 个；N=3 有 6 个...
    // signs 全 +1 表示玻色统计（费米会 +/- 交替）
    PermutationSet perms = PermutationSet::generate(N);

    // ---------------- 5. 求解器配置 ----------------
    // SolverConfig 控制 TDVP 细节:
    //   lambda_C:  C 矩阵的 Tikhonov 正则化（防病态）
    //   resolve_u: 每步后用广义特征值重新求 u（更稳）
    //   dtao_max:  虚时间最大步长（自适应增长的上限）
    SolverConfig config;
    config.lambda_C    = 1e-8;
    config.rcond       = 1e-4;
    config.resolve_u   = true;
    config.dtao_grow   = 1.5;
    config.dtao_max    = 10.0;
    config.energy_tol  = 1e-12;
    config.optimize_A  = true;
    config.optimize_B  = true;
    config.optimize_R  = true;

    // ---------------- 6. 初始能量（SVM 之后）----------------
    // build_HS 构造 H 和 S 矩阵 (大小 K×K)
    // lowest_energy 解广义特征值问题 H·v = λ·S·v，返回最低 λ
    // 这比 compute_total_energy 更严格（后者只算当前 u 对应的 ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩）
    double E_init = lowest_energy(svm.H, svm.S);
    set_u_from_eigenvector(basis, svm.H, svm.S);  // 把 u 设成最佳线性组合

    std::cout << "\nSVM 后 E = " << std::setprecision(10) << E_init << "\n";
    std::cout << "目标 E   = " << E_exact << "  (N × ℏω/2)\n";

    // ---------------- 7. Phase 2: TDVP 虚时间微调 ----------------
    // 传入 &perms 让 evolution 用 resolve_u 路径（每步后重新求 u）
    std::cout << "\n--- Phase 2: TDVP 虚时间微调 ---\n";
    evolution(alpha_z_list, basis, /*dtao=*/1e-3,
              /*max_steps=*/2000, /*tol=*/1e-12,
              terms, config, &perms);

    // ---------------- 8. 最终能量（用最佳 u 的广义特征值结果）----------------
    auto [H_final, S_final] = build_HS(basis, perms, terms);
    double E_final = lowest_energy(H_final, S_final);

    std::cout << "\n最终 E = " << std::setprecision(12) << E_final << "\n";
    std::cout << "误差   = " << std::scientific << (E_final - E_exact)
              << std::defaultfloat << "\n";
}


#pragma endregion


#pragma region 任务 B: 网格法精确被踢谐振子 (ground truth)
// ============================================================================
// 【任务 B】网格法精确演化被踢谐振子
// ============================================================================
// 物理:
//   H(t) = H₀ + V_kick(x, t)
//   H₀   = -1/(2m) d²/dx² + 1/2 m ω² x²     ← 自由哈密顿量（动能+谐振子阱）
//   V_kick = κ cos(2 k_L x) · Σ_n δ(t - nT)  ← 周期性瞬时踢
//
// 每周期 T 的演化:
//   ψ(t+T) = U_kick · U_free(T) · ψ(t)
//   U_free(T) = exp(-i H₀ T)                 ← 自由演化
//   U_kick   = exp(-i κ cos(2 k_L x))        ← 瞬时相位旋转（位置算符对角）
//
// 数值方法:
//   (1) 把 x 离散成 N_grid 个等距格点 z_j ∈ [-z_max, +z_max]
//   (2) -d²/dx² 用 3 点有限差分 → 三对角矩阵
//   (3) H₀ 在格点上写成矩阵，对角化 → (λ_k, u_k) 本征对
//   (4) U_free(T) = U · diag(exp(-iλ_k T)) · U†  （本征基上对角）
//   (5) U_kick 在位置空间是对角，逐点乘相位即可
//
// 这是"精确"的数值方法，只受限于:
//   - 格点密度 dz（影响高动量态的精度）
//   - 盒子宽度 z_max（波函数要衰减到边界）
// ============================================================================
static void task2_grid_kicked_rotor(int n_kicks) {
    std::cout << "\n=====================================================\n";
    std::cout << " 任务 B: 网格法精确演化被踢谐振子 (N=1)\n";
    std::cout << "=====================================================\n";

    // ---------------- 参数（都从 physical_constants.hpp 读，或在此改）----------------
    // 物理常数: ℏ = m = ω = 1（谐振子单位）
    const double omega_val = omega;     // 谐振子频率
    const double mass_val  = mass;      // 粒子质量
    const double k_L_val   = k_L;       // 光晶格半波数
    const double kappa_val = kappa;     // kick 强度
    const double T_period  = 1.0;       // kick 周期（一个谐振子周期）

    // 网格参数: 格点越多越精确但越慢
    const int    N_grid = 512;
    const double z_max  = 15.0;          // 盒子 [-15, +15]（对 ω=1 基态宽度 σ=1，够了）

    std::cout << "物理参数:  ω=" << omega_val << ", m=" << mass_val
              << ", k_L=" << k_L_val << ", κ=" << kappa_val
              << ", T=" << T_period << "\n";
    std::cout << "数值参数:  N_grid=" << N_grid << ", z_max=" << z_max
              << ", n_kicks=" << n_kicks << "\n\n";

    // ---------------- 1. 建立格点 ----------------
    // dz = 两个相邻格点的距离
    // 格点位置 z[j] = -z_max + j·dz, j = 0, 1, ..., N_grid-1
    // Eigen::VectorXd 是一个实数向量
    const double dz = 2.0 * z_max / (N_grid - 1);
    Eigen::VectorXd z(N_grid);
    for (int j = 0; j < N_grid; j++) {
        z(j) = -z_max + j * dz;
    }

    // ---------------- 2. 构建 H₀ = T + V_ho 的有限差分矩阵 ----------------
    // 动能: -1/(2m) d²ψ/dx² ≈ -1/(2m) · (ψ_{j+1} - 2ψ_j + ψ_{j-1}) / dz²
    //   → 主对角 +1/(m·dz²),  次对角 -1/(2m·dz²)
    // 谐振子: V(z_j) = 1/2 · m · ω² · z_j²  → 加到主对角
    //
    // Eigen::MatrixXd::Zero(N, N) 建一个 N×N 全 0 实矩阵
    const double coeff = 1.0 / (2.0 * mass_val * dz * dz);
    Eigen::MatrixXd H0 = Eigen::MatrixXd::Zero(N_grid, N_grid);
    for (int j = 0; j < N_grid; j++) {
        // 对角项: 动能贡献 2·coeff + 谐振子势 1/2·mω²·z²
        H0(j, j) = 2.0 * coeff + 0.5 * mass_val * omega_val * omega_val * z(j) * z(j);
        // 次对角项（动能的左右邻居）
        if (j > 0)          H0(j, j - 1) = -coeff;
        if (j < N_grid - 1) H0(j, j + 1) = -coeff;
    }

    // ---------------- 3. 对角化 H₀，得到本征能量和本征态 ----------------
    // SelfAdjointEigenSolver: 专门用于实对称矩阵的对角化，比一般方法快且稳定
    // eigenvalues():  按升序排列的本征值向量 (λ_0, λ_1, ...)
    // eigenvectors(): 对应的本征向量作为列，第 k 列就是 u_k
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H0);
    Eigen::VectorXd lambda = es.eigenvalues();   // 本征值
    Eigen::MatrixXd U      = es.eigenvectors();  // 本征矢（列向量）

    std::cout << "谐振子本征能量对比（期待 E_n = n + 1/2）:\n";
    for (int k = 0; k < 5; k++) {
        std::cout << "  n=" << k
                  << "  E_grid=" << std::setprecision(10) << lambda(k)
                  << "  E_exact=" << (k + 0.5)
                  << "  误差=" << std::scientific << (lambda(k) - (k + 0.5))
                  << std::defaultfloat << "\n";
    }
    std::cout << "\n";

    // ---------------- 4. 预先算好自由演化相位 exp(-iλ_k T) ----------------
    // 这些相位永远不变，存下来每次 kick 循环直接用
    Eigen::VectorXcd free_phase(N_grid);
    for (int k = 0; k < N_grid; k++) {
        // std::exp(Cd(0, x)) = cos(x) + i·sin(x) = e^(ix)
        free_phase(k) = std::exp(Cd(0.0, -lambda(k) * T_period));
    }

    // ---------------- 5. 预先算好 kick 相位 exp(-iκ cos(2k_L z_j)) ----------------
    // 位置算符 x̂ 在格点基里是对角的（对 ψ(z_j) 就是乘以 z_j）
    // 所以 cos(2k_L x̂) 也是对角: 乘 cos(2k_L z_j)
    // 整个 kick 算符对每个格点都是一个相位因子
    Eigen::VectorXcd kick_phase(N_grid);
    for (int j = 0; j < N_grid; j++) {
        kick_phase(j) = std::exp(Cd(0.0, -kappa_val * std::cos(2.0 * k_L_val * z(j))));
    }

    // ---------------- 6. 初始波函数 = 数值基态 ----------------
    // 最低本征值对应的本征矢 u_0 就是网格上的基态
    // .cast<Cd>() 把实向量转成复向量（后续演化会乘复相位）
    Eigen::VectorXcd psi = U.col(0).cast<Cd>();
    // 归一化到 ∫|ψ|²dz = 1 的连续意义
    double norm0 = (psi.adjoint() * psi)(0).real() * dz;
    psi /= std::sqrt(norm0);

    // ---------------- 7. 打印初始态信息 ----------------
    auto compute_energy = [&](const Eigen::VectorXcd& p) -> double {
        // lambda 函数（匿名函数）: 捕获 H0 和 dz，返回 ⟨p|H₀|p⟩·dz
        return (p.adjoint() * (H0.cast<Cd>() * p))(0).real() * dz;
    };
    auto compute_norm = [&](const Eigen::VectorXcd& p) -> double {
        return (p.adjoint() * p)(0).real() * dz;
    };

    std::cout << "初始态（网格基态）:\n";
    std::cout << "  E = " << std::setprecision(10) << compute_energy(psi)
              << "  (期待 0.5)\n";
    std::cout << "  ‖ψ‖² = " << compute_norm(psi) << "\n\n";

    // ---------------- 8. 主循环: 每周期做 (自由演化 T) → (瞬时 kick) ----------------
    std::cout << std::setw(6) << "kick"
              << std::setw(18) << "E"
              << std::setw(18) << "ΔE"
              << std::setw(16) << "‖ψ‖²\n";
    std::cout << std::string(58, '-') << "\n";

    const double E_init = compute_energy(psi);
    std::cout << std::setw(6) << 0
              << std::setw(18) << std::setprecision(10) << E_init
              << std::setw(18) << 0.0
              << std::setw(16) << compute_norm(psi) << "\n";

    for (int n = 1; n <= n_kicks; n++) {
        // ---- (a) 自由演化 T: ψ ← exp(-i H₀ T) ψ = U · diag(e^{-iλT}) · U† · ψ ----
        //   U.transpose() 给出 U† (对实矩阵 U†=Uᵀ)
        //   .cast<Cd>() 保证矩阵-向量相乘时类型匹配
        Eigen::VectorXcd c = U.transpose().cast<Cd>() * psi;   // 转到本征基: c = U† ψ
        for (int k = 0; k < N_grid; k++) {
            c(k) *= free_phase(k);                              // 本征基上对角演化
        }
        psi = U.cast<Cd>() * c;                                 // 回到位置基: ψ = U · c

        // ---- (b) 瞬时 kick: ψ(z_j) ← exp(-iκ cos(2k_L z_j)) · ψ(z_j) ----
        //   位置基上对角算符，逐点乘相位
        for (int j = 0; j < N_grid; j++) {
            psi(j) *= kick_phase(j);
        }

        // ---- (c) 观测: 能量、归一化 ----
        double E_now    = compute_energy(psi);
        double norm_now = compute_norm(psi);
        std::cout << std::setw(6) << n
                  << std::setw(18) << std::setprecision(10) << E_now
                  << std::setw(18) << (E_now - E_init)
                  << std::setw(16) << norm_now << "\n";
    }

    std::cout << "\n(‖ψ‖² 应始终 ≈ 1，E 会随着 kick 次数先涨后局域化)\n";
}


#pragma endregion


#pragma region 任务 C: N=1 ECG 基态被踢一次 (解析 kick + 观测量)
// ============================================================================
// 【任务 C】N=1 谐振子基态被踢一次，计算踢后的态和观测量
// ============================================================================
// 流程:
//   Step 1: 用 SVM + TDVP 制备 N=1 谐振子基态（同任务 A，但固定 N=1）
//   Step 2: 记录踢之前的观测量: E, ⟨x²⟩, ⟨p²⟩
//   Step 3: 调用 apply_analytic_kick 做一次 kick
//           (内部是 Jacobi-Anger 展开 + u' = S⁻¹ K u)
//   Step 4: 记录踢之后的观测量，并报告 fidelity
//
// 物理预期（ψ₀(x) = π^(-1/4) exp(-x²/2)，κ=1, k_L=0.5）:
//   踢之前: E=0.5,  ⟨x²⟩=0.5,  ⟨p²⟩=0.5
//   踢之后: kick 是纯相位 → |ψ(x)|² 不变 → ⟨x²⟩ 不变 = 0.5
//          而 ⟨p²⟩ 增加了 4κ²k_L²·⟨sin²(2k_Lx)⟩
//          对 k_L=0.5: ⟨sin²(2k_Lx)⟩ = (1-e^(-1))/2 ≈ 0.316
//          所以 ⟨p²⟩_after ≈ 0.5 + 4·1·0.25·0.316 ≈ 0.816
//          E_after = (⟨p²⟩ + ⟨x²⟩)/2 ≈ (0.816 + 0.5)/2 ≈ 0.658
//
// Kick fidelity:
//   fidelity = ⟨ψ_kicked|S|ψ_kicked⟩ / ⟨ψ_原|S|ψ_原⟩
//   理论上 kick 是幺正的，fidelity = 1。若基底不够张成被踢态，fidelity < 1
//   (跟 Part 3 的 "ECG 流形对 kick 不闭合" 直接对应)
// ============================================================================
static void task3_single_kick(int K) {
    std::cout << "\n=====================================================\n";
    std::cout << " 任务 C: N=1 谐振子基态被踢一次 (K=" << K << ")\n";
    std::cout << "=====================================================\n";

    const int N = 1;

    #pragma region Step 1: 制备 N=1 基态 (SVM + TDVP)
    // ---------------- Step 1: 制备 N=1 谐振子基态（复用任务 A 的流程）----------------
    auto terms = HamiltonianTerms::kinetic_harmonic();
    PermutationSet perms = PermutationSet::generate(N);

    std::cout << "\n--- Step 1: 准备基态 (SVM + TDVP) ---\n";
    SvmResult svm = svm_build_basis(N, K, 2000, terms, 42, 0.0);
    std::vector<BasisParams> basis = svm.basis;
    set_u_from_eigenvector(basis, svm.H, svm.S);

    // 构建 alpha_z_list (u, B, R, A 全打开)
    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K; i++) alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) alpha_z_list.push_back({2, i, j, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) alpha_z_list.push_back({3, i, j, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
            for (int k = j; k < N; k++) alpha_z_list.push_back({4, i, j, k});

    SolverConfig config;
    config.lambda_C = 1e-8;
    config.rcond = 1e-4;
    config.resolve_u = true;
    config.dtao_grow = 1.5;
    config.dtao_max = 10.0;
    config.energy_tol = 1e-12;
    config.optimize_A = true;
    config.optimize_B = true;
    config.optimize_R = true;

    evolution(alpha_z_list, basis, 1e-3, 2000, 1e-12, terms, config, &perms);

    // 用最佳广义特征向量重设 u
    auto [H_gs, S_gs] = build_HS(basis, perms, terms);
    set_u_from_eigenvector(basis, H_gs, S_gs);

    #pragma endregion

    #pragma region Step 1.5: 打印基态基函数
    // ---------------- Step 1.5: 打印基态的基函数参数 ----------------
    std::cout << "\n--- Step 1.5: 基态的 ECG 基函数 ---\n";
    print_basis("基态 basis (踢之前)", basis);
    std::cout << "  注: 精确基态 exp(-x²/2) 对应 A+B = 0.5, R = 0\n";
    std::cout << "      若 K>1, 多个基函数会线性叠加逼近这个单高斯\n";
    #pragma endregion

    #pragma region Step 2: 踢之前的观测量
    // ---------------- Step 2: 踢之前的观测量 ----------------
    std::cout << "\n--- Step 2: 踢之前的观测量 ---\n";
    StateObs before = compute_observables(basis);
    print_observables("踢之前", before);
    std::cout << "  (理论期待: E=0.5, ⟨x²⟩=0.5, ⟨p²⟩=0.5)\n";
    #pragma endregion

    #pragma region Step 3: 应用 analytic kick
    // ---------------- Step 3: 应用 analytic kick ----------------
    //   U_kick = exp(-i κ cos(2 k_L x))
    //   内部通过 Jacobi-Anger 展开构造 K 矩阵，然后 u' = S⁻¹ · K · u
    //   basis 的 A, B, R 不动（高斯中心/宽度不变），只改 u
    //   这相当于把"被踢的真实态"投影到 ECG 基上
    std::cout << "\n--- Step 3: 应用一次 kick ---\n";
    std::cout << "  参数: κ = " << kappa << ", k_L = " << k_L << "\n";
    double fidelity = apply_analytic_kick(basis, perms,
                                           /*kappa=*/kappa,
                                           /*k_L=*/k_L,
                                           /*n_bessel=*/20,
                                           /*print_diag=*/true);

    #pragma endregion

    #pragma region Step 3.5: 打印 kick 后基函数
    // ---------------- Step 3.5: 打印 kick 后的基函数参数 ----------------
    std::cout << "\n--- Step 3.5: Kick 后的 ECG 基函数 ---\n";
    print_basis("kick 后 basis", basis);
    std::cout << "  注: A, B, R 不变 (analytic kick 不改这些), 只有 u 被更新\n";
    std::cout << "      u 现在带虚部, 编码了 kick 注入的动量相位\n";
    #pragma endregion

    #pragma region Step 4: 踢之后的观测量 + 理论对比
    // ---------------- Step 4: 踢之后的观测量 ----------------
    std::cout << "\n--- Step 4: 踢之后的观测量 ---\n";
    StateObs after = compute_observables(basis);
    print_observables("踢之后", after);

    // 理论值（见文件头注释）
    double x2_expected = 0.5;                              // kick 不改位置分布
    double sin2_avg    = 0.5 * (1.0 - std::exp(-4.0 * k_L * k_L));
    double p2_expected = 0.5 + 4.0 * kappa * kappa * k_L * k_L * sin2_avg;
    double E_expected  = 0.5 * (p2_expected + x2_expected);

    std::cout << "\n--- 对比理论值 ---\n";
    std::cout << std::setw(12) << "量" << std::setw(18) << "ECG"
              << std::setw(18) << "理论" << std::setw(14) << "误差\n";
    std::cout << std::string(62, '-') << "\n";
    std::cout << std::setprecision(8);
    std::cout << std::setw(12) << "E"
              << std::setw(18) << after.E
              << std::setw(18) << E_expected
              << std::setw(14) << std::scientific << (after.E - E_expected)
              << std::defaultfloat << "\n";
    std::cout << std::setw(12) << "⟨x²⟩"
              << std::setw(18) << after.x2
              << std::setw(18) << x2_expected
              << std::setw(14) << std::scientific << (after.x2 - x2_expected)
              << std::defaultfloat << "\n";
    std::cout << std::setw(12) << "⟨p²⟩"
              << std::setw(18) << after.p2
              << std::setw(18) << p2_expected
              << std::setw(14) << std::scientific << (after.p2 - p2_expected)
              << std::defaultfloat << "\n";

    std::cout << "\n--- 分析 ---\n";
    std::cout << "  Kick fidelity = " << std::setprecision(8) << fidelity
              << "  (1.0 = 基底完美表达, <1 = 基底不够)\n";
    std::cout << "  ΔE            = " << std::setprecision(6) << (after.E - before.E) << "\n";
    std::cout << "  Δ⟨x²⟩         = " << (after.x2 - before.x2)
              << "  (kick 是纯相位, 理论上 = 0)\n";
    std::cout << "  Δ⟨p²⟩         = " << (after.p2 - before.p2)
              << "  (kick 往动量里注入能量)\n";
    #pragma endregion
}
#pragma endregion


#pragma region 入口 main: 命令行解析 & 任务分发
// ============================================================================
// 入口函数
// ============================================================================
// argc = 命令行参数个数（含程序名本身）
// argv = 每个参数的字符串数组
// 例如 ./ecg1d_simple --task2 --n-kicks 20
//   argc=4, argv[0]="./ecg1d_simple", argv[1]="--task2", ...
int main(int argc, char** argv) {
    bool run_task1 = false;
    bool run_task2 = false;
    bool run_task3 = false;
    int  task1_N   = 1;   // 任务 A: 粒子数
    int  task1_K   = 1;   // 任务 A: 基函数数量
    int  n_kicks   = 10;  // 任务 B: kick 次数
    int  task3_K   = 5;   // 任务 C: 基函数数量（N=1）

    // 解析命令行
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if      (arg == "--task1")   run_task1 = true;
        else if (arg == "--task2")   run_task2 = true;
        else if (arg == "--task3")   run_task3 = true;
        else if (arg == "--N"        && i + 1 < argc) task1_N  = std::stoi(argv[++i]);
        else if (arg == "--K"        && i + 1 < argc) task1_K  = std::stoi(argv[++i]);
        else if (arg == "--K3"       && i + 1 < argc) task3_K  = std::stoi(argv[++i]);
        else if (arg == "--n-kicks"  && i + 1 < argc) n_kicks  = std::stoi(argv[++i]);
        else {
            std::cerr << "未知参数: " << arg << "\n";
            std::cerr << "用法: " << argv[0]
                      << " [--task1] [--task2] [--task3]"
                      << " [--N 粒子数] [--K 基函数数(任务A)]"
                      << " [--K3 基函数数(任务C)] [--n-kicks N]\n";
            return 1;
        }
    }

    // 不指定则全部都跑
    if (!run_task1 && !run_task2 && !run_task3) {
        run_task1 = true;
        run_task2 = true;
        run_task3 = true;
    }

    if (run_task1) task1_harmonic_ground_state(task1_N, task1_K);
    if (run_task2) task2_grid_kicked_rotor(n_kicks);
    if (run_task3) task3_single_kick(task3_K);

    return 0;
}
#pragma endregion
