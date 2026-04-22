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
#include "target_fitting.hpp"     // Route 1: target-fitting TDVP
#include "realtime_tdvp.hpp"      // 实时 TDVP 动基底演化

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
    std::cout << "  注: A, B, R 看起来不变, 但这只是代码实现的选择 ——\n";
    std::cout << "      物理上 kick 应该把 R 虚部按 Jacobi-Anger 展开移位:\n";
    std::cout << "        U_kick·φ_R = Σ_n (-i)^n J_n(κ)·φ_{R+i·n·k_L/B}\n";
    std::cout << "      也就是每个原 ECG 派生出 (2n_max+1) 个新 R 的 ECG.\n";
    std::cout << "      apply_analytic_kick 没真的加新基函数, 而是把这个\n";
    std::cout << "      \"应在新 R 上的态\" 投影到原 R 的基底里, 只改 u.\n";
    std::cout << "      → fidelity < 1 时, 说明原基底没覆盖这些移位 ECG,\n";
    std::cout << "        投影丢了信息. (策略 HTML 里 Strategy A vs C 的核心)\n";
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


#pragma region 任务 D: Route 1 target-fitting TDVP + finite-diff 验证
// ============================================================================
// 【任务 D】Route 1 目标拟合 TDVP（单次 kick 用虚拟 H_eff target-fit）
//
// 数学：|ψ_tar⟩ = e^{-iκcos(2k_L x)}|ψ_0⟩ 用 Jacobi-Anger 截断到 |n|≤n_max
//       H_eff = -|ψ_tar⟩⟨ψ_tar|/N_tar
//       g_α = -(A/N_tar)·b_α, A = ⟨ψ_tar|ψ⟩, b_α = ⟨∂_α ψ|ψ_tar⟩
//       TDVP: C·ż = -g
//
// 两个子命令：
//   --task4           跑完整 target-fitting 循环，对比 task3 的 apply_analytic_kick
//   --verify-g-alpha  用中心差分对拍 assemble_grad_target_fit 的 g_α 公式
// ============================================================================

// Compute E_unnorm(z) = <psi(z)|H_eff|psi(z)> = -|A(z)|^2 / N_tar (unnormalized
// by <psi|psi>). Analytic g_alpha matches ∂E_unnorm/∂z_alpha^* (P⊥ absorbs the S
// normalization in the full energy -|A|^2/(N_tar·S)).
static Cd compute_E_unnorm(const std::vector<BasisParams>& trial,
                           const TargetBasis& target,
                           Cd N_tar,
                           const PermutationSet& perms) {
    Cd A = cross_overlap(target, trial, perms);
    return -A * std::conj(A) / N_tar;
}

// Perturb basis along a single AlphaIndex by complex amount `delta` using
// update_basis_function. Returns a perturbed copy.
static std::vector<BasisParams> perturb_basis(const std::vector<BasisParams>& basis,
                                              const AlphaIndex& alpha,
                                              Cd delta) {
    std::vector<BasisParams> perturbed = basis;
    std::vector<AlphaIndex> single{alpha};
    VectorXcd dz(1);
    dz(0) = delta;
    update_basis_function(perturbed, dz, 1.0, single, /*updata_constant=*/0);
    return perturbed;
}

// Central-difference Wirtinger derivative: ∂E/∂z_α^* ≈ 0.5·(FD_x + i·FD_y)
// where FD_x = (E(+δ)-E(-δ))/(2δ) with real δ perturbation of u_α (or A/B/R)
// and FD_y = (E(+iδ)-E(-iδ))/(2δ) with imaginary perturbation.
static Cd finite_diff_g_alpha(const std::vector<BasisParams>& basis,
                              const AlphaIndex& alpha,
                              const TargetBasis& target,
                              Cd N_tar,
                              const PermutationSet& perms,
                              double delta) {
    Cd E_px = compute_E_unnorm(perturb_basis(basis, alpha, Cd(+delta, 0)),
                               target, N_tar, perms);
    Cd E_mx = compute_E_unnorm(perturb_basis(basis, alpha, Cd(-delta, 0)),
                               target, N_tar, perms);
    Cd FD_x = (E_px - E_mx) / (2.0 * delta);

    Cd E_py = compute_E_unnorm(perturb_basis(basis, alpha, Cd(0, +delta)),
                               target, N_tar, perms);
    Cd E_my = compute_E_unnorm(perturb_basis(basis, alpha, Cd(0, -delta)),
                               target, N_tar, perms);
    Cd FD_y = (E_py - E_my) / (2.0 * delta);

    return 0.5 * (FD_x + Cd(0, 1) * FD_y);
}

// Prepare a K=1 ground-state ECG basis (mirrors task3 Step 1, trimmed to essentials).
// For K=1 harmonic oscillator, SVM finds the exact ground state (single Gaussian),
// so we skip TDVP iff the initial energy is already at 0.5 (within 1e-8). Otherwise
// run a short TDVP pass (max 50 steps) as a safety net.
static std::vector<BasisParams> prepare_ground_state_K1(int K) {
    const int N = 1;
    auto terms = HamiltonianTerms::kinetic_harmonic();
    PermutationSet perms = PermutationSet::generate(N);

    SvmResult svm = svm_build_basis(N, K, 2000, terms, 42, 0.0);
    std::vector<BasisParams> basis = svm.basis;
    set_u_from_eigenvector(basis, svm.H, svm.S);

    Cd E_init = compute_total_energy(basis, terms);
    if (std::abs(E_init.real() - 0.5 * N) < 1e-8) {
        return basis;   // exactly at ground state; no need for TDVP
    }

    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K; i++) alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) alpha_z_list.push_back({2, i, j, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) alpha_z_list.push_back({3, i, j, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
            for (int k = j; k < N; k++) alpha_z_list.push_back({4, i, j, k});

    SolverConfig cfg;
    cfg.lambda_C = 1e-8;
    cfg.rcond = 1e-4;
    cfg.resolve_u = true;
    cfg.dtao_grow = 1.5;
    cfg.dtao_max = 10.0;
    cfg.energy_tol = 1e-12;
    cfg.optimize_A = true;
    cfg.optimize_B = true;
    cfg.optimize_R = true;

    evolution(alpha_z_list, basis, 1e-3, /*max_steps=*/50, 1e-12, terms, cfg, &perms);
    auto [H_gs, S_gs] = build_HS(basis, perms, terms);
    set_u_from_eigenvector(basis, H_gs, S_gs);
    return basis;
}

static void verify_g_alpha(int K, int n_max, double delta) {
    std::cout << "\n=====================================================\n";
    std::cout << " 验证 g_α：有限差分对拍 (K=" << K << ", n_max=" << n_max
              << ", δ=" << delta << ")\n";
    std::cout << "=====================================================\n";

    const int N = 1;
    PermutationSet perms = PermutationSet::generate(N);

    // Step 1: K=1 pre-kick 基态（rebalance 使 B>0 便于 kick 编码）
    std::cout << "\n--- Step 1: 准备 K=" << K << " 基态 ---\n";
    std::vector<BasisParams> basis_raw = prepare_ground_state_K1(K);
    std::vector<BasisParams> basis = rebalance_AB_for_kick(basis_raw, /*b_min=*/0.5);
    std::cout << "  rebalanced basis (A+B 不变, B>=0.5 使 kick 编码良定):\n";
    print_basis("rebalanced basis", basis);

    // Step 2: 故意扰动使 z ≠ 最优（否则 g_α 全是 0，看不出 bug）
    std::cout << "\n--- Step 2: 扰动 basis 使其离最优点 ---\n";
    for (int i = 0; i < K; i++) {
        basis[i].A(0, 0) += Cd(0.05 * (i + 1), 0.02);
        basis[i].B(0, 0) += Cd(0.03, -0.01);
        basis[i].R(0)    += Cd(0.1, -0.05);
        basis[i].u       += Cd(0.01, 0.02);
    }
    print_basis("扰动后 trial basis", basis);

    // Step 3: target（同样从 rebalanced 基底构造）
    std::cout << "\n--- Step 3: 构造目标态 ---\n";
    std::vector<BasisParams> pre_kick_raw = prepare_ground_state_K1(K);
    std::vector<BasisParams> pre_kick_basis = rebalance_AB_for_kick(pre_kick_raw, 0.5);
    TargetBasis target = build_kicked_target(pre_kick_basis, kappa, k_L, n_max);
    std::cout << "  |target| = " << target.size()
              << " (K=" << K << " × (2·" << n_max << "+1))\n";
    Cd N_tar = target_norm(target, perms);
    std::cout << "  N_tar = " << N_tar << " (虚部应 ≈ 0)\n";

    // Step 4: alpha_z_list (全部参数)
    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K; i++) alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K; i++) alpha_z_list.push_back({2, i, 0, 0});
    for (int i = 0; i < K; i++) alpha_z_list.push_back({3, i, 0, 0});
    for (int i = 0; i < K; i++) alpha_z_list.push_back({4, i, 0, 0});

    // Step 5: analytic g
    std::cout << "\n--- Step 5: 解析 g_α ---\n";
    Cd A_val, N_chk;
    VectorXcd g_analytic = assemble_grad_target_fit(alpha_z_list, basis, target,
                                                    perms, A_val, N_chk);
    std::cout << "  A(z) = " << A_val << "\n";
    std::cout << "  |A|^2 / N_tar = " << std::norm(A_val) / N_tar.real() << "\n";

    // Step 6: finite diff 对每个 alpha
    std::cout << "\n--- Step 6: 有限差分 g_α 对拍 ---\n";
    std::cout << std::setw(4) << "idx"
              << std::setw(4) << "a1"
              << std::setw(4) << "a2"
              << std::setw(4) << "a3"
              << std::setw(4) << "a4"
              << std::setw(36) << "g_analytic"
              << std::setw(36) << "g_finite_diff"
              << std::setw(14) << "|rel_err|\n";
    std::cout << std::string(102, '-') << "\n";

    double max_rel_err = 0.0;
    int n_pass = 0, n_fail = 0;
    const double tol = 1e-5;

    for (int a = 0; a < static_cast<int>(alpha_z_list.size()); a++) {
        Cd g_num = finite_diff_g_alpha(basis, alpha_z_list[a], target,
                                       N_tar, perms, delta);
        Cd g_ana = g_analytic(a);
        double denom = std::max(std::abs(g_ana), 1e-10);
        double rel_err = std::abs(g_num - g_ana) / denom;
        if (rel_err > max_rel_err) max_rel_err = rel_err;
        const char* mark = (rel_err < tol) ? "PASS" : "FAIL";
        if (rel_err < tol) n_pass++; else n_fail++;

        std::cout << std::setw(4) << a
                  << std::setw(4) << alpha_z_list[a].a1
                  << std::setw(4) << alpha_z_list[a].a2
                  << std::setw(4) << alpha_z_list[a].a3
                  << std::setw(4) << alpha_z_list[a].a4
                  << "  " << std::setw(34) << fmt_cd(g_ana, 6)
                  << "  " << std::setw(34) << fmt_cd(g_num, 6)
                  << "  " << std::scientific << std::setprecision(3)
                  << std::setw(11) << rel_err << " " << mark
                  << std::defaultfloat << "\n";
    }

    std::cout << "\n--- 结果 ---\n";
    std::cout << "  PASS: " << n_pass << " / " << alpha_z_list.size() << "\n";
    std::cout << "  FAIL: " << n_fail << " / " << alpha_z_list.size() << "\n";
    std::cout << "  max |rel_err| = " << std::scientific << max_rel_err
              << std::defaultfloat << "   (tol=" << tol << ")\n";
    std::cout << "  "
              << (n_fail == 0 ? "✓ g_α 公式验证通过" : "✗ g_α 公式有偏差！")
              << "\n";
}

static void task4_target_fitting_tdvp(int K_pre, int n_max, int n_mom,
                                       TargetFitConfig::UMode u_mode) {
    std::cout << "\n=====================================================\n";
    std::cout << " 任务 D: Route 1 target-fitting TDVP (K_pre=" << K_pre
              << ", n_max=" << n_max
              << ", n_mom=" << n_mom
              << ", u_mode="
              << (u_mode == TargetFitConfig::UPDATE_LINEAR ? "LINEAR" : "RESOLVE")
              << ")\n";
    std::cout << "=====================================================\n";

    const int N = 1;
    PermutationSet perms = PermutationSet::generate(N);

    // Step 1: K=1 pre-kick 基态（作为 trial 初值 + target 构造源）
    std::cout << "\n--- Step 1: 准备 pre-kick 基态 ---\n";
    std::vector<BasisParams> pre_kick_raw = prepare_ground_state_K1(K_pre);
    std::vector<BasisParams> pre_kick_basis = rebalance_AB_for_kick(pre_kick_raw, 0.5);
    print_basis("pre-kick basis (rebalanced, B≥0.5)", pre_kick_basis);

    // Step 2: 踢之前观测量（确认基态正确）
    std::cout << "\n--- Step 2: 踢之前观测量 ---\n";
    StateObs before = compute_observables(pre_kick_basis);
    print_observables("踢之前", before);

    // Step 3: 目标态
    std::cout << "\n--- Step 3: 构造目标态 |ψ_tar⟩ ---\n";
    TargetBasis target = build_kicked_target(pre_kick_basis, kappa, k_L, n_max);
    Cd N_tar = target_norm(target, perms);
    std::cout << "  |target| = " << target.size() << "\n";
    std::cout << "  N_tar = " << N_tar << "\n";

    // Step 3.5: 动量扩基 —— K_pre=1 的单高斯装不下多动量扇区，
    // 所以在 target-fitting 之前先用 augment_basis_with_momentum 把 n=±1,...,±n_mom
    // 动量副本加进 trial basis。TDVP 再微调每个分量的 (u, A, B, R)。
    std::cout << "\n--- Step 3.5: 动量扩基 (augment_basis_with_momentum, n_mom=" << n_mom
              << ") ---\n";
    std::vector<BasisParams> trial_init = augment_basis_with_momentum(
        pre_kick_basis, k_L, /*n_mom=*/n_mom, /*b_val=*/0.5, /*max_cond=*/1e6);
    const int K_trial = static_cast<int>(trial_init.size());
    std::cout << "  K_pre=" << pre_kick_basis.size()
              << " → K_trial=" << K_trial << "\n";
    print_basis("扩基后 trial basis", trial_init);

    // 初始 fidelity（trial = 扩基后, target = kicked state）
    // 新加入的动量副本 u=0, 先做一次 optimal-u 再看 F，否则数字会非常难看
    {
        VectorXcd u_init = solve_optimal_u(trial_init, target, perms);
        for (int i = 0; i < K_trial; i++) trial_init[i].u = u_init(i);
    }
    Cd A0 = cross_overlap(target, trial_init, perms);
    Cd S0 = overlap(trial_init);
    double F0 = std::norm(A0) / (N_tar.real() * S0.real());
    std::cout << "  初始 F(augmented trial | target) = " << std::setprecision(10) << F0
              << ", 1-F = " << std::scientific << (1.0 - F0)
              << std::defaultfloat << "\n";

    // Step 4: alpha_z_list（对 K_trial 个基函数的 u, B, R, A 全放开）
    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({2, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({3, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({4, i, 0, 0});

    // Step 5: target-fitting evolution
    std::cout << "\n--- Step 5: Target-fitting TDVP 主循环 ---\n";
    TargetFitConfig fit_cfg;
    fit_cfg.u_mode = u_mode;
    fit_cfg.max_steps = 2000;
    fit_cfg.dtao_init = 1e-3;
    fit_cfg.fidelity_tol = 1e-10;
    fit_cfg.verbose = true;
    fit_cfg.print_every = 50;

    TargetFitResult result = target_fitting_evolution(alpha_z_list,
                                                      trial_init,
                                                      target, fit_cfg);

    // Step 6: 结果
    std::cout << "\n--- Step 6: 收敛后基底 & 观测量 ---\n";
    print_basis("收敛后 basis", result.basis);
    StateObs after = compute_observables(result.basis);
    print_observables("Route 1 后", after);

    // 理论值（与 task3 Step 4 相同）
    double x2_expected = 0.5;
    double sin2_avg    = 0.5 * (1.0 - std::exp(-4.0 * k_L * k_L));
    double p2_expected = 0.5 + 4.0 * kappa * kappa * k_L * k_L * sin2_avg;
    double E_expected  = 0.5 * (p2_expected + x2_expected);

    std::cout << "\n--- 与理论对比 ---\n";
    std::cout << std::setw(12) << "量" << std::setw(18) << "Route 1"
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
              << std::setw(18) << x2_expected << "\n";
    std::cout << std::setw(12) << "⟨p²⟩"
              << std::setw(18) << after.p2
              << std::setw(18) << p2_expected
              << std::setw(14) << std::scientific << (after.p2 - p2_expected)
              << std::defaultfloat << "\n";

    std::cout << "\n--- 摘要 ---\n";
    std::cout << "  F_final        = " << std::setprecision(12) << result.F_final << "\n";
    std::cout << "  1 - F_final    = " << std::scientific << (1.0 - result.F_final) << std::defaultfloat << "\n";
    std::cout << "  Steps          = " << result.n_steps << "\n";
    std::cout << "  Line-search fails = " << result.line_search_fails << "\n";
}
#pragma endregion


#pragma region 任务 E: Route 1 踢后态的实时 TDVP 演化（T+V_ho）
// ============================================================================
// 【任务 E】对任务 D 得到的 kick 后态跑实时 TDVP，哈密顿量 H = T + V_ho
//
// 物理：
//   单粒子谐振子，经典轨迹周期 T_ω = 2π/ω = 2π（ω=1）
//   ⟨x²⟩(t), ⟨p²⟩(t) 应该做周期 π 的振荡（因为 x² ~ cos²(ωt+φ) 周期 π）
//   回归：|⟨ψ(0)|ψ(T_ω)⟩|² 应该 = 1（模全局相位）
//   能量 E(t) 严格守恒（H 不含 t，TDVP 保能量）
//
// 实现：
//   - 共用任务 D 的踢后态准备流程（rebalance + augment + Route 1 target-fit）
//   - 主循环 realtime_tdvp_evolution 默认用 RK4 + dt=1e-3
//   - 采样 E/norm/⟨x²⟩/⟨p²⟩/⟨ψ(0)|ψ(t)⟩ 轨迹
// ============================================================================

// 复用任务 D 的踢后态准备流程，返回 Route 1 收敛后的 basis（K_trial 维）
static std::vector<BasisParams> prepare_kicked_state_route1(
    int K_pre, int n_max, int n_mom,
    TargetFitConfig::UMode u_mode = TargetFitConfig::UPDATE_RESOLVE,
    bool verbose = false) {

    const int N = 1;
    PermutationSet perms = PermutationSet::generate(N);

    std::vector<BasisParams> pre_kick_raw = prepare_ground_state_K1(K_pre);
    std::vector<BasisParams> pre_kick_basis = rebalance_AB_for_kick(pre_kick_raw, 0.5);

    TargetBasis target = build_kicked_target(pre_kick_basis, kappa, k_L, n_max);
    std::vector<BasisParams> trial_init = augment_basis_with_momentum(
        pre_kick_basis, k_L, n_mom, 0.5, 1e6);
    const int K_trial = static_cast<int>(trial_init.size());

    // 初始 u 解析最优化
    {
        VectorXcd u_init = solve_optimal_u(trial_init, target, perms);
        for (int i = 0; i < K_trial; i++) trial_init[i].u = u_init(i);
    }

    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({2, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({3, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({4, i, 0, 0});

    TargetFitConfig fit_cfg;
    fit_cfg.u_mode = u_mode;
    fit_cfg.max_steps = 2000;
    fit_cfg.dtao_init = 1e-3;
    fit_cfg.fidelity_tol = 1e-12;
    fit_cfg.verbose = verbose;

    TargetFitResult result = target_fitting_evolution(alpha_z_list, trial_init,
                                                      target, fit_cfg);
    if (verbose) {
        std::cout << "  [prepare_kicked_state_route1] F_final=" << std::setprecision(10)
                  << result.F_final << "\n";
    }
    return result.basis;
}

static void task5_realtime_evolve_kicked(int n_max, int n_mom,
                                          double T_total, double dt,
                                          RtIntegrator integrator) {
    std::cout << "\n=====================================================\n";
    std::cout << " 任务 E: Route 1 踢后态的实时 TDVP 演化 (H = T + V_ho)\n";
    std::cout << "   n_max=" << n_max << ", n_mom=" << n_mom
              << ", T_total=" << T_total << ", dt=" << dt
              << ", integrator=" << (integrator == RtIntegrator::RK4 ? "RK4" : "Euler")
              << "\n";
    std::cout << "=====================================================\n";

    const int N = 1;

    // Step 1: 踢后态
    std::cout << "\n--- Step 1: 准备 kick 后的 ECG 态 (Route 1) ---\n";
    std::vector<BasisParams> kicked_basis = prepare_kicked_state_route1(
        /*K_pre=*/1, n_max, n_mom, TargetFitConfig::UPDATE_RESOLVE, /*verbose=*/true);
    const int K_trial = static_cast<int>(kicked_basis.size());
    std::cout << "  K_trial = " << K_trial << "\n";

    // 初始观测
    StateObs init = compute_observables(kicked_basis);
    print_observables("t=0 (踢后)", init);

    // Step 2: alpha_z_list 全开（u, B, R, A 都演化）
    std::vector<AlphaIndex> alpha_z_list;
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({1, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({2, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({3, i, 0, 0});
    for (int i = 0; i < K_trial; i++) alpha_z_list.push_back({4, i, 0, 0});

    // Step 3: 实时 TDVP 演化
    std::cout << "\n--- Step 2: 实时 TDVP 演化 t ∈ [0, " << T_total << "] ---\n";
    HamiltonianTerms terms = HamiltonianTerms::kinetic_harmonic();  // 只含 T + V_ho
    SolverConfig solver_cfg;
    solver_cfg.lambda_C = 0.0;   // 实时演化，去掉 Tikhonov 以保幺正性
    solver_cfg.rcond = 1e-4;

    RealtimeEvolutionConfig rt_cfg;
    rt_cfg.dt = dt;
    rt_cfg.integrator = integrator;
    rt_cfg.sample_every = std::max(1, static_cast<int>(std::round(T_total / dt / 40)));
    rt_cfg.verbose = true;
    rt_cfg.print_every = std::max(50, rt_cfg.sample_every * 5);

    RealtimeEvolutionResult rt = realtime_tdvp_evolution(
        alpha_z_list, kicked_basis, T_total, terms, solver_cfg, rt_cfg);

    // Step 4: 轨迹打印
    std::cout << "\n--- Step 3: 轨迹 E(t), ⟨x²⟩(t), ⟨p²⟩(t), norm(t), |F(t)|² ---\n";
    std::cout << std::setw(10) << "t"
              << std::setw(14) << "E"
              << std::setw(14) << "norm"
              << std::setw(14) << "<x²>"
              << std::setw(14) << "<p²>"
              << std::setw(16) << "|<psi0|psi(t)>|²\n";
    std::cout << std::string(82, '-') << "\n";
    double norm0 = rt.trace.norm[0];
    for (size_t k = 0; k < rt.trace.t.size(); k++) {
        double fid = std::norm(rt.trace.overlap0[k]) / (norm0 * rt.trace.norm[k]);
        std::cout << std::setw(10) << std::fixed << std::setprecision(4) << rt.trace.t[k]
                  << std::setw(14) << std::setprecision(8) << rt.trace.E[k]
                  << std::setw(14) << rt.trace.norm[k]
                  << std::setw(14) << rt.trace.x2[k]
                  << std::setw(14) << rt.trace.p2[k]
                  << std::setw(16) << fid << "\n";
    }

    // 守恒量偏差
    double E_init = rt.trace.E.front();
    double E_final = rt.trace.E.back();
    double norm_final = rt.trace.norm.back();
    double dE_max = 0.0, dnorm_max = 0.0;
    for (size_t k = 0; k < rt.trace.t.size(); k++) {
        dE_max = std::max(dE_max, std::abs(rt.trace.E[k] - E_init));
        dnorm_max = std::max(dnorm_max, std::abs(rt.trace.norm[k] - norm0) / norm0);
    }
    double fid_T = std::norm(rt.trace.overlap0.back()) / (norm0 * norm_final);

    std::cout << "\n--- 守恒量 & 周期性 ---\n";
    std::cout << "  E(0)           = " << std::setprecision(10) << E_init << "\n";
    std::cout << "  E(T_total)     = " << E_final << "\n";
    std::cout << "  max |ΔE|       = " << std::scientific << dE_max
              << std::defaultfloat << "\n";
    std::cout << "  max |Δnorm|/n0 = " << std::scientific << dnorm_max
              << std::defaultfloat << "\n";
    std::cout << "  |⟨ψ₀|ψ(T)⟩|²  = " << std::setprecision(10) << fid_T
              << "  (T=2π 理论 → 1.0)\n";
    std::cout << "  收敛步数       = " << rt.n_steps << "\n";
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
    bool run_task4 = false;
    bool run_task5 = false;
    bool run_verify = false;
    int  task1_N   = 1;   // 任务 A: 粒子数
    int  task1_K   = 1;   // 任务 A: 基函数数量
    int  n_kicks   = 40;  // 任务 B: kick 次数
    int  task3_K   = 5;   // 任务 C: 基函数数量（N=1）
    int  task4_K   = 1;   // 任务 D: K_pre (pre-kick 基底数, 默认 1 够 harmonic ground state)
    int  task4_nmax = 4;  // 任务 D: Jacobi-Anger 截断
    int  task4_nmom = 2;  // 任务 D: 动量扩基每侧个数 (trial K' = K_pre + 2·n_mom 个 momentum 副本)
    std::string task4_umode = "resolve"; // "resolve" 或 "linear"
    double verify_delta = 1e-5;
    double task5_T      = 2.0 * M_PI;     // 任务 E: 演化时长 (一个谐振周期)
    double task5_dt     = 1e-3;           // 任务 E: 积分步长
    std::string task5_integrator = "rk4"; // "rk4" 或 "euler"

    // 解析命令行
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if      (arg == "--task1")   run_task1 = true;
        else if (arg == "--task2")   run_task2 = true;
        else if (arg == "--task3")   run_task3 = true;
        else if (arg == "--task4")   run_task4 = true;
        else if (arg == "--task5")   run_task5 = true;
        else if (arg == "--verify-g-alpha") run_verify = true;
        else if (arg == "--N"        && i + 1 < argc) task1_N  = std::stoi(argv[++i]);
        else if (arg == "--K"        && i + 1 < argc) task1_K  = std::stoi(argv[++i]);
        else if (arg == "--K3"       && i + 1 < argc) task3_K  = std::stoi(argv[++i]);
        else if (arg == "--K4"       && i + 1 < argc) task4_K  = std::stoi(argv[++i]);
        else if (arg == "--n-max"    && i + 1 < argc) task4_nmax = std::stoi(argv[++i]);
        else if (arg == "--n-mom"    && i + 1 < argc) task4_nmom = std::stoi(argv[++i]);
        else if (arg == "--u-mode"   && i + 1 < argc) task4_umode = argv[++i];
        else if (arg == "--delta"    && i + 1 < argc) verify_delta = std::stod(argv[++i]);
        else if (arg == "--T-free"   && i + 1 < argc) task5_T = std::stod(argv[++i]);
        else if (arg == "--dt"       && i + 1 < argc) task5_dt = std::stod(argv[++i]);
        else if (arg == "--integrator" && i + 1 < argc) task5_integrator = argv[++i];
        else if (arg == "--n-kicks"  && i + 1 < argc) n_kicks  = std::stoi(argv[++i]);
        else {
            std::cerr << "未知参数: " << arg << "\n";
            std::cerr << "用法: " << argv[0]
                      << " [--task1] [--task2] [--task3] [--task4] [--verify-g-alpha]"
                      << " [--N 粒子数] [--K 基函数数(任务A)]"
                      << " [--K3 基函数数(C)] [--K4 K_pre(D,默认1)]"
                      << " [--n-max 4] [--n-mom 2] [--u-mode resolve|linear]"
                      << " [--delta 1e-5] [--n-kicks N]\n";
            return 1;
        }
    }

    // 不指定则全部都跑（只含 A/B/C；task4/5/verify 要显式触发）
    if (!run_task1 && !run_task2 && !run_task3 && !run_task4 && !run_task5 && !run_verify) {
        run_task1 = true;
        run_task2 = true;
        run_task3 = true;
    }

    if (run_task1) task1_harmonic_ground_state(task1_N, task1_K);
    if (run_task2) task2_grid_kicked_rotor(n_kicks);
    if (run_task3) task3_single_kick(task3_K);

    if (run_verify) verify_g_alpha(task4_K, task4_nmax, verify_delta);
    if (run_task4) {
        TargetFitConfig::UMode mode = TargetFitConfig::UPDATE_RESOLVE;
        if (task4_umode == "linear") mode = TargetFitConfig::UPDATE_LINEAR;
        task4_target_fitting_tdvp(task4_K, task4_nmax, task4_nmom, mode);
    }
    if (run_task5) {
        RtIntegrator itg = (task5_integrator == "euler")
                           ? RtIntegrator::Euler
                           : RtIntegrator::RK4;
        task5_realtime_evolve_kicked(task4_nmax, task4_nmom, task5_T, task5_dt, itg);
    }

    return 0;
}
#pragma endregion
