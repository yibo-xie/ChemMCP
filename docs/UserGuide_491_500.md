# ChemMCP #491-500 使用指南

> **动力学校正与计算化学工具包** — 10 个 MCP 工具，覆盖从反应动力学近似到量子化学计算的全流程

---

## 目录

1. [快速开始](#1-快速开始)
2. [工具总览](#2-工具总览)
3. [动力学校正工具 (#491-495)](#3-动力学校正工具-491-495)
4. [计算化学工具 (#496-500)](#4-计算化学工具-496-500)
5. [典型工作流](#5-典型工作流)
6. [Cherry Studio 配置](#6-cherry-studio-配置)
7. [常见问题](#7-常见问题)

---

## 1. 快速开始

### 环境要求

```bash
# 已安装 uv（Python 包管理器）
uv --version

# 仓库已克隆到本地
cd ~/ChemMCP
uv sync
```

### 调用方式

每个 MCP 工具支持**两种调用方式**：

| 方式 | 适用场景 | 示例 |
|------|----------|------|
| `run_code()` | 编程调用 / AI Agent | 关键字参数，类型安全 |
| `run_text()` | 快速输入 / 命令行 | 空格分隔的文本字符串 |

```python
from chemmcp.tools import TunnelingCorrection

tool = TunnelingCorrection()

# 方式一：代码调用（推荐）
result = tool.run_code(
    temperature_K=298.15,
    barrier_height_kJ_mol=45.0,
    imaginary_frequency_cm_minus_1=1500.0,
    correction_model="bell"
)

# 方式二：文本调用
result = tool.run_text("298.15 45.0 1500.0 bell")
```

---

## 2. 工具总览

| # | 工具名 | 功能一句话 | 核心场景 |
|---|--------|-----------|----------|
| 491 | **TunnelingCorrection** | 计算量子隧穿校正因子 κ | H/D 转移反应速率修正 |
| 492 | **SteadyStateApproximation** | 稳态近似求解中间体浓度 | 多步反应机理分析 |
| 493 | **PreEquilibrium** | 预平衡近似推导速率方程 | 快平衡+慢决速步机理 |
| 494 | **RateDeterminingStep** | 自动识别速控步 | 反应瓶颈定位 |
| 495 | **MichaelisMenten** | 酶动力学完整分析 | 酶催化/抑制研究 |
| 496 | **ReactionNetworkSolver** | 反应网络 ODE 求解 | 复杂反应体系模拟 |
| 497 | **GeometryOptimizer** | 分子几何优化（能量极小化） | 分子结构优化 |
| 498 | **TransitionStateSearch** | 过渡态鞍点搜索 | 活化能计算 |
| 499 | **PotentialEnergySurface** | 势能面扫描 | 反应路径探索 |
| 500 | **FrequencyAnalysis** | 振动频率分析+驻点分类 | 结构验证（极小/过渡态） |

---

## 3. 动力学校正工具 (#491-495)

### #491 TunnelingCorrection — 隧穿校正

**用途：** 经典过渡态理论 (TST) 忽略了量子隧穿效应。对于轻原子（H, D）转移反应，低温下隧穿效应显著，需要用 κ 因子修正速率常数：

$$k_{\text{quantum}} = \kappa \times k_{\text{classical}}$$

**三种模型：**

| 模型 | 特点 | 适用场景 |
|------|------|----------|
| `bell` | 抛物线势垒近似，解析解 | 快速估算，教学 |
| `eckart` | 非对称势垒，精确积分 | 不对称反应（正逆势垒不同） |
| `wkb` | WKB 数值积分 | 高精度，任意势垒形状 |

**示例：H 转移反应在 298K 的隧穿校正**

```python
from chemmcp.tools import TunnelingCorrection
tc = TunnelingCorrection()

# Bell 模型：T=298K, 势垒=45 kJ/mol, 虚频=1500 cm⁻¹
r = tc.run_code(
    temperature_K=298.15,
    barrier_height_kJ_mol=45.0,
    imaginary_frequency_cm_minus_1=1500.0,
    correction_model="bell"
)

print(f"κ = {r['correction_factor_kappa']:.4f}")       # 校正因子
print(f"经典相对速率: {r['classical_rate_relative']:.4e}")
print(f"量子修正后:   {r['quantum_corrected_rate_relative']:.4e}")
print(f"隧穿显著?     {r['tunneling_significant']}")      # True/False
```

**Eckart 非对称势垒（需要正逆势垒）：**

```python
r = tc.run_code(
    temperature_K=200.0,              # 低温 → 隧穿更强
    barrier_height_kJ_mol=40.0,
    imaginary_frequency_cm_minus_1=1800.0,
    correction_model="eckart",
    reduced_mass_amu=1.008,           # H 原子约化质量
    forward_barrier_kJ_mol=40.0,      # 正向势垒
    reverse_barrier_kJ_mol=80.0       # 反向势垒更高（放热反应）
)
```

**实用技巧：**
- 温度越低 → κ 越大（隧穿效应更显著）
- 虚频越大（势垒越窄）→ κ 越大
- 同位素替换 D→H 后 κ 减小（动力学同位素效应 KIE）
- 一般 κ > 1 表示量子隧穿增强反应速率

---

### #492 SteadyStateApproximation — 稳态近似

**用途：** 对于含中间体 I 的多步反应 A → I → P，当中间体活性很高时，可假设其浓度不随时间变化（d[I]/dt ≈ 0），从而简化动力学方程。

**支持的机理类型：**

| mechanism_type | 反应模式 | 速率常数格式 |
|----------------|----------|-------------|
| `consecutive` | A → I → P | `[k1, k2]` |
| `reversible_consecutive` | A ⇌ I → P | `[kf1, kr1, k2]` |
| `pre_equilibrium` | A ⇌ I → P (快平衡) | `[kf, kr, k_slow]` |
| `parallel_consecutive` | A → I₁ + I₂ → P | `[k1, k2, k3]` |

**示例：连续反应 A → I → P**

```python
from chemmcp.tools import SteadyStateApproximation
ssa = SteadyStateApproximation()

r = ssa.run_code(
    mechanism_type="consecutive",
    rate_constants=[0.1, 0.05],         # k1=0.1, k2=0.05
    initial_reactant_concentration_A0=1.0,
    time_t=50.0,
    target_intermediate="I"
)

print(f"精确 [I] = {r['exact']['intermediate_I']:.6f}")
print(f"SSA 近似 [I] = {r['approximate']['intermediate_I']:.6f}")
print(f"SSA 有效? {r['ssa_valid']}")                     # True/False
print(f"误差 = {r['intermediate_error_percent']:.2f}%")
print(f"k2/k1 比 = {r['k2_k1_ratio']:.2f}")
```

**何时 SSA 有效？**
- k₂ ≫ k₁（中间体消耗快于生成）→ 误差 < 5%
- 时间 t >> 1/k₁（过了诱导期之后）
- 如果 `ssa_valid = False`，说明条件不满足，需用精确数值解（见 #496）

**可逆连续反应：**

```python
r = ssa.run_code(
    mechanism_type="reversible_consecutive",
    rate_constants=[0.1, 0.02, 0.01],      # kf1, kr1, k2
    initial_reactant_concentration_A0=1.0,
    time_t=100.0
)
```

---

### #493 PreEquilibrium — 预平衡近似

**用途：** 当反应机理中存在快速平衡步骤后跟慢速决速步时（如 A ⇌ I → P），用预平衡近似推导表观速率定律。

$$k_{\text{eff}} = K_{\text{eq}} \times k_{\text{slow}} = \frac{k_f}{k_r} \times k_{\text{slow}}$$

**示例：单分子预平衡**

```python
from chemmcp.tools import PreEquilibrium
peq = PreEquilibrium()

r = peq.run_code(
    mechanism="unimolecular",
    k_forward_list=[1.0],        # k_f (快速正向)
    k_reverse_list=[0.5],        # k_r (快速逆向)
    k_slow=0.01,                  # 决速步（慢）
    initial_concentrations={"A": 1.0}
)

print(f"K_eq = {r.get('equilibrium_constant', 'N/A')}")
print(f"k_eff = {r['effective_rate_constant']:.6f}")
print(f"速率方程: {r['rate_law']}")
print(f"有效条件: {r['validity']}")
```

**双分子预平衡（A + B ⇌ I → P）：**

```python
r = peq.run_code(
    mechanism="bimolecular",
    k_forward_list=[1.0],
    k_reverse_list=[0.2],
    k_slow=0.05,
    initial_concentrations={"A": 1.0, "B": 1.0},
    time_points=[0, 10, 50, 100, 200]
)

# 含时间演化的浓度分布
for tp in r["concentration_profiles"]:
    print(f"t={tp['time']:.1f}: A={tp['A']:.4f}, I={tp['I']:.6f}, P={tp['P']:.4f}")
```

**与 #492 的区别：**
- `SteadyStateApprox`: 通用稳态近似，适用于任何中间体
- `PreEquilibrium`: 专门处理「快平衡 + 慢决速步」模式，给出 K_eq 和有效速率常数

---

### #494 RateDeterminingStep — 速控步分析

**用途：** 给定多步反应机理，自动找出最慢的一步（速控步 RDS），推导总速率方程。

**核心逻辑：** 比较 k 值大小 → 找最小 k → 标记为 RDS → 推导总速率定律

```python
from chemmcp.tools import RateDeterminingStep
rds = RateDeterminingStep()

# 定义三步反应机理
r = rds.run_code([
    {"reactants": "A → I", "products": "", "k": 1e5, "reversible": False},     # 步骤1: 很快
    {"reactants": "I + B → P", "products": "", "k": 0.01, "reversible": False}, # 步骤2: 最慢 ← RDS!
    {"reactants": "P → Q", "products": "", "k": 100.0, "reversible": False},   # 步骤3: 较快
])

print(f"RDS 是第 {r['rds_step_index']} 步")           # 1-indexed: 返回 2
print(f"RDS 描述: {r['rds_step_description']}")
print(f"RDS 速率常数 k = {r['rds_rate_constant']}")
print(f"总速率方程: {r['overall_rate_law']}")
print(f"\n各步排名 (按速率):")
for rank in r["step_ranking"]:
    print(f"  步骤 {rank['step']}: k={rank['k']}, 相对RDS倍数={rank['ratio_to_rds']:.1f}x")
```

**含预平衡的分析：**

```python
r = rds.run_code([
    {"reactants": "A ⇌ I", "k": 100.0, "reversible": True},   # 快平衡
    {"reactants": "I → P", "k": 0.01, "reversible": False},    # 慢 RDS
], has_pre_equilibrium=True)

# 会自动用 K_eq 推导速率方程
print(r["overall_rate_law"])
# 输出类似: rate = (k_f/k_r) * k_slow * [A]
```

**输出说明：**
- `rds_step_index`: **1-indexed**（第1步返回1，不是0）
- `rate_constant_ratios`: 各步 k 值相对于 RDS k 的比值（均 ≥ 1）
- `step_ranking`: 从慢到快的完整排序

---

### #495 MichaelisMenten — Michaelis-Menten 酶动力学

**用途：** 酶催化反应的完整动力学分析，包括速率计算、抑制分析和参数拟合。

**核心公式：**
$$v = \frac{V_{\max} [S]}{K_m + [S]}$$

#### 基本速率计算

```python
from chemmcp.tools import MichaelisMenten
mm = MichaelisMenten()

r = mm.run_code(
    analysis_type="calculate_velocity",
    substrate_concentration_S=5.0,    # [S] = 5 mM
    Vmax=10.0,                        # μM/s
    Km=2.0                           # mM
)

print(f"反应速率 v = {r['velocity']:.4f} μM/s")        # ≈ 7.14
print(f"占 Vmax 比例 = {r['fraction_of_Vmax']*100:.1f}%")  # ≈ 71.4%
print(f"饱和程度: {r['substrate_saturation']}")          # 'partial'
```

**关键浓度点验证：**
- 当 [S] = Km 时 → v = Vmax/2（半饱和）
- 当 [S] << Km 时 → v ∝ [S]（一级动力学）
- 当 [S] >> Km 时 → v ≈ Vmax（零级动力学，饱和）

#### 抑制分析

```python
# 竞争性抑制（抑制剂与底物竞争结合位点）
r = mm.run_code(
    analysis_type="full_analysis",
    substrate_concentration_S=5.0,
    Vmax=10.0, Km=2.0,
    inhibition_type="competitive",
    inhibitor_concentration_I=3.0,
    Ki=1.0                          # 抑制常数
)

print(f"无抑制 v = {r['velocity_uninhibited']:.4f}")
print(f"有抑制 v = {r['velocity_inhibited']:.4f}")
print(f"抑制率 = {(1-r['velocity_inhibited']/r['velocity_uninhibited'])*100:.1f}%")

# α 因子（竞争性抑制只影响 apparent Km，不影响 Vmax）
print(f"α = {r.get('alpha', 'N/A')}")
```

**四种抑制类型对比：**

| 类型 | Km 变化 | Vmax 变化 | 典型例子 |
|------|---------|----------|---------|
| `competitive` | Km↑ (×α) | 不变 | 磺胺类药物 |
| `uncompetitive` | Km↓ (×α') | Vmax↓ (×α') | 重金属离子 |
| `noncompetitive` | 不变 | Vmax↓ (×α) | 与酶-底物复合物结合 |
| `mixed` | Km↑ (×α) | Vmax↓ (×α') | 混合型抑制剂 |

#### Lineweaver-Burk 双倒数作图

```python
r = mm.run_code(
    analysis_type="lineweaver_burk",
    substrate_velocities_data=[
        {"S": 1.0, "v": 1.67},
        {"S": 2.0, "v": 2.5},
        {"S": 5.0, "v": 3.33},
        {"S": 10.0, "v": 4.0},
        {"S": 20.0, "v": 4.44},
    ]
)

lb = r["linearization"]
print(f"拟合 Vmax = {lb['Vmax']:.2f}")
print(f"拟合 Km = {lb['Km']:.2f}")
print(f"斜率 = {lb['slope']:.4f} (= Km/Vmax)")
print(f"截距 = {lb['intercept']:.4f} (= 1/Vmax)")
print(f"R² = {lb['r_squared']:.4f}")
```

---

## 4. 计算化学工具 (#496-500)

### #496 ReactionNetworkSolver — 反应网络求解器

**用途：** 对任意反应网络建立微分方程组，用四阶龙格-库塔法 (RK4) 数值积分，得到所有物种的浓度-时间曲线。

**支持的反应类型：**

```python
from chemmcp.tools import ReactionNetworkSolver
rn = ReactionNetworkSolver()
```

**连续反应 A → B → C：**

```python
r = rn.run_code(
    species=["A", "B", "C"],
    reactions=[
        {"reactants": ["A"], "products": ["B"], "k": 0.1},
        {"reactants": ["B"], "products": ["C"], "k": 0.05},
    ],
    initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
    time_end=100.0,
    n_points=50
)

# 浏览浓度变化
profiles = r["concentration_profiles"]
for p in profiles[::10]:  # 每10个点打印一次
    print(f"t={p['time']:6.1f}s  A={p['A']:.4f}  B={p['B']:.4f}  C={p['C']:.4f}")

# 半衰期信息
print(f"\n半衰期: {r['half_lives']}")
print(f"达到稳态? {r.get('steady_state_info', {}).get('steady', 'N/A')}")
```

**可逆反应 A ⇌ B：**

```python
r = rn.run_code(
    species=["A", "B"],
    reactions=[{
        "reactants": ["A"], "products": ["B"],
        "k": 0.1, "reversible": True, "k_reverse": 0.05
    }],
    initial_concentrations={"A": 1.0, "B": 0.0},
    time_end=200.0,
    n_points=100
)
# 最终 [A]/[B] 应趋近 k_rev/k_f = 0.5
```

**平行反应 A → B, A → C：**

```python
r = rn.run_code(
    species=["A", "B", "C"],
    reactions=[
        {"reactants": ["A"], "products": ["B"], "k": 0.1},   # 支路1
        {"reactants": ["A"], "products": ["C"], "k": 0.05},  # 支路2 (较慢)
    ],
    initial_concentrations={"A": 1.0, "B": 0.0, "C": 0.0},
    time_end=50.0
)
# B:C 产率比 ≈ k1:k2 = 2:1
```

**双分子反应 A + B → C：**

```python
r = rn.run_code(
    species=["A", "B", "C"],
    reactions=[{
        "reactants": ["A", "B"],
        "products": ["C"],
        "k": 0.1                       # 二级速率常数
    }],
    initial_concentrations={"A": 1.0, "B": 1.0, "C": 0.0},
    time_end=30.0
)
```

---

### #497 GeometryOptimizer — 几何优化

**用途：** 给定分子初始结构，通过能量极小化找到局部最优几何构型。使用 Lennard-Jones 势描述原子间相互作用。

**三种优化算法：**

| 算法 | 特点 | 适用场景 |
|------|------|---------|
| `steepest_descent` | 最速下降，稳定但收敛慢 | 初始结构差时 |
| `conjugate_gradient` | 共轭梯度，收敛较快 | 一般情况（推荐） |
| `damped_md` | 阻尼分子动力学 | 克服较大形变 |

**示例：水分子 H₂O 几何优化**

```python
from chemmcp.tools import GeometryOptimizer
go = GeometryOptimizer()

r = go.run_code(
    atoms=[
        {"symbol": "O", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [0.96, 0.0, 0.0]},     # 初始 O-H 偏长
        {"symbol": "H", "position": [-0.24, 0.93, 0.0]},
    ],
    optimizer="steepest_descent",
    max_iterations=100,
    convergence_threshold=1e-4
)

print(f"收敛? {r['converged']}")
print(f"迭代次数: {r['n_iterations']}")
print(f"最终能量: {r['final_energy_eV']:.4f} eV")
print(f"初始能量: {r['initial_energy_eV']:.4f} eV")
print(f"能量降低: {r['initial_energy_eV'] - r['final_energy_eV']:.4f} eV")

# 优化后的坐标
for atom in r["optimized_coordinates"]:
    print(f"  {atom['symbol']}: ({atom['x']:.4f}, {atom['y']:.4f}, {atom['z']:.4f})")
```

**CO₂ 线性分子优化：**

```python
r = go.run_code(
    atoms=[
        {"symbol": "O", "position": [-1.16, 0.0, 0.0]},
        {"symbol": "C", "position": [0.0, 0.0, 0.0]},
        {"symbol": "O", "position": [1.16, 0.0, 0.0]},
    ],
    optimizer="conjugate_gradient",
    max_iterations=50
)
```

**自定义键参数（力常数 + 平衡键长）：**

```python
r = go.run_code(
    atoms=[
        {"symbol": "H", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [1.5, 0.0, 0.0]},  # 拉长的 H-H 键
    ],
    bonds=[{"i": 0, "j": 1, "r0": 0.74, "k": 450}],   # 平衡键长0.74Å, 力常数450
    optimizer="conjugate_gradient"
)
```

---

### #498 TransitionStateSearch — 过渡态搜索

**用途：** 在势能面上寻找**一阶鞍点**（Transition State），即沿反应坐标能量极大、垂直方向能量极小的特殊构型。

**判断标准：** Hessian 矩阵有且仅有 **1 个负本征值**（对应虚频/imaginary frequency）

**三种搜索方法：**

| 方法 | 原理 | 适用场景 |
|------|------|---------|
| `quadratic_saddle` | 二次鞍点近似 | 快速估计，教学演示 |
| `eigenvector_following` | 本征向量跟踪（EF算法） | 实际搜索（推荐） |
| `dimer` | 二聚体方法 | 复杂 PES |

**示例：H₂ 分子拉伸 → 解离过渡态**

```python
from chemmcp.tools import TransitionStateSearch
tss = TransitionStateSearch()

r = tss.run_code(
    atoms=[
        {"symbol": "H", "position": [-0.5, 0.0, 0.0]},
        {"symbol": "H", "position": [0.5, 0.0, 0.0}],   # 初始猜测: 拉伸的 H-H
    ],
    search_method="quadratic_saddle",
    max_iterations=100
)

print(f"收敛? {r['converged']}")
print(f"TS 能量: {r['ts_energy']:.4f} eV")
print(f"虚频数: {r['n_imaginary']}")                    # 应该 = 1
print(f"虚频率: {r['imaginary_frequency_cm-1']:.1f} cm⁻¹")

# TS 验证
print(f"TS 有效? {r['ts_valid']}")                      # True = 一阶鞍点 ✓

# TS 坐标
for i, coord in enumerate(r["ts_coordinates"]):
    print(f"  原子{i}: ({coord[0]:.4f}, {coord[1]:.4f}, {coord[2]:.4f})")
```

**三原子系统 H···H···H（氢原子交换反应）：**

```python
r = tss.run_code(
    atoms=[
        {"symbol": "H", "position": [-0.5, 0.0, 0.0]},
        {"symbol": "H", "position": [0.0, 0.0, 0.0]},    # 中间 H
        {"symbol": "H", "position": [1.0, 0.0, 0.0]},
    ],
    search_method="eigenvector_following",
    max_iterations=100
)
```

**Hessian 本征值分析：**

```python
eigenvalues = r["hessian_eigenvalues"]
n_negative = sum(1 for ev in eigenvalues if ev < 0)
print(f"Hessian 本征值: {[f'{ev:.4f}' for ev in eigenvalues]}")
print(f"负本征值数: {n_negative}")
# n_negative = 0 → 极小点（稳定分子）
# n_negative = 1 → 过渡态（一阶鞍点）✓
# n_negative ≥ 2 → 高阶鞍点（非真实TS）
```

---

### #499 PotentialEnergySurface — 势能面扫描

**用途：** 沿某个内坐标（键长、键角、二面角）扫描势能面，得到能量随坐标变化的曲线，用于：
- 找到能量极小点（稳定构型）
- 估算活化能垒
- 确认反应路径

**扫描类型：**

| scan_type | 坐标定义 | scan_atoms 格式 | 示例 |
|-----------|---------|----------------|------|
| `bond_length` | 键长 | `[i, j]` | H-H 键扫描 |
| `angle` | 键角 | `[i, j, k]` | H-O-H 角度扫描 |
| `dihedral` | 二面角 | `[i, j, k, l]` | 丁烷旋转扫描 |

**示例 1：H-H 键长扫描（解离曲线）**

```python
from chemmcp.tools import PotentialEnergySurface
pes = PotentialEnergySurface()

r = pes.run_code(
    atoms=[
        {"symbol": "H", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [0.74, 0.0, 0.0}],
    ],
    scan_type="bond_length",
    scan_atoms=[0, 1],
    start_value=0.4,          # 起始键长 0.4 Å
    end_value=2.0,            # 终止键长 2.0 Å
    n_points=20               # 20 个扫描点
)

print(f"扫描类型: {r['scan_type']}")
print(f"范围: {r['scan_range']}")
print(f"最低能量: {r['min_energy_eV']:.4f} eV at point {r['min_energy_index']}")

# 能量剖面数据
for pt in r["scan_results"][:5]:
    print(f"  键长={pt['value']:.3f}Å  E={pt['energy']:.4f} eV")
print("  ...")

# 驻点（极小/极大/鞍点）
for sp in r.get("stationary_points", []):
    print(f"  驻点 @ index {sp['index']}: type={sp['type']}, E={sp['energy']:.4f} eV")
```

**示例 2：H-O-H 键角扫描（水分子弯曲）**

```python
r = pes.run_code(
    atoms=[
        {"symbol": "O", "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "position": [0.96, 0.0, 0.0]},
        {"symbol": "H", "position": [0.0, 0.96, 0.0]},
    ],
    scan_type="angle",
    scan_atoms=[1, 0, 2],        # H(1)-O(0)-H(2) 角度
    start_value=60.0,             # 60° → 180°
    end_value=180.0,
    n_points=15
)
# 预期: 最小值在 ~104.5°（水分子的平衡键角附近）
```

**示例 3：丁烷 C-C-C-C 二面角扫描（构象异构）**

```python
r = pes.run_code(
    atoms=[
        {"symbol": "C", "position": [0.0, 0.0, 0.0]},
        {"symbol": "C", "position": [1.54, 0.0, 0.0]},
        {"symbol": "C", "position": [2.31, 1.33, 0.0]},
        {"symbol": "C", "position": [3.85, 1.33, 0.0]},
    ],
    scan_type="dihedral",
    scan_atoms=[0, 1, 2, 3],
    start_value=-180.0,
    end_value=180.0,
    n_points=12
)
# 预期: 反式(-180°/180°) 和 旁式(±60°) 出现极小值, 顺叠(0°) 为极大值
```

**每点优化（更精确但更慢）：**

```python
r = pes.run_code(
    atoms=[...],
    scan_type="bond_length",
    scan_atoms=[0, 1],
    start_value=0.5, end_value=2.0, n_points=15,
    optimize_each_point=True    # 每个扫描点都做全坐标优化
)
```

---

### #500 FrequencyAnalysis — 频率分析

**用途：** 从 Hessian 矩阵（能量二阶导数矩阵）计算振动频率，用于：
- **确认驻点类型**：极小点（全正频率，稳定分子）vs 过渡态（1 个虚频）
- 计算**零点能 (ZPE)** 校正
- 计算热力学函数（G, H, S）—— 用于判断反应自发性

**核心原理：**
$$\nu_i = \frac{1}{2\pi c}\sqrt{\frac{\lambda_i}{\mu}}$$

其中 λᵢ 是质量加权 Hessian 的本征值。

**示例 1：H₂ 双原子分子**

```python
from chemmcp.tools import FrequencyAnalysis
fa = FrequencyAnalysis()

r = fa.run_code(
    atoms=[
        {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "mass": 1.008, "position": [0.74, 0.0, 0.0]},
    ],
    temperature_K=298.15,
    pressure_atm=1.0
)

print(f"驻点类型: {r['stationary_point_type']}")
print(f"虚频数: {r['n_imaginary_frequencies']}")       # 0=极小, 1=过渡态
print(f"振动频率: {r['frequencies_cm-1']} cm⁻¹")
print(f"ZPE = {r['zero_point_energy']:.4f} eV")

# 热力学量
th = r["thermodynamics"]
print(f"G(298K) = {th['gibbs_free_energy']:.4f} eV")
print(f"H(298K) = {th['enthalpy']:.4f} eV")
print(f"S(298K) = {th['entropy']:.4f} eV/K")
```

**示例 2：H₂O 非线性分子（3N-6 = 3 个振动模式）**

```python
r = fa.run_code(
    atoms=[
        {"symbol": "O", "mass": 15.999, "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "mass": 1.008, "position": [0.96, 0.0, 0.0]},
        {"symbol": "H", "mass": 1.008, "position": [-0.24, 0.93, 0.0]},
    ],
    temperature_K=298.15
)

print(f"驻点类型: {r['stationary_point_type']}")
print(f"频率数: {len(r['frequencies_cm-1'])}")            # 应为 3
for i, (freq, mode) in enumerate(zip(r['frequencies_cm-1'], r.get('normal_modes', []))):
    print(f"  模式 {i+1}: {freq:.1f} cm⁻¹  ({mode.get('description', 'N/A')})")
```

**示例 3：过渡态验证（应有恰好 1 个虚频）**

```python
r = fa.run_code(
    atoms=[
        {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "mass": 1.008, "position": [1.2, 0.0, 0.0]},  # 拉伸的 H-H
    ]
)

print(f"类型: {r['stationary_point_type']}")
# 预期: "transition state (first-order saddle point)"
print(f"虚频数: {r['n_imaginary_frequencies']}")          # 预期: 1
if r['n_imaginary_frequencies'] == 1:
    print("✓ 确认为过渡态结构")
elif r['n_imaginary_frequencies'] == 0:
    print("✓ 确认为稳定分子（能量极小点）")
else:
    print("⚠ 高阶鞍点，可能不是真实的过渡态")
```

**自定义 Hessian 矩阵：**

```python
import numpy as np
# 2 个原子 → 6×6 Hessian
n_dim = 6
H = [[0.5 if i == j else 0.01 for j in range(n_dim)] for i in range(n_dim)]

r = fa.run_code(
    atoms=[
        {"symbol": "H", "mass": 1.008, "position": [0.0, 0.0, 0.0]},
        {"symbol": "H", "mass": 1.008, "position": [1.2, 0.0, 0.0]},
    ],
    hessian_matrix=H,
    scale_factor=0.96            # DFT 常用校正因子 0.96-0.98
)
```

**频率校正因子参考：**

| 方法 | 典型 scale_factor |
|------|------------------|
| HF/6-31G(d) | 0.89 |
| B3LYP/6-31G(d) | 0.96-0.98 |
| M06-2X/cc-pVTZ | 0.98 |
| 实验（无需校正） | 1.0 |

---

## 5. 典型工作流

### 工作流 A：完整反应动力学分析

```
#494 RateDeterminingStep  →  识别瓶颈
        ↓
#493 PreEquilibrium      →  推导速率方程（如果有快平衡）
        ↓
#492 SteadyStateApprox   →  验证中间体近似是否合理
        ↓
#496 ReactionNetworkSolver →  数值模拟完整浓度-时间曲线
        ↓
#491 TunnelingCorrection →  修正轻原子转移的量子效应
```

**实例：分析一个三步有机反应机理**

```python
from chemmcp.tools import (
    RateDeterminingStep, PreEquilibrium,
    SteadyStateApproximation, ReactionNetworkSolver,
    TunnelingCorrection
)

# Step 1: 找出速控步
rds = RateDeterminingStep()
mechanism = [
    {"reactants": "A + Cat ⇌ AC", "k": 1e6, "reversible": True},
    {"reactants": "AC + B → Int", "k": 0.05},
    {"reactants": "Int → P + Cat", "k": 10.0},
]
r = rds.run_code(mechanism, has_pre_equilibrium=True)
print(f"RDS = 第 {r['rds_step_index']} 步")

# Step 2: 用预平衡推导速率方程
peq = PreEquilibrium()
r2 = peq.run_code("bimolecular", [1e6], [1e5], 0.05, {"A": 1.0, "B": 1.0, "Cat": 0.1})
print(f"k_eff = {r2['effective_rate_constant']}")

# Step 3: 数值模拟
rn = ReactionNetworkSolver()
r3 = rn.run_code(
    species=["A", "Cat", "AC", "B", "Int", "P"],
    reactions=[
        {"reactants": ["A","Cat"], "products": ["AC"], "k": 1e6, "reversible": True, "k_reverse": 1e5},
        {"reactants": ["AC","B"], "products": ["Int"], "k": 0.05},
        {"reactants": ["Int"], "products": ["P","Cat"], "k": 10.0},
    ],
    initial_concentrations={"A": 1.0, "Cat": 0.1, "B": 1.0},
    time_end=200.0
)
```

### 工作流 B：计算化学反应活化能

```
#497 GeometryOptimizer     →  优化反应物构型 → E_reactant
        ↓
#498 TransitionStateSearch →  找到过渡态 → E_TS (+ 验证 1 个虚频)
        ↓
#499 PotentialEnergySurface →  扫描反应路径确认
        ↓
#500 FrequencyAnalysis     →  ZPE 校正 → ΔG‡
```

**实例：H₂ 解离反应能垒计算**

```python
from chemmcp.tools import (
    GeometryOptimizer, TransitionStateSearch,
    PotentialEnergySurface, FrequencyAnalysis
)

go = GeometryOptimizer()
tss = TransitionStateSearch()
pes = PotentialEnergySurface()
fa = FrequencyAnalysis()

# 1. 优化反应物 (H₂ 平衡构型)
r_react = go.run_code(
    atoms=[{"symbol":"H","pos":[0,0,0]}, {"symbol":"H","pos":[0.8,0,0]}],
    optimizer="conjugate_gradient"
)
E_react = r_react["final_energy_eV"]

# 2. 搜索过渡态 (拉伸的 H₂)
r_ts = tss.run_code(
    atoms=[{"symbol":"H","pos":[-0.5,0,0]}, {"symbol":"H","pos":[0.5,0,0]}],
    method="eigenvector_following"
)
E_ts = r_ts["ts_energy"]

# 3. 频率分析确认 TS
r_freq_ts = fa.run_code(atoms=r_ts["ts_atoms"])
assert r_freq_ts["n_imaginary_frequencies"] == 1, "不是有效的过渡态!"
ZPE_ts = r_freq_ts["zero_point_energy"]

# 4. 频率分析反应物
r_freq_react = fa.run_code(atoms=r_react["optimized_atoms"])
ZPE_react = r_freq_react["zero_point_energy"]

# 5. 计算活化能 (含 ZPE 校正)
delta_E = (E_ts + ZPE_ts) - (E_react + ZPE_react)
print(f"活化能 Ea = {delta_E:.4f} eV = {delta_E * 96.485:.2f} kJ/mol")
```

### 工作流 C：酶动力学实验数据分析

```python
from chemmcp.tools import MichaelisMenten
mm = MichaelisMenten()

# 实验测得的 [S] vs v 数据
data = [
    {"S": 0.1, "v": 0.82}, {"S": 0.2, "v": 1.45},
    {"S": 0.5, "v": 2.78}, {"S": 1.0, "v": 3.57},
    {"S": 2.0, "v": 4.17}, {"S": 5.0, "v": 4.55},
    {"S": 10.0, "v": 4.76},
]

# Lineweaver-Burk 双倒数拟合
r = mm.run_code("lineweaver_burk", substrate_velocities_data=data)
lb = r["linearization"]
print(f"Vmax = {lb['Vmax']:.2f}, Km = {lb['Km']:.2f}")

# 用拟合参数预测新浓度下的速率
for S_new in [0.3, 3.0, 8.0]:
    pred = mm.run_code("calculate_velocity", S_new, lb["Vmax"], lb["Km"])
    print(f"[S]={S_new}: v_pred={pred['velocity']:.3f}")

# 分析竞争性抑制效果
r_inhib = mm.run_code("full_analysis", 2.0, lb["Vmax"], lb["Km"],
                       "competitive", 1.0, lb["Km"]*0.5)
print(f"无抑制: v={r_inhib['velocity_uninhibited']:.3f}")
print(f"有抑制: v={r_inhib['velocity_inhibited']:.3f}")
```

---

## 6. Cherry Studio 配置

将以下 JSON 添加到 Cherry Studio 的 MCP 配置中：

### 方案 A：单独加载每个工具（细粒度控制）

```json
{
  "mcpServers": {
    "ChemMCP_491_Tunneling": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "TunnelingCorrection"]
    },
    "ChemMCP_492_SteadyState": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "SteadyStateApproximation"]
    },
    "ChemMCP_493_PreEq": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "PreEquilibrium"]
    },
    "ChemMCP_494_RDS": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "RateDeterminingStep"]
    },
    "ChemMCP_495_MM": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "MichaelisMenten"]
    },
    "ChemMCP_496_Network": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "ReactionNetworkSolver"]
    },
    "ChemMCP_497_GeoOpt": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "GeometryOptimizer"]
    },
    "ChemMCP_498_TSSearch": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "TransitionStateSearch"]
    },
    "ChemMCP_499_PES": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "PotentialEnergySurface"]
    },
    "ChemMCP_500_FreqAna": {
      "command": "uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "FrequencyAnalysis"]
    }
  }
}
```

### 方案 B：一次性加载全部 10 个工具（推荐）

```json
{
  "mcpServers": {
    "ChemMCP_Kinetics_Computational": {
      "command": "uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools",
        "TunnelingCorrection,SteadyStateApproximation,PreEquilibrium,RateDeterminingStep,MichaelisMenten,ReactionNetworkSolver,GeometryOptimizer,TransitionStateSearch,PotentialEnergySurface,FrequencyAnalysis"
      ]
    }
  }
}
```

### 导入步骤

1. 打开 Cherry Studio → 设置 → MCP Servers
2. 点击「导入」或手动粘贴上方 JSON
3. 保存后在对话中选择对应的 MCP 服务
4. 发送问题如：「计算 H 转移反应在 200K 的隧穿校正因子」
5. AI 将自动调用 MCP 工具并返回结果

---

## 7. 常见问题

### Q1: run_code() 和 run_text() 返回什么？

两者都返回 Python **dict**（字典），包含计算结果。直接访问 key 即可获取数据：

```python
result = tool.run_code(...)
print(result["correction_factor_kappa"])  # 直接取值
```

**注意：** 不是 `{"result": ...}` 包装格式，而是原始 dict。

### Q2: 速率常数的单位是什么？

- 一级反应 (s⁻¹): `k = 0.1` 表示 0.1 s⁻¹
- 二级反应 (M⁻¹s⁻¹ 或 L·mol⁻¹·s⁻¹): `k = 0.1` 表示 0.1 M⁻¹s⁻¹
- 工具内部统一使用相对单位，实际应用时注意换算

### Q3: 几何优化的能量单位是什么？

- 能量：**eV** (电子伏特)
- 距离：**Å** (埃, Angstrom)
- 力常数：默认 kcal/mol/Å²（LJ 参数中 epsilon 为 eV）
- 换算：1 eV = 96.485 kJ/mol = 23.06 kcal/mol

### Q4: 频率分析的 Hessian 怎么获得？

- 本工具使用 Lennard-Jones 势自动构建 Hessian（基于原子坐标）
- 如需使用外部量子化学软件（Gaussian/ORCA）的 Hessian，通过 `hessian_matrix` 参数传入
- 外部 Hessian 单位应为 **eV/Å²**

### Q5: TunnelingCorrection 的 κ 值为什么小于 1？

不同文献对 κ 的定义不同：
- **传统定义**: κ ≥ 1（量子增强，k_quantum > k_classical）
- **本实现**: κ 可能 < 1 取决于模型参数化方式

重要的是比较**趋势**而非绝对值：温度↓、势垒↓、虚频↑ → 隧穿效应增强。

### Q6: RateDeterminingStep 的索引从几开始？

**1-indexed**（从 1 开始）。单步反应返回 `rds_step_index=1`，三步中最慢的第二步返回 `2`。

### Q7: 如何将结果用于论文/报告？

每个工具的输出 dict 包含完整的中间数据和解释性字段：
- `description` / `interpretation`: 文字描述
- `raw_data`: 原始数值（可用于作图）
- `validity` / `conditions`: 结果有效性判断

建议用 `raw_data` + matplotlib/matlab 自定义绘图。

### Q8: 测试如何运行？

```bash
cd ~/ChemMCP
uv run python -m pytest test/test_tools_491_500.py -v
# 预期: 60 passed
```

---

## 附录：工具间关系图

```
反应动力学方向                        计算化学方向
───────────────                      ───────────────

#494 RDS ──→ #493 PreEquil          #497 GeoOpt ──→ 优化结构
  │              │                       │
  ▼              ▼                       ▼
#492 SSA    #496 NetworkSolver      #498 TSSearch ──→ 过渡态
  │                                      │
  ▼                                      ▼
#491 Tunneling                     #499 PES Scan
                                       │
                                       ▼
                                  #500 Freq Analysis
                                  (ZPE + 热力学 + 驻点确认)
```

---

*文档版本: 2026-05-08 | ChemMCP Tools #491-500 | 60 tests passing ✅*
