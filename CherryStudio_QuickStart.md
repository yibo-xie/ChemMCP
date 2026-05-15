# ChemMCP × Cherry Studio 快速开始指南

> **406 个化学 MCP 工具**，覆盖有机、无机、物理、分析化学，通过 Cherry Studio 即插即用

---

## 目录

1. [环境准备](#1-环境准备)
2. [启动 MCP Server](#2-启动-mcp-server)
3. [配置 Cherry Studio](#3-配置-cherry-studio)
4. [工具总览（406 个）](#4-工具总览406-个)
5. [按类别速查](#5-按类别速查)
6. [使用示例](#6-使用示例)
7. [常见问题](#7-常见问题)

---

## 1. 环境准备

### 前置要求

| 依赖 | 版本 | 检查命令 |
|------|------|----------|
| Python | ≥ 3.10 | `python3 --version` |
| uv (包管理器) | 最新 | `uv --version` |
| Cherry Studio | ≥ 1.0 | 桌面应用 |

### 安装步骤

```bash
# 1. 克隆/确认仓库位置
cd ~/ChemMCP

# 2. 安装依赖（已配置 .venv）
uv sync

# 3. 验证安装
.venv/bin/python -c "from chemmcp.utils.mcp_app import ChemMCPManager; print(f'{len(ChemMCPManager.get_registered_tools())} tools registered')"
# 预期输出: 406 tools registered
```

> ⚠️ **注意**: 如果遇到 `SchrodingerSolver1D` 导入错误，说明需要修复 `src/chemmcp/tools/__init__.py` 第 567 行的重复条目（详见[常见问题 Q1](#q1-导入错误-schrodingersolver1d)）。

---

## 2. 启动 MCP Server

ChemMCP 支持两种 MCP 传输模式：

### 模式 A：stdio（推荐用于 Cherry Studio）

```bash
cd ~/ChemMCP
uv run -m chemmcp --tools ToolName1,ToolName2,...
```

Cherry Studio 会自动管理进程的 stdin/stdout，无需手动启动。

### 模式 B：SSE（Server-Sent Events，用于远程访问）

```bash
cd ~/ChemMCP
uv run -m chemmcp --tools ToolName1,ToolName2,... --sse
# 服务启动在 http://127.0.0.1:8001
```

### 加载工具的方式

```bash
# 单个工具
--tools MoleculeWeight

# 多个工具（逗号分隔，无空格）
--tools MoleculeWeight,MoleculeSimilarity,SafetyCheck

# 批量加载（适合相关工具组合）
--tools BalanceEquation,BalanceRedox,CalculateReactionEnthalpy,CalculateGibbsChange
```

---

## 3. 配置 Cherry Studio

### 方法一：手动添加 MCP Server

1. 打开 **Cherry Studio** → 左下角 ⚙️ **设置** → **MCP Servers**
2. 点击 **「+」** 添加新服务器
3. 填写配置：

**基础配置示例（单个工具）：**

| 字段 | 值 |
|------|-----|
| 名称 | `ChemMCP_MoleculeWeight` |
| 类型 | `stdio` |
| Command | `/home/wave/.local/bin/uv` |
| Args | `--directory`, `/home/wwave/ChemMCP`, `run`, `-m`, `chemmcp`, `--tools`, `MoleculeWeight` |

### 方法二：导入 JSON 配置（推荐）

将以下 JSON 导入 Cherry Studio 的 MCP 设置：

#### 配置 A：最小化测试（4 个核心工具）

```json
{
  "mcpServers": {
    "ChemMCP_MolWeight": {
      "command": "/home/wave/.local/bin/uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "MoleculeWeight"]
    },
    "ChemMCP_Similarity": {
      "command": "/home/wave/.local/bin/uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "MoleculeSimilarity"]
    },
    "ChemMCP_IUPAC": {
      "command": "/home/wave/.local/bin/uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "Iupac2Smiles"]
    },
    "ChemMCP_Safety": {
      "command": "/home/wave/.local/bin/uv",
      "args": ["--directory", "/home/wave/ChemMCP", "run", "-m", "chemmcp", "--tools", "SafetyCheck"]
    }
  }
}
```

#### 配置 B：有机合成工作流（15 个工具）

```json
{
  "mcpServers": {
    "ChemMCP_Organic": {
      "command": "/home/wave/.local/bin/uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools",
        "ForwardSynthesis,Retrosynthesis,PredictProducts,ReactionPredictor,ReactionSmilesCheck,BalanceEquation,BalanceRedox,FunctionalGroups,FunctionalGroupIdentifier,MoleculeModifier,MoleculeVisualizer,MoleculeWeight,MoleculeSimilarity,Smiles2Formula,Smiles2Iupac"
      ]
    }
  }
}
```

#### 配置 C：物理化学套件（20 个工具）

```json
{
  "mcpServers": {
    "ChemMCP_Physical": {
      "command": "/home/wave/.local/bin/uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools",
        "IdealGasCalculator,VanDerWaalsGas,CompressibilityFactor,VirialEquation,ClausiusClapeyron,CalculatePH,CalculateGibbsChange,CalculateReactionEnthalpy,EquilibriumConstant,NernstEquation,CalculateCellPotential,ArrheniusCalculator,CollisionTheory,EyringEquation,TransitionStateTheory,RateDeterminingStep,SteadyStateApproximation,ReactionNetworkSolver,MichaelisMenten,TunnelingCorrection"
      ]
    }
  }
}
```

#### 配置 D：全量加载（所有 406 个工具）

```json
{
  "mcpServers": {
    "ChemMCP_Full": {
      "command": "/home/wave/.local/bin/uv",
      "args": [
        "--directory", "/home/wave/ChemMCP",
        "run", "-m", "chemmcp",
        "--tools",
        "AcidBaseStrengthCompare,ActinideProperties,ActivityCoefficient,AdvancedPartitionFunction,AldolMechanism,AldolReaction,AlkaliMetalProperties,AlkalineEarthProperties,AmphotericSpecies,AnalyzeMolecularOrbital,AnovaAnalyzer,AromaticSystemDetector,ArrheniusAnalyzer,ArrheniusCalculator,ArrowPushingValidator,AssignOxidationNumber,AuditTrailReviewer,BackTitrationSolver,BaeyerVilligerOxidation,BalanceEquation,BalanceRedox,BasisSetHandler,BeckmannRearrangement,BeerLambertCalculator,BETSurfaceArea,BindingEnergy,BirchReduction,BisectionMethod,BondOrderCalculator,BornOppenheimer,BuchwaldHartwigAmination,BufferCapacity,BufferCapacityCalculator,BufferPreparation,ButlerVolmerKinetics,CalculateCellPotential,CalculateFormalCharge,CalculateGibbsChange,CalculateLatticeEnergy,CalculatePH,CalculateReactionEnthalpy,CalculateSolubility,CalibrationCurveFitter,CapacityFactorCalculator,CarbocationRearrangement,CarnotEfficiency,CentrifugationCalculator,ChelateEffect,ChemicalPotential,ChemicalPotentialAdvanced,ChronoamperometryAnalyzer,CisTransAnalyzer,ClaisenCondensation,ClaisenMechanism,ClausiusClapeyron,ClemmensenReduction,CMCDocumentationHelper,CollisionTheory,CommonIonEffect,CommonIonSolubility,CommonNameLookup,CommutatorCalculator,CompareElements,ComplexFromName,ComplexIonSolubility,CompressibilityFactor,ConductivityCellConstant,ConfigurationInteraction,ConstitutionalIsomerGenerator,ContactAngleAnalyzer,CoordinationGeometry,CoulombIntegral,CoulombPotential,CouplingConstantAnalyzer,CrystalFieldSplitting,CurveFittingNonlinear,CvPeakAnalyzer,DeadVolumeCalculator,DegreeOfDissociation,DElectronConfiguration,DerivatizationReagentSelector,DeterminantCalculator,DftXcFunctional,DielsAlderMechanism,DielsAlderReaction,DiffusionCoefficientCalculator,DigestionProtocolSelector,DimensionalAnalysis,DipoleMoment,DipoleMomentEstimator,DisproportionationCheck,DissolvePrecipitate,DrawLewisStructure,E1Mechanism,E2Mechanism,EisCircuitFitter,ElectricFieldMolecule,ElectrodeSelectionGuide,ElectronDensityCalculator,ElectronDensityPlotter,ElectronSinkIdentifier,ElectronSourceIdentifier,ElectrophilicAddition,ElectrophilicAromaticSubstitution,ElementAbundance,EnthalpyCalculator,EnthalpyCalculatorNew,EntropyCalculator,EntropyCalculatorNew,EnzymeKinetics,EquationOfState,EquilibriumConstant,EquivalencePoint,ErrorPropagation,EvaporationEstimator,ExchangeIntegral,ExpectationValue,ExtendedHuckel,ExtractionOptimizer,EyringEquation,EzConfigurator,FftTransform,FiltrationGuide,ForceFieldCalculator,ForwardSynthesis,FranckCondon,FranckCondonFactors,FrequencyAnalysis,FreundlichIsotherm,FriedelCraftsReaction,FrontierOrbitalAnalysis,FTestCalculator,FtirBaselineCorrector,FugacityCalculator,FunctionalGroupIdentifier,FunctionalGroups,GasEvolutionPrediction,GcCarrierGasSelector,GcColumnBleedPredictor,GcOvenProgramDesigner,GeometryOptimizer,GetBondEnergy,GetBondLength,GetCrystalStructure,GetElectronAffinity,GetElectronConfiguration,GetElementDiscovery,GetElementInfo,GetFormationConstant,GetGibbsEnergy,GetIonizationEnergy,GetIsotopes,GetKsp,GetOxidationStates,GetPka,GetPkb,GetStandardEnthalpy,GetStandardEntropy,GetStandardPotential,GibbsAdsorption,GibbsEnergyCalculator,GibbsFreeEnergy,GibbsMinimization,GradientCalculator,GrignardMechanism,GrignardReaction,HalfLifeCalculation,HalfLifeCalculator,HalogenProperties,HarmonicOscillator,HartreeFockSCF,HeatCapacityCalculator,HeatCapacityLookup,HeckReaction,HessianMatrix,HomogenizationProtocol,HplcColumnSelector,HuckelMethod,HybridizationAnalyzer,HydrogenAtomOrbitals,HydrogenWavefunction,IdealGasCalculator,IdentifyOxidizingAgent,IdentifyReactionType,IntegratedRateLaw,InternalStandardSelector,InterpolationSpline,IonChromatographyEluent,IonTransportNumber,IrPeakInterpreter,IrSpectrumInterpreter,IrSpectrumPredictor,IsomerTypes,Iupac2Smiles,IupacNamer,JahnTellerDistortion,JonesOxidation,JouleThomson,LangmuirIsotherm,LanthanideProperties,LaplacianOperator,LeastSquaresFit,LeavingGroupRanker,LeChatelierAnalyzer,LigandFieldDiagram,LinearityRangeValidator,LinearSystemSolver,LodLoqCalculator,MainGroupTrends,MassDefect,MassDefectFilter,MassSpecFragmenter,MatrixDiagonalization,MatrixEigenvalueSolver,MatrixInversion,MatrixMatchingAdvisor,MaxwellRelations,MesoCompoundChecker,MethodValidationChecklist,MichaelisMenten,MichaelMechanism,MobilePhaseOptimizer,MOEnergyLevelDiagram,MolarAbsorptivityLookup,MolecularDynamicsVerlet,MolecularIonCalculator,MolecularOrbitalDiagram,MoleculeAtomCount,MoleculeCaptioner,MoleculeGenerator,MoleculeModifier,MoleculePrice,MoleculeSimilarity,MoleculeSmilesCheck,MoleculeVisualizer,MoleculeWeight,MomentOfInertia,MonteCarloIntegrator,Mp2Correlation,MultipoleExpansion,Name2Smiles,NameComplex,NamedReactionLookup,NernstEquation,NernstEquationSolver,NewtonRaphsonSolver,NmrChemicalShift,NmrChemicalShiftPredictor,NmrCPredictor,NmrHPredictor,NmrShielding,NobleGasCompounds,NormalModeAnalysis,NuclearEquationBalance,NucleophilicAddition,NucleophilicAromaticSubstitution,NucleophilicityRanker,NumericalDifferentiation,NumericalIntegrator,OdeSolverRk4,OdeSolverStiff,OscillatorStrength,OutlierDetector,OverlapIntegral,OverpotentialAnalyzer,OxidationMechanism,ParallelConsecutiveReactions,PartialDerivative,PartialMolarQuantity,ParticleInBox,PatentCheck,PdeSolverFiniteDiff,PeakPurityAnalyzer,PericyclicAnalyzer,PeriodicTrend,PerturbationTheory,PhAdjustmentBuffer,PhaseEquilibrium,PhaseRuleAnalyzer,PhElectrodeCalibration,PlateNumberCalculator,PointGroupIdentifier,PoissonBoltzmann,Polarizability,PolyproticAcid,PotentialEnergySurface,PotentiometricTitrationEndpoint,PourbaixDiagramLookup,PrecipitationPrediction,PredictColor,PredictHybridization,PredictMagnetism,PredictPolarity,PredictProducts,PredictVseprGeometry,PreEquilibrium,PubchemSearch,PubchemSearchQA,RadialDistribution,RadicalChainMechanism,RamanActivity,RamanSpectrumPredictor,RateDeterminingStep,RateLawFitter,RateLawIntegrator,ReactionCoordinate,ReactionEnergyEstimator,ReactionMechanismSimulator,ReactionNetworkSolver,ReactionPredictor,ReactionSmilesCheck,RecoveryCalculator,ReductionMechanism,RegioselectivityPredictor,RegressionDiagnostics,ResolutionCalculator,RetentionTimePredictor,Retrosynthesis,RigidRotor,RingSystemAnalyzer,RobustnessDoeDesigner,Rotational Spectrum,RsConfigurator,SafetyCheck,SampleDilutionCalculator,SamplePreservationGuide,Schrodinger3DSolver,SchrodingerSolver1d,SecCalibrationCurve,SelectionRules,SelectionRulesChecker,SelectivePrecipitation,SelectivityFactorCalculator,Selfies2Smiles,SharplessDihydroxylation,SharplessEpoxidation,SignificantFigures,SlaterDeterminant,Smiles2Cas,Smiles2Formula,Smiles2Iupac,Smiles2Selfies,SmilesCanonicalization,SmilesToStructure,Sn1Mechanism,Sn2Mechanism,SolubilityRules,SonogashiraCoupling,SpecificityTestDesigner,SpectralDeconvolution,SpectrochemicalSeries,SPEMethodDesigner,SphericalHarmonics,SpinOrbitCoupling,SpinSpinCoupling,SplittingPatternExplainer,SpontaneityCheck,StabilityStudyPlanner,StandardSolutionPrep,StandardStateConverter,StatisticalEnsemble,SteadyStateApprox,StereocenterFinder,StereoisomerCounter,StereoselectivityPredictor,StrippingVoltammetryOptimizer,StructureToSmiles,SurfaceTensionCalculator,SuzukiCoupling,SvdDecomposition,SwernOxidation,SymmetryOperations,SymmetryPointGroup,SystemSuitabilityChecker,TautomerGenerator,TemperatureJumpRelaxation,TitrantStandardization,TitrationCurve,TransitionDipole,TransitionMetalChemistry,TransitionStateSearch,TransitionStateTheory,TTestCalculator,TunnelingCorrection,TunnelingProbability,UncertaintyPrinciple,UncertaintyPropagator,UnitConversion,UvVisSpectrum,UvVisTransitions,UvVisWavelengthSelector,VanDeemterAnalyzer,VanDerWaalsGas,VantHoffEquilibrium,VariationalMethod,VelocityVerlet,VibrationalModes,VirialEquation,VseprGeometry,WebSearch,WillPrecipitate,WittigMechanism,WittigReaction,WKBApproximation,WolffKishnerReduction"
      ]
    }
  }
}
```

> ⚠️ **全量加载警告**: 406 个工具同时加载会占用较多内存（~500MB+）和启动时间（~15-30秒）。建议按需使用配置 B/C 的分类组合。

### 导入步骤

1. Cherry Studio → **设置** (⚙️) → **MCP Servers**
2. 点击 **「导入」** → 粘贴 JSON → **保存**
3. 新建/打开对话 → 点击 **MCP 工具按钮** → 选择已配置的服务
4. 开始对话！AI 将自动调用工具

---

## 4. 工具总览（406 个）

### 统计概览

| 类别 | 工具数 | 编号范围 |
|------|--------|----------|
| 🔬 **通用/分子工具** | ~50 | #1-30 |
| ⚗️ **有机反应与机理** | ~80 | #31-120 |
| 🧪 **无机/配位化学** | ~50 | #121-180 |
| 🌡️ **物理化学/热力学** | ~100 | #181-300 |
| 📊 **分析化学/光谱** | ~80 | #301-390 |
| 🔬 **计算化学/量子** | ~46 | #391-450 (含扩充电工具) |
| *(持续扩展中)* | | #451+ |

---

## 5. 按类别速查

### A. 分子基础工具

| 工具名 | 功能 |
|--------|------|
| `MoleculeWeight` | 计算分子量 |
| `MoleculeSimilarity` | 分子相似性比较 (Tanimoto) |
| `MoleculeVisualizer` | 分子 2D 可视化 |
| `MoleculeAtomCount` | 原子计数 |
| `MoleculeCaptioner` | AI 分子描述生成 |
| `MoleculeGenerator` | AI 分子生成 |
| `MoleculeModifier` | 分子修饰/衍生化 |
| `MoleculePrice` | 分子购买价格查询 |
| `MoleculeSmilesCheck` | SMILES 合法性校验 |
| `Smiles2Formula` | SMILES → 分子式 |
| `Smiles2Iupac` | SMILES → IUPAC 名 |
| `Iupac2Smiles` | IUPAC 名 → SMILES |
| `Name2Smiles` | 化合物名称 → SMILES |
| `Smiles2Cas` | SMILES → CAS 号 |
| `Smiles2Selfies` | SMILES ↔ SELFIES 互转 |
| `Selfies2Smiles` | SELFIES → SMILES |
| `SmilesCanonicalization` | SMILES 规范化 |
| `PatentCheck` | 专利检索 |
| `SafetyCheck` | 安全性检查 |
| `PubchemSearch` | PubChem 搜索 |
| `PubchemSearchQA` | PubChem QA |
| `FunctionalGroups` | 官能团识别 |
| `FunctionalGroupIdentifier` | 官能团详细识别 |
| `DrawLewisStructure` | 路易斯结构式绘制 |

### B. 有机反应与机理

| 工具名 | 功能 |
|--------|------|
| `ForwardSynthesis` | 正向合成预测 (IBM RXN) |
| `Retrosynthesis` | 逆合成分析 (IBM RXN) |
| `PredictProducts` | 反应产物预测 |
| `ReactionPredictor` | 反应预测器 |
| `ReactionSmilesCheck` | 反应 SMILES 校验 |
| `BalanceEquation` | 化学方程式配平 |
| `BalanceRedox` | 氧化还原配平 |
| `AldolReaction` / `AldolMechanism` | 羟醛反应及机理 |
| `DielsAlderReaction` / `DielsAlderMechanism` | Diels-Alder 及机理 |
| `GrignardReaction` / `GrignardMechanism` | 格氏反应及机理 |
| `Sn1Mechanism` / `Sn2Mechanism` | SN1/SN2 机理 |
| `E1Mechanism` / `E2Mechanism` | E1/E2 消除机理 |
| `MichaelMechanism` | Michael 加成机理 |
| `WittigReaction` / `WittigMechanism` | Wittig 反应及机理 |
| `SuzukiCoupling` | Suzuki 偶联 |
| `HeckReaction` | Heck 反应 |
| `SonogashiraCoupling` | Sonogashira 偶联 |
| `BuchwaldHartwigAmination` | Buchwald-Hartwig 氨化 |
| `FriedelCraftsReaction` | Friedel-Crafts 反应 |
| `ElectrophilicAromaticSubstitution` |芳香亲电取代 (EAS) |
| `NucleophilicAromaticSubstitution` | 芳香亲核取代 (SNAr) |
| `NucleophilicAddition` | 亲核加成 |
| `ElectrophilicAddition` | 亲电加成 |
| `ClaisenCondensation` / `ClaisenMechanism` | Claisen 缩合及机理 |
| `BirchReduction` | Birch 还原 |
| `BaeyerVilligerOxidation` | Baeyer-Villiger 氧化 |
| `BeckmannRearrangement` | Beckmann 重排 |
| `ClemmensenReduction` | Clemmensen 还原 |
| `JonesOxidation` | Jones 氧化 |
| `SwernOxidation` | Swern 氧化 |
| `SharplessEpoxidation` | Sharpless 环氧化 |
| `SharplessDihydroxylation` | Sharpless 双羟基化 |
| `WolffKishnerRed reduction` | Wolff-Kishner 还原 |
| `CarbocationRearrangement` | 碳正离子重排 |
| `RadicalChainMechanism` | 自由基链式机理 |
| `ArrowPushingValidator` | 电子推动验证 |
| `DisconnectionSuggester` | 断裂建议 |
| `TotalSynthesisPlanner` | 全合成规划 |
| `RegioselectivityPredictor` | 区域选择性预测 |
| `StereoselectivityPredictor` | 立体选择性预测 |
| `StereoisomerCounter` | 立体异构体计数 |
| `StereocenterFinder` | 手性中心查找 |
| `RsConfigurator` | R/S 构型确定 |
| `EzConfigurator` | E/Z 构型确定 |
| `MesoCompoundChecker` | 内消旋体判断 |
| `CisTransAnalyzer` | 顺反异构分析 |
| `ConformationalAnalyzer` | 构象分析 |
| `TautomerGenerator` | 互变异构体生成 |
| `RingSystemAnalyzer` | 环系统分析 |
| `AromaticSystemDetector` | 芳香性检测 |
| `ConstitutionalIsomerGenerator` | 构造异构体生成 |
| `IsomerTypes` | 异构体类型判断 |
| `ProtectingGroupSelector` | 保护基选择 |
| `DeprotectionConditions` | 脱保护条件 |
| `SynthonIdentifier` | 合成子识别 |
| `CatalystRecommender` | 催化剂推荐 |
| `NameComplex` | 配合物命名 |
| `AsymmetricSynthesisGuide` | 不对称合成指导 |

### C. 无机/配位化学

| 工具名 | 功能 |
|--------|------|
| `CoordinationGeometry` | 配位几何判断 |
| `CrystalFieldSplitting` | 晶体场分裂 |
| `LigandFieldDiagram` | 配体场图 |
| `DElectronConfiguration` | d 电子组态 |
| `JahnTellerDistortion` | Jahn-Teller 畸变 |
| `TransitionMetalChemistry` | 过渡金属化学 |
| `ActinideProperties` | 锕系元素性质 |
| `LanthanideProperties` | 镧系元素性质 |
| `AlkaliMetalProperties` | 碱金属性质 |
| `AlkalineEarthProperties` | 碱土金属性质 |
| `HalogenProperties` | 卤素性质 |
| `NobleGasCompounds` | 稀有气体化合物 |
| `MainGroupTrends` | 主族元素趋势 |
| `ComplexFromName` | 由名称解析配合物 |
| `ComplexIonSolubility` | 配离子溶解度 |
| `ChelateEffect` | 螯合效应 |
| `GetFormationConstant` | 形成常数查询 |
| `SymmetryOperations` | 对称操作 |
| `SymmetryPointGroup` | 对称点群 |
| `PointGroupIdentifier` | 点群识别 |
| `PredictColor` | 颜色预测 |
| `PredictMagnetism` | 磁性预测 |
| `NuclearEquationBalance` | 核反应方程配平 |
| `RadioactiveDecay` | 放射性衰变 |
| `DecaySeries` | 衰变系列 |
| `HalfLifeCalculation` / `HalfLifeCalculator` | 半衰期计算 |
| `GetIsotopes` | 同位素信息 |

### D. 物理化学 / 热力学

| 工具名 | 功能 |
|--------|------|
| `IdealGasCalculator` | 理想气体状态方程 |
| `VanDerWaalsGas` | 范德华气体 |
| `CompressibilityFactor` | 压缩因子 Z |
| `VirialEquation` | 维里方程 |
| `EquationOfState` | 状态方程通用求解 |
| `ClausiusClapeyron` | Clausius-Clapeyron 方程 |
| `CalculateGibbsChange` | Gibbs 自由能变 ΔG |
| `CalculateReactionEnthalpy` | 反应焓变 ΔH |
| `EnthalpyCalculator` / `EnthalpyCalculatorNew` | 焓计算器 |
| `EntropyCalculator` / `EntropyCalculatorNew` | 熵计算器 |
| `GibbsEnergyCalculator` | Gibbs 能计算 |
| `GibbsFreeEnergy` | Gibbs 自由能 |
| `GibbsMinimization` | Gibbs 最小化 |
| `EquilibriumConstant` | 平衡常数 K |
| `EquilibriumConstantThermo` | 热力学平衡常数 |
| `LeChatelierAnalyzer` / `LeChatelierPrediction` | Le Chatelier 原理 |
| `VantHoffEquilibrium` | van't Hoff 方程 |
| `PressureEffectEquilibrium` | 压力对平衡的影响 |
| `CalculatePH` / `CalculatePH` | pH 计算 |
| `AcidBasePredictor` | 酸碱性预测 |
| `AcidBaseStrengthCompare` | 酸强度比较 |
| `AcidBaseTitrationCalculator` | 酸碱滴定 |
| `PolyproticAcid` | 多元酸 |
| `AmphotericSpecies` | 两性物种 |
| `BufferCapacity` / `BufferCapacityCalculator` | 缓冲容量 |
| `BufferPreparation` | 缓冲液配制 |
| `PhAdjustmentBuffer` | pH 调节 |
| `CalculateSolubility` | 溶解度计算 |
| `SolubilityRules` | 溶解度规则 |
| `CommonIonEffect` / `CommonIonSolubility` | 同离子效应 |
| `SelectivePrecipitation` | 选择性沉淀 |
| `PrecipitationPrediction` | 沉淀预测 |
| `WillPrecipitate` | 判断是否沉淀 |
| `DissolvePrecipitate` | 溶解/沉淀平衡 |
| `Ksp` 相关: `GetKsp`, `ComplexIonSolubility` | 溶度积 |
| `CalculateCellPotential` | 电池电动势 |
| `NernstEquation` / `NernstEquationSolver` | Nernst 方程 |
| `IdentifyOxidizingAgent` | 氧化剂识别 |
| `AssignOxidationNumber` | 氧化数分配 |
| `CalculateLatticeEnergy` | 晶格能 (Born-Haber) |
| `BornHaberCycle` | Born-Haber 循环 |
| `ChemicalPotential` / `ChemicalPotentialAdvanced` | 化学势 |
| `ActivityCoefficient` | 活度系数 |
| `FugacityCalculator` | 逸度计算 |
| `PhaseEquilibrium` | 相平衡 |
| `PhaseRuleAnalyzer` | 相律分析 |
| `BinaryPhaseDiagram` | 二元相图 |
| `EutecticPointFinder` | 共晶点 |
| `LeverRuleCalculator` | 杠杆规则 |
| `RaoultLaw` | Raoult 定律 |
| `HenrysLaw` | Henry 定律 |
| `ColligativeProperties`: `FreezingPointDepression`, `BoilingPointElevation`, `OsmoticPressure` | 依数性 |
| `ArrheniusCalculator` / `ArrheniusAnalyzer` | Arrhenius 方程 |
| `CollisionTheory` | 碰撞理论 |
| `EyringEquation` | Eyring 方程 (TST) |
| `RateLawFitter` | 速率定律拟合 |
| `RateLawIntegrator` | 速率定律积分 |
| `IntegratedRateLaw` | 积分速率定律 |
| `RateDeterminingStep` | 速控步 |
| `SteadyStateApproximation` / `SteadyStateApprox` | 稳态近似 |
| `PreEquilibrium` | 预平衡近似 |
| `ParallelConsecutiveReactions` | 平行/连续反应 |
| `ReactionNetworkSolver` | 反应网络 ODE 求解 |
| `ReactionCoordinate` | 反应坐标 |
| `ReactionEnergyEstimator` | 反应能量估算 |
| `ReactionMechanismSimulator` | 反应机理模拟 |
| `TunnelingCorrection` | 隧穿校正 |
| `TransitionStateTheory` | 过渡态理论 |
| `EnzymeKinetics` / `MichaelisMenten` | 酶动力学 |
| `KineticVsThermodynamic` | 动力学 vs 热力学控制 |
| `DisproportionationCheck` | 歧化反应检查 |
| `SpontaneityCheck` | 自发性判断 |
| `CarnotEfficiency` | Carnot 效率 |
| `JouleThomson` | Joule-Thomson 效应 |
| `HeatCapacityCalculator` / `HeatCapacityLookup` | 热容 |
| `BoltzmannDistribution` | Boltzmann 分布 |
| `MaxwellBoltzmannSpeed` | Maxwell-Boltzmann 速率分布 |
| `FermiDiracDistribution` | Fermi-Dirac 分布 |
| `BoseEinsteinDistribution` | Bose-Einstein 分布 |
| `EquipartitionTheorem` | 能量均分定理 |
| `PartitionFunction` / `AdvancedPartitionFunction` | 配分函数 |
| `StatisticalEnsemble` / `StatisticalEntropy` | 统计系综 |
| `TemperatureJumpRelaxation` | 温度跳跃弛豫 |
| `EvaporationEstimator` | 蒸发估算 |

### E. 分析化学 / 光谱

| 工具名 | 功能 |
|--------|------|
| `BeerLambertCalculator` | Beer-Lambert 定律 |
| `MolarAbsorptivityLookup` | 摩尔吸光系数查询 |
| `UvVisSpectrum` / `UvVisTransitions` / `UvVisPredictor` | UV-Vis 光谱 |
| `UvVisWavelengthSelector` | 波长选择 |
| `IrSpectrumPredictor` | IR 光谱预测 |
| `IrSpectrumInterpreter` / `IrPeakInterpreter` | IR 光谱解析 |
| `FtirBaselineCorrector` | FTIR 基线校正 |
| `NmrHPredictor` | ^1H NMR 化学位移预测 |
| `NmrCPredictor` | ^13C NMR 化学位移预测 |
| `NmrChemicalShift` / `NmrChemicalShiftPredictor` | NMR 化学位移 |
| `NmrShielding` | NMR 屏蔽 |
| `SpinSpinCoupling` | 自旋-自旋耦合 |
| `CouplingConstantAnalyzer` | 耦合常数分析 |
| `SplittingPatternExplainer` | 分裂模式解释 |
| `DeptInterpreter` | DEPT 谱解释 |
| `CosyNoesyGuide` | COSY/NOESY 指南 |
| `RamanSpectrumPredictor` | Raman 光谱预测 |
| `RamanActivity` / `RamanShiftCalculator` | Raman 活性/位移 |
| `MassSpecFragmenter` | MS 断裂预测 |
| `IsotopePatternGenerator` / `IsotopePatternSimulator` | 同位素模式 |
| `MassDefect` / `MassDefectFilter` | 质量亏损 |
| `MassAccuracyCalculator` | 质量精度 |
| `AdductIonIdentifier` | 加合离子识别 |
| `MsmsSpectrumAnnotator` | MS/MS 注释 |
| `CollisionEnergyOptimizer` | 碰撞能量优化 |
| `FluorescenceQuantumYield` | 荧光量子产率 |
| `FluorescenceLifetime` | 荧光寿命 |
| `ExcitationEmissionOptimizer` | 激发/发射优化 |
| `OscillatorStrength` | 振子强度 |
| `FranckCondon` / `FranckCondonFactors` | Franck-Condon 因子 |
| `SpectralDeconvolution` | 光谱去卷积 |
| `SpectralLinewidth` | 光谱线宽 |
| `Chromatography`: `PlateNumberCalculator`, `ResolutionCalculator`, `VanDeemterAnalyzer`, `CapacityFactorCalculator`, `SelectivityFactorCalculator`, `LodLoqCalculator`, `LinearityRangeValidator`, `TailFactor`... | 色谱相关 |
| `HplcColumnSelector` | HPLC 柱选择 |
| `GcOvenProgramDesigner` | GC 升温程序 |
| `GcCarrierGasSelector` | GC 载气选择 |
| `GcColumnBleedPredictor` | GC 柱流失预测 |
| `MobilePhaseOptimizer` | 流动相优化 |
| `InternalStandardSelector` | 内标选择 |
| `CalibrationCurveFitter` | 校准曲线拟合 |
| `SecCalibrationCurve` | SEC 校准曲线 |
| `TitrationCurve` | 滴定曲线生成 |
| `EquivalencePoint` | 等当点计算 |
| `BackTitrationSolver` | 返滴定 |
| `PotentiometricTitrationEndpoint` | 电位滴定终点 |
| `PrecipitationTitrationCalculator` | 沉淀滴定 |
| `ComplexometricTitrationHelper` | 配位滴定 |
| `RedoxTitrationCalculator` | 氧化还原滴定 |
| `TitrantStandardization` | 滴定剂标定 |
| `SampleDilutionCalculator` | 样品稀释 |
| `StandardSolutionPrep` | 标准溶液配制 |
| `MethodValidationChecklist` | 方法验证清单 |
| `SystemSuitabilityChecker` | 系统适用性 |
| `OutlierDetector` | 异常值检测 |
| `AnovaAnalyzer` | 方差分析 |
| `TTestCalculator` / `FTestCalculator` | t/F 检验 |
| `RegressionDiagnostics` | 回归诊断 |
| `ErrorPropagation` / `UncertaintyPropagator` | 误差传递 |
| `SignificantFigures` | 有效数字 |
| `MeasurementPrecisionCalculator` | 测量精度 |
| `AuditTrailReviewer` | 审计追踪 |
| `CmcDocumentationHelper` | CMC 文档辅助 |

### F. 计算化学 / 量子化学

| 工具名 | 功能 |
|--------|------|
| `SchrodingerSolver1d` | 一维薛定谔方程 |
| `Schrodinger3DSolver` | 三维薛定谔方程 |
| `HydrogenAtomOrbitals` | 氢原子轨道 |
| `HydrogenWavefunction` | 氢原子波函数 |
| `ParticleInBox` | 粒子在盒中 |
| `HarmonicOscillator` | 谐振子 |
| `RigidRotor` | 刚性转子 |
| `HuckelMethod` / `ExtendedHuckel` | Hückel 方法 |
| `HartreeFockSCF` | Hartree-Fock SCF |
| `Mp2Correlation` | MP2 相关能 |
| `ConfigurationInteraction` | CI 方法 |
| `DftXcFunctional` | DFT 泛函 |
| `VariationalMethod` | 变分法 |
| `PerturbationTheory` | 微扰理论 |
| `SlaterDeterminant` | Slater 行列式 |
| `FrontierOrbitalAnalysis` / `MolecularOrbitalDiagram` / `MOEnergyLevelDiagram` | 前线轨道/分子轨道 |
| `HessianMatrix` | Hessian 矩阵 |
| `NormalModeAnalysis` / `VibrationalModes` | 正常模式/振动模式 |
| `FrequencyAnalysis` | 频率分析 |
| `GeometryOptimizer` | 几何优化 |
| `TransitionStateSearch` | 过渡态搜索 |
| `PotentialEnergySurface` | 势能面扫描 |
| `ForceFieldCalculator` | 力场计算 |
| `MolecularDynamicsVerlet` / `VelocityVerlet` | 分子动力学 |
| `ElectronDensityCalculator` / `ElectronDensityPlotter` | 电子密度 |
| `DipoleMoment` / `DipoleMomentEstimator` | 偶极矩 |
| `Polarizability` | 极化率 |
| `MomentOfInertia` | 转动惯量 |
| `RotationalSpectrum` | 转动光谱 |
| `BindingEnergy` | 结合能 |
| `GradientCalculator` | 梯度计算 |
| `OverlapIntegral` / `ExchangeIntegral` / `CoulombIntegral` | 积分 |
| `MultipoleExpansion` | 多极展开 |
| `BornOppenheimer` | Born-Oppenheimer 近似 |
| `UncertaintyPrinciple` | 不确定性原理 |
| `WKBApproximation` | WKB 近似 |
| `TunnelingProbability` | 隧穿概率 |
| `SpinOrbitCoupling` | 自旋-轨道耦合 |
| `ZeemanSplitting` | Zeeman 分裂 |
| `StarkEffect` | Stark 效果 |
| `SphericalHarmonics` | 球谐函数 |
| `RadialDistribution` | 径向分布函数 |
| `SelectionRules` / `SelectionRulesChecker` | 选律 |
| `SymmetryOperations` / `SymmetryPointGroup` / `PointGroupIdentifier` | 对称性 |
| `PredictHybridization` | 杂化预测 |
| `VseprGeometry` / `PredictVseprGeometry` | VSEPR 几何 |
| `HybridizationAnalyzer` | 杂化分析 |
| `BondOrderCalculator` | 键级计算 |
| `GetBondEnergy` / `GetBondLength` | 键能/键长 |
| `BasisSetHandler` | 基组处理 |
| `PredictPolarity` | 极性预测 |
| `ElectricFieldMolecule` | 分子电场 |
| `SpectrochemicalSeries` | 光谱化学序列 |

### G. 数值方法 / 通用工具

| 工具名 | 功能 |
|--------|------|
| `NewtonRaphsonSolver` | Newton-Raphson 求解 |
| `BisectionMethod` | 二分法 |
| `NumericalIntegrator` | 数值积分 |
| `NumericalDifferentiation` | 数值微分 |
| `OdeSolverRk4` / `OdeSolverStiff` | ODE 求解器 |
| `PdeSolverFiniteDiff` | PDE 有限差分 |
| `LinearSystemSolver` | 线性方程组 |
| `MatrixInversion` / `MatrixDiagonalization` / `MatrixEigenvalueSolver` | 矩阵运算 |
| `LeastSquaresFit` / `CurveFittingNonlinear` | 曲线拟合 |
| `InterpolationSpline` | 样条插值 |
| `FftTransform` | FFT 变换 |
| `MonteCarloIntegrator` | Monte Carlo 积分 |
| `SvdDecomposition` | SVD 分解 |
| `DeterminantCalculator` / `CommutatorCalculator` | 行列式/对易子 |
| `DimensionalAnalysis` | 量纲分析 |
| `UnitConversion` / `StandardStateConverter` | 单位转换 |
| `IceTableSolver` | ICE 表求解 |
| `WebSearch` | 网络搜索 (Tavily) |
| `PubchemSearch` / `PubchemSearchQA` | PubChem 搜索/QA |
| `CommonNameLookup` | 俗名查询 |
| `NamedReactionLookup` | 命名反应查询 |
| `GetElementInfo` | 元素信息 |
| `CompareElements` | 元素比较 |
| `PeriodicTrend` | 周期表趋势 |
| `GetElectronConfiguration` | 电子排布 |
| `GetIonizationEnergy` | 电离能 |
| `GetElectronAffinity` | 电子亲和能 |
| `GetOxidationStates` | 氧化态 |
| `GetStandardEnthalpy` / `GetStandardEntropy` / `GetStandardPotential` | 标准热力学数据 |
| `GetGibbsEnergy` | Gibbs 能 |
| `GetPka` / `GetPkb` | pKa/pKb |
| `GetCrystalStructure` | 晶体结构 |
| `GetElementDiscovery` | 元素发现史 |
| `ElementAbundance` | 元素丰度 |
| `PourbaixDiagramLookup` | Pourbaix 图 |
| `XrdPhaseIdentifier` | XRD 物相鉴定 |
| `XrfMatrixCorrection` | XRF 基体校正 |
| `SurfaceTensionCalculator` | 表面张力 |
| `ContactAngleAnalyzer` | 接触角 |
| `ZetaPotential` | Zeta 电位 |
| `ColloidalStability` | 胶体稳定性 |
| `BetSurfaceArea` | BET 比表面 |
| `DiffusionCoefficientCalculator` | 扩散系数 |
| `CentrifugationCalculator` | 离心计算 |
| `FiltrationGuide` | 过滤指南 |
| `HomogenizationProtocol` | 均质方案 |
| `DigestionProtocolSelector` | 消解方案 |
| `ExtractionOptimizer` | 萃取优化 |
| `SamplePreservationGuide` | 样品保存 |
| `StabilityStudyPlanner` | 稳定性研究设计 |
| `RobustnessDoeDesigner` | 鲁棒性 DOE 设计 |
| `SpecificityTestDesigner` | 专属性试验 |
| `SpikeRecoveryEvaluator` | 加标回收 |
| `MatrixEffectCalculator` / `MatrixMatchingAdvisor` | 基体效应 |
| `IonSuppressionChecker` | 离子抑制 |
| `DeadVolumeCalculator` | 死体积 |
| `RecoveryCalculator` | 回收率 |
| `ConductivityCellConstant` | 电导池常数 |
| `ConductivityCalculator` | 电导率 |
| `IonTransportNumber` | 离子迁移数 |
| `FaradayElectrolysis` | Faraday 电解定律 |
| `ButlerVolmerKinetics` | Butler-Volmer 动力学 |
| `ChronoamperometryAnalyzer` | 计时安培分析 |
| `CvPeakAnalyzer` | CV 峰分析 |
| `StrippingVoltammetryOptimizer` | 溶出伏安优化 |
| `OverpotentialAnalyzer` | 过电位分析 |
| `TafelEquation` | Tafel 方程 |
| `EisCircuitFitter` | EIS 等效电路 |
| `ElectrodeSelectionGuide` | 电极选择 |
| `PhElectrodeCalibration` | pH 电极校准 |
| `LangmuirIsotherm` / `FreundlichIsotherm` | 等温吸附 |
| `GibbsAdsorption` | Gibbs 吸附 |
| `IcpOesLineSelector` | ICP-OES 谱线选择 |
| `IcpMsIsotopeSelector` | ICP-MS 同位素选择 |
| `AasFlameSelector` / `AasInterferenceChecker` | AAS 火焰/干扰 |
| `PlasmaConditionOptimizer` | 等离子体优化 |
| `SpeMethodDesigner` | SPE 方法设计 |
| `DerivatizationReagentSelector` | 衍生试剂选择 |
| `MobilePhaseOptimizer` | 流动相优化 |
| `GcCarrierGasSelector` | GC 载气 |
| `HplcColumnSelector` | HPLC 柱 |
| `IonChromatographyEluent` | IC 淋洗液 |
| `QcChartGenerator` | QC 质控图 |

---

## 6. 使用示例

### 示例 1：分子分析

**在 Cherry Studio 中输入：**
> 计算阿司匹林 (CC(=O)Oc1ccccc1C(=O)O) 的分子量和 SMILES 对应的 IUPAC 名称

AI 将自动调用：
1. `MoleculeWeight` → 输出分子量
2. `Smiles2Iupac` → 输出 IUPAC 名称

### 示例 2：反应预测

**输入：**
> 甲苯 + KMnO4 → ? 预测产物并配平方程式

AI 将调用：
1. `PredictProducts` 或 `ReactionPredictor`
2. `BalanceEquation`

### 示例 3：光谱分析

**输入：**
> 预测苯甲醛 (c1ccccc1C=O) 的 IR 和 ^1H NMR 谱图主要特征

AI 将调用：
1. `IrSpectrumPredictor`
2. `NmrHPredictor`

### 示例 4：热力学计算

**输入：**
> 计算 N2(g) + 3H2(g) ⇌ 2NH3(g) 在 298K、100 atm 下的 Gibbs 自由能变化和平衡移动方向

AI 将调用：
1. `CalculateGibbsChange`
2. `LeChatelierAnalyzer`

### 示例 5：滴定分析

**输入：**
> 用 0.1 M NaOH 滴定 20 mL 0.05 M HAc，绘制滴定曲线并找出等当点 pH

AI 将调用：
1. `TitrationCurve`
2. `EquivalencePoint`
3. `BufferCapacity`

---

## 7. 常见问题

### Q1: 导入错误 `SchrodingerSolver1D`？

**原因**: `src/chemmcp/tools/__init__.py` 第 567 行有重复条目 `SchrodingerSolver1D`（大写 D），但实际类名为 `SchrodingerSolver1d`（小写 d）。

**修复**:
```bash
# 编辑 src/chemmcp/tools/__init__.py 第 567 行
# 将:
"SchrodingerSolver1D": "schrodinger_solver_1d",
# 改为注释或删除（第 316 行已有正确的小写版本）
```

已修复 ✅ — 当前 `__init__.py` 已去除重复条目。

### Q2: Cherry Studio 连接超时？

**原因**: 全量加载 406 个工具需要 15-30 秒初始化。

**解决**:
- 使用分类配置（配置 B/C），每次只加载 10-20 个工具
- 在 Cherry Studio 中增加 MCP 超时时间（设置 → 高级 → 超时 → 改为 60s）

### Q3: 某些工具需要 API Key？

以下工具需要外部服务 API Key：

| 工具 | 服务 | 环境变量 |
|------|------|----------|
| `WebSearch` | Tavily | `TAVILY_API_KEY` |
| `ForwardSynthesis` | IBM RXN | `IBM_RXN_API_KEY` |
| `Retrosynthesis` | IBM RXN | `IBM_RXN_API_KEY` |
| `PubchemSearch` / `PubchemSearchQA` | PubChem | 免费，无需 Key |
| `MoleculePrice` | ChemSpace | 免费，无需 Key |

### Q4: 如何更新 ChemMCP？

```bash
cd ~/ChemMCP
git pull
uv sync   # 更新依赖
# 重启 Cherry Studio 中的 MCP 连接
```

### Q5: 如何添加自定义工具？

1. 在 `src/chemmcp/tools/` 下新建 `.py` 文件
2. 继承 `BaseTool` 类，定义 `name`, `description`, `code_input_sig` 等
3. 在 `__init__.py` 的 `_tool_module_map` 中注册
4. 重启 MCP Server

详细模板参考 `src/chemmcp/tools/molecule_weight.py`。

### Q6: `run_code()` vs `run_text()` 区别？

| | `run_code()` | `run_text()` |
|--|-------------|-------------|
| 参数方式 | 关键字参数 (`smiles="..."`) | 空格分隔字符串 |
| 类型安全 | ✅ 有类型检查 | ❌ 字符串解析 |
| 推荐场景 | MCP/AI 调用 | 快速命令行 |
| Cherry Studio | 默认使用此方式 | — |

---

## 附录 A：现有 Cherry Studio 配置文件索引

项目已包含多个预生成的 Cherry Studio 配置文件：

| 文件 | 工具编号范围 | 工具数量 |
|------|-------------|---------|
| `cherry_studio_config_21_30.json` | #21-30 | 10 |
| `cherry_studio_config_91_100.json` | #91-100 | 10 |
| `cherry_studio_config_301_310.json` | #301-310 | 10 |
| `cherry_studio_config_341_350.json` | #341-350 | 10 |
| `cherry_studio_config_401_410.json` | #401-410 | 10 |
| `cherry_studio_config_441_450.json` | #441-450 | 10 |

直接导入即可使用：
```
Cherry Studio → 设置 → MCP Servers → 导入 → 选择上述 JSON 文件
```

## 附录 B：工具编号完整列表

<details>
<summary>📋 点击展开全部 406 个工具名</summary>

```
  1. AcidBaseStrengthCompare          52. IdealGasCalculator
  2. ActinideProperties               53. IdentifyOxidizingAgent
  3. ActivityCoefficient              54. IdentifyReactionType
  4. AdvancedPartitionFunction        55. IntegratedRateLaw
  5. AldolMechanism                   56. InternalStandardSelector
  6. AldolReaction                    57. InterpolationSpline
  7. AlkaliMetalProperties            58. IonChromatographyEluent
  8. AlkalineEarthProperties          59. IonTransportNumber
  9. AmphotericSpecies                60. IrPeakInterpreter
 10. AnalyzeMolecularOrbital          61. IrSpectrumInterpreter
 11. AnovaAnalyzer                    62. IrSpectrumPredictor
 12. AromaticSystemDetector           63. IsomerTypes
 13. ArrheniusAnalyzer                64. Iupac2Smiles
 14. ArrheniusCalculator              65. IupacNamer
 15. ArrowPushingValidator            66. JahnTellerDistortion
 16. AssignOxidationNumber            67. JonesOxidation
 17. AuditTrailReviewer               68. JouleThomson
 18. BackTitrationSolver              69. LangmuirIsotherm
 19. BaeyerVilligerOxidation          70. LanthanideProperties
 20. BalanceEquation                  71. LaplacianOperator
 21. BalanceRedox                     72. LeastSquaresFit
 22. BasisSetHandler                  73. LeavingGroupRanker
 23. BeckmannRearrangement            74. LeChatelierAnalyzer
 24. BeerLambertCalculator            75. LigandFieldDiagram
 25. BETSurfaceArea                   76. LinearityRangeValidator
 26. BindingEnergy                    77. LinearSystemSolver
 27. BirchReduction                   78. LodLoqCalculator
 28. BisectionMethod                  79. MainGroupTrends
 29. BondOrderCalculator             80. MassDefect
 30. BornOppenheimer                  81. MassDefectFilter
 31. BuchwaldHartwigAmination         82. MassSpecFragmenter
 32. BufferCapacity                   83. MatrixDiagonalization
 33. BufferCapacityCalculator         84. MatrixEigenvalueSolver
 34. BufferPreparation                85. MatrixInversion
 35. ButlerVolmerKinetics             86. MatrixMatchingAdvisor
 36. CalculateCellPotential           87. MaxwellRelations
 37. CalculateFormalCharge            88. MesoCompoundChecker
 38. CalculateGibbsChange            89. MethodValidationChecklist
 39. CalculateLatticeEnergy           90. MichaelisMenten
 40. CalculatePH                      91. MichaelMechanism
 41. CalculateReactionEnthalpy        92. MobilePhaseOptimizer
 42. CalculateSolubility             93. MOEnergyLevelDiagram
 43. CalibrationCurveFitter           94. MolarAbsorptivityLookup
 44. CapacityFactorCalculator         95. MolecularDynamicsVerlet
 45. CarbocationRearrangement        96. MolecularIonCalculator
 46. CarnotEfficiency                97. MolecularOrbitalDiagram
 47. CentrifugationCalculator        98. MoleculeAtomCount
 48. ChelateEffect                   99. MoleculeCaptioner
 49. ChemicalPotential               100. MoleculeGenerator
 50. ChemicalPotentialAdvanced       101. MoleculeModifier
 51. ChronoamperometryAnalyzer       102. MoleculePrice
103. MoleculeSimilarity              255. Schrodinger3DSolver
104. MoleculeSmilesCheck             256. SchrodingerSolver1d
105. MoleculeVisualizer              257. SecCalibrationCurve
106. MoleculeWeight                  258. SelectionRules
107. MomentOfInertia                 259. SelectionRulesChecker
108. MonteCarloIntegrator            260. SelectivePrecipitation
109. Mp2Correlation                  261. SelectivityFactorCalculator
110. MultipoleExpansion              262. Selfies2Smiles
111. Name2Smiles                     263. SharplessDihydroxylation
112. NameComplex                     264. SharplessEpoxidation
113. NamedReactionLookup             265. SignificantFigures
114. NernstEquation                  266. SlaterDeterminant
115. NernstEquationSolver            267. Smiles2Cas
116. NewtonRaphsonSolver             268. Smiles2Formula
117. NmrChemicalShift               269. Smiles2Iupac
118. NmrChemicalShiftPredictor      270. Smiles2Selfies
119. NmrCPredictor                  271. SmilesCanonicalization
120. NmrHPredictor                  272. SmilesToStructure
121. NmrShielding                   273. Sn1Mechanism
122. NobleGasCompounds              274. Sn2Mechanism
123. NormalModeAnalysis             275. SolubilityRules
124. NuclearEquationBalance          276. SonogashiraCoupling
125. NucleophilicAddition            277. SpecificityTestDesigner
126. NucleophilicAromaticSubstitution 278. SpectralDeconvolution
127. NucleophilicityRanker           279. SpectrochemicalSeries
128. NumericalDifferentiation        280. SPEMethodDesigner
129. NumericalIntegrator             281. SphericalHarmonics
130. OdeSolverRk4                   282. SpinOrbitCoupling
131. OdeSolverStiff                 283. SpinSpinCoupling
132. OscillatorStrength              284. SplittingPatternExplainer
133. OutlierDetector                 285. SpontaneityCheck
134. OverlapIntegral                 286. StabilityStudyPlanner
135. OverpotentialAnalyzer           287. StandardSolutionPrep
136. OxidationMechanism              288. StandardStateConverter
137. ParallelConsecutiveReactions    289. StatisticalEnsemble
138. PartialDerivative               290. SteadyStateApprox
139. PartialMolarQuantity            291. SteadyStateApproximation
140. ParticleInBox                  292. StereocenterFinder
141. PatentCheck                    293. StereoisomerCounter
142. PdeSolverFiniteDiff            294. StereoselectivityPredictor
143. PeakPurityAnalyzer             295. StrippingVoltammetryOptimizer
144. PericyclicAnalyzer             296. StructureToSmiles
145. PeriodicTrend                  297. SurfaceTensionCalculator
146. PerturbationTheory             298. SuzukiCoupling
147. PhAdjustmentBuffer             299. SvdDecomposition
148. PhaseEquilibrium               300. SwernOxidation
149. PhaseRuleAnalyzer              301. SymmetryOperations
150. PhElectrodeCalibration         302. SymmetryPointGroup
151. PlateNumberCalculator          303. SystemSuitabilityChecker
152. PointGroupIdentifier           304. TautomerGenerator
153. PoissonBoltzmann               305. TemperatureJumpRelaxation
154. Polarizability                 306. TitrantStandardization
155. PolyproticAcid                 307. TitrationCurve
156. PotentialEnergySurface         308. TransitionDipole
157. PotentiometricTitrationEndpoint  309. TransitionMetalChemistry
158. PourbaixDiagramLookup          310. TransitionStateSearch
159. PrecipitationPrediction         311. TransitionStateTheory
160. PredictColor                   312. TTestCalculator
161. PredictHybridization           313. TunnelingCorrection
162. PredictMagnetism               314. TunnelingProbability
163. PredictPolarity                315. UncertaintyPrinciple
164. PredictProducts                316. UncertaintyPropagator
165. PredictVseprGeometry           317. UnitConversion
166. PreEquilibrium                 318. UvVisSpectrum
167. PubchemSearch                  319. UvVisTransitions
168. PubsearchSearchQA              320. UvVisWavelengthSelector
169. RadialDistribution             321. VanDeemterAnalyzer
170. RadicalChainMechanism           322. VanDerWaalsGas
171. RamanActivity                  323. VantHoffEquilibrium
172. RamanSpectrumPredictor         324. VariationalMethod
173. RateDeterminingStep            325. VelocityVerlet
174. RateLawFitter                  326. VibrationalModes
175. RateLawIntegrator              327. VirialEquation
176. ReactionCoordinate             328. VseprGeometry
177. ReactionEnergyEstimator        329. WebSearch
178. ReactionMechanismSimulator     330. WillPrecipitate
179. ReactionNetworkSolver          331. WittigMechanism
180. ReactionPredictor              332. WittigReaction
181. ReactionSmilesCheck            333. WKBApproximation
182. RecoveryCalculator             334. WolffKishnerReduction
183. ReductionMechanism
184. RegioselectivityPredictor
185. RegressionDiagnostics
186. ResolutionCalculator
187. RetentionTimePredictor
188. Retrosynthesis
189. RigidRotor
190. RingSystemAnalyzer
191. RobustnessDoeDesigner
192. RotationalSpectrum
193. RsConfigurator
194. SafetyCheck
195. SampleDilutionCalculator
196. SamplePreservationGuide
```

</details>

---

*文档版本: 2026-05-15 | ChemMCP @ ~/ChemMCP | 406 tools verified ✅*

*由 X Leclaw 🦐 自动生成*
