# ChemMCP

## 项目简介

ChemMCP 是一个**全面、易用且可扩展的化学工具包，专为大型语言模型和AI助手设计**，兼容[模型上下文协议（MCP）](https://modelcontextprotocol.org/)。通过集成**涵盖有机化学、无机化学、物理化学和分析化学的强大工具**，ChemMCP**赋能通用AI模型具备专业化学能力**，使其能够执行分子分析、性质预测、反应合成、热力学计算和光谱模拟等任务，而无需进行特定领域的训练。ChemMCP 可以无缝集成到您的研究工作流中，用于数据处理、智能体应用和计算实验。

[MCP（模型上下文协议）](https://modelcontextprotocol.io/introduction) 是一个允许AI模型通过标准化接口访问外部工具和资源的框架。ChemMCP 利用此架构，弥合了通用AI模型与专业化学工具之间的鸿沟，实现了化学专业知识与AI工作流的无缝集成。

具体而言，ChemMCP 提供以下关键特性：

- **🔌 为AI助手提供即插即用的化学工具**：ChemMCP工具可在几分钟内集成到任何[支持MCP的LLM客户端](https://github.com/punkpeye/awesome-mcp-clients)中，让有机、无机、物化、分析化学领域的研究者无需额外训练即可为LLM增强专业能力。
- **🛠️ 支持自定义工作流的独立工具包**：凭借其解耦设计和统一接口，ChemMCP工具可以轻松导入到您的科学工作流中，用于处理数据、组装流水线步骤、进行计算实验或构建定制的智能体应用——通过MCP或Python皆可，任您选择。
- **📦 模块化与可扩展设计**：添加新工具就像编写一个Python文件一样简单。所有工具都遵循一致的架构，确保接口清晰、易于维护。
- **🧪 跨学科覆盖**：旨在满足有机合成、无机化合物分析、物理化学计算和分析光谱学的需求。

我们将持续在ChemMCP中添加和维护工具。**非常欢迎您的贡献，例如反馈意见、维护现有工具或添加新工具！**

## 工具列表

根据功能和化学子学科，工具分类如下：

- **通用工具**：提供跨所有化学领域的广泛信息检索和网络搜索。
- **分子工具**：提供与化学化合物及其性质相关的各种分析、预测和转换（适用于有机和无机分子）。
- **反应工具**：预测化学反应产物，并为合成给定产物建议潜在反应物（有机和无机合成）。
- **物理化学工具**：执行热力学计算、动力学模拟和量子化学性质预测。
- **分析化学工具**：模拟光谱数据（核磁共振、红外、质谱）并辅助色谱分析。

## 伦理与负责任使用声明

ChemMCP 是一个开源工具包，它将语言模型和智能体与化学工具及公开可用的化学数据相结合，以支持跨所有化学子学科的AI for Science研究。虽然ChemMCP提供了强大的功能，但必须承认其使用可能带来的潜在风险。

1.  **安全与责任**

    ChemMCP包含旨在帮助识别危险分子和反应的安全检查工具；但是，工具包本身并不强制使用这些工具。由于ChemMCP不是独立的智能体，而是一个开源资源，我们无法保证每个用户都会采用安全措施。

    与ChemMCP配对的大型语言模型通常包含其自身的安全机制，通常会拒绝涉及非法或不道德应用的请求。尽管如此，用户仍须全权负责确保所有活动符合适用的安全规程、机构规定以及使用地所有司法管辖区的法律要求。

    由于ChemMCP仅访问公开可用的工具和数据，我们对任何危险或非法用途不承担责任。用户必须验证其工作流程是安全、合法且符合伦理的，特别是在处理反应性化合物、含能材料或受控物质时。

2.  **预期用途**

    ChemMCP 仅提供给学术、工业和政府机构用于合法研究、教育和调查目的。

    严禁将ChemMCP用于设计、制造或推荐有害物质（例如化学毒素、武器或非法药物）。任何试图利用该工具包进行恶意活动的行为都违反了我们的条款和伦理准则。

3.  **限制与免责声明**

    ChemMCP不保证其输出的准确性、完整性或安全性。所有计算、预测和推断均“按原样”提供，不作任何形式的保证。用户应运用专业判断，并在适当时进行实验验证，特别是对于：
    *   可能涉及危险中间体的合成路线
    *   用于安全关键应用的物理性质预测
    *   为监管决策提供信息的数据解析

    对于因使用ChemMCP而产生的任何直接或间接后果（财务、法律或其他方面），我们概不负责。

4.  **贡献与安全保障**

    我们鼓励社区贡献，以增强ChemMCP的安全功能并促进负责任的使用。贡献者在引入新工具或数据源时，必须记录任何已识别的风险、潜在故障模式和缓解策略，特别是涉及以下内容的工具：
    *   高能化合物或反应中间体
    *   毒理学或环境影响预测
    *   受控物质数据库

    在合并新功能之前，维护者应审查拟议的更改是否存在可能的滥用场景，并相应更新文档。

5.  **用户协议**

    通过安装或调用ChemMCP，您同意：

    1.  在使用工具包时遵守所有适用的法律、法规、机构政策和职业道德准则。
    2.  在处理潜在危险数据、化合物或反应条件时，应用合理的安全检查（包括人工和自动化检查）。
    3.  避免任何可能助长有害物质或非法材料的创造、分发或使用的活动。
    4.  承认在关键应用中，计算预测不能替代实验验证。
    5.  对您使用ChemMCP及其部署产生的任何结果承担全部责任。

## 许可证

ChemMCP 基于 [Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/) 分发。该工具包的实现严重依赖于开源项目——最引人注目的是 [RDKit](https://github.com/rdkit/rdkit)（BSD 3-Clause）和 [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk)（MIT）。其他主要开源依赖项（其代码被全部或部分使用）以及任何必需的托管服务或软件均在下表中列出。使用ChemMCP及其工具，即表示您同意遵守所有引用的许可证和服务条款。

| **工具名称** | **主要开源依赖** | **托管服务/软件** | **化学领域** |
| :--- | :--- | :--- | :--- |
| BbbpPredictor | [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) (MIT), [Uni-Core](https://github.com/dptech-corp/Uni-Core) (MIT) | - | 药物/有机化学 |
| ForwardSynthesis | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT), [rxn4chemistry](https://github.com/rxn4chemistry/rxn4chemistry) (MIT) | [IBM RXN for Chemistry](https://rxn.app.accelerate.science/) | 有机合成 |
| FunctionalGroups | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT) | - | 有机化学 |
| HivInhibitorPredictor | [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) (MIT), [Uni-Core](https://github.com/dptech-corp/Uni-Core) (MIT) | - | 药物/有机化学 |
| Iupac2Smiles | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT), [molbloom](https://github.com/whitead/molbloom) (MIT), [PubChemPy](https://github.com/mcs07/PubChemPy) (MIT) | [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [ChemSpace](https://chem-space.com/) | 通用化学 |
| LogDPredictor | [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) (MIT), [Uni-Core](https://github.com/dptech-corp/Uni-Core) (MIT) | - | 物理/药物化学 |
| MoleculeAtomCount | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT) | - | 通用化学 |
| MoleculeCaptioner | [MolT5](https://github.com/blender-nlp/MolT5) (BSD 3-Clause) | - | 通用化学 |
| MoleculeGenerator | [MolT5](https://github.com/blender-nlp/MolT5) (BSD 3-Clause) | - | 通用化学 |
| MoleculeModifier | [synspace](https://github.com/whitead/synspace) (MIT) | - | 有机/合成化学 |
| MoleculePrice | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT) | [ChemSpace](https://chem-space.com/) | 通用化学 |
| MoleculeSimilarity | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT) | - | 通用化学 |
| MoleculeSmilesCheck | - | - | 通用化学 |
| MoleculeVisualizer | - | - | 通用化学 |
| MoleculeWeight | - | - | 通用化学 |
| Name2Smiles | [PubChemPy](https://github.com/mcs07/PubChemPy) (MIT) | - | 通用化学 |
| PatentCheck | [molbloom](https://github.com/whitead/molbloom) (MIT) | - | 通用化学 |
| PubchemSearch | - | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | 通用化学 |
| PubchemSearchQA | - | [PubChem](https://pubchem.ncbi.nlm.nih.gov/), 自定义LLMs | 通用化学 |
| ReactionSmilesCheck | - | - | 通用化学 |
| Retrosynthesis | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT), [rxn4chemistry](https://github.com/rxn4chemistry/rxn4chemistry) (MIT) | [IBM RXN for Chemistry](https://rxn.app.accelerate.science/) | 有机合成 |
| SafetyCheck | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT) | [PubChem](https://pubchem.ncbi.nlm.nih.gov/), 自定义LLMs | 通用化学 |
| Selfies2Smiles | [selfies](https://github.com/aspuru-guzik-group/selfies) (Apache License 2.0) | - | 计算化学 |
| SideEffectPredictor | [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) (MIT), [Uni-Core](https://github.com/dptech-corp/Uni-Core) (MIT) | - | 药物/有机化学 |
| Smiles2Cas | [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) (MIT) | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | 通用化学 |
| Smiles2Formula | - | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | 通用化学 |
| Smiles2Iupac | [PubChemPy](https://github.com/mcs07/PubChemPy) (MIT) | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | 通用化学 |
| Smiles2Selfies | [selfies](https://github.com/aspuru-guzik-group/selfies) (Apache License 2.0) | - | 计算化学 |
| SmilesCanonicalization | [LlaSMol](https://github.com/OSU-NLP-Group/LLM4Chem) (MIT) | - | 计算化学 |
| SolubilityPredictor | [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) (MIT), [Uni-Core](https://github.com/dptech-corp/Uni-Core) (MIT) | - | 物理化学 |
| ToxicityPredictor | [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) (MIT), [Uni-Core](https://github.com/dptech-corp/Uni-Core) (MIT) | - | 药物/环境化学 |
| WebSearch | [tavily-python](https://github.com/tavily-ai/tavily-python) (MIT) | [Tavily](https://www.tavily.com/) | 通用 |

**说明**：表中标记为 `-` 的单元格表示该工具由我们原创创建，除了RDKit和MCP之外不直接依赖其他开源软件，或不使用外部托管服务和软件。

**未来工具开发计划**：我们计划扩展ChemMCP，增加专门用于以下领域的工具：
*   **无机化学**：配合物分析、晶体场理论计算、对称操作
*   **物理化学**：热力学计算器、动力学模拟器、量子化学性质预测器
*   **分析化学**：核磁/红外/质谱预测器、色谱保留时间估算器

**免责声明**：

*   上表中未明确指明的任何开源软件依赖项（包括由诸如RDKit、PyTorch或科学计算库等软件包引入的任何间接或传递依赖项）仍受其自身许可证条款的约束。对于所有依赖项及其相应的许可证义务的完整清单，用户应查阅项目的requirements.txt文件或使用许可证合规性工具（例如pip-licenses）。
*   本文引用的托管服务和应用程序编程接口（API）——例如IBM RXN for Chemistry、PubChem或外部托管的语言模型——均受其自身的服务条款、可接受使用政策或等效合同协议的约束。用户有责任审查并遵守这些外部服务提供商规定的所有适用条款和条件。
*   我们的软件基于开源代码和数据构建，我们尊重其创造者的所有权和知识产权。此外，许多工具也基于某些托管服务和软件，我们相信它们的使用条款与我们的研究目的兼容。在上表中，我们已尽最大努力列出了它们的代码库/网站并提供其许可证。如有需要，我们欢迎原作者或开发者提出修改或删除相关工具的请求。