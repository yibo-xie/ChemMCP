import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ConformationalAnalyzer(BaseTool):
    """
    构象稳定性分析工具 - 分析分子的构象稳定性，生成Newman投影式、椅式构象、能量图和扭转张力分析。
    支持烷烃、环己烷衍生物、取代环己烷、丁烷等常见体系。
    """
    __version__ = "0.1.0"
    name             = "ConformationalAnalyzer"
    func_name        = "analyze_conformation"
    description      = "Analyze molecular conformational stability with Newman projections (staggered/eclipsed), chair conformations, energy diagrams, and torsional strain analysis for acyclic and cyclic molecules."
    implementation_description = "Knowledge-based conformational analysis system covering Newman projections (ethane, butane, substituted systems), cyclohexane chair/boat/twist-boat conformations, ring-flip analysis, and 1,3-diaxial interaction calculations."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Conformation", "Newman Projection", "Chair Conformation", "Cyclohexane", "Torsional Strain", "Stereochemistry"]
    required_envs    = []

    code_input_sig   = [
        ("molecule", "str", "N/A", "Molecule to analyze: 'ethane', 'butane', 'neopentane', 'cyclohexane', 'methylcyclohexane', 'trans-1,2-dimethylcyclohexane', or general SMILES/name."),
        ("analysis_type", "str", "newman", "Analysis type: 'newman' (default), 'chair', 'general', or 'all'."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'molecule [analysis_type]'. Example: 'butane newman' or 'cyclohexane chair'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed conformational analysis with ASCII projections, energy values, strain breakdown, and stability ranking."),
    ]

    examples         = [
        {
            "code_input": {"molecule": "butane", "analysis_type": "newman"},
            "text_input": {"input_params": "butane newman"},
            "output": {"result": "Newman projection of butane... anti most stable... gauche +0.9 kcal/mol..."},
        },
        {
            "code_input": {"molecule": "methylcyclohexane", "analysis_type": "chair"},
            "text_input": {"input_params": "methylcyclohexane chair"},
            "output": {"result": "Chair conformation: equatorial methyl favored by 1.8 kcal/mol..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_database()

    def _build_database(self):
        """Build conformational analysis database."""
        
        # === Torsional Energies (kcal/mol) ===
        # Energy cost for eclipsing interactions
        self.eclipsing_energies = {
            "H-H": 1.0,       # H eclipsing H in ethane
            "H-CH3": 1.4,     # H eclipsing CH3 in propane/ethane-like
            "CH3-CH3": 3.6,   # Two methyl groups eclipsing each other
            "H-F": 1.5,      # Estimated
            "H-OH": 1.0,
            "H-NH2": 1.2,
            "CH3-F": 3.5,
            "CH3-OH": 3.0,
            "CH3-Ph": 3.0,   # Methyl eclipsing phenyl
            "Ph-Ph": 6.0,    # Two phenyl groups (very large)
        }

        # Gauche interactions (staggered but not anti)
        self.gauche_energies = {
            "CH3-CH3": 0.9,   # Butane gauche penalty
            "CH3-OH": 0.5,    # Moderate
            "CH3-F": 0.5,
            "CH3-Cl": 0.5,
            "OH-OH": 0.5,     # If applicable
            "large-large": 1.5-3.0,  # Two large groups
        }

        # === Newman Projections Database ===
        self.newman_systems = {
            "ethane": {
                "formula": "CH3-CH3",
                "conformations": [
                    {
                        "name": "Staggered",
                        "dihedral_angle": "60 degrees",
                        "energy": 0.0,  # Reference
                        "description": "Most stable. C-H bonds offset by 60deg.",
                        "projection": """
     H       H         H
      \\     |        /
       C --- C       (view along C-C bond)
      /     |        \\
     H       H        H
""",
                        "torsional_strain": "~0",
                    },
                    {
                        "name": "Eclipsed",
                        "dihedral_angle":"0 degrees",
                        "energy": 3.0,  # ~3 kcal/mol barrier
                        "description": "Transition state / least stable. C-H bonds aligned.",
                        "projection": """
     H       H
      \\   /
       C = C       (eclipsed)
      /   \\
     H       H
""",
                        "torsional_strain": "3.0 kcal/mol (3 × H/H eclipsing)",
                    },
                ],
                "barrier_height": "3.0 kcal/mol (rotation barrier)",
                "note": "Ethane rotation is essentially free at room temperature (barrier easily overcome).",
            },
            "butane": {
                "formula": "CH3-CH2-CH2-CH3",
                "conformations": [
                    {
                        "name": "Anti (anti-periplanar)",
                        "dihedral_angle":"180 degrees",
                        "energy": 0.0,
                        "description": "**Most stable conformation.** The two methyl groups are 180deg apart — maximum separation.",
                        "projection": """
      CH3          H
        |          /
    H3C-C ------ C-H    (looking down C2-C3 bond)
        |          \\
        H         CH3
""",
                        "strain_breakdown": {"torsional": "minimal (staggered)", "steric": "none (methyls far apart)", "total": "0.0 (reference)"},
                    },
                    {
                        "name": "Gauche (+synclinal)",
                        "dihedral_angle":"+-60 degrees",
                        "energy": 0.9,
                        "description": "Local minimum. Methyl groups 60deg apart — some steric repulsion.",
                        "projection": """
      CH3         H
        |        /
    H3C-C ----- C-H
        |       \\
        H       CH3     (dihedral ≈ +60deg)
""",
                        "strain_breakdown": {"torsional": "low (staggered)", "steric": "+0.9 (gauche CH3/CH3)", "total": "+0.9 kcal/mol"},
                    },
                    {
                        "name": "Eclipsed (partially eclipsed, CH3/H)",
                        "dihedral_angle":"+-120 degrees",
                        "energy": 3.4,
                        "description": "Transition state region. One CH3 eclipses an H; the other CH3 is gauche to H.",
                        "strain_breakdown": {"torsional": "3.0 (eclipsing)", "steric": "~0.4", "total": "~3.4 kcal/mol"},
                    },
                    {
                        "name": "Fully Eclipsed (syn-periplanar)",
                        "dihedral_angle":"0 degrees",
                        "energy": 4.5-5.0,
                        "description": "**Highest energy.** Both methyl groups eclipse each other directly.",
                        "projection": """
      CH3       CH3
        \\     /
         C = C        (fully eclipsed)
        /     \\
       H        H
""",
                        "strain_breakdown": {"torsional": "3.0 (eclipsing H/H)", "steric": "1.5-2.0 (CH3/CH3 eclipsing)", "total": "~4.5-5.0 kcal/mol"},
                    },
                ],
                "energy_profile": "Anti(0) → Gauche(+0.9) → Eclipsed CH3/H(+3.4) → Syn(+4.9) → Eclipsed CH3/H(+3.4) → Gauche(+0.9) → Anti(0)",
                "population_at_25C": "Anti: ~70%, Gauche: ~30% (each), Eclipsed: negligible",
                "note": "The gauche population (~30%) is significant! This explains why some reactions proceed through gauche conformers.",
            },
            "propane": {
                "formula": "CH3-CH2-CH3",
                "note": "Similar to ethane but with slightly higher barrier (~3.3 kcal/mol) due to CH3/H eclipsing (vs H/H in ethane).",
                "barrier": "~3.3 kcal/mol",
            },
        }

        # === Cyclohexane Conformations ===
        self.cyclohexane_conformations = {
            "chair": {
                "energy": 0.0,  # Reference
                "strain_energy": "0.0 kcal/mol (strain-free!)",
                "description": "**The ideal conformation.** All bonds staggered, no angle strain (109.5deg ideal), no torsional strain. This is why cyclohexane is the most common carbocycle in nature.",
                "features": [
                    "All C-C-C bond angles ≈ 111deg (close to tetrahedral 109.5deg)",
                    "All adjacent C-H bonds are staggered (no torsional strain)",
                    "12 H atoms: 6 axial (alternating up/down), 6 equatorial (outward)",
                    "Two distinct chair forms interconvert via ring flip (ΔG‡ ≈ 10-12 kcal/mol)",
                ],
                "axial_equatorial": "During ring flip: all axial ↔ equatorial and vice versa",
                "chair_ascii": """
        H(ax,up)   H(eq)
             \\   /
              C1 --- C6
             /   \\     /   \\
    H(eq) C2       C5   H(ax,down)
             \\   /     \\   /
              C3 --- C4
             /   \\
    H(ax,down) H(eq)
    
    Chair cyclohexane (axial = vertical, equatorial = outward)
""",
            },
            "boat": {
                "energy": 6.5,  # kcal/mol above chair
                "strain_energy": "~6.5 kcal/mol (mainly flagpole steric + torsional)",
                "description": "Higher energy conformation. Has eclipsed bonds and 'flagpole' H-H interaction.",
                "problems": ["Torsional strain from eclipsed bonds along the 'sides'", "Flagpole interaction: two H atoms at the bow/stern are forced close together (~1.83 Å)"],
                "not_a_minimum": "Boat is actually a transition state between twist-boats, not a true minimum.",
            },
            "twist_boat": {
                "energy": 5.5,  # kcal/mol above chair
                "strain_energy": "~5.5 kcal/mol",
                "description": "Local minimum (slightly lower than boat). Relieves some flagpole and torsional strain by twisting.",
                "population": "~<1% at room temperature (but important for ring inversion pathway)",
            },
            "half_chair": {
                "energy": 10.5,  # Transition state for ring flip (approx)
                "description": "Transition state between chair and twist-boat during ring inversion. Highest point on the conformational landscape.",
                "role": "TS for chair-chair interconversion (ring flip)",
            },
        }

        # === Disubstituted Cyclohexane Analysis ===
        self.disubstituted_cyclohexane = {
            "cis-1,2-dimethyl": {
                "dieq_deia": "one eq, one ax (always)",
                "more_stable": "diaxial is very unstable (1,3-diaxial interactions); diequatorial impossible for cis-1,2",
                "energy_difference": "N/A (only one reasonable chair form: one-up-one-down)",
                "a_value_cost": "~1.8 kcal/mol (one axial Me)",
                "ring_flip_product": "Same arrangement (one eq, one ax) — just which one is axial flips",
            },
            "trans-1,2-dimethyl": {
                "forms": ["diequatorial (most stable)", "diaxial (very unstable)"],
                "more_stable": "diequatorial",
                "energy_difference": "~2 × 1.8 = ~3.6 kcal/mol (two A-values)",
                "population": ">99% diequatorial at 25degC",
            },
            "cis-1,3-dimethyl": {
                "forms": ["diequatorial (most stable)", "diaxial (unstable)"],
                "more_stable": "diequatorial",
                "energy_difference": "~2 × 1.8 = ~3.6 kcal/mol",
            },
            "trans-1,3-dimethyl": {
                "dieq_deia": "one eq, one ax (always)",
                "energy_penalty": "~1.8 kcal/mol (one axial Me)",
            },
            "cis-1,4-dimethyl": {
                "forms": ["diequatorial (most stable)", "diaxial"],
                "more_stable": "diequatorial",
                "energy_difference": "~3.6 kcal/mol",
            },
            "trans-1,4-dimethyl": {
                "dieq_deia": "one eq, one ax",
                "energy_penalty": "~1.8 kcal/mol",
            },
        }

        # 1,3-Diaxial interaction energies
        self.diaxial_interactions = {
            "H/H": 0.40,      # Reference: axial H ... axial H (1,3)
            "H/CH3": 0.90,    # Axial H ... axial CH3
            "CH3/CH3": 1.80,  # Axial CH3 ... axial CH3 (= A-value of methyl roughly)
            "H/OH": 0.45,
            "H/F": 0.30,
            "H/Cl": 0.35,
            "H/Br": 0.30,
            "CH3/OH": 1.0,
            "CH3/F": 0.85,
            "CH3/Cl": 0.95,
            "CH3/t-Bu": 3.5,   # Very severe!
            "t-Bu/t-Bu": ">10", # Extremely unfavorable
        }

    def _run_base(self, molecule: str, analysis_type: str = "newman") -> str:
        """Perform conformational analysis."""
        mol_lower = molecule.lower().strip()
        atype = analysis_type.lower().strip()
        
        parts = [f"## Conformational Analysis: `{molecule}`\n"]
        parts.append(f"**Analysis Type:** {atype}\n")

        # Determine what kind of analysis to do
        is_cyclic = any(w in mol_lower for w in ["cyclo", "ring", "chair", "hexane", "pyranose"])
        is_newman_target = any(w in mol_lower for w in ["ethane", "butane", "propane", "pentane", "neopentane"])

        if atype == "all" or (not is_cyclic and not is_newman_target):
            # Do comprehensive analysis
            result = []
            if is_cyclic:
                result.extend(self._analyze_cyclic(mol_lower))
            else:
                result.extend(self._analyze_acyclic(mol_lower))
            parts += result
        elif atype in ("chair", "cyclo", "cyclic") or is_cyclic:
            parts += self._analyze_cyclic(mol_lower)
        else:
            parts += self._analyze_acyclic(mol_lower)

        return "\n".join(parts)

    def _analyze_acyclic(self, mol):
        """Analyze acyclic molecule (Newman projections)."""
        parts = []
        
        if "butane" in mol:
            sys_data = self.newman_systems["butane"]
            parts.append(f"### 🔄 Newman Projection Analysis: **Butane** ({sys_data['formula']})\n")
            
            for conf in sys_data["conformations"]:
                parts.append(f"#### {conf['name']} (dihedral: {conf['dihedral_angle']})\n")
                parts.append(f"- **Energy:** ΔG = **+{conf['energy']} kcal/mol**\n")
                parts.append(f"- {conf['description']}\n")
                if conf.get("projection"):
                    parts.append(f"```\n{conf['projection']}\n```\n")
                if conf.get("strain_breakdown"):
                    sb = conf["strain_breakdown"]
                    parts.append(f"- **Strain Breakdown:**\n")
                    for k, v in sb.items():
                        parts.append(f"  - {k}: {v}\n")
            
            parts.append(f"\n### 📈 Energy Profile\n")
            parts.append(f"{sys_data.get('energy_profile', 'N/A')}\n")
            parts.append(f"\n**Population at 25degC:** {sys_data.get('population_at_25C', 'N/A')}\n")
            parts.append(f"> {sys_data.get('note', '')}\n")
            
        elif "ethane" in mol:
            sys_data = self.newman_systems["ethane"]
            parts.append(f"### 🔄 Newman Projection: **Ethane**\n")
            for conf in sys_data["conformations"]:
                parts.append(f"#### {conf['name']}\n")
                parts.append(f"- **Energy:** ΔG = **+{conf['energy']} kcal/mol**\n")
                parts.append(f"- {conf['description']}\n")
                if conf.get("projection"):
                    parts.append(f"```\n{conf['projection']}\n```\n")
            parts.append(f"\n**Rotation Barrier:** {sys_data['barrier_height']}\n")
            parts.append(f"> {sys_data.get('note', '')}\n")
            
        else:
            # Generic acyclic analysis
            parts.append("### 🔄 General Acyclic Conformation Analysis\n")
            parts.append("""
#### Key Principles for Acyclic Conformations

| Principle | Description |
|---|---|
| **Staggered > Eclipsed** | Staggered is always favored by ~3 kcal/mol per eclipsed interaction |
| **Anti > Gauche** | Anti (180deg dihedral) favored over gauche (±60deg) for larger groups |
| **Gauche Penalty** | CH3/CH3: +0.9 kcal/mol; larger groups: +1-3 kcal/mol |
| **Eclipsing Cost** | H/H: 1.0; H/CH3: 1.4; CH3/CH3: 3.6 kcal/mol |
| **1,3-Syn-axial** | Not applicable to acyclic, but analogous to gauche |

#### Common Systems

**Propane (CH3-CH2-CH3):**
- Barrier: ~3.3 kcal/mol (CH3/H eclipsing)
- Similar to ethane but slightly higher

**Neopentane ((CH3)4C):**
- Very crowded central carbon
- Three methyls on each side create significant gauche interactions
- Prefers all-staggered arrangements

**1,2-Dichloroethane:**
- Anti favored over gauche (Cl/Cl gauche penalty)
- In gas phase: anti dominant
- In solution: can have more gauche due to solvation effects

""")
            parts.append(f"> For specific molecules like **butane**, use `analysis_type='newman'` with explicit name.\n")

        return parts

    def _analyze_cyclic(self, mol):
        """Analyze cyclic molecule (cyclohexane-based)."""
        parts = []

        # Chair conformation basics
        parts.append("### 🪑 Cyclohexane Chair Conformation\n")
        chair = self.cyclohexane_conformations["chair"]
        parts.append(f"- **Strain Energy:** {chair['strain_energy']}")
        parts.append(f"- {chair['description']}\n")
        parts.append("**Key Features:**\n")
        for f in chair["features"]:
            parts.append(f"- {f}\n")
        if chair.get("chair_ascii"):
            parts.append(f"```\n{chair['chair_ascii']}\n```\n")

        # Other conformations
        parts.append("\n### 📊 Cyclohexane Conformation Energy Landscape\n")
        parts.append("| Conformation | Energy (kcal/mol) | Type |")
        parts.append("|---|---|---|")
        for key, data in self.cyclohexane_conformations.items():
            role = data.get("role", "Minimum" if data["energy"] < 7 else "Transition State")
            emoji = "🪑" if key == "chair" else "⛵" if "boat" in key else "🔀"
            parts.append(f"| {emoji} {key.replace('_', ' ').title()} | +{data['energy']:.1f} | {role} |")
        
        parts.append("\n> **Ring flip barrier:** ΔG‡ ≈ **10-12 kcal/mol** (fast at room temperature)\n")

        # Check for disubstituted patterns
        disub_match = None
        for key, data in self.disubstituted_cyclohexane.items():
            if key.replace("-", "").replace(",", "") in mol.replace("-", "").replace(",", ""):
                disub_match = (key, data)
                break
        
        if disub_match or "dimethyl" in mol or "substituted" in mol:
            parts.append("\n### 🔬 Disubstituted Cyclohexane Analysis\n")
            if disub_match:
                key, data = disub_match
                parts.append(f"#### **{key.title()}**\n")
                for k, v in data.items():
                    parts.append(f"- **{k.replace('_', ' ').title()}:** {v}\n")
            else:
                parts.append("| Isomer | Most Stable Form | Energy Penalty |")
                parts.append("|---|---|---|")
                for key, data in self.disubstituted_cyclohexane.items():
                    penalty = data.get("energy_difference", data.get("energy_penalty", "N/A"))
                    stable = data.get("more_stable", data.get("forms", ["N/A"])[0] if isinstance(data.get("forms"), list) else "N/A")
                    parts.append(f"| {key} | {stable} | {penalty} |")

        # 1,3-Diaxial interactions reference
        parts.append("\n### ⚡ 1,3-Diaxial Interaction Energies\n")
        parts.append("| Interaction | Energy (kcal/mol) |")
        parts.append("|---|---|")
        sorted_da = sorted(self.diaxial_interactions.items(), key=lambda x: float(str(x[1]).replace('~','').replace('>','')) if isinstance(x[1], (int, float)) else 0, reverse=True)
        for pair, energy in sorted_da[:12]:
            parts.append(f"| Ax-{pair} | {energy} |")

        return parts

    def _run_text(self, input_params: str) -> str:
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Please provide a molecule name. Example: 'butane newman' or 'cyclohexane chair'")
        molecule = parts[0]
        atype = parts[1] if len(parts) > 1 else "newman"
        return self._run_base(molecule, atype)
