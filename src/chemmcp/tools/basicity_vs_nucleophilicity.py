import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive data: common bases/nucleophiles with pKa (conjugate acid) and nucleophilicity trends
_NUCLEOPHILE_DATA = {
    # (pKa of conjugate acid, protic solvent ranking, aprotic solvent ranking, notes)
    "I-": {"pka": -10, "protic_rank": 1, "aprotic_rank": 4, "type": "weak_base_strong_nucleophile", "notes": "Excellent nucleophile, very weak base. Large size, polarizable."},
    "Br-": {"pka": -9, "protic_rank": 2, "aprotic_rank": 3, "type": "weak_base_strong_nucleophile", "notes": "Very good nucleophile, weak base. More polarizable than Cl-."},
    "Cl-": {"pka": -7, "protic_rank": 3, "aprotic_rank": 2, "type": "weak_base_moderate_nucleophile", "notes": "Good nucleophile in polar aprotic solvents; moderate in protic."},
    "F-": {"pka": 3.2, "protic_rank": 8, "aprotic_rank": 1, "type": "strong_base_strong_nucleophile_aprotic", "notes": "Weak nucleophile in protic solvents (H-bonding); excellent in aprotic."},
    "HO-": {"pka": 15.7, "protic_rank": 5, "aprotic_rank": 6, "type": "strong_base_good_nucleophile", "notes": "Strong base and good nucleophile. Favors E2 over SN2 with hindered substrates."},
    "CH3O-": {"pka": 15.5, "protic_rank": 6, "aprotic_rank": 7, "type": "strong_base_good_nucleophile", "notes": "Similar to HO- but more basic and sterically similar."},
    "t-BuO-": {"pka": 17, "protic_rank": 9, "aprotic_rank": 9, "type": "strong_base_poor_nucleophile", "notes": "Strong base but very poor nucleophile (steric bulk). Promotes E2."},
    "CH3COO-": {"pka": 4.76, "protic_rank": 7, "aprotic_rank": 5, "type": "weak_base_weak_nucleophile", "notes": "Weak base and weak nucleophile. Rarely used as nucleophile."},
    "CN-": {"pka": 9.2, "protic_rank": 4, "aprotic_rank": 8, "type": "moderate_base_good_nucleophile", "notes": "Good nucleophile, moderate base. Small and not highly solvated."},
    "N3-": {"pka": 4.7, "protic_rank": None, "aprotic_rank": None, "type": "moderate_base_good_nucleophile", "notes": "Good nucleophile, linear shape reduces steric hindrance."},
    "NH3": {"pka": 9.2, "protic_rank": None, "aprotic_rank": None, "type": "weak_base_moderate_nucleophile", "notes": "Neutral nucleophile. Moderate nucleophilicity, weak base."},
    "CH3NH2": {"pka": 10.6, "protic_rank": None, "aprotic_rank": None, "type": "moderate_base_moderate_nucleophile", "notes": "Stronger base than NH3, good nucleophile for SN2."},
    "(CH3)2NH": {"pka": 10.7, "protic_rank": None, "aprotic_rank": None, "type": "moderate_base_moderate_nucleophile", "notes": "More basic/nucleophilic than NH3 but with steric effects."},
    "(CH3)3N": {"pka": 9.8, "protic_rank": None, "aprotic_rank": None, "type": "moderate_base_poor_nucleophile", "notes": "Sterically hindered; poor nucleophile despite moderate basicity."},
    "H2O": {"pka": 15.7, "protic_rank": None, "aprotic_rank": None, "type": "very_weak_base_very_weak_nucleophile", "notes": "Very weak base and nucleophile. Only reacts with very reactive substrates."},
    "CH3S-": {"pka": 6.5, "protic_rank": None, "aprotic_rank": None, "type": "weak_base_excellent_nucleophile", "notes": "Excellent nucleophile (large, polarizable S). Weak base (pKa of CH3SH)."},
    "PhS-": {"pka": 6.5, "protic_rank": None, "aprotic_rank": None, "type": "weak_base_excellent_nucleophile", "notes": "Soft nucleophile, excellent for SN2 with primary/secondary alkyl halides."},
    "RS-": {"pka": ~8-11, "protic_rank": None, "aprotic_rank": None, "type": "weak_base_excellent_nucleophile", "notes": "Thiolates are excellent soft nucleophiles, weak bases."},
}

_SOLVENT_EFFECTS = {
    "protic": {
        "description": "Protic solvents (water, alcohols) form hydrogen bonds with anions.",
        "effect": "Strong H-bonding solvates small anions (F-, HO-) heavily, reducing their nucleophilicity. Large/polarizable ions (I-, Br-) are less affected.",
        "nucleophilicity_order": "I- > Br- > Cl- > F- (in protic solvents)",
        "basicity_order": "F- > HO- > CH3O- > CN- > N3- > NH3 > I-",
    },
    "aprotic_polar": {
        "description": "Polar aprotic solvents (DMSO, acetone, DMF, acetonitrile) do NOT H-bond to anions.",
        "effect": "Small 'hard' anions are poorly solvated and become VERY nucleophilic. Order reverses from protic.",
        "nucleophilicity_order": "F- > Cl- > Br- > I- (in polar aprotic solvents)",
    },
}


@ChemMCPManager.register_tool
class BasicityVsNucleophilicity(BaseTool):
    """
    分析碱性与亲核性关系的工具。
    基于溶剂效应、结构因素、硬度/软度等化学原理，分析给定碱/亲核体的碱性强度与亲核性强度的关系。
    """
    __version__      = "0.1.0"
    name             = "BasicityVsNucleophilicity"
    func_name        = "analyze_basicity_nucleophilicity"
    description      = "Analyze the relationship between basicity and nucleophilicity for a given species, considering solvent effects, steric factors, and HSAB principles."
    implementation_description = "Uses embedded chemical data tables of pKa values, nucleophilicity rankings in different solvents, HSAB theory, and steric factors to provide comprehensive analysis."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Reaction"]
    tags             = ["Basicity", "Nucleophilicity", "Organic Chemistry", "HSAB Theory", "Solvent Effects"]
    required_envs    = []

    code_input_sig   = [
        ("species", "str", "N/A", "The chemical species to analyze (e.g., 'HO-', 'CH3O-', 'I-', 't-BuO-', 'CN-')."),
        ("solvent", "str", "protic", "Solvent environment: 'protic', 'aprotic_polar', or 'nonpolar'."),
        ("compare_with", "str", "None", "Optional second species to compare with (e.g., 'Br-'). Use 'None' for no comparison."),
    ]

    text_input_sig   = [
        ("input_text", "str", "N/A", "Input string: 'species [solvent] [compare_species]'. Example: 'HO- protic CH3O-'"),
    ]

    output_sig       = [
        ("result", "str", "Detailed analysis string covering basicity, nucleophilicity, solvent effects, and comparison."),
    ]

    examples         = [
        {
            "code_input": {
                "species": "HO-",
                "solvent": "protic",
                "compare_with": "t-BuO-"
            },
            "text_input": {
                "input_text": "HO- protic t-BuO-"
            },
            "output": {
                "result": """## Basicity vs Nucleophilicity Analysis

### Species: HO- (Hydroxide)

**Basicity:** Strong base (pKa of conjugate acid H₂O = 15.7)
**Nucleophilicity in protic solvent:** Good (ranked ~5th among common nucleophiles)
**Type:** Strong base / Good nucleophile

---

### Key Factors:

1. **Basicity**: HO- is a strong base — it eagerly accepts protons.
2. **Nucleophilicity**: In protic solvents, HO- is moderately nucleophilic because it is strongly solvated by H-bonds, which partially neutralizes its reactivity.
3. **Sterics**: Small size → minimal steric hindrance → can attack crowded electrophilic centers.

### Solvent Effect:
In **protic solvents**, HO- is strongly H-bonded/solvated, reducing its nucleophilicity relative to its basicity.

### Comparison: HO- vs t-BuO-

| Property | HO- | t-BuO- |
|----------|-----|--------|
| Basicity | Strong (pKa 15.7) | Stronger (pKa ~17) |
| Nucleophilicity | Good | Very Poor |
| Steric Bulk | Small | Very Bulky |
| Typical Role | SN2 + E2 competitor | E2 promoter |

**Conclusion:** t-BuO- is MORE basic but LESS nucleophilic than HO- due to extreme steric hindrance. t-BuO- promotes elimination (E2); HO- can do both SN2 and E2."""
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, species: str, solvent: str = "protic", compare_with: str = "None") -> str:
        """Core logic: analyze basicity vs nucleophilicity."""
        species = species.strip()
        solvent = solvent.strip().lower()
        cmp = compare_with.strip() if compare_with and compare_with.upper() != "NONE" else None

        lines = []
        lines.append("## Basicity vs Nucleophilicity Analysis\n")

        # Look up primary species
        data = _NUCLEOPHILE_DATA.get(species)
        if data is None:
            # Try case-insensitive match
            for k, v in _NUCLEOPHILE_DATA.items():
                if k.lower() == species.lower():
                    species = k
                    data = v
                    break

        if data is None:
            return f"## Analysis for '{species}'\n\n⚠️ Species '{species}' not found in database. Available species:\n" + \
                   ", ".join(sorted(_NUCLEOPHILE_DATA.keys()))

        # Primary analysis
        type_desc = {
            "weak_base_strong_nucleophile": "Weak base / Strong nucleophile",
            "weak_base_moderate_nucleophile": "Weak base / Moderate nucleophile",
            "weak_base_weak_nucleophile": "Weak base / Weak nucleophile",
            "strong_base_good_nucleophile": "Strong base / Good nucleophile",
            "strong_base_poor_nucleophile": "Strong base / Poor nucleophile",
            "strong_base_strong_nucleophile_aprotic": "Strong base / Strong nucleophile (in aprotic only)",
            "moderate_base_good_nucleophile": "Moderate base / Good nucleophile",
            "moderate_base_moderate_nucleophile": "Moderate base / Moderate nucleophile",
            "moderate_base_poor_nucleophile": "Moderate base / Poor nucleophile",
            "very_weak_base_very_weak_nucleophile": "Very weak base / Very weak nucleophile",
            "weak_base_excellent_nucleophile": "Weak base / Excellent nucleophile (soft/ polarizable)",
        }

        lines.append(f"### Species: {species}")
        lines.append(f"**Basicity:** {self._basicity_desc(data['pka'])}")
        if data.get('protic_rank') is not None:
            lines.append(f"**Nucleophilicity (protic):** Ranked #{data['protic_rank']} among common nucleophiles")
        if data.get('aprotic_rank') is not None:
            lines.append(f"**Nucleophilicity (aprotic):** Ranked #{data['aprotic_rank']} among common nucleophiles")
        lines.append(f"**Classification:** {type_desc.get(data['type'], data['type'])}")
        lines.append(f"**Notes:** {data['notes']}")
        lines.append("")

        # Key factors
        lines.append("---\n### Key Factors:")
        factors = self._analyze_factors(species, data)
        for f in factors:
            lines.append(f"- **{f[0]}:** {f[1]}")
        lines.append("")

        # Solvent effect
        lines.append("### Solvent Effect:")
        if solvent in ("protic", "water", "alcohol"):
            lines.append(_SOLVENT_EFFECTS["protic"]["effect"])
        elif solvent in ("aprotic_polar", "aprotic", "dmso", "dmf", "acetone"):
            lines.append(_SOLVENT_EFFECTS["aprotic_polar"]["effect"])
        else:
            lines.append("Nonpolar solvents have minimal differential effect on nucleophilicity vs basicity.")
        lines.append("")

        # Comparison
        if cmp:
            cmp_data = _NUCLEOPHILE_DATA.get(cmp)
            if cmp_data is None:
                for k, v in _NUCLEOPHILE_DATA.items():
                    if k.lower() == cmp.lower():
                        cmp = k
                        cmp_data = v
                        break

            lines.append(f"---\n### Comparison: {species} vs {cmp}\n")
            if cmp_data:
                lines.append("| Property | {} | {} |".format(species, cmp))
                lines.append("|----------|------|------|")
                lines.append("| Basicity (pKa of conj. acid) | {} | {} |".format(data['pka'], cmp_data['pka']))
                b1 = "Strong" if data['pka'] > 7 else ("Moderate" if data['pka'] > 0 else "Weak")
                b2 = "Strong" if cmp_data['pka'] > 7 else ("Moderate" if cmp_data['pka'] > 0 else "Weak")
                lines.append("| Basic Strength | {} | {} |".format(b1, b2))
                r1 = data.get('protic_rank', 'N/A')
                r2 = cmp_data.get('protic_rank', 'N/A')
                lines.append("| Nucl. Rank (protic) | {} | {} |".format(r1, r2))
                t1 = type_desc.get(data['type'], data['type'])
                t2 = type_desc.get(cmp_data['type'], cmp_data['type'])
                lines.append("| Classification | {} | {} |".format(t1, t2))
                lines.append("")
                lines.append("**Conclusion:** " + self._conclude_comparison(species, data, cmp, cmp_data))
            else:
                lines.append(f"⚠️ Comparison species '{cmp}' not found in database.")

        lines.append("\n---\n*Analysis based on standard organic chemistry principles (HSAB theory, solvent effects, steric factors).*")
        return "\n".join(lines)

    def _run_text(self, input_text: str) -> str:
        parts = input_text.strip().split()
        species = parts[0] if len(parts) > 0 else "HO-"
        solvent = parts[1] if len(parts) > 1 else "protic"
        cmp = parts[2] if len(parts) > 2 else "None"
        return self._run_base(species, solvent, cmp)

    def _basicity_desc(self, pka):
        if pka > 12:
            return f"Strong base (pKa of conjugate acid = {pka})"
        elif pka > 5:
            return f"Moderate base (pKa of conjugate acid = {pka})"
        elif pka > 0:
            return f"Weak base (pKa of conjugate acid = {pka})"
        else:
            return f"Very weak base (pKa of conjugate acid = {pka}, conjugate acid is strong)"

    def _analyze_factors(self, species, data):
        factors = []
        pka = data['pka']
        # Charge factor
        if species.endswith('-'):
            factors.append(("Charge", "Anionic species — generally more basic and nucleophilic than neutral counterparts."))
        elif species.endswith(')') or '(' in species:
            factors.append(("Charge", "Neutral species — typically less basic/nucleophilic than anions."))

        # Size / Polarizability
        if any(x in species for x in ['I', 'Br', 'S']):
            factors.append(("Size/Polarizability", "Large, polarizable atom present — enhances nucleophilicity (soft nucleophile) without increasing basicity proportionally."))
        elif any(x in species for x in ['F', 'O', 'N']) and 't-Bu' not in species:
            factors.append(("Size/Polarizability", "Small, hard atom — basicity correlates better with nucleophilicity in aprotic solvents."))

        # Sterics
        if 't-Bu' in species or '(CH3)3' in species:
            factors.append(("Steric Hindrance", "Highly branched/bulky structure — severely hinders SN2 attack despite strong basicity. Favors elimination."))
        elif '(CH3)2' in species or 'CH3' in species:
            factors.append(("Steric Hindrance", "Some steric bulk present — may reduce nucleophilicity relative to smaller analogs."))

        # HSAB
        if pka < 0 and any(x in species for x in ['I', 'Br', 'Cl', 'S']):
            factors.append(("HSAB Character", "Soft nucleophile — prefers soft electrophiles (e.g., alkyl iodides, allylic/benzylic positions)."))
        elif pka > 10 and any(x in species for x in ['O', 'F', 'N']):
            factors.append(("HSAB Character", "Hard nucleophile — prefers hard electrophiles (e.g., methyl halides, carbonyl carbons)."))

        return factors

    def _conclude_comparison(self, s1, d1, s2, d2):
        b_diff = d1['pka'] - d2['pka']
        r1 = d1.get('protic_rank', 99)
        r2 = d2.get('protic_rank', 99)

        conclusions = []
        if abs(b_diff) > 3:
            stronger = s1 if b_diff > 0 else s2
            weaker = s2 if b_diff > 0 else s1
            conclusions.append(f"{stronger} is STRONGER BASE than {weaker}")

        if r1 != 99 and r2 != 99 and r1 != r2:
            better_nu = s1 if r1 < r2 else s2
            worse_nu = s2 if r1 < r2 else s1
            conclusions.append(f"{better_nu} is BETTER NUCLEOPHILE (protic) than {worse_nu}")

        # Check for inversion
        if b_diff > 2 and r1 != 99 and r2 != 99 and r1 > r2:
            conclusions.append("⚠️ INVERSION: The stronger base is the WORSE nucleophile (likely due to sterics or solvation)")
        elif b_diff < -2 and r1 != 99 and r2 != 99 and r1 < r2:
            conclusions.append("⚠️ INVERSION: The weaker base is the BETTER nucleophile (likely due to polarizability)")

        if not conclusions:
            return f"{s1} and {s2} show similar basicity-nucleophilicity profiles."
        return ". ".join(conclusions) + "."
