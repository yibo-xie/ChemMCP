import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class InductiveEffectAnalyzer(BaseTool):
    """
    诱导效应分析工具 - 分析分子中取代基的诱导效应对酸性、碱性、反应活性和其他性质的影响。
    通过σ键传递电子效应，随距离增加而衰减。
    """
    __version__ = "0.1.0"
    name             = "InductiveEffectAnalyzer"
    func_name        = "analyze_inductive_effect"
    description      = "Analyze inductive effects (electron withdrawal/donation through σ-bonds) on molecular properties: acidity, basicity, reactivity, and spectroscopic shifts."
    implementation_description = "Knowledge-based system using electronegativity differences, field effects, and distance-dependent attenuation to predict how substituents affect molecular properties through σ-bond inductive effects."
    oss_dependencies = []
    services_and_software = []
    categories       = ["Molecule"]
    tags             = ["Inductive Effect", "Substituent Effects", "Physical Organic Chemistry", "Acidity", "Reactivity"]
    required_envs    = []

    code_input_sig   = [
        ("molecule", "str", "N/A", "Molecule name, SMILES, or structural description (e.g., 'chloroacetic acid', 'trifluoroethanol', '4-nitrobenzoic acid')."),
        ("focus_property", "str", "general", "Property to focus on: 'acidity', 'basicity', 'reactivity', 'nmr_shift', 'dipole', or 'general' for overall analysis."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'molecule [focus_property]'. Example: 'chloroacetic acid acidity'."),
    ]

    output_sig       = [
        ("result", "str", "Detailed inductive effect analysis with predicted property changes, quantitative estimates where possible."),
    ]

    examples         = [
        {
            "code_input": {"molecule": "chloroacetic acid", "focus_property": "acidity"},
            "text_input": {"input_params": "chloroacetic acid acidity"},
            "output": {"result": "Analysis of -I effect of Cl on acetic acid acidity..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_database()

    def _build_database(self):
        """Build inductive effect database."""
        
        # Electronegativity values (Pauling scale)
        self.en = {
            "F": 3.98, "O": 3.44, "Cl": 3.16, "N": 3.04, "Br": 2.96,
            "S": 2.58, "C": 2.55, "I": 2.66, "H": 2.20,
            "P": 2.19, "B": 2.04, "Si": 1.90, "Ge": 2.01,
            "As": 2.18, "Se": 2.55, "Te": 2.10,
            "metal_na": 0.93, "metal_k": 0.82, "metal_li": 0.98,
        }

        # Inductive effect strength ranking (relative to H)
        # Positive value = electron-withdrawing (-I), negative = electron-donating (+I)
        self.i_effects = {
            # Strong -I groups (complete)
            "NR3+": +1.30, "NH3+": +1.20, "NO2": +0.67, "CN": +0.58,
            "SO2R": +0.72, "SO3R": +0.80, "CF3": +0.54, "CCl3": +0.50,
            "F": (+0.54, "direct"), "Br": (+0.47, "direct"), "Cl": (+0.47, "direct"),
            "I": (+0.40, "direct"),
            "COOR": +0.45, "COR": +0.50, "CHO": +0.42, "CONH2": +0.28,
            "OH": (+0.33, "direct"), "OR": (+0.30, "direct"),
            "NHCOR": +0.25, "NHR": (+0.08, "direct"),
            "C6H5": +0.06, "CH=CH2": +0.05, "C≡CH": +0.18,
            "H": 0.0,
            # +I groups (negative values)
            "D": -0.02, "SiR3": -0.15, "GeR3": -0.10,
            "C(CH3)3": -0.10, "CH(CH3)2": -0.07, "CH2R": -0.04, "CH3": -0.04,
            "CR3" : -0.01,  # Alkyl general
        }

        # Distance attenuation factor (roughly 1/2^n per bond for simple -I)
        # More accurate: exponential decay with distance constant ~2-3 bonds
        
        # Known pKa shift data for calibration
        self.pka_shift_data = {
            # (parent, substituent, position): ΔpKa
            ("acetic acid", "F", "alpha"): -2.05,   # Fluoroacetic acid 2.66 vs 4.76
            ("acetic acid", "Cl", "alpha"): -1.90,  # Chloroacetic acid 2.86 vs 4.76
            ("acetic acid", "Cl", "beta"): -0.40,   # β-Chloropropionic acid ≈ 4.4
            ("acetic acid", "Cl", "gamma"): -0.20,  # γ-Chlorobutyric acid ≈ 4.5
            ("acetic acid", "Cl2", "alpha"): -3.86, # Dichloroacetic acid 1.26 vs 4.76 (approx additive)
            ("acetic acid", "Cl3", "alpha"): -4.06, # Trichloroacetic acid 0.70 vs 4.76
            ("phenol", "NO2", "para"): -2.85,       # p-Nitrophenol 7.15 vs 10.00
            ("phenol", "NO2", "meta"): -1.61,       # m-Nitrophenol 8.39 vs 10.00
            ("phenol", "NO2", "ortho"): -2.83,      # o-Nitrophenol 7.17 vs 10.00
            ("phenol", "Cl", "para"): -0.62,        # p-Chlorophenol 9.38 vs 10.00
            ("phenol", "Cl", "meta"): -0.98,        # m-Chlorophenol 9.02 vs 10.00
            ("phenol", "CH3", "para"): +0.26,       # p-Cresol 10.26 vs 10.00
            ("phenol", "OCH3", "para"): -0.73,      # Actually net effect is complex (pKa of guaiacol ≈ 9.98)
            ("benzoic acid", "NO2", "para"): -0.79, # p-Nitrobenzoic acid 3.41 vs 4.20
            ("benzoic acid", "NO2", "meta"): -0.71, # m-Nitrobenzoic acid 3.49 vs 4.20
            ("benzoic acid", "NH2", "para"): +0.79, # p-Aminobenzoic acid 4.99 vs 4.20
            ("benzoic acid", "OH", "para"): -0.16, # p-Hydroxybenzoic acid 4.54 vs 4.20 (complex)
            ("benzoic acid", "CF3", "para"): -0.43, # p-CF3-benzoic acid ≈ 3.77
            ("ethanol", "F", "alpha"): -1.5,        # 2-Fluoroethanol pKa ≈ 14.0 vs 15.9
            ("ethanol", "Cl", "alpha"): -1.0,       # 2-Chloroethanol pKa ≈ 14.9
            ("methanol", "CF3", "alpha"): -5.0,     # TFE pKa ≈ 12.5 vs 15.5
        }

        # NMR chemical shift effects (¹H, ppm relative to parent)
        self.nmr_inductive_shifts = {
            "alpha_to_EWG": "+2.0 to +4.0 ppm (deshielded)",
            "alpha_to_EDG": "-0.2 to -0.5 ppm (slightly shielded)",
            "beta_to_EWG": "+0.3 to +1.0 ppm",
            "gamma_to_EWG": "-0.1 to +0.2 ppm (often negligible or slight γ-gauche effect)",
        }

    def _identify_substituents(self, molecule: str) -> list:
        """Identify substituents and their positions from molecule description."""
        mol_lower = molecule.lower()
        substituents = []
        
        # Pattern matching for common substitution patterns
        patterns = [
            (r"fluoro|fluorine|\bf\b", "F", "unknown"),
            (r"chloro|chlorine|\bcl\b", "Cl", "unknown"),
            (r"bromo|bromine|\bbr\b", "Br", "unknown"),
            (r"iodo|iodine|\bi\b(?!pr)", "I", "unknown"),
            (r"trifluoro|cf3|tfa", "CF3", "unknown"),
            (r"trichloro|ccl3", "CCl3", "unknown"),
            (r"dichloro|cl2ch|chcl2", "CCl2H", "unknown"),
            (r"cyano|cyano.*group|cn|nitrile", "CN", "unknown"),
            (r"nitro|no2", "NO2", "unknown"),
            (r"methoxy|och3|meo", "OCH3", "unknown"),
            (r"hydroxy|oh", "OH", "unknown"),
            (r"amino|nh2", "NH2", "unknown"),
            (r"methyl|ch3", "CH3", "unknown"),
            (r"tert.butyl|t.bu|c\(ch3\)3", "C(CH3)3", "unknown"),
            (r"phenyl|ph", "Ph", "unknown"),
            (r"carboxyl|cooh", "COOH", "unknown"),
            (r"carbonyl|cor|c=o|ketone", "COR", "unknown"),
            (r"ester|coor|coome|cooet", "COOR", "unknown"),
            (r"sulfonyl|so2", "SO2R", "unknown"),
            (r"dimethylamino|nme2|n\(ch3\)2", "NMe2", "unknown"),
            (r"trimethylammonium|n\(ch3\)3\+", "N(CH3)3+", "unknown"),
            (r"para.*|p[- ]|4-", None, "para"),  # Position only
            (r"meta.*|m[- ]|3-", None, "meta"),
            (r"ortho.*|o[- ]|2-", None, "ortho"),
            (r"alpha|α|2-|β.position", None, "alpha"),
            (r"beta|β|3-", None, "beta"),
            (r"gamma|γ|4-", None, "gamma"),
        ]
        
        import re
        for pattern, sub, pos in patterns:
            if re.search(pattern, mol_lower, re.IGNORECASE):
                if sub:
                    substituents.append({"substituent": sub, "position": pos})
                else:
                    # Update position of last added substituent
                    if substituents:
                        substituents[-1]["position"] = pos

        return substituents

    def _estimate_pka_shift(self, base_pka: float, substituent: str, position: str) -> float:
        """Estimate pKa shift based on inductive effect."""
        # Get base -I strength
        i_strength = self.i_effects.get(substituent, 0)
        if isinstance(i_strength, tuple):
            i_strength = i_strength[0]
        
        if i_strength == 0:
            return 0.0
        
        # Distance attenuation
        distance_factors = {
            "alpha": 1.0,
            "beta": 0.35,
            "gamma": 0.12,
            "delta": 0.04,
            "ortho": 0.4,     # Aromatic ortho — both inductive + resonance proximity
            "meta": 0.25,     # Aromatic meta — mostly inductive
            "para": 0.15,     # Aromatic para — resonance dominates over induction
            "unknown": 0.3,   # Default guess
        }
        
        d_factor = distance_factors.get(position, 0.3)
        
        # For EWGs: lower pKa (more acidic) → negative shift
        # Scale factor calibrated roughly against experimental data
        # α-halogenated acids show ~2 pKa units per Cl at alpha position
        shift = -i_strength * d_factor * 4.0  # Empirical scaling
        
        return round(shift, 2)

    def _run_base(self, molecule: str, focus_property: str = "general") -> str:
        """Analyze inductive effects."""
        subs = self._identify_substituents(molecule)
        
        parts = [f"## Inductive Effect Analysis: `{molecule}`\n"]
        parts.append(f"**Focus Property:** {focus_property}\n")

        if not subs:
            parts.append("### ⚠️ No Substituents Detected\n")
            parts.append("Could not identify specific substituents from the molecule name.\n")
            parts.append("**Try more specific names:**\n")
            parts.append("- `chloroacetic acid` / `2-chloroethanol` / `trifluoroacetic acid`\n")
            parts.append("- `p-nitrophenol` / `m-chlorobenzoic acid` / `4-cyanobutanoic acid`\n")
            parts.append("- `α-fluoro ketone` / `β-keto ester` / `γ-amino alcohol`\n")
            return "\n".join(parts)

        parts += [f"### 🔍 Identified Substituents\n"]
        for s in subs:
            sub_name = s["substituent"]
            pos = s["position"]
            
            # Get EN difference info
            en_val = None
            atom_en = self.en.get(sub_name[0] if len(sub_name) > 0 else "")
            
            # Get -I/+I classification
            i_eff = self.i_effects.get(sub_name, 0)
            if isinstance(i_eff, tuple):
                i_eff = i_eff[0]
            
            if i_eff > 0.1:
                classification = f"**Electron-Withdrawing (-I)** (strength: {i_eff:+.2f})"
                direction = "pulls electron density away"
                color = "🔴"
            elif i_eff < -0.05:
                classification = f"**Electron-Donating (+I)** (strength: {i_eff:+.2f})"
                direction = "pushes electron density toward"
                color = "🟢"
            else:
                classification = "**Neutral/Very Weak** inductive effect"
                direction = "minimal electronic effect"
                color = "⚪"

            parts.append(f"- **{sub_name}** ({pos} position): {classification} {color}")
            parts.append(f"  - This group *{direction}* the reaction center via σ-bonds\n")

        # Property-specific analysis
        parts.append(f"\n### 📊 Effect on `{focus_property}`\n")

        if focus_property in ("acidity", "general"):
            parts.append("#### Acidity Analysis\n")
            total_shift = 0.0
            for s in subs:
                shift = self._estimate_pka_shift(7.0, s["substituent"], s["position"])  # Use dummy base
                total_shift += shift
                direction_text = "increases acidity (lowers pKa)" if shift < 0 else "decreases acidity (raises pKa)"
                parts.append(f"- **{s['substituent']}** ({s['position']}): estimated ΔpKa ≈ {shift:+.2f} → {direction_text}\n")
            
            if total_shift != 0:
                parts.append(f"> **Net inductive effect:** ΔpKa ≈ {total_shift:+.2f}\n")

        if focus_property in ("basicity", "general"):
            parts.append("\n#### Basicity Analysis\n")
            for s in subs:
                i_eff = self.i_effects.get(s["substituent"], 0)
                if isinstance(i_eff, tuple):
                    i_eff = i_eff[0]
                if i_eff > 0:
                    parts.append(f"- **{s['substituent']}**: Decreases basicity (EWG destabilizes conjugate acid)\n")
                elif i_eff < -0.05:
                    parts.append(f"- **{s['substituent']}**: Increases basicity slightly (EDG stabilizes cation)\n")

        if focus_property in ("reactivity", "general"):
            parts.append("\n#### Reactivity Analysis\n")
            for s in subs:
                i_eff = self.i_effects.get(s["substituent"], 0)
                if isinstance(i_eff, tuple):
                    i_eff = i_eff[0]
                if i_eff > 0.3:
                    parts.append(f"- **{s['substituent']}**: \n")
                    parts.append(f"  - **SN1/SN2:** Accelerates SN1 near the site (stabilizes carbocation); may slow SN2 (reduces nucleophilicity of adjacent atoms)\n")
                    parts.append(f"  - **EAS:** Deactivates ring; meta-directing (if aromatic)\n")
                    parts.append(f"  - **Elimination:** May accelerate E1/E2 by destabilizing adjacent C-H bonds\n")
                elif i_eff < -0.05:
                    parts.append(f"- **{s['substituent']}**: \n")
                    parts.append(f"  - **SN1:** Slightly destabilizes carbocation (weak EDG doesn't help much)\n")
                    parts.append(f"  - **EAS:** Activates ring; ortho/para-directing (if aromatic)\n")

        if focus_property in ("nmr_shift", "general"):
            parts.append("\n#### NMR Chemical Shift Prediction\n")
            for s in subs:
                pos = s["position"]
                if pos in ("alpha", "ortho"):
                    i_eff = self.i_effects.get(s["substituent"], 0)
                    if isinstance(i_eff, tuple):
                        i_eff = i_eff[0]
                    if i_eff > 0.2:
                        parts.append(f"- **{s['substituent']}** ({pos}): Deshields nearby protons by **+1.5 to +4.0 ppm** (downfield shift)\n")
                    elif i_eff < -0.05:
                        parts.append(f"- **{s['substituent']}** ({pos}): Slight shielding (**-0.2 to -0.5 ppm**) or minimal effect\n")
                elif pos in ("beta", "meta"):
                    parts.append(f"- **{s['substituent']}** ({pos}): Moderate deshielding **+0.3 to +1.0 ppm**\n")
                else:
                    parts.append(f"- **{s['substituent']}** ({pos}): Small or negligible effect on NMR shifts\n")

        # Dipole moment analysis
        if focus_property in ("dipole", "general"):
            parts.append("\n#### Dipole Moment Contribution\n")
            for s in subs:
                i_eff = self.i_effects.get(s["substituent"], 0)
                if isinstance(i_eff, tuple):
                    i_eff = i_eff[0]
                if abs(i_eff) > 0.2:
                    direction = "toward the substituent (δ+ on carbon side)" if i_eff > 0 else "away from the substituent"
                    magnitude = min(abs(i_eff) * 2.0, 3.0)  # Rough estimate in Debye
                    parts.append(f"- **{s['substituent']}**: Contributes ~{magnitude:.1f} D to dipole, pointing {direction}\n")

        # General rules summary
        parts.append("""
---

### 📐 Key Principles of Inductive Effects

| Feature | Description |
|---|---|
| **Transmission** | Through σ-bonds only; decays with distance (~1/3 per bond) |
| **Direction** | Determined by electronegativity difference |
| **Saturation** | Effect becomes negligible after 3-4 bonds |
| **Additivity** | Multiple substituents have approximately additive effects |
| **Field effect** | Long-range component operating through space/solvent |

> 💡 **Remember:** Inductive effects are permanent (don't depend on mechanism). They always operate regardless of the reaction type. Resonance effects can reinforce or oppose them depending on context.
""")

        return "\n".join(parts)

    def _run_text(self, input_params: str) -> str:
        parts = input_params.strip().split()
        if len(parts) < 1:
            raise ChemMCPError("Input must include molecule name. Format: 'molecule [property]'")
        molecule = parts[0]
        prop = parts[1] if len(parts) > 1 else "general"
        return self._run_base(molecule, prop)
