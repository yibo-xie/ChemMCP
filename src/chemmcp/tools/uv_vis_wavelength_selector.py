"""
UV-Vis Wavelength Selector — UV-Vis最佳检测波长选择
基于发色团、共轭体系、溶剂效应推荐最佳检测波长
"""
import logging
from typing import Optional, List, Dict, Any

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见发色团 UV-Vis 数据 (λmax in nm, ε in M⁻¹cm⁻¹) ────────────
CHROMOPHORE_DATA: Dict[str, dict] = {
    # ── Simple chromophores ──
    "alkene_cc": {"name": "Isolated C=C", "lambda_max_nm": 170, "epsilon": 15000,
                  "solvent_effect": "minor", "notes": "π→π* transition; below typical instrument range."},
    "alkyne_cc": {"name": "C≡C (isolated)", "lambda_max_nm": 178, "epsilon": 10000,
                  "solvent_effect": "minor", "notes": "π→π*; often obscured by solvent cutoff."},
    "carbonyl_c_o": {"name": "C=O (ketone/aldehyde)", "lambda_max_nm": 280, "epsilon": 15,
                     "solvent_effect": "+5-10nm in polar solvent", "notes": "n→π*; weak (forbidden)."},
    "carboxyl": {"name": "COOH / COOR", "lambda_max_nm": 204, "epsilon": 60,
                 "solvent_effect": "moderate", "notes": "n→π*; very weak absorption."},
    "ester": {"name": "Ester C=O", "lambda_max_nm": 205, "epsilon": 60,
              "solvent_effect": "moderate", "notes": "n→π* similar to acid."},
    "amide": {"name": "Amide", "lambda_max_nm": 220, "epsilon": "~63 (weak)",
              "solvent_effect": "significant", "notes": "Weak n→π* ~220nm; stronger at ~190nm (π→π*)."},
    "nitro": {"name": "Nitro (-NO₂)", "lambda_max_nm": 270, "epsilon": 20,
              "solvent_effect": "moderate", "notes": "n→π* weak; also π→π* ~200nm (ε>5000)."},
    "nitrite": {"name": "Nitrite (-ONO)", "lambda_max_nm": 280, "epsilon": 40,
                "solvent_effect": "moderate", "notes": "n→π* and π→π* transitions."},
    "azo": {"name": "Azo (-N=N-)", "lambda_max_nm": 340, "epsilon": "~20 (n→π*)",
            "solvent_effect": "strong", "notes": "Intense color if conjugated; λmax shifts with substitution."},
    "imine_cn": {"name": "C=N (imine)", "lambda_max_nm": 235, "epsilon": 100,
                 "solvent_effect": "moderate", "notes": "n→π* weak."},

    # ── Conjugated systems ──
    "diene_isolated": {"name": "Diene (isolated)", "lambda_max_nm": "~175+175",
                       "epsilon": "~10000", "solvent_effect": "minor",
                       "notes": "Essentially same as isolated alkene ×2."},
    "diene_conjugated": {"name": "Conjugated diene", "lambda_max_nm": 217, "epsilon": 21000,
                         "solvent_effect": "minor", "notes": "Butadiene-type; strong π→π*."},
    "triene_conjugated": {"name": "Conjugated triene", "lambda_max_nm": 258, "epsilon": 35000,
                          "solvent_effect": "minor", "notes": "Woodward-Fieser applicable."},
    "benzene": {"name": "Benzene ring", "lambda_max_nm": 255, "epsilon": 215,
                "solvent_effect": "minor", "notes": "Fine structure: 255, 200, 180 nm bands. B-band 255nm (forbidden)."},
    "toluene": {"name": "Toluene (alkyl benzene)", "lambda_max_nm": 261, "epsilon": 300,
                "solvent_effect": "minor", "notes": "Alkyl substituent red-shifts ~5nm."},
    "phenol": {"name": "Phenol", "lambda_max_nm": 270, "epsilon": 1450,
               "solvent_effect": "strong (pH dependent)", "notes": "OH is auxochrome; pH shift: acidic → neutral λmax change."},
    "phenolate": {"name": "Phenolate anion", "lambda_max_nm": 287, "epsilon": 2600,
                  "solvent_effect": "N/A (aqueous base)", "notes": "Deprotonation greatly enhances intensity + red shift."},
    "aniline": {"name": "Aniline", "lambda_max_nm": 230, "epsilon": 8600,
                "solvent_effect": "strong (pH dependent)", "notes": "NH₂ auxochrome; protonation destroys conjugation."},
    "styrene": {"name": "Styrene (vinyl benzene)", "lambda_max_nm": 248, "epsilon": 14500,
                "solvent_effect": "minor", "notes": "Extended conjugation through vinyl group."},
    "naphthalene": {"name": "Naphthalene", "lambda_max_nm": 275, "epsilon": 5600,
                    "solvent_effect": "minor", "notes": "Fused rings; multiple bands: 275, 312 (forbidden), 220."},
    "anthracene": {"name": "Anthracene", "lambda_max_nm": 356, "epsilon": 7900,
                   "solvent_effect": "minor", "notes": "3 fused rings; visible fluorescence possible."},
    "biphenyl": {"name": "Biphenyl", "lambda_max_nm": 248, "epsilon": 17000,
                 "solvent_effect": "minor", "notes": "Two phenyl rings; planarity affects conjugation."},

    # ── Heteroaromatics ──
    "pyridine": {"name": "Pyridine", "lambda_max_nm": 256, "epsilon": 2750,
                 "solvent_effect": "moderate", "notes": "N-heteroaromatic; fine structure similar to benzene."},
    "pyrrole": {"name": "Pyrrole", "lambda_max_nm": 240, "epsilon": "unknown",
                "solvent_effect": "moderate", "notes": "Five-membered N-heterocycle."},
    "quinoline": {"name": "Quinoline", "lambda_max_nm": 275, "epsilon": 4500,
                  "solvent_effect": "minor", "notes": "Fused benzene-pyridine."},
    "isoquinoline": {"name": "Isoquinoline", "lambda_max_nm": 266, "epsilon": 6000,
                     "solvent_effect": "minor", "notes": "Isomer of quinoline."},

    # ── Common organic molecules with well-known spectra ──
    "acetone": {"name": "Acetone", "lambda_max_nm": 279, "epsilon": 15,
                "solvent_effect": "+5nm in water vs hexane", "notes": "n→π* of carbonyl; very weak."},
    "ethanol": {"name": "Ethanol", "lambda_max_nm": "<185", "epsilon": "N/A",
                "solvent_effect": "N/A", "notes": "No useful UV chromophore above 185nm; common solvent."},
    "methanol": {"name": "Methanol", "lambda_max_nm": "<177", "epsilon": "N/A",
                 "solvent_effect": "N/A", "notes": "UV-transparent down to ~205nm (HPLC grade); common solvent."},
    "acetonitrile": {"name": "Acetonitrile", "lambda_max_nm": "<190", "epsilon": "N/A",
                     "solvent_effect": "N/A", "notes": "HPLC-grade ACN transparent to ~190nm."},
    "hexane": {"name": "Hexane", "lambda_max_nm": "<195", "epsilon": "N/A",
               "solvent_effect": "N/A", "notes": "Excellent UV transparency for nonpolar analysis."},
    "water_hplc": {"name": "Water (HPLC grade)", "lambda_max_nm": "<191", "epsilon": "N/A",
                   "solvent_effect": "N/A", "notes": "Lowest UV cutoff of common solvents (~167nm theoretical)."},

    # ── Biomolecules ──
    "protein_280": {"name": "Proteins (Trp/Tyr)", "lambda_max_nm": 280, "epsilon": "variable",
                    "solvent_effect": "minor", "notes": "Trp ε~5500, Tyr ε~1400 M⁻¹cm⁻¹ at 280nm. Standard for protein quantification."},
    "dna_260": {"name": "DNA/RNA", "lambda_max_nm": 260, "epsilon": "variable",
                "solvent_effect": "minor", "notes": "Nucleic acid bases absorb strongly at 260nm. A260/280 ratio checks purity."},
    "nad_nadh": {"name": "NADH/NADPH", "lambda_max_nm": 340, "epsilon": 6220,
                 "solvent_effect": "minor", "notes": "Reduced form absorbs at 340nm; oxidized form does not. Key for enzyme assays."},
    "heme_soret": {"name": "Heme (Soret band)", "lambda_max_nm": 400, "epsilon": "~100000",
                   "solvent_effect": "moderate", "notes": "Extremely intense Soret band; also Q-bands at 500-600nm."},
    "chlorophyll_a": {"name": "Chlorophyll a", "lambda_max_nm": 663, "epsilon": "~90000",
                      "solvent_effect": "strong (solvent-dependent)", "notes": "Blue: ~430nm, Red: ~663nm (in ether). Green appearance."},
    "flavin_adenine": {"name": "FAD/FMN (flavins)", "lambda_max_nm": 450, "epsilon": 12200,
                       "solvent_effect": "minor", "notes": "Yellow color; 450nm (oxidized) → reduced loses this band."},
    "retinal": {"name": "Retinal (visual pigment)", "lambda_max_nm": 498, "epsilon": 42000,
                "solvent_effect": "strong", "notes": "Visible absorption; basis of vision. Protonated Schiff base shifts λmax."},
    "caffeine": {"name": "Caffeine", "lambda_max_nm": 273, "epsilon": 9700,
                 "solvent_effect": "minor", "notes": "Common HPLC-UV IS; good λmax for quantification."},
    "ibuprofen": {"name": "Ibuprofen", "lambda_max_nm": 222, "epsilon": 14500,
                  "solvent_effect": "minor", "notes": "Carboxylic acid NSAID; aromatic π→π*."},
    "paracetamol": {"name": "Paracetamol (Acetaminophen)", "lambda_max_nm": 245, "epsilon": "~9000",
                    "solvent_effect": "minor", "notes": "Phenolic acetanilide; good UV activity."},
}

# ── 溶剂紫外截止波长 ──────────────────────────────────────────────
SOLVENT_CUTOFF: Dict[str, float] = {
    "water": 191, "acetonitrile": 190, "methanol": 205, "ethanol": 210,
    "hexane": 195, "heptane": 200, "cyclohexane": 210, "dichloromethane": 235,
    "chloroform": 245, "thf": 230, "dioxane": 215, "dmf": 270, "dmso": 268,
    "acetone": 330, "ethyl_acetate": 260, "toluene": 286, "pyridine": 305,
    "isopropanol": 210, "butanol": 210, "pentane": 203,
}


@ChemMCPManager.register_tool
class UvVisWavelengthSelector(BaseTool):
    """
    UV-Vis 最佳波长选择器：根据化合物结构、发色团、溶剂等信息，
    推荐最佳检测波长，并给出光谱特征和注意事项。
    """
    __version__ = "0.1.0"
    name = "UvVisWavelengthSelector"
    func_name = "select_uvvis_wavelength"
    description = "Recommend optimal UV-Vis detection wavelength(s) based on compound structure, chromophores, conjugation, and solvent effects."
    implementation_description = "Uses a built-in database of chromophore λmax/ε values, Woodward-Fieser rules for dienes/enones, and solvent cutoff data to recommend optimal analytical wavelengths with sensitivity estimates."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["UV-Vis", "Wavelength Selection", "Chromophore", "Spectroscopy", "Analytical Chemistry", "Detection"]
    required_envs = []

    code_input_sig = [
        ("compound_name", "str", "", "Compound name or key to look up in database."),
        ("chromophores", "list", "[]", "List of chromophore keys (e.g., ['benzene', 'carbonyl_c_o'])."),
        ("conjugation_level", "int", "0", "Number of conjugated double bonds (for Woodward-Fieser estimation)."),
        ("auxochromes", "list", "[]", "List of auxochrome groups attached (e.g., ['-OH', '-NH₂', '-OCH₃'])."),
        ("solvent", "str", "", "Solvent used (e.g., 'methanol', 'water', 'hexane')."),
        ("detection_goal", "str", "sensitivity", "Goal: 'sensitivity' (max ε), 'selectivity', 'avoid_interference', or 'general'."),
        ("min_wavelength_nm", "float", "190", "Instrument minimum wavelength."),
        ("max_wavelength_nm", "str", "800", "Instrument maximum wavelength."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "E.g., 'benzene methanol' or 'diene_conjugated 3 -OH -OCH₃ ethanol'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with recommended wavelength(s), molar absorptivity, solvent compatibility, and practical notes."),
    ]

    examples = [
        {
            "code_input": {
                "compound_name": "benzene",
                "chromophores": [],
                "conjugation_level": 0,
                "auxochromes": [],
                "solvent": "methanol",
                "detection_goal": "sensitivity",
                "min_wavelength_nm": 190,
                "max_wavelength_nm": 800,
            },
            "text_input": {
                "input_params": "benzene methanol",
            },
            "output": {
                "result": {
                    "mode": "wavelength_selection",
                    "compound": "benzene",
                    "note": "Recommended wavelengths based on chromophore data.",
                }
            }
        },
        {
            "code_input": {
                "compound_name": "",
                "chromophores": ["diene_conjugated"],
                "conjugation_level": 3,
                "auxochromes": ["-OH"],
                "solvent": "ethanol",
                "detection_goal": "sensitivity",
                "min_wavelength_nm": 200,
                "max_wavelength_nm": 800,
            },
            "text_input": {
                "input_params": "diene_conjugated 3 -OH ethanol",
            },
            "output": {
                "result": {
                    "mode": "wavelength_selection",
                    "note": "Woodward-Fieser estimated wavelength.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _lookup_compound(self, name: str) -> Optional[dict]:
        key = name.lower().strip().replace(" ", "_")
        # Direct match
        if key in CHROMOPHORE_DATA:
            return CHROMOPHORE_DATA[key]
        # Fuzzy match
        for k, v in CHROMOPHORE_DATA.items():
            if key in k or k in key or name.lower() in v["name"].lower():
                return v
        return None

    @staticmethod
    def _woodward_fieser_diene(base: int, n_ext: int, n_ring: int = 0,
                               subs: list = None) -> int:
        """Woodward-Fieser rules for conjugated dienes (base values in nm)."""
        value = base  # parent acyclic=214, homoannular=253, heteroannular=214
        # Extending double bonds
        value += n_ext * 30
        # Ring residues
        ring_residue = min(n_ring, 1) * 39  # simplified
        value += ring_residue
        # Substituent increments
        sub_increments = {
            "-alkyl": 5, "-ring residue": 5, "-OR": 6, "-SR": 30,
            "-NR₂": 60, "-OCOR": 0, "-Cl, -Br": 5, "-R (allylic)": 5,
        }
        if subs:
            for s in subs:
                for pattern, inc in sub_increments.items():
                    if s.replace(" ", "").lower() in pattern.lower().replace(" ", "") \
                       or pattern.lower().replace(" ", "") in s.replace(" ", "").lower():
                        value += inc
                        break
        return value

    @staticmethod
    def _woodward_fieser_enone(base: int, alpha_subs: list = None,
                               beta_subs: list = None, gamma_delta_subs: list = None,
                               ex_cycles: list = None) -> int:
        """Woodward-Fieser rules for α,β-unsaturated ketones (base: acyclic/open=215, hexacyclic=214, pentacyclic=225)."""
        value = base
        alpha_inc = {"α-alkyl": 10, "α-OH, α-OR": 35, "α-OAc": 10}
        beta_inc = {"β-alkyl": 12, "β-OMe, β-OR": 35}
        gamma_delta_inc = {"γ, δ-alkyl (extended)": 30}
        ex_inc = {"exocyclic double bond": 5, "homodiene component": 39}

        for lst, inc_dict in [(alpha_subs or [], alpha_inc),
                               (beta_subs or [], beta_inc),
                               (gamma_delta_subs or [], gamma_delta_inc)]:
            for s in lst:
                for pat, val in inc_dict.items():
                    if pat.split(" ")[0].replace(",", "") in s.upper() or s.upper() in pat.upper():
                        value += val
                        break

        return value

    def _calc_auxochrome_shift(self, base_lambda: float, auxochromes: list) -> float:
        """Estimate bathochromic shift from auxochromes."""
        shift = 0
        for aux in auxochromes:
            a = aux.upper().strip()
            if a in ("-OH", "OH", "HYDROXY"):
                shift += 7  # phenolic OH: +5-7nm
            elif a in ("-NH2", "NH2", "AMINO"):
                shift += 20  # amino: +15-20nm
            elif a in ("-OCH3", "OCH3", "METHOXY"):
                shift += 7  # methoxy: +5-7nm
            elif a in ("-Cl", "-BR", "HALOGEN"):
                shift += 2  # halogen: small shift
            elif a in ("-NO2", "NITRO"):
                shift += 15  # nitro: significant shift
            else:
                shift += 5  # generic auxochrome
        return shift

    def _run_base(self, compound_name: str = "", chromophores: list = None,
                  conjugation_level: int = 0, auxochromes: list = None,
                  solvent: str = "", detection_goal: str = "sensitivity",
                  min_wavelength_nm: float = 190.0,
                  max_wavelength_nm: float = 800.0) -> dict:

        if chromophores is None:
            chromophores = []
        if auxochromes is None:
            auxochromes = []

        # Try direct lookup first
        looked_up = None
        if compound_name:
            looked_up = self._lookup_compound(compound_name)

        recommendations = []
        all_chromophores_used = []

        if looked_up:
            lam = looked_up["lambda_max_nm"]
            eps_str = str(looked_up["epsilon"])
            rec = {
                "source": f"Database lookup: {looked_up['name']}",
                "wavelength_nm": lam,
                "molar_absorptivity": eps_str,
                "transition_type": looked_up.get("notes", ""),
                "confidence": "high" if isinstance(lam, (int, float)) else "estimated",
            }
            recommendations.append(rec)
            all_chromophores_used.append(compound_name)

        # Process each specified chromophore
        for ch_key in chromophores:
            if ch_key in CHROMOPHORE_DATA:
                ch = CHROMOPHORE_DATA[ch_key]
                lam = ch["lambda_max_nm"]
                # Apply auxochrome shift
                if auxochromes:
                    shift = self._calc_auxochrome_shift(lam if isinstance(lam, (int, float)) else 250, auxochromes)
                    shifted_lam = lam + shift if isinstance(lam, (int, float)) else lam
                else:
                    shifted_lam = lam
                    shift = 0

                recommendations.append({
                    "source": f"Chromophore: {ch['name']}",
                    "wavelength_nm": shifted_lam,
                    "original_lambda_nm": lam,
                    "auxochrome_shift_nm": round(shift, 1) if shift else 0,
                    "molar_absorptivity": str(ch["epsilon"]),
                    "transition_type": ch.get("notes", ""),
                    "confidence": "high",
                })
                all_chromophores_used.append(ch_key)

        # Woodward-Fieser estimation for conjugated systems
        if conjugation_level > 0:
            if any(d in str(chromophores) for d in ["diene"]):
                wf_lambda = self._woodward_fieser_diene(
                    214, max(0, conjugation_level - 2), 0, auxochromes)
                recommendations.append({
                    "source": f"Woodward-Fieser (diene, n={conjugation_level})",
                    "wavelength_nm": wf_lambda,
                    "molar_absorptivity": "estimated 20000-50000",
                    "transition_type": "π→π* (allowed)",
                    "confidence": "estimated",
                })
            elif any(k in str(chromophores) for k in ["carbonyl", "enone"]):
                wf_lambda = self._woodward_fieser_enone(215, auxochromes, auxochromes)
                recommendations.append({
                    "source": f"Woodward-Fieser (enone, n={conjugation_level})",
                    "wavelength_nm": wf_lambda,
                    "molar_absorptivity": "estimated 10000-30000",
                    "transition_type": "π→π* + n→π*",
                    "confidence": "estimated",
                })

        # Solvent check
        solvent_info = {}
        if solvent:
            skey = solvent.lower().strip().replace(" ", "_")
            cutoff = SOLVENT_CUTOFF.get(skey) or SOLVENT_CUTOFF.get(solvent.lower())
            if cutoff:
                solvent_info = {"solvent": solvent, "uv_cutoff_nm": cutoff,
                                 "compatible": [r["wavelength_nm"] for r in recommendations
                                                if isinstance(r.get("wavelength_nm"), (int, float))
                                                and r["wavelength_nm"] > cutoff + 10]}

        # Filter by instrument range and select best
        valid_recs = [r for r in recommendations
                      if isinstance(r.get("wavelength_nm"), (int, float))
                      and min_wavelength_nm <= r["wavelength_nm"] <= max_wavelength_nm]

        best = valid_recs[0] if valid_recs else (recommendations[0] if recommendations else None)

        return {"result": {
            "mode": "wavelength_selection",
            "compound_name": compound_name or "not specified",
            "chromophores_specified": chromophores or ["auto-detected from name"],
            "auxochromes": auxochromes or [],
            "conjugation_level": conjugation_level,
            "solvent": solvent or "not specified",
            "detection_goal": detection_goal,
            "instrument_range_nm": [min_wavelength_nm, max_wavelength_nm],
            "all_candidates": recommendations,
            "best_recommendation": best,
            "solvent_compatibility": solvent_info,
            "practical_notes": self._generate_notes(best, solvent, detection_goal),
            "measurement_tips": [
                "Always run a blank (solvent-only) scan under identical conditions.",
                "For quantitative work, verify linearity at the chosen wavelength (Beer-Lambert range A=0.2-0.8).",
                "Consider scanning 200-800nm first to identify unexpected peaks or interferences.",
                "Temperature can affect λmax by ±1-3nm for some compounds.",
                "pH-sensitive compounds (phenols, anilines): buffer the mobile phase/solvent.",
            ],
        }}

    @staticmethod
    def _generate_notes(best: dict, solvent: str, goal: str) -> List[str]:
        notes = []
        if not best:
            notes.append("⚠ No specific recommendation available. Provide more structural information.")
            return notes

        wl = best.get("wavelength_nm", 0)
        eps = best.get("molar_absorptivity", "?")
        conf = best.get("confidence", "")

        if goal == "sensitivity":
            notes.append(f"🎯 For maximum sensitivity: use λ = {wl}nm where ε = {eps}.")
        elif goal == "selectivity":
            notes.append(f"🎯 For selectivity: consider secondary peaks or derivative spectroscopy.")
        elif goal == "avoid_interference":
            notes.append(f"🎯 To avoid interference: choose a less common wavelength region if matrix absorbs near {wl}nm.")

        if conf == "estimated":
            notes.append("⚠ This wavelength is an ESTIMATE from empirical rules. Experimental verification recommended.")

        if solvent:
            notes.append(f"Verify that your chosen wavelength is above the solvent cutoff for {solvent}.")

        if isinstance(wl, (int, float)) and wl < 210:
            notes.append("⚠ Low-wavelength region (<210nm): ensure high-purity solvents and clean optics.")
        elif isinstance(wl, (int, float)) and wl > 350:
            notes.append("Visible-range detection: colored compounds may be visible to the eye.")

        return notes

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            compound = parts[0] if parts else ""
            solvent = parts[1] if len(parts) > 1 else ""
            conj = int(parts[2]) if len(parts) > 2 else 0
            aux = parts[3:] if len(parts) > 3 else []
            return self._run_base(compound_name=compound, solvent=solvent,
                                   conjugation_level=conj, auxochromes=aux)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input '{input_params}': {e}")
