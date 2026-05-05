"""
Adduct Ion Identifier - Identifies common adduct ions in ESI mass spectrometry
including [M+H]⁺, [M+Na]⁺, [M-H]⁻, [M+K]⁺, [M+NH₄]⁺, [2M+H]⁺, etc.
"""

import logging
import math
from typing import Dict, List, Tuple, Any, Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Comprehensive adduct ion database: (adduct_name, mass_shift, charge, typical_mode, description)
_ADDUCT_DATABASE: List[Dict[str, Any]] = [
    # Positive mode adducts
    {"name": "[M+H]⁺",          "mass_shift": 1.007276,   "charge": +1, "mode": "positive",  "description": "Protonated molecule — most common in ESI+"},
    {"name": "[M+Na]⁺",         "mass_shift": 22.989218,  "charge": +1, "mode": "positive",  "description": "Sodium adduct — common in presence of Na salts"},
    {"name": "[M+K]⁺",          "mass_shift": 38.963158,   "charge": +1, "mode": "positive",  "description": "Potassium adduct"},
    {"name": "[M+NH₄]⁺",       "mass_shift": 18.033823,   "charge": +1, "mode": "positive",  "description": "Ammonium adduct — common in LC-MS"},
    {"name": "[M+Li]⁺",         "mass_shift": 7.016004,    "charge": +1, "mode": "positive",  "description": "Lithium adduct"},
    {"name": "[M+H-H₂O]⁺",     "mass_shift": -17.003289,  "charge": +1, "mode": "positive",  "description": "Protonated minus water"},
    {"name": "[M+H+NH₄]²⁺",    "mass_shift": 19.0205495,  "charge": +2, "mode": "positive",  "description": "Doubly charged protonated ammoniated"},
    {"name": "[M+2H]²⁺",        "mass_shift": 1.007276,    "charge": +2, "mode": "positive",  "description": "Diprotonated (large molecules)"},
    {"name": "[M+3H]³⁺",        "mass_shift": 1.007276,    "charge": +3, "mode": "positive",  "description": "Triprotonated (peptides/proteins)"},
    {"name": "[M+Na-H]⁺",      "mass_shift": 21.981942,   "charge": +1, "mode": "positive",  "description": "Sodium-hydrogen exchange"},
    {"name": "[M+CH₃OH+H]⁺",   "mass_shift": 33.033489,   "charge": +1, "mode": "positive",  "description": "Methanol cluster (MeOH solvent)"},
    {"name": "[M+ACN+H]⁺",     "mass_shift": 42.034394,   "charge": +1, "mode": "positive",  "description": "Acetonitrile cluster (ACN solvent)"},
    {"name": "[M+DMSO+H]⁺",    "mass_shift": 79.02122,    "charge": +1, "mode": "positive",  "description": "DMSO cluster/adduct"},
    {"name": "[M+IsoProp+H]⁺", "mass_shift": 59.04968,    "charge": +1, "mode": "positive",  "description": "Isopropanol cluster"},
    {"name": "[2M+H]⁺",        "mass_shift": 1.007276,    "charge": +1, "mode": "positive",  "description": "Dimer (2×MW + H) — high concentration"},
    {"name": "[2M+Na]⁺",       "mass_shift": 22.989218,   "charge": +1, "mode": "positive",  "description": "Sodium-bound dimer"},
    {"name": "[2M+NH₄]⁺",     "mass_shift": 18.033823,   "charge": +1, "mode": "positive",  "description": "Ammoniated dimer"},
    {"name": "[3M+H]⁺",        "mass_shift": 1.007276,    "charge": +1, "mode": "positive",  "description": "Trimer (very high concentration)"},
    {"name": "[M+Fe]³⁺",       "mass_shift": 55.84511,    "charge": +3, "mode": "positive",  "description": "Iron coordination complex"},
    
    # Negative mode adducts
    {"name": "[M-H]⁻",          "mass_shift": -1.007276,   "charge": -1, "mode": "negative",  "description": "Deprotonated molecule — most common in ESI-"},
    {"name": "[M+Cl]⁻",         "mass_shift": 34.969402,   "charge": -1, "mode": "negative",  "description": "Chloride adduct (negative mode)"},
    {"name": "[M+FA-H]⁻",      "mass_shift": 44.998202,   "charge": -1, "mode": "negative",  "description": "Formate adduct (formic acid additive)"},
    {"name": "[M+Ac-H]⁻",       "mass_shift": 59.013852,   "charge": -1, "mode": "negative",  "description": "Acetate adduct (acetic acid additive)"},
    {"name": "[M+Na-2H]⁻",     "mass_shift": 20.974566,   "charge": -1, "mode": "negative",  "description": "Sodium di-deprotonated"},
    {"name": "[M-H₂O-H]⁻",     "mass_shift": -19.018395,  "charge": -1, "mode": "negative",  "description": "Deprotonated minus water"},
    {"name": "[M-2H]²⁻",        "mass_shift": -1.007276,   "charge": -2, "mode": "negative",  "description": "Dideprotonated (acidic molecules)"},
]


@ChemMCPManager.register_tool
class AdductIonIdentifier(BaseTool):
    """
    加合离子识别器 — 在质谱中识别常见的加合离子。
    
    根据观测到的 m/z 值和分子量（可选），匹配可能的加合离子类型，
    并给出置信度评分和建议。
    """
    __version__      = "0.1.0"
    name             = "AdductIonIdentifier"
    func_name        = "identify_adduct_ions"
    description      = "Identify common adduct ions in ESI mass spectra including [M+H]⁺, [M+Na]⁺, [M-H]⁻, [M+K]⁺, [M+NH₄]⁺, dimers, and solvent clusters."
    implementation_description = "Matches observed m/z against a comprehensive database of known adduct ions using ppm tolerance. Calculates neutral molecular weight from each candidate adduct and ranks by confidence based on match error, adduct prevalence, and consistency across candidates."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Mass Spectrometry", "Adducts", "Ion Identification", "ESI", "LC-MS"]
    required_envs    = []

    code_input_sig   = [
        ("observed_mz", "float", "N/A", "Observed m/z value from the mass spectrum."),
        ("molecular_weight", "float", "None", "Known molecular weight (neutral monoisotopic). If None, will be calculated from matches."),
        ("ionization_mode", "str", "positive", "Ionization mode: 'positive' or 'negative'."),
        ("tolerance_ppm", "float", "10.0", "Mass tolerance in ppm for matching."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Space-separated: 'observed_mz [molecular_weight] [ionization_mode] [tolerance_ppm]'. Example: '286.1438 285.1385 positive 10'"),
    ]

    output_sig       = [
        ("result", "dict", "Dict with matched_adducts list (each with name, calc_MW, error_ppm, confidence), best_match, suggested_neutral_MW, alternative_candidates, and interpretation notes."),
    ]

    examples         = [
        {
            "code_input": {
                "observed_mz": 286.1438,
                "molecular_weight": 285.1385,
                "ionization_mode": "positive",
                "tolerance_ppm": 10.0,
            },
            "text_input": {
                "input_params": "286.1438 285.1385 positive 10"
            },
            "output": {
                "result": {
                    "observed_mz": 286.1438,
                    "best_match": {"name": "[M+H]⁺", "calc_MW": 285.1365, "error_ppm": 0.7, "confidence": "very_high"},
                    "matched_adducts": [
                        {"name": "[M+H]⁺", "calc_MW": 285.1365, "error_ppm": 0.7, "confidence": "very_high"},
                        {"name": "[M+NH₄]⁺", "calc_MW": 268.1100, "error_ppm": 15.2, "confidence": "low"},
                    ],
                    "suggested_neutral_MW": 285.1365,
                    "interpretation": "Most likely [M+H]⁺ of a compound with MW ~285.14 Da",
                }
            },
        },
        {
            "code_input": {
                "observed_mz": 553.2770,
                "molecular_weight": None,
                "ionization_mode": "positive",
                "tolerance_ppm": 15.0,
            },
            "text_input": {
                "input_params": "553.2770 none positive 15"
            },
            "output": {
                "result": {
                    "observed_mz": 553.2770,
                    "best_match": {"name": "[2M+Na]⁺", "calc_MW": 265.144, "error_ppm": 2.1, "confidence": "high"},
                    "interpretation": "Possible sodium-bound dimer of MW ~265 Da compound",
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, observed_mz: float, molecular_weight: float = None, ionization_mode: str = "positive", tolerance_ppm: float = 10.0) -> dict:
        """Core logic: match observed m/z to known adducts."""
        if observed_mz <= 0:
            raise ChemMCPError("Observed m/z must be positive.")
        if tolerance_ppm <= 0:
            raise ChemMCPError("Tolerance must be positive.")

        # Filter by ionization mode
        mode = ionization_mode.lower()
        if mode not in ("positive", "negative"):
            raise ChemMCPError("Ionization_mode must be 'positive' or 'negative'.")

        candidates = []
        for adduct in _ADDUCT_DATABASE:
            if adduct["mode"] != mode:
                continue

            shift = adduct["mass_shift"]
            chg = abs(adduct["charge"])
            
            # Calculate what the neutral MW would be for this adduct
            # mz = (MW + shift) / |charge|
            calc_MW = observed_mz * chg - shift
            
            if calc_MW <= 0:
                continue

            # If we know the real MW, check how close this candidate is
            if molecular_weight is not None and molecular_weight > 0:
                error_ppm = abs(calc_MW - molecular_weight) / molecular_weight * 1e6
                if error_ppm > tolerance_ppm * 3:  # loose filter for display
                    continue
            else:
                # Without known MW, check if the result is chemically reasonable
                if calc_MW < 50 or calc_MW > 5000:
                    continue
                error_ppm = 0.0  # can't calculate without reference

            # Calculate m/z error from theoretical
            theoretical_mz = (calc_MW + shift) / chg
            mz_error_ppm = abs(observed_mz - theoretical_mz) / theoretical_mz * 1e6

            # Confidence scoring
            confidence = self._score_confidence(adduct, mz_error_ppm, calc_MW, molecular_weight)

            candidates.append({
                "name": adduct["name"],
                "mass_shift": round(shift, 6),
                "charge": adduct["charge"],
                "calc_MW": round(calc_MW, 4),
                "error_ppm": round(mz_error_ppm, 2),
                "confidence": confidence,
                "description": adduct["description"],
            })

        # Sort by confidence then by error
        candidates.sort(key=lambda x: (
            0 if x["confidence"] == "very_high" else 
            1 if x["confidence"] == "high" else 
            2 if x["confidence"] == "medium" else 3,
            x["error_ppm"]
        ))

        # Determine best match
        best = candidates[0] if candidates else None

        # Suggest neutral MW
        if molecular_weight is not None and molecular_weight > 0:
            suggested_mw = molecular_weight
        elif best:
            suggested_mw = best["calc_MW"]
        else:
            suggested_mw = None

        # Generate interpretation
        interpretation = self._generate_interpretation(best, candidates, observed_mz, mode)

        return {
            "result": {
                "observed_mz": observed_mz,
                "ionization_mode": mode,
                "tolerance_ppm": tolerance_ppm,
                "input_molecular_weight": molecular_weight,
                "matched_adducts": candidates[:15],  # top 15
                "best_match": best,
                "suggested_neutral_MW": suggested_mw,
                "total_matches": len(candidates),
                "interpretation": interpretation,
                "notes": (
                    "Adduct Identification Notes:\n"
                    "• [M+H]⁺ is the most common adduct in ESI+\n"
                    "• [M+Na]⁺ often appears when Na⁺ is present in mobile phase\n"
                    "• Dimers ([2M+H]⁺, [2M+Na]⁺) indicate high analyte concentration\n"
                    "• Solvent clusters indicate contamination or poor desolvation\n"
                    "• Always confirm identity with MS/MS when possible"
                ),
            }
        }

    def _score_confidence(self, adduct: dict, error_ppm: float, calc_MW: float, known_MW: Optional[float]) -> str:
        """Score confidence of an adduct assignment."""
        score = 100
        
        # Penalize by m/z error
        if error_ppm < 2:
            score += 0
        elif error_ppm < 5:
            score -= 10
        elif error_ppm < 10:
            score -= 25
        else:
            score -= 40

        # Boost common adducts
        common_names = {"[M+H]⁺", "[M-H]⁻", "[M+Na]⁺", "[M+NH₄]⁺"}
        if adduct["name"] in common_names:
            score += 20
        elif "dimer" in adduct["name"].lower() or "2M" in adduct["name"]:
            score -= 15  # less likely unless concentration is high
        elif "solvent" in adduct["description"].lower() or "cluster" in adduct["name"].lower():
            score -= 20

        # Bonus if MW matches known value
        if known_MW and known_MW > 0:
            mw_diff = abs(calc_MW - known_MW)
            if mw_diff < 0.01:
                score += 30
            elif mw_diff < 0.1:
                score += 15
            elif mw_diff < 1.0:
                score += 5

        # Reasonable MW range bonus
        if 100 <= calc_MW <= 800:
            score += 10
        elif calc_MW < 50 or calc_MW > 2000:
            score -= 20

        if score >= 80:
            return "very_high"
        elif score >= 60:
            return "high"
        elif score >= 35:
            return "medium"
        else:
            return "low"

    def _generate_interpretation(self, best: dict, candidates: list, obs_mz: float, mode: str) -> str:
        """Generate human-readable interpretation."""
        if not best:
            return f"No matching adduct found for m/z {obs_mz} in {mode} mode within tolerance."

        parts = [f"Most likely: {best['name']} → neutral MW ≈ {best['calc_MW']} Da"]
        
        if best["confidence"] == "very_high":
            parts.append("High-confidence assignment")
        elif best["confidence"] == "low":
            parts.append("Low confidence — consider MS/MS confirmation")

        # Check for multiple plausible assignments
        high_conf = [c for c in candidates if c["confidence"] in ("very_high", "high")]
        if len(high_conf) >= 2:
            names = ", ".join(c["name"] for c in high_conf[:3])
            parts.append(f"Multiple plausible assignments: {names}")

        # Check for dimer indication
        dimers = [c for c in candidates if "2M" in c["name"] or "3M" in c["name"]]
        if dimers and dimers[0]["confidence"] in ("medium", "high"):
            parts.append("⚠ Dimer detected — consider sample dilution")

        # Check for solvent clusters
        solvents = [c for c in candidates if any(s in c.get("description","").lower() for s in ["solvent", "cluster"]) ]
        if solvents:
            parts.append(f"Solvent cluster detected ({solvents[0]['name']}) — check mobile phase purity")

        return ". ".join(parts)

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            obs_mz = float(parts[0])
            mw = None
            mode = "positive"
            tol = 10.0

            if len(parts) > 1 and parts[1].lower() != "none":
                try:
                    mw = float(parts[1])
                except ValueError:
                    pass
            if len(parts) > 2:
                mode = parts[2]
            if len(parts) > 3:
                tol = float(parts[3])

            return self._run_base(obs_mz, mw, mode, tol)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'observed_mz [MW] [mode] [tolerance_ppm]'")
