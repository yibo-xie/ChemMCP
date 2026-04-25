import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class AmphotericSpecies(BaseTool):
    """
    两性物质pH计算。
    支持两性阴离子（如HCO3-、HPO4^2-）、两性氢氧化物（如Al(OH)3、Zn(OH)2）、氨基酸等。
    """
    __version__ = "0.1.0"
    name = "AmphotericSpecies"
    func_name = "calculate_amphoteric_ph"
    description = "Calculate the pH of amphoteric species solutions (amphoteric anions, amphoterichydroxides, amino acids)."
    implementation_description = "Uses appropriate approximation formulas for different types of amphoteric substances: [H+]≈√(Ka1·Ka2) for acid salts like HCO3-, pI=(pKa1+pKa2)/2 for amino acids, equilibrium analysis for amphoteric hydroxides."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Amphoteric", "pH", "Equilibrium", "Acid-Base"]
    required_envs = []

    code_input_sig = [
        ("substance", "str", "N/A", "Chemical formula or name of the amphoteric substance, e.g., 'HCO3-', 'Al(OH)3', 'glycine'."),
        ("concentration", "float", "0.1", "Concentration in mol/L."),
        ("ka", "float", "None", "Optional: custom Ka value. If provided, overrides built-in data."),
        ("kb", "float", "None", "Optional: custom Kb value. If provided, overrides built-in data."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'substance concentration [ka] [kb]'. Example: 'HCO3- 0.1' or 'glycine 0.05'."),
    ]

    output_sig = [
        ("ph", "float", "Calculated pH value."),
        ("h_conc", "float", "[H+] concentration (mol/L)."),
        ("dominant_species", "str", "The dominant species present at this pH."),
        ("substance_type", "str", "Type classification: 'acid_salt', 'amino_acid', 'amphoteric_hydroxide', or 'water'."),
        ("explanation", "str", "Detailed explanation of the calculation method and result."),
    ]

    examples = [
        {
            "code_input": {
                "substance": "HCO3-",
                "concentration": 0.1,
                "ka": None,
                "kb": None,
            },
            "text_input": {
                "input_params": "HCO3- 0.1",
            },
            "output": {
                "ph": 8.31,
                "h_conc": 4.9e-9,
                "dominant_species": "HCO3- and H2CO3/CO3^2- mixture",
                "substance_type": "acid_salt",
                "explanation": "For NaHCO3 solution: [H+] ≈ √(Ka1·Ka2) = √(4.45×10⁻⁷ × 4.69×10⁻¹¹) ≈ 4.57×10⁻⁹ M, pH ≈ 8.34.",
            },
        },
        {
            "code_input": {
                "substance": "glycine",
                "concentration": 0.05,
                "ka": None,
                "kb": None,
            },
            "text_input": {
                "input_params": "glycine 0.05",
            },
            "output": {
                "ph": 6.06,
                "h_conc": 8.7e-7,
                "dominant_species": "zwitterion H3N+CH2COO-",
                "substance_type": "amino_acid",
                "explanation": "For glycine (isoelectric point): pH = pI = (pKa1 + pKa2)/2 = (2.34 + 9.60)/2 = 5.97.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Initialize database of amphoteric substances with Ka/Kb values."""
        # Format: {name: {"type": ..., "Ka": [...], "Kb": [...], "pKa": [...], "pKb": [...]}}
        # For acid salts (like HCO3-): Ka = second dissociation constant, Kb = Kw/first_Ka
        # For amphoteric hydroxides: Ka (as acid), Kb (as base)
        self._db = {
            # --- Acid salts (anions of polyprotic acids) ---
            "HCO3-": {
                "type": "acid_salt",
                "Ka": 4.69e-11,      # Ka2 of H2CO3
                "Ka_prev": 4.45e-7,   # Ka1 of H2CO3
                "pKa": 10.33,
                "pKa_prev": 6.35,
            },
            "HS-": {
                "type": "acid_salt",
                "Ka": 1.0e-19,        # Ka2 of H2S (approximate)
                "Ka_prev": 9.5e-8,     # Ka1 of H2S
                "pKa": 19.0,
                "pKa_prev": 7.02,
            },
            "HPO4^2-": {
                "type": "acid_salt",
                "Ka": 4.2e-13,         # Ka3 of H3PO4
                "Ka_prev": 6.23e-8,    # Ka2 of H3PO4
                "pKa": 12.38,
                "pKa_prev": 7.21,
            },
            "H2PO4-": {
                "type": "acid_salt",
                "Ka": 6.23e-8,         # Ka2 of H3PO4
                "Ka_prev": 7.52e-3,    # Ka1 of H3PO4
                "pKa": 7.21,
                "pKa_prev": 2.12,
            },

            # --- Amphoteric hydroxides ---
            "Al(OH)3": {
                "type": "amphoteric_hydroxide",
                "Ka": 1.0e-11,          # As acid: Al(OH)3 ⇌ Al(OH)2+ + OH-
                "Kb": 5.0e-5,           # As base: Al(OH)3 + 3H+ ⇌ Al^3+ + 3H2O (effective)
                "pKa": 11.0,
                "pKb": 4.3,
            },
            "Zn(OH)2": {
                "type": "amphoteric_hydroxide",
                "Ka": 1.0e-12,
                "Kb": 1.0e-4,
                "pKa": 12.0,
                "pKb": 4.0,
            },
            "Pb(OH)2": {
                "type": "amphoteric_hydroxide",
                "Ka": 1.0e-14,
                "Kb": 1.0e-3,
                "pKa": 14.0,
                "pKb": 3.0,
            },
            "Cr(OH)3": {
                "type": "amphoteric_hydroxide",
                "Ka": 1.0e-12,
                "Kb": 1.0e-5,
                "pKa": 12.0,
                "pKb": 5.0,
            },
            "Sn(OH)2": {
                "type": "amphoteric_hydroxide",
                "Ka": 1.0e-13,
                "Kb": 1.0e-4,
                "pKa": 13.0,
                "pKb": 4.0,
            },

            # --- Amino acids ---
            "glycine": {
                "type": "amino_acid",
                "pKa1": 2.34,   # -COOH group
                "pKa2": 9.60,   # -NH3+ group
                "pI": 5.97,
            },
            "alanine": {
                "type": "amino_acid",
                "pKa1": 2.34,
                "pKa2": 9.69,
                "pI": 6.01,
            },
            "valine": {
                "type": "amino_acid",
                "pKa1": 2.32,
                "pKa2": 9.62,
                "pI": 5.97,
            },
            "aspartic_acid": {
                "type": "amino_acid",
                "pKa1": 1.99,
                "pKa2": 3.90,
                "pKa3": 9.90,
                "pI": 2.95,
            },
            "lysine": {
                "type": "amino_acid",
                "pKa1": 2.18,
                "pKa2": 8.95,
                "pKa3": 10.53,
                "pI": 9.74,
            },

            # --- Water (self-ionization) ---
            "H2O": {
                "type": "water",
                "Kw": 1.0e-14,
            },
        }

        self._aliases = {
            "bicarbonate": "HCO3-", "sodium bicarbonate": "HCO3-", "nahco3": "HCO3-",
            "hydrogen sulfide ion": "HS-", "bisulfide": "HS-",
            "hydrogen phosphate": "HPO4^2-", "hpo4(2-)": "HPO4^2-",
            "dihydrogen phosphate": "H2PO4-", "h2po4-": "H2PO4-",
            "aluminum hydroxide": "Al(OH)3", "al(oh)3": "Al(OH)3",
            "zinc hydroxide": "Zn(OH)2", "zn(oh)2": "Zn(OH)2",
            "lead hydroxide": "Pb(OH)2", "pb(oh)2": "Pb(OH)2",
            "chromium hydroxide": "Cr(OH)3", "cr(oh)3": "Cr(OH)3",
            "tin(ii) hydroxide": "Sn(OH)2", "sn(oh)2": "Sn(OH)2",
            "water": "H2O",
        }

    def _run_base(self, substance: str, concentration: float = 0.1,
                  ka: Optional[float] = None, kb: Optional[float] = None) -> dict:
        """Core logic: calculate pH of amphoteric substance."""
        if concentration <= 0:
            raise ChemMCPError("Concentration must be positive.")

        kw = 1.0e-14
        name_key = self._resolve_name(substance)

        if name_key not in self._db:
            raise ChemMCPError(
                f"Unknown amphoteric substance: '{substance}'. "
                f"Available: {list(self._db.keys())}"
            )

        data = self._db[name_key]
        stype = data["type"]

        if stype == "acid_salt":
            ph_result = self._calc_acid_salt(data, concentration, ka, kb, kw)
        elif stype == "amino_acid":
            ph_result = self._calc_amino_acid(data, concentration, ka, kb, kw)
        elif stype == "amphoteric_hydroxide":
            ph_result = self._calc_ampho_hydroxide(data, concentration, ka, kb, kw)
        elif stype == "water":
            ph_result = {"ph": 7.0, "h_conc": 1e-7}
        else:
            raise ChemMCPError(f"Unknown substance type: {stype}")

        ph_result["substance_type"] = stype
        logger.info(f"AmphotericSpecies: {substance} C={concentration} → pH={ph_result['ph']:.2f}")
        return ph_result

    def _calc_acid_salt(self, data: dict, conc: float,
                        ka_override=None, kb_override=None, kw=1e-14) -> dict:
        """Calculate pH for acid salt solutions like HCO3-, H2PO4-."""
        Ka = ka_override if ka_override is not None else data.get("Ka", 0)
        Ka_prev = data.get("Ka_prev", 0)

        # Approximation: [H+] ≈ sqrt(Ka_prev * Ka)
        h_conc = math.sqrt(Ka_prev * Ka)
        ph = -math.log10(h_conc) if h_conc > 0 else 7.0

        dominant = f"{data.get('name', 'Amphoteric anion')} with small amounts of protonated/deprotonated forms"

        explanation = (
            f"For {self._current_substance}: [H⁺] ≈ √(Ka₁·Ka₂) "
            f"= √({Ka_prev:.2e} × {Ka:.2e}) = {h_conc:.2e} M\n"
            f"pH = {ph:.2f}. This is the approximate pH of a solution of the "
            f"amphoteric intermediate species."
        )

        return {
            "ph": round(ph, 2),
            "h_conc": round(h_conc, 15),
            "dominant_species": dominant,
            "explanation": explanation,
        }

    def _calc_amino_acid(self, data: dict, conc: float,
                         ka_override=None, kb_override=None, kw=1e-14) -> dict:
        """Calculate pH for amino acids using isoelectric point."""
        pKa1 = data.get("pKa1", 2.34)
        pKa2 = data.get("pKa2", 9.60)

        if "pI" in data:
            pI = data["pI"]
        else:
            pI = (pKa1 + pKa2) / 2.0

        h_conc = 10 ** (-pI)
        ph = pI

        explanation = (
            f"For {self._current_substance}: isoelectric point pI = (pKa₁ + pKa₂)/2 "
            f"= ({pKa1} + {pKa2}) / 2 = {pI:.2f}\n"
            f"At pI, the zwitterion form dominates. pH = {ph:.2f}, [H⁺] = {h_conc:.2e} M."
        )

        return {
            "ph": round(ph, 2),
            "h_conc": round(h_conc, 20),
            "dominant_species": "zwitterion (net neutral)",
            "explanation": explanation,
        }

    def _calc_ampho_hydroxide(self, data: dict, conc: float,
                              ka_override=None, kb_override=None, kw=1e-14) -> dict:
        """Calculate pH for amphoteric hydroxides (rough estimate)."""
        Ka = ka_override if ka_override is not None else data.get("Ka", 1e-12)
        Kb = kb_override if kb_override is not None else data.get("Kb", 1e-5)

        # Rough estimate: pH where both equilibria are considered
        # As base: MOH + H2O ⇌ [M(OH)2]+ + OH-
        # As acid: MOH ⇌ [MO]- + H+
        # Approximate pH from balance
        if Kb > Ka:
            # More basic than acidic
            oh_conc = math.sqrt(Kb * conc)
            if oh_conc > 1e-7:
                poh = -math.log10(oh_conc)
                ph = 14 - poh
            else:
                ph = 7.0
            h_conc = 10 ** (-ph)
            dominant = "slightly basic, undissolved solid dominates"
        else:
            # More acidic than basic
            h_conc = math.sqrt(Ka * conc)
            ph = -math.log10(h_conc) if h_conc > 0 else 7.0
            dominant = "slightly acidic, undissolved solid dominates"

        explanation = (
            f"For {self._current_substance}: Ka={Ka:.2e}, Kb={Kb:.2e}\n"
            f"The substance shows amphoteric behavior. At C={conc} M, "
            f"estimated pH ≈ {ph:.2f}. Note: actual pH depends on "
            f"suspension/solubility details."
        )

        return {
            "ph": round(ph, 2),
            "h_conc": round(h_conc, 20),
            "dominant_species": dominant,
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if len(parts) < 2:
                raise ValueError("Need substance and concentration. Format: 'substance concentration [ka] [kb]'")
            substance = parts[0]
            conc = float(parts[1])
            ka = float(parts[2]) if len(parts) > 2 else None
            kb = float(parts[3]) if len(parts) > 3 else None
            return self._run_base(substance, conc, ka, kb)
        except ValueError as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'substance concentration [ka] [kb]'")

    def _resolve_name(self, name: str) -> str:
        """Resolve substance name/alias to canonical key."""
        self._current_substance = name
        name_lower = name.strip()
        if name_lower in self._db:
            return name_lower
        alias_lower = name_lower.lower()
        if alias_lower in self._aliases:
            return self._aliases[alias_lower]
        # Try case-insensitive match on db keys
        for k in self._db:
            if k.lower() == alias_lower:
                return k
        return name_lower
