import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Common radionuclides and their primary decay modes
# Format: {symbol: (decay_mode, daughter, half_life_years)}
DECAY_DATA = {
    # Alpha emitters
    "U-238": ("α", "Th-234", 4.468e9),
    "U-235": ("α", "Th-231", 7.04e8),
    "U-234": ("α", "Th-230", 2.455e5),
    "Th-232": ("α", "Ra-228", 1.405e10),
    "Ra-226": ("α", "Rn-222", 1600),
    "Rn-222": ("α", "Po-218", 3.8235),
    "Po-210": ("α", "Pb-206", 0.3795),
    "Pu-239": ("α", "U-235", 24110),
    "Am-241": ("α", "Np-237", 432.2),
    # Beta-minus emitters
    "Sr-90": ("β⁻", "Y-90", 28.79),
    "Cs-137": ("β⁻", "Ba-137m", 30.08),
    "C-14": ("β⁻", "N-14", 5730),
    "H-3": ("β⁻", "He-3", 12.32),
    "P-32": ("β⁻", "S-32", 14.29),
    "Co-60": ("β⁻", "Ni-60*", 5.2714),
    "I-131": ("β⁻", "Xe-131", 8.02),
    "Y-90": ("β⁻", "Zr-90", 64.0 / (365.25)),
    "Kr-85": ("β⁻", "Rb-85", 10.72),
    "Tc-99": ("β⁻", "Ru-99", 211000),
    # Beta-plus / Positron emitters
    "Na-22": ("β⁺", "Ne-22", 2.6019),
    "F-18": ("β⁺", "O-18", 109.77 / 365.25 * 24 * 3600),  # ~110 min in years
    "C-11": ("β⁺", "B-11", 20.36 / 365.25 * 24 * 3600),
    # Gamma emitters (isomeric transition / gamma decay)
    "Ba-137m": ("γ(IT)", "Ba-137", 2.552 / (365.25 * 24 * 3600)),  # 2.55 min
    "Tc-99m": ("γ(IT)", "Tc-99", 6.0067 / (24)),  # ~6 hours in years
    "Co-60*": ("γ", "Ni-60", 0.0),  # excited state decays by gamma
}


@ChemMCPManager.register_tool
class RadioactiveDecay(BaseTool):
    """
    放射性衰变类型判断工具。
    根据核素名称或符号判断其主要衰变类型（α/β⁻/β⁺/γ）。
    """
    __version__      = "0.1.0"
    name             = "RadioactiveDecay"
    func_name        = "identify_decay_type"
    description      = "Identify the primary radioactive decay type (alpha, beta-minus, beta-plus, or gamma) for a given radionuclide."
    implementation_description = "Uses a built-in database of common radionuclides and their primary decay modes, half-lives, and daughter products."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Nuclear Chemistry", "Radioactive Decay", "Isotopes"]
    required_envs    = []

    code_input_sig   = [
        ("nuclide", "str", "N/A", "The nuclide identifier (e.g., 'U-238', 'C-14', 'Co-60')."),
    ]

    text_input_sig   = [
        ("nuclide_str", "str", "N/A", "The nuclide identifier as a string."),
    ]

    output_sig       = [
        ("nuclide", "str", "The input nuclide identifier."),
        ("decay_type", "str", "Primary decay mode: α (alpha), β⁻ (beta-minus), β⁺ (beta-plus/positron), γ (gamma/IT)."),
        ("daughter", "str", "The daughter product nuclide."),
        ("half_life", "str", "Half-life with appropriate units."),
        ("description", "str", "Human-readable description of the decay process."),
    ]

    examples         = [
        {
            "code_input": {"nuclide": "U-238"},
            "text_input": {"nuclide_str": "U-238"},
            "output": {
                "nuclide": "U-238",
                "decay_type": "α",
                "daughter": "Th-234",
                "half_life": "4.47×10⁹ years",
                "description": "Uranium-238 decays via alpha emission to Thorium-234.",
            }
        },
        {
            "code_input": {"nuclide": "C-14"},
            "text_input": {"nuclide_str": "C-14"},
            "output": {
                "nuclide": "C-14",
                "decay_type": "β⁻",
                "daughter": "N-14",
                "half_life": "5730 years",
                "description": "Carbon-14 decays via beta-minus emission to Nitrogen-14.",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _format_half_life(self, years: float) -> str:
        """Format half-life into human-readable string."""
        if years == 0:
            return "prompt (excited state)"
        if years < 1 / (365.25 * 24 * 3600):
            seconds = years * 365.25 * 24 * 3600
            return f"{seconds:.2g} s"
        elif years < 1 / (365.25 * 24):
            hours = years * 365.25 * 24
            return f"{hours:.2g} h"
        elif years < 1 / 12:
            days = years * 365.25
            return f"{days:.2g} d"
        elif years < 100:
            if years >= 1:
                return f"{years:.2g} years"
            months = years * 12
            return f"{months:.1f} months"
        elif years < 1e6:
            return f"{years:.3g} years"
        else:
            return f"{years:.3e} years"

    def _run_base(self, nuclide: str) -> dict:
        """
        Identify the decay type of a given radionuclide.
        """
        # Normalize input
        key = nuclide.strip()
        # Try exact match first
        if key not in DECAY_DATA:
            # Try case-insensitive
            for k in DECAY_DATA:
                if k.lower() == key.lower():
                    key = k
                    break
            else:
                raise ChemMCPError(
                    f"Nuclide '{nuclide}' not found in database. "
                    f"Available nuclides: {', '.join(sorted(DECAY_DATA.keys()))}"
                )

        decay_type, daughter, half_life_years = DECAY_DATA[key]
        hl_str = self._format_half_life(half_life_years)

        # Generate description
        type_names = {
            "α": "alpha emission",
            "β⁻": "beta-minus emission",
            "β⁺": "beta-plus (positron) emission",
            "γ(IT)": "gamma decay (isomeric transition)",
            "γ": "gamma emission",
        }
        desc = (
            f"{key} decays via {type_names.get(decay_type, decay_type)} "
            f"to {daughter}. Half-life: {hl_str}."
        )

        logger.info(f"Decay identification for {key}: {decay_type} → {daughter}")

        return {
            "nuclide": key,
            "decay_type": decay_type,
            "daughter": daughter,
            "half_life": hl_str,
            "description": desc,
        }

    def _run_text(self, nuclide_str: str) -> dict:
        return self._run_base(nuclide_str)
