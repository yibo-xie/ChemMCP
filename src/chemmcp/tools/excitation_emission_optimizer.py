"""
Excitation/Emission Wavelength Optimizer — 激发/发射波长优化选择工具 (#323)

功能：
  根据荧光物质的光谱数据和仪器约束，选择最优激发/发射波长组合。
  内置常见荧光物质参考数据库。
"""

import logging
import math
from typing import Optional, List, Dict, Any, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# ── 常见荧光物质光谱数据库 (λ_ex_max, λ_em_max, Stokes shift) ──
FLUOROPHORE_DB: Dict[str, Dict[str, Any]] = {
    "fluorescein": {
        "ex_max": 494, "ex_range": [440, 520],
        "em_max": 514, "em_range": [500, 650],
        "quantum_yield": 0.93,
        "extinction_coefficient": 93000,  # M⁻¹cm⁻¹ at λ_ex_max
        "solvent": "0.1M NaOH (aq)",
        "notes": "High quantum yield; pH sensitive",
    },
    "rhodamine_b": {
        "ex_max": 554, "ex_range": [520, 580],
        "em_max": 577, "em_range": [570, 650],
        "quantum_yield": 0.49,
        "extinction_coefficient": 116000,
        "solvent": "ethanol",
        "notes": "Photostable; widely used",
    },
    "rhodamine_6g": {
        "ex_max": 488, "ex_range": [450, 520],
        "em_max": 550, "em_range": [530, 620],
        "quantum_yield": 0.95,
        "extinction_coefficient": 116000,
        "solvent": "ethanol",
        "notes": "Very high QY; laser dye",
    },
    "quinine_sulfate": {
        "ex_max": 348, "ex_range": [320, 380],
        "em_max": 450, "em_range": [400, 550],
        "quantum_yield": 0.54,
        "extinction_coefficient": 10000,
        "solvent": "0.5M H2SO4",
        "notes": "Quantum yield standard",
    },
    "gfp": {   # Green Fluorescent Protein
        "ex_max": 395, "ex_range": [370, 420],
        "em_max": 509, "em_range": [480, 560],
        "quantum_yield": 0.79,
        "extinction_coefficient": 83000,
        "solvent": "aqueous buffer pH 7-8",
        "notes": "Also excitable at 475 nm (minor band)",
    },
    "egfp": {  # Enhanced GFP
        "ex_max": 488, "ex_range": [470, 510],
        "em_max": 509, "em_range": [480, 560],
        "quantum_yield": 0.60,
        "extinction_coefficient": 55000,
        "solvent": "aqueous buffer pH 7-8",
        "notes": "Optimized for 488 nm argon laser",
    },
    "dapi": {
        "ex_max": 358, "ex_range": [340, 380],
        "em_max": 461, "em_range": [430, 520],
        "quantum_yield": 0.46,
        "extinction_coefficient": 35000,
        "solvent": "buffer with DNA",
        "notes": "DNA stain; fluorescence increases ~20x when bound to DNA",
    },
    "cy5": {
        "ex_max": 649, "ex_range": [620, 670],
        "em_max": 670, "em_range": [655, 720],
        "quantum_yield": 0.28,
        "extinction_coefficient": 250000,
        "solvent": "aqueous/pH neutral",
        "notes": "Far-red; low background",
    },
    "alexa_fluor_488": {
        "ex_max": 495, "ex_range": [465, 525],
        "em_max": 519, "em_range": [505, 600],
        "quantum_yield": 0.92,
        "extinction_coefficient": 73000,
        "solvent": "aqueous buffer",
        "notes": "Photostable alternative to fluorescein",
    },
    "tryptophan": {
        "ex_max": 280, "ex_range": [260, 300],
        "em_max": 348, "em_range": [320, 400],
        "quantum_yield": 0.13,
        "extinction_coefficient": 5600,
        "solvent": "water pH 7",
        "notes": "Intrinsic protein fluorophore",
    },
    "anthracene": {
        "ex_max": 356, "ex_range": [340, 380],
        "em_max": 402, "em_range": [380, 450],
        "quantum_yield": 0.27,
        "extinction_coefficient": 8000,
        "solvent": "ethanol",
        "notes": "PAH standard",
    },
}


@ChemMCPManager.register_tool
class ExcitationEmissionOptimizer(BaseTool):
    """
    激发/发射波长优化选择工具。
    根据荧光物质光谱数据、仪器约束和Stokes位移要求，推荐最优波长。
    """
    __version__                = "0.1.0"
    name                       = "ExcitationEmissionOptimizer"
    func_name                  = "optimize_wavelengths"
    description                = ("Optimize selection of excitation and emission wavelengths "
                                 "for fluorescence measurements based on fluorophore properties and instrument constraints.")
    implementation_description = ("Uses built-in spectral database of common fluorophores or accepts user-provided "
                                 "spectral data to find optimal (λ_ex, λ_em) pairs maximizing signal-to-noise ratio.")
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Fluorescence", "Spectroscopy", "Wavelength Selection",
                                   "Analytical Chemistry", "Optimization"]
    required_envs              = []

    code_input_sig = [
        ("compound_name",               "str",   "N/A",         "Name of the fluorescent compound (e.g., 'fluorescein', 'rhodamine_b')."),
        ("instrument_min_ex_nm",        "float", "200.0",       "Instrument minimum excitation wavelength (nm)."),
        ("instrument_max_ex_nm",        "float", "800.0",       "Instrument maximum excitation wavelength (nm)."),
        ("instrument_min_em_nm",        "float", "200.0",       "Instrument minimum emission wavelength (nm)."),
        ("instrument_max_em_nm",        "float", "900.0",       "Instrument maximum emission wavelength (nm)."),
        ("stokes_shift_min_nm",         "float", "20.0",        "Minimum acceptable Stokes shift (nm) to avoid scatter interference."),
        ("priority",                    "str",   "balanced",     "Optimization priority: 'sensitivity', 'stability', or 'balanced'."),
    ]

    text_input_sig = [
        ("input_params",                "str",   "N/A",
         "Space-separated: compound_name [min_ex max_ex min_em max_em stokes_min priority]"),
    ]

    output_sig = [
        ("optimal_ex_wavelength",       "float", "Recommended excitation wavelength (nm)."),
        ("optimal_em_wavelength",       "float", "Recommended emission wavelength (nm)."),
        ("stokes_shift_nm",             "float", "Stokes shift between optimal ex/em (nm)."),
        ("fluorophore_info",            "dict",  "Database info for the selected compound."),
        ("recommendations",             "str",   "Practical recommendations for measurement setup."),
        ("warnings",                    "list",  "Any warnings about constraints or trade-offs."),
    ]

    examples = [
        {
            "code_input": {
                "compound_name": "fluorescein",
            },
            "text_input": {"input_params": "fluorescein"},
            "output": {
                "optimal_ex_wavelength": 494.0,
                "optimal_em_wavelength": 514.0,
                "stokes_shift_nm": 20.0,
                "recommendations": "",
            },
        },
        {
            "code_input": {
                "compound_name": "quinine_sulfate",
                "stokes_shift_min_nm": 80.0,
            },
            "text_input": {"input_params": "quinine_sulfate 200 800 200 900 80 balanced"},
            "output": {
                "optimal_ex_wavelength": 348.0,
                "optimal_em_wavelength": 450.0,
                "stokes_shift_nm": 102.0,
                "recommendations": "",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Load fluorophore database."""
        self.db = FLUOROPHORE_DB

    def _find_compound(self, name: str) -> Dict[str, Any]:
        """Case-insensitive lookup."""
        key = name.strip().lower().replace(" ", "_")
        if key in self.db:
            return self.db[key]
        # fuzzy match
        matches = [k for k in self.db if key in k.lower()]
        if len(matches) == 1:
            return self.db[matches[0]]
        elif matches:
            raise ChemMCPError(f"Ambiguous compound '{name}'. Did you mean: {matches}?")
        available = ", ".join(sorted(self.db.keys()))
        raise ChemMCPError(f"Compound '{name}' not found. Available: {available}")

    def _run_base(
        self,
        compound_name: str,
        instrument_min_ex_nm: float = 200.0,
        instrument_max_ex_nm: float = 800.0,
        instrument_min_em_nm: float = 200.0,
        instrument_max_em_nm: float = 900.0,
        stokes_shift_min_nm: float = 20.0,
        priority: str = "balanced",
    ) -> Dict[str, Any]:
        """Core optimization logic."""
        # Look up compound
        info = self._find_compound(compound_name)
        ex_max = info["ex_max"]
        em_max = info["em_max"]
        ex_range = info.get("ex_range", [ex_max - 20, ex_max + 20])
        em_range = info.get("em_range", [em_max - 20, em_max + 20])

        warnings = []

        # Check instrument compatibility
        if ex_max < instrument_min_ex_nm:
            warnings.append(f"Excitation peak {ex_max} nm is below instrument minimum ({instrument_min_ex_nm} nm).")
        if ex_max > instrument_max_ex_nm:
            warnings.append(f"Excitation peak {ex_max} nm exceeds instrument maximum ({instrument_max_ex_nm} nm).")
        if em_max < instrument_min_em_nm:
            warnings.append(f"Emission peak {em_max} nm is below instrument minimum ({instrument_min_em_nm} nm).")
        if em_max > instrument_max_em_nm:
            warnings.append(f"Emission peak {em_max} nm exceeds instrument maximum ({instrument_max_em_nm} nm).")

        # Clamp to instrument range
        opt_ex = max(instrument_min_ex_nm, min(instrument_max_ex_nm, ex_max))
        opt_em = max(instrument_min_em_nm, min(instrument_max_em_nm, em_max))

        # Adjust for minimum Stokes shift
        stokes = opt_em - opt_ex
        if stokes < stokes_shift_min_nm:
            # Shift emission higher to meet Stokes requirement
            opt_em = opt_ex + stokes_shift_min_nm
            if opt_em > instrument_max_em_nm:
                # Try shifting excitation lower instead
                opt_em = em_max
                opt_ex = opt_em - stokes_shift_min_nm
                if opt_ex < instrument_min_ex_nm:
                    warnings.append(
                        f"Cannot achieve minimum Stokes shift of {stokes_shift_min_nm} nm within instrument range."
                    )
            stokes = opt_em - opt_ex

        # Priority-based fine-tuning
        rec_parts = []
        if priority == "sensitivity":
            rec_parts.append("Use λ_ex at absorption maximum for highest sensitivity.")
        elif priority == "stability":
            offset = min(5.0, (ex_range[1] - ex_max) / 2)
            opt_ex = ex_max + offset  # Slightly off-peak reduces photobleaching
            rec_parts.append("Excite slightly off-peak to reduce photobleaching.")
        else:  # balanced
            rec_parts.append("Use peak wavelengths for balanced performance.")

        # General recommendations
        qy = info.get("quantum_yield", "N/A")
        if isinstance(qy, float) and qy > 0.8:
            rec_parts.append(f"High quantum yield ({qy}) — excellent for trace analysis.")
        elif isinstance(qy, float) and qy < 0.2:
            rec_parts.append(f"Low quantum yield ({qy}) — consider signal averaging or more sensitive detector.")

        notes = info.get("notes", "")
        if notes:
            rec_parts.append(f"Note: {notes}")

        logger.info(f"Optimized for '{compound_name}': λ_ex={opt_ex:.1f} nm, λ_em={opt_em:.1f} nm, Δλ={stokes:.1f} nm")
        return {
            "optimal_ex_wavelength": round(opt_ex, 2),
            "optimal_em_wavelength": round(opt_em, 2),
            "stokes_shift_nm": round(stokes, 2),
            "fluorophore_info": info,
            "recommendations": " ".join(rec_parts),
            "warnings": warnings,
        }

    def _run_text(self, input_params: str) -> Dict[str, Any]:
        """Parse text input."""
        try:
            parts = input_params.split()
            if not parts:
                raise ValueError("Empty input.")

            compound_name = parts[0].replace("_", " ")
            kwargs = {}
            if len(parts) > 1:
                kwargs["instrument_min_ex_nm"] = float(parts[1])
            if len(parts) > 2:
                kwargs["instrument_max_ex_nm"] = float(parts[2])
            if len(parts) > 3:
                kwargs["instrument_min_em_nm"] = float(parts[3])
            if len(parts) > 4:
                kwargs["instrument_max_em_nm"] = float(parts[4])
            if len(parts) > 5:
                kwargs["stokes_shift_min_nm"] = float(parts[5])
            if len(parts) > 6:
                kwargs["priority"] = parts[6]

            return self._run_base(compound_name, **kwargs)
        except ChemMCPError:
            raise
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}")
