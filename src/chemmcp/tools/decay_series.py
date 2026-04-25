import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

# Major radioactive decay series data
# Each series: list of (nuclide, decay_mode, half_life_str)
DECAY_SERIES = {
    "uranium-238": {
        "name": "Uranium-238 Series (4n+2)",
        "parent": "U-238",
        "stable_end": "Pb-206",
        "chain": [
            ("U-238",    "α",   "4.47×10⁹ y"),
            ("Th-234",   "β⁻",  "24.1 d"),
            ("Pa-234m",  "β⁻/IT", "1.17 min"),
            ("Pa-234",   "β⁻",  "6.70 h"),
            ("U-234",    "α",   "2.45×10⁵ y"),
            ("Th-230",   "α",   "7.54×10⁴ y"),
            ("Ra-226",   "α",   1600),
            ("Rn-222",   "α",   "3.82 d"),
            ("Po-218",   "α/β⁻", "3.10 min"),
            ("Pb-214",   "β⁻",  "26.8 min"),
            ("Bi-214",   "α/β⁻", "19.9 min"),
            ("Po-214",   "α",   "164 μs"),
            ("Pb-210",   "β⁻",  "22.3 y"),
            ("Bi-210",   "β⁻",  "5.01 d"),
            ("Po-210",   "α",   "138 d"),
            ("Pb-206",   "stable", "stable"),
        ],
    },
    "thorium-232": {
        "name": "Thorium-232 Series (4n)",
        "parent": "Th-232",
        "stable_end": "Pb-208",
        "chain": [
            ("Th-232",   "α",   "1.40×10¹⁰ y"),
            ("Ra-228",   "β⁻",  "5.75 y"),
            ("Ac-228",   "β⁻",  "6.15 h"),
            ("Th-228",   "α",   "1.91 y"),
            ("Ra-224",   "α",   "3.66 d"),
            ("Rn-220",   "α",   "55.6 s"),
            ("Po-216",   "α",   "0.145 s"),
            ("Pb-212",   "β⁻/α", "10.64 h"),
            ("Bi-212",   "α/β⁻", "60.55 min"),
            ("Po-212",   "α",   "0.30 μs") if False else ("Tl-208", "β⁻", "3.05 min"),  # branch
            ("Pb-208",   "stable", "stable"),
        ],
    },
    "uranium-235": {
        "name": "Uranium-235 / Actinium Series (4n+3)",
        "parent": "U-235",
        "stable_end": "Pb-207",
        "chain": [
            ("U-235",    "α",   "7.04×10⁸ y"),
            ("Th-231",   "β⁻",  "25.5 h"),
            ("Pa-231",   "α",   "3.28×10⁴ y"),
            ("Ac-227",   "β⁻/α", "21.77 y"),
            ("Th-227",   "α",   "18.72 d"),
            ("Ra-223",   "α",   "11.43 d"),
            ("Rn-219",   "α",   "3.96 s"),
            ("Po-215",   "α",   "1.78 ms"),
            ("Pb-211",   "β⁻",  "36.1 min"),
            ("Bi-211",   "α/β⁻", "2.14 min"),
            ("Tl-207",   "β⁻",  "4.77 min"),
            ("Pb-207",   "stable", "stable"),
        ],
    },
}


@ChemMCPManager.register_tool
class DecaySeries(BaseTool):
    """
    放射性衰变系列查询工具。
    查询三大天然放射系（铀系、锕系、钍系）的完整衰变链。
    """
    __version__      = "0.1.0"
    name             = "DecaySeries"
    func_name        = "query_decay_series"
    description      = "Query radioactive decay series (Uranium-238, Thorium-232, Uranium-235) showing complete decay chains with nuclides, decay modes, and half-lives."
    implementation_description = "Uses built-in data for the three natural radioactive decay series (4n, 4n+2, 4n+3), returning full decay chains from parent to stable lead isotope."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Nuclear Chemistry", "Decay Series", "Radioactivity"]
    required_envs    = []

    code_input_sig   = [
        ("series_name", "str", "N/A", "Name of the decay series: 'uranium-238', 'thorium-232', 'uranium-235', or 'all'."),
        ("format_type", "str", "detailed", "Output format: 'summary' (brief overview) or 'detailed' (full chain)."),
    ]

    text_input_sig   = [
        ("query_str", "str", "N/A", "Series name and optional format, e.g., 'uranium-238 detailed' or 'all summary'."),
    ]

    output_sig       = [
        ("series_name", "str", "The queried series name(s)."),
        ("series_info", "dict", "Detailed information about the decay series including chain steps."),
    ]

    examples         = [
        {
            "code_input": {"series_name": "uranium-238", "format_type": "summary"},
            "text_input": {"query_str": "uranium-238 summary"},
            "output": {
                "series_name": "Uranium-238 Series (4n+2)",
                "series_info": {
                    "parent": "U-238",
                    "stable_end_product": "Pb-206",
                    "total_steps": 14,
                    "mass_number_rule": "4n + 2",
                    "alpha_decays": 8,
                    "beta_decays": 6,
                }
            }
        },
        {
            "code_input": {"series_name": "all", "format_type": "summary"},
            "text_input": {"query_str": "all summary"},
            "output": {
                "series_name": "All Three Natural Decay Series",
                "series_info": {
                    "uranium-238": {"parent": "U-238", "end": "Pb-206", "rule": "4n+2"},
                    "uranium-235": {"parent": "U-235", "end": "Pb-207", "rule": "4n+3"},
                    "thorium-232": {"parent": "Th-232", "end": "Pb-208", "rule": "4n"},
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, series_name: str, format_type: str = "detailed") -> dict:
        """Query decay series information."""
        key = series_name.strip().lower()

        if key == "all":
            result = {}
            for sk, sv in DECAY_SERIES.items():
                result[sk] = self._format_series(sv, format_type)
            return {
                "series_name": "All Three Natural Decay Series",
                "series_info": result,
            }

        # Find matching series
        matched = None
        for sk in DECAY_SERIES:
            if sk.startswith(key) or key in sk:
                matched = sk
                break

        if matched is None:
            available = ", ".join(DECAY_SERIES.keys())
            raise ChemMCPError(f"Unknown series '{series_name}'. Available: {available}")

        return {
            "series_name": DECAY_SERIES[matched]["name"],
            "series_info": self._format_series(DECAY_SERIES[matched], format_type),
        }

    def _format_series(self, series_data: dict, fmt: str) -> dict:
        """Format series data according to requested format."""
        chain = series_data["chain"]

        if fmt == "summary":
            alpha_count = sum(1 for _, mode, _ in chain if "α" in mode)
            beta_count = sum(1 for _, mode, _ in chain if "β" in mode)
            return {
                "parent": series_data["parent"],
                "stable_end_product": series_data["stable_end"],
                "total_steps": len([c for c in chain if c[1] != "stable"]),
                "mass_number_rule": series_data["name"].split("(")[-1].rstrip(")"),
                "alpha_decays": alpha_count,
                "beta_decays": beta_count,
            }
        else:
            # detailed format
            steps = []
            for i, (nuclide, mode, hl) in enumerate(chain):
                steps.append({
                    "step": i + 1,
                    "nuclide": nuclide,
                    "decay_mode": mode,
                    "half_life": str(hl),
                })
            return {
                "parent": series_data["parent"],
                "stable_end_product": series_data["stable_end"],
                "series_full_name": series_data["name"],
                "decay_chain": steps,
            }

    def _run_text(self, query_str: str) -> dict:
        parts = query_str.strip().split()
        series = parts[0] if parts else "all"
        fmt = parts[1] if len(parts) > 1 else "detailed"
        return self._run_base(series, fmt)
