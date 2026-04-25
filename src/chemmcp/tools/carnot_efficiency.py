import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

R = 8.314  # J/(mol·K)


@ChemMCPManager.register_tool
class CarnotEfficiency(BaseTool):
    """
    计算卡诺循环效率及各过程功热。
    
    卡诺循环四个过程：
    1. 等温膨胀 (1→2): T=Th恒定, 吸热 Qh = n·R·Th·ln(V2/V1)
    2. 绝热膨胀 (2→3): Q=0, T从Th降到Tc
    3. 等温压缩 (3→4): T=Tc恒定, 放热 Qc = n·R·Tc·ln(V4/V3)
    4. 绝热压缩 (4→1): Q=0, T从Tc升到Th
    
    效率: η = 1 - Tc/Th
    """
    __version__ = "0.1.0"
    name = "CarnotEfficiency"
    func_name = "calculate_carnot_efficiency"
    description = "Calculate Carnot cycle efficiency and work/heat for each of the four processes."
    implementation_description = "Computes efficiency η = 1 - Tc/Th, heat absorbed (Qh), heat released (Qc), net work (W_net), and details for all four processes using ideal gas relations."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Thermodynamics", "Carnot Cycle", "Heat Engine", "Efficiency", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("th", "float", "N/A", "Hot reservoir temperature in Kelvin."),
        ("tc", "float", "N/A", "Cold reservoir temperature in Kelvin."),
        ("n_moles", "float", "1.0", "Amount of gas in moles."),
        ("v1", "float", "0.01", "Initial volume V1 in liters (for process calculations)."),
        ("compression_ratio", "float", "5.0", "Compression ratio r = V2/V1 for isothermal expansion."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'th tc [n_moles] [v1] [compression_ratio]'. Example: '600 300 1.0 0.01 5'"),
    ]

    output_sig = [
        ("efficiency", "float", "Carnot cycle efficiency η (dimensionless, 0-1)."),
        ("efficiency_percent", "float", "Efficiency as percentage (%)."),
        ("q_h", "float", "Heat absorbed from hot reservoir Qh in J."),
        ("q_c", "float", "Heat released to cold reservoir Qc in J (absolute value)."),
        ("w_net", "float", "Net work output W_net in J."),
        ("processes", "str", "Details of the four processes."),
        ("explanation", "str", "Summary explanation."),
    ]

    examples = [
        {
            "code_input": {
                "th": 600.0,
                "tc": 300.0,
                "n_moles": 1.0,
                "v1": 0.01,
                "compression_ratio": 5.0,
            },
            "text_input": {
                "input_params": "600 300",
            },
            "output": {
                "efficiency": 0.5,
                "efficiency_percent": 50.0,
                "q_h": 806.2,
                "q_c": 403.1,
                "w_net": 403.1,
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, th: float, tc: float, n_moles: float = 1.0, v1: float = 0.01, compression_ratio: float = 5.0) -> dict:
        """Core logic: calculate Carnot cycle parameters."""
        if th <= 0 or tc <= 0:
            raise ChemMCPError("Temperatures must be positive in Kelvin.")
        if tc >= th:
            raise ChemMCPError(f"Hot temperature ({th} K) must be greater than cold temperature ({tc} K).")
        if compression_ratio <= 0:
            raise ChemMCPError("Compression ratio must be positive.")
        if v1 <= 0:
            raise ChemMCPError("Initial volume must be positive.")

        # Efficiency
        eta = 1.0 - tc / th

        # Process calculations (using R in L·kPa/(K·mol) for volume in L → J conversion)
        # Use R = 8.314 J/(mol·K) directly; volume needs to be consistent
        # For ideal gas: PV = nRT, so V in m³ gives P in Pa
        # Let's use V in L, convert: 1 L = 0.001 m³
        V1 = v1
        r = compression_ratio

        # Process 1→2: Isothermal expansion at Th
        V2 = V1 * r
        Qh = n_moles * R * th * math.log(r)  # J (positive, heat absorbed)

        # Process 2→3: Adiabatic expansion, Th → Tc
        # TV^(γ-1) = const, but we need γ. Assume diatomic gas: γ = 7/5 = 1.4
        gamma = 1.4  # diatomic ideal gas
        # From adiabatic relation: Th*V2^(γ-1) = Tc*V3^(γ-1)
        V3 = V2 * (th / tc) ** (1.0 / (gamma - 1))

        # Process 3→4: Isothermal compression at Tc
        # Adiabatic 4→1: Tc*V4^(γ-1) = Th*V1^(γ-1)
        V4 = V1 * (th / tc) ** (1.0 / (gamma - 1))

        Qc = abs(n_moles * R * tc * math.log(V4 / V3))  # J (absolute value, heat released)

        # Net work
        W_net = Qh - Qc

        # Process details
        processes = (
            f"Process 1→2 (Isothermal expansion @ {th} K):\n"
            f"  V1={V1:.4f} L → V2={V2:.4f} L\n"
            f"  Qh = nRT·ln(V2/V1) = {n_moles}×{R}×{th}×ln({r}) = {Qh:.2f} J (absorbed)\n\n"
            f"Process 2→3 (Adiabatic expansion):\n"
            f"  T: {th} K → {tc} K, V2={V2:.4f} L → V3={V3:.4f} L\n"
            f"  Q = 0 (adiabatic)\n\n"
            f"Process 3→4 (Isothermal compression @ {tc} K):\n"
            f"  V3={V3:.4f} L → V4={V4:.4f} L\n"
            f"  Qc = |nRT·ln(V4/V3)| = {Qc:.2f} J (released)\n\n"
            f"Process 4→1 (Adiabatic compression):\n"
            f"  T: {tc} K → {th} K, V4={V4:.4f} L → V1={V1:.4f} L\n"
            f"  Q = 0 (adiabatic)"
        )

        explanation = (
            f"Carnot Cycle Summary:\n"
            f"Th = {th} K, Tc = {tc} K\n"
            f"η = 1 - Tc/Th = 1 - {tc}/{th} = {eta:.4f} ({eta*100:.2f}%)\n"
            f"Qh = {Qh:.2f} J, Qc = {Qc:.2f} J, W_net = {W_net:.2f} J\n"
            f"W_net/Qh = {W_net/Qh:.4f} (= η ✓)"
        )

        logger.info(f"Carnot: Th={th}K, Tc={tc}K → η={eta:.4f}, W_net={W_net:.2f}J")
        return {
            "efficiency": round(eta, 6),
            "efficiency_percent": round(eta * 100, 4),
            "q_h": round(Qh, 2),
            "q_c": round(Qc, 2),
            "w_net": round(W_net, 2),
            "processes": processes,
            "explanation": explanation,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            th = float(parts[0])
            tc = float(parts[1])
            n = float(parts[2]) if len(parts) > 2 else 1.0
            v1 = float(parts[3]) if len(parts) > 3 else 0.01
            cr = float(parts[4]) if len(parts) > 4 else 5.0
            return self._run_base(th, tc, n, v1, cr)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'th tc [n_moles] [v1] [compression_ratio]'")
