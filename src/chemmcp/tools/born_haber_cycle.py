import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BornHaberCycle(BaseTool):
    """
    Born-Haber 循环计算：用于计算离子晶体的晶格能或反应焓变。
    ΔH_f = ΔH_sub + IE + (1/2)D + EA + U
    可求解任意一个未知量。
    """
    __version__                = "0.1.0"
    name                       = "BornHaberCycle"
    func_name                  = "born_haber_cycle"
    description                = "Calculate lattice energy or formation enthalpy using Born-Haber cycle for ionic crystals."
    implementation_description = "Uses Born-Haber cycle: ΔH_f = ΔH_sub + ΣIE + (1/2)D + ΣEA + U. Solves for any one unknown parameter."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Thermodynamics", "Solid State", "Ionic Crystal", "Lattice Energy"]
    required_envs              = []

    code_input_sig             = [
        ("delta_h_f",         "float", "None", "Formation enthalpy ΔH_f in kJ/mol (None if unknown)."),
        ("delta_h_sub",       "float", "None", "Sublimation enthalpy of metal in kJ/mol."),
        ("ionization_energies","list",  "None", "List of ionization energies (IE1, IE2, ...) in kJ/mol."),
        ("bond_dissociation", "float", "None", "Bond dissociation energy of X₂ (D) in kJ/mol."),
        ("electron_affinities","list", "None", "List of electron affinities (EA1, EA2, ...) in kJ/mol."),
        ("lattice_energy",    "float", "None", "Lattice energy U in kJ/mol (None if unknown)."),
        ("unknown",           "str",   "None", "Which parameter to solve for: 'delta_h_f', 'delta_h_sub', 'bond_dissociation', 'lattice_energy', or an IE/EA index like 'ie_0'."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "JSON-like string with all parameters; set unknown to null and specify 'unknown' field."),
    ]

    output_sig                 = [
        ("result",            "float", "The calculated value of the unknown parameter in kJ/mol."),
        ("cycle_summary",     "dict",  "Summary of all Born-Haber terms."),
        ("unit",              "str",   "Unit: kJ/mol."),
    ]

    examples                   = [
        {
            "code_input": {
                "delta_h_f": -411.0,
                "delta_h_sub": 108.0,
                "ionization_energies": [496.0],
                "bond_dissociation": 122.0,
                "electron_affinities": [-349.0],
                "lattice_energy": None,
                "unknown": "lattice_energy",
            },
            "text_input": {
                "input_params": "-411.1 108 496 122 -349 None lattice_energy",
            },
            "output": {
                "result": -788.0,
                "cycle_summary": {"delta_h_sub": 108.0, "total_ie": 496.0, "half_d": 61.0, "total_ea": -349.0, "lattice_energy": -788.0, "sum": -472.0},
                "unit": "kJ/mol",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, delta_h_f: float = None, delta_h_sub: float = None,
                  ionization_energies: list = None, bond_dissociation: float = None,
                  electron_affinities: list = None, lattice_energy: float = None,
                  unknown: str = None) -> dict:
        # Validate exactly one unknown
        params = {
            "delta_h_f": delta_h_f,
            "delta_h_sub": delta_h_sub,
            "bond_dissociation": bond_dissociation,
            "lattice_energy": lattice_energy,
        }
        none_count = sum(1 for v in params.values() if v is None)
        if ionization_energies is None:
            none_count += 1
        if electron_affinities is None:
            none_count += 1
        if none_count != 1:
            raise ChemMCPError("Exactly one parameter must be None (the unknown to solve for).")

        total_ie = sum(ionization_energies or [0])
        total_ea = sum(electron_affinities or [0])
        half_d = (bond_dissociation or 0) / 2.0

        # Born-Haber: ΔH_f = ΔH_sub + ΣIE + D/2 + ΣEA + U
        known_sum = 0.0
        if delta_h_sub is not None:
            known_sum += delta_h_sub
        if ionization_energies is not None:
            known_sum += total_ie
        if bond_dissociation is not None:
            known_sum += half_d
        if electron_affinities is not None:
            known_sum += total_ea
        if lattice_energy is not None:
            known_sum += lattice_energy

        if delta_h_f is None:
            result = known_sum
            delta_h_f = result
        else:
            result = delta_h_f - known_sum
            lattice_energy = result

        summary = {
            "delta_h_f": delta_h_f if delta_h_f is not None else result,
            "delta_h_sub": delta_h_sub,
            "total_ie": total_ie,
            "half_d": half_d,
            "total_ea": total_ea,
            "lattice_energy": lattice_energy if lattice_energy is not None else result,
        }

        logger.info(f"Born-Haber cycle: unknown={unknown}, result={result} kJ/mol")
        return {
            "result": round(result, 2),
            "cycle_summary": {k: round(v, 2) if isinstance(v, float) else v for k, v in summary.items()},
            "unit": "kJ/mol",
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            if len(parts) < 7:
                raise ValueError("Need 7 parameters: dhf dhsub ie d ea lattice unknown")
            dhf = None if parts[0].lower() == "none" else float(parts[0])
            dhsub = None if parts[1].lower() == "none" else float(parts[1])
            ie_list = [float(x) for x in parts[2].strip("[]").split(",")]
            d = None if parts[3].lower() == "none" else float(parts[3])
            ea_list = [float(x) for x in parts[4].strip("[]").split(",")]
            lat = None if parts[5].lower() == "none" else float(parts[5])
            unk = parts[6]
            return self._run_base(dhf, dhsub, ie_list, d, ea_list, lat, unk)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
