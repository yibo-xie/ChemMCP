import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BondEnergyCalculation(BaseTool):
    """
    用键能估算反应焓变。
    ΔH = Σ(断裂反应物中的键能) - Σ(形成产物中的键能)
    """
    __version__                = "0.1.0"
    name                       = "BondEnergyCalculation"
    func_name                  = "bond_energy_calculation"
    description                = "Estimate reaction enthalpy change using bond dissociation energies."
    implementation_description = "ΔH_rxn = Σ(bonds broken in reactants) - Σ(bonds formed in products). All energies in kJ/mol."
    oss_dependencies           = []
    services_and_software      = []
    categories                 = ["General"]
    tags                       = ["Thermodynamics", "Bond Energy", "Enthalpy", "Reaction"]
    required_envs              = []

    code_input_sig             = [
        ("bonds_broken",      "list",  "N/A",   "List of (bond_name, bond_energy_kj_per_mol, count) for bonds broken in reactants."),
        ("bonds_formed",      "list",  "N/A",   "List of (bond_name, bond_energy_kj_per_mol, count) for bonds formed in products."),
    ]

    text_input_sig             = [
        ("input_params",      "str",   "N/A",   "Semi-structured string describing broken and formed bonds with energies."),
    ]

    output_sig                 = [
        ("delta_h",           "float", "Estimated reaction enthalpy change ΔH in kJ/mol."),
        ("total_broken",      "float", "Total energy of bonds broken (kJ/mol)."),
        ("total_formed",      "float", "Total energy of bonds formed (kJ/mol)."),
        ("reaction_type",     "str",   "'exothermic' if ΔH < 0, 'endothermic' if ΔH > 0, 'thermoneutral' if ≈0."),
    ]

    examples                   = [
        {
            "code_input": {
                "bonds_broken": [("N≡N", 941, 1), ("H-H", 436, 3)],
                "bonds_formed": [("N-H", 391, 6)],
            },
            "text_input": {
                "input_params": "broken: N≡N=941x1 H-H=436x3; formed: N-H=391x6",
            },
            "output": {
                "delta_h": -93.0,
                "total_broken": 2249.0,
                "total_formed": 2346.0,
                "reaction_type": "exothermic",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, bonds_broken: list, bonds_formed: list) -> dict:
        total_broken = 0.0
        for name, energy, count in bonds_broken:
            total_broken += energy * count

        total_formed = 0.0
        for name, energy, count in bonds_formed:
            total_formed += energy * count

        delta_h = total_broken - total_formed

        if abs(delta_h) < 1e-6:
            reaction_type = "thermoneutral"
        elif delta_h < 0:
            reaction_type = "exothermic"
        else:
            reaction_type = "endothermic"

        logger.info(f"Bond energy: broken={total_broken}, formed={total_formed}, ΔH={delta_h} ({reaction_type})")
        return {
            "delta_h": round(delta_h, 2),
            "total_broken": round(total_broken, 2),
            "total_formed": round(total_formed, 2),
            "reaction_type": reaction_type,
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            # Parse format: "broken: A=B=100x2 C-D=200x1; formed: E-F=300x3"
            broken = []
            formed = []
            sections = input_params.split(";")
            for sec in sections:
                sec = sec.strip()
                if sec.startswith("broken:"):
                    part = sec[len("broken:"):].strip()
                    for item in part.split():
                        # Format: BondName=EnergyxCount
                        eq_idx = item.rfind("=")
                        x_idx = item.find("x")
                        name = item[:eq_idx]
                        energy = float(item[eq_idx+1:x_idx])
                        count = int(item[x_idx+1:])
                        broken.append((name, energy, count))
                elif sec.startswith("formed:"):
                    part = sec[len("formed:"):].strip()
                    for item in part.split():
                        eq_idx = item.rfind("=")
                        x_idx = item.find("x")
                        name = item[:eq_idx]
                        energy = float(item[eq_idx+1:x_idx])
                        count = int(item[x_idx+1:])
                        formed.append((name, energy, count))
            return self._run_base(broken, formed)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Expected format: 'broken: Bond=E*x ...; formed: Bond=E*x ...'")
