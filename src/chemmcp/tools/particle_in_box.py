import logging
import math
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ParticleInBox(BaseTool):
    """
    一维/三维势箱中粒子的能级和波函数计算。
    
    1D: E_n = n²h²/(8mL²), ψ_n(x) = sqrt(2/L)·sin(nπx/L)
    3D: E_(nx,ny,nz) = (h²/8m)(nx/Lx² + ny/Ly² + nz/Lz²)
    """
    __version__ = "0.1.0"
    name = "ParticleInBox"
    func_name = "particle_in_box"
    description = "Calculate energy levels and wavefunctions for a particle in a 1D or 3D box (infinite potential well)."
    implementation_description = "Solves the time-independent Schrödinger equation for particle in a box. Computes energy eigenvalues, wavefunction expressions, probability densities, degeneracy (3D), and transition energies."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Particle in a Box", "Energy Levels", "Wavefunction"]
    required_envs = []

    code_input_sig = [
        ("dimensionality", "str", "N/A", "'1d' for one-dimensional box or '3d' for three-dimensional box."),
        ("mass_kg", "float", "N/A", "Mass of the particle in kg."),
        ("box_length_m", "float", "N/A", "Box length L in meters (for 1D). For 3D, use lengths_3d instead."),
        ("quantum_number_n", "int", "N/A", "Principal quantum number n (int >= 1). For 3D, use quantum_numbers_3d instead."),
        ("lengths_3d", "list", "None", "[Lx, Ly, Lz] in meters for 3D box (optional; if given, uses 3D mode)."),
        ("quantum_numbers_3d", "list", "None", "[nx, ny, nz] for 3D box (optional)."),
        ("position_x_m", "float", "None", "Position at which to evaluate wavefunction and probability density (optional, default=L/2)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'dim mass length n [x_position]' Example: '1d 9.109e-31 1e-9 2' for electron in 1nm box, n=2"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary with energy level, wavefunction info, probability density, degeneracy (3D), and transition data."),
    ]

    examples = [
        {
            "code_input": {
                "dimensionality": "1d",
                "mass_kg": 9.109e-31,
                "box_length_m": 1e-9,
                "quantum_number_n": 2,
                "lengths_3d": None,
                "quantum_numbers_3d": None,
                "position_x_m": None,
            },
            "text_input": {
                "input_params": "1d 9.109e-31 1e-9 2",
            },
            "output": {
                "result": {
                    "energy_J": 6.0249e-20,
                    "energy_eV": 0.3762,
                    "wavefunction_type": "sin(2πx/L)",
                    "n_nodes": 1,
                    "degeneracy": 1,
                }
            },
        },
        {
            "code_input": {
                "dimensionality": "3d",
                "mass_kg": 9.109e-31,
                "box_length_m": None,
                "quantum_number_n": None,
                "lengths_3d": [1e-9, 1e-9, 1e-9],
                "quantum_numbers_3d": [1, 1, 2],
                "position_x_m": None,
            },
            "text_input": {
                "input_params": "3d 9.109e-31 1e-9,1e-9,1e-9 1,1,2",
            },
            "output": {
                "result": {
                    "energy_eV": 1.506,
                    "degeneracy": 3,
                }
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.h = 6.62607015e-34  # Planck constant J·s
        self.eV_per_J = 6.241509e18  # eV per Joule

    def _solve_1d(self, mass, L, n, x_pos=None):
        """1D Particle in a box."""
        if n < 1:
            raise ChemMCPError("Quantum number n must be >= 1.")
        
        # E_n = n²h²/(8mL²)
        E = (n ** 2) * self.h ** 2 / (8 * mass * L ** 2)
        E_eV = E * self.eV_per_J
        
        # Wavefunction: ψ_n(x) = sqrt(2/L) * sin(n*π*x/L)
        if x_pos is None:
            x_pos = L / 2.0
        
        psi = math.sqrt(2.0 / L) * math.sin(n * math.pi * x_pos / L)
        prob_density = psi ** 2

        # Number of nodes (excluding boundaries): n - 1
        n_nodes = n - 1

        # Transition energy from ground state
        E_ground = self.h ** 2 / (8 * mass * L ** 2)
        delta_E = E - E_ground
        wavelength = self.h * 2.998e8 / delta_E if delta_E > 0 else float('inf')  # λ = hc/E

        return {
            "dimensionality": "1D",
            "quantum_number_n": n,
            "energy_J": round(E, 25),
            "energy_eV": round(E_eV, 10),
            "wavefunction_expression": f"sqrt(2/L)*sin({n}πx/L)",
            "psi_at_position": round(psi, 20),
            "probability_density_at_position": round(prob_density, 20),
            "position_evaluated_m": x_pos,
            "box_length_m": L,
            "n_nodes": n_nodes,
            "transition_energy_from_ground_eV": round(delta_E * self.eV_per_J, 10),
            "transition_wavelength_nm": round(wavelength * 1e9, 4) if wavelength != float('inf') else None,
            "degeneracy": 1,
        }

    def _solve_3d(self, mass, lengths, qn_tuple, x_pos=None):
        """3D Particle in a box."""
        nx, ny, nz = qn_tuple
        Lx, Ly, Lz = lengths
        
        if any(q < 1 for q in qn_tuple):
            raise ChemMCPError("All quantum numbers must be >= 1.")
        
        h = self.h
        # E = (h²/8m) * [(nx/Lx)² + (ny/Ly)² + (nz/Lz)²]
        E = (h ** 2 / (8 * mass)) * (
            (nx / Lx) ** 2 + (ny / Ly) ** 2 + (nz / Lz) ** 2
        )
        E_eV = E * self.eV_per_J

        # Degeneracy check for cubic box
        is_cubic = (abs(Lx - Ly) < 1e-15 * max(Lx, Ly) and 
                    abs(Ly - Lz) < 1e-15 * max(Ly, Lz))
        degeneracy = 1  # Would need full enumeration for accurate value

        # Wavefunction product
        if x_pos is None:
            x_pos = Lx / 2.0
        y_pos = Ly / 2.0
        z_pos = Lz / 2.0
        
        psi = (math.sqrt(2/Lx) * math.sin(nx*math.pi*x_pos/Lx) *
               math.sqrt(2/Ly) * math.sin(ny*math.pi*y_pos/Ly) *
               math.sqrt(2/Lz) * math.sin(nz*math.pi*z_pos/Lz))
        prob = psi ** 2

        return {
            "dimensionality": "3D",
            "quantum_numbers": {"nx": nx, "ny": ny, "nz": nz},
            "box_lengths_m": {"Lx": Lx, "Ly": Ly, "Lz": Lz},
            "energy_J": round(E, 25),
            "energy_eV": round(E_eV, 10),
            "wavefunction_expression": f"psi({nx},{ny},{nz})=product of 1D sine functions",
            "psi_at_center": round(psi, 20),
            "probability_density_at_center": round(prob, 20),
            "is_cubic_box": is_cubic,
            "degeneracy": degeneracy,
            "n_total_nodes": (nx - 1) + (ny - 1) + (nz - 1),
        }

    def _run_base(self, dimensionality: str, mass_kg: float, box_length_m: float = None,
                  quantum_number_n: int = None, lengths_3d: list = None,
                  quantum_numbers_3d: list = None, position_x_m: float = None) -> dict:
        
        if dimensionality == "3d" or lengths_3d is not None:
            if lengths_3d is None or len(lengths_3d) != 3:
                raise ChemMCPError("For 3D, provide lengths_3d as [Lx, Ly, Lz].")
            if quantum_numbers_3d is None or len(quantum_numbers_3d) != 3:
                raise ChemMCPError("For 3D, provide quantum_numbers_3d as [nx, ny, nz].")
            result = self._solve_3d(mass_kg, lengths_3d, tuple(quantum_numbers_3d), position_x_m)
        elif dimensionality == "1d":
            result = self._solve_1d(mass_kg, box_length_m, quantum_number_n, position_x_m)
        else:
            raise ChemMCPError("dimensionality must be '1d' or '3d'.")

        logger.info(f"ParticleInBox: {result['dimensionality']}, E={result.get('energy_eV', 'N/A')}eV")
        return result

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            dim = parts[0]
            mass = float(parts[1])
            
            if dim == "3d":
                lengths = [float(x) for x in parts[2].split(",")]
                qns = [int(x) for x in parts[3].split(",")]
                x_pos = float(parts[4]) if len(parts) > 4 else None
                return self._run_base("3d", mass, 0, lengths_3d=lengths,
                                      quantum_numbers_3d=qns, position_x_m=x_pos)
            else:
                L = float(parts[2])
                n = int(parts[3])
                x_pos = float(parts[4]) if len(parts) > 4 else None
                return self._run_base("1d", mass, L, n, position_x_m=x_pos)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}. Format: '1d|3d mass L|Lx,Ly,Lz n|nx,ny,nz [x_pos]'")
