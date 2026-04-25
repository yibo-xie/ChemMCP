import logging
import math
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Common quantum mechanical operators and their expectation value formulas
# for standard model systems (particle in a box, harmonic oscillator, hydrogen atom)

@ChemMCPManager.register_tool
class ExpectationValue(BaseTool):
    """
    力学量期望值计算工具。
    支持一维势箱、谐振子、氢原子等标准量子系统的各种力学量期望值计算。
    """
    __version__ = "0.1.0"
    name = "ExpectationValue"
    func_name = "calculate_expectation_value"
    description = "Calculate expectation values of various mechanical observables for standard quantum systems."
    implementation_description = "Uses analytical solutions of the Schrödinger equation for particle-in-a-box, harmonic oscillator, and hydrogen atom to compute <x>, <p>, <x²>, <p²>, <E>, <L>, <L_z> etc."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Quantum Mechanics", "Expectation Value", "Wave Function", "Physical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("system", "str", "N/A", "Quantum system: 'particle_in_box', 'harmonic_oscillator', or 'hydrogen_atom'."),
        ("observable", "str", "N/A", "Observable to compute: 'x', 'x2', 'p', 'p2', 'E', 'T', 'V', 'L', 'Lz', 'L2', 'r', 'r2', '1/r'."),
        ("n", "int", "1", "Principal quantum number (or state index)."),
        ("l", "int", "0", "Angular momentum quantum number (for hydrogen atom)."),
        ("m", "int", "0", "Magnetic quantum number (for hydrogen atom)."),
        ("L", "float", "1.0", "Box length (m) for particle_in_box, or oscillator length scale (m)."),
        ("mass_kg", "float", "9.109e-31", "Particle mass in kg."),
        ("omega", "float", "1e14", "Angular frequency (rad/s) for harmonic oscillator."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: system observable n [l m] L [mass_kg omega]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing the expectation value, units, formula used, and physical interpretation."),
    ]

    examples = [
        {
            "code_input": {
                "system": "particle_in_box",
                "observable": "x",
                "n": 1,
                "l": 0,
                "m": 0,
                "L": 1e-9,
                "mass_kg": 9.109e-31,
                "omega": 0,
            },
            "text_input": {
                "input_params": "particle_in_box x 1 0 0 1e-9",
            },
            "output": {
                "result": {
                    "system": "particle_in_box",
                    "observable": "x",
                    "quantum_state_n": 1,
                    "expectation_value": 5e-10,
                    "units": "m",
                    "formula": "<x> = L/2",
                    "interpretation": "For any stationary state of particle in a box, <x> is at the center.",
                }
            }
        },
        {
            "code_input": {
                "system": "harmonic_oscillator",
                "observable": "E",
                "n": 2,
                "l": 0,
                "m": 0,
                "L": 0,
                "mass_kg": 9.109e-31,
                "omega": 1e14,
            },
            "text_input": {
                "input_params": "harmonic_oscillator E 2 0 0 0 9.109e-31 1e14",
            },
            "output": {
                "result": {
                    "system": "harmonic_oscillator",
                    "observable": "E",
                    "quantum_state_n": 2,
                    "expectation_value": 2.5 * 1.054571817e-34 * 1e14,
                    "units": "J",
                    "formula": "<E> = (n + 1/2)ℏω",
                    "interpretation": "Energy expectation equals eigenvalue for stationary states.",
                }
            }
        },
        {
            "code_input": {
                "system": "hydrogen_atom",
                "observable": "r",
                "n": 1,
                "l": 0,
                "m": 0,
                "L": 0,
                "mass_kg": 0,
                "omega": 0,
            },
            "text_input": {
                "input_params": "hydrogen_atom r 1 0 0",
            },
            "output": {
                "result": {
                    "system": "hydrogen_atom",
                    "observable": "r",
                    "quantum_state_n": 1,
                    "l": 0,
                    "expectation_value_m": 7.937e-11,
                    "units": "m",
                    "formula": "<r> = a₀[3n² - l(l+1)] / 2",
                    "interpretation": "Ground state mean electron-proton distance is 1.5·a₀ ≈ 79.4 pm.",
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self.hbar = 1.054571817e-34   # J·s
        self.a0 = 5.29177210903e-11   # Bohr radius, m
        self.eV = 1.602176634e-19      # J per eV
        self.Ry = 13.605693122994      # Rydberg constant in eV

    def _run_base(self, system: str, observable: str, n: int = 1, l: int = 0, m: int = 0,
                  L: float = 1.0, mass_kg: float = 9.109e-31, omega: float = 1e14) -> dict:
        """Core logic: compute expectation value."""
        sys = system.lower().replace(" ", "_")
        obs = observable.lower().strip()

        if sys == "particle_in_box":
            return self._pib(obs, n, L, mass_kg)
        elif sys == "harmonic_oscillator":
            return self._ho(obs, n, mass_kg, omega)
        elif sys == "hydrogen_atom":
            return self._hydrogen(obs, n, l, m)
        else:
            raise ChemMCPError(
                f"Unknown system '{system}'. Use 'particle_in_box', 'harmonic_oscillator', or 'hydrogen_atom'."
            )

    # ── Particle in 1D Box (0 to L) ──────────────────────────────────
    def _pib(self, obs: str, n: int, L: float, mass: float) -> dict:
        if n < 1:
            raise ChemMCPError("Quantum number n must be >= 1.")
        if L <= 0:
            raise ChemMCPError("Box length L must be positive.")

        En = (n ** 2 * math.pi ** 2 * self.hbar ** 2) / (2 * mass * L ** 2)

        if obs == "x":
            val = L / 2.0
            fmt = "<x> = L/2"
        elif obs == "x2":
            val = L ** 2 * (1 / 3 - 1 / (2 * math.pi ** 2 * n ** 2))
            fmt = "<x²> = L²(1/3 - 1/(2π²n²))"
        elif obs == "p":
            val = 0.0
            fmt = "<p> = 0"
        elif obs == "p2":
            val = (n * math.pi * self.hbar / L) ** 2
            fmt = "<p²> = (nπℏ/L)²"
        elif obs in ("e", "energy"):
            En = (n ** 2 * math.pi ** 2 * self.hbar ** 2) / (2 * mass * L ** 2)
            val = En
            fmt = f"<E> = Eₙ = n²π²ℏ²/(2mL²)"
        elif obs in ("t", "kinetic"):
            En = (n ** 2 * math.pi ** 2 * self.hbar ** 2) / (2 * mass * L ** 2)
            val = En  # V=0 inside box
            fmt = "<T> = Eₙ"
        elif obs == "delta_x":
            dx_sq = L ** 2 * (1 / 12 - 1 / (2 * math.pi ** 2 * n ** 2))
            val = math.sqrt(dx_sq)
            fmt = "Δx = √(<x²> - <x>²)"
        elif obs == "delta_p":
            dp_sq = (n * math.pi * self.hbar / L) ** 2
            val = math.sqrt(dp_sq)
            fmt = "Δp = |<p²>|^(1/2)"
        else:
            raise ChemMCPError(f"Unknown observable '{obs}' for particle_in_box.")

        return self._make_result("particle_in_box", obs, val, n=n, L=L, formula=fmt)

    # ── 1D Harmonic Oscillator ───────────────────────────────────────
    def _ho(self, obs: str, n: int, mass: float, omega: float) -> dict:
        if n < 0:
            raise ChemMCPError("Quantum number n must be >= 0.")
        if omega <= 0:
            raise ChemMCPError("Angular frequency omega must be positive.")
        if mass <= 0:
            raise ChemMCPError("Mass must be positive.")

        En = (n + 0.5) * self.hbar * omega

        if obs == "x":
            val = 0.0
            fmt = "<x> = 0"
        elif obs == "x2":
            val = (n + 0.5) * self.hbar / (mass * omega)
            fmt = "<x²> = (n+½)ℏ/(mω)"
        elif obs == "p":
            val = 0.0
            fmt = "<p> = 0"
        elif obs == "p2":
            val = (n + 0.5) * self.hbar * mass * omega
            fmt = "<p²> = (n+½)ℏmω"
        elif obs in ("e", "energy"):
            val = En
            fmt = "<E> = (n+½)ℏω"
        elif obs in ("t", "kinetic"):
            val = En / 2.0
            fmt = "<T> = <V> = E/2 (virial theorem)"
        elif obs in ("v", "potential"):
            val = En / 2.0
            fmt = "<V> = <T> = E/2 (virial theorem)"
        elif obs == "delta_x":
            val = math.sqrt((n + 0.5) * self.hbar / (mass * omega))
            fmt = "Δx = √(<x²>)"
        elif obs == "delta_p":
            val = math.sqrt((n + 0.5) * self.hbar * mass * omega)
            fmt = "Δp = √(<p²>)"
        else:
            raise ChemMCPError(f"Unknown observable '{obs}' for harmonic_oscillator.")

        return self._make_result("harmonic_oscillator", obs, val, n=n, omega=omega, formula=fmt)

    # ── Hydrogen Atom ────────────────────────────────────────────────
    def _hydrogen(self, obs: str, n: int, l: int, m: int) -> dict:
        if n < 1:
            raise ChemMCPError("Principal quantum number n must be >= 1.")
        if l < 0 or l >= n:
            raise ChemMCPError(f"Angular momentum l must satisfy 0 ≤ l < n (got l={l}, n={n}).")
        if abs(m) > l:
            raise ChemMCPError(f"|m| must be ≤ l (got m={m}, l={l}).")

        En_eV = -self.Ry / (n * n)

        if obs in ("e", "energy"):
            val = En_eV * self.eV
            fmt = f"<E> = -Ry/n² = {En_eV:.6f} eV"
            return self._make_result("hydrogen_atom", obs, val, n=n, l=l, m=m,
                                     formula=fmt, units="J", extra={"value_eV": round(En_eV, 8)})
        elif obs == "r":
            # <r> = (a0/2)[3n² - l(l+1)]
            val = (self.a0 / 2) * (3 * n * n - l * (l + 1))
            fmt = f"<r> = a₀/2 · [3n² - l(l+1)]"
        elif obs == "r2":
            # <r²> = (a0² n² / 2)[5n² + 1 - 3l(l+1)]
            val = (self.a0 ** 2 * n * n / 2) * (5 * n * n + 1 - 3 * l * (l + 1))
            fmt = "<r²> = a₀²n²/2 · [5n²+1-3l(l+1)]"
        elif obs == "1_over_r":
            val = 1 / (self.a0 * n * n)
            fmt = "<1/r> = 1/(a₀n²)"
        elif obs in ("l2", "l_squared"):
            val = l * (l + 1) * self.hbar ** 2
            fmt = "<L²> = l(l+1)ℏ²"
        elif obs in ("lz", "l_z"):
            val = m * self.hbar
            fmt = "<L_z> = mℏ"
        elif obs in ("l", "angular_momentum_mag"):
            val = math.sqrt(l * (l + 1)) * self.hbar
            fmt = "|L| = √(l(l+1)) · ℏ"
        else:
            raise ChemMCPError(f"Unknown observable '{obs}' for hydrogen_atom.")

        return self._make_result("hydrogen_atom", obs, val, n=n, l=l, m=m, formula=fmt)

    # ── Helper ───────────────────────────────────────────────────────
    def _make_result(self, system: str, obs: str, value: float, **kw) -> dict:
        result = {
            "system": system,
            "observable": obs,
            "expectation_value": f"{value:.6e}" if abs(value) < 1e-10 or abs(value) > 1e10 else round(value, 15),
            "units": kw.get("units", self._default_units(system, obs)),
            "formula": kw.get("formula", ""),
        }
        # Add quantum numbers if present
        for qk in ("n", "l", "m"):
            if qk in kw and kw[qk] is not None:
                result[f"quantum_{qk}"] = kw[qk]
        if kw.get("extra"):
            result.update(kw["extra"])
        result["interpretation"] = self._interpret(system, obs, value, **kw)
        return {"result": result}

    @staticmethod
    def _default_units(sys: str, obs: str) -> str:
        unit_map = {
            "x": "m", "x2": "m²", "p": "kg·m/s", "p2": "kg²·m²/s²",
            "E": "J", "energy": "J", "T": "J", "kinetic": "J", "V": "J", "potential": "J",
            "delta_x": "m", "delta_p": "kg·m/s",
            "r": "m", "r2": "m²", "1_over_r": "1/m",
            "L2": "J²·s²", "l_squared": "J²·s²",
            "Lz": "J·s", "l_z": "J·s", "L": "J·s", "angular_momentum_mag": "J·s",
        }
        return unit_map.get(obs, "unknown")

    def _interpret(self, sys: str, obs: str, val: float, **kw) -> str:
        """Generate human-readable interpretation."""
        if sys == "particle_in_box":
            if obs == "x":
                return f"For state n={kw.get('n')}, the particle is equally likely found anywhere; average position is at center L/2 = {val:.3e} m."
            elif obs == "E":
                return f"Ground state energy of particle in box: {val:.3e} J ({val/self.eV:.3f} eV)."
        elif sys == "harmonic_oscillator":
            if obs == "E":
                ev_val = val / self.eV
                return f"Zero-point energy included: E_{kw.get('n')} = ({kw.get('n')}+½)ℏω = {ev_val:.6f} eV."
            elif obs in ("T", "V"):
                return f"By virial theorem, <T> = <V> = E/2 = {val:.3e} J."
        elif sys == "hydrogen_atom":
            if obs == "r":
                pm = val * 1e12
                return f"Mean electron-nucleus distance for n={kw.get('n')}, l={kw.get('l')}: {pm:.3f} pm ({val/self.a0:.3f} a₀)."
            elif obs == "E":
                eV_val = val / self.eV
                return f"Binding energy for n={kw.get('n')}: {eV_val:.6f} eV."
        return f"Expectation value of {obs} in {sys}: {val:.6e}"

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            system = parts[0]
            observable = parts[1]
            n = int(parts[2]) if len(parts) > 2 else 1
            l = int(parts[3]) if len(parts) > 3 else 0
            m = int(parts[4]) if len(parts) > 4 else 0
            L = float(parts[5]) if len(parts) > 5 else 1.0
            mass = float(parts[6]) if len(parts) > 6 else 9.109e-31
            omega = float(parts[7]) if len(parts) > 7 else 1e14
            return self._run_base(system, observable, n, l, m, L, mass, omega)
        except (IndexError, ValueError) as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
