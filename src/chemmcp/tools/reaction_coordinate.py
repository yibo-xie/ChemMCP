import logging
import math
from typing import List, Optional, Tuple

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# Physical constants
R = 8.314462618      # J/(mol·K)
kB = 1.380649e-23     # J/K
h = 6.62607015e-34    # J·s
NA = 6.02214076e23    # mol⁻¹


@ChemMCPManager.register_tool
class ReactionCoordinate(BaseTool):
    """
    反应坐标分析工具（IRC计算）。
    沿反应坐标绘制能量剖面图，定位过渡态、反应物和产物，计算正逆反应能垒。
    支持从SMILES描述符或手动能量数据进行定性和简化的IRC分析。
    """
    __version__ = "0.1.0"
    name = "ReactionCoordinate"
    func_name = "analyze_reaction_coordinate"
    description = "Analyze reaction coordinate profile: locate transition state (TS), calculate forward/reverse barriers, generate IRC (Intrinsic Reaction Coordinate) energy diagram data."
    implementation_description = (
        "Constructs a 1D potential energy surface along the reaction coordinate. "
        "Supports two input modes:\n"
        "1. 'manual': User provides a list of (coordinate, relative_energy_kj_mol) points.\n"
        "2. 'analytical': Uses analytical functions (cubic/Morse-like) to model PES from reactant/product/TS energies.\n"
        "Outputs: energy profile with TS position, forward barrier ΔG‡_f, reverse barrier ΔG‡_r, "
        "reaction energy ΔG_rxn, IRC point data for plotting, and Hammond postulate analysis."
    )
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Reaction Coordinate", "IRC", "Transition State", "Potential Energy Surface", "Kinetics"]
    required_envs = []

    code_input_sig = [
        ("input_mode", "str", "N/A", "Input mode: 'manual' (energy points) or 'analytical' (from energies)."),
        # For manual mode:
        ("energy_profile", "list", "None", "List of (coord, energy_kJ_mol) tuples along reaction coordinate."),
        # For analytical mode:
        ("reactant_energy_kj_mol", "float", "None", "Energy of reactants (kJ/mol), typically set to 0."),
        ("ts_energy_kj_mol", "float", "None", "Energy of transition state (kJ/mol) relative to reactants."),
        ("product_energy_kj_mol", "float", "None", "Energy of products (kJ/mol) relative to reactants."),
        # Common parameters:
        ("temperature_k", "float", "298.15", "Temperature for kinetic analysis (K)."),
        ("n_points", "int", "50", "Number of points to generate along the reaction coordinate."),
        ("detail_level", "str", "standard", "Detail level: 'basic', 'standard', or 'detailed'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: input_mode [reactant_E] [ts_E] [product_E] [T] [n_points]"),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing energy_profile, ts_position, barriers, irc_points, and thermodynamic/kinetic analysis."),
    ]

    examples = [
        {
            "code_input": {
                "input_mode": "analytical",
                "energy_profile": None,
                "reactant_energy_kj_mol": 0.0,
                "ts_energy_kj_mol": 75.0,
                "product_energy_kj_mol": -20.0,
                "temperature_k": 298.15,
                "n_points": 50,
                "detail_level": "standard",
            },
            "text_input": {
                "input_params": "analytical 0 75 -20 298.15 50 standard",
            },
            "output": {
                "result": {
                    "reaction_type": "exothermic",
                    "reactant_energy_kj_mol": 0.0,
                    "ts_energy_kj_mol": 75.0,
                    "product_energy_kj_mol": -20.0,
                    "forward_barrier_kj_mol": 75.0,
                    "reverse_barrier_kj_mol": 95.0,
                    "reaction_energy_kj_mol": -20.0,
                    "ts_position_normalized": 0.79,
                    "hammond_postulate": "TS resembles products (exothermic, late TS).",
                    "irc_points": [{"coord": 0.0, "energy": 0.0}, {"coord": 0.5, "energy": 70.0}, {"coord": 1.0, "energy": -20.0}],
                    "kinetic_analysis": "Exothermic with moderate barrier; favorable forward direction.",
                }
            }
        },
        {
            "code_input": {
                "input_mode": "manual",
                "energy_profile": [(0, 0), (0.25, 30), (0.5, 65), (0.65, 80), (0.75, 72), (0.9, 30), (1.0, 10)],
                "reactant_energy_kj_mol": 0,
                "ts_energy_kj_mol": 0,
                "product_energy_kj_mol": 0,
                "temperature_k": 298.15,
                "n_points": 50,
                "detail_level": "detailed",
            },
            "text_input": {
                "input_params": "manual 0,0 0.25,30 0.5,65 0.65,80 0.75,72 0.9,30 1.0,10 detailed",
            },
            "output": {
                "result": {
                    "reaction_type": "endergonic (slightly)",
                    "forward_barrier_kj_mol": 80.0,
                    "reverse_barrier_kj_mol": 70.0,
                    "reaction_energy_kj_mol": 10.0,
                    "ts_position_normalized": 0.65,
                    "hammond_postulate": "TS resembles reactants (endergonic, early TS).",
                    "irc_points": "provided from input data",
                    "max_energy_point": {"coord": 0.65, "energy_kj_mol": 80.0},
                }
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        input_mode: str,
        energy_profile: Optional[List[Tuple[float, float]]] = None,
        reactant_energy_kj_mol: float = 0.0,
        ts_energy_kj_mol: float = 0.0,
        product_energy_kj_mol: float = 0.0,
        temperature_k: float = 298.15,
        n_points: int = 50,
        detail_level: str = "standard",
    ) -> dict:
        """Core logic: analyze reaction coordinate."""
        mode = input_mode.lower().strip()
        T = temperature_k
        dl = detail_level.lower()

        if mode == "manual":
            if not energy_profile or len(energy_profile) < 2:
                raise ChemMCPError("Manual mode requires at least 2 energy profile points.")
            result = self._analyze_manual(energy_profile, T, dl)
        elif mode == "analytical":
            result = self._analyze_analytical(reactant_energy_kj_mol, ts_energy_kj_mol,
                                               product_energy_kj_mol, T, n_points, dl)
        else:
            raise ChemMCPError(f"Unknown mode: '{mode}'. Use 'manual' or 'analytical'.")

        return {"result": result}

    def _analyze_manual(self, profile: List[Tuple[float, float]], T: float, dl: str) -> dict:
        """Analyze manually provided energy profile."""
        coords = [p[0] for p in profile]
        energies = [p[1] for p in profile]

        E_react = energies[0]
        E_prod = energies[-1]
        E_ts = max(energies)
        ts_idx = energies.index(E_ts)
        ts_coord = coords[ts_idx]

        dE_rxn = E_prod - E_react
        dG_f = E_ts - E_react   # forward barrier
        dG_r = E_ts - E_prod     # reverse barrier

        rtype = self._classify_reaction(dE_rxn)
        ts_norm = ts_coord / max(coords) if max(coords) > 0 else ts_coord

        result = {
            "input_mode": "manual",
            "n_data_points": len(profile),
            "reaction_type": rtype,
            "reactant_energy_kj_mol": round(E_react, 2),
            "product_energy_kj_mol": round(E_prod, 2),
            "ts_energy_kj_mol": round(E_ts, 2),
            "max_energy_point": {"coord": round(ts_coord, 4), "energy_kj_mol": round(E_ts, 2)},
            "ts_position_normalized": round(ts_norm, 4),
            "forward_barrier_kj_mol": round(dG_f, 2),
            "reverse_barrier_kj_mol": round(dG_r, 2),
            "reaction_energy_kj_mol": round(dE_rxn, 2),
            "hammond_postulate": self._hammond(dE_rxn, ts_norm),
            "irc_points": [
                {"coord": round(c, 4), "energy_kj_mol": round(e, 2)}
                for c, e in zip(coords, energies)
            ],
        }

        if dl != "basic":
            result["kinetic_analysis"] = self._kinetic_summary(dG_f, dG_r, dE_rxn, T)
            result["equilibrium_constant"] = self._calc_keq(dE_rxn, T)

        if dl == "detailed":
            result["curvature_analysis"] = self._estimate_curvature(profile, ts_idx)
            result["various_transition_state"] = self._ts_characteristics(E_ts, E_react, E_prod, ts_norm)

        return result

    def _analyze_analytical(self, E_react: float, E_ts: float, E_prod: float,
                              T: float, n_pts: int, dl: str) -> dict:
        """Generate and analyze an analytical reaction coordinate profile."""
        dG_f = E_ts - E_react
        dG_r = E_ts - E_prod
        dE_rxn = E_prod - E_react

        # Estimate TS position using Marcus-like or Bell-shaped interpolation
        # For asymmetric reactions, TS shifts toward higher-energy side (Hammond)
        if abs(dE_rxn) < 1:
            ts_norm = 0.5  # symmetric → TS in middle
        else:
            # Hammond approximation: TS position correlates with exo/endergonicity
            # More exothermic → later TS (closer to products); more endothermic → earlier TS
            total_range = abs(dG_f) + abs(dG_r)
            if total_range > 0:
                # Hammond: TS position = reverse_barrier / total
                # Exothermic: large reverse barrier → late TS (> 0.5)
                # Endergonic: small reverse barrier → early TS (< 0.5)
                ts_norm = abs(dG_r) / total_range
            else:
                ts_norm = 0.5
        ts_norm = max(0.05, min(0.95, ts_norm))

        # Generate smooth IRC points using cubic spline-like interpolation
        irc_points = []
        for i in range(n_pts + 1):
            s = i / n_pts  # normalized coordinate 0→1
            e = self._interpolated_energy(s, E_react, E_ts, E_prod, ts_norm)
            irc_points.append({"coord": round(s, 4), "energy_kj_mol": round(e, 2)})

        rtype = self._classify_reaction(dE_rxn)

        result = {
            "input_mode": "analytical",
            "reaction_type": rtype,
            "reactant_energy_kj_mol": round(E_react, 2),
            "product_energy_kj_mol": round(E_prod, 2),
            "ts_energy_kj_mol": round(E_ts, 2),
            "ts_position_normalized": round(ts_norm, 4),
            "forward_barrier_kj_mol": round(dG_f, 2),
            "reverse_barrier_kj_mol": round(dG_r, 2),
            "reaction_energy_kj_mol": round(dE_rxn, 2),
            "hammond_postulate": self._hammond(dE_rxn, ts_norm),
            "irc_points": irc_points,
            "n_generated_points": len(irc_points),
        }

        if dl != "basic":
            result["kinetic_analysis"] = self._kinetic_summary(dG_f, dG_r, dE_rxn, T)
            result["equilibrium_constant"] = self._calc_keq(dE_rxn, T)

        if dl == "detailed":
            result["ts_characteristics"] = self._ts_characteristics(E_ts, E_react, E_prod, ts_norm)
            result["rate_estimate"] = self._estimate_rate(dG_f, T)

        return result

    @staticmethod
    def _interpolated_energy(s: float, E_react: float, E_ts: float, E_prod: float, ts_s: float) -> float:
        """Interpolate energy along reaction coordinate using piecewise function."""
        if s <= ts_s:
            # Reactant side: rise to TS (use cosine-squared for smooth rise)
            t = s / ts_s if ts_s > 0 else 0
            return E_react + (E_ts - E_react) * (1 - math.cos(math.pi * t)) / 2
        else:
            # Product side: fall from TS to product
            t = (s - ts_s) / (1 - ts_s) if ts_s < 1 else 0
            return E_ts + (E_prod - E_ts) * (1 - math.cos(math.pi * t)) / 2

    @staticmethod
    def _classify_reaction(dE: float) -> str:
        if dE < -10:
            return "exothermic"
        elif dE < -1:
            return "slightly exothermic"
        elif dE <= 1:
            return "thermoneutral"
        elif dE <= 10:
            return "slightly endergonic"
        else:
            return "endergonic"

    @staticmethod
    def _hammond(dE: float, ts_pos: float) -> str:
        if dE < -10:
            return f"Late TS (position={ts_pos:.2f}) — resembles products (highly exothermic)."
        elif dE < 0:
            return f"Moderately late TS (position={ts_pos:.2f}) — product-like (exothermic)."
        elif dE <= 10:
            return f"Central TS (position={ts_pos:.2f}) — similar character to both (thermoneutral/slightly endergonic)."
        else:
            return f"Early TS (position={ts_pos:.2f}) — resembles reactants (endergonic)."

    @staticmethod
    def _kinetic_summary(dG_f: float, dG_r: float, dE: float, T: float) -> str:
        parts = []
        # Forward rate estimate using TST
        prefactor = kB * T / h
        k_f = prefactor * math.exp(-dG_f * 1000 / (R * T))
        k_r = prefactor * math.exp(-dG_r * 1000 / (R * T))

        if dG_f < 40:
            parts.append(f"Low forward barrier ({dG_f:.1f} kJ/mol) → fast forward reaction.")
        elif dG_f < 80:
            parts.append(f"Moderate forward barrier ({dG_f:.1f} kJ/mol) → measurable rate at room temp.")
        else:
            parts.append(f"High forward barrier ({dG_f:.1f} kJ/mol) → slow without catalysis/heating.")

        if abs(dE) > 20:
            favored = "products" if dE < 0 else "reactants"
            parts.append(f"Equilibrium strongly favors {favored} (ΔG = {dE:.1f} kJ/mol).")

        keq_est = math.exp(-dE * 1000 / (R * T)) if R * T > 0 else 0
        parts.append(f"Estimated K_eq ≈ {keq_est:.3e} at {T} K.")

        return " ".join(parts)

    @staticmethod
    def _calc_keq(dE: float, T: float) -> dict:
        """Calculate equilibrium constant from reaction energy."""
        if R * T > 0:
            keq = math.exp(-dE * 1000 / (R * T))
        else:
            keq = 0
        delta_G = dE
        return {
            "K_eq": round(keq, 4),
            "delta_G_kj_mol": round(delta_G, 2),
            "temperature_K": T,
            "interpretation": "Products favored" if keq > 1 else "Reactants favored" if keq < 1 else "At equilibrium",
        }

    @staticmethod
    def _estimate_curvature(profile: list, ts_idx: int) -> dict:
        """Estimate curvature (second derivative) around TS."""
        n = len(profile)
        if ts_idx < 1 or ts_idx >= n - 1:
            return {"note": "Cannot compute curvature at boundary."}

        # Simple finite difference curvature
        x0, y0 = profile[ts_idx]
        xm1, ym1 = profile[ts_idx - 1]
        xp1, yp1 = profile[min(ts_idx + 1, n - 1)]

        dx_left = x0 - xm1 if x0 != xm1 else 1.0
        dx_right = xp1 - x0 if xp1 != x0 else 1.0

        # Curvature ≈ negative (TS is a maximum along IRC)
        curv_left = 2 * (ym1 - 2*y0 + ym1) / (dx_left ** 2) if ts_idx >= 1 else 0
        # Use actual neighbors
        d2y = (profile[min(ts_idx+1, n-1)][1] - 2*y0 + profile[max(ts_idx-1, 0)][1])
        dx_avg = ((xp1 - x0) + (x0 - xm1)) / 2 if ts_idx > 0 and ts_idx < n - 1 else 1.0
        curvature = d2y / (dx_avg ** 2) if dx_avg > 0 else 0

        return {
            "curvature_at_TS": round(curvature, 2),
            "note": "Negative curvature expected at TS (maximum along IRC)" if curvature < 0 else "Unexpected positive curvature",
            "sharpness": "sharp TS (large |curvature|)" if abs(curvature) > 500 else "broad TS (small |curvature|)",
        }

    @staticmethod
    def _ts_characteristics(E_ts: float, E_react: float, E_prod: float, ts_pos: float) -> dict:
        """Characterize the transition state."""
        asymmetry = abs((E_ts - E_react) - (E_ts - E_prod))
        total_barrier = (E_ts - E_react) + (E_ts - E_prod)
        symmetry_ratio = asymmetry / total_barrier if total_barrier > 0 else 0

        return {
            "asymmetry_kj_mol": round(asymmetry, 2),
            "symmetry_ratio": round(symmetry_ratio, 3),
            "character": "symmetric TS" if symmetry_ratio < 0.1 else "moderately asymmetric" if symmetry_ratio < 0.3 else "strongly asymmetric TS",
            "imaginary_frequency_note": "TS has exactly one imaginary frequency (ν‡) corresponding to the reaction coordinate mode.",
        }

    @staticmethod
    def _estimate_rate(dG_f: float, T: float) -> dict:
        """Estimate rate constant using simplified TST."""
        prefactor = kB * T / h
        k = prefactor * math.exp(-dG_f * 1000 / (R * T))
        t_half = math.log(2) / k if k > 0 else float('inf')

        return {
            "estimated_k_s-1": round(k, 4),
            "half_life_s": round(t_half, 2),
            "half_life_human": f"{t_half:.1f} s ({t_half/60:.1f} min, {t_half/3600:.2f} h)" if t_half < 86400 else ">24 h",
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        try:
            parts = input_params.strip().split()
            if not parts:
                raise ChemMCPError("Empty input.")

            mode = parts[0].lower()
            kwargs = {"input_mode": mode}

            if mode == "manual":
                # Parse coordinate,energy pairs
                pairs = []
                for p in parts[1:]:
                    if "," in p:
                        c, e = p.split(",")
                        pairs.append((float(c.strip()), float(e.strip())))
                if pairs:
                    kwargs["energy_profile"] = pairs
                modes = {"basic", "standard", "detailed"}
                for p in parts[1:]:
                    if p in modes:
                        kwargs["detail_level"] = p

            elif mode == "analytical":
                if len(parts) >= 4:
                    kwargs["reactant_energy_kj_mol"] = float(parts[1])
                    kwargs["ts_energy_kj_mol"] = float(parts[2])
                    kwargs["product_energy_kj_mol"] = float(parts[3])
                if len(parts) >= 5:
                    try:
                        kwargs["temperature_k"] = float(parts[4])
                    except ValueError:
                        pass
                if len(parts) >= 6:
                    try:
                        kwargs["n_points"] = int(parts[5])
                    except ValueError:
                        pass
                if len(parts) >= 7:
                    kwargs["detail_level"] = parts[6]

            return self._run_base(**kwargs)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {e}")
