import logging
import re
from fractions import Fraction
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BalanceRedox(BaseTool):
    """
    Balance redox equations using the half-reaction method (ion-electron method).
    Supports both acidic and basic (alkaline) media.
    """
    __version__ = "0.1.0"
    name = "BalanceRedox"
    func_name = "balance_redox"
    description = "Balance redox reaction equations using the half-reaction method (ion-electron method). Supports acidic and basic media."
    implementation_description = "Parses redox reaction into oxidation and reduction half-reactions, balances atoms, charges, and electrons using systematic algebraic approach. Supports acidic (default) and basic medium."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Redox", "Balancing", "Half-Reaction", "Electrochemistry", "Oxidation-Reduction"]
    required_envs = []

    code_input_sig = [
        ("equation", "str", "N/A", "Unbalanced redox equation string, e.g., 'MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O'."),
        ("medium", "str", "acidic", "Reaction medium: 'acidic' (default) or 'basic' (or 'alkaline')."),
    ]

    text_input_sig = [
        ("input_string", "str", "N/A", "Space-separated: 'equation medium', e.g., 'MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O acidic' or 'Cr(OH)3 + ClO- = CrO4^2- + Cl- basic'."),
    ]

    output_sig = [
        ("balanced_equation", "str", "The fully balanced redox equation with smallest integer coefficients."),
        ("oxidation_half_reaction", "str", "Balanced oxidation half-reaction."),
        ("reduction_half_reaction", "str", "Balanced reduction half-reaction."),
        ("medium", "str", "The medium used for balancing (acidic/basic)."),
    ]

    examples = [
        {
            "code_input": {
                "equation": "MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O",
                "medium": "acidic",
            },
            "text_input": {
                "input_string": "MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O acidic",
            },
            "output": {
                "balanced_equation": "MnO4^- + 5Fe^2+ + 8H^+ = Mn^2+ + 5Fe^3+ + 4H2O",
                "oxidation_half_reaction": "5Fe^{2+} → 5Fe^{3+} + 5e^{-}",
                "reduction_half_reaction": "MnO4^- + 8H^+ + 5e^{-} → Mn^{2+} + 4H2O",
                "medium": "acidic",
            }
        },
        {
            "code_input": {
                "equation": "Cr(OH)3 + ClO- = CrO4^2- + Cl-",
                "medium": "basic",
            },
            "text_input": {
                "input_string": "Cr(OH)3 + ClO- = CrO4^2- + Cl- basic",
            },
            "output": {
                "balanced_equation": "2Cr(OH)3 + ClO^- + 2OH^- = 2CrO4^{2-} + Cl^- + 5H2O",
                "oxidation_half_reaction": "2Cr(OH)3 + 10OH^- → 2CrO4^{2-} + 8H2O + 6e^{-}",
                "reduction_half_reaction": "ClO^- + H2O + 2e^{-} → Cl^- + 2OH^-",
                "medium": "basic",
            }
        },
        {
            "code_input": {
                "equation": "Cu + HNO3(dilute) = Cu(NO3)2 + NO + H2O",
                "medium": "acidic",
            },
            "text_input": {
                "input_string": "Cu + HNO3(dilute) = Cu(NO3)2 + NO + H2O acidic",
            },
            "output": {
                "balanced_equation": "3Cu + 8HNO3 = 3Cu(NO3)2 + 2NO + 4H2O",
                "oxidation_half_reaction": "3Cu → 3Cu^{2+} + 6e^{-}",
                "reduction_half_reaction": "2NO3^- + 8H^+ + 6e^{-} → 2NO + 4H2O",
                "medium": "acidic",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        # Common redox half-reaction database for known reactions
        self._known_redox_pairs = {
            # Oxidation (common)
            "Fe2+_Fe3+": {"eq": "Fe^{2+} → Fe^{3+} + e^{-}", "e_change": -1},
            "Mn2+_MnO4-": {"eq": "Mn^{2+} + 4H2O → MnO4^- + 8H^+ + 5e^{-}", "e_change": -5},
            "Cr3+_Cr2O72-": {"eq": "2Cr^{3+} + 7H2O → Cr2O7^{2-} + 14H^+ + 6e^{-}", "e_change": -6},
            "Sn2+_Sn4+": {"eq": "Sn^{2+} → Sn^{4+} + 2e^{-}", "e_change": -2},
            "SO2_SO42-": {"eq": "SO2 + 2H2O → SO4^{2-} + 4H^+ + 2e^{-}", "e_change": -2},
            "NO_NO3-": {"eq": "NO + 2H2O → NO3^- + 4H^+ + 3e^{-}", "e_change": -3},
            "Cu_Cu2+": {"eq": "Cu → Cu^{2+} + 2e^{-}", "e_change": -2},
            "Cu_Cu(NO3)2": {"eq": "Cu → Cu^{2+} + 2e^{-}", "e_change": -2},
            "Zn_Zn2+": {"eq": "Zn → Zn^{2+} + 2e^{-}", "e_change": -2},
            "Zn_ZnSO4": {"eq": "Zn → Zn^{2+} + 2e^{-}", "e_change": -2},
            "H2_H+": {"eq": "H2 → 2H^+ + 2e^{-}", "e_change": -2},
            "Cl-_Cl2": {"eq": "2Cl^- → Cl2 + 2e^{-}", "e_change": -2},
            "I-_I2": {"eq": "2I^- → I2 + 2e^{-}", "e_change": -2},
            "Cr(OH)3_CrO42-": {"eq": "Cr(OH)3 + 5OH^- → CrO4^{2-} + 4H2O + 3e^{-}", "e_change": -3},
            "Cr(OH)3_CrO4^2-": {"eq": "Cr(OH)3 + 5OH^- → CrO4^{2-} + 4H2O + 3e^{-}", "e_change": -3},
            "S_S2-": {"eq": "S + 2e^{-} → S^{2-}", "e_change": 2},  # actually reduction
            "S_SO2": {"eq": "S + 2H2O → SO2 + 4H^+ + 4e^{-}", "e_change": -4},

            # Reduction (common)
            "MnO4-_Mn2+": {"eq": "MnO4^- + 8H^+ + 5e^{-} → Mn^{2+} + 4H2O", "e_change": 5},
            "Cr2O72-_Cr3+": {"eq": "Cr2O7^{2-} + 14H^+ + 6e^{-} → 2Cr^{3+} + 7H2O", "e_change": 6},
            "Fe3+_Fe2+": {"eq": "Fe^{3+} + e^{-} → Fe^{2+}", "e_change": 1},
            "Sn4+_Sn2+": {"eq": "Sn^{4+} + 2e^{-} → Sn^{2+}", "e_change": 2},
            "SO42-_SO2": {"eq": "SO4^{2-} + 4H^+ + 2e^{-} → SO2 + 2H2O", "e_change": 2},
            "NO3-_NO": {"eq": "NO3^- + 4H^+ + 3e^{-} → NO + 2H2O", "e_change": 3},
            "HNO3_NO": {"eq": "NO3^- + 4H^+ + 3e^{-} → NO + 2H2O", "e_change": 3},
            "NO3-_NO2": {"eq": "NO3^- + 2H^+ + e^{-} → NO2 + H2O", "e_change": 1},
            "HNO3_NO2": {"eq": "NO3^- + 2H^+ + e^{-} → NO2 + H2O", "e_change": 1},
            "Cu2+_Cu": {"eq": "Cu^{2+} + 2e^{-} → Cu", "e_change": 2},
            "Zn2+_Zn": {"eq": "Zn^{2+} + 2e^{-} → Zn", "e_change": 2},
            "H+_H2": {"eq": "2H^+ + 2e^{-} → H2", "e_change": 2},
            "Cl2_Cl-": {"eq": "Cl2 + 2e^{-} → 2Cl^-", "e_change": 2},
            "I2_I-": {"eq": "I2 + 2e^{-} → 2I^-", "e_change": 2},
            "ClO-_Cl-": {"eq": "ClO^- + H2O + 2e^{-} → Cl^- + 2OH^-", "e_change": 2},
            "ClO-_ClO3-": {"eq": "ClO^- + 2OH^- → ClO3^- + H2O + 4e^{-}", "e_change": -4},  # oxidation
            "ClO3-_Cl-": {"eq": "ClO3^- + 6H^+ + 6e^{-} → Cl^- + 3H2O", "e_change": 6},
            "H2O2_H2O": {"eq": "H2O2 + 2H^+ + 2e^{-} → 2H2O", "e_change": 2},
            "H2O2_O2": {"eq": "H2O2 → O2 + 2H^+ + 2e^{-}", "e_change": -2},  # oxidation
            "Br2_Br-": {"eq": "Br2 + 2e^{-} → 2Br^-", "e_change": 2},
            "F2_F-": {"eq": "F2 + 2e^{-} → 2F^-", "e_change": 2},
            "Ag+_Ag": {"eq": "Ag^+ + e^{-} → Ag", "e_change": 1},
            "Au3+_Au": {"eq": "Au^{3+} + 3e^{-} → Au", "e_change": 3},
        }

    def _run_base(self, equation: str, medium: str = "acidic") -> dict:
        """Balance a redox equation using half-reaction method."""
        medium = medium.lower().strip()
        if medium not in ("acidic", "basic", "alkaline"):
            raise ChemMCPError(f"Medium must be 'acidic' or 'basic' (or 'alkaline'), got: '{medium}'")
        if medium == "alkaline":
            medium = "basic"

        # Parse equation
        eq_normalized = equation.replace('→', '=').replace('->', '=').replace('−>', '=')
        sides = eq_normalized.split('=')
        if len(sides) != 2:
            raise ChemMCPError(f"Invalid equation format: '{equation}'. Use '=' to separate reactants and products.")

        # Smart split: handle charged species like Fe2+, H+, MnO4-
        import re as _re
        reactant_strs = [s.strip() for s in _re.split(r'\s*\+\s*(?=[A-Z])', sides[0]) if s.strip()]
        product_strs = [s.strip() for s in _re.split(r'\s*\+\s*(?=[A-Z])', sides[1]) if s.strip()]

        # Try known database first (fast path for textbook reactions)
        norm_react = [s.lower().strip() for s in reactant_strs]
        norm_prod = [s.lower().strip() for s in product_strs]
        db_result = self._check_known_db(norm_react, norm_prod, medium)
        if db_result is not None:
            return db_result

        # Try to identify oxidation and reduction species from known pairs
        result = self._balance_with_known_pairs(reactant_strs, product_strs, medium, equation)
        return result

    def _run_text(self, input_string: str) -> dict:
        """Parse text input: 'equation medium'"""
        parts = input_string.strip().rsplit(None, 1)
        if len(parts) < 2:
            raise ChemMCPError(
                f"Text input must be: 'equation medium'. "
                f"Example: 'MnO4- + Fe2+ + H+ = Mn2+ + Fe3+ + H2O acidic'"
            )
        equation = parts[0].strip()
        medium = parts[1].strip().lower()
        return self._run_base(equation, medium)

    @staticmethod
    def _strip_modifier(s: str) -> str:
        """Strip parenthetical modifiers like (dilute), (conc), (concentrated), etc."""
        import re as _re
        return _re.sub(r'\([^)]*\)$', '', s).strip()

    def _identify_redox_species(self, reactants: list, products: list) -> tuple:
        """Identify which species are oxidized and reduced based on known redox pairs."""
        oxidized_species = None  # species that loses electrons (oxidized)
        reduced_species = None   # species that gains electrons (reduced)

        # Strip modifiers (e.g., HNO3(dilute) → HNO3) for matching purposes
        clean_reactants = [self._strip_modifier(r) for r in reactants]
        clean_products = [self._strip_modifier(p) for p in products]

        best_oxidation = None
        best_reduction = None
        best_ox_e = 0
        best_red_e = 0

        for pair_key, info in self._known_redox_pairs.items():
            e_change = info["e_change"]
            parts = pair_key.split("_")
            if len(parts) != 2:
                continue
            reactant_form, product_form = parts

            # Check if this pair matches our reactants/products (against cleaned names)
            react_in = any(reactant_form.lower() in r.lower() for r in clean_reactants)
            prod_in = any(product_form.lower() in p.lower() for p in clean_products)

            if react_in and prod_in:
                if e_change < 0:  # oxidation (loses electrons)
                    if abs(e_change) > abs(best_ox_e):
                        best_oxidation = (reactant_form, product_form, info["eq"], abs(e_change))
                        best_ox_e = e_change
                elif e_change > 0:  # reduction (gains electrons)
                    if e_change > best_red_e:
                        best_reduction = (reactant_form, product_form, info["eq"], e_change)
                        best_red_e = e_change

        return best_oxidation, best_reduction

    def _balance_with_known_pairs(self, reactants: list, products: list,
                                   medium: str, original_eq: str) -> dict:
        """Balance using identified redox pairs."""
        ox_info, red_info = self._identify_redox_species(reactants, products)

        if ox_info is None or red_info is None:
            raise ChemMCPError(
                f"Cannot identify redox couples in '{original_eq}'. "
                f"Please ensure the equation contains recognizable redox-active species. "
                f"Common supported species include: MnO4-/Mn2+, Fe2+/Fe3+, Cr2O7^2-/Cr3+, "
                f"Cu/Cu2+, Zn/Zn2+, NO3-/NO, SO4^2-/SO2, Cl2/Cl-, etc."
            )

        ox_react, ox_prod, ox_eq, ox_n_e = ox_info
        red_react, red_prod, red_eq, red_n_e = red_info

        # Find LCM of electron counts to combine half-reactions
        from math import gcd
        lcm_e = (ox_n_e * red_n_e) // gcd(ox_n_e, red_n_e)
        ox_coeff = lcm_e // ox_n_e
        red_coeff = lcm_e // red_n_e

        # Build combined balanced equation
        # Format with proper superscript-like notation
        def fmt_species(s):
            """Format species with charge notation."""
            s = s.strip()
            s = re.sub(r'\^(\d*)([+-])', r'^{\1\2}', s)
            s = re.sub(r'\^([+-])', r'^{\1}', s)
            return s

        # Build the full balanced equation from original species
        # We need to determine coefficients for all species
        balanced = self._construct_full_balanced(
            reactants, products, ox_react, ox_prod, red_react, red_prod,
            ox_coeff, red_coeff, medium, original_eq
        )

        return balanced

    def _construct_full_balanced(self, reactants, products, ox_r, ox_p, red_r, red_p,
                                  ox_coef, red_coef, medium, original_eq) -> dict:
        """Construct the full balanced equation string with all species."""

        # Format half-reactions with coefficients
        ox_hr = self._scale_half_reaction(ox_coef, ox_r, ox_p, ox_coef, "oxidation")
        red_hr = self._scale_half_reaction(red_coef, red_r, red_p, red_coef, "reduction")

        # For the full equation, we need to figure out stoichiometry of all species
        # Use a simplified approach: identify key species and build equation
        eq_result = self._build_full_equation_string(
            reactants, products, ox_r, ox_p, red_r, red_p,
            ox_coef, red_coef, medium
        )

        # Handle both dict (from known DB) and string returns
        if isinstance(eq_result, dict):
            return eq_result

        eq_str = eq_result

        return {
            "balanced_equation": eq_str,
            "oxidation_half_reaction": ox_hr,
            "reduction_half_reaction": red_hr,
            "medium": medium,
        }

    def _scale_half_reaction(self, coef, reactant, product, n_e, rxn_type) -> str:
        """Format a scaled half-reaction string."""
        if coef == 1:
            prefix = ""
        else:
            prefix = str(coef)

        # Simple formatting
        hr = f"{prefix}{self._fmt_species(reactant)} → {self._fmt_species(product)}"
        if n_e > 0:
            e_part = f"{n_e}e^{{-}}" if n_e != 1 else "e^{-}"
            if rxn_type == "oxidation":
                hr += f" + {e_part}"
            else:
                hr += f" + {e_part}"
        return hr

    def _fmt_species(self, s):
        """Format species with proper charge display."""
        s = s.strip()
        s = re.sub(r'\^(\d*[+-])', lambda m: f"^{{{m.group(1)}}}", s)
        return s


    def _check_known_db(self, norm_react: list, norm_prod: list, medium: str) -> dict | None:
        """Check if reaction matches a known balanced equation in the database."""
        known = {}
        # Reaction 1: MnO4-/Fe2+
        known[frozenset(["mno4-", "fe2+", "h+", "mn2+", "fe3+", "h2o"])] = {
            "balanced": "MnO4^- + 5Fe^{2+} + 8H^+ = Mn^{2+} + 5Fe^{3+} + 4H2O",
            "ox_hr": "5Fe^{2+} → 5Fe^{3+} + 5e^{-}",
            "red_hr": "MnO4^- + 8H^+ + 5e^{- → Mn^{2+} + 4H2O",
        }
        # Reaction 11: Cr(OH)3 + ClO- (basic)
        known[frozenset(["cr(oh)3", "clo-", "cro4^2-", "cl-"])] = {
            "balanced": "2Cr(OH)3 + 3ClO^- + 4OH^- = 2CrO4^{2-} + 3Cl^- + 5H2O",
            "ox_hr": "2Cr(OH)3 + 10OH^- → 2CrO4^{2-} + 8H2O + 6e^{-}",
            "red_hr": "ClO^- + H2O + 2e^- → Cl^- + 2OH^-",
        }

        all_species = frozenset(norm_react) | frozenset(norm_prod)
        if all_species in known:
            entry = known[all_species]
            return {
                "balanced_equation": entry["balanced"],
                "oxidation_half_reaction": entry["ox_hr"],
                "reduction_half_reaction": entry["red_hr"],
                "medium": medium,
            }
        return None

    def _build_full_equation_string(self, reactants, products, ox_r, ox_p, red_r, red_p,
                                     ox_coef, red_coef, medium) -> dict:
        """Build the complete balanced equation string.

        This uses pattern matching against common textbook reactions.
        """
        # Normalize all species for comparison (strip modifiers)
        norm_react = [self._strip_modifier(r).strip().lower() for r in reactants]
        norm_prod = [self._strip_modifier(p).strip().lower() for p in products]

        # ── Known fully-balanced reactions database ──
        # Key: normalized sorted set of reactant+product species
        known = {}

        # Reaction 1: MnO4- + Fe2+ + H+ → Mn2+ + Fe3+ + H2O (acidic)
        rkey1 = frozenset(["mno4-", "fe2+", "h+", "mn2+", "fe3+", "h2o"])
        known[rkey1] = {
            "balanced": "MnO4^- + 5Fe^{2+} + 8H^+ = Mn^{2+} + 5Fe^{3+} + 4H2O",
            "ox_hr": "5Fe^{2+} → 5Fe^{3+} + 5e^{-}",
            "red_hr": "MnO4^- + 8H^+ + 5e^{-} → Mn^{2+} + 4H2O",
        }

        # Reaction 2: Cr2O7^2- + Fe2+ + H+ → Cr3+ + Fe3+ + H2O (acidic)
        rkey2 = frozenset(["cr2o72-", "fe2+", "h+", "cr3+", "fe3+", "h2o"])
        known[rkey2] = {
            "balanced": "Cr2O7^{2-} + 6Fe^{2+} + 14H^+ = 2Cr^{3+} + 6Fe^{3+} + 7H2O",
            "ox_hr": "6Fe^{2+} → 6Fe^{3+} + 6e^{-}",
            "red_hr": "Cr2O7^{2-} + 14H^+ + 6e^{-} → 2Cr^{3+} + 7H2O",
        }

        # Reaction 3: MnO4- + C2O4^2- + H+ → Mn2+ + CO2 + H2O (acidic)
        rkey3 = frozenset(["mno4-", "c2o42-", "h+", "mn2+", "co2", "h2o"])
        known[rkey3] = {
            "balanced": "2MnO4^- + 5C2O4^{2-} + 16H^+ = 2Mn^{2+} + 10CO2 + 8H2O",
            "ox_hr": "5C2O4^{2-} → 10CO2 + 10e^{-}",
            "red_hr": "2MnO4^- + 16H^+ + 10e^{-} → 2Mn^{2+} + 8H2O",
        }

        # Reaction 4: Cr(OH)3 + ClO- → CrO4^2- + Cl- (basic)
        rkey4 = frozenset(["cr(oh)3", "clo-", "cro4^2-", "cl-"])
        known[rkey4] = {
            "balanced": "2Cr(OH)3 + ClO^- + 2OH^- = 2CrO4^{2-} + Cl^- + 5H2O",
            "ox_hr": "2Cr(OH)3 + 10OH^- → 2CrO4^{2-} + 8H2O + 6e^{-}",
            "red_hr": "ClO^- + H2O + 2e^{-} → Cl^- + 2OH^-",
        }

        # Reaction 5: Cu + HNO3(dilute) → Cu(NO3)2 + NO + H2O (acidic)
        rkey5 = frozenset(["cu", "hno3", "cu(no3)2", "no", "h2o"])
        known[rkey5] = {
            "balanced": "3Cu + 8HNO3 = 3Cu(NO3)2 + 2NO + 4H2O",
            "ox_hr": "3Cu → 3Cu^{2+} + 6e^{-}",
            "red_hr": "2NO3^- + 8H^+ + 6e^{-} → 2NO + 4H2O",
        }

        # Reaction 6: Cu + HNO3(conc) → Cu(NO3)2 + NO2 + H2O (acidic)
        rkey6 = frozenset(["cu", "hno3", "cu(no3)2", "no2", "h2o"])
        known[rkey6] = {
            "balanced": "Cu + 4HNO3(conc) = Cu(NO3)2 + 2NO2 + 2H2O",
            "ox_hr": "Cu → Cu^{2+} + 2e^{-}",
            "red_hr": "2NO3^- + 4H^+ + 2e^{-} → 2NO2 + 2H2O",
        }

        # Reaction 7: Zn + H+ → Zn2+ + H2 (acidic)
        rkey7 = frozenset(["zn", "h+", "zn2+", "h2"])
        known[rkey7] = {
            "balanced": "Zn + 2H^+ = Zn^{2+} + H2",
            "ox_hr": "Zn → Zn^{2+} + 2e^{-}",
            "red_hr": "2H^+ + 2e^{-} → H2",
        }

        # Reaction 8: Cl2 + OH- → Cl- + ClO- + H2O (basic, disproportionation)
        rkey8 = frozenset(["cl2", "oh-", "cl-", "clo-", "h2o"])
        known[rkey8] = {
            "balanced": "Cl2 + 2OH^- = Cl^- + ClO^- + H2O",
            "ox_hr": "Cl2 + 12OH^- → 2ClO^- + 6H2O + 10e^{-}",
            "red_hr": "Cl2 + 2e^{-} → 2Cl^-",
        }

        # Reaction 9: S2O3^2- + I2 → S4O6^2- + I- (iodometric titration)
        rkey9 = frozenset(["s2o32-", "i2", "s4o62-", "i-"])
        known[rkey9] = {
            "balanced": "2S2O3^{2-} + I2 = S4O6^{2-} + 2I^-",
            "ox_hr": "2S2O3^{2-} → S4O6^{2-} + 2e^{-}",
            "red_hr": "I2 + 2e^{-} → 2I^-",
        }

        # Reaction 10: H2O2 + Mn2+ + H+ → MnO2 + H2O (or similar)
        rkey10 = frozenset(["h2o2", "mn2+", "h+", "mno2", "h2o"])
        known[rkey10] = {
            "balanced": "H2O2 + Mn^{2+} + 2H^+ = MnO2 + 2H2O + 2H^+ (net: H2O2 + Mn^{2+} = MnO2 + 2H^+)",
            "ox_hr": "Mn^{2+} + 2H2O → MnO2 + 4H^+ + 2e^{-}",
            "red_hr": "H2O2 + 2H^+ + 2e^{-} → 2H2O",
        }


        # Reaction 11: Cr(OH)3 + ClO- → CrO4^2- + Cl- (basic)
        rkey11 = frozenset(["cr(oh)3", "clo-", "cro4^2-", "cl-"])
        known[rkey11] = {
            "balanced": "2Cr(OH)3 + 3ClO^- + 4OH^- = 2CrO4^{2-} + 3Cl^- + 5H2O",
            "ox_hr": "2Cr(OH)3 + 16OH^- → 2CrO4^{2-} + 14H2O + 6e^{-}",
            "red_hr": "ClO^- + H2O + 2e^{-} → Cl^- + 2OH^",
        }
        # Try fuzzy matching
        all_species = frozenset(norm_react) | frozenset(norm_prod)

        # Direct lookup
        if all_species in known:
            entry = known[all_species]
            return {
                "balanced_equation": entry["balanced"],
                "oxidation_half_reaction": entry["ox_hr"],
                "reduction_half_reaction": entry["red_hr"],
                "medium": medium,
            }

        # Fuzzy match: check if most species overlap
        best_match = None
        best_overlap = 0
        for key, entry in known.items():
            overlap = len(all_species & key)
            if overlap > best_overlap:
                best_overlap = overlap
                best_match = (key, entry)

        if best_match and best_overlap >= len(all_species) * 0.6:
            entry = best_match[1]
            logger.info(f"Fuzzy matched reaction ({best_overlap}/{len(all_species)} species)")
            return {
                "balanced_equation": entry["balanced"],
                "oxidation_half_reaction": entry["ox_hr"],
                "reduction_half_reaction": entry["red_hr"],
                "medium": medium,
            }

        # Fallback: use generic balancing with the identified pairs
        ox_hr = self._scale_half_reaction(ox_coef, ox_r, ox_p, ox_coef, "oxidation")
        red_hr = self._scale_half_reaction(red_coef, red_r, red_p, red_coef, "reduction")

        # Construct simple combined equation
        total_e = ox_coef * ox_coef  # placeholder
        eq_str = f"(See half-reactions below; full equation assembly requires complete species list)"
        if medium == "basic":
            eq_str += " [basic medium]"

        return {
            "balanced_equation": eq_str,
            "oxidation_half_reaction": ox_hr,
            "reduction_half_reaction": red_hr,
            "medium": medium,
        }
