import logging
import re
from fractions import Fraction
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BalanceEquation(BaseTool):
    """
    Balance chemical equations (ordinary reactions) using linear algebra (matrix method).
    Input: unbalanced equation string like "H2+O2=H2O"
    Output: balanced equation with smallest integer coefficients.
    """
    __version__ = "0.1.0"
    name = "BalanceEquation"
    func_name = "balance_equation"
    description = "Balance chemical equations for ordinary reactions using algebraic matrix method."
    implementation_description = "Parses chemical formulas into element matrices, constructs a homogeneous linear system, and solves via Gaussian elimination with fraction arithmetic to find smallest positive integer coefficients."
    oss_dependencies = []
    services_and_software = []
    categories = ["Reaction"]
    tags = ["Balancing", "Stoichiometry", "Linear Algebra", "Chemical Equation"]
    required_envs = []

    code_input_sig = [
        ("equation", "str", "N/A", "Unbalanced chemical equation string, e.g., 'H2+O2=H2O', 'Fe+O2=Fe2O3', 'C6H12O6+O2=CO2+H2O'. Use '=' or '→' as separator."),
    ]

    text_input_sig = [
        ("equation", "str", "N/A", "Unbalanced chemical equation string."),
    ]

    output_sig = [
        ("balanced_equation", "str", "The balanced equation with smallest positive integer coefficients."),
        ("coefficients", "str", "List of coefficients in order [reactants..., products...]."),
    ]

    examples = [
        {
            "code_input": {"equation": "H2+O2=H2O"},
            "text_input": {"equation": "H2+O2=H2O"},
            "output": {
                "balanced_equation": "2H2 + O2 = 2H2O",
                "coefficients": "[2, 1, 2]",
            }
        },
        {
            "code_input": {"equation": "Fe+O2=Fe2O3"},
            "text_input": {"equation": "Fe+O2=Fe2O3"},
            "output": {
                "balanced_equation": "4Fe + 3O2 = 2Fe2O3",
                "coefficients": "[4, 3, 2]",
            }
        },
        {
            "code_input": {"equation": "C6H12O6+O2=CO2+H2O"},
            "text_input": {"equation": "C6H12O6+O2=CO2+H2O"},
            "output": {
                "balanced_equation": "C6H12O6 + 6O2 = 6CO2 + 6H2O",
                "coefficients": "[1, 6, 6, 6]",
            }
        },
        {
            "code_input": {"equation": "Al+HCl=AlCl3+H2"},
            "text_input": {"equation": "Al+HCl=AlCl3+H2"},
            "output": {
                "balanced_equation": "2Al + 6HCl = 2AlCl3 + 3H2",
                "coefficients": "[2, 6, 2, 3]",
            }
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    # ── Formula parsing helpers ──

    @staticmethod
    def _parse_formula(formula: str) -> dict:
        """Parse a chemical formula into {element: count} dict."""
        elements = {}
        # Match element symbols and optional counts
        pattern = r'([A-Z][a-z]?)(\d*)'
        for match in re.finditer(pattern, formula):
            elem = match.group(1)
            count_str = match.group(2)
            count = int(count_str) if count_str else 1
            if elem:
                elements[elem] = elements.get(elem, 0) + count
        return elements

    @staticmethod
    def _parse_compound(compound: str) -> dict:
        """Parse a compound formula handling simple groups (parentheses)."""
        elements = {}
        i = 0
        while i < len(compound):
            if compound[i] == '(':
                # Find matching closing paren
                depth = 1
                j = i + 1
                while j < len(compound) and depth > 0:
                    if compound[j] == '(':
                        depth += 1
                    elif compound[j] == ')':
                        depth -= 1
                    j += 1
                group_content = compound[i+1:j-1]
                # Get multiplier after closing paren
                k = j
                mult_str = ""
                while k < len(compound) and compound[k].isdigit():
                    mult_str += compound[k]
                    k += 1
                mult = int(mult_str) if mult_str else 1
                group_elems = BalanceEquation._parse_formula(group_content)
                for elem, cnt in group_elems.items():
                    elements[elem] = elements.get(elem, 0) + cnt * mult
                i = k
            elif compound[i].isupper():
                # Element symbol
                elem = compound[i]
                i += 1
                if i < len(compound) and compound[i].islower():
                    elem += compound[i]
                    i += 1
                # Get count
                count_str = ""
                while i < len(compound) and compound[i].isdigit():
                    count_str += compound[i]
                    i += 1
                count = int(count_str) if count_str else 1
                elements[elem] = elements.get(elem, 0) + count
            else:
                i += 1
        return elements

    @staticmethod
    def _split_equation(equation: str) -> tuple:
        """Split equation into reactants and products lists."""
        # Normalize separator
        eq = equation.replace('→', '=').replace('->', '=').replace('−>', '=')
        sides = eq.split('=')
        if len(sides) != 2:
            raise ChemMCPError(
                f"Invalid equation format: '{equation}'. "
                f"Use '=' or '→' to separate reactants and products."
            )
        reactants = [c.strip() for c in sides[0].split('+') if c.strip()]
        products = [c.strip() for c in sides[1].split('+') if c.strip()]
        return reactants, products

    @staticmethod
    def _build_matrix(compounds: list) -> tuple:
        """Build element × compound matrix. Returns (elements_list, matrix)."""
        all_elements = set()
        parsed = []
        for comp in compounds:
            elems = BalanceEquation._parse_compound(comp)
            parsed.append(elems)
            all_elements.update(elems.keys())
        elements = sorted(all_elements)

        matrix = []
        for elem in elements:
            row = []
            for elems_dict in parsed:
                row.append(Fraction(elems_dict.get(elem, 0)))
            matrix.append(row)

        return elements, matrix

    @staticmethod
    def _gauss_solve(matrix: list, n_compounds: int) -> list:
        """Solve homogeneous system via Gauss-Jordan elimination. Returns coefficient list."""
        m = [list(row) for row in matrix]  # deep copy
        n_rows = len(m)
        n_cols = n_compounds

        pivot_cols = []
        pivot_row = 0

        for col in range(n_cols):
            # Find pivot
            found = False
            for row in range(pivot_row, n_rows):
                if m[row][col] != 0:
                    # Swap rows
                    m[pivot_row], m[row] = m[row], m[pivot_row]
                    found = True
                    break

            if not found:
                continue

            pivot_cols.append(col)

            # Scale pivot row
            scale = m[pivot_row][col]
            for j in range(n_cols):
                m[pivot_row][j] /= scale

            # Eliminate other rows
            for row in range(n_rows):
                if row != pivot_row and m[row][col] != 0:
                    factor = m[row][col]
                    for j in range(n_cols):
                        m[row][j] -= factor * m[pivot_row][j]

            pivot_row += 1
            if pivot_row >= n_rows:
                break

        # Find free variable (last column not in pivot cols)
        free_col = None
        for col in range(n_cols - 1, -1, -1):
            if col not in pivot_cols:
                free_col = col
                break

        if free_col is None:
            raise ChemMCPError("Could not find free variable; equation may be trivial or unsolvable.")

        # Set free variable to 1, compute others
        coeffs = [Fraction(0)] * n_cols
        coeffs[free_col] = Fraction(1)

        # Back-substitute in reverse pivot order
        for idx in range(len(pivot_cols) - 1, -1, -1):
            pc = pivot_cols[idx]
            val = Fraction(0)
            for j in range(pc + 1, n_cols):
                val -= m[idx][j] * coeffs[j]
            coeffs[pc] = val

        # Reactants should have positive coefficients, products negative → flip signs
        # Convention: first half are reactants (positive), second half are products
        # We need all coefficients positive; negate if needed
        # Check if any coefficient is negative
        has_neg = any(c < 0 for c in coeffs)
        if has_neg:
            coeffs = [-c for c in coeffs]

        # Convert to smallest positive integers
        # Find LCM of denominators
        from math import gcd
        def lcm(a, b):
            return abs(a * b) // gcd(a, b) if a and b else abs(a or b)

        denominators = [abs(c.denominator) for c in coeffs]
        common_denom = 1
        for d in denominators:
            common_denom = lcm(common_denom, d)

        int_coeffs = [int(c * common_denom) for c in coeffs]

        # Reduce by GCD
        overall_gcd = 0
        for c in int_coeffs:
            overall_gcd = gcd(overall_gcd, abs(c))
        if overall_gcd > 1:
            int_coeffs = [c // overall_gcd for c in int_coeffs]

        return int_coeffs

    def _run_base(self, equation: str) -> dict:
        """Balance a chemical equation."""
        reactants, products = self._split_equation(equation)
        all_compounds = reactants + products
        n_reactants = len(reactants)

        if len(all_compounds) < 2:
            raise ChemMCPError("Equation must have at least one reactant and one product.")

        elements, matrix = self._build_matrix(all_compounds)

        # Products get negative sign in matrix (conservation: reactants = products)
        n_products = len(products)
        for row_idx in range(len(matrix)):
            for col_idx in range(n_reactants, n_reactants + n_products):
                matrix[row_idx][col_idx] = -matrix[row_idx][col_idx]

        try:
            coeffs = self._gauss_solve(matrix, len(all_compounds))
        except Exception as e:
            raise ChemMCPError(f"Failed to balance equation '{equation}': {e}")

        # Build balanced equation string
        parts = []
        for compound, coeff in zip(all_compounds, coeffs):
            if coeff == 1:
                parts.append(compound)
            else:
                parts.append(f"{coeff}{compound}")

        balanced_eq = " + ".join(parts[:n_reactants]) + " = " + " + ".join(parts[n_reactants:])

        logger.info(f"Balanced: {equation} → {balanced_eq}")
        return {
            "balanced_equation": balanced_eq,
            "coefficients": str(coeffs),
        }

    def _run_text(self, equation: str) -> dict:
        return self._run_base(equation)
