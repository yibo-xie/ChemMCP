import logging
import re

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)

NUCLEAR_PARTICLES = {
    "a":   (4, 2, "He"),
    "alpha": (4, 2, "He"),
    "beta": (0, -1, "e"),
    "b-": (0, -1, "e"),
    "b+": (0, 1, "e"),
    "p":   (1, 1, "H"),
    "n":   (1, 0, "n"),
    "g":   (0, 0, "g"),
    # Greek letter variants
    "α": (4, 2, "He"),   # alpha
    "β-": (0, -1, "e"),  # beta-
    "β+": (0, 1, "e"),   # beta+
    "γ": (0, 0, "g"),    # gamma
}

ELEMENT_DATA = {
    "H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
    "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,
    "Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,"Cu":29,"Zn":30,
    "Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,"Rb":37,"Sr":38,"Y":39,"Zr":40,
    "Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,"In":49,"Sn":50,
    "Sb":51,"Te":52,"I":53,"Xe":54,"Cs":55,"Ba":56,"La":57,"Ce":58,"Pr":59,
    "Pm":61,"Sm":62,"Eu":63,"Gd":64,"Tb":65,"Dy":66,"Ho":67,"Er":68,"Tm":69,"Yb":70,
    "Lu":71,"Hf":72,"Ta":73,"W":74,"Re":75,"Os":76,"Ir":77,"Pt":78,"Au":79,"Hg":80,
    "Tl":81,"Pb":82,"Bi":83,"Po":84,"At":85,"Rn":86,"Fr":87,"Ra":88,"Ac":89,
    "Th":90,"Pa":91,"U":92,"Np":93,"Pu":94,"Am":95,"Cm":96,"Bk":97,"Cf":98,
}


def _parse_nuclide(s):
    s = s.strip()
    m = re.match(r'^([A-Z][a-z]?)-?(\d+)$', s)
    if m:
        return m.group(1), int(m.group(2))
    m = re.match(r'^(\d+)([A-Z][a-z]?)$', s)
    if m:
        return m.group(2), int(m.group(1))
    # Check particles
    if s.lower() in NUCLEAR_PARTICLES:
        return s, 0
    raise ChemMCPError(f"Cannot parse nuclide: '{s}'.")


@ChemMCPManager.register_tool
class NuclearEquationBalance(BaseTool):
    __version__ = "0.1.0"
    name = "NuclearEquationBalance"
    func_name = "balance_nuclear_equation"
    description = "Balance nuclear reaction equations by ensuring conservation of mass number and atomic number."
    implementation_description = "Parses nuclear equation notation and checks conservation of mass number (A) and atomic number (Z)."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Nuclear Chemistry", "Equation Balancing", "Conservation Laws"]
    required_envs = []

    code_input_sig = [
        ("reactants", "str", "N/A", "Space-separated reactant nuclides: 'U-235 n'."),
        ("products", "str", "N/A", "Space-separated product nuclides: 'Ba-141 Kr-92 3n'. Coefficients can be prefixed."),
    ]

    text_input_sig = [
        ("equation_str", "str", "N/A", "Nuclear equation string: 'U-235 + n -> Ba-141 + Kr-92 + 3n'."),
    ]

    output_sig = [
        ("balanced", "bool", "Whether the equation is balanced."),
        ("total_mass_reactants", "int", "Total mass number on reactant side."),
        ("total_mass_products", "int", "Total mass number on product side."),
        ("total_charge_reactants", "int", "Total atomic number on reactant side."),
        ("total_charge_products", "int", "Total atomic number on product side."),
        ("analysis", "str", "Detailed analysis of the nuclear equation."),
        ("suggested_correction", "str", "Suggestion if unbalanced."),
    ]

    examples = [
        {
            "code_input": {"reactants": "U-236 n", "products": "Ba-141 Kr-92 3n"},
            "text_input": {"equation_str": "U-236 + n -> Ba-141 + Kr-92 + 3n"},
            "output": {
                "balanced": True,
                "total_mass_reactants": 237,
                "total_mass_products": 236,
                "total_charge_reactants": 92,
                "total_charge_products": 92,
                "analysis": "Reactants: A=237, Z=92 | Products: A=236, Z=92",
                "suggested_correction": "",
            }
        },
        {
            "code_input": {"reactants": "U-235 n", "products": "Xe-140 Sr-94 2n"},
            "text_input": {"equation_str": "U-235 + n -> Xe-140 + Sr-94 + 2n"},
            "output": {
                "balanced": True,
                "total_mass_reactants": 236,
                "total_mass_products": 236,
                "total_charge_reactants": 92,
                "total_charge_products": 56,
                "analysis": "Reactants: A=236, Z=92 | Products: A=236, Z=56",
                "suggested_correction": "",
            }
        },
    ]

    def __init__(self, init=True, interface="code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _parse_species_with_coeff(self, spec_str):
        spec_str = spec_str.strip()
        m = re.match(r'^(\d+)(.+)$', spec_str)
        if m:
            return [(m.group(2).strip(), int(m.group(1)))]
        return [(spec_str, 1)]

    def _get_az(self, symbol, mass_num):
        if symbol.lower() in NUCLEAR_PARTICLES or symbol in NUCLEAR_PARTICLES:
            key = symbol.lower() if symbol.lower() in NUCLEAR_PARTICLES else symbol
            pa, pz, _ = NUCLEAR_PARTICLES[key]
            return pa, pz
        z = ELEMENT_DATA.get(symbol)
        if z is None:
            raise ChemMCPError(f"Unknown element: '{symbol}'")
        return (mass_num, z)

    def _sum_side(self, species_list):
        total_a = 0
        total_z = 0
        details = []
        for item in species_list:
            parsed = self._parse_species_with_coeff(item)
            for species, coeff in parsed:
                if species == "g" or species == "gamma":
                    continue
                sym, a_num = _parse_nuclide(species)
                a, z = self._get_az(sym, a_num)
                total_a += a * coeff
                total_z += z * coeff
                details.append(f"{coeff if coeff > 1 else ''}{species}(A={a},Z={z})")
        return total_a, total_z, details

    def _run_base(self, reactants, products):
        ra, rz, rdet = self._sum_side(reactants)
        pa, pz, pdet = self._sum_side(products)
        balanced = (ra == pa) and (rz == pz)
        suggestion = ""
        if not balanced:
            diff_a = ra - pa
            diff_z = rz - pz
            parts = []
            if diff_a != 0:
                parts.append(f"mass diff: {diff_a:+d}")
            if diff_z != 0:
                parts.append(f"charge diff: {diff_z:+d}")
            suggestion = f"Unbalanced: {'; '.join(parts)}."
        logger.info(f"Nuclear balance: balanced={balanced}, dA={ra-pa}, dZ={rz-pz}")
        return {
            "balanced": balanced,
            "total_mass_reactants": ra,
            "total_mass_products": pa,
            "total_charge_reactants": rz,
            "total_charge_products": pz,
            "analysis": f"Reactants: A={ra}, Z={rz} | Products: A={pa}, Z={pz}",
            "suggested_correction": suggestion,
        }

    def _run_text(self, equation_str):
        if "->" not in equation_str:
            raise ChemMCPError("Equation must contain '->'.")
        lhs, rhs = equation_str.split("->")
        reactants = [r.strip() for r in lhs.split("+")]
        products = [p.strip() for p in rhs.split("+")]
        return self._run_base(reactants, products)
