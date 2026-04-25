import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MainGroupTrends(BaseTool):
    """
    主族元素性质递变规律查询工具（s区和p区元素）。
    覆盖原子半径、电离能、电负性、原子体积、金属性/非金属性、
    氧化物酸碱性、氢化物稳定性等周期性递变规律，包括对角线规则。
    """
    __version__ = "0.1.0"
    name = "MainGroupTrends"
    func_name = "get_main_group_trends"
    description = "Query main group (s-block and p-block) periodic trends: atomic radius, ionization energy, electronegativity, metallic character, oxide acidity/basicity, hydride stability, diagonal relationships, and group-specific patterns."
    implementation_description = "Comprehensive database of periodic trends for main group elements (Groups 1-2, 13-18). Covers horizontal and vertical trends with explanations based on effective nuclear charge, shielding, orbital penetration, and quantum mechanical effects. Includes diagonal rule (Li-Mg, Be-Al, B-Si) and key anomalies explained."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Periodic Trends", "Main Group", "s-Block", "p-Block", "Periodic Table", "Diagonal Rule"]
    required_envs = []

    code_input_sig = [
        ("trend_name", "str", "all", "Name of trend to query: 'atomic_radius', 'ionization_energy', 'electronegativity', 'metallic_character', 'oxide_acidity', 'hydride_stability', 'diagonal_rule', 'group_summary', or 'all'."),
        ("group_number", "int", "0", "Optional specific group number (1-18) to focus on. 0 = all groups."),
        ("period", "int", "0", "Optional period number (1-7). 0 = all periods."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'trend_name [group] [period]'. Example: 'ionization_energy' or 'diagonal_rule' or 'oxide_acidity 15'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing trend data, explanation, examples, and exceptions."),
    ]

    examples = [
        {
            "code_input": {"trend_name": "atomic_radius", "group_number": 0, "period": 0},
            "text_input": {"input_params": "atomic_radius"},
            "output": {"result": {"trend": "atomic radius trend data", "across_period": "decreases", "down_group": "increases"}}
        },
        {
            "code_input": {"trend_name": "diagonal_rule", "group_number": 0, "period": 0},
            "text_input": {"input_params": "diagonal_rule"},
            "output": {"result": {"pairs": ["Li-Mg, Be-Al, B-Si"], "explanation": "diagonal relationship due to similar charge density"}}
        },
    ]

    TRENDS_DB = {
        "atomic_radius": {
            "definition": "Half the distance between nuclei of two adjacent atoms of the same element (metallic radius for metals, covalent/van der Waals radius for nonmetals).",
            "across_period_LEFT_to_RIGHT": "DECREASES significantly",
            "reason_across": "Increasing nuclear charge (+1 proton per element) pulls electrons closer; electrons added to same principal shell → poor mutual shielding → Z_eff increases steadily from ~2.5 to ~10+ across a period.",
            "down_group_TOP_to_BOTTOM": "INCREASES",
            "reason_down": "Each new electron shell is farther from nucleus; shielding by inner electrons outweighs increased nuclear charge.",
            "key_data_pm": {
                "Period_2": "Li(152) > Be(112) > B(85) > C(77) > N(75) > O(73) > F(72) > Ne(71)",
                "Period_3": "Na(186) > Mg(160) > Al(143) > Si(118) > P(110) > S(103) > Cl(99) > Ar(98)",
                "Group_1": "Li(152) < Na(186) < K(227) < Rb(248) < Cs(265)",
                "Group_17": "F(42) < Cl(79) < Br(94) < I(114)",
            },
            "anomalies": [
                "Oxygen radius slightly larger than expected (electron-electron repulsion in small p-shell)",
                "Gallium (Ga, 122 pm) has similar radius to Al (143 pm) despite being one period lower — d-block contraction effect",
            ],
        },

        "ionization_energy": {
            "definition": "Minimum energy required to remove the most loosely bound electron from a gaseous atom at ground state. Units: kJ/mol (or eV).",
            "general_trend": "INCREASES across a period (left→right), DECREASES down a group (top→bottom)",
            "reason_across": "Z_eff increases → valence electrons held more tightly",
            "reason_down": "Outer electrons farther from nucleus + more inner-shell shielding → easier to remove",
            "periodic_dips": {
                "Be(899) < B(801)": "B loses a p-electron (higher energy, shielded by s²); Be removes stable s² pair",
                "N(1402) < O(1314)": "O loses paired p-electron (e⁻-e⁻ repulsion makes removal easier); N has stable half-filled p³",
                "Mg(738) < Al(578)": "Al loses p-electron vs Mg removing stable s²",
                "P(1012) < S(1000)": "S loses paired e⁻ from p⁴; P has stable half-filled p³",
            },
            "special_values_kjmol": {
                "HIGHEST": "He (2372), He⁺ still holds e⁻ very tightly due to tiny size + no shielding",
                "LOWEST (non-radioactive)": "Cs (375.7), Fr (~380)",
                "notable": "F (1681) — high but NOT highest because of e⁻-e⁻ repulsion in small 2p shell",
            },
            "successive_IEs": "IE₁ << IE₂ << IE₃ ... (each successive electron harder to remove; large jumps indicate valence shell completion)",
        },

        "electronegativity": {
            "definition": "Ability of an atom in a chemical bond to attract shared electrons toward itself. Most common scale: Pauling (dimensionless, 0.7–4.0 range).",
            "trend": "INCREASES left→right across period, DECREASES top→bottom down group",
            "pauling_values": {
                "H": 2.20,
                "Period_2": "Li(0.98) < Be(1.57) < B(2.04) < C(2.55) < N(3.04) < O(3.44) < F(3.98)",
                "Period_3": "Na(0.93) < Mg(1.31) < Al(1.61) < Si(1.90) < P(2.19) < S(2.58) < Cl(3.16)",
                "Group_1": "Li(0.98) > Na(0.93) > K(0.82) > Rb(0.82) > Cs(0.79)",
                "Group_17": "F(3.98) > Cl(3.16) > Br(2.96) > I(2.66) > At(2.20)",
            },
            "most_electronegative": "Fluorine (3.98) — defines the scale (χ_F was set to 4.0 originally, refined to 3.98)",
            "least_electronegative": "Fr/Cs (~0.79) — most electropositive elements",
            "applications": "Predict bond polarity (Δχ > 1.7 ≈ ionic; Δχ < 1.7 ≈ covalent), oxidation direction in redox reactions, acid strength trends",
        },

        "metallic_character": {
            "definition": "Extent to which an element exhibits properties typical of metals: luster, malleability, good electrical/thermal conductivity, tendency to lose electrons (form cations), basic oxides.",
            "trend": "INCREASES down a group, DECREASES across a period (left→right)",
            "metalloid_line": "Stair-step line from B-Si-As-Te-At separates metals (left/below) from nonmetals (right/above); elements on the line are metalloids",
            "metallic_character_by_region": {
                "s-block (Grp 1-2)": "Strong metals (except H)",
                "d-block": "All metals (except some say Hg is borderline)",
                "p-block_left (Grp 13-14)": "Metals (Al, Ga, In, Tl, Sn, Pb) + metalloids (B, Si, Ge)",
                "p-block_right (Grp 15-18)": "Nonmetals (N, O, F, Cl + noble gases) + metalloids (As, Sb, Te) + some metals (Bi, Po)",
            },
            "physical_correlates": "Metallic character correlates with: low IE, low EN, large atomic radius, thermal/electrical conductivity, basic oxide character",
        },

        "oxide_acidity_basicity": {
            "definition": "Character of the highest oxide formed by each main group element when combined with oxygen.",
            "trend_across_period": "Basic → Amphoteric → Acidic (moving right)",
            "trend_down_group": "Acidity decreases / Basicity increases (within same group)",
            "examples_by_group": {
                "Group_1": "Li₂O, Na₂O, K₂O, Rb₂O, Cs₂O — ALL strongly BASIC (alkali metal oxides form strong bases MOH with water)",
                "Group_2": "BeO (amphoteric!) < MgO (weakly basic) < CaO < SrO < BaO (strongly basic)",
                "Group_13": "Al₂O₃ (amphoteric) — reacts with both acids AND bases; B₂O₃ (acidic); Ga₂O₃, In₂O₃, Tl₂O₃ amphoteric/basic",
                "Group_14": "CO₂, SiO₂ (ACIDIC — acidic oxides); GeO₂ (amphoteric); SnO, PbO (amphoteric/basic)",
                "Group_15": "N₂O₅, P₄O₁₀ (ACIDIC — form oxyacids); As₄O₆ (weakly acidic); Sb₂O₃ (amphoteric); Bi₂O₃ (basic! — inert pair effect)",
                "Group_16": "SO₃ (acidic); SeO₃ (acidic); TeO₃ (less acidic/amphoteric); PoO₂ (basic/amphoteric)",
                "Group_17": "Cl₂O₇, SO₃-type (ACIDIC — form strong oxyacids like HClO₄)",
            },
            "rule_of_thumb": "Oxide character reflects the element's position: metallic elements → basic oxides; nonmetallic → acidic oxides; intermediate → amphoteric",
        },

        "hydride_stability_and_type": {
            "classification": {
                "ionic_hydrides": "Group 1 & 2 (s-block): M-H ionic/salt-like; high mp; react violently with water to give H₂; strong bases (H⁻ is powerful base/proton acceptor)",
                "covalent_hydrides": "Group 13-17 (p-block): molecular compounds; volatile; acidity increases down Group 14-17",
                "metallic_interstitial": "d-block: H occupies interstitial sites in metal lattice; non-stoichiometric; variable composition",
            },
            "thermal_stability_trend": "DECREASES down each group for covalent hydrides",
            "stability_order_Grp16": "H₂O >> H₂S > H₂Se > H₂Te > H₂Po (bond energy: O-H 463 > S-H 363 > Se-H 276 > Te-H ~238 kJ/mol)",
            "acidity_trend_Grp16": "H₂O (neutral, pKa=14) < H₂S (weak, pKa≈7) < H₂Se (pKa≈3.9) < H₂Te (pKa≈2.6) — weaker H-X bond = stronger acid",
            "reducing_power_trend": "Increases down group (weaker E-H bond breaks more easily, donates H⁻/H• more readily)",
            "special_cases": [
                "HF is weak acid (strong H-bonding) while other HX are strong acids",
                "H₂O has anomalously high bp (100°C vs H₂S -60°C) due to hydrogen bonding",
                "NH₃ is basic (lone pair on N); PH₃ much less basic; AsH₃/SbH₃ essentially neutral",
                "SnH₄, PbH₄ highly unstable (inert pair effect — heavy p-block reluctant to use s-electrons in bonding)",
            ],
        },

        "diagonal_rule": {
            "definition": "Elements that are diagonally adjacent in the periodic table often show similar properties due to similar charge-to-size ratios (charge density).",
            "pairs": [
                {
                    "pair": "Li — Mg",
                    "similarities": [
                        "Both form nitrides directly with N₂ (Li₃N, Mg₃N₂)",
                        "Both form normal oxides (Li₂O, MgO) rather than peroxides/superoxides",
                        "Carbonates decompose on heating (Li₂CO₃ → Li₂O + CO₂; MgCO₃ similarly)",
                        "Both form covalent organometallic compounds (RLi, RMgX Grignard reagents)",
                        "Hydroxides/carbonates have limited solubility in water (unlike other G1/G2)",
                        "Both have diagonal relationship in their salts (e.g., both form hydrated salts)",
                    ]
                },
                {
                    "pair": "Be — Al",
                    "similarities": [
                        "Both oxides are AMPHOTERIC (react with acids AND bases)",
                        "Both form covalent compounds (high polarizing power from small size/high charge)",
                        "Both hydroxides are amphoteric (Be(OH)₂, Al(OH)₃)",
                        "Both form carbides that release methane on hydrolysis (Be₂C, Al₄C₃)",
                        "Both are rendered passive by concentrated HNO₃ (protective oxide layer)",
                        "Chlorides act as Lewis acids / bridge polymers (BeCl₂ chain structure, AlCl₃ dimeric)",
                    ]
                },
                {
                    "pair": "B — Si",
                    "similarities": [
                        "Both form semiconducting elements (B: semiconductor, Si: quintessential semiconductor)",
                        "Both oxides are acidic (B₂O₃, SiO₂) and form glassy solids",
                        "Both form covalent halides that hydrolyze readily (BCl₃, SiCl₄)",
                        "Both form hydrides that are volatile and spontaneously flammable (B₂H₆ boranes, SiH₄ silane)",
                        "Both form oxygen-containing acids (H₃BO₃ boric acid, H₄SiO₄ silicic acid — both WEAK acids)",
                        "Both exist in allotropes (multiple structural forms)",
                    ]
                },
            ],
            "underlying_reason": "Moving one step right-down in the periodic table simultaneously INCREASES atomic radius (like going down a group) and INCREASES nuclear charge/Z_eff (like going across a period). These two effects partially cancel, giving similar charge density → similar chemistry.",
        },

        "group_summaries": {
            "Group_1_Alkali_Metals": "ns¹ configuration; strongest reducing agents; form +1 ions exclusively; oxides/hydroxides strongly basic; reactivity increases dramatically down group; all react with water (rate: Li << Na < K < Rb < Cs); stored under oil (except Li which floats on oil — stored in petroleum jelly)",
            "Group_2_Alkaline_Earth": "ns² configuration; form +2 ions; harder than G1; Be shows anomalous amphoteric/covalent behavior; reactivity with water increases down group (Be none, Mg slow/steam, Ca moderate, Sr/Ba vigorous); sulfate solubility DECREASES down group (opposite of most trends); hydroxide solubility INCREASES down group",
            "Group_13_Boron_Group": "ns²np¹; B is metalloid (forms covalent network/molecular compounds); Al and below are metals; +3 dominant oxidation state; Tl also shows +1 (inert pair effect); B forms unique multicenter-bonded compounds (boranes); Al is protected by oxide layer (amphoteric)",
            "Group_14_Carbon_Group": "ns²np²; contains nonmetal (C), metalloids (Si, Ge), metals (Sn, Pb); +4 and +2 oxidation states; +2 becomes more stable down group (inert pair effect: Pb(II) > Pb(IV)); catenation ability greatest for C (C-C bonds strong); Si also shows significant catenation (silicones); CO/CO₂ (molecular) vs SiO₂ (network solid) — fundamental difference",
            "Group_15_Pnictogens": "ns²np³; +5, +3, -3 oxidation states; metallic character increases down group (N, P: nonmetals; As, Sb: metalloids; Bi: metal); basicity of hydrides NH₃ >> PH₃ > AsH₃; acidity of oxoacids increases with oxidation state (HNO₃ > HNO₂); Bi(V) is strong oxidizer (inert pair makes Bi(III) stable); N₂ extremely inert (triple bond 941 kJ/mol); P can form P≡P only with extreme stabilization",
            "Group_16_Calcogens": "ns²np⁴; -2, +2, +4, +6 oxidation states; O is VERY different (small, high EN, forms H-bonds, paramagnetic O₂); Po is radioactive metal; trend: nonmetal (O, S) → metalloid (Se, Te) → metal (Po); H-bonding important for O compounds (H₂O anomaly); allotropy: O₂/O₃, S (many allotropes: S₈ rings, chains, etc.)",
            "Group_17_Halogens": "ns²np⁵; -1 dominant; positive oxidation states up to +7 (Cl, Br, I); reactivity F₂ >> Cl₂ > Br₂ > I₂; hydrogen halide acidity HF(weak) < HCl < HBr < HI(strongest); interhalogen compounds exist (ClF, BrF₃, IF₅, IF₇); oxidizing power X₂ decreases down group; reducing power X⁻ increases down group",
            "Group_18_Noble_Gases": "ns²np⁶ full valence; historically called 'inert gases'; He has highest IE (2372 kJ/mol); Xe and Kr form compounds with F and O (most notably XeF₂, XeF₄, XeF₆, XeO₃); He used in diving gas (prevents nitrogen narcosis), cryogenics (lowest bp), arc welding; Radon is radioactive health hazard indoors; Ne used in lighting/signs",
        },
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, trend_name: str, group_number: int = 0, period: int = 0) -> dict:
        trend = trend_name.lower().strip().replace(" ", "_")

        if trend == "all":
            return {"result": self.TRENDS_DB}

        if trend not in self.TRENDS_DB:
            valid = list(self.TRENDS_DB.keys())
            raise ChemMCPError(f"Trend '{trend_name}' not found. Options: {valid}")

        data = self.TRENDS_DB[trend]

        # Filter by group if specified
        if group_number > 0:
            grp_key = f"Group_{group_number}"
            filtered = {k: v for k, v in data.items() if grp_key.lower() in k.lower() or grp_key in str(v)}
            if filtered:
                return {"result": {"trend": trend_name, **filtered}}
            return {"result": {"trend": trend_name, "note": f"No group-specific data for Group {group_number} in this trend.", **data}}

        return {"result": {"trend": trend_name, **data}}

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            trend = parts[0] if parts else "all"
            grp = int(parts[1]) if len(parts) > 1 else 0
            per = int(parts[2]) if len(parts) > 2 else 0
            return self._run_base(trend, grp, per)
        except ValueError:
            raise ChemMCPError(f"Invalid format. Use: 'trend_name [group_number] [period]'")
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}")
