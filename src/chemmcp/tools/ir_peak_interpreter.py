import logging
from typing import List, Optional, Union

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class IrPeakInterpreter(BaseTool):
    """
    IR 光谱特征峰解释工具。
    根据波数（cm⁻¹）归属官能团。
    覆盖范围：4000-400 cm⁻¹
    """
    __version__ = "0.1.0"
    name = "IrPeakInterpreter"
    func_name = "interpret_ir_peak"
    description = "Interpret infrared (IR) absorption peaks by wavenumber. Assigns functional groups to given peak positions."
    implementation_description = "Uses a comprehensive database of IR absorption frequencies (4000-400 cm⁻¹) covering O-H, N-H, C-H, C≡C, C≡N, C=O, C=C, NO₂, and fingerprint region absorptions. Supports single or multiple peak input."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Spectroscopy", "IR", "Functional Group", "Vibrational Spectroscopy", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("wavenumber", "float", "N/A", "Wavenumber in cm⁻¹ to interpret. Can also pass a list of floats."),
        ("tolerance", "float", "15.0", "Tolerance in cm⁻¹ for matching (default ±15)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated wavenumbers. Example: '3400 2950 1710 1600'"),
    ]

    output_sig = [
        ("interpretations", "list", "List of interpretations for each peak, including functional group assignment, vibration type, intensity, and notes."),
    ]

    examples = [
        {
            "code_input": {"wavenumber": 1715.0, "tolerance": 15.0},
            "text_input": {"input_params": "1715"},
            "output": {
                "interpretations": [{
                    "wavenumber": 1715,
                    "assignment": "Saturated aliphatic ketone C=O stretch",
                    "region": "Carbonyl (1750-1650)",
                    "intensity": "strong",
                }]
            },
        },
        {
            "code_input": {"wavenumber": [3400, 2950, 1710, 1600, 1450], "tolerance": 20},
            "text_input": {"input_params": "3400 2950 1710 1600 1450"},
            "output": {
                "interpretations": "list of assignments for each peak",
            },
        },
    ]

    # ========== IR ABSORPTION DATABASE ==========
    # Format: (center_wavenumber, range_low, range_high, assignment, intensity, notes)
    _IR_DATABASE = [
        # === O-H stretching region ===
        (3600, 3650, 3590, "Free O-H stretch (alcohol/phenol)", "sharp", "Sharp; H-bonded shifts to lower frequency and broadens"),
        (3450, 3550, 3200, "Hydrogen-bonded O-H stretch (alcohol/phenol)", "broad, strong", "Broad due to hydrogen bonding; position depends on H-bond strength"),
        (3300, 3400, 3200, "O-H stretch (carboxylic acid dimer)", "very broad, very strong", "Very broad (300-600 cm⁻¹ wide); centers around 3000"),
        (3620, 3640, 3600, "Free O-H stretch (intramolecular H-bond)", "sharp", "Internal H-bond gives intermediate sharpness"),

        # === N-H stretching region ===
        (3500, 3550, 3450, "N-H stretch (primary amine, asymmetric)", "medium", "Two bands: asym ~3450, sym ~3370"),
        (3370, 3420, 3320, "N-H stretch (primary amine, symmetric)", "medium", "Second band of primary amine"),
        (3460, 3500, 3420, "N-H stretch (secondary amine)", "weak-medium", "Single band; can be absent if weakly polar"),
        (3350, 3400, 3250, "N-H stretch (amide, primary)", "medium", "Amide I + II bands; N-H appears here"),
        (3300, 3350, 3180, "N-H stretch (amide, secondary)", "strong", "One N-H band; sharper than alcohol O-H"),
        (3330, 3380, 3280, "N-H stretch (pyrrole / indole)", "weak", "N-H in aromatic heterocycle; sharp"),
        (3040, 3100, 3030, "=C-H stretch (aromatic)", "medium", "Multiple fine peaks just above 3000"),
        (3020, 3080, 3000, "=C-H stretch (alkene)", "medium", "Just above 3000 cm⁻¹"),

        # === C-H stretching region ===
        (2960, 2980, 2940, "C-H stretch (alkane, asymmetric CH₃)", "strong", "Asymmetric stretch of methyl group"),
        (2870, 2900, 2850, "C-H stretch (alkane, symmetric CH₃)", "medium", "Symmetric stretch of methyl group"),
        (2930, 2950, 2900, "C-H stretch (alkane, asymmetric CH₂)", "strong", "Asymmetric stretch of methylene"),
        (2850, 2870, 2830, "C-H stretch (alkane, symmetric CH₂)", "medium", "Symmetric stretch of methylene"),
        (2890, 2910, 2870, "C-H stretch (aldehyde C-H)", "medium", "Fermi doublet with ~2820; characteristic of -CHO"),
        (2820, 2840, 2800, "C-H stretch (aldehyde C-H, Fermi resonance)", "medium", "Paired with ~2890; diagnostic for aldehydes"),
        (2780, 2830, 2680, "C-H stretch (methoxy / OCH₃)", "medium-strong", "Characteristic of ethers with O-CH₃ group"),

        # === Triple bond region ===
        (2260, 2260, 2220, "C≡N stretch (nitrile)", "medium", "Sharp; conjugation lowers to ~2220"),
        (2190, 2220, 2180, "C≡N stretch (conjugated nitrile)", "medium", "Lowered by conjugation with double bonds/aromatics"),
        (2140, 2160, 2120, "C≡C stretch (terminal alkyne)", "weak-medium", "Terminal alkynes show this band; ≡C-H at ~3300"),
        (2100, 2140, 2060, "Cumulene C=C=C stretch", "variable", "Allene/ketene type structures"),
        (2100, 2150, 2050, "-N₃⁺ (azide) asymmetric stretch", "strong", "Azide compounds; also ~1340 cm⁻¹ symmetrical"),
        (2130, 2170, 2090, "-N=C=O (isocyanate) stretch", "strong", "Isocyanates (R-N=C=O)"),

        # === Carbonyl region (most important!) ===
        (1745, 1760, 1730, "C=O stretch (saturated ester)", "strong", "Ester carbonyl; α,β-unsaturation lowers by ~20"),
        (1740, 1755, 1725, "C=O stretch (aldehyde)", "strong", "Aliphatic aldehyde; conjugation lowers"),
        (1740, 1755, 1725, "C=O stretch (γ-lactone, 5-membered ring)", "strong", "Ring strain raises frequency vs acyclic ester"),
        (1735, 1745, 1720, "C=O stretch (saturated carboxylic acid)", "strong", "Often broad due to H-bonding"),
        (1725, 1740, 1710, "C=O stretch (saturated ketone)", "strong", "The 'standard' ketone reference point (~1715)"),
        (1720, 1735, 1705, "C=O stretch (α,β-unsaturated aldehyde)", "strong", "Lowered by ~20 from saturated by conjugation"),
        (1715, 1728, 1700, "C=O stretch (α,β-unsaturated ketone)", "strong", "Conjugation lowers from ~1715 to ~1685"),
        (1710, 1720, 1695, "C=O stretch (carboxylic acid dimer)", "strong", "Dimeric acids; broader than ketones"),
        (1705, 1720, 1690, "C=O stretch (α-diketone / β-ketoester)", "strong", "Two C=O groups may split into two bands"),
        (1700, 1715, 1685, "C=O stretch (α,β-unsaturated ester)", "strong", "Lowered by conjugation"),
        (1695, 1710, 1680, "C=O stretch (conjugated aldehyde)", "strong", "Further lowered by extended conjugation"),
        (1690, 1705, 1675, "C=O stretch (α,β-unsaturated acid)", "strong", "Acid + unsaturation = lower frequency"),
        (1685, 1700, 1670, "C=O stretch (aryl ketone / diaryl ketone)", "strong", "Aromatic conjugation lowers frequency"),
        (1680, 1695, 1665, "C=O stretch (α,β,γ,δ-unsaturated ketone)", "strong", "Extended conjugation further lowers"),
        (1675, 1685, 1660, "C=O stretch (amide, primary)", "strong", "Amide I band; N-H bend (Amide II) at ~1550"),
        (1670, 1680, 1655, "C=O stretch (amide, secondary)", "strong", "Secondary amide; often broader"),
        (1660, 1675, 1645, "C=O stretch (amide, tertiary)", "strong", "Tertiary amide; no N-H"),
        (1655, 1670, 1635, "C=O stretch (β-lactam, 4-membered ring)", "strong", "Ring strain raises frequency significantly"),
        (1660, 1680, 1640, "C=O stretch (urea / urethane)", "strong", "N-C(=O)-N type structure"),
        (1650, 1670, 1630, "C=O stretch (hindered cyclic ketone)", "strong", "Steric effects can raise or lower"),
        (1745, 1770, 1720, "C=O stretch (anhydride, two bands)", "two strong", "Anhydrides show TWO C=O bands: ~1820 & ~1760"),
        (1820, 1840, 1790, "C=O stretch (anhydride, higher freq band)", "strong", "Coupled asymmetric C=O stretch of anhydride"),
        (1765, 1780, 1750, "C=O stretch (anhydride, lower freq band)", "strong", "Coupled symmetric C=O stretch of anhydride"),
        (1815, 1830, 1795, "C=O stretch (acid chloride / acyl halide)", "very strong", "Inductive effect of halogen raises frequency"),
        (1805, 1820, 1790, "C=O stretch (β-lactone, 4-ring ester)", "strong", "High ring strain → high frequency"),
        (1740, 1760, 1720, "C=O stretch (δ-lactone, 6-ring ester)", "strong", "Moderate ring strain effect"),
        (1730, 1745, 1715, "C=O stretch (formate ester HCOOR)", "strong", "Formates slightly higher than other esters"),

        # === C=C and C=N region ===
        (1650, 1670, 1620, "C=C stretch (alkene)", "variable", "Stronger when symmetry is reduced; conjugation lowers"),
        (1645, 1655, 1635, "C=C stretch (terminal alkene =CH₂)", "medium-strong", "Terminal alkenes stronger than internal"),
        (1655, 1670, 1640, "C=C stretch (cis-disubstituted alkene)", "medium", "cis usually weaker than trans"),
        (1670, 1680, 1660, "C=C stretch (trans-disubstituted alkene)", "weak/absent", "Trans alkenes often IR-inactive or very weak"),
        (1620, 1640, 1590, "C=C stretch (conjugated diene)", "medium", "Conjugation splits into multiple bands"),
        (1600, 1610, 1580, "C=C stretch (aromatic ring, quadrant stretch)", "variable", "One of the 4 characteristic aromatic bands"),
        (1580, 1595, 1570, "C=C stretch (aromatic ring, semicircle stretch)", "variable", "Second aromatic band; sensitive to substitution"),
        (1510, 1520, 1490, "C=C stretch (aromatic ring)", "medium", "Third aromatic band; very characteristic"),
        (1450, 1470, 1430, "C=C stretch (aromatic ring)", "medium", "Fourth aromatic band; overlaps with CH₂ bend"),
        (1690, 1700, 1640, "C=N stretch (imine / oxime / Schiff base)", "variable", "Lower than C=O; oximes near 1660"),
        (1660, 1680, 1640, "C=N stretch (oxime C=N)", "medium", "Oxime functional group"),
        (1650, 1670, 1630, "C=N stretch (azine / hydrazone)", "medium", "N=N-C= systems"),

        # === Fingerprint region (important functional groups) ===
        (1465, 1480, 1440, "CH₂ bending (scissoring)", "medium", "Methylene scissoring; nearly all organic compounds"),
        (1450, 1465, 1435, "CH₃ asymmetric bending (asym deformation)", "medium", "Methyl umbrella mode"),
        (1385, 1395, 1370, "CH₃ symmetric bending (umbrella mode)", "medium-strong", "Characteristic methyl signal; splitting indicates isopropyl/t-butyl"),
        (1365, 1390, 1350, "CH₃ bending (gem-dimethyl / isopropyl split)", "doublet", "Isopropyl: equal intensity doublet ~1385/1365"),
        (1360, 1370, 1340, "CH₃ bending (tert-butyl split)", "doublet", "t-Butyl: unequal doublet, low-freq more intense"),
        (1375, 1385, 1365, "CH₃ umbrella (isolated methyl on heteroatom)", "medium", "O-CH₃, N-CH₃ often shifted"),
        (1440, 1460, 1415, "O-H bend (in-plane, carboxylic acids)", "medium", "Carboxylic acid OH deformation"),
        (1420, 1440, 1400, "C-O-H bend (alcohol)", "medium", "Associated with OH-bearing carbons"),
        (1400, 1420, 1380, "=C-H bend (in-plane, alkene/aromatic)", "medium", "Alkene/aromatic CH in-plane bends"),
        (1330, 1350, 1310, "NO₂ symmetric stretch (nitro)", "strong", "Nitro compounds: asym ~1550, sym ~1350"),
        (1555, 1570, 1540, "NO₂ asymmetric stretch (nitro)", "strong", "Nitro compounds always show both asym and sym bands"),
        (1320, 1340, 1300, "S=O asymmetric stretch (sulfone)", "strong", "Sulfones: asym ~1350, sym ~1150"),
        (1160, 1180, 1140, "S=O symmetric stretch (sulfone)", "strong", "Second sulfone SO₂ band"),
        (1310, 1330, 1290, "C-H rock (wagging, CH₂)", "medium", "Long-chain polymethylene shows regular series"),
        (1300, 1320, 1270, "C-N stretch (aromatic amine)", "medium", "Aromatic amines"),
        (1280, 1300, 1240, "C-N stretch (aliphatic amine)", "medium", "Aliphatic amines"),
        (1280, 1300, 1240, "Ar-O stretch (aryl ether / phenol C-O)", "strong", "Aromatic ether linkage"),
        (1250, 1265, 1230, "C-O stretch (aryl alkyl ether Ar-O-R)", "strong", "Phenolic ether"),
        (1240, 1260, 1210, "C-O stretch (saturated aliphatic ether)", "strong", "Aliphatic ether (R-O-R')"),
        (1220, 1240, 1200, "C-O stretch (tertiary alcohol)", "strong", "Tertiary C-OH stretch"),
        (1170, 1190, 1140, "C-O stretch (secondary alcohol)", "strong", "Secondary C-OH stretch"),
        (1060, 1080, 1030, "C-O stretch (primary alcohol)", "strong", "Primary C-OH stretch"),
        (1120, 1140, 1100, "C-O-C asymmetric stretch (ether)", "strong", "Ether antisymmetric stretch"),
        (1100, 1120, 1070, "C-O-C symmetric stretch (cyclic ether)", "strong", "THF-type cyclic ether"),
        (1150, 1170, 1130, "C-O stretch (ester C-O-C)", "strong", "Ester C-O single bond stretch"),
        (1040, 1060, 1020, "C-O stretch (primary alcohol, branched)", "strong", "Branched primary alcohol"),
        (1030, 1050, 1010, "C-O stretch (cyclohexanol type)", "strong", "Cyclic/alicyclic alcohol"),
        (980, 1000, 960, "=C-H bend (out-of-plane, trans alkene)", "strong", "Trans RCH=CHR': strong band ~965"),
        (910, 930, 890, "=C-H bend (out-of-plane, terminal =CH₂)", "strong", "Terminal vinyl: two bands ~990 & ~910"),
        (850, 870, 830, "=C-H bend (oop, trisubstituted alkene)", "strong", "Trisubstituted R₂C=CHR"),
        (800, 820, 780, "=C-H bend (oop, cis alkene)", "strong", "Cis RCH=CHR: ~700"),
        (720, 730, 710, "CH₂ rocking (long-chain polymethylene)", "medium", "(CH₂)n where n≥4: solid indicator of long chain"),
        (770, 780, 740, "C-H oop bend (monosubstituted benzene)", "strong", "Two bands: ~690 & ~750 for monosubstituted benzene"),
        (695, 710, 685, "C-H oop bend (monosubstituted benzene, second band)", "strong", "Companion to ~750 band"),
        (780, 800, 760, "C-H oop bend (ortho-disubstituted benzene)", "strong", "One strong band for ortho substitution"),
        (860, 880, 840, "C-H oop bend (para-disubstituted benzene)", "strong", "One band for para substitution"),
        (780, 800, 760, "C-H oop bend (meta-disubstituted benzene)", "strong", "Two bands for meta substitution (main ~780, minor ~880)"),
        (830, 850, 800, "C-H oop bend (1,2,3-trisubstituted benzene)", "strong", "Pattern changes with substitution"),
        (720, 740, 700, "C-H oop bend (1,2,4-trisubstituted benzene)", "strong", "Complex patterns for tri/pentasubstitution"),
        (690, 710, 670, "C-H oop bend (penta/hexasubstituted benzene)", "weak", "Highly substituted aromatics have weak/absent oop bands"),
        (650, 670, 630, "C-Cl stretch (chloro compound)", "strong", "C-Cl bond vibration"),
        (600, 620, 580, "C-Br stretch (bromo compound)", "strong", "C-Br bond vibration"),
        (550, 570, 530, "C-I stretch (iodo compound)", "medium", "C-I bond vibration (low frequency)"),
        (500, 520, 480, "C-S stretch (thiol/thioether)", "weak", "C-S bond; thiols also S-H at ~2550"),
        (2550, 2580, 2520, "S-H stretch (thiol)", "weak", "Thiol S-H; much weaker than O-H/N-H"),
    ]

    def __init__(
        self,
        init: bool = True,
        interface: str = "code"
    ):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        """Pre-sort database for efficient lookup."""
        self._sorted_db = sorted(self._IR_DATABASE, key=lambda x: x[0])

    def _find_matches(self, wavenumber: float, tolerance: float = 15.0) -> list:
        """Find all IR assignments matching within tolerance."""
        matches = []
        for center, lo, hi, assignment, intensity, notes in self._IR_DATABASE:
            if min(lo, hi) - tolerance <= wavenumber <= max(lo, hi) + tolerance:
                deviation = abs(wavenumber - center)
                matches.append({
                    "wavenumber_input": round(wavenumber, 1),
                    "nearest_reference": center,
                    "deviation_cm": round(deviation, 1),
                    "assignment": assignment,
                    "range": f"{lo}-{hi}",
                    "intensity": intensity,
                    "notes": notes,
                })

        # Sort by deviation (closest match first)
        matches.sort(key=lambda x: x["deviation_cm"])
        return matches

    def _run_base(self, wavenumber: Union[float, List[float]], tolerance: float = 15.0) -> dict:
        """
        Interpret IR peak(s).

        Args:
            wavenumber: Single wavenumber (float) or list of wavenumbers
            tolerance: Matching tolerance in cm⁻¹

        Returns:
            Dict with interpretation results
        """
        if isinstance(wavenumber, (int, float)):
            wavenumbers = [float(wavenumber)]
        elif isinstance(wavenumber, list):
            wavenumbers = [float(w) for w in wavenumber]
        else:
            raise ChemMCPError(f"wavenumber must be float or list, got {type(wavenumber)}")

        if not wavenumbers:
            raise ChemMCPError("At least one wavenumber value is required.")

        results = []
        for wn in wavenumbers:
            matches = self._find_matches(wn, tolerance)
            if matches:
                best = matches[0]  # closest match
                best["all_matches"] = matches[:5]  # include up to 5 alternatives
                results.append(best)
            else:
                results.append({
                    "wavenumber_input": round(wn, 1),
                    "assignment": "No common functional group absorption found in this region.",
                    "note": "This region may correspond to fingerprint region overtones or uncommon vibrations.",
                    "all_matches": [],
                })

        return {
            "interpretations": results,
            "total_peaks": len(results),
            "tolerance_used": tolerance,
        }

    def _run_text(self, input_params: str) -> dict:
        """Parse text input."""
        parts = input_params.strip().split()
        if not parts:
            raise ChemMCPError("Input required. Format: 'wavenumber1 wavenumber2 ...'")

        try:
            wavenumbers = [float(p) for p in parts]
        except ValueError as e:
            raise ChemMCPError(f"Invalid wavenumber value: {e}")

        return self._run_base(wavenumbers)
