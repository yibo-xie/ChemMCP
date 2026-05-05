import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


# SPE sorbent database
SORBENT_DATABASE = {
    "C18": {
        "full_name": "Octadecyl (C18)",
        "mechanism": "Reversed-phase, hydrophobic interaction",
        "retains": ["non_polar", "moderately_non_polar"],
        "typical_analytes": ["PAHs", "pesticides", "PCBs", "drugs", "steroids", "fatty_acids", "aromatics"],
        "weak_retention": ["very_polar", "ionic"],
        "condition_solvent": "Methanol or acetonitrile",
        "equilibration_solvent": "Water or buffer (pH 2-8)",
        "wash_solvent": "5-20% methanol in water",
        "elution_solvent": "Methanol, acetonitrile, or mixture with ethyl acetate",
        "ph_range": "2-8",
        "capacity_mg_g": "5-15%",
        "best_for": "General non-polar to moderately polar compounds from aqueous matrices.",
        "notes": "Most widely used SPE sorbent. Start here if unsure.",
    },
    "C8": {
        "full_name": "Octyl (C8)",
        "mechanism": "Reversed-phase, hydrophobic interaction (weaker than C18)",
        "retains": ["non_polar", "moderately_non_polar"],
        "typical_analytes": ["pesticides", "drugs", "moderately_hydrophobic_compounds"],
        "weak_retention": ["very_polar", "ionic"],
        "condition_solvent": "Methanol or acetonitrile",
        "equilibration_solvent": "Water or buffer",
        "wash_solvent": "5-15% methanol in water",
        "elution_solvent": "Methanol or acetonitrile",
        "ph_range": "2-8",
        "capacity_mg_g": "3-10%",
        "best_for": "Compounds that are too strongly retained on C18; easier elution.",
        "notes": "Shorter chain than C18 — less retention, easier elution of strong hydrophobics.",
    },
    "PS-DVB": {
        "full_name": "Polystyrene-Divinylbenzene (PS-DVB)",
        "mechanism": "Reversed-phase + π-π interactions",
        "retains": ["non_polar", "aromatic", "polar_aromatic"],
        "typical_analytes": ["aromatic_pesticides", "phenols", "explosives", "dyes", "aromatic_amines"],
        "weak_retention": ["aliphatic", "very_polar_aliphatic"],
        "condition_solvent": "Methanol or acetonitrile",
        "equilibration_solvent": "Water or buffer (wide pH range)",
        "wash_solvent": "Water or low organic",
        "elution_solvent": "Acetonitrile, THF, methanol/acetone mixtures",
        "ph_range": "0-14",
        "capacity_mg_g": "High",
        "best_for": "Aromatic compounds, wide pH stability needed.",
        "notes": "Excellent pH stability (0-14). Good for basic compounds at high pH.",
    },
    "SAX": {
        "full_name": "Strong Anion Exchange (SAX, quaternary amine)",
        "mechanism": "Ion exchange — anionic attraction",
        "retains": ["acids", "anions", "acidic_compounds"],
        "typical_analytes": ["carboxylic_acids", "sulfonic_acids", "organic_anions", "nucleotides", "phosphates"],
        "weak_retention": ["neutrals", "cations", "bases"],
        "condition_solvent": "Methanol then water/pH buffer",
        "equilibration_solvent": "Water or buffer (pH > pKa of analyte)",
        "wash_solvent": "Buffer or water/methanol mixture",
        "elution_solvent": "Acidic buffer (pH < 2) or high-salt solution",
        "ph_range": "0-14",
        "capacity_mg_g": "0.7-1.0 meq/g",
        "best_for": "Strong acids and anions requiring ion-exchange retention.",
        "notes": "Always ionized — retains acidic compounds across all pH ranges.",
    },
    "SCX": {
        "full_name": "Strong Cation Exchange (SCX, sulfonic acid)",
        "mechanism": "Ion exchange — cationic attraction",
        "retains": ["bases", "cations", "basic_compounds"],
        "typical_analytes": ["amines", "amino_acids", "basic_drugs", "catecholamines", "metals_as_cations"],
        "weak_retention": ["neutrals", "anions", "acids"],
        "condition_solvent": "Methanol then water/pH buffer",
        "equilibration_solvent": "Water or buffer (pH < pKa of analyte)",
        "wash_solvent": "Buffer or water/methanol mixture",
        "elution_solvent": "Basic buffer (pH > 10) or high-salt solution",
        "ph_range": "0-14",
        "capacity_mg_g": "0.7-1.0 meq/g",
        "best_for": "Strong bases and cations requiring ion-exchange retention.",
        "notes": "Always ionized — retains basic compounds across all pH ranges.",
    },
    "WAX": {
        "full_name": "Weak Anion Exchange (WAX, tertiary amine)",
        "mechanism": "Ion exchange (pH-dependent) — anionic attraction when deprotonated",
        "retains": ["weak_acids", "anions_at_high_pH"],
        "typical_analytes": ["weak_carboxylic_acids", "phenols", "weakly_acidic_drugs"],
        "weak_retention": ["strong_acids_at_low_pH", "neutrals"],
        "condition_solvent": "Methanol then buffer (pH ~8)",
        "equilibration_solvent": "Buffer pH 6-8",
        "wash_solvent": "Buffer or weak organic wash",
        "elution_solvent": "Acidic buffer (pH 2-4) — protonates the sorbent, releasing analytes",
        "ph_range": "0-12",
        "capacity_mg_g": "0.5-0.8 meq/g",
        "best_for": "Weak acids where milder elution conditions are desired.",
        "notes": "Can be eluted under mild acidic conditions by protonating the sorbent.",
    },
    "WCX": {
        "full_name": "Weak Cation Exchange (WCX, carboxylic acid)",
        "mechanism": "Ion exchange (pH-dependent) — cationic attraction when deprotonated",
        "retains": ["weak_bases", "cations_at_low_pH"],
        "typical_analytes": ["weak_bases", "amines_with_pKa_5-8", "some_drugs"],
        "weak_retention": ["strong_bases_at_high_pH", "neutrals"],
        "condition_solvent": "Methanol then buffer (pH ~6)",
        "equilibration_solvent": "Buffer pH 4-6",
        "wash_solvent": "Buffer or weak organic wash",
        "elution_solvent": "Basic buffer (pH 9-11) — deprotonates sorbent, releases analytes",
        "ph_range": "2-12",
        "capacity_mg_g": "0.5-0.8 meq/g",
        "best_for": "Weak bases where milder elution is acceptable.",
        "notes": "Can be eluted under mild basic conditions.",
    },
    "HLB": {
        "full_name": "Hydrophilic-Lipophilic Balance (N-vinylpyrrolidone-divinylbenzene copolymer)",
        "mechanism": "Balanced reversed-phase with hydrophilic character",
        "retains": ["wide_range", "polar_to_nonpolar"],
        "typical_analytes": ["drugs", "pesticides", "antibiotics", "mycotoxins", "polar_pesticides", "multiclass"],
        "weak_retention": [],
        "condition_solvent": "Methanol or acetonitrile",
        "equilibration_solvent": "Water or buffer (pH 0-14)",
        "wash_solvent": "5% methanol/water (removes interferences without losing analytes)",
        "elution_solvent": "Methanol, acetonitrile, or mixture with MTBE/ethyl acetate",
        "ph_range": "0-14",
        "capacity_mg_g": "High (~20-30%)",
        "best_for": "Multiclass multi-residue analysis; broad-spectrum extraction from various matrices.",
        "notes": "Universal sorbent for multi-residue methods. Excellent water wettability — won't dry out if run dry.",
    },
    "Florisil": {
        "full_name": "Florisil (magnesium silicate)",
        "mechanism": "Normal-phase adsorption (polar interactions)",
        "retains": ["polar_compounds_from_nonpolar"],
        "typical_analytes": ["pesticides_from_petroleum_ether", "lipids", "pigments", "organochlorines"],
        "weak_retention": ["nonpolar_in_polar_solvent"],
        "condition_solvent": "Dry solvent (hexane, DCM)",
        "equilibration_solvent": "Non-polar solvent (hexane)",
        "wash_solvent": "Hexane or hexane with small % polar modifier",
        "elution_solvent": "Acetone, ethyl acetate, or ethyl acetate/hexane mixtures",
        "ph_range": "N/A (normal phase)",
        "capacity_mg_g": "Moderate",
        "best_for": "Cleanup of lipid-rich samples; pesticide residue analysis from fatty matrices.",
        "notes": "Must keep sorbent dry during use. Commonly used for pesticide cleanup after extraction.",
    },
    "Silica": {
        "full_name": "Silica gel (SiO2)",
        "mechanism": "Normal-phase adsorption (hydrogen bonding, dipole-dipole)",
        "retains": ["polar_compounds_from_nonpolar_solvent"],
        "typical_analytes": ["surfactants", "glycerides", "polar_pesticides", "pigments"],
        "weak_retention": ["nonpolar_compounds"],
        "condition_solvent": "Non-polar solvent (hexane)",
        "equilibration_solvent": "Hexane or non-polar solvent",
        "wash_solvent": "Hexane or hexane/dichloromethane",
        "elution_solvent": "Methanol, acetone, or polar solvent in hexane",
        "ph_range": "N/A (normal phase)",
        "capacity_mg_g": "Moderate (5-10%)",
        "best_for": "Normal-phase cleanup; fractionation based on polarity.",
        "notes": "Keep dry during normal-phase operation. Deactivate with water for modified polarity.",
    },
    "GraphitizedCarbon_Black": {
        "full_name": "Graphitized Carbon Black (GCB)",
        "mechanism": "Planar adsorption + specific interactions",
        "retains": ["planar_molecules", "pigments", "chlorophyll", "steroids"],
        "typical_analytes": ["planar_pesticides", "pigments", "chlorophyll_cleanup", "PAHs"],
        "weak_retention": ["non_planar", "aliphatic"],
        "condition_solvent": "Methanol, then water/buffer",
        "equilibration_solvent": "Water or buffer",
        "wash_solvent": "Water/methanol mixture",
        "elution_solvent": "Toluene/MeOH mixture, dichloromethane/methanol, or ACN/toluene",
        "ph_range": "1-13",
        "capacity_mg_g": "Moderate-High",
        "best_for": "Removal of pigments (chlorophyll) and planar interferences from plant extracts.",
        "notes": "Excellent for removing chlorophyll from plant matrix. Can also retain planar pesticides too strongly — optimize amount carefully.",
    },
}


@ChemMCPManager.register_tool
class SPEMethodDesigner(BaseTool):
    """
    固相萃取方法设计：选择填料和洗脱条件。
    根据目标分析物性质、样品基质和检测要求设计完整SPE方案。
    """
    __version__ = "0.1.0"
    name = "SPEMethodDesigner"
    func_name = "design_spe_method"
    description = "Design a complete solid-phase extraction (SPE) method including sorbent selection, conditioning, loading, washing, and elution protocols."
    implementation_description = "Uses a comprehensive SPE sorbent database with rule-based matching of analyte properties, matrix type, and detection requirements."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["SPE", "Solid Phase Extraction", "Sample Preparation", "Analytical Chemistry"]
    required_envs = []

    code_input_sig = [
        ("analytes", "str", "N/A", "Target analytes, comma-separated (e.g., 'pesticides, PAHs, drugs')."),
        ("analyte_properties", "str", "nonpolar", "Analyte properties: 'nonpolar', 'polar', 'acidic', 'basic', 'amphoteric', 'mixed'."),
        ("sample_matrix", "str", "water", "Sample matrix: 'water', 'urine', 'blood', 'plasma', 'food', 'soil_extract', 'environmental_water'."),
        ("sample_volume_ml", "float", "100.0", "Sample volume to load (mL)."),
        ("detection_method", "str", "LC-MS", "Detection method: 'LC-MS', 'GC-MS', 'UV', 'fluorescence', 'general'."),
        ("ph_adjustment", "float", "None", "Desired sample pH before loading (None for auto-recommend)."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'analytes analyte_properties sample_matrix [sample_volume_ml] [detection_method] [ph_adjustment]'"),
    ]

    output_sig = [
        ("recommended_sorbent", "str", "Recommended SPE sorbent/cartridge type."),
        ("method_protocol", "dict", "Complete SPE protocol with step-by-step instructions."),
        ("rationale", "str", "Explanation of sorbent selection."),
        ("alternatives", "list", "Alternative sorbents if primary not available."),
        ("tips", "list", "Practical tips and troubleshooting."),
    ]

    examples = [
        {
            "code_input": {
                "analytes": "pesticides, drugs",
                "analyte_properties": "mixed",
                "sample_matrix": "water",
                "sample_volume_ml": 500.0,
                "detection_method": "LC-MS",
            },
            "text_input": {
                "input_params": "pesticides,drugs mixed water 500.0 LC-MS",
            },
            "output": {
                "recommended_sorbent": "Hydrophilic-Lipophilic Balance (HLB)",
                "rationale": "HLB provides balanced retention for both polar and non-polar analytes with high capacity.",
            },
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(
        self,
        analytes: str,
        analyte_properties: str = "nonpolar",
        sample_matrix: str = "water",
        sample_volume_ml: float = 100.0,
        detection_method: str = "LC-MS",
        ph_adjustment: float = None,
    ) -> dict:
        """Core logic: design complete SPE method."""
        props = analyte_properties.lower().strip()
        matrix = sample_matrix.lower().strip()
        analytes_lower = analytes.lower()

        # Sorbent selection logic
        sorbent_key = self._select_sorbent(props, analytes_lower, matrix)
        sorb = SORBENT_DATABASE[sorbent_key]

        # Auto pH recommendation
        ph_rec = ph_adjustment
        ph_rationale = ""
        if ph_rec is None:
            ph_rec, ph_rationale = self._recommend_ph(props, sorbent_key)

        # Build protocol
        protocol = self._build_protocol(sorb, sorbent_key, matrix, sample_volume_ml, ph_rec, detection_method)

        # Alternatives
        alternatives = []
        for key, s in SORBENT_DATABASE.items():
            if key != sorbent_key:
                alt_score = self._score_sorbent_for(key, props, analytes_lower, matrix)
                if alt_score >= 0.4:
                    alternatives.append({
                        "sorbent": s["full_name"],
                        "key": key,
                        "score": round(alt_score, 2),
                        "why_choose_this": s["best_for"],
                    })
        alternatives.sort(key=lambda x: x["score"], reverse=True)

        # Tips
        tips = self._generate_tips(sorbent_key, matrix, sample_volume_ml)

        logger.info(f"SPE method designed: {sorbent_key} for {analytes} from {matrix}")
        return {
            "recommended_sorbent": sorb["full_name"],
            "method_protocol": protocol,
            "rationale": f"{sorb['full_name']} selected: {sorb['best_for']}",
            "ph_recommendation": {"target_ph": ph_rec, "reasoning": ph_rationale} if ph_rationale else None,
            "alternatives": alternatives[:4],
            "tips": tips,
        }

    def _select_sorbent(self, props, analytes, matrix):
        """Select best sorbent based on analyte properties."""
        if props == "mixed" or props == "amphoteric":
            return "HLB"
        elif props == "acidic":
            return "WAX"
        elif props == "basic":
            return "WCX"
        elif props == "polar":
            return "PS-DVB"
        elif props == "nonpolar":
            # Check for specific analyte hints
            if any(a in analytes for a in ["pesticide", "drug", "pa"]):
                return "C18"
            elif any(a in analytes for a in ["aromatic", "phenol", "dye"]):
                return "PS-DVB"
            else:
                return "C18"
        return "C18"

    def _recommend_ph(self, props, sorbent_key):
        """Recommend sample pH."""
        if props == "acidic":
            return 6.0, "pH 6-7 ensures acidic analytes are deprotonated (ionized) for anion exchange retention on WAX/SAX."
        elif props == "basic":
            return 8.0, "pH 8-9 ensures basic analytes are protonated (ionized) for cation exchange retention on WCX/SCX."
        elif props == "polar":
            return 7.0, "Neutral pH minimizes ionization effects for RP retention of polar compounds."
        else:
            return 7.0, "Neutral pH suitable for general reversed-phase SPE."

    def _build_protocol(self, sorb, key, matrix, vol_ml, ph, detection):
        """Build step-by-step SPE protocol."""
        # Condition volume depends on cartridge size
        cart_size = self._recommend_cartridge_size(vol_ml)

        protocol = {
            "cartridge_recommendation": cart_size,
            "steps": [
                {
                    "step": 1,
                    "name": "Conditioning (活化)",
                    "action": f"Pass {cart_size.split('(')[0].strip()} mL of {sorb['condition_solvent']} through cartridge. Do NOT let sorbent dry completely.",
                    "solvent": sorb["condition_solvent"],
                    "volume_ml": float("".join(c for c in cart_size.split("(")[0].strip() if c.isdigit() or c == ".")) or 5.0,
                },
                {
                    "step": 2,
                    "name": "Equilibration (平衡)",
                    "action": f"Pass same volume of {sorb['equilibration_solvent']} through cartridge.",
                    "solvent": sorb["equilibration_solvent"],
                    "volume_ml": float("".join(c for c in cart_size.split("(")[0].strip() if c.isdigit() or c == ".")) or 5.0,
                },
                {
                    "step": 3,
                    "name": "Sample Loading (上样)",
                    "action": f"Load {vol_ml} mL of sample at flow rate 1-5 mL/min. Adjust pH to {ph} before loading if necessary.",
                    "flow_rate_ml_min": "1-5",
                    "note": "For large volumes (>250mL), consider using vacuum manifold or positive pressure processor.",
                },
                {
                    "step": 4,
                    "name": "Washing (淋洗)",
                    "action": f"Pass {sorb['wash_solvent']}. Purpose: remove matrix interferences while retaining target analytes.",
                    "solvent": sorb["wash_solvent"],
                    "volume_ml": 5.0,
                },
                {
                    "step": 5,
                    "name": "Drying (可选干燥)",
                    "action": "Dry cartridge under vacuum for 5 min (or pass air/N2) to remove residual water before elution with organic solvent.",
                    "optional": True,
                },
                {
                    "step": 6,
                    "name": "Elution (洗脱)",
                    "action": f"Elute with {sorb['elution_solvent']}. Collect eluent in clean tube. Use 2 × 2 mL for complete recovery.",
                    "solvent": sorb["elution_solvent"],
                    "volume_ml": 4.0,
                    "note": "Elute slowly (<1 mL/min) for maximum recovery.",
                },
            ],
            "post_elution": {
                "evaporation": "Evaporate eluent to near-dryness under gentle nitrogen stream at 35-40°C." if detection in ("GC-MS", "LC-MS") else "Dilute to final volume as needed.",
                "reconstitution": "Reconstitute in 200 μL initial mobile phase or compatible solvent." if detection in ("LC-MS", "GC-MS") else "Ready for analysis or dilute as needed.",
            },
        }
        return protocol

    def _recommend_cartridge_size(self, vol_ml):
        if vol_ml <= 50:
            return "1 mL (30 mg)"
        elif vol_ml <= 250:
            return "3 mL (60 mg)"
        elif vol_ml <= 1000:
            return "6 mL (200 mg)"
        else:
            return "12-15 mL (500 mg - 1 g)"

    def _score_sorbent_for(self, key, props, analytes, matrix):
        score = 0.0
        s = SORBENT_DATABASE[key]
        if props in s.get("retains", []) or props.replace("polar", "") in str(s.get("retains", [])):
            score += 2.0
        if any(a in str(s.get("typical_analytes", [])) for a in analytes.split(",")):
            score += 2.0
        if props == "mixed" and key == "HLB":
            score += 3.0
        if props == "acidic" and "Anion" in s["full_name"]:
            score += 2.0
        if props == "basic" and "Cation" in s["full_name"]:
            score += 2.0
        if props == "nonpolar" and key in ("C18", "C8"):
            score += 1.5
        if props == "polar" and key in ("PS-DVB", "HLB"):
            score += 1.5
        return max(0.0, min(score / 4.0, 1.0))

    def _generate_tips(self, key, matrix, vol_ml):
        tips = [
            "Never let sorbent dry completely between conditioning and sample loading (except GCB/silica normal-phase).",
            "Flow rates: loading ≤5 mL/min, washing ≤5 mL/min, elution ≤1 mL/min.",
            "Pre-filter or centrifuge samples containing particulate matter before SPE.",
        ]
        if vol_ml > 500:
            tips.append(f"For large sample volumes ({vol_ml} mL), consider using disk format SPE instead of cartridge for faster throughput.")
        if matrix in ("blood", "plasma"):
            tips.append("For biological fluids: protein precipitation (e.g., 3x ACN) before SPE can extend cartridge life.")
        if key in ("C18", "C8"):
            tips.append("C18/C8 cartridges should be stored wetted with methanol to prevent degradation of bonded phase.")
        return tips

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.split()
            analytes = parts[0]
            props = parts[1] if len(parts) > 1 else "nonpolar"
            matrix = parts[2] if len(parts) > 2 else "water"
            vol = float(parts[3]) if len(parts) > 3 else 100.0
            det = parts[4] if len(parts) > 4 else "LC-MS"
            ph = float(parts[5]) if len(parts) > 5 else None
            return self._run_base(analytes, props, matrix, vol, det, ph)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse text input: {str(e)}. Format: 'analytes props [matrix] [vol] [det] [pH]'")
