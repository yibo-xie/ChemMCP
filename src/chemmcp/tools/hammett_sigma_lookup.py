import logging
from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HammettSigmaLookup(BaseTool):
    """
    Hammett σ 常数查询工具 - 查询取代基的 Hammett σm、σp、σp+ 值，并解释电子效应。
    覆盖常见取代基的完整 Hammett 参数数据库。
    """
    __version__ = "0.1.0"
    name             = "HammettSigmaLookup"
    func_name        = "hammett_sigma_lookup"
    description      = "Look up Hammett σ constants (σm, σp, σp+) for substituents and interpret their electronic effects on aromatic reactivity and acidity."
    implementation_description = "Comprehensive database of Hammett substituent constants with electronic effect interpretation, including σm (meta), σp (para), σp+ (para for cationic systems), and Swain-Lupton F/R parameters."
    oss_dependencies = []
    services_and_software = []
    categories       = ["General"]
    tags             = ["Hammett", "Substituent Effects", "Physical Organic Chemistry", "Linear Free Energy Relationships"]
    required_envs    = []

    code_input_sig   = [
        ("substituent", "str", "N/A", "Substituent name or notation: e.g., 'p-NO2', 'm-Cl', 'p-OCH3', 'NMe2', 'CF3', 'H' for reference."),
        ("include_swain_lupton", "bool", "False", "Whether to include Swain-Lupton field (F) and resonance (R) parameters."),
    ]

    text_input_sig   = [
        ("input_params", "str", "N/A", "Substituent query. Example: 'p-nitro' or 'm-chloro' or 'para-methoxy'."),
    ]

    output_sig       = [
        ("result", "str", "Hammett constants table, electronic effect interpretation, and application guidance."),
    ]

    examples         = [
        {
            "code_input": {"substituent": "p-NO2", "include_swain_lupton": False},
            "text_input": {"input_params": "p-NO2"},
            "output": {"result": "σm=0.71, σp=0.78... Strong electron-withdrawing..."},
        },
        {
            "code_input": {"substituent": "p-OCH3", "include_swain_lupton": False},
            "text_input": {"input_params": "p-OCH3"},
            "output": {"result": "σm=0.10, σp=-0.27... Electron-donating via resonance..."},
        },
    ]

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        self._build_database()

    def _build_database(self):
        """Build comprehensive Hammett constant database."""
        # Format: {canonical_name: {
        #   display_names: [aliases],
        #   sigma_m: float, sigma_p: float, sigma_p_plus: float (or None),
        #   sigma_p_minus: float (or None),
        #   swain_F: float, swain_R: float,
        #   classification: str,
        #   interpretation: str,
        # }}
        self.database = {
            "h": {
                "display_names": ["H", "h", "hydrogen"],
                "sigma_m": 0.00, "sigma_p": 0.00, "sigma_p_plus": 0.00, "sigma_p_minus": 0.00,
                "swain_F": 0.00, "swain_R": 0.00,
                "classification": "Reference",
                "interpretation": "Reference substituent. No electronic effect by definition.",
            },
            # === Strong Electron Donors (+R dominant) ===
            "nme2": {
                "display_names": ["NMe2", "nme2", "N(CH3)2", "dimethylamino", "p-NMe2", "p-dimethylamino", "para-NMe2", "para-dimethylamino"],
                "sigma_m": -0.15, "sigma_p": -0.83, "sigma_p_plus": -1.7, "sigma_p_minus": -0.15,
                "swain_F": 0.10, "swain_R": -0.83,
                "classification": "Strong +R donor",
                "interpretation": "Very strong electron donor via resonance (+R). The lone pair on N can delocalize into the ring at para position. At meta position, only weak -I from N electronegativity is felt (σm slightly positive). In σp+ (cationic transition states), donation is dramatically enhanced (σp+ = -1.7). Dominates in electrophilic aromatic substitution (strongly activating/ortho,para-directing).",
            },
            "nh2": {
                "display_names": ["NH2", "nh2", "amino", "p-NH2", "para-amino"],
                "sigma_m": -0.16, "sigma_p": -0.66, "sigma_p_plus": -1.3, "sigma_p_minus": -0.12,
                "swain_F": 0.02, "swain_R": -0.68,
                "classification": "Strong +R donor",
                "interpretation": "Strong electron donor via resonance. Similar to NMe2 but slightly less donating due to N-H bond polarity. Strongly activates benzene toward EAS. Protonation in acid converts it to -EWG (NH3+).",
            },
            "oh": {
                "display_names": ["OH", "oh", "hydroxy", "hydroxyl", "p-OH", "para-hydroxy", "phenol-type"],
                "sigma_m": 0.12, "sigma_p": -0.37, "sigma_p_plus": -0.92, "sigma_p_minus": -0.25,
                "swain_F": 0.29, "swain_R": -0.64,
                "classification": "+R donor / -I acceptor",
                "interpretation": "Net electron donor at para position due to strong +R (lone pair donation into ring) outweighing -I (O electronegativity). At meta, only -I operates → σm > 0. Phenols are much more acidic than alcohols because phenoxide is stabilized by this +R delocalization.",
            },
            "och3": {
                "display_names": ["OMe", "ome", "OCH3", "och3", "methoxy", "p-OMe", "p-OCH3", "para-methoxy", "para-methoxy"],
                "sigma_m": 0.10, "sigma_p": -0.27, "sigma_p_plus": -0.78, "sigma_p_minus": -0.20,
                "swain_F": 0.27, "swain_R": -0.51,
                "classification": "+R donor / -I acceptor",
                "interpretation": "Electron-donating at para via resonance (+R), similar to OH but weaker (the methyl reduces O's electron density availability). Common ortho/para director in EAS. Anisole is activated toward nitration, halogenation etc.",
            },
            "oc2h5": {
                "display_names": ["OEt", "oet", "OC2H5", "oc2h5", "ethoxy", "p-OEt", "para-ethoxy"],
                "sigma_m": 0.10, "sigma_p": -0.24, "sigma_p_plus": -0.70,
                "swain_F": 0.27, "swain_R": -0.48,
                "classification": "+R donor",
                "interpretation": "Similar to OCH3 but slightly less donating due to larger alkyl group.",
            },
            
            # === Weak Electron Donors (alkyl / hyperconjugation) ===
            "ch3": {
                "display_names": ["Me", "me", "CH3", "ch3", "methyl", "p-Me", "p-CH3", "para-methyl", "p-cresol-type"],
                "sigma_m": -0.07, "sigma_p": -0.17, "sigma_p_plus": -0.31, "sigma_p_minus": -0.10,
                "swain_F": -0.04, "swain_R": -0.13,
                "classification": "Weak +I donor (hyperconjugation)",
                "interpretation": "Weak electron donor via hyperconjugation (σC-H bond donation) and inductive (+I) effect. Slightly activates ring toward EAS. Ortho/para-directing. Effect is small but consistent across many reactions.",
            },
            "c2h5": {
                "display_names": ["Et", "et", "C2H5", "c2h5", "ethyl", "p-Et", "para-ethyl"],
                "sigma_m": -0.07, "sigma_p": -0.15, "sigma_p_plus": -0.30,
                "swain_F": -0.05, "swain_R": -0.10,
                "classification": "Weak +I donor",
                "interpretation": "Similar to CH3; very weak hyperconjugative donor.",
            },
            "ch(ch3)2": {
                "display_names": ["i-Pr", "ipr", "CH(CH3)2", "isopropyl", "p-i-Pr", "para-isopropyl"],
                "sigma_m": -0.05, "sigma_p": -0.15, "sigma_p_plus": -0.28,
                "swain_F": -0.03, "swain_R": -0.12,
                "classification": "Weak +I donor",
                "interpretation": "Weakly donating; steric bulk can influence reactions beyond electronic effects.",
            },
            "c(ch3)3": {
                "display_names": ["t-Bu", "tbu", "C(CH3)3", "tert-butyl", "p-t-Bu", "para-tert-butyl"],
                "sigma_m": -0.10, "sigma_p": -0.20, "sigma_p_plus": -0.26,
                "swain_F": -0.08, "swain_R": -0.12,
                "classification": "Weak +I donor (with steric effects)",
                "interpretation": "Slightly more donating than CH3 due to three methyl groups' cumulative hyperconjugation. Significant steric bulk can block approach of reagents (steric inhibition).",
            },
            "ph": {
                "display_names": ["Ph", "ph", "phenyl", "C6H5", "p-Ph", "para-phenyl"],
                "sigma_m": 0.06, "sigma_p": -0.01, "sigma_p_plus": -0.18,
                "swain_F": 0.08, "swain_R": -0.08,
                "classification": "Very weak donor (π-system)",
                "interpretation": "Nearly neutral electronically. The phenyl group can act as a very weak π-donor or π-acceptor depending on context. Often considered 'spectator' electronically but adds conjugation length.",
            },

            # === Halogens (-I > +R, unique case) ===
            "f": {
                "display_names": ["F", "f", "fluoro", "fluorine", "p-F", "p-F", "para-fluoro", "para-fluorine", "m-F", "meta-fluoro", "meta-fluorine"],
                "sigma_m": 0.34, "sigma_p": 0.06, "sigma_p_plus": -0.07, "sigma_p_minus": 0.34,
                "swain_F": 0.45, "swain_R": -0.39,
                "classification": "Strong -I, strong +R (unique)",
                "interpretation": "Most electronegative element → strongest -I effect. But also has lone pairs that can donate into ring via +R (especially para). Net result: deactivating (σm > 0, σp > 0) but ortho/para-directing in EAS (due to +R control of orientation). This apparent contradiction (deactivating yet o/p-directing) is classic organic chemistry!",
            },
            "cl": {
                "display_names": ["Cl", "cl", "chloro", "chlorine", "p-Cl", "para-chloro", "m-Cl", "meta-chloro"],
                "sigma_m": 0.37, "sigma_p": 0.23, "sigma_p_plus": 0.11, "sigma_p_minus": 0.23,
                "swain_F": 0.42, "swain_R": -0.19,
                "classification": "Moderate -I, weak +R",
                "interpretation": "Deactivating (net EWG) but ortho/para-directing. -I dominates over +R. Larger size means poorer orbital overlap for +R than F. Chlorobenzene undergoes EAS slower than benzene but gives ortho+para products.",
            },
            "br": {
                "display_names": ["Br", "br", "bromo", "bromine", "p-Br", "para-bromo", "m-Br", "meta-bromo"],
                "sigma_m": 0.39, "sigma_p": 0.23, "sigma_p_plus": 0.15, "sigma_p_minus": 0.23,
                "swain_F": 0.44, "swain_R": -0.21,
                "classification": "Moderate -I, weak +R",
                "interpretation": "Similar to Cl but slightly stronger -I (better polarizability). Same deactivating/o,p-directing pattern.",
            },
            "i": {
                "display_names": ["I", "i", "iodo", "iodine", "p-I", "para-iodo", "m-I", "meta-iodo"],
                "sigma_m": 0.35, "sigma_p": 0.18, "sigma_p_plus": 0.14, "sigma_p_minus": 0.24,
                "swain_F": 0.40, "swain_R": -0.22,
                "classification": "Moderate -I, weak +R",
                "interpretation": "Least electronegative halogen but most polarizable. Weakest deactivating halogen in EAS. Still ortho/para-directing.",
            },

            # === Strong Electron Withdrawers (-I, -R) ===
            "cf3": {
                "display_names": ["CF3", "cf3", "trifluoromethyl", "p-CF3", "para-CF3", "para-trifluoromethyl", "m-CF3", "meta-CF3"],
                "sigma_m": 0.43, "sigma_p": 0.54, "sigma_p_plus": 0.61,
                "swain_F": 0.38, "swain_R": 0.19,
                "classification": "Strong -I (dominant)",
                "interpretation": "Very strong electron-withdrawing group primarily through -I induction (three highly electronegative F atoms pull electron density through σ-bonds). Minimal +R ability (no π-lone pairs). Strongly deactivating, meta-directing in EAS. Dramatically increases acidity when attached to phenols/carboxylic acids.",
            },
            "cn": {
                "display_names": ["CN", "cn", "cyano", "cyano group", "p-CN", "para-CN", "para-cyano", "m-CN", "meta-cyano"],
                "sigma_m": 0.56, "sigma_p": 0.66, "sigma_p_plus": 0.66, "sigma_p_minus": 0.90,
                "swain_F": 0.51, "swain_R": 0.19,
                "classification": "Strong -I, moderate -R",
                "interpretation": "Strong electron-withdrawing via both -I (sp N is electronegative) and -R (C≡N π* accepts electrons). One of the strongest neutral EWGs. Meta-directing, strongly deactivating. Lowers pKa of benzoic acid from 4.20 to 3.55 (p-cyanobenzoic acid).",
            },
            "no2": {
                "display_names": ["NO2", "no2", "nitro", "p-NO2", "para-NO2", "para-nitro", "m-NO2", "meta-nitro", "nitro group"],
                "sigma_m": 0.71, "sigma_p": 0.78, "sigma_p_plus": 1.27, "sigma_p_minus": 1.24,
                "swain_F": 0.65, "swain_R": 0.16,
                "classification": "Very strong -I, strong -R",
                "interpretation": "The strongest common neutral EWG. Powerful -I (N+ electronegativity) combined with strong -R (π* orbital accepts electron density). Extremely deactivating, meta-directing. p-Nitrobenzoic acid pKa = 3.41 (vs 4.20 for PhCOOH). Picric acid (2,4,6-trinitrophenol) has pKa ≈ 0.38 — as strong as some mineral acids! σp+ >> σp shows enhanced withdrawal in cationic TS (e.g., SN1 rate acceleration).",
            },
            "so2ch3": {
                "display_names": ["SO2Me", "so2ch3", "methylsulfonyl", "mesyl", "p-SO2CH3", "para-SO2Me"],
                "sigma_m": 0.60, "sigma_p": 0.72, "sigma_p_plus": 1.00,
                "swain_F": 0.59, "swain_R": 0.18,
                "classification": "Strong -I, strong -R",
                "interpretation": "Similar to NO2 in strength. S(VI) is highly oxidized and electron-deficient. Very strongly deactivating and meta-directing.",
            },
            "so2r_generic": {
                "display_names": ["SO2R", "SO2Ar", "sulfonyl", "sulfone", "p-SO2Ph"],
                "sigma_m": 0.68, "sigma_p": 0.80,
                "swain_F": 0.62, "swain_R": 0.22,
                "classification": "Strong -I, strong -R",
                "interpretation": "General sulfonyl group. Strong EWG similar to NO2.",
            },
            "sor": {
                "display_names": ["SOR", "sulfinyl", "p-SOMe"],
                "sigma_m": 0.40, "sigma_p": 0.35,
                "swain_F": 0.40, "swain_R": 0.0,
                "classification": "Moderate -I",
                "interpretation": "Sulfoxide (S=O). Moderate EWG, weaker than sulfone.",
            },
            "sr": {
                "display_names": ["SR", "thioalkyl", "p-SMe", "para-SMe", "methylthio"],
                "sigma_m": 0.23, "sigma_p": 0.03, "sigma_p_plus": -0.21, "sigma_p_minus": 0.36,
                "swain_F": 0.28, "swain_R": -0.24,
                "classification": "Weak -I, moderate +R",
                "interpretation": "Thioether analog of OR. More polarizable than O, so +R is weaker relative to -I. Slightly deactivating overall but can be o/p-directing in some contexts.",
            },

            # === Carbonyl-containing groups ===
            "cho": {
                "display_names": ["CHO", "cho", "formyl", "aldehyde", "p-CHO", "para-formyl", "p-aldehyde"],
                "sigma_m": 0.35, "sigma_p": 0.42, "sigma_p_plus": 0.73, "sigma_p_minus": 0.73,
                "swain_F": 0.31, "swain_R": 0.14,
                "classification": "Moderate -I, moderate -R",
                "interpretation": "Aldehyde carbonyl withdraws electrons via -I (polar C=O) and -R (π* acceptance). Meta-directing, deactivating. Benzaldehyde is ~100× less reactive than benzene toward EAS (nitration requires harsher conditions).",
            },
            "coch3": {
                "display_names": ["COMe", "come", "COCH3", "acetyl", "p-COMe", "p-acetyl", "para-acetyl", "acetyl group"],
                "sigma_m": 0.38, "sigma_p": 0.50, "sigma_p_plus": 0.87, "sigma_p_minus": 0.87,
                "swain_F": 0.32, "swain_R": 0.20,
                "classification": "Moderate -I, moderate -R",
                "interpretation": "Ketone carbonyl. Similar to CHO but slightly stronger due to additional methyl electron donation to carbonyl carbon making it more electron-deficient? Actually, alkyl ketones are generally similar or slightly stronger EWGs than aldehydes in practice. Acetophenone is meta-directing.",
            },
            "cooh": {
                "display_names": ["COOH", "cooh", "carboxy", "carboxylic acid", "p-COOH", "para-carboxy", "p-carboxylic acid"],
                "sigma_m": 0.37, "sigma_p": 0.45, "sigma_p_plus": 0.85, "sigma_p_minus": 0.85,
                "swain_F": 0.33, "swain_R": 0.17,
                "classification": "Moderate -I, moderate -R",
                "interpretation": "Carboxylic acid. Behaves similarly to aldehyde/ketone (carbonyl EWG). Benzoic acid is the reference compound for σp determination! At low pH (protonated COOH), it's even more withdrawing.",
            },
            "coor": {
                "display_names": ["COOR", "coor", "ester", "p-COOMe", "p-COOEt", "para-ester"],
                "sigma_m": 0.37, "sigma_p": 0.45, "sigma_p_plus": 0.85, "sigma_p_minus": 0.85,
                "swain_F": 0.33, "swain_R": 0.17,
                "classification": "Moderate -I, moderate -R",
                "interpretation": "Ester carbonyl. Nearly identical to COOH in electronic effect (the OR vs OH difference is minimal for Hammett parameters). Methyl benzoate is meta-directing.",
            },
            "conhr": {
                "display_names": ["CONH2", "conh2", "carbamoyl", "amide", "p-CONH2", "para-amide"],
                "sigma_m": 0.28, "sigma_p": 0.36, "sigma_p_plus": 0.61,
                "swain_F": 0.24, "swain_R": 0.14,
                "classification": "Weak-moderate -I, weak -R",
                "interpretation": "Amide carbonyl. Weaker EWG than ester/aldehyde because N donates electrons into carbonyl (reducing its electron deficiency). Benzamide is still meta-directing but less deactivating than acetophenone.",
            },
            "cor": {
                "display_names": ["COR_general", "acyl"],
                "sigma_m": 0.38, "sigma_p": 0.50,
                "swain_F": 0.33, "swain_R": 0.19,
                "classification": "Acyl (general)",
                "interpretation": "General acyl group. All acyl groups are meta-directing deactivators.",
            },

            # === Other groups ===
            "nh3+": {
                "display_names": ["NH3+", "nh3+", "anilinium", "protonated amino"],
                "sigma_m": 0.86, "sigma_p": 0.86,
                "swain_F": 0.86, "swain_R": 0.0,
                "classification": "Very strong -I (full positive charge)",
                "interpretation": "Protonated amine — formally positively charged. Extremely electron-withdrawing (stronger than NO2!). Anilinium ion is ~10,000× less reactive than benzene toward EAS. Meta-directing (only position not directly adjacent to positive charge).",
            },
            "nr3+": {
                "display_names": ["NR3+", "trimethylanilinium", "quaternary ammonium"],
                "sigma_m": 0.88, "sigma_p": 0.82,
                "swain_F": 0.88, "swain_R": 0.0,
                "classification": "Very strong -I (positive charge)",
                "interpretation": "Quaternary ammonium — permanent positive charge. Similar to NH3+. Extremely deactivating.",
            },
            "sihme3": {
                "display_names": ["SiMe3", "sime3", "trimethylsilyl", "TMS", "p-SiMe3"],
                "sigma_m": -0.07, "sigma_p": -0.07, "sigma_p_plus": -0.21,
                "swain_F": -0.13, "swain_R": 0.07,
                "classification": "Weak +I donor (β-silicon effect)",
                "interpretation": "Silicon is electropositive relative to carbon (β-silicon effect). Acts as an electron donor (opposite to carbon analogs!). Unique element where Si is less EN than C.",
            },
            "b(oh)2": {
                "display_names": ["B(OH)2", "b(oh)2", "boronic acid", "p-B(OH)2"],
                "sigma_m": 0.10, "sigma_p": 0.45,
                "swain_F": 0.28, "swain_R": 0.19,
                "classification": "Para-stronger EWG",
                "interpretation": "Boron is electron-deficient (empty p-orbital). Stronger effect at para where π-interaction with empty B p-orbital is possible. Used extensively in Suzuki coupling.",
            },
            "ch=ch2": {
                "display_names": ["CH=CH2", "vinyl", "ethenyl", "p-vinyl"],
                "sigma_m": 0.06, "sigma_p": -0.02, "sigma_p_plus": -0.16,
                "swain_F": 0.08, "swain_R": -0.08,
                "classification": "Very weak donor/neutral",
                "interpretation": "Vinyl group. Nearly neutral, similar to Ph. Can extend conjugation.",
            },
            "cch": {
                "display_names": ["C≡CH", "ethynyl", "p-ethynyl"],
                "sigma_m": 0.21, "sigma_p": 0.23,
                "swain_F": 0.22, "swain_R": 0.03,
                "classification": "Weak EWG",
                "interpretation": "Ethynyl group. sp-carbon is relatively electronegative. Slightly withdrawing.",
            },
            "no": {
                "display_names": ["NO", "nitroso", "p-NO"],
                "sigma_m": 0.12, "sigma_p": -0.12, "sigma_p_plus": 0.0, "sigma_p_minus": 0.63,
                "swain_F": 0.24, "swain_R": -0.33,
                "classification": "+R donor (like O)",
                "interpretation": "Nitroso group. Has lone pair on N that can donate (+R), similar to OH/NH2 behavior. Unusual in having large σp- value.",
            },
            "n2+": {
                "display_names": ["N2+", "diazonium"],
                "sigma_m": 1.8, "sigma_p": 1.9, "sigma_p_plus": 2.0,
                "swain_F": 1.9, "swain_R": 0.1,
                "classification": "Extreme -I (positive charge + good leaving group)",
                "interpretation": "Diazonium ion — extremely electron-withdrawing due to formal positive charge on N. Arenediazonium ions are very reactive toward substitution (Sandmeyer reaction).",
            },
        }

    def _find_substituent(self, query: str):
        """Find substituent in database by any alias."""
        q = query.strip().lower()
        
        # Direct match
        if q in self.database:
            return self.database[q], q

        # Search all display names
        for key, data in self.database.items():
            for name in data["display_names"]:
                if q == name.lower():
                    return data, key
        
        # Partial/fuzzy match
        best_match = None
        best_score = 0
        for key, data in self.database.items():
            for name in data["display_names"]:
                if q in name.lower() or name.lower() in q:
                    score = len(q) if q in name.lower() else len(name)
                    if score > best_score:
                        best_match = (data, key)
                        best_score = score
        
        return best_match

    def _run_base(self, substituent: str, include_swain_lupton: bool = False) -> str:
        """Look up Hammett constants."""
        result = self._find_substituent(substituent)
        
        if result is None:
            # Show available options
            available = sorted([k for k in self.database.keys()])
            return f"## Hammett σ Constant Lookup: `{substituent}`\n\n### ⚠️ Not Found\n\nSubstituent '{substituent}' not found in database.\n\n**Available substituents:**\n" + ", ".join(available[:40]) + f"\n\n...and {len(available)-40} more. Try using format like 'p-NO2', 'm-Cl', 'OCH3', 'CF3', etc."

        data, key = result

        parts = [f"## Hammett σ Constants: {substituent}\n"]
        
        # Also list other names
        aliases = [n for n in data["display_names"] if n.lower() != key]
        if aliases:
            parts.append(f"**Also known as:** {', '.join(aliases[:8])}\n")

        parts.append("### 📊 Hammett σ Values\n")
        parts.append("| Constant | Value | Interpretation |")
        parts.append("|---|---|---|")
        parts.append(f"| **σm (meta)** | **{data['sigma_m']:+.2f}** | {'EWG (deactivating)' if data['sigma_m'] > 0 else 'EDG (activating)' if data['sigma_m'] < 0 else 'Reference'} |")
        parts.append(f"| **σp (para)** | **{data['sigma_p']:+.2f}** | {'EWG (deactivating)' if data['sigma_p'] > 0 else 'EDG (activating)' if data['sigma_p'] < 0 else 'Reference'} |")
        
        if data.get("sigma_p_plus") is not None:
            sp_plus = data["sigma_p_plus"]
            diff_sp = sp_plus - data["sigma_p"]
            parts.append(f"| **σp⁺ (cationic)** | **{sp_plus:+.2f}** | Enhanced {'withdrawal' if diff_sp > 0 else 'donation'} in cationic TS (Δ={diff_sp:+.2f}) |")
        
        if data.get("sigma_p_minus") is not None:
            sp_min = data["sigma_p_minus"]
            diff_sm = sp_min - data["sigma_p"]
            parts.append(f"| **σp⁻ (anionic)** | **{sp_min:+.2f}** | Enhanced {'withdrawal' if diff_sm > 0 else 'donation'} in anionic TS (Δ={diff_sm:+.2f}) |")

        if include_swain_lupton:
            parts.append("\n### 📐 Swain-Lupton Parameters\n")
            parts.append("| Parameter | Value | Meaning |")
            parts.append("|---|---|---|")
            parts.append(f"| **F (field)** | **{data['swain_F']:+.2f}** | Pure inductive effect ({'-I' if data['swain_F'] > 0 else '+I'}) |")
            parts.append(f"| **R (resonance)** | **{data['swain_R']:+.2f}** | Pure resonance effect ({'-R (withdrawing)' if data['swain_R'] > 0 else '+R (donating)'}) |")
            
            total = abs(data['swain_F']) + abs(data['swain_R'])
            if total > 0:
                f_pct = abs(data['swain_F']) / total * 100
                r_pct = abs(data['swain_R']) / total * 100
                parts.append(f"\nEffect breakdown: **{f_pct:.0f}% field (inductive)**, **{r_pct:.0f}% resonance**")

        parts.append(f"\n### 🔬 Classification & Interpretation\n")
        parts.append(f"**Classification:** {data['classification']}\n")
        parts.append(f"{data['interpretation']}\n")

        # Application guidance
        parts.append("### 🧪 Practical Applications\n")
        sm = data["sigma_m"]
        sp = data["sigma_p"]

        # EAS direction
        if sp < sm:
            eas_dir = "**Ortho/Para-directing** (despite possibly being deactivating)"
        elif sp >= sm:
            eas_dir = "**Meta-directing**"
        
        if sm < 0:
            eas_act = "**Activating** (faster than benzene)"
        elif sm < 0.3:
            eas_act = "**Weakly deactivating** (slower than benzene)"
        else:
            eas_act = "**Strongly deactivating** (much slower than benzene)"

        parts.append(f"- **Electrophilic Aromatic Substitution:** {eas_dir}, {eas_act}")
        
        # Acidity effect
        if sp > 0:
            acid_eff = f"Increases acidity of phenols/benzoic acids (ΔpKa ≈ -{abs(sp)*2:.1f} for para-substituted phenol)"
        else:
            acid_eff = f"Decreases acidity of phenols/benzoic acids (ΔpKa ≈ +{abs(sp)*2:.1f} for para-substituted phenol)"
        parts.append(f"- **Acidity Effect:** {acid_eff}")
        
        # SN1/SN2
        if sp > 0.3:
            sn1 = "Accelerates SN1 (stabilizes carbocation)"
        elif sp < -0.2:
            sn1 = "Decelerates SN1 (destabilizes carbocation)"
        else:
            sn1 = "Modest effect on SN1 rates"
        parts.append(f"- **SN1 Reactivity:** {sn1}")

        # Comparison table
        parts.append("\n### 📊 Comparison with Selected Substituents\n")
        parts.append("| Substituent | σm | σp | Classification |")
        parts.append("|---|---|---|---|")
        comparisons = ["h", "ch3", "och3", "f", "cl", "cn", "no2", "cf3", "nh2", "nme2"]
        for cmp_key in comparisons:
            if cmp_key in self.database:
                d = self.database[cmp_key]
                marker = " ← **current**" if cmp_key == key else ""
                parts.append(f"| {cmp_key.upper()} | {d['sigma_m']:+.2f} | {d['sigma_p']:+.2f} | {d['classification']} |{marker}")

        return "\n".join(parts)

    def _run_text(self, input_params: str) -> str:
        input_params = input_params.strip()
        if not input_params:
            raise ChemMCPError("Please provide a substituent name. Example: 'p-NO2', 'm-Cl'")
        return self._run_base(input_params)
