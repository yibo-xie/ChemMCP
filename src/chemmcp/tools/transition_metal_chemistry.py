import logging
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPError
from ..utils.mcp_app import ChemMCPManager

logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class TransitionMetalChemistry(BaseTool):
    """
    过渡金属特征化学查询工具。
    覆盖第4-7周期过渡金属（d区元素）的典型性质：变价性、有色化合物、配合物形成、
    催化活性、磁性、d电子组态与颜色/磁性关系等。
    """
    __version__ = "0.1.0"
    name = "TransitionMetalChemistry"
    func_name = "get_transition_metal_chemistry"
    description = "Query transition metal (d-block) chemistry: variable oxidation states, colored compounds, coordination complexes, catalytic activity, magnetic properties, d-electron configurations, and characteristic reactions."
    implementation_description = "Built-in database of transition metal properties organized by element and by concept (oxidation states, colors, magnetism, catalysis). Includes d-d transition color explanation, common coordination geometries, and industrial applications."
    oss_dependencies = []
    services_and_software = []
    categories = ["General"]
    tags = ["Transition Metals", "d-Block", "Coordination Chemistry", "Oxidation States", "Catalysis", "Magnetism"]
    required_envs = []

    code_input_sig = [
        ("element", "str", "N/A", "Element symbol or name (e.g., 'Fe', 'iron', 'Cu', or 'all' for all)."),
        ("property_type", "str", "all", "'oxidation_states', 'colors', 'complexes', 'catalysis', 'magnetism', 'reactions', 'trends', or 'all'."),
    ]

    text_input_sig = [
        ("input_params", "str", "N/A", "Space-separated: 'element [property_type]'. Example: 'Fe oxidation_states' or 'all trends'."),
    ]

    output_sig = [
        ("result", "dict", "Dictionary containing requested data."),
    ]

    examples = [
        {
            "code_input": {"element": "Fe", "property_type": "oxidation_states"},
            "text_input": {"input_params": "Fe oxidation_states"},
            "output": {"result": {"element": "Iron", "common_ox_states": [+2, +3], "stable_ox_states": [] }}
        },
    ]

    DATABASE = {
        "Ti": {
            "name": "Titanium", "period": 4, "group": 4,
            "electron_config": "[Ar] 3d² 4s²",
            "common_ox_states": [+2, +3, +4], "most_stable": +4,
            "characteristic_color": "Ti(IV) compounds are white/colorless; Ti(III) is violet/red (e.g., TiCl₃ is purple)",
            "key_compounds": ["TiO₂ (rutile/anatase — white pigment, sunscreen)", "TiCl₄ (colorless liquid, precursor to Ti metal via Kroll process)", "BaTiO₃ (ferroelectric)"],
            "key_complexes": ["[Ti(H₂O)₆]³⁺ (violet — d¹, one d-d transition)", "[TiF₆]²⁻"],
            "catalysis": ["Ziegler-Natta catalysts (TiCl₄ + Al(C₂H₅)₃ for polypropylene polymerization)"],
            "biology": None,
            "applications": ["Aerospace alloys (lightweight, strong, corrosion-resistant)", "Pigments (TiO₂ — most used white pigment)", "Medical implants (biocompatible)"],
        },
        "V": {
            "name": "Vanadium", "period": 4, "group": 5,
            "electron_config": "[Ar] 3d³ 4s²",
            "common_ox_states": [+2, +3, +4, +5], "most_stable": "+5 (as VO₂⁺ / VO₄³⁻)",
            "characteristic_color": "V(II) purple; V(III) green; V(IV) blue (VO²⁺ vanadyl); V(V) yellow (VO₃⁻ / VO₄³⁻)",
            "key_compounds": ["V₂O₅ (orange solid — important industrial catalyst)", "NH₄VO₃ (ammonium metavanadate)", "NaVO₃"],
            "key_complexes": ["[V(H₂O)₆]²⁺ (violet), [V(H₂O)₆]³⁺ (green), [VO(H₂O)₅]²⁺ (blue), [VO₄]³⁻ (yellow)"],
            "catalysis": ["V₂O₅ in Contact Process for H₂SO₄ production (oxidizes SO₂ → SO₃)"],
            "biology": "Some enzymes contain V (vanadium nitrogenases in certain bacteria)",
            "applications": ["Steel alloying (tool steels — increases strength/hardness)", "Redox flow batteries (vanadium redox battery)", "Catalyst (sulfuric acid production)"],
        },
        "Cr": {
            "name": "Chromium", "period": 4, "group": 6,
            "electron_config": "[Ar] 3d⁵ 4s¹",
            "common_ox_states": [+1, +2, +3, +6], "most_stable": "+3 (kinetically inert); +6 (thermodynamically stable but strong oxidizer)",
            "characteristic_color": "Cr(III) green/violet; Cr(VI) yellow/orange (CrO₄²⁻ yellow, Cr₂O₇²⁻ orange); Cr(II) blue",
            "key_compounds": ["K₂Cr₂O₇ (potassium dichromate — orange, strong oxidizer in acid)", "Na₂CrO₄ (chromate — yellow)", "CrCl₃ (green)", "CrO₃ (red, dangerous)"],
            "key_complexes": ["[Cr(H₂O)₆]³⁺ (violet), [Cr(NH₃)₆]³⁺ (yellow), [Cr(CN)₆]³⁻", "[Cr(edta)]⁻"],
            "catalysis": ["Various Cr-catalyzed organic oxidations"],
            "biology": "Glucose tolerance factor contains Cr(III); essential trace element for sugar metabolism",
            "applications": ["Stainless steel plating (chrome plating — shiny, corrosion-resistant)", "Leather tanning (Cr(III) salts cross-link proteins)", "Pigments (chrome green/yellow)", "Refractory materials (chromite bricks)"],
            "toxicity": "Cr(VI) compounds are HIGHLY TOXIC (carcinogenic, mutagenic); Cr(III) is an essential nutrient",
        },
        "Mn": {
            "name": "Manganese", "period": 4, "group": 7,
            "electron_config": "[Ar] 3d⁵ 4s²",
            "common_ox_states": [+2, +3, +4, +6, +7], "most_stable": "+2 (aqueous); +7 (permanganate — strong oxidizer)",
            "characteristic_color": "Mn(II) pale pink (very pale!); Mn(IV) black/brown (MnO₂); Mn(VII) deep purple (MnO₄⁻)",
            "key_compounds": ["KMnO₄ (potassium permanganate — deep purple, powerful oxidizer/disinfectant)", "MnO₂ (black, cathode in dry-cell batteries)", "KMnO₄ (purple)",
                              "MnSO₄ (pink solution)"],
            "key_complexes": ["[Mn(H₂O)₆]²⁺ (pale pink — d⁵ high-spin, spin-forbidden d-d transitions = very weak color)"],
            "catalysis": ["MnO₂ as catalyst in decomposition of H₂O₂ (oxygen generation in labs)"],
            "biology": "Essential in photosystem II water-splitting complex (Mn₄CaO₅ cluster); also essential enzyme cofactor (Mn-SOD)",
            "applications": ["Steel making (ferromanganese — deoxidizes and desulfurizes steel)", "Batteries (alkaline batteries use MnO₂ cathode)", "Water treatment (KMnO₄ disinfectant/oxidizer)", "Fertilizer (micronutrient)"],
        },
        "Fe": {
            "name": "Iron", "period": 4, "group": 8,
            "electron_config": "[Ar] 3d⁶ 4s²",
            "common_ox_states": [+2, +3], "less_common": [+6],
            "most_stable": "+3 (aerobic); +2 (anaerobic/reducing conditions)",
            "characteristic_color": "Fe(II) pale green (oxidizes to Fe(III)); Fe(III) yellow/brown (often rust-colored from hydrolysis); Fe(VI) dark purple (ferrate)",
            "key_compounds": ["Fe₂O₃ (hematite — red, iron ore/rust)", "Fe₃O₄ (magnetite — black, magnetic)", "FeSO₄·7H₂O (green vitriol)", "FeCl₃ (brown/yellow)",
                              "K₂FeO₄ (potassium ferrate — purple, strong oxidizer)"],
            "catalysis": ["Fe in Haber-Bosch process (promoter)", "Fenton's reagent (Fe²⁺ + H₂O₂ → •OH radical for organic degradation)", "Heme enzymes (cytochrome P450, catalase, peroxidase)"],
            "biology": "CRITICAL: hemoglobin (O₂ transport — Fe in heme), myoglobin (O₂ storage), cytochromes (ETC), ferritin (Fe storage), iron-sulfur clusters",
            "applications": ["Steel production (98% of Fe consumption goes to steel)", "Catalysis (various industrial processes)", "Medicine (iron supplements, contrast agents)", "Magnetic materials (Fe₃O₄, Fe₃C in magnets)"],
        },
        "Co": {
            "name": "Cobalt", "period": 4, "group": 9,
            "electron_config": "[Ar] 3d⁷ 4s²",
            "common_ox_states": [+2, +3], "rare": [+4],
            "most_stable": "+2 (aqueous); +3 (in complexes like [Co(NH₃)₆]³⁺)",
            "characteristic_color": "Co(II) pink (hydrated); Co(III) complexes vary widely (depends on ligands)",
            "key_compounds": ["CoCl₂ (anhydrous = blue; hydrated = pink — humidity indicator!)", "Co₃O₄ (black)", "Co(NO₃)₂"],
            "key_complexes": ["[Co(H₂O)₆]²⁺ (pink), [Co(NH₃)₆]²⁺ (rapidly oxidized to Co(III)), [Co(NH₃)₆]³⁺ (yellow-orange, low-spin d⁶)", "[Co(CN)₆]³⁻ (yellow), [Co(en)₃]³⁺ (optical isomers! — first resolved chiral complex)"],
            "catalysis": ["Hydroformylation (Co carbonyl catalysts — historical)", "Pincer complexes for hydrogenation"],
            "biology": "Vitamin B₁₂ (cobalamin — only known biological role of Co; coenzyme for methyl transfer reactions)",
            "applications": ["Superalloys (jet engines, gas turbines — high temperature strength)", "Li-ion battery cathodes (LCO, NCA)", "Magnets (Alnico, Sm-Co)", "Pigments (cobalt blue — Thenard's blue)", "Radiotherapy (Co-60 gamma source)"],
        },
        "Ni": {
            "name": "Nickel", "period": 4, "group": 10,
            "electron_config": "[Ar] 3d⁸ 4s²",
            "common_ox_states": [+2], "less_common": [+1, +3, +4],
            "most_stable": "+2",
            "characteristic_color": "Ni(II) green (hydrated); many Ni(II) complexes are colored (d⁸ has allowed transitions)",
            "key_compounds": ["NiSO₄ (green)", "Ni(OH)₂ (green)", "NiO (green)", "Ni(CO)₄ (volatile liquid — Mond process for Ni purification)"],
            "key_complexes": ["[Ni(H₂O)₆]²⁺ (green), [Ni(NH₃)₆]²⁺ (blue/purple), [Ni(CN)₄]²⁻ (YELLOW — square planar! unusual for d⁸)", "[Ni(en)₃]²⁺ (purple — octahedral)"],
            "catalysis": ["Raney Ni (hydrogenation catalyst)", "Ni in steam reforming (methane → syngas)", "Cross-coupling (Negishi, Kumada couplings with Ni catalysts)"],
            "biology": "Several Ni-containing enzymes (urease, hydrogenase, methyl-coenzyme M reductase in methanogens)",
            "applications": ["Stainless steel (austenitic SS contains 8-10% Ni)", "Coins (US nickel actually 25% Ni)", "Batteries (Ni-Cd, Ni-MH)", "Electroplating (corrosion-resistant coating)", "Catalyst (hydrogenation, reforming)"],
        },
        "Cu": {
            "name": "Copper", "period": 4, "group": 11,
            "electron_config": "[Ar] 3d¹⁰ 4s¹",
            "common_ox_states": [+1, +2], "rare": [+3],
            "most_stable": "+2 (aqueous); +1 (solid/copper(I) compounds)",
            "characteristic_color": "Cu metal: reddish-brown; Cu(I): colorless/white; Cu(II): BLUE (characteristic! — Jahn-Teller distorted d⁹)",
            "key_compounds": ["CuSO₄·5H₂O (blue vitriol — blue crystals, algicide)", "Cu₂O (red — cuprous oxide)", "CuO (black)", "CuCl₂ (green-brown)",
                              "Cu(NO₃)₂ (blue)"],
            "key_complexes": ["[Cu(H₂O)₆]²⁺ (blue — classic JT elongated octahedron), [Cu(NH₃)₄(H₂O)₂]²⁺ (deep royal blue — Tollens'-like test)", "[CuCl₄]²⁻ (yellow-green), [Cu(edta)]²⁻ (deep blue)"],
            "catalysis": ["Cu in Wacker process (ethene → acetaldehyde)", "Cu in Ullmann coupling", "Glucose oxidase (contains Cu)"],
            "biology": "Cytochrome c oxidase (terminal ETC enzyme — reduces O₂ to H₂O), superoxide dismutase (Cu,Zn-SOD), hemocyanin (arthropod/mollusk O₂ transport)",
            "applications": ["Electrical wiring (second best conductor after Ag; much cheaper)", "Plumbing (pipes, fittings)", "Coinage (US penny = Cu-plated Zn)", "Antimicrobial surfaces (Cu alloys kill bacteria/viruses)", "Roofs/architectural (patina formation)"],
        },
        "Zn": {
            "name": "Zinc", "period": 4, "group": 12,
            "electron_config": "[Ar] 3d¹⁰ 4s²",
            "common_ox_states": ["+2 ONLY"], "note": "Full d-subshell — NOT a true transition metal by IUPAC definition (no partially filled d orbitals in any common state)",
            "characteristic_color": "ALL Zn(II) compounds are WHITE/COLORLESS (d¹⁰ — no d-d transitions possible!)",
            "key_compounds": ["ZnO (white powder — sunblock, rubber additive)", "ZnS (zinc blende — phosphor, luminescent when doped)", "ZnCl₂ (deliquescent white solid — flux/galvanizing)"],
            "key_complexes": ["[Zn(H₂O)₆]²⁺ (colorless), [Zn(NH₃)₄]²⁺ (colorless), [Zn(CN)₄]²⁻ (colorless)"],
            "catalysis": ["Zn enzymes (carbonic anhydrase, alcohol dehydrogenase, zinc finger proteins)"],
            "biology": "ESSENTIAL: >300 enzymes require Zn (carbonic anhydrase, DNA/RNA polymerases, zinc finger transcription factors, alcohol dehydrogenase)",
            "applications": ["Galvanization (Zn coating on steel prevents rusting)", "Die-casting alloys (brass = Cu-Zn)", "Batteries (Zn-carbon, Zn-air)", "Sunscreen (ZnO nanoparticles — physical UV blocker)", "Nutritional supplements"],
        },
        "Ag": {
            "name": "Silver", "period": 5, "group": 11,
            "electron_config": "[Kr] 4d¹⁰ 5s¹",
            "common_ox_states": [+1], "less_common": [+2, +3],
            "most_stable": "+1",
            "characteristic_color": "Ag metal: white metallic luster; Ag(I) compounds usually colorless (d¹⁰); AgBr pale yellow; AgI yellow",
            "key_compounds": ["AgNO₃ (silver nitrate — colorless, soluble, used in photography/halide detection)", "AgCl (white precipitate — darkens on exposure to light)",
                              "AgBr (pale yellow — photographic film emulsion)", "Ag₂O (dark brown — decomposes to Ag + O₂ on heating)"],
            "key_complexes": ["[Ag(NH₃)₂]⁺ (Tollens' reagent — diammine silver(I), colorless, mild oxidizer for aldehydes)", "[Ag(CN)₂]⁻ (used in silver electroplating), [Ag(S₂O₃)₂]³⁻ (photographic fixing)"],
            "catalysis": ["Ag catalyst in ethylene oxide production (historical)", "Silver nanoparticles (antibacterial)"],
            "biology": "No known biological role; Ag⁺ is antimicrobial (disrupts bacterial enzymes/membranes)",
            "applications": ["Jewelry/coinage (sterling silver = 92.5% Ag)", "Photography (historical — Ag halides light-sensitive)", "Electronics (high-performance conductors, contacts)", "Antimicrobial (wound dressings, coatings, water filters)", "Solar panels (Ag paste in PV cells)"],
            "tarnishing": "Ag + H₂S → Ag₂S (black) + H₂ — tarnish reaction with atmospheric H₂S",
        },
        "Pt": {
            "name": "Platinum", "period": 6, "group": 10,
            "electron_config": "[Xe] 4f¹⁴ 5d⁹ 6s¹",
            "common_ox_states": [+2, +4], "also": ["0 (in complexes)"],
            "most_stable": "+4 (thermodynamically); kinetically very inert (both Pt(II) and Pt(IV))",
            "characteristic_color": "Pt metal: silvery-white; Pt(II) complexes often colored (square planar d⁸); Pt(IV) usually less colored",
            "key_compounds": ["H₂PtCl₆ (chloroplatinic acid — yellow, precursor to most Pt compounds)", "PtCl₂", "PtO₂"],
            "key_complexes": ["[PtCl₄]²⁻ (square planar)", "[Pt(NH₃)₄]²⁺ (colorless)", "cis-[PtCl₂(NH₃)₂] (cisplatin — anticancer drug! YELLOW)",
                           "trans-[PtCl₂(NH₃)₂] (transplatin — inactive, shows importance of geometry)"],
            "catalysis": ["Automotive catalytic converters (Pt/Pd/Rh — oxidize CO, reduce NOx)", "Pt in fuel cell electrodes (ORR catalyst)",
                          "Hydroation/silane addition (Speier's catalyst — H₂PtCl₆)", "PtO₂ (Adams' catalyst — hydrogenations)"],
            "biology": "None naturally; but cisplatin is a landmark anticancer drug (cross-links DNA)",
            "applications": ["Catalytic converters (largest use of Pt)", "Jewelry", "Laboratory equipment (electrodes, crucibles)", "Anticancer drugs (cisplatin, carboplatin)", "Thermocouples (Pt-Rh)"],
        },
        "Au": {
            "name": "Gold", "period": 6, "group": 11,
            "electron_config": "[Xe] 4f¹⁴ 5d¹⁰ 6s¹",
            "common_ox_states": [+1, +3], "rare": ["-1 (aurides like CsAu)"], "most_common_in_nature": "0 (native gold)",
            "most_stable": "+3 (in compounds); +1 (disproportionates unless stabilized)",
            "characteristic_color": "Au metal: distinctive YELLOW (relativistic effects contract 6s orbital, raise 5d energy → absorb blue light)",
            "key_compounds": ["HAuCl₄ (chloroauric acid — yellow)", "AuCl₃ (red, dimeric Au₂Cl₆)", "Au₂O₃ (dark brown)"],
            "key_complexes": ["[Au(CN)₂]⁻ (cyanidation of gold ore extraction — dissolves Au)", "[AuCl₄]⁻ (square planar)"],
            "catalysis": ["Au nanoparticles (selective oxidation, CO oxidation at low T)", "Au-catalyzed alkyne hydration"],
            "biology": "None; Au is biologically inert (one reason it's used in dentistry/implants)",
            "applications": ["Monetary reserves/jewelry", "Electronics (corrosion-free contacts, bond wires in chips)", "Dentistry (crowns, inlays — biocompatible)", "Auric cyanidation (gold extraction from ore)", "Nanomedicine (AuNP drug delivery, photothermal therapy)", "Space (gold foil on spacecraft IR reflector)"],
            "special_note": "Gold is the MOST NOBLE metal — does not react with O₂, H₂O, acids (including non-oxidizing acids); only dissolved by aqua regia (HCl:HNO₃ 3:1)",
        },
    }

    CONCEPT_DATA = {
        "why_colored": "Transition metal compounds are colored due to **d-d electronic transitions**. In a free ion, all 5 d orbitals are degenerate. In a ligand field (crystal field splitting), they split into different energy levels. Electrons can absorb visible light to jump between these split levels. The complementary color of absorbed light is what we see.\n\nKey factors:\n- d⁰ (Sc³⁺, Ti⁴⁺) and d¹⁰ (Zn²⁺, Ag⁺): NO d-d transitions → COLORLESS\n- d¹-d⁹: colored (except some symmetric cases)\n- Intensity depends on whether transitions are Laporte-allowed (weak in centrosymmetric complexes)\n- Color also affected by charge-transfer bands (often more intense than d-d)",
        "why_variable_oxidation": "Transition metals have similar energies for (n-1)d and ns electrons → both can participate in bonding → multiple oxidation states accessible. The range of oxidation states is widest in the middle of each transition series (e.g., Mn: +2 to +7) and narrows toward the ends.",
        "why_form_complexes": "Small, highly charged metal ions have: (1) vacant low-energy orbitals to accept lone pairs, (2) high charge density to polarize/attract electron donors. This makes them excellent Lewis acids that form coordinate covalent bonds with Lewis bases (ligands).",
        "catalytic_activity_reasons": "(1) Variable oxidation states allow redox cycles (metal shuttles between states during reaction), (2) Ability to adsorb reactants onto surface or into coordination sphere, (3) Large surface area (nanoparticles), (4) Moderate M-L bond strengths (strong enough to bind reactants, weak enough to release products).",
    }

    def __init__(self, init: bool = True, interface: str = "code"):
        super().__init__(init=init, interface=interface)

    def _init_modules(self):
        pass

    def _run_base(self, element: str, property_type: str = "all") -> dict:
        element = element.strip().capitalize()
        prop_type = property_type.lower().strip()

        if prop_type == "trends":
            return {"result": self._get_trends()}

        if prop_type in self.CONCEPT_DATA:
            return {"result": {prop_type: self.CONCEPT_DATA[prop_type]}}

        if element == "All":
            result = {}
            for sym in self.DATABASE:
                result[sym] = self._filter_data(self.DATABASE[sym], prop_type)
            return {"result": result}

        if element not in self.DATABASE:
            raise ChemMCPError(f"Element '{element}' not found. Options: {list(self.DATABASE.keys()) + ['All']}")

        data = self.DATABASE[element]
        return {"result": {**{"element": element}, **self._filter_data(data, prop_type)}}

    def _filter_data(self, data: dict, prop_type: str) -> dict:
        if prop_type == "all":
            return {k: v for k, v in data.items() if k not in ("period", "group")}
        elif prop_type in data:
            return {prop_type: data.get(prop_type)}
        else:
            keys = ("name", "electron_config", "common_ox_states", "most_stable", "characteristic_color",
                     "key_compounds", "key_complexes", "catalysis", "biology", "applications")
            return {k: data.get(k) for k in keys if k in data}

    def _get_trends(self) -> dict:
        return {
            "across_period_4": {
                "atomic_radius": "Decreases then levels off (Sc > Ti > V > Cr ≈ Mn < Fe < Co < Ni < Cu > Zn — d-block contraction)",
                "ionization_energy": "Generally increases across (with dips at Cr [3d⁵4s¹] and Cu [3d¹⁰4s¹])",
                "number_of_ox_states": "Increases to middle (Mn has most: +2 to +7), then decreases toward end (Zn only +2)",
                "ionic_radius_M2p": "Decreases: Mn²⁺(83) > Fe²⁺(78) > Co²⁺(74.5) > Ni²⁺(69) > Cu²⁺(73) > Zn²⁺(74) pm (anomaly at Cu/Zn due to Jahn-Teller/d¹⁰)",
                "stability_of_higher_ox_states": "Increases across period (Ti(+4) stable, Cr(+6) strong oxidizer, Fe(+6) very strong oxidizer)",
                "color_intensity": "Generally increases toward middle of series (more unpaired electrons → more transitions)",
                "complex_stability": "Irregular trend ( Irving-Williams series for M²⁺: Mn²⁺ < Fe²⁺ < Co²⁺ < Ni²⁺ < Cu²⁺ >> Zn²⁺ )",
                "catalytic_importance": "Fe, Co, Ni, Pd, Pt are most industrially important",
            },
            "irving_williams_series": "Complex stability (for identical ligands): Mn²⁺ < Fe²⁺ < Co²⁺ < Ni²⁺ < Cu²⁺ >> Zn²⁺\nReason: combination of ionic radius (decreases left→right) and LFSE (varies irregularly)",
            "color_rules": [
                "d⁰ and d¹⁰ ions form colorless compounds (no d-d transitions): Sc³⁺, Ti⁴⁺, Zn²⁺, Ag⁺, Cd²⁺, Hg²⁺",
                "d⁵ high-spin (e.g., Mn²⁺, Fe³⁺) are VERY pale (spin-forbidden transitions): Mn²⁺ barely pink, Fe³⁺ pale violet",
                "d⁹ (Cu²⁺) is ALWAYS blue/green (Jahn-Teller distortion guarantees asymmetry)",
                "Charge-transfer (LMCT/MLCT) bands are often MORE intense than d-d bands (e.g., permanganate MnO₄⁻ purple from LMCT, not d-d)",
            ],
            "special_notes": [
                "Sc-Y-La-Lu group: more like main-group metals (dominant +3 ox state, few complexes, no colored compounds)",
                "Group 12 (Zn, Cd, Hg): d¹⁰ configuration — NOT true transition metals by IUPAC (no partial d shell in common states)",
                "Relativistic effects become significant for 5d/6d elements (Au yellow color, Hg liquid at RT, Pt/Au nobility)",
            ]
        }

    def _run_text(self, input_params: str) -> dict:
        try:
            parts = input_params.strip().split()
            elem = parts[0] if parts else "All"
            prop = parts[1] if len(parts) > 1 else "all"
            return self._run_base(elem, prop)
        except Exception as e:
            raise ChemMCPError(f"Failed to parse: {str(e)}. Format: 'element [property_type]'")
