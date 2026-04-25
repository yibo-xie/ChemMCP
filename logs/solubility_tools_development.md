# ChemMCP Solubility Tools (#66-70) Development Log

## Overview
Developed 5 new MCP tools for chemical equilibrium and solubility calculations,
completing the #61-70 tool set in the ChemMCP project.

## Tool Summary

| # | Tool Name | Function | Purpose |
|---|-----------|----------|---------|
| 66 | SelectivePrecipitation | `selective_precipitation` | Design selective precipitation separation schemes based on classic qualitative analysis (Groups I-V) |
| 67 | DissolvePrecipitate | `dissolve_precipitate` | Analyze dissolution conditions for precipitates (acid/base/complexation/redox) |
| 68 | CommonIonSolubility | `common_ion_solubility` | Calculate common ion effect on solubility (50+ Ksp database) |
| 69 | ComplexIonSolubility | `complex_ion_solubility` | Calculate complex ion formation effect on solubility (20+ Kf database) |
| 70 | SolubilityRules | `solubility_rules` | Query standard solubility rules for qualitative prediction (9 rule groups, 100+ compounds) |

## Core Logic Verification

### #66 SelectivePrecipitation
- **Algorithm**: Classic cation group separation (H₂S system)
  - Group I: HCl → Ag⁺, Pb²⁺, Hg₂²⁺ (chloride insoluble)
  - Group II: H₂S (0.3M H⁺) → Cu²⁺, Cd²⁺, Pb²⁺ (acidic sulfide)
  - Group III: (NH₄)S (basic) → Fe³⁺, Al³⁺, Cr³⁺ (basic sulfide/hydroxide)
  - Group IV: (NH₄)₂CO₃ → Ba²⁺, Sr²⁺, Ca²⁺ (carbonate)
  - Group V: soluble (Na⁺, K⁺, NH₄⁺, Mg²⁺)
- **Test**: Ag⁺/Ba²⁺/Cu²⁺/Mg²⁺ → 3-step scheme ✅
- **pH tracking**: Each step includes pH range for optimal precipitation

### #67 DissolvePrecipitate
- **4 dissolution mechanisms**:
  1. Acid dissolution (protonation of anion): carbonates, sulfides, hydroxides
  2. Base dissolution (amphoteric): Al(OH)₃, Zn(OH)₂, Cr(OH)₃, Pb(OH)₂
  3. Complexation: AgCl + 2NH₃ → [Ag(NH₃)₂]⁺ + Cl⁻
  4. Redox: CuS + oxidants (HNO₃/H₂O₂/aqua regia)
- **Database**: 30+ precipitates with Ksp, Kf, minimum concentrations
- **Test**: AgCl→3 methods (NH₃/Na₂S₂O₃/Na₂S₂O₃+KCN), CuS→redox ✅

### #68 CommonIonSolubility
- **Formula**: For AmBn(s) ⇌ mAⁿ⁺ + nBᵐ⁻, Ksp = [A]ᵐ[B]ⁿ
  - With common ion [B] = C: S = (Ksp/Cⁿ)^(1/m) for m=1; or solve polynomial
- **Database**: 50+ compounds with Ksp, stoichiometry (m,n), Mw
- **Test**: AgCl + 0.1M Cl⁻ → S reduced 7516× (100% decrease) ✅
- **Output**: Pure water solubility, ion-effect solubility, reduction factor, percent decrease, ion concentrations

### #69 ComplexIonSolubility
- **Formula**: Overall K = Ksp × Kf (formation constant of complex)
  - Enhanced solubility: S' ≈ (K × [L]^n)^(1/m) where n = coordination number
- **Database**: 20+ metal-ligand pairs with Kf values
  - Metals: Ag⁺, Cu²⁺, Zn²⁺, Al³⁺, Fe³⁺, Hg²⁺, Ni²⁺, Co²⁺, Pb²⁺, Cd²⁺, etc.
  - Ligands: NH₃, CN⁻, OH⁻, S₂O₃²⁻, Cl⁻, I⁻, SCN⁻, C₂O₄²⁻, en, EDTA, CH₃COO⁻, F⁻
- **Test**: AgCl + 1M NH₃ → 4123× enhancement → [Ag(NH₃)₂]⁺ ✅
- **Test**: Cu(OH)₂ + 2M NH₃ → 110357× enhancement → [Cu(NH₃)₄]²⁺ ✅

### #70 SolubilityRules
- **9 rule groups**: Nitrates/Acetates, Halides, Sulfates, Carbonates/Phosphates, Hydroxides, Sulfides, Chromates, Oxalates, Fluorides
- **Compound parser**: Extracts cation/anion from formula (case-insensitive matching)
- **Rule engine**: Always-soluble ions → Anion-specific rules → Exception checking → Slightly-soluble classification
- **Test**: NaCl=soluble, BaSO4=insoluble, AgF=soluble(exception), CaCO3=insoluble ✅

## Files Created/Modified

### New Files:
- `src/chemmcp/tools/selective_precipitation.py` (#66)
- `src/chemmcp/tools/dissolve_precipitate.py` (#67)
- `src/chemmcp/tools/common_ion_solubility.py` (#68)
- `src/chemmcp/tools/complex_ion_solubility.py` (#69)
- `src/chemmcp/tools/solubility_rules.py` (#70)
- `tests/test_solubility_tools.py` (16 test cases)

### Modified Files:
- `src/chemmcp/tools/__init__.py` (added #61-70 registration entries)

## Test Results
```
✅ 16/16 tests passed
✅ All 10 tools (#61-70) import successfully
✅ 21/21 functional smoke tests passed
```

## Cherry Studio Import Configuration
See: `logs/cherry_studio_config.json`
