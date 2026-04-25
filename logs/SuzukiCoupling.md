# SuzukiCoupling (Tool #160)

## Test Input

```python
SuzukiCoupling(**{"organoboron_smiles": "phenylboronic acid", "halide_smiles": "bromobenzene", "ligand": "PPh3", "base": "K2CO3", "solvent": "dioxane/H2O"})
```

## Output (4890 chars)

```json
{
  "result": {
    "reaction": "Suzuki coupling: phenylboronic acid + bromobenzene → cross-coupled product",
    "organoboron_analysis": {
      "type": "boronic acid",
      "input": "phenylboronic acid",
      "note": "Moderately stable; protodeboronation possible"
    },
    "halide_analysis": {
      "type": "aryl bromide",
      "input": "bromobenzene",
      "note": "Standard substrate; good cost/reactivity balance"
    },
    "catalytic_cycle": [
      [
        "1",
        "Oxidative addition",
        "Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X) complex. Rate: RI > RBr >> RCl (activated) >> RF. Vinyl/Ar halides faster than alkyl."
      ],
      [
        "2",
        "Transmetalation",
        "Organoboron (R²-B(OH)2 or ester) activated by base → transfers R² to Pd, displacing X. Base forms boronate [R²-B(OH)3]⁻ which transmetalates more readily."
      ],
      [
        "3",
        "Reductive elimination",
        "Pd(II)(R¹)(R²) eliminates R¹-R² coupled product + regenerates Pd(0). Stereoretentive for stereochemical centers on both partners."
      ]
    ],
    "mechanism_steps": [
      {
        "step": "1",
        "name": "Oxidative addition",
        "detail": "Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X) complex. Rate: RI > RBr >> RCl (activated) >> RF. Vinyl/Ar halides faster than alkyl."
      },
      {
        "step": "2",
        "name": "Transmetalation",
        "detail": "Organoboron (R²-B(OH)2 or ester) activated by base → transfers R² to Pd, displacing X. Base forms boronate [R²-B(OH)3]⁻ which transmetalates more readily."
      },
      {
        "step": "3",
        "name": "Reductive elimination",
        "detail": "Pd(II)(R¹)(R²) eliminates R¹-R² coupled product + regenerates Pd(0). Stereoretentive for stereochemical centers on both partners."
      }
    ],
    "catalyst_recommendation": {
      "primary": "Pd(PPh3)4 (2-3 mol%) or Pd(OAc)2/PPh3 (3 mol%/6 mol%)",
      "reason": "Standard catalyst system for aryl bromides/iodides"
    },
    "optimal_conditions": {
      "catalyst_loading": "1-3 mol% Pd (standard)",
      "ligand_loading": "2-4 mol% (relative to Pd)",
      "base": "K2CO3 (2-3 eq)",
      "solvent": "dioxane/H2O",
      "concentration": "0.1-0.5 M",
      "temperature": "60-80°C",
      "time": "2-16 hours",
      "atmosphere": "N2 or Ar (DEGASSED solvents — oxygen causes homocoupling!)",
      "monitoring": "TLC or GC/MS for consumption of limiting reagent",
      "workup": "Dilute with water, extract with EtOAc (×3), wash with brine, dry (Na2SO4), concentrate, purify by column chromatography",
      "palladium_removal": "For pharmaceutical: pass through silica-bound thiol scavenger resin or aqueous KF wash to reduce Pd to <10 ppm"
    },
    "scope": [
      "Aryl-aryl (biaryl) couplings: EXCELLENT — the flagship transformation",
      "Aryl-vinyl (styrene) couplings: Excellent — stereochemistry retained",
      "Vinyl-vinyl (diene) couplings: Good — stereochemistry retained",
      "Heteroaryl couplings: Good — pyridine, thiophene, furan, indole all work (may need optimization)",
      "Alkyl-aryl couplings: CHALLENGING — primary alkyl halides work with modern catalysts; secondary/tertiary difficult",
      "Couplings with ortho-substituted partners: Works with bulky ligands (SPhos, XPhos)",
      "Large-scale: Industrially proven (e.g., boscalid, valsartan synthesis)"
    ],
    "applicable_limitations": [
      {
        "issue": "Protodeboronation",
        "problem": "Boronic acids can decompose under basic conditions → Ar-H byproduct",
        "solution": "Use boronic ester (Bpin) instead; lower temperature; minimize reaction time"
      },
      {
        "issue": "Oxygen sensitivity",
        "problem": "O2 causes oxidative homocoupling of boronic species",
        "solution": "Degas solvents thoroughly; sparge with inert gas; use sealed tube"
      }
    ],
    "ligand_effects": {
      "electron_rich_bulky_phosphines": [
        "Accelerate oxidative addition (especially for Ar-Cl)",
        "Prevent β-hydride elimination",
        "Examples: SPhos, XPhos, RuPhos, DavePhos, BrettPhos, JohnPhos"
      ],
      "bidentate_phosphines": [
        "Stabilize Pd intermediates; prevent Pd black formation",
        "Examples: dppf, dppe, dppp, Xantphos"
      ],
      "N_heterocyclic_carbenes_NHC": [
        "Very strong σ-donors; highly active; air-stable complexes",
        "Examples: IPr, IMes, PEPPSI series"
      ],
      "water_soluble_ligands": [
        "TPPTS (triphenylphosphine trisulfonate); enable aqueous-phase Suzuki",
        "Green chemistry advantage"
      ]
    },
    "typical_yields": "80-95%",
    "base_role": "Activates organoboron compound via formation of tetracoordinate boronate [R-B(OH)3]⁻ which undergoes transmetalation readily",
    "summary": "Suzuki coupling predicted: 80-95%. Key concern: Protodeboronation"
  }
}
```

## Summary

- **reaction**: Suzuki coupling: phenylboronic acid + bromobenzene → cross-coupled product
- **organoboron_analysis**: {'type': 'boronic acid', 'input': 'phenylboronic acid', 'note': 'Moderately stable; protodeboronation possible'}
- **halide_analysis**: {'type': 'aryl bromide', 'input': 'bromobenzene', 'note': 'Standard substrate; good cost/reactivity balance'}
- **catalytic_cycle**: [('1', 'Oxidative addition', 'Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X) complex. Rate: RI > RBr >> RCl (activated) >>...
- **mechanism_steps**: [{'step': '1', 'name': 'Oxidative addition', 'detail': 'Pd(0) inserts into R¹-X bond → Pd(II)(R¹)(X) complex. Rate: RI >...
