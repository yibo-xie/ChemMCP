# WittigReaction (Tool #158)

## Test Input

```python
WittigReaction(**{"carbonyl_smiles": "benzaldehyde", "ylide_type": "non-stabilized", "phosphonium_salt": "Ph3P=CH2", "base": "n-BuLi"})
```

## Output (4017 chars)

```json
{
  "result": {
    "reaction": "benzaldehyde + Ph3P=CH2 (non-stabilized ylide) → benzaldehyde-CH2 alkene (disubstituted terminal/trisubstituted)",
    "carbonyl_analysis": {
      "type": "aldehyde",
      "name": "benzaldehyde",
      "reactivity": "high (good electrophile)",
      "steric_demand": "low"
    },
    "ylide_analysis": {
      "description": "R = alkyl, H (no electron-withdrawing group on ylide carbon)",
      "reactivity": "very high (reacts at -78°C to 0°C)",
      "e_z_selectivity": "Z-alkene favored (kinetic product via betaine pathway)",
      "base_needed": "strong base (n-BuLi, NaHMDS)",
      "examples": [
        "CH2=PPh3 (methylene ylide)",
        "MeCH=PPh3 (ethylidene ylide)"
      ],
      "specified_as": "non-stabilized",
      "salt": "Ph3P=CH2"
    },
    "e_z_selectivity": {
      "prediction": "Terminal alkene (no E/Z isomerism) or Z-favored if disubstituted",
      "rule": "Non-stabilized ylides give kinetic Z-product via betaine pathway"
    },
    "mechanism": [
      {
        "step": "1",
        "name": "Ylide formation",
        "desc": "Ph3P=CR2 formed by deprotonation of phosphonium salt with strong base (n-BuLi, NaH, KHMDS, etc.)"
      },
      {
        "step": "2",
        "name": "[2+2] cycloaddition",
        "desc": "Ylide carbonyl carbon attacks carbonyl C → oxaphosphetane (4-membered ring)"
      },
      {
        "step": "3",
        "name": "Collapse",
        "desc": "Oxaphosphetane collapses: P-O bond forms, C=C bond forms simultaneously"
      },
      {
        "step": "4",
        "name": "Products",
        "desc": "Triphenylphosphine oxide (Ph3P=O) + alkene (R2C=CR'2) — driving force is strong P=O bond (ΔH ~544 kJ/mol)"
      }
    ],
    "product_prediction": {
      "name": "benzaldehyde-CH2 alkene (disubstituted terminal/trisubstituted)",
      "class": "alkene (C=C double bond)",
      "byproduct": "triphenylphosphine oxide (Ph3P=O)",
      "yield": "80-95%"
    },
    "variant_comparison": [
      {
        "variant": "Horner-Wadsworth-Emmons (HWE)",
        "reason": "Better E-selectivity, easier purification (water-soluble byproduct), milder base"
      },
      {
        "variant": "Still-Gennari",
        "reason": "Specifically designed for Z-selective olefination"
      }
    ],
    "scope": [
      "Aldehydes: all types work well (aromatic, aliphatic, α,β-unsaturated)",
      "Ketones: work but slower than aldehydes (steric hindrance); stabilized ylides preferred",
      "Esters: NOT reactive in classical Wittig (use Tebbe or Petasis reagent instead)",
      "Intramolecular Wittig: excellent for cyclic alkene synthesis (macrocycles, medium rings)",
      "α,β-Unsaturated carbonyls: can give 1,2-addition (normal Wittig) or conjugate addition depending on conditions"
    ],
    "limitations": [
      "Sterically hindered ketones react poorly (use Tebbe/Petasis or Julia olefination)",
      "E/Z selectivity can be hard to control for semi-stabilized ylides",
      "Classical Wittig does not work with esters, amides, acid chlorides as electrophiles",
      "Phosphonium salt synthesis requires alkyl halide (may not be readily available)",
      "Ph3P=O byproduct removal can be tedious (column chromatography usually needed)",
      "Base-sensitive functional groups incompatible with ylide formation conditions",
      "Ylides are air/moisture sensitive and malodorous (phosphine smell)"
    ],
    "optimal_conditions": {
      "solvent": "anhydrous THF (standard)",
      "base": "n-BuLi",
      "temperature": "-78°C → RT",
      "atmosphere": "N2/Ar (strictly anhydrous and oxygen-free)",
      "order_of_addition": "Add base to phosphonium salt (form ylide), then add carbonyl solution",
      "workup": "Quench with sat. NH4Cl, extract with Et2O/EtOAc, purify by column chromatography"
    },
    "summary": "Wittig reaction: benzaldehyde-CH2 alkene (disubstituted terminal/trisubstituted). E/Z: Terminal alkene (no E/Z isomerism) or Z-favored if disubstituted. Yield: 80-95%."
  }
}
```

## Summary

- **reaction**: benzaldehyde + Ph3P=CH2 (non-stabilized ylide) → benzaldehyde-CH2 alkene (disubstituted terminal/trisubstituted)
- **carbonyl_analysis**: {'type': 'aldehyde', 'name': 'benzaldehyde', 'reactivity': 'high (good electrophile)', 'steric_demand': 'low'}
- **ylide_analysis**: {'description': 'R = alkyl, H (no electron-withdrawing group on ylide carbon)', 'reactivity': 'very high (reacts at -78°...
- **e_z_selectivity**: {'prediction': 'Terminal alkene (no E/Z isomerism) or Z-favored if disubstituted', 'rule': 'Non-stabilized ylides give k...
- **mechanism**: [{'step': '1', 'name': 'Ylide formation', 'desc': 'Ph3P=CR2 formed by deprotonation of phosphonium salt with strong base...
