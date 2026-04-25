# DielsAlderReaction (Tool #156)

## Test Input

```python
DielsAlderReaction(**{"diene_smiles": "C=CC=C", "dienophile_smiles": "maleic anhydride", "temperature_c": 80, "solvent": ""})
```

## Output (3310 chars)

```json
{
  "result": {
    "reaction": "1,3-butadiene + maleic anhydride → D-A adduct",
    "diene_analysis": {
      "patterns": [
        "^c=cc=c$",
        "butadiene",
        "1,3-butadiene"
      ],
      "s_cis": true,
      "planar": true,
      "reactivity": "moderate",
      "name": "1,3-butadiene"
    },
    "dienophile_analysis": {
      "patterns": [
        "maleic.anhydride",
        "O=C1OC(=O)C=C1"
      ],
      "ewg": "two carbonyls + anhydride",
      "reactivity": "excellent",
      "name": "maleic anhydride"
    },
    "stereochemistry": {
      "cycloaddition_mode": "[4s+2s] suprafacial (concerted)",
      "stereospecificity": "Yes — configuration of both components is preserved in product",
      "endo_product": "Kinetically favored (Alder endo rule; secondary orbital interactions)",
      "exo_product": "Thermodynamically favored (less steric strain)",
      "predicted_endo_exo_ratio": "Highly dependent on substituents; typically endo-major for cyclic dienophiles (anhydrides, quinones)",
      "chiral_induction_possible": "Yes — via chiral auxiliaries, chiral Lewis acids, or organocatalysis",
      "relative_config": "cis-dienophile → cis-substituted in cyclohexene ring; trans → trans"
    },
    "regiochemistry": {
      "regioselectivity": "N/A (symmetrical components → single product)"
    },
    "product_prediction": {
      "name": "Cycloadduct from 1,3-butadiene + maleic anhydride",
      "ring_system": "cyclohexene derivative",
      "bicyclic": false,
      "endobicyclic": false,
      "new_bonds_formed": "2 σ bonds (C-C) + 1 new π bond (if dienophile was alkene) or 2 π bonds (if alkyne)",
      "degree_of_unsaturation_change": "-1 (two π bonds consumed, one new π formed for alkene dienophile)"
    },
    "feasibility_assessment": {
      "rating": "excellent",
      "score": 85,
      "factors": [
        "Excellent dienophile"
      ]
    },
    "optimal_conditions": {
      "temperature": "80°C (adjust based on reactivity: -78°C to 150°C typical range)",
      "solvent": "neat/toluene",
      "time": "minutes (reactive pairs) to days (unreactive pairs)",
      "atmosphere": "N2/Ar (air-sensitive dienes/dienophiles should be protected)",
      "workup": "Cool, concentrate, purify by column chromatography or recrystallization"
    },
    "substituent_effects": [
      [
        "diene EDG (+OMe, +Me, +NR2)",
        "Rate ↑↑ (raises HOMO)",
        "Normal demand faster"
      ],
      [
        "diene EWG (+COMe, +CN, +NO2)",
        "Rate ↓↓ (lowers HOMO)",
        "May switch to inverse demand"
      ],
      [
        "dienophile EWG (+COMe, +CN, +CHO, +NO2)",
        "Rate ↑↑ (lowers LUMO)",
        "Normal demand faster"
      ],
      [
        "dienophile EDG (+OMe, +NR2)",
        "Rate ↓ (raises LUMO)",
        "Normal demand slower; good for inverse demand"
      ],
      [
        "dienophile = alkyne",
        "Rate slower than alkene",
        "Product: cyclohexadiene (aromatizable)"
      ],
      [
        "dienophile = C=O",
        "Possible but slow",
        "Product: dihydropyran (hetero-D-A)"
      ],
      [
        "dienophile = C≡N",
        "Moderate reactivity",
        "Product: dihydropyridine precursor"
      ]
    ],
    "retro_da_notes": null,
    "summary": "D-A reaction: excellent. "
  }
}
```

## Summary

- **reaction**: 1,3-butadiene + maleic anhydride → D-A adduct
- **diene_analysis**: {'patterns': ['^c=cc=c$', 'butadiene', '1,3-butadiene'], 's_cis': True, 'planar': True, 'reactivity': 'moderate', 'name'...
- **dienophile_analysis**: {'patterns': ['maleic.anhydride', 'O=C1OC(=O)C=C1'], 'ewg': 'two carbonyls + anhydride', 'reactivity': 'excellent', 'nam...
- **stereochemistry**: {'cycloaddition_mode': '[4s+2s] suprafacial (concerted)', 'stereospecificity': 'Yes — configuration of both components i...
- **regiochemistry**: {'regioselectivity': 'N/A (symmetrical components → single product)'}
