# AldolReaction (Tool #154)

## Test Input

```python
AldolReaction(**{"substrate1_smiles": "CC=O", "substrate2_smiles": "", "catalyst_type": "base", "solvent": "EtOH"})
```

## Output (3076 chars)

```json
{
  "result": {
    "reaction_type": "Self-aldol of CC=O",
    "substrate1_analysis": {
      "alpha_h_pKa": 17,
      "enolization": "easy",
      "electrophile_reactivity": "high",
      "self_aldol_tendency": "high",
      "examples": [
        "acetaldehyde (CH3CHO)",
        "propionaldehyde (CH3CH2CHO)",
        "butyraldehyde (CH3(CH2)2CHO)"
      ],
      "classified_as": "aliphatic_aldehyde",
      "input": "CC=O"
    },
    "substrate2_analysis": null,
    "catalyst_analysis": {
      "selected": "generic base",
      "type": "base",
      "notes": "Base-type catalyst specified: base. See NaOH/KOH or NaOEt for typical conditions."
    },
    "product_prediction": {
      "aldol_adduct": "β-hydroxy carbonyl dimer from CC=O",
      "condensation_product": "α,β-unsaturated carbonyl from CC=O",
      "functional_groups_adduct": [
        "carbonyl",
        "secondary alcohol"
      ],
      "functional_groups_condensed": [
        "carbonyl",
        "conjugated alkene"
      ]
    },
    "optimal_conditions": {
      "catalyst_loading": "10-20 mol%",
      "solvent": "EtOH",
      "temperature": "0°C → RT (addition), then reflux (dehydration)",
      "atmosphere": "N2 or Ar (especially for LDA, Grignard-type bases)",
      "workup": "Quench with sat. NH4Cl or dilute acid; extract with EtOAc",
      "purification": "Column chromatography or distillation (if volatile)"
    },
    "scope_and_limitations": {
      "scope": [
        "Works well with aromatic aldehydes → cinnamaldehyde derivatives"
      ],
      "limitations": [
        "Aliphatic aldehydes tend toward multiple condensations (product retains α-H)"
      ]
    },
    "yield_expectation": "30-40% (estimated)",
    "predicted_side_reactions": [
      {
        "reaction": "Multiple condensation",
        "description": "Product still has α-H → can undergo further aldol reactions",
        "prevention": "Control stoichiometry, use kinetic enolate, low T"
      },
      {
        "reaction": "Dehydration to conjugated system",
        "description": "β-hydroxy carbonyl eliminates water under basic/thermal conditions",
        "prevention": "Desired in aldol condensation; control with T and time"
      },
      {
        "reaction": "Racemization",
        "description": "If chiral center adjacent to carbonyl, enolization causes racemization",
        "prevention": "Use mild conditions or alternative disconnection"
      }
    ],
    "recommendations": [
      "Self-aldol: control stoichiometry to minimize multiple additions",
      "Consider using aromatic aldehyde for cleaner crossed variant",
      "Upgrade to LDA/NaOEt for improved selectivity if needed"
    ],
    "stereochemical_notes": [
      "New stereocenters: up to 2 (α and β positions of new bond)",
      "Zimmerman-Traxler model: 6-membered cyclic TS determines syn/anti diastereoselectivity",
      "E-enolate → anti aldol; Z-enolate → syn aldol (Zimmerman-Traxler)"
    ],
    "summary": "Self-aldol of CC=O. Product: α,β-unsaturated carbonyl from CC=O. Expected yield: 30-40% (estimated)."
  }
}
```

## Summary

- **reaction_type**: Self-aldol of CC=O
- **substrate1_analysis**: {'alpha_h_pKa': 17, 'enolization': 'easy', 'electrophile_reactivity': 'high', 'self_aldol_tendency': 'high', 'examples':...
- **substrate2_analysis**: None
- **catalyst_analysis**: {'selected': 'generic base', 'type': 'base', 'notes': 'Base-type catalyst specified: base. See NaOH/KOH or NaOEt for typ...
- **product_prediction**: {'aldol_adduct': 'β-hydroxy carbonyl dimer from CC=O', 'condensation_product': 'α,β-unsaturated carbonyl from CC=O', 'fu...
