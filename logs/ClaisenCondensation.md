# ClaisenCondensation (Tool #155)

## Test Input

```python
ClaisenCondensation(**{"ester1_smiles": "CC(=O)OCC", "ester2_smiles": "", "base": "NaOEt", "solvent": "EtOH"})
```

## Output (4403 chars)

```json
{
  "result": {
    "reaction_type": "Self-Claisen condensation of CC(=O)OCC",
    "ester1_analysis": {
      "classified_as": "generic_aliphatic_ester",
      "input": "CC(=O)OCC"
    },
    "ester2_analysis": null,
    "base_analysis": {
      "selected": "NaOEt/EtOH",
      "type": "classic",
      "strength": "pKa EtOH ≈ 16",
      "conditions": "reflux in absolute EtOH",
      "notes": "Classic Claisen conditions; base must match ester alkoxy group to prevent mixed esters via transesterification",
      "catalytic_vs_stoich": "1 equiv minimum (consumed to form β-ketoester salt)"
    },
    "mechanism_summary": [
      {
        "step": "1",
        "name": "Enolate formation",
        "description": "Base (e.g., OEt⁻) abstracts α-proton from ester → resonance-stabilized enolate"
      },
      {
        "step": "2",
        "name": "Nucleophilic acyl substitution",
        "description": "Enolate attacks carbonyl carbon of another ester molecule → tetrahedral intermediate"
      },
      {
        "step": "3",
        "name": "Elimination of alkoxide",
        "description": "Tetrahedral intermediate collapses, expelling -OR (e.g., -OEt) → β-ketoester"
      },
      {
        "step": "4",
        "name": "Deprotonation of product",
        "description": "The β-ketoester is more acidic (pKa ~11) than starting ester → deprotonated by remaining base → stable enolate salt"
      },
      {
        "step": "5",
        "name": "Acidic workup",
        "description": "Add dilute acid → protonate enolate → neutral β-ketoester product"
      }
    ],
    "product_prediction": {
      "name": "β-ketoester from CC(=O)OCC self-condensation",
      "class": "β-ketoester (1,3-dicarbonyl)"
    },
    "optimal_conditions": {
      "base_amount": "1 equivalent (minimum; consumed in reaction)",
      "solvent": "EtOH (absolute/anhydrous)",
      "temperature": "reflux temperature of solvent",
      "time": "2-12 hours",
      "atmosphere": "N2/Ar (anhydrous)",
      "workup": "Cool, acidify with dilute HCl (pH ~4), extract with organic solvent",
      "critical_note": "MUST be anhydrous — water destroys base and stops reaction"
    },
    "scope": [
      "Aliphatic esters with α-H work well",
      "Ethyl acetate is the classic substrate → ethyl acetoacetate",
      "Cyclic diesters → Dieckmann cyclization (excellent for 5/6 rings)",
      "Crossed Claisen works cleanly when one ester lacks α-H",
      "β-Ketoesters are versatile synthetic intermediates (alkylation, Knoevenagel, etc.)"
    ],
    "limitations": [
      "Esters without α-H cannot undergo self-Claisen (only as electrophiles)",
      "Both esters having α-H in crossed Claisen → messy mixture",
      "Sterically hindered α-positions (e.g., t-butyl ester α-position) react slowly",
      "Transesterification can compete if base/ester alkoxy groups do not match",
      "Requires anhydrous conditions (base consumed by any water present)",
      "Product must be acidic enough to be fully deprotonated (drives equilibrium)"
    ],
    "variants": {
      "Crossed_Claisen": {
        "key_requirement": "One ester MUST lack α-H (e.g., ethyl benzoate, ethyl formate, aromatic ester)",
        "example": "EtOAc + EtOBz → benzoylacetate",
        "risk": "If both have α-H → mixture of 4 products"
      },
      "Dieckmann_condensation": {
        "ring_sizes": "5- and 6-membered rings favored (entropic advantage); 7-membered possible but slower",
        "example": "Diethyl heptanedioate (diethyl pimelate) → ethyl 2-oxocyclopentanecarboxylate (5-membered ring)",
        "ring_preference": "5 > 6 >> 7 membered rings"
      },
      "Claisen-Schmidt": {
        "note": "This is actually aldol (ketone/aldehyde + aldehyde), not a true Claisen. Named similarly but different mechanism."
      },
      "Acetoacetic_ester_synthesis": {
        "sequence": "(1) Claisen self-condensation of EtOAc → acetoacetic ester (2) Alkylation (3) Hydrolysis + decarboxylation → methyl ketone"
      },
      "Malonic_ester_synthesis": {
        "relation": "Parallel to acetoacetic ester synthesis"
      },
      "Thorpie-Ziegler_modification": {
        "notes": "Broader scope beyond simple esters"
      }
    },
    "yield_expectation": "42-58%",
    "summary": "Self-Claisen condensation of CC(=O)OCC. Product: β-ketoester from CC(=O)OCC self-condensation. Driving force: β-ketoester stabilization."
  }
}
```

## Summary

- **reaction_type**: Self-Claisen condensation of CC(=O)OCC
- **ester1_analysis**: {'classified_as': 'generic_aliphatic_ester', 'input': 'CC(=O)OCC'}
- **ester2_analysis**: None
- **base_analysis**: {'selected': 'NaOEt/EtOH', 'type': 'classic', 'strength': 'pKa EtOH ≈ 16', 'conditions': 'reflux in absolute EtOH', 'not...
- **mechanism_summary**: [{'step': '1', 'name': 'Enolate formation', 'description': 'Base (e.g., OEt⁻) abstracts α-proton from ester → resonance-...
