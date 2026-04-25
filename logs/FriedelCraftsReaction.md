# FriedelCraftsReaction (Tool #159)

## Test Input

```python
FriedelCraftsReaction(**{"arene_smiles": "benzene", "electrophile_type": "alkyl", "electrophile_spec": "CH3CH2Cl", "catalyst": "AlCl3"})
```

## Output (2716 chars)

```json
{
  "result": {
    "reaction": "FC Alkylation: benzene + CH3CH2Cl → alkylated arene",
    "arene_analysis": {
      "patterns": [
        "benzene",
        "c1ccccc1"
      ],
      "substituent": null,
      "activation": "neutral",
      "directing": "none",
      "reactivity": "baseline"
    },
    "electrophile_analysis": {
      "type": "alkyl",
      "spec": "CH3CH2Cl",
      "active_species": "carbocation (R⁺)",
      "rearrangement_possible": true,
      "stability": "depends on structure (1° rearranges, 2°/3° stable)"
    },
    "catalyst_analysis": {
      "selected": "AlCl3",
      "strength": "very strong",
      "stoich": "catalytic (alkylation) or ≥1 eq (acylation)",
      "scope": "universal for FC",
      "handling": "moisture sensitive, exothermic on addition",
      "notes": "Gold standard FC catalyst"
    },
    "mechanism": [],
    "orientation": {
      "pattern": "monosubstitution (all positions equivalent for benzene)",
      "isomer_distribution": "single product"
    },
    "product_prediction": {
      "name": "alkylated arene",
      "class": "alkyl arene",
      "yield": "50-80%"
    },
    "applicable_limitations": [
      {
        "issue": "Polyalkylation",
        "problem": "Product more activated than starting material",
        "solution": "Use large excess of arene (10+ eq)"
      },
      {
        "issue": "Rearrangement",
        "problem": "Carbocation may rearrange (especially 1° R groups)",
        "solution": "Use acylation + reduction (Clemmensen/WK) instead"
      }
    ],
    "fc_alkylation_vs_acylation": {
      "rearrangement": "Alkylation: YES (carbocation); Acylation: NO (acylium stable)",
      "poly_substitution": "Alkylation: YES (product activated); Acylation: NO (product deactivated)",
      "catalyst_amount": "Alkylation: catalytic; Acylation: stoichiometric (complexes product)",
      "product_type": "Alkyl arene vs aryl ketone",
      "indirect_route": "Acylation → Clemmensen/WK = rearrangement-free alkylation"
    },
    "optimal_conditions": {
      "catalyst_loading": "catalytic (0.1-0.2 eq)",
      "solvent": "CS2, CH2Cl2, CH3NO2, or neat arene (excess as solvent)",
      "temperature": "RT to gentle reflux (0-80°C typical)",
      "atmosphere": "anhydrous (exclude moisture — deactivates Lewis acid)",
      "addition_order": "Add catalyst to arene solution, then add electrophile slowly with cooling",
      "workup": "Pour onto ice/water carefully (EXOTHERMIC!); extract with organic solvent; wash, dry, purify",
      "safety": "HIGHLY EXOTHERMIC reaction initiation — add slowly with good cooling and stirring!"
    },
    "summary": "FC alkyllation: alkylated arene. Yield: 50-80%. Key: Polyalkylation"
  }
}
```

## Summary

- **reaction**: FC Alkylation: benzene + CH3CH2Cl → alkylated arene
- **arene_analysis**: {'patterns': ['benzene', 'c1ccccc1'], 'substituent': None, 'activation': 'neutral', 'directing': 'none', 'reactivity': '...
- **electrophile_analysis**: {'type': 'alkyl', 'spec': 'CH3CH2Cl', 'active_species': 'carbocation (R⁺)', 'rearrangement_possible': True, 'stability':...
- **catalyst_analysis**: {'selected': 'AlCl3', 'strength': 'very strong', 'stoich': 'catalytic (alkylation) or ≥1 eq (acylation)', 'scope': 'univ...
- **mechanism**: []
