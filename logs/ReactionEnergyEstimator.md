# ReactionEnergyEstimator (Tool #153)

## Test Input

```python
ReactionEnergyEstimator(**{"reactants_smiles": "C=C + H2", "products_smiles": "CC", "temperature_k": 298.15})
```

## Output (1050 chars)

```json
{
  "result": {
    "reaction": "C=C + H2 → CC",
    "delta_h": 0,
    "delta_h_unit": "kJ/mol",
    "delta_s_estimate": 130,
    "delta_s_unit": "J/(mol·K)",
    "delta_g_estimate": -38.8,
    "delta_g_unit": "kJ/mol",
    "equilibrium_constant": "6.18e+06",
    "log_k": 6.79,
    "spontaneous": true,
    "spontaneity_description": "Spontaneous (negative ΔG)",
    "feasibility_rating": "very good — high conversion expected (>99%)",
    "method_used": "no_data_available",
    "confidence": "very_low",
    "entropy_method": "gas_molecule_count_heuristic + reaction_type_pattern",
    "temperature_k": 298.15,
    "thermodynamic_notes": [
      "Near-thermoneutral: entropy dominates the spontaneity.",
      "Large positive ΔS: favored by increasing temperature."
    ],
    "temperature_effect": "Endothermic-like (ΔS>0): increasing T makes reaction MORE favorable (larger K).",
    "practical_considerations": [
      "Consider kinetic factors: even thermodynamically favorable reactions may require catalyst, heat, or activation."
    ]
  }
}
```

## Summary

- **reaction**: C=C + H2 → CC
- **delta_h**: 0
- **delta_h_unit**: kJ/mol
- **delta_s_estimate**: 130
- **delta_s_unit**: J/(mol·K)
