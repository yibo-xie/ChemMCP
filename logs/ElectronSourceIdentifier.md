# ElectronSourceIdentifier (Tool #152)

## Test Input

```python
ElectronSourceIdentifier(**{"reaction_input": "CH3CHO + NaBH4 → CH3CH2OH", "analysis_mode": "detailed"})
```

## Output (1540 chars)

```json
{
  "result": {
    "reaction_input": "CH3CHO + NaBH4 → CH3CH2OH",
    "primary_electron_source": "Na",
    "source_species": [
      "Na"
    ],
    "source_reasoning": "Sodium: strong reducing agent, Na → Na⁺ + e⁻",
    "matched_reducing_agents": [
      {
        "agent": "Na",
        "oxidation_change": "0 → +1",
        "description": "Sodium: strong reducing agent, Na → Na⁺ + e⁻"
      },
      {
        "agent": "H2",
        "oxidation_change": "0 → +1",
        "description": "Hydrogen gas: H₂ → 2H⁺ + 2e⁻"
      },
      {
        "agent": "C",
        "oxidation_change": "0 → +II/+IV",
        "description": "Carbon: C → CO/CO₂, reducing in metallurgy"
      },
      {
        "agent": "NaBH4",
        "oxidation_change": "H(-I) → +I",
        "description": "Sodium borohydride: hydride donor, H⁻ → H⁺"
      }
    ],
    "organic_reduction_patterns": [],
    "oxidation_state_changes": [],
    "electrons_transferred": null,
    "electrons_transferred_detail": {
      "note": "Insufficient data.",
      "total": null
    },
    "electron_flow_description": "Primary electron source: **Na**. Sodium: strong reducing agent, Na → Na⁺ + e⁻",
    "confidence": "high",
    "analysis_mode": "detailed",
    "additional_notes": [
      "ℹ️ Multiple potential reductants detected (4). Strongest identified as primary source."
    ],
    "half_reactions": [
      {
        "type": "oxidation (electron source)",
        "equation": "Na → Na⁺ + e⁻",
        "electrons": 1,
        "note": "alkali metal"
      }
    ]
  }
}
```

## Summary

- **reaction_input**: CH3CHO + NaBH4 → CH3CH2OH
- **primary_electron_source**: Na
- **source_species**: ['Na']
- **source_reasoning**: Sodium: strong reducing agent, Na → Na⁺ + e⁻
- **matched_reducing_agents**: [{'agent': 'Na', 'oxidation_change': '0 → +1', 'description': 'Sodium: strong reducing agent, Na → Na⁺ + e⁻'}, {'agent':...
