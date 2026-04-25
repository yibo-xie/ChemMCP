# ElectronSinkIdentifier (Tool #151)

## Test Input

```python
ElectronSinkIdentifier(**{"reaction_input": "2KMnO4 + 16HCl → 2KCl + 2MnCl2 + 5Cl2 + 8H2O", "analysis_mode": "detailed"})
```

## Output (2650 chars)

```json
{
  "result": {
    "reaction_input": "2KMnO4 + 16HCl → 2KCl + 2MnCl2 + 5Cl2 + 8H2O",
    "primary_electron_sink": "KMnO4",
    "sink_species": [
      "KMnO4"
    ],
    "sink_reasoning": "Permanganate: Mn(VII) reduced, color purple→colorless/brown",
    "matched_oxidizing_agents": [
      {
        "agent": "KMnO4",
        "oxidation_change": "Mn(+7) → +2/+4/+6",
        "description": "Permanganate: Mn(VII) reduced, color purple→colorless/brown"
      },
      {
        "agent": "Cl2",
        "oxidation_change": "0 → -1",
        "description": "Chlorine: Cl₂ → 2Cl⁻"
      },
      {
        "agent": "MnO4⁻",
        "oxidation_change": "varies",
        "description": "Permanganate: Strong oxidant, Mn(VII)→lower"
      }
    ],
    "organic_oxidation_patterns": [],
    "oxidation_state_changes": [
      {
        "element": "Mn",
        "from_oxidation_state": "+7",
        "to_oxidation_state": "+2",
        "change": -5,
        "is_reduction": true,
        "is_oxidation": false
      },
      {
        "element": "Cl",
        "from_oxidation_state": "-1",
        "to_oxidation_state": "0",
        "change": 1,
        "is_reduction": false,
        "is_oxidation": true
      },
      {
        "element": "H",
        "from_oxidation_state": "0",
        "to_oxidation_state": "+1",
        "change": 1,
        "is_reduction": false,
        "is_oxidation": true
      }
    ],
    "electrons_transferred": 5,
    "electrons_transferred_detail": {
      "electrons_gained_by_sink": 5,
      "electrons_lost_by_source": 2,
      "total": 5
    },
    "electron_flow_description": "Primary electron sink: **KMnO4**. Permanganate: Mn(VII) reduced, color purple→colorless/brown Electron source(s): Cl, H — these species are oxidized (lose electrons). Electron sink(s): Mn (+7 → +2) — these species are reduced (gain electrons).",
    "confidence": "high",
    "analysis_mode": "detailed",
    "additional_notes": [
      "ℹ️ Multiple potential oxidants detected (3). The strongest/most specific one is identified as the primary sink."
    ],
    "half_reactions": [
      {
        "type": "reduction (electron sink)",
        "equation": "MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O  (acidic)",
        "alternative": "MnO4⁻ + 2H2O + 3e⁻ → MnO2 + 4OH⁻  (basic)",
        "electrons": 5
      },
      {
        "type": "oxidation (electron source)",
        "equation": "Cl (oxidation state -1) → Cl (oxidation state 0) + 1e⁻",
        "electrons": 1
      },
      {
        "type": "oxidation (electron source)",
        "equation": "H (oxidation state 0) → H (oxidation state +1) + 1e⁻",
        "electrons": 1
      }
    ]
  }
}
```

## Summary

- **reaction_input**: 2KMnO4 + 16HCl → 2KCl + 2MnCl2 + 5Cl2 + 8H2O
- **primary_electron_sink**: KMnO4
- **sink_species**: ['KMnO4']
- **sink_reasoning**: Permanganate: Mn(VII) reduced, color purple→colorless/brown
- **matched_oxidizing_agents**: [{'agent': 'KMnO4', 'oxidation_change': 'Mn(+7) → +2/+4/+6', 'description': 'Permanganate: Mn(VII) reduced, color purple...
