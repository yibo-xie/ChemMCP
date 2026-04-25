# GrignardReaction (Tool #157)

## Test Input

```python
GrignardReaction(**{"grignard_reagent": "CH3MgBr", "electrophile_smiles": "benzaldehyde", "solvent": "dry Et2O"})
```

## Output (5208 chars)

```json
{
  "result": {
    "reaction": "CH3MgBr + benzaldehyde → CH3CH(elec.get('name','R'))OH (secondary alcohol)",
    "reagent_analysis": {
      "formula": "CH3MgBr",
      "type": "primary alkyl (methyl)",
      "nucleophilicity": "very high (carbanion character)",
      "basicity": "strong base"
    },
    "electrophile_analysis": {
      "type": "aldehyde",
      "name": "benzaldehyde",
      "electrophile": "formaldehyde (HCHO)",
      "product": "primary alcohol (+1 carbon)",
      "notes": "R-MgX + HCHO → RCH2OH"
    },
    "mechanism": [
      "1. Nucleophilic attack: CH3MgBr transfers R⁻ (with carbanion character) to electrophilic center of benzaldehyde",
      "2. Tetrahedral intermediate formed (alkoxide if carbonyl electrophile)",
      "3. Acidic workup (sat. NH4Cl or dilute HCl/H2SO4): protonation of intermediate",
      "4. Product isolation: CH3CH(elec.get('name','R'))OH (secondary alcohol)"
    ],
    "product_prediction": {
      "name": "CH3CH(elec.get('name','R'))OH (secondary alcohol)",
      "class": "secondary alcohol",
      "yield_estimate": "75-92%"
    },
    "scope_table": [
      {
        "electrophile": "formaldehyde (HCHO)",
        "product": "primary alcohol (+1 carbon)",
        "notes": "R-MgX + HCHO → RCH2OH"
      },
      {
        "electrophile": "aldehyde (R'CHO)",
        "product": "secondary alcohol",
        "notes": "R-MgX + R'CHO → RCH(R')OH"
      },
      {
        "electrophile": "ketone (R'2CO)",
        "product": "tertiary alcohol",
        "notes": "R-MgX + R'2CO → CR'3(R)OH"
      },
      {
        "electrophile": "ester (R'COOR\")",
        "product": "tertiary alcohol (after 2 equiv)",
        "notes": "First: ketone intermediate; second: tertiary alcohol. 2 eq R-MgX needed."
      },
      {
        "electrophile": "acid chloride (R'COCl)",
        "product": "tertiary alcohol (after 2 equiv)",
        "notes": "Reactive — often over-addition. Can stop at ketone at low T with CuI catalysis."
      },
      {
        "electrophile": "epoxide",
        "product": "alcohol (ring-opened at less substituted C)",
        "notes": "Regioselective SN2-like opening; CuI catalysis can invert selectivity"
      },
      {
        "electrophile": "CO2 (dry ice)",
        "product": "carboxylic acid (+1 carbon)",
        "notes": "Classic carboxylation: R-MgX + CO2 → RCOOH after acidic workup"
      },
      {
        "electrophile": "nitrile (R'CN)",
        "product": "ketone after hydrolysis",
        "notes": "R-MgX + R'CN → imine intermediate → hydrolyze → ketone"
      },
      {
        "electrophile": "DMF, DMSO, etc.",
        "product": "aldehyde after hydrolysis",
        "notes": "Formyl equivalent: R-MgX + HC(OR)3 → aldehyde; DMF → aldehyde"
      },
      {
        "electrophile": "alkyl halide (R'X)",
        "product": "coupling (C-C bond formation)",
        "notes": "Possible but competes with Wurtz-type coupling; transition metal catalysis preferred nowadays"
      },
      {
        "electrophile": "oxygen (O2)",
        "product": "alcohol (after hydrolysis)",
        "notes": "Side reaction to avoid — keep under inert atmosphere"
      }
    ],
    "applicable_limitations": [
      {
        "issue": "Enolizable protons α to carbonyl",
        "problem": "Can be deprotonated instead of addition",
        "solution": "Use lower T and controlled addition"
      }
    ],
    "safety_notes": [
      [
        "Pyrophoric",
        "Grignard reagents ignite spontaneously in air",
        "Keep under N2/Ar at all times; use septum techniques"
      ],
      [
        "Exothermic",
        "Formation and carbonyl addition are highly exothermic",
        "Cool in ice bath; add slowly with stirring"
      ],
      [
        "Ether solvent fire hazard",
        "Et2O is highly flammable, forms explosive peroxides",
        "Test for peroxides; use fresh or properly stored ether"
      ],
      [
        "Mg dust",
        "Fine Mg is flammable",
        "Handle carefully; avoid ignition sources"
      ],
      [
        "Quenching",
        "Always quench excess Grignard carefully with sat. NH4Cl or dilute acid",
        "Never add water directly to concentrated Grignard solution"
      ],
      [
        "Pressure buildup",
        "Reaction produces gas during quenching",
        "Quench slowly with cooling and venting"
      ]
    ],
    "workup_procedure": "After reaction completion (monitored by TLC or GC), cool to 0°C.\nSlowly pour onto saturated aqueous NH4Cl solution (or dilute HCl) with vigorous stirring.\nExtract with Et2O or EtOAc (×3).\nWash combined organics with brine, dry over Na2SO4 or MgSO4.\nFilter, concentrate, purify by column chromatography or distillation.",
    "optimal_conditions": {
      "solvent": "dry Et2O",
      "temperature": "0°C → RT (addition exothermic!)",
      "atmosphere": "N2 or Ar (rigorous exclusion of air/moisture)",
      "addition_order": "Add Grignard reagent dropwise to electrophile solution (minimizes enolization/wurtz side reactions)",
      "monitoring": "TLC or GC for consumption of electrophile"
    },
    "summary": "Grignard addition: CH3CH(elec.get('name','R'))OH (secondary alcohol). Yield: 75-92%."
  }
}
```

## Summary

- **reaction**: CH3MgBr + benzaldehyde → CH3CH(elec.get('name','R'))OH (secondary alcohol)
- **reagent_analysis**: {'formula': 'CH3MgBr', 'type': 'primary alkyl (methyl)', 'nucleophilicity': 'very high (carbanion character)', 'basicity...
- **electrophile_analysis**: {'type': 'aldehyde', 'name': 'benzaldehyde', 'electrophile': 'formaldehyde (HCHO)', 'product': 'primary alcohol (+1 carb...
- **mechanism**: ['1. Nucleophilic attack: CH3MgBr transfers R⁻ (with carbanion character) to electrophilic center of benzaldehyde', '2. ...
- **product_prediction**: {'name': "CH3CH(elec.get('name','R'))OH (secondary alcohol)", 'class': 'secondary alcohol', 'yield_estimate': '75-92%'}
