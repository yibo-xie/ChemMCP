# CalculateLatticeEnergy

> **Version:** 0.1.0 | **Function:** `calculate_lattice_energy`

## Description
Calculate or query lattice energy of ionic solids using Born-Landé equation and Kapustinskii approximation.

## Implementation Details
Calculates lattice energy via Born-Landé equation: U = -(N_A · M · z⁺ · z⁻ · e²)/(4πε₀ r₀)(1 - 1/n), where M is Madelung constant, z are ion charges, r₀ is interionic distance, n is Born exponent. Also provides Kapustinskii approximation as fallback and experimental values for comparison.

## Input Format
### Code Input
| Param | Type | Default | Desc |
|-------|------|--------|------|
| `('compound` | str | N/A | Ionic compound formula (e.g., NaCl, MgO, CaF2)'), |
| `('method` | str | born_landé | Calculation method: born_lande, kapustinskii, or both'), |

## Output Format
| Field | Type | Description |
|-------|------|-------------|
| `compound` | str | 'Compound formula' |
| `lattice_energy_kj_mol` | float | 'Calculated lattice energy (kJ/mol, negative = exothermic |
| `method` | str | 'Calculation method used' |
| `parameters_used` | dict | 'Parameters: Madelung constant, charges, interionic distance, Born exponent' |
| `experimental_value` | float | 'Experimental/literature value if available' |
| `description` | str | 'Explanation of calculation' |

## Examples
```python
examples = [
{'code_input': {'compound': 'NaCl', 'method': 'born-lande'}, 'text_input': {'query': 'NaCl'}, 'output': {'compound': 'NaCl', 'lattice_energy_kj_mol': -786, 'description': '...', 'method': '...', 'parameters_used': {...}, 'experimental_value': -787}},
        {'code_input': {'compound': 'MgO', 'method': 'born-lande'}, 'text_input': {'query': 'MgO'}, 'output': {'compound': 'MgO', 'lattice_energy_kj_mol': -3795, 'description': '...', 'method': '...', 'parameters_used': {...}, 'experimental_value': -3795}},
        {'code_input': {'compound': 'CaF2', 'method': 'kapustinskii'}, 'text_input': {'query': 'CaF2'}, 'output': {'compound': 'CaF2', 'lattice_energy_kj_mol': -2838, 'description': '...', 'method': '...', 'parameters_used': {...}, 'experimental_value': -2838}},
]
```

## Cherry Studio JSON Config
```json
{
  "name": "CalculateLatticeEnergy",
  "description": "Calculate or query lattice energy of ionic solids using Born-Landé equation and Kapustinskii approximation.",
  "function": "calculate_lattice_energy",
  "parameters": {
    "('compound": {
      "type": "str",
      "description": "Ionic compound formula (e.g., NaCl, MgO, CaF2)'),"
    },
    "('method": {
      "type": "str",
      "description": "Calculation method: born_lande, kapustinskii, or both'),"
    }
  }
}
```
