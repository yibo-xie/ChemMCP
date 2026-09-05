"""Generate Cherry Studio JSON config for ChemMCP tools #21-30."""
import importlib
import json

tool_names = [
    ('SymmetryPointGroup', 'symmetry_point_group'),
    ('CoordinationGeometry', 'coordination_geometry'),
    ('BalanceEquation', 'balance_equation'),
    ('BalanceRedox', 'balance_redox'),
    ('IdentifyReactionType', 'identify_reaction_type'),
    ('PredictProducts', 'predict_products'),
    ('GetStandardPotential', 'get_standard_potential'),
    ('CalculateCellPotential', 'calculate_cell_potential'),
    ('NernstEquation', 'nernst_equation'),
    ('IdentifyOxidizingAgent', 'identify_oxidizing_agent'),
]

mcp_tools = []
for class_name, module_name in tool_names:
    mod = importlib.import_module(f'chemmcp.tools.{module_name}')
    cls = getattr(mod, class_name)
    
    params = []
    for item in cls.code_input_sig:
        pname = item[0]
        ptype = item[1]
        pdesc = item[-1] if len(item) > 2 else ''
        params.append({
            'name': pname,
            'type': ptype,
            'description': pdesc,
        })
    
    mcp_tools.append({
        'name': cls.func_name,
        'description': cls.description,
        'parameters': {
            'type': 'object',
            'properties': {p['name']: {'type': p['type'].lower(), 'description': p['description']} for p in params},
            'required': [p['name'] for p in params],
        }
    })

config = {
    "mcpServers": {
        "ChemMCP": {
            "type": "stdio",
            "command": "python3",
            "args": ["-m", "chemmcp"],
            "env": {"PYTHONPATH": "/home/wave/ChemMCP/src"},
        }
    },
    "tools": mcp_tools,
}

output_path = "/home/wave/ChemMCP/cherry_studio_config_21_30.json"
with open(output_path, 'w') as f:
    json.dump(config, f, indent=2, ensure_ascii=False)

print(f"✅ Generated {output_path} with {len(mcp_tools)} tools")
for t in mcp_tools:
    print(f"  - {t['name']}: {t['description'][:70]}...")
