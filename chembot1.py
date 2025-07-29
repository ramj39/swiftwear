import requests
import rdkit
from rdkit import Chem

# ✅ Validate SMILES using RDKit
def validate_smiles(smiles):
    return Chem.MolFromSmiles(smiles) is not None if smiles else False

# 🔹 PubChem Resolver
def resolve_with_pubchem(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES,MolecularWeight/JSON"
        response = requests.get(url, timeout=5)
        data = response.json()
        props = data['PropertyTable']['Properties'][0]
        return {
            'source': 'PubChem',
            'name': name,
            'smiles': props.get('IsomericSMILES'),
            'mw': props.get('MolecularWeight')
        }
    except Exception as e:
        return None

# 🔹 NIH CIR Resolver
def resolve_with_cir(name):
    try:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
        response = requests.get(url, timeout=5)
        smiles = response.text.strip()
        return {
            'source': 'CIR',
            'name': name,
            'smiles': smiles,
            'mw': None
        } if validate_smiles(smiles) else None
    except Exception:
        return None

# 🔹 OPSIN Resolver
def resolve_with_opsin(name):
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        response = requests.get(url, timeout=5)
        data = response.json()
        smiles = data.get('smiles')
        return {
            'source': 'OPSIN',
            'name': name,
            'smiles': smiles,
            'mw': data.get('molWeight')
        } if validate_smiles(smiles) else None
    except Exception:
        return None

# 🧠 Unified Resolver Chain
def resolve_name_to_smiles(name):
    for resolver in [resolve_with_pubchem, resolve_with_cir, resolve_with_opsin]:
        result = resolver(name)
        if result and validate_smiles(result['smiles']):
            return result
    return {
        'source': 'None',
        'name': name,
        'smiles': None,
        'mw': None
    }
