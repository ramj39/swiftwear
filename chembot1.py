import sys
import csv
import time
import requests
from rdkit import Chem
from requests.utils import quote

# ✅ Validate SMILES using RDKit
def validate_smiles(smiles):
    if not smiles:
        return False
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False

# 🔹 PubChem Resolver
def resolve_with_pubchem(name):
    try:
        encoded_name = quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/IsomericSMILES,MolecularWeight/JSON"
        headers = {'User-Agent': 'ChemicalResolver/1.0'}
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        # Handle case where no properties are found
        if not data.get('PropertyTable', {}).get('Properties'):
            return None
            
        props = data['PropertyTable']['Properties'][0]
        return {
            'source': 'PubChem',
            'name': name,
            'smiles': props.get('IsomericSMILES'),
            'mw': props.get('MolecularWeight')
        }
    except Exception as e:
        print(f"PubChem error for '{name}': {str(e)}", file=sys.stderr)
        return None

# 🔹 NIH CIR Resolver
def resolve_with_cir(name):
    try:
        encoded_name = quote(name)
        url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded_name}/smiles"
        headers = {'User-Agent': 'ChemicalResolver/1.0'}
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        smiles = response.text.strip()
        
        if not validate_smiles(smiles):
            return None
            
        return {
            'source': 'CIR',
            'name': name,
            'smiles': smiles,
            'mw': None
        }
    except Exception as e:
        print(f"CIR error for '{name}': {str(e)}", file=sys.stderr)
        return None

# 🔹 OPSIN Resolver
def resolve_with_opsin(name):
    try:
        encoded_name = quote(name)
        url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.json"
        headers = {'User-Agent': 'ChemicalResolver/1.0'}
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()
        smiles = data.get('smiles')
        
        if not validate_smiles(smiles):
            return None
            
        return {
            'source': 'OPSIN',
            'name': name,
            'smiles': smiles,
            'mw': data.get('molWeight')
        }
    except Exception as e:
        print(f"OPSIN error for '{name}': {str(e)}", file=sys.stderr)
        return None

# 🧠 Unified Resolver
def resolve_name_to_smiles(name):
    # Try all resolvers in sequence
    for resolver in [resolve_with_pubchem, resolve_with_opsin, resolve_with_cir]:
        result = resolver(name)
        if result:
            return result
    
    # Fallback: Try name variations
    variations = [
        name,
        name.lower(),
        name.title(),
        name.upper()
    ]
    
    # Try variations with different resolvers
    for variant in variations:
        for resolver in [resolve_with_pubchem, resolve_with_opsin, resolve_with_cir]:
            result = resolver(variant)
            if result:
                result['name'] = name  # Keep original name
                return result
    
    return {
        'source': 'None',
        'name': name,
        'smiles': None,
        'mw': None
    }

def main():
    if len(sys.argv) != 3:
        print("Usage: python resolver.py input.txt output.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        writer.writerow(['name', 'smiles', 'source', 'mw'])
        
        for i, row in enumerate(reader):
            if not row:
                continue
            name = row[0].strip()
            if not name:
                continue
                
            # Skip header row
            if i == 0 and name.lower() in ['name', 'compound', 'chemical']:
                continue
                
            try:
                result = resolve_name_to_smiles(name)
                writer.writerow([
                    result['name'],
                    result['smiles'] or '',
                    result['source'],
                    result['mw'] or ''
                ])
            except Exception as e:
                print(f"Error processing '{name}': {str(e)}", file=sys.stderr)
                writer.writerow([name, '', 'Error', ''])
            
            outfile.flush()  # Write after each row
            time.sleep(0.3)  # Be polite to servers

if __name__ == '__main__':
    main()
