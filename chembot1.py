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

# 🔹 PubChem Resolver with enhanced error handling
def resolve_with_pubchem(name):
    try:
        encoded_name = quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/IsomericSMILES,MolecularWeight/JSON"
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'}
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()  # Raise exception for HTTP errors
        data = response.json()
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

# 🔹 NIH CIR Resolver with multiple endpoints
def resolve_with_cir(name):
    endpoints = [
        "smiles",
        "stdinchikey",
        "iupac_name"
    ]
    
    for endpoint in endpoints:
        try:
            encoded_name = quote(name)
            url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded_name}/{endpoint}"
            headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'}
            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()
            
            content = response.text.strip()
            if endpoint == "smiles" and validate_smiles(content):
                return {
                    'source': 'CIR',
                    'name': name,
                    'smiles': content,
                    'mw': None
                }
            elif endpoint == "stdinchikey" and content:
                # Use InChIKey to get SMILES from PubChem
                inchikey = content
                pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/IsomericSMILES/JSON"
                pub_response = requests.get(pubchem_url, headers=headers, timeout=10)
                pub_response.raise_for_status()
                pub_data = pub_response.json()
                smiles = pub_data['PropertyTable']['Properties'][0]['IsomericSMILES']
                if validate_smiles(smiles):
                    return {
                        'source': 'CIR+PubChem',
                        'name': name,
                        'smiles': smiles,
                        'mw': None
                    }
        except Exception as e:
            print(f"CIR error for '{name}' ({endpoint}): {str(e)}", file=sys.stderr)
            continue
    
    return None

# 🔹 OPSIN Resolver with improved error handling
def resolve_with_opsin(name):
    try:
        encoded_name = quote(name)
        url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.json"
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'}
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()
        smiles = data.get('smiles')
        if validate_smiles(smiles):
            return {
                'source': 'OPSIN',
                'name': name,
                'smiles': smiles,
                'mw': data.get('molWeight')
            }
    except Exception as e:
        print(f"OPSIN error for '{name}': {str(e)}", file=sys.stderr)
    
    return None

# 🧠 Unified Resolver Chain with fallback strategies
def resolve_name_to_smiles(name):
    # Try all resolvers in sequence
    for resolver in [resolve_with_pubchem, resolve_with_opsin, resolve_with_cir]:
        result = resolver(name)
        if result and result.get('smiles'):
            return result
    
    # Fallback 1: Try name without special characters
    clean_name = ''.join(e for e in name if e.isalnum() or e in " -_")
    if clean_name != name:
        for resolver in [resolve_with_pubchem, resolve_with_opsin, resolve_with_cir]:
            result = resolver(clean_name)
            if result and result.get('smiles'):
                result['name'] = name  # Keep original name in results
                return result
    
    # Fallback 2: Try common name variations
    variations = [
        name.lower(),
        name.title(),
        name.upper()
    ]
    
    for variant in variations:
        if variant == name:
            continue
        for resolver in [resolve_with_pubchem, resolve_with_opsin, resolve_with_cir]:
            result = resolver(variant)
            if result and result.get('smiles'):
                result['name'] = name  # Keep original name in results
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
                
            # Skip header row if exists
            if i == 0 and name.lower() in ['name', 'compound', 'chemical']:
                continue
                
            result = resolve_name_to_smiles(name)
            writer.writerow([
                result['name'],
                result['smiles'] if result['smiles'] else '',
                result['source'],
                result['mw'] if result['mw'] is not None else ''
            ])
            outfile.flush()  # Write after each row
            time.sleep(0.2)   # Be polite to servers

if _name_ == '_main_':
    main()
