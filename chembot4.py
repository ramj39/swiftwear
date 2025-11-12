import streamlit as st
import requests
import concurrent.futures
from rdkit import Chem
from requests.utils import quote
st.markdown(
    """
    <style>
    body, .stApp {
        background: linear-gradient(45deg, #ff9a9e 0%, #fad0c4 99%,#fad0c4 100%);
        min-height: 100vh;
        background-attachment: fixed;
    }
    </style>
    """,
    unsafe_allow_html=True
)
# Constants
PUBCHEM_FIELDS = [
    "MolecularFormula", "MolecularWeight", "XLogP", "InChI", "IUPACName", "CanonicalSMILES", "IsomericSMILES"
]
OUTPUT_FIELDS = [
    'source', 'formula', 'mw', 'xlogp', 'inchi', 'iupac'
]

# Utility
def validate_smiles(smiles):
    if not smiles:
        return False
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False

# Resolver Functions
def resolve_with_pubchem(name):
    try:
        encoded = quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/{','.join(PUBCHEM_FIELDS)}/JSON"
        r = requests.get(url, headers={'User-Agent': 'ChemicalResolver/1.0'}, timeout=5)
        r.raise_for_status()
        data = r.json()
        props = data.get('PropertyTable', {}).get('Properties', [])
        if not props:
            return None
        p = props[0]
        return {
            'source': 'PubChem',
            'CID': p.get('cid'),
            'smiles': p.get('IsomericSMILES'),
            'canonical_smiles': p.get('CanonicalSMILES'),
            'formula': p.get('MolecularFormula'),
            'mw': p.get('MolecularWeight'),
            'xlogp': p.get('XLogP'),
            'inchi': p.get('InChI'),
            'iupac': p.get('IUPACName'),
        }
    except Exception:
        return None

def resolve_with_cir(name):
    try:
        encoded = quote(name)
        url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/smiles"
        r = requests.get(url, headers={'User-Agent': 'ChemicalResolver/1.0'}, timeout=5)
        r.raise_for_status()
        smiles = r.text.strip()
        if not validate_smiles(smiles):
            return None
        return {
            'source': 'CIR',
            'smiles': smiles,
            'canonical_smiles': None,
            'formula': None,
            'mw': None,
            'xlogp': None,
            'inchi': None,
            'iupac': None
        }
    except Exception:
        return None

def resolve_with_opsin(name):
    try:
        encoded = quote(name)
        url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded}.json"
        r = requests.get(url, headers={'User-Agent': 'ChemicalResolver/1.0'}, timeout=5)
        r.raise_for_status()
        data = r.json()
        smiles = data.get("smiles")
        if not validate_smiles(smiles):
            return None
        return {
            'source': 'OPSIN',
            'smiles': smiles,
            'canonical_smiles': None,
            'formula': data.get("formula"),
            'mw': data.get("molWeight"),
            'xlogp': None,
            'inchi': None,
            'iupac': None
        }
    except Exception:
        return None

def resolve_name_to_all_props(name):
    result = resolve_with_pubchem(name)
    if result:
        return result
    with concurrent.futures.ThreadPoolExecutor() as ex:
        futures = [ex.submit(resolve_with_cir, name), ex.submit(resolve_with_opsin, name)]
        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            if res:
                return res
    if " " in name:
        result = resolve_with_pubchem(name.replace(" ", ""))
        if result:
            return result
    return {field: '' for field in OUTPUT_FIELDS} | {'source': 'None'}

def result_for_table(result):
    return {
        'Molecular Formula': result.get('formula', '-') or "-",
        'Molecular Weight': result.get('mw', '-') or "-",
        'XLogP': result.get('xlogp', '-') or "-",
        'InChI': result.get('inchi', '-') or "-",
        'IUPAC Name': result.get('iupac', '-') or "-",
        'Source': result.get('source', '-') or "-"
    }

# Streamlit UI
st.title("Chemical Name to Properties Resolver")

with st.form("resolver_form"):
    name = st.text_input("Enter chemical name", "")
    submitted = st.form_submit_button("Resolve")
    if submitted and name.strip():
        st.info(f"Resolving properties for: `{name}`")
        result = resolve_name_to_all_props(name.strip())
        st.write(f"**Source used:** `{result.get('source', '')}`")
        st.table(result_for_table(result))

        if result.get('smiles'):
            st.success("SMILES found!")
            st.write("SMILES:", result['smiles'])
        else:
            st.warning("No SMILES found for this name.")

# Optional: show PubChem response for a fixed CID
if st.checkbox("Show PubChem response for CID 2346"):
    cid = "2346"
    api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,InChI/JSON"
    r = requests.get(api_url)
    pubchem_response = r.json()
    st.json(pubchem_response)
    smiles = pubchem_response["PropertyTable"]["Properties"][0].get("CanonicalSMILES", "Not found")
    st.write("Canonical SMILES:", smiles)
