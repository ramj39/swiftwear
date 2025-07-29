import requests
import pandas as pd
import streamlit as st

# 🧪 Streamlit app setup
st.set_page_config(page_title="Compound Properties Resolver", layout="wide")
st.title("🔍 PubChem Compound Properties Resolver")
isomeric = props.get("IsomericSMILES")
canonical = props.get("CanonicalSMILES")

smiles = isomeric if isomeric else canonical if canonical else "Not Found"

# 🌟 Lookup by compound name (direct query)
def fetch_properties_by_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/MolecularWeight,IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return {
                "Name": name,
                "MolecularWeight": props[0].get("MolecularWeight"),
                "IsomericSMILES": props[0].get("IsomericSMILES", "Not Found")
            }
    return None

# 🔁 Fallback: Resolve to CID first, then fetch properties
def fetch_properties_by_cid(name):
    cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    cid_response = requests.get(cid_url)
    if cid_response.status_code == 200:
        cid_data = cid_response.json()
        cids = cid_data.get("IdentifierList", {}).get("CID", [])
        if cids:
            cid = cids[0]
            prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,IsomericSMILES/JSON"
            prop_response = requests.get(prop_url)
            if prop_response.status_code == 200:
                prop_data = prop_response.json()
                props = prop_data.get("PropertyTable", {}).get("Properties", [])
                if props:
                    return {
                        "Name": name,
                        "MolecularWeight": props[0].get("MolecularWeight"),
                        "IsomericSMILES": props[0].get("IsomericSMILES", "Not Found")
                    }
    # ❗ If everything fails
    return {
        "Name": name,
        "MolecularWeight": "Not Found",
        "IsomericSMILES": "Not Found"
    }

# 🧠 Main input logic
user_input = st.text_input("Enter compound names (comma-separated):")

if user_input:
    compound_names = [name.strip() for name in user_input.split(",")]
    compound_data = []

    for name in compound_names:
        result = fetch_properties_by_name(name)
        if result:
            compound_data.append(result)
        else:
            st.warning(f"Primary lookup failed for '{name}'. Trying fallback via CID...")
            fallback_result = fetch_properties_by_cid(name)
            compound_data.append(fallback_result)

    # 📊 Results
    df = pd.DataFrame(compound_data)
    st.subheader("🧬 Compound Lookup Results")
    st.dataframe(df)

    # 📥 CSV download option
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 Download CSV",
        data=csv,
        file_name="compound_properties.csv",
        mime="text/csv"
    )
