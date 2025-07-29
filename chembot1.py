import requests
import pandas as pd
import streamlit as st

# 🧬 Streamlit config
st.set_page_config(page_title="Compound Properties Resolver", layout="wide")
st.title("🔍 PubChem Compound Properties Resolver")

# 🌟 Primary lookup by compound name
def fetch_properties_by_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/MolecularWeight,IsomericSMILES,CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        properties = data.get("PropertyTable", {}).get("Properties", [])
        if properties:
            props = properties[0]
            return {
                "Name": name,
                "MolecularWeight": props.get("MolecularWeight", "Not Found"),
                "IsomericSMILES": props.get("IsomericSMILES", "Not Found"),
                "CanonicalSMILES": props.get("CanonicalSMILES", "Not Found")
            }
    return None

# 🔁 Fallback via CID
def fetch_properties_by_cid(name):
    cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    cid_response = requests.get(cid_url)
    if cid_response.status_code == 200:
        cid_data = cid_response.json()
        cids = cid_data.get("IdentifierList", {}).get("CID", [])
        if cids:
            cid = cids[0]
            prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,IsomericSMILES,CanonicalSMILES/JSON"
            prop_response = requests.get(prop_url)
            if prop_response.status_code == 200:
                prop_data = prop_response.json()
                properties = prop_data.get("PropertyTable", {}).get("Properties", [])
                if properties:
                    props = properties[0]
                    return {
                        "Name": name,
                        "MolecularWeight": props.get("MolecularWeight", "Not Found"),
                        "IsomericSMILES": props.get("IsomericSMILES", props.get("CanonicalSMILES", "Not Found")),
                        "CanonicalSMILES": props.get("CanonicalSMILES", "Not Found")
                    }
    # 🛑 Final fallback
    return {
        "Name": name,
        "MolecularWeight": "Not Found",
        "IsomericSMILES": "Not Found",
        "CanonicalSMILES": "Not Found"
    }

# 🧠 Input field
user_input = st.text_input("Enter compound names (comma-separated):")

# 🚀 Main logic
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

    # 📊 Show results
    df = pd.DataFrame(compound_data)
    st.subheader("🧬 Compound Lookup Results")
    st.dataframe(df)

    # 📥 CSV export
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 Download CSV",
        data=csv,
        file_name="compound_properties.csv",
        mime="text/csv"
    )
