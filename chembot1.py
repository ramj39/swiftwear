import requests
import pandas as pd
import streamlit as st

st.title("PubChem Compound Properties Explorer")

# User input via Streamlit text box
user_input = st.text_input("Enter compound names (comma-separated):")

# Only proceed if user has entered something
if user_input:
    compound_names = [name.strip() for name in user_input.split(",")]
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/MolecularWeight,CanonicalSMILES/JSON"
    compound_data = []

    for name in compound_names:
        url = base_url.format(name)
        response = requests.get(url)
        if response.status_code == 200:
            try:
                data = response.json()
                if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                    properties = data["PropertyTable"]["Properties"][0]
                    compound_data.append({
                        "Name": name,
                        "MolecularWeight": properties.get("MolecularWeight"),
                        "CanonicalSMILES": properties.get("CanonicalSMILES")
                    })
                else:
                    st.warning(f"No valid properties found for {name}.")
            except KeyError:
                st.error(f"Key Error: Missing properties for {name}.")
        else:
            st.error(f"Error {response.status_code}: Unable to fetch data for {name}.")

    # Display result table
    if compound_data:
        df = pd.DataFrame(compound_data)
        st.dataframe(df)

        # Option to download CSV
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("Download CSV", csv, "compound_properties.csv", "text/csv")

