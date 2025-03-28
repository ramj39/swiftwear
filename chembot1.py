import requests
import pandas as pd

# Prompt the user to input compound names (comma-separated)
user_input = input("Enter compound names (comma-separated): ")
compound_names = [name.strip() for name in user_input.split(",")]

# Base URL for PubChem API to query by name
base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/MolecularWeight,CanonicalSMILES/JSON"

# Data storage for results
compound_data = []

# Query each compound name and retrieve properties
for name in compound_names:
    url = base_url.format(name)
    response = requests.get(url)
    print(f"Querying: {url}")  # Debug: Print the API URL being queried
    if response.status_code == 200:  # Only handle valid responses
        try:
            data = response.json()
            print(f"Response Data for {name}: {data}")  # Debug: Print response content
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                properties = data["PropertyTable"]["Properties"][0]
                compound_data.append({
                    "Name": name,
                    "MolecularWeight": properties.get("MolecularWeight"),
                    "CanonicalSMILES": properties.get("CanonicalSMILES")
                })
            else:
                print(f"Warning: No valid properties found for {name}.")
        except KeyError:
            print(f"Key Error: Missing properties for {name}.")
    else:  # Handle non-successful status codes here
        print(f"Error {response.status_code}: Unable to fetch data for {name}.")

# Convert the results to a DataFrame for tabulation
df = pd.DataFrame(compound_data)

# Display the table
print(df)

# Save the table to a CSV file
df.to_csv("compound_properties_by_name.csv", index=False)

# Uncomment to save as Excel if needed
# df.to_excel("compound_properties_by_name.xlsx", index=False)

input("Press Enter to exit...")
