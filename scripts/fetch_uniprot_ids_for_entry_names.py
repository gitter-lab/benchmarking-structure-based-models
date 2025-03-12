import requests
import pandas as pd


def fetch_uniprot_id(entry_name):
    print(entry_name)
    url = f"https://rest.uniprot.org/uniprotkb/{entry_name}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        # Extract the UniProt ID
        # print(data)
        uniprot_id = data.get('primaryAccession', None)
        return uniprot_id
    else:
        print(f"Failed to fetch data for {entry_name}. Status code: {response.status_code}")
        return None

def main():
    entry_name_to_uniprot_id_map = {}
    failed_assay_names = []
    # List of UniProt entry names
    entry_names = ['P12345', 'Q8NEE5', 'A0A0B4J3P0']  # Replace with your list of entry names

    # Fetch and print the UniProt IDs for each entry
    for entry_name in pd.read_csv('assets/DMS_substitutions_ref_file.csv')['UniProt_ID'].tolist():
        uniprot_id = fetch_uniprot_id(entry_name)
        if uniprot_id:
            entry_name_to_uniprot_id_map[entry_name] = uniprot_id
            print(f"UniProt ID for {entry_name}: {uniprot_id}")
        else:
            failed_assay_names.append(entry_name)
    
    pd.DataFrame(data={'entry_name': list(entry_name_to_uniprot_id_map.keys()), 
                       'uniprot_id': list(entry_name_to_uniprot_id_map.values())}
    ).to_csv('assets/entry_name_to_uniprot_id_map.csv', index=False)

if __name__ == '__main__':
    main()    