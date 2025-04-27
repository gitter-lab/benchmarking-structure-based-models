import argparse
import requests, sys, json
import pandas as pd

def compute_disorder_region_in_assay(disprot_proteins_df):
    combined_df_data = {}
    for i, r in disprot_proteins_df.iterrows():
        combined_df_data[r['entry_name']] = r['uniprot_id']

    combined_df = pd.DataFrame(data={'entry_name': list(combined_df_data.keys()), 
                                        'uniprot_id': list(combined_df_data.values())})
    combined_df['mutated_sequences'] = ""
    combined_df['disorder_regions'] = ""
    for i, r in combined_df.iterrows():
        combined_df.loc[i, 'mutated_sequences'] = ",".join(
            row['region_sequence'] 
        for i, row in disprot_proteins_df[disprot_proteins_df['entry_name'] == r['entry_name']].iterrows())
        combined_df.loc[i, 'disorder_regions'] = ",".join(
            str(row['start']) + "-" + str(row['end']) 
        for i, row in disprot_proteins_df[disprot_proteins_df['entry_name'] == r['entry_name']].iterrows())
    return combined_df

def longest_common_substring(string1, string2):
    len1, len2 = len(string1), len(string2)
    longest_substring = ""
    for i in range(len1):
        for j in range(len2):
            current_substring = ""
            k = 0
            while i + k < len1 and j + k < len2 and string1[i + k] == string2[j + k]:
                current_substring += string1[i + k]
                k += 1
            if len(current_substring) > len(longest_substring):
                longest_substring = current_substring
    return longest_substring

def fetch_uniprot_sequence(uniprot_acc):
    headers = {
    "accept": "application/json"
    }
    base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}"

    response = requests.get(base_url, headers=headers)
    if not response.ok:
        response.raise_for_status()
        sys.exit()

    data = response.json()
    return data['sequence']['value']


def main(args):
    disprot_df = pd.read_csv(args['disprot_data_path'], 
                              sep='\t').drop(columns=['obsolete'])
    proteingym_uniprot_id_df = pd.read_csv(args['unprotID_to_entry_name_file_path'])

    disprot_proteins_df = proteingym_uniprot_id_df.merge(
        disprot_df[['acc', 'start', 'end', 'region_sequence']], 
        right_on='acc', left_on='uniprot_id', )

    combined_disorder_df = compute_disorder_region_in_assay(disprot_proteins_df)

    ref_file_df = pd.read_csv('assets/DMS_substitutions_ref_file.csv')
    ref_file_df = ref_file_df.merge(
        combined_disorder_df[['entry_name', 'mutated_sequences', 'disorder_regions', 'uniprot_id']], 
        left_on='UniProt_ID', right_on='entry_name', how='inner'
    ).drop(columns=['entry_name'])

    '''
        Filter disordered regions to reflect regions targeted by assays, i.e.
        throw away disordered regions not in the target assay.
    '''
    ref_file_df['disordered_regions'] = ""

    get_region_start_end = lambda x: (int(x.split('-')[0]), int(x.split('-')[1]))

    for i in ref_file_df.index:
        uniprot_seq = fetch_uniprot_sequence(ref_file_df.loc[i, 'uniprot_id'])
        longest_seq_match = longest_common_substring(uniprot_seq, ref_file_df.loc[i, 'target_seq'])

        if len(uniprot_seq) == len(ref_file_df.loc[i, 'target_seq']):
            target_seq_start = 0
        else:
            target_seq_start = uniprot_seq.find(longest_seq_match)
        target_seq_end = target_seq_start + len(ref_file_df.loc[i, 'target_seq']) - 1

        try:
            disordered_regions = ref_file_df.loc[i, 'disorder_regions'].split(',')
        except:
            disordered_regions = [ref_file_df.loc[i, 'disorder_regions']]
        filtered_regions = []
        for r in disordered_regions:
            r_start, r_end = get_region_start_end(r)
            
            # Convert region start/end to 0 indexed
            r_start -= 1
            r_end -= 1

            if (r_end <  target_seq_start) or (r_start >  target_seq_end):
                continue
            reg = str(max(r_start,  target_seq_start)) +"-" + \
                        str(min(target_seq_end, r_end))
            
            found_region = False
            for j in range(len(filtered_regions)):

                if (int(filtered_regions[j].split('-')[0]) > (int(reg.split('-')[1]) - target_seq_start)) or (int(filtered_regions[j].split('-')[1]) < (int(reg.split('-')[0]) - target_seq_start)):
                    continue

                filtered_regions[j] = str(min(int(filtered_regions[j].split('-')[0]), int(reg.split('-')[0]) - target_seq_start)) + "-" + \
                    str(max(int(filtered_regions[j].split('-')[1]), int(reg.split('-')[1]) - target_seq_start))
                #print(min(int(filtered_regions[j].split('-')[0]), int(reg.split('-')[0])))
                #print(max(int(filtered_regions[j].split('-')[1]), int(reg.split('-')[1])) - target_seq_start)
                #print(filtered_regions[j])
                found_region = True
                break
            if not found_region:
                filtered_regions.append(str(int(reg.split('-')[0]) - target_seq_start) + "-" + str(int(reg.split('-')[1]) - target_seq_start))
        if len(filtered_regions) > 0:    
            ref_file_df.loc[i, 'disordered_regions'] = ','.join(filtered_regions)

    ref_file_df.drop(columns=['disorder_regions'], inplace=True)
    ref_file_df.to_csv('assets/disordered_DMS_substitutions_ref_file.csv')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--disprot_data_path", type=str, 
                        default='assets/DisProt_release_2024_12.tsv')
    parser.add_argument("--unprotID_to_entry_name_file_path", type=str, 
                        default='assets/entry_name_to_uniprot_id_map.csv')
    args = vars(parser.parse_args())

    main(args)