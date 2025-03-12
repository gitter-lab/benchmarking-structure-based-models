import argparse


def compute_disorder_region_in_assay(disprot_proteins_df):
    disordered_regions = []
    for i, row in disprot_proteins_df.iterrows():
        # print(row['consensus'])
        regions = []
        start_idx = None
        end_idx = None
        for i, c in enumerate(row['consensus']):
            if c == 'D' or c == 'T':
                if start_idx is None:
                    start_idx = i
                    end_idx = i
                else:
                    end_idx = i
            else:
                if start_idx is not None:
                    regions.append(f"{start_idx}-{end_idx}")
                    start_idx = None
        if start_idx is not None:
            regions.append(f"{start_idx}-{end_idx}")
        disordered_regions.append(','.join(regions))
    return disordered_regions



def main(args):
    disprot_df = pd.read_csv(args['disprot_data_path'], 
                              sep='\t').drop(columns=['obsolete'])
    proteingym_uniprot_id_df = pd.read_csv(args['unprotID_to_entry_name_file_path'])

    disprot_proteins_df = proteingym_uniprot_id_df.merge(
        disprot_df[['acc', 'consensus']],
        right_on='acc', left_on='uniprot_id'
    )

    disordered_regions = compute_disorder_region_in_assay(disprot_proteins_df)

    ref_file_df = pd.read_csv('assets/DMS_substitutions_ref_file.csv')

    ref_file_df = ref_file_df.merge(
        disprot_proteins_df[['entry_name', 'disordered_regions', 'uniprot_id']],
        left_on='UniProt_ID',
        right_on='entry_name',
        how='inner'
    ).drop(columns=['entry_name'])

    '''
        Filter disordered regions to reflect regions targeted by assays, i.e.
        throw away disordered regions not in the target assay.
    '''
    for i in ref_file_df.index:
        disordered_regions = ref_file_df.iloc[i]['disordered_regions'].split(',')
        filtered_regions = []
        for r in disordered_regions:
            if int(r.split('-')[0]) >= ref_file_df.iloc[i]['seq_len']:
                continue
            elif int(r.split('-')[1]) >= ref_file_df.iloc[i]['seq_len']:
                filtered_regions.append(f"{r.split('-')[0]}-{ref_file_df.iloc[i]['seq_len']-1}")
            else:
                filtered_regions.append(r)
        ref_file_df.loc[i, 'disordered_regions'] = ','.join(filtered_regions)
    ref_file_df.to_csv('assets/disordered_DMS_substitutions_ref_file.csv')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--disprot_data_path", type=str, 
                        default='assets/DisProt_release_2024_12_consensus_regions.tsv')
    parser.add_argument("--unprotID_to_entry_name_file_path", type=str, 
                        default='assets/entry_name_to_uniprot_id_map.csv')
    args = parser.parse_args()

    main(args)