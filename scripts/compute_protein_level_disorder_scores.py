import argparse
import requests, sys, json
import pandas as pd
import math
from copy import deepcopy
from scipy.stats import spearmanr
import numpy as np

from matplotlib import pyplot as plt

models_2_path_map = {'ESM2': ('ESM2/650M/', 'esm2_t33_650M_UR50D'),
          'ESM-IF1': ('ESM-IF1/', 'esmif1_ll'),
          'ProtSSN': ('ProtSSN', 'ProtSSN_ensemble'),
          'TranceptEVE': ('TranceptEVE/TranceptEVE_L/', 'avg_score'),
          }
def parse_disorder_scores_for_models(model_name, ref_file_df, args):
    model_scores = {}
    dms_scores = {}
    for i,row in ref_file_df.iterrows():
        model_scores[row['DMS_id']] = pd.read_csv(f"{args['zero_shot_scores']}/{models_2_path_map[model_name][0]}/{row['DMS_filename']}")
        dms_scores[row['DMS_id']] = pd.read_csv(f"{args['DMS_substitutions_assays_path']}/{row['DMS_filename']}")[['mutant', 'DMS_score']]



    model_disordered_scores = {}
    for i,row in ref_file_df.iterrows():
        if row['DMS_id'] not in model_scores:
            continue

        model_disordered_scores[row['DMS_id']] = model_scores[row['DMS_id']].copy()
        model_disordered_scores[row['DMS_id']]['is_disordered'] = [False]* len(model_disordered_scores[row['DMS_id']])

        if 'DMS_score' not in model_disordered_scores[row['DMS_id']].columns:
            model_disordered_scores[row['DMS_id']] = model_disordered_scores[row['DMS_id']].merge(dms_scores[row['DMS_id']][['mutant', 'DMS_score']], how='inner', on='mutant')
        if isinstance(row['disordered_regions'], float):
            continue
        for d_range in row['disordered_regions'].split(','):
            upper_limit = int(d_range.split('-')[1])
            low_limit = int(d_range.split('-')[0])
            model_disordered_scores[row['DMS_id']]['is_disordered'] = model_disordered_scores[row['DMS_id']]['is_disordered'] | \
            model_disordered_scores[row['DMS_id']]['mutant'].str.split(":").apply(lambda x: any([int(s[1:-1])-1 <= upper_limit and int(s[1:-1])-1 >= low_limit for s in x]))
    
    scores = []
    for i in range(len(ref_file_df)):
        scores.append(spearmanr(model_disordered_scores[ref_file_df.iloc[i]['DMS_id']]['DMS_score'],
            model_disordered_scores[ref_file_df.iloc[i]['DMS_id']][models_2_path_map[model_name][1]]).statistic)
    ref_file_df['spearman score'] = scores
    return ref_file_df

def main(args):
    disorder_ref_file_df = pd.read_csv(
        args['disorder_ref_file']
    )
    combined_spearman_df = disorder_ref_file_df
    model_scores = {}

    for mname in models_2_path_map.keys():
        model_scores[mname] = parse_disorder_scores_for_models(mname, combined_spearman_df.copy(), args)

    barWidth = 1/(len(models_2_path_map) + 1)

    brs = {}
    mnames = list(models_2_path_map.keys())
    print()
    brs[mnames[0]] = np.arange(len(model_scores[mnames[0]].groupby(['UniProt_ID', 'coarse_selection_type'])['spearman score'].mean().groupby('coarse_selection_type').mean().index)) 
    for i in range(1, len(mnames)):
        brs[mnames[i]] = [x + i*barWidth for x in brs[mnames[0]]] 

    for mname in models_2_path_map.keys():
        plt.bar(brs[mname], 
        model_scores[mname].groupby(['UniProt_ID', 'coarse_selection_type'])['spearman score'].mean().groupby('coarse_selection_type').mean().tolist(),
        width = barWidth, edgecolor ='grey', label =mname, alpha=0.8)


    plt.xticks([r + barWidth for r in range(len(model_scores[mnames[0]].groupby(['UniProt_ID', 'coarse_selection_type'])['spearman score'].mean().groupby('coarse_selection_type').mean().index))], 
            model_scores[mnames[0]].groupby(['UniProt_ID', 'coarse_selection_type'])['spearman score'].mean().groupby('coarse_selection_type').mean().index.tolist(),
            rotation=45)

    plt.ylabel('Spearman Correlation for disordered proteins', fontweight ='bold', fontsize = 10) 
    plt.xlabel('Function type', fontweight ='bold', fontsize = 15) 

    plt.legend()
    plt.show()
    plt.show()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--disorder_ref_file", type=str, 
                        default='assets/disordered_DMS_substitutions_ref_file.csv')
    parser.add_argument("--zero_shot_scores", type=str, required=True)
    parser.add_argument("--DMS_substitutions_assays_path", type=str, required=True)
    #parser.add_argument("--unprotID_to_entry_name_file_path", type=str, 
    #                    default='assets/entry_name_to_uniprot_id_map.csv')
    args = vars(parser.parse_args())

    main(args)