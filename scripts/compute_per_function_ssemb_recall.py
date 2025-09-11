import pandas as pd
import matplotlib.pyplot as plt
import sys, os
from scipy.stats import spearmanr
import yaml
import numpy as np

def calc_toprecall(true_scores, model_scores, top_true=10, top_model=10):
    '''
        Code Taken from:
        https://github.com/OATML-Markslab/ProteinGym/blob/main/proteingym/performance_DMS_benchmarks.py
        For consistency with ProteinGym scores

        The code is available under an MIT License, copyrighted under:

        MIT License
        Copyright (c) 2023 OATML-Markslab, Pascal Notin, Aaron Kollasch, Daniel Ritter, Lood van Niekerk
    '''
    top_true = (true_scores >= np.percentile(true_scores, 100-top_true))
    top_model = (model_scores >= np.percentile(model_scores, 100-top_model))
    
    TP = (top_true) & (top_model)
    recall = TP.sum() / (top_true.sum()) if top_true.sum() > 0 else 0
    
    return (recall)


def main():
    with open("config/ssemb.yaml") as yaml_file:
        args = yaml.safe_load(yaml_file)

    proteingym_ref_file = args["proteingym_ref_file"]

    if not os.path.isfile(proteingym_ref_file):
        print(f"ProteinGym repository not found at path: {proteingym_ref_file}\n\
                Please run scripts/fetch_proteingym_assets.sh from the root directory of the project.")
        sys.exit(1)

    reference_df = pd.read_csv(proteingym_ref_file).rename(
        columns={"coarse_selection_type": "Selection Type"}
    )

    ssemb_score_file = (
        pd.read_csv(args["ssemb_score_file_path"])
        .dropna()
        .rename(columns={"dms_id": "DMS_id"})
    )

    dms_assays = ssemb_score_file["DMS_id"].unique().tolist()

    score_list = []
    for assay in dms_assays:
        score_list.append(
            calc_toprecall(
                ssemb_score_file[ssemb_score_file["DMS_id"] == assay]["score_dms"],
                ssemb_score_file[ssemb_score_file["DMS_id"] == assay]["score_ml"]
            )
        )

    ssemb_score_file = pd.DataFrame(
        data={
            "DMS_id": dms_assays,
            "ndcg": score_list,
        }
    )
    """
    ssemb_score_file["Spearman"] = spearmanr(
        ssemb_score_file["score_ml"], ssemb_score_file["score_dms"]
    ).statistic
    """

    model_score_df = ssemb_score_file[["DMS_id", "ndcg"]].merge(
        reference_df[["DMS_id", "UniProt_ID", "Selection Type"]],
        how="left",
        on="DMS_id",
    )

    average_scores_by_function = (
        model_score_df.groupby(["UniProt_ID", "Selection Type"])["ndcg"]
        .mean()
        .groupby("Selection Type")
        .mean()
        .round(4)
    )

    print(average_scores_by_function)

    print(f"Average spearman rho: {average_scores_by_function.mean()}")


if __name__ == "__main__":
    main()
