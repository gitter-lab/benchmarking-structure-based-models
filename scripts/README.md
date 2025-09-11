# Scripts

- `compute_per_function_model_scores.py`: Script to compute spearman rho scores for models already in proteinGym
- `compute_per_function_model_recall.py`: Script to compute Top 10 recall scores for models already in proteinGym.
- `compute_per_function_ssemb_scores.py`: Script to compute spearman rho scores for SSEmb model
- `compute_per_function_ssemb_recall.py`: Script to compute Top 10 recall scores for SSEmb model.
- `disordered_proteins.py`: Computes the UniProt IDs of proteins in ProteinGym and DisProt release
- `generate_box_plots.py`: Generate box plots for scores from all the models across 216 ProteinGym assays
- `compute_per_function_ensemble_scores.py`: Script to compute Top 10 recall/ Spearman scores for all 4 Ensemble models
- `fetch_proteingym_assets.sh`: Script to download DMS data from ProteinGym.
- `fetch_uniprot_ids_for_entry_names.py`: Map Uniprot ID for each protein in DisPort to its ID found in ProteinGym.
- `compute_disorder_in_proteingym.py`: Compute disordered regions in target proteins from ProteinGym using DisProt database annotations.
- `compute_protein_level_disorder.py`: Compute zero-shot spearman correlation scores for proteins that contain any mutation in discovered disordered regions.
- `setup_ssemb.sh`: Script to download model and data for SSEmb.
- `setup_proteinGym.sh`: Script to download and setup assets needed to compute spearman correlation scores using experimental structures (instead of AF-2 provided structures from ProteinGym)
