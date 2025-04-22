"""
The functions in this module handle running multiple models and providers on an adata_dict.
"""
# pylint: disable=line-too-long
import time

import anndict as adt


def process_adata_dict(
    adata_dict: adt.AdataDict,
    llm_config: dict,
    initial_resolution: int = 0.5, # pylint: disable=unused-argument
    subcluster_resolution: int = 0.2, # pylint: disable=unused-argument
    expected_cell_type_col: str = None,
) -> tuple[dict, dict, dict, dict]:
    """
    Process a dictionary of AnnData objects through a series of AI-powered cell type annotation steps.
    This function performs multiple steps of cell type annotation and analysis using AI:
    1. Configures the LLM backend with provided settings
    2. Performs initial cell type annotation
    3. Simplifies cell type labels by removing redundancy

    Parameters
    ----------
    adata_dict
        Dictionary of AnnData objects to process.
        Keys are identifiers and values are AnnData objects.

    llm_config
        Configuration dictionary for the Language Learning Model (LLM) backend.
        Must contain 'model' key specifying which model to use.

    initial_resolution
        Initial resolution parameter for leiden clustering, by default 0.5

    subcluster_resolution
        Resolution parameter for subclustering, by default 0.2

    expected_cell_type_col
        Column name in the AnnData object that contains expected cell types.

    Returns
    -------
    A dictionary containing:
    - obs_dict (dict): Dictionary of observation DataFrames from processed AnnData objects
    - appropriate_resolution_dict (dict or None): Dictionary of appropriate resolutions for each dataset
    - label_results (dict): Results from initial AI cell type annotation
    - simplified_mappings (dict): Mapping of original to simplified cell type labels

    Notes
    -----
    The function currently has some commented-out functionality for subclustering and 
    further annotation refinement which can be uncommented if needed.

    """

    #Configure the backend to work with the specific provider and model
    adt.configure_llm_backend(**llm_config)

    #Extract model name from the config
    model = llm_config['model']

    #Determine appropriate cluster resolutions using AI
    #This will leave the final column as 'leiden' in the .obs of each anndata
    #appropriate_resolution_dict = adt.ai_determine_leiden_resolution_adata_dict(adata_dict, initial_resolution=initial_resolution)
    appropriate_resolution_dict = None

    #Get the expected cell types for each adata
    def get_expected_cell_types(adata):
        return list(adata.obs[expected_cell_type_col].unique())
    expected_cell_types = adata_dict.fapply(get_expected_cell_types)

    #Use AI to do automatic interpretation of Differentially Expressed Genes
    #This will be a first-draft annotation that we will subcluster, reannotate, then re-merge
    # label_results = adt.wrappers.ai_annotate_cell_type_adata_dict(adata_dict, groupby='leiden', n_top_genes=10, new_label_column=f'{model}_ai_cell_type', tissue_of_origin_col="tissue")
    label_results = adata_dict.fapply(adt.ai_annotate_from_expected_cell_types, groupby='leiden', n_top_genes=10, expected_cell_types=expected_cell_types, new_label_column=f'{model}_ai_cell_type', tissue_of_origin_col="tissue")
    # print(label_results)
    print('ai cell type found')

    #These labels seem to have some redundancy, let's merge them with AI
    # simplified_mappings = adt.wrappers.simplify_obs_column_adata_dict(adata_dict, f'{model}_ai_cell_type', f'{model}_simplified_ai_cell_type', simplification_level='redundancy-removed')
    # print('simplified ai cell type found')

    obs_dict = {key: adata.obs for key, adata in adata_dict.items()}

    return obs_dict, appropriate_resolution_dict, label_results

def run_multiple_providers_models(adata_dict, provider_endpoint_dict, provider_config, expected_cell_type_col):
    """
    Run multiple language model providers and models to process annotated data.
    This function processes a dictionary of AnnData objects through multiple LLM providers
    and their respective models, applying the specified configurations for each provider.

    Parameters
    ----------

    adata_dict (dict): Dictionary containing AnnData objects to be processed
    provider_endpoint_dict (dict): Dictionary mapping providers to their available models
    provider_config (dict): Configuration settings for each provider
    expected_cell_type_col (str): Column name in the .obs of each AnnData in ``adata_dict`` that contains expected cell types

    Returns
    -------
    dict: Dictionary containing results for each provider-model combination with the following structure:
        {
            'provider_model': {
                'obs_dict': Dictionary of observations,
                'appropriate_resolution': Dictionary of appropriate resolutions,
                'label_results': Results of labeling process,
                'simplified_mappings': Simplified cluster mappings
    Raises
    ------
    ValueError: If configuration for a specified provider is not found in provider_config

    Notes
    -----
    - Includes a 60-second sleep between model runs to prevent API rate limiting
    - Model specifications in provider_config will be overridden by the model parameter
        from provider_endpoint_dict if they differ

    """
    results = {}
    for provider, models in provider_endpoint_dict.items():
        if provider not in provider_config:
            raise ValueError(f"Configuration for provider '{provider}' not found in provider_config")

        for model in models:
            key = f'{provider}_{model}'
            print(f"Key: {key}")

            llm_config = provider_config[provider].copy()

            # Ensure the model is correctly set in the configuration
            if 'model' in llm_config and llm_config['model'] != model:
                print(f"Warning: Overriding model in config from '{llm_config['model']}' to '{model}'")
            llm_config['model'] = model


            # obs_dict, appropriate_resolution_dict, label_results, simplified_mappings, sub_cluster_mappings, simplified_sub_cluster_mappings = process_adata_dict(adata_dict, llm_config)
            obs_dict, appropriate_resolution_dict, label_results = process_adata_dict(adata_dict, llm_config, expected_cell_type_col=expected_cell_type_col)
            results[key] = {
                'obs_dict' : obs_dict,
                'appropriate_resolution': appropriate_resolution_dict,
                'label_results': label_results,
            }

            #sleep 1 min in between models to prevent API overload (rate limiter is reset each time a new model is initiated)
            time.sleep(60)

    return results
