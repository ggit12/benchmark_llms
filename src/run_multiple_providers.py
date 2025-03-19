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
    subcluster_resolution: int = 0.2 # pylint: disable=unused-argument
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

    #Use AI to do automatic interpretation of Differentially Expressed Genes
    #This will be a first-draft annotation that we will subcluster, reannotate, then re-merge
    label_results = adt.wrappers.ai_annotate_cell_type_adata_dict(adata_dict, groupby='leiden', n_top_genes=10, new_label_column=f'{model}_ai_cell_type', tissue_of_origin_col="tissue")
    print(label_results)
    print('ai cell type found')

    #These labels seem to have some redundancy, let's merge them with AI
    simplified_mappings = adt.wrappers.simplify_obs_column_adata_dict(adata_dict, f'{model}_ai_cell_type', f'{model}_simplified_ai_cell_type', simplification_level='redundancy-removed')
    print('simplified ai cell type found')

    #Now get cell type labels with more granularity
    #First, subcluster within each cell type
    #Resolution is passed as a single number here, but could be a dict with keys adata_dict.keys() with values of resolution, or dictionary with keys adata_dict['some_key'].obs[groupby].unique() and values as numbers
    # adata_dict = adt.leiden_sub_cluster_adata_dict(adata_dict, groupby=f'{model}_simplified_ai_cell_type', key_added=f'{model}_leiden_subcluster', resolution=subcluster_resolution)
    # print('adata_dict has been subclustered')

    #Now, annotate each subcluster, passing tissue info as well
    # adata_dict, sub_cluster_mappings = adt.ai_annotate_cell_sub_type_adata_dict(adata_dict, n_top_genes=15, cell_type_col=[f'{model}_simplified_ai_cell_type'], sub_cluster_col=f'{model}_leiden_subcluster', new_label_column=f'{model}_ai_cell_sub_type', tissue_of_origin_col='tissue')
    # print("adata_dict has had subclusters annotated")

    #Use AI to remove redundancy again
    # simplified_sub_cluster_mappings = adt.simplify_obs_column_adata_dict(adata_dict, f'{model}_ai_cell_sub_type', f'{model}_simplified_ai_cell_sub_type', simplification_level='redundancy-removed, keep periods in keys')
    # print("adata_dict has had subcluster annotations simplified")

    obs_dict = {key: adata.obs for key, adata in adata_dict.items()}

    return obs_dict, appropriate_resolution_dict, label_results, simplified_mappings#, sub_cluster_mappings#, simplified_sub_cluster_mappings

def run_multiple_providers_models(adata_dict, provider_endpoint_dict, provider_config):
    """
    Run multiple language model providers and models to process annotated data.
    This function processes a dictionary of AnnData objects through multiple LLM providers
    and their respective models, applying the specified configurations for each provider.
    
    Parameters
    ----------

    adata_dict (dict): Dictionary containing AnnData objects to be processed
    provider_endpoint_dict (dict): Dictionary mapping providers to their available models
    provider_config (dict): Configuration settings for each provider
    
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
            obs_dict, appropriate_resolution_dict, label_results, simplified_mappings = process_adata_dict(adata_dict, llm_config)
            results[key] = {
                'obs_dict' : obs_dict,
                'appropriate_resolution': appropriate_resolution_dict,
                'label_results': label_results,
                'simplified_mappings': simplified_mappings,
                # 'sub_cluster_mappings': sub_cluster_mappings,
                # 'simplified_sub_cluster_mappings': simplified_sub_cluster_mappings
            }

            #sleep 1 min in between models to prevent API overload (rate limiter is reset each time a new model is initiated)
            time.sleep(60)

    return results
