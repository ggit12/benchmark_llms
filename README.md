# Benchmarking LLMs Pipeline

This repository contains a pipeline for benchmarking Large Language Models (LLMs) introduced in: 
> #### Benchmarking Cell Type Annotation by Large Language Models with AnnDictionary  
> **George Crowley, Tabula Sapiens Consortium, Stephen R. Quake**  
> *bioRxiv* 2024.10.10.617605  
> [doi: https://doi.org/10.1101/2024.10.10.617605](https://doi.org/10.1101/2024.10.10.617605)

## Repository Structure

- `src/`: Contains core functions for
    - Plotting utilities
    - Benchmarking functions

- `benchmark_pipeline/`: Main pipeline directory
    - `Snakefile`: Defines the pipeline workflow and dependencies
    - `.env`: Configuration file for paths and API keys (create your own)
    - `scripts/`: Processing scripts in numerical order
    - `res/`: Output directory for results and figures

## Setup

1. Clone the repository
2. Create your own `.env` file in `benchmark_pipeline/` with:
     - Data paths
     - API keys
     - Other necessary configurations

## Pipeline Execution

Files are numbered sequentially to indicate processing order. Input/output relationships are defined in the Snakefile.

## Note

To run this pipeline, you must configure your own `.env` file with appropriate credentials and paths.
We provide a template file `.env_template`. To get started:

```bash
cp benchmark_pipeline/.env_template benchmark_pipeline/.env
```

Then modify the `.env` file in your preferred text editor with your specific configurations.
