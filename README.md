# Benchmarking LLMs Pipeline

This repository contains a pipeline for benchmarking Large Language Models (LLMs) at cell type annotation introduced in: 
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
    - `scripts/`: Processing scripts in labeled in the order that they are run
    - `res/`: Output directory for results and figures

## Setup

1. Clone the repository
```bash
git clone https://github.com/ggit12/benchmark_llms.git
```
2. Create a local `.env` file by copying the template:
```bash
cd benchmark_llms
cp benchmark_pipeline/.env_template benchmark_pipeline/.env
```
3. Edit your `.env` file in `benchmark_pipeline/` as needed, including:
    - LLM provider API keys
    - Data path
4. Create a conda environment with `anndict`, `snakemake`, and `selenium`
```bash
conda create -n benchmark_llms python=3.12
conda activate benchmark_llms
conda install -c conda-forge tbb numba
pip install -r requirements.txt
```

## Pipeline Execution

Run the pipeline from within `benchmark_pipeline/`:
```bash
cd benchmark_pipeline/
snakemake --snakefile Snakefile
```

Files are numbered sequentially to indicate processing order. Input and output of each script are defined in the Snakefile.
