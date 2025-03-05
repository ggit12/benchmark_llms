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
    - LLM used for downstream analysis
4. Create a conda environment with `anndict`, `snakemake`, and `selenium`
```bash
conda create -n benchmark_llms python=3.12
conda activate benchmark_llms
conda install -c conda-forge tbb numba
pip install -r requirements.txt
```

## Pipeline Execution
Files are numbered sequentially to indicate processing order. Input and output of each script are defined in the Snakefile.

To run the pipeline:  
First `cd` into `benchmark_pipeline/`:
```bash
cd benchmark_pipeline/
```

Then run the pipeline with snakemake:
```bash
snakemake --snakefile Snakefile
```

Or, if you are running on a computing cluster with slurm, you can run the pipeline like this:
```bash
tmux new -s snakemake -d
tmux send-keys -t snakemake "cd $(pwd) && \
conda activate benchmark_llms && \
snakemake --snakefile Snakefile --profile slurm_profile" C-m
```

Note that this specific command might require some debugging, depending on your specific system configurations. For example, if conda is nto initialized in your `.bashrc`, this might fail. 
In this case, you can try manually by following the example below.
```bash
tmux new -s snakemake
cd /path/to/benchmark_llms/benchmark_pipeline

<<<Initialize conda here (however you normally would)>>>

conda activate benchmark_llms
snakemake --snakefile Snakefile --profile slurm_profile
```
Then `Ctrl+B`, `D` to detach from the session and leave it running in the background.
