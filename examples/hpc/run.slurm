#!/usr/bin/env bash
#SBATCH --job-name=funlr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=funlr-%j.log

# Example allocation, not a measured production resource requirement.
# Add your site's account/partition directives and tune memory/time to your data.
# Submit with FUNLR_PREFIX and FUNLR_CONFIG set to absolute paths.
set -euo pipefail
: "${FUNLR_PREFIX:?Set FUNLR_PREFIX to the installed environment path}"
: "${FUNLR_CONFIG:?Set FUNLR_CONFIG to the YAML configuration path}"
export PATH="${FUNLR_PREFIX}/bin:${PATH}"
funlr run --config "${FUNLR_CONFIG}" --threads "${SLURM_CPUS_PER_TASK:-1}"
