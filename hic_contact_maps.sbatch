#!/bin/bash
#SBATCH --job-name=hic_maps
#SBATCH --partition=research
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=hic_maps.%j.log
#
# hic_contact_maps.sbatch -- SLURM wrapper for hic_contact_maps.sh
#
#   sbatch hic_contact_maps.sbatch ASSEMBLY.fa[.gz] HIC_READ_DIR [OUTDIR]
#
# Only scheduler resources and environment set-up live here; the pipeline and
# its options (BINSIZE, MAPQ, PLOT_OPTS, ...) are documented in hic_contact_maps.sh.
# Pass options with:  sbatch --export=ALL,BINSIZE=10000,PLOT_OPTS="--cmap fall" ...
# Override resources on the command line, e.g.  sbatch --partition=research_long --mem=200G ...
#
#   CONDA_ENV   env (name or prefix path) activated if hic_contact_maps.sh is not
#               already on PATH  [/groups/assembly/software/conda/hic_contact_maps
#               if it exists, else hic_contact_maps]

set -euo pipefail

if ! command -v hic_contact_maps.sh >/dev/null 2>&1; then
    GROUP_ENV=/groups/assembly/software/conda/hic_contact_maps
    if [[ -z ${CONDA_ENV:-} ]]; then
        if [[ -d $GROUP_ENV/conda-meta ]]; then CONDA_ENV=$GROUP_ENV; else CONDA_ENV=hic_contact_maps; fi
    fi
    CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "ERROR: conda not found" >&2; exit 1; }
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    set +u
    conda activate "$CONDA_ENV" || { echo "ERROR: cannot activate conda env '$CONDA_ENV' (set CONDA_ENV to an env name or prefix path)" >&2; exit 1; }
    set -u
fi

PIPELINE=$(command -v hic_contact_maps.sh || true)
# not installed in the env: fall back to a copy in the submit dir
[[ -n $PIPELINE ]] || { [[ -x ${SLURM_SUBMIT_DIR:-.}/hic_contact_maps.sh ]] && PIPELINE=${SLURM_SUBMIT_DIR:-.}/hic_contact_maps.sh; }
[[ -n $PIPELINE ]] || { echo "ERROR: hic_contact_maps.sh not found in env '${CONDA_ENV:-}' or submit dir (run install_hic_contact_maps.sh)" >&2; exit 1; }

echo "[$(date '+%F %T')] job ${SLURM_JOB_ID:-local} on $(hostname): $PIPELINE $*"
exec "$PIPELINE" "$@"
