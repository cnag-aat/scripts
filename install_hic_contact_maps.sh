#!/bin/bash
#
# install_hic_contact_maps.sh -- create/update the hic_contact_maps conda env and
# install the pipeline scripts into it, so they are on PATH once it is activated.
#
#   install_hic_contact_maps.sh [ENV_PREFIX]
#
#   e.g. install_hic_contact_maps.sh /groups/assembly/software/conda/hic_contact_maps
#        install_hic_contact_maps.sh                  # named env "hic_contact_maps"
#
# Re-run after pulling changes: the env is only updated if hic_contact_maps.yml
# changed, and the scripts are always re-copied.
#
# Installs into <env>/bin:   hic_contact_maps.sh  plot_hic_maps.py  simplify_fasta_headers.py
#        into <env>/share/hic_contact_maps:  hic_contact_maps.sbatch (copy it and
#        adjust #SBATCH resources for your cluster)

set -euo pipefail

SRC=$(dirname "$(realpath "$0")")
YML=$SRC/hic_contact_maps.yml
SCRIPTS=(hic_contact_maps.sh plot_hic_maps.py simplify_fasta_headers.py)
for f in "$YML" "${SCRIPTS[@]/#/$SRC/}" "$SRC/hic_contact_maps.sbatch"; do
    [[ -f $f ]] || { echo "ERROR: $f not found" >&2; exit 1; }
done

CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "ERROR: conda not found" >&2; exit 1; }
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
CONDA=conda; command -v mamba >/dev/null 2>&1 && CONDA=mamba

if [[ $# -ge 1 ]]; then
    PREFIX=$(realpath -m "$1" 2>/dev/null || python3 -c 'import os,sys;print(os.path.abspath(sys.argv[1]))' "$1")
    TARGET=(-p "$PREFIX")
else
    PREFIX=$CONDA_BASE/envs/hic_contact_maps
    TARGET=(-n hic_contact_maps)
fi

if [[ ! -d $PREFIX/conda-meta ]]; then
    echo "creating env $PREFIX"
    $CONDA env create -y "${TARGET[@]}" -f "$YML"
elif ! cmp -s "$YML" "$PREFIX/share/hic_contact_maps/hic_contact_maps.yml"; then
    echo "updating env $PREFIX (hic_contact_maps.yml changed)"
    $CONDA env update "${TARGET[@]}" -f "$YML" --prune
else
    echo "env $PREFIX is up to date"
fi

mkdir -p "$PREFIX/share/hic_contact_maps"
for s in "${SCRIPTS[@]}"; do
    install -m 755 "$SRC/$s" "$PREFIX/bin/$s"
    echo "installed $PREFIX/bin/$s"
done
install -m 644 "$SRC/hic_contact_maps.sbatch" "$YML" "$PREFIX/share/hic_contact_maps/"
echo "installed $PREFIX/share/hic_contact_maps/hic_contact_maps.sbatch"

GIT=$(git -C "$SRC" describe --always --dirty 2>/dev/null || echo unknown)
echo "$GIT" > "$PREFIX/share/hic_contact_maps/VERSION"
echo
echo "done (scripts from $SRC @ $GIT). Use:"
echo "  conda activate ${TARGET[1]}"
echo "  sbatch $PREFIX/share/hic_contact_maps/hic_contact_maps.sbatch ASSEMBLY.fa.gz HIC_READ_DIR [OUTDIR]"
