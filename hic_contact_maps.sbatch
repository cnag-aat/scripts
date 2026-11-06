#!/bin/bash
#SBATCH --job-name=hic_maps
#SBATCH --partition=research
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=hic_maps.%j.log
#
# hic_contact_maps.sbatch -- Hi-C contact maps from an assembly + paired Hi-C reads
#
#   chromap (--preset hic, PCR-dup removal) -> .pairs.gz
#     -> cooler cload/zoomify/balance     -> <prefix>.mcool  (5 kb base)
#     -> PretextMap / PretextSnapshot     -> <prefix>.pretext + png
#     -> plot_hic_maps.py                 -> labelled whole-genome + per-chromosome png/pdf
#
# Usage
#   sbatch hic_contact_maps.sbatch ASSEMBLY.fa[.gz] HIC_READ_DIR [OUTDIR]
#
#   e.g. sbatch hic_contact_maps.sbatch ibLobDeci3.2.asm.fa.gz hic_reads/ hic_maps
#
# Read pairs are found in HIC_READ_DIR by mate suffix:
#   X.1.fastq.gz/X.2.fastq.gz  X_1.fq.gz/X_2.fq.gz  X_R1.fastq.gz  X_R1_001.fastq.gz ...
#
# Options (environment variables; pass with  sbatch --export=ALL,BINSIZE=10000 ...)
#   BINSIZE     base (finest) resolution of the mcool, bp           [5000]
#   MAPQ        min MAPQ for chromap pairs                          [1]
#   RENAME      1 = simplify ENA/NCBI headers while decompressing   [1]
#               (CHR_1_OZ212012.1, MT_OZ212029.1, scaffold_68_CBAAFG020000022.1;
#                old->new table in <prefix>.names.tsv; see simplify_fasta_headers.py)
#   LABEL_MIN   size fallback: label seqs >= this many bp, used only if no
#               sequence is named like a chromosome (CHR_/Chr_/chr_/SUPER_)  [1000000]
#   PLOT_OPTS   extra options for plot_hic_maps.py                  [""]
#               e.g. PLOT_OPTS="--chroms-only --cmap fall"
#   PRETEXT_SORT  PretextMap --sortby: nosort|length|name           [nosort]
#   CONDA_ENV   conda env (name or prefix path) to activate if chromap is not on PATH
#               [/groups/assembly/software/conda/hic_contact_maps if it exists,
#                else hic_contact_maps]
#   THREADS     threads                                [SLURM_CPUS_PER_TASK or 8]
#   KEEP_TMP    1 = keep the uncompressed reference and chromap index  [0]
#
# Runs outside SLURM too:  bash hic_contact_maps.sbatch asm.fa.gz hic_reads/
# Conda env: conda env create -f hic_contact_maps.yml   (or -p /path/to/prefix)

set -euo pipefail

usage() { sed -n '/^# Usage/,/^# Read pairs/p' "$0" | grep -v '^# Read pairs' | sed 's/^# \{0,1\}//'; exit 1; }
[[ $# -ge 2 ]] || usage

[[ -s $1 ]] || { echo "ERROR: assembly '$1' not found (relative paths are from $(pwd))" >&2; exit 1; }
[[ -d $2 ]] || { echo "ERROR: read directory '$2' not found (relative paths are from $(pwd))" >&2; exit 1; }
ASM=$(realpath "$1")
READDIR=$(realpath "$2")
OUTDIR=${3:-hic_maps}

BINSIZE=${BINSIZE:-5000}
MAPQ=${MAPQ:-1}
LABEL_MIN=${LABEL_MIN:-1000000}
RENAME=${RENAME:-1}
PLOT_OPTS=${PLOT_OPTS:-}
PRETEXT_SORT=${PRETEXT_SORT:-nosort}
GROUP_ENV=/groups/assembly/software/conda/hic_contact_maps
if [[ -z ${CONDA_ENV:-} ]]; then
    if [[ -d $GROUP_ENV/conda-meta ]]; then CONDA_ENV=$GROUP_ENV; else CONDA_ENV=hic_contact_maps; fi
fi
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
KEEP_TMP=${KEEP_TMP:-0}

# plot_hic_maps.py and simplify_fasta_headers.py live next to this script. Under
# SLURM, $0 is a copy in the spool dir, so ask SLURM for the original path; also
# try the submit dir and PATH.
ORIG_SCRIPT=$0
if [[ -n ${SLURM_JOB_ID:-} ]] && command -v scontrol >/dev/null; then
    ORIG_SCRIPT=$(scontrol show job -o "$SLURM_JOB_ID" | tr ' ' '\n' | sed -n 's/^Command=//p' | head -n1)
    ORIG_SCRIPT=${ORIG_SCRIPT:-$0}
fi
SCRIPTDIR=
for d in "$(dirname "$(realpath "$ORIG_SCRIPT")")" "$(dirname "$(realpath "$0")")" \
         "$(realpath "${SLURM_SUBMIT_DIR:-.}")" \
         "$(dirname "$(command -v plot_hic_maps.py 2>/dev/null || echo /nonexistent/x)")"; do
    if [[ -x "$d/plot_hic_maps.py" && -x "$d/simplify_fasta_headers.py" ]]; then SCRIPTDIR=$d; break; fi
done
[[ -n "$SCRIPTDIR" ]] || { echo "ERROR: plot_hic_maps.py and simplify_fasta_headers.py not found next to $ORIG_SCRIPT, in ${SLURM_SUBMIT_DIR:-.}, or on PATH" >&2; exit 1; }

# --- environment -------------------------------------------------------------
if ! command -v chromap >/dev/null 2>&1; then
    CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "ERROR: chromap not on PATH and conda not found" >&2; exit 1; }
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    set +u; conda activate "$CONDA_ENV" || { echo "ERROR: cannot activate conda env '$CONDA_ENV' (set CONDA_ENV to an env name or prefix path)" >&2; exit 1; }; set -u
fi
for t in chromap samtools bgzip pigz cooler PretextMap PretextSnapshot python3; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done

log() { echo "[$(date '+%F %T')] $*"; }

# --- inputs ------------------------------------------------------------------
PREFIX=$(basename "$ASM"); PREFIX=${PREFIX%.gz}; PREFIX=${PREFIX%.fa}; PREFIX=${PREFIX%.fasta}; PREFIX=${PREFIX%.fna}
mkdir -p "$OUTDIR"/{ref,plots}
cd "$OUTDIR"

R1=(); R2=()
shopt -s nullglob
for f in "$READDIR"/*; do
    # <stem><sep>[R]1[_001].(fastq|fq).gz  ->  mate 2
    if [[ $(basename "$f") =~ ^(.*[._])(R?)1((_001)?\.f(ast)?q(\.gz)?)$ ]]; then
        m="$READDIR/${BASH_REMATCH[1]}${BASH_REMATCH[2]}2${BASH_REMATCH[3]}"
        [[ -e "$m" ]] || { echo "ERROR: no mate for $f (expected $m)" >&2; exit 1; }
        R1+=("$f"); R2+=("$m")
    fi
done
shopt -u nullglob
[[ ${#R1[@]} -gt 0 ]] || { echo "ERROR: no read-1 files found in $READDIR" >&2; exit 1; }

log "assembly : $ASM"
log "prefix   : $PREFIX"
log "read pairs (${#R1[@]}):"
for i in "${!R1[@]}"; do echo "    $(basename "${R1[$i]}")  $(basename "${R2[$i]}")"; done
log "outdir   : $(pwd)   threads: $THREADS   bin: $BINSIZE   MAPQ>=$MAPQ   label>=$LABEL_MIN bp"

# --- reference ---------------------------------------------------------------
# Steps are skipped if their output exists, so a failed job can simply be resubmitted.
REF=ref/$PREFIX.fa
PAIRS=$PREFIX.pairs.gz
get_ref() {
    if [[ ! -s $REF ]]; then
        log "decompressing reference"
        if [[ $RENAME == 1 ]]; then
            log "  renaming headers (ENA/NCBI -> CHR_<n>_<acc>, MT_<acc>, <scaffold>_<acc>); map: $PREFIX.names.tsv"
            pigz -dc -f -p 4 "$ASM" | python3 "$SCRIPTDIR/simplify_fasta_headers.py" - --map "$PREFIX.names.tsv" > "$REF"
        else
            pigz -dc -f -p 4 "$ASM" > "$REF"
        fi
    fi
    [[ -s $REF.fai ]] || samtools faidx "$REF"
}
[[ -s $REF.fai ]] || get_ref
cut -f1,2 "$REF.fai" > "$PREFIX.chrom.sizes"     # assembly order is kept everywhere

# --- mapping -> pairs --------------------------------------------------------
if [[ ! -s $PAIRS ]]; then
    get_ref
    if [[ ! -s ref/$PREFIX.chromap.idx ]]; then
        log "chromap index"
        chromap -i -t "$THREADS" -r "$REF" -o "ref/$PREFIX.chromap.idx"
    fi
    log "chromap mapping"
    chromap --preset hic --remove-pcr-duplicates -q "$MAPQ" -t "$THREADS" \
        -x "ref/$PREFIX.chromap.idx" -r "$REF" \
        -1 "$(IFS=,; echo "${R1[*]}")" -2 "$(IFS=,; echo "${R2[*]}")" \
        -o "$PREFIX.pairs" 2> chromap.log
    tail -n 30 chromap.log
    bgzip -@ "$THREADS" "$PREFIX.pairs"
fi

# --- cooler: 5 kb base, zoomify, balance ------------------------------------
# nice multiples of the base bin up to 10 Mb
RES=$(python3 -c "
b=$BINSIZE; out=[]
for e in range(0,5):
    for m in (1,2,5):
        r=b*m*10**e
        if r<=10_000_000 and r not in out: out.append(r)
print(','.join(map(str,sorted(out))))")

COOL=$PREFIX.${BINSIZE}.cool
MCOOL=$PREFIX.mcool
if [[ ! -s $MCOOL ]]; then
    log "cooler cload pairs @ $BINSIZE"
    pigz -dc -p 4 "$PAIRS" | \
        cooler cload pairs --assembly "$PREFIX" -c1 2 -p1 3 -c2 4 -p2 5 \
            "$PREFIX.chrom.sizes:$BINSIZE" - "$COOL"
    log "cooler zoomify + balance -> $RES"
    cooler zoomify --nproc "$THREADS" --balance --resolutions "$RES" -o "$MCOOL" "$COOL"
    rm -f "$COOL"
fi

# --- PretextMap --------------------------------------------------------------
# Rewrite the pairs header from the .fai so Pretext keeps assembly order
# (chromap writes #chromsize lines in its own order).
PRETEXT=$PREFIX.pretext
if [[ ! -s $PRETEXT ]]; then
    log "PretextMap"
    { echo "## pairs format v1.0"
      awk '{printf "#chromsize:\t%s\t%s\n", $1, $2}' "$REF.fai"
      printf "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2\n"
      pigz -dc -p 4 "$PAIRS" | grep -v '^#'
    } | PretextMap -o "$PRETEXT" --sortby "$PRETEXT_SORT" --mapq 0 > pretextmap.log 2>&1
fi
log "PretextSnapshot (full map)"
mkdir -p plots/pretext
PretextSnapshot -m "$PRETEXT" --sequences "=full" -r 4096 -f png -o plots/pretext > pretextsnapshot.log 2>&1 \
    || echo "[warn] PretextSnapshot failed (map is still usable in PretextView)" >&2

# --- publication plots -------------------------------------------------------
log "plot_hic_maps.py"
# shellcheck disable=SC2086
python3 "$SCRIPTDIR/plot_hic_maps.py" "$MCOOL" -o plots --prefix "$PREFIX" --min-size "$LABEL_MIN" $PLOT_OPTS

if [[ $KEEP_TMP != 1 ]]; then rm -f "ref/$PREFIX.chromap.idx" "$REF"; fi

log "done"
echo "  pairs   : $(pwd)/$PAIRS"
echo "  mcool   : $(pwd)/$MCOOL   (resolutions: $RES)"
echo "  pretext : $(pwd)/$PRETEXT"
echo "  plots   : $(pwd)/plots/"
