#!/bin/sh
# Assign input arguments to variables
# Sample command line: ./run_hifiasm_ont.sh hifiasm marathon_assembly \
# /scratch_isilon/groups/assembly/data/projects/BGE/qqOdiTrog/assembly/s02.1_p01.1_Filtlong/ont.filtlong.fastq.gz \
# 128 qqOdiTrog_hifiasm_ont /scratch_isilon/groups/assembly/data/projects/BGE/qqOdiTrog/asm2_hifiasm/s03.3_p02.1_hifiasm 24:00:00 qqOdiTrog

PROGRAM=${1:-hifiasm}
QUEUE=${2:-marathon_assembly}
ONT=${3:-/scratch_isilon/groups/assembly/data/projects/BGE/qqOdiTrog/assembly/s02.1_p01.1_Filtlong/ont.filtlong.fastq.gz}
THREADS_TOTAL=${4:-128}
PREFIXOUT=${5:-qqOdiTrog_hifiasm_ont}
INPUTDIR=${6:-/scratch_isilon/groups/assembly/data/projects/BGE/qqOdiTrog/asm2_hifiasm/s03.3_p02.1_hifiasm}
TIME_TOTAL=${7:-24:00:00}
TOLID=${8:-qqOdiTrog}
HIC1=${9:-}
HIC2=${10:-}

# Check if READ1 file exists
if [ -n "$HIC1" ] && [ ! -f "$HIC1" ]; then
  HIC1=""
fi

# Check if READ2 file exists
if [ -n "$HIC2" ] && [ ! -f "$HIC2" ]; then
  HIC2=""
fi


# Echo the variables to check their values
echo "PROGRAM: $PROGRAM"
echo "TOLID: $TOLID"
echo "Threads: $THREADS_TOTAL"
echo "Time: $TIME_TOTAL"

# Create a temporary SLURM script
TEMP_SLURM_SCRIPT="$INPUTDIR/${PROGRAM}_${THREADS_TOTAL}_$(date +%s).slurm"

cat << EOF > $TEMP_SLURM_SCRIPT
#!/bin/bash
#SBATCH --job-name=${PROGRAM}.${TOLID}.${THREADS_TOTAL}
#SBATCH --qos=$QUEUE
#SBATCH --partition=general
#SBATCH --requeue
#SBATCH --output=${PROGRAM}.${TOLID}.${THREADS_TOTAL}.%j.out
#SBATCH --error=${PROGRAM}.${TOLID}.${THREADS_TOTAL}.%j.err
#SBATCH --mem=950GB
#SBATCH -c $THREADS_TOTAL
#SBATCH --time=$TIME_TOTAL


#2. Activate environment:

source ~/.bashrc;

conda activate /software/assembly/conda/hifiasm0.24.0-r702/

#3. Run with HiC Phasing

cd $INPUTDIR

EOF

# Add the appropriate hifiasm command based on the availability of HiC reads

if [ -n "$HIC1" ] && [ -n "$HIC2" ]; then
cat << EOF >> $TEMP_SLURM_SCRIPT
hifiasm -o $PREFIXOUT.asm -t $THREADS_TOTAL --ont --h1 $HIC1 --h2 $HIC2 $ONT 2> $PREFIXOUT.asm.log
EOF
else
cat << EOF >> $TEMP_SLURM_SCRIPT
hifiasm -o $PREFIXOUT.asm -t $THREADS_TOTAL --ont $ONT 2> $PREFIXOUT.asm.log
EOF
fi

sbatch $TEMP_SLURM_SCRIPT


