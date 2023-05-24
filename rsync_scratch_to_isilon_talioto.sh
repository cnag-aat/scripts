#!/bin/bash
#SBATCH --mem-per-cpu=6gb           # Memory per processor
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1           # Time limit hrs:min:sec
#SBATCH -q long

#### sbatch --array=0-3 ~/repos/scripts/rsync_scratch_to_isilon.sh
jobs=()
for s in `cat final_dirs_to_sync.txt`
do
jobs+=(${s})
done
d=${jobs[SLURM_ARRAY_TASK_ID]}
echo $d
if [[ -d /freezer/scratch/devel/talioto/denovo_assemblies/${d} ]]
then
    rsync -vrt --exclude='*.snakemake' /freezer/scratch/devel/talioto/denovo_assemblies/${d} /scratch_isilon/groups/assembly/data/archived_projects/assembly2/
fi
rsync -vrt --exclude='*.snakemake' $d /scratch_isilon/groups/assembly/data/archived_projects/assembly2/
