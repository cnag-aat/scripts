#!/bin/bash
#SBATCH --mem-per-cpu=6gb           # Memory per processor
#SBATCH --time=6:00:00              # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1           # Time limit hrs:min:sec

#### sbatch --array=0-60 ~/repos/scripts/rsync_scratch_to_isilon_annotationdir.sh

jobs=()
for s in `cat dir_list_to_sync.txt`
do
jobs+=(${s})
done
d=${jobs[SLURM_ARRAY_TASK_ID]}
echo $d
if [[ -d /freezer/scratch/devel/talioto/de_novo_annotation/${d} ]]
then
    rsync -av --exclude='*.snakemake' /freezer/scratch/devel/talioto/de_novo_annotation/${d} /scratch_isilon/groups/assembly/data/archived_projects/annotation/
fi
rsync -av --exclude='*.snakemake' $d /scratch_isilon/groups/assembly/data/archived_projects/annotation/