#! /bin/bash


# User: fcruz
# Date: 2024-07-18

# Brief description/notes:
#   1. This script will create an sbatch per fastq to perform the ENA upload using our credentials 
#   USAGE: go to the <tolid>/ENA/reads/ont directory and run:

#       ena_submit_lftp_fastQ_upload_sbatch.sh tolid

# Writing the .sbatch file

tolid=$1;

for file in *.fastq.gz; 

  do     
    #basename per fastq
    b=$(basename $file .gz);
    sbatchfile=$tolid"_upload_"$b".sbatch";
    #echo "$sbatchfile";
    echo -e '#!/usr/bin/env bash' > $sbatchfile;
    echo "">> $sbatchfile;
    echo -e '#SBATCH --cpus-per-task=4'>>$sbatchfile;
    echo -e '#SBATCH --job-name lftp_transfer_fastq_'$tolid'_'$b >> $sbatchfile; 
    echo -e '#SBATCH --output=lftp_transfer_fastq_'$tolid'_'$b'_%j.out' >> $sbatchfile;
    echo -e '#SBATCH --error=lftp_transfer_fastq_'$tolid'_'$b'_%j.err'>> $sbatchfile;
    echo -e '#SBATCH --time=12:00:00'>> $sbatchfile;
    echo -e '#SBATCH --partition research'>> $sbatchfile;
    echo -e '#SBATCH --mem 5G'>> $sbatchfile;

    echo "";>>$sbatchfile;

    echo "#upload commands:">>$sbatchfile;

    echo "";>>$sbatchfile;

    echo -e 'module load lftp; ' >> $sbatchfile;

    echo "";>>$sbatchfile;

    echo lftp -c '"'open -u Webin-1543,BGWSEi8y webin.ebi.ac.uk";" put -c -O . $file'"' >>$sbatchfile; #-c option allows continue, reput. It requires permission to overwrite remote files
    
    chmod ug+x $sbatchfile;

    echo submitting $sbatchfile; 
    
    sbatch $sbatchfile;

  done


