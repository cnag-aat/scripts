#! /bin/bash


# User: fcruz
# Date: 2025-03-05

# Brief description/notes:
#   creates basic directory structure to keep organized the curation folder
#
#

# Create subdirectories
mkdir -p curation/1_yahs_pretext/

mkdir -p curation/2_checks/in
mkdir -p curation/2_checks/blobtoolkit/btk_out
mkdir -p curation/2_checks/WGA

mkdir -p curation/3_decontaminate/in
mkdir -p curation/3_decontaminate/out

mkdir -p curation/4_rapid/in
mkdir -p curation/4_rapid/out

mkdir -p curation/5_evaluations/in
mkdir -p curation/5_evaluations/blobtoolkit/btk_out
mkdir -p curation/5_evaluations/claws/out
mkdir -p curation/5_evaluations/ear/out

mkdir -p curation/6_review/in
mkdir -p curation/6_review/rapid/out
mkdir -p curation/6_review/claws/
mkdir -p curation/5_evaluations/ear/

mkdir -p curation/7_ENA/

# Grant Permissions to Team
find curation -type d -exec chmod 1770 '{}' \;
