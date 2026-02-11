#!/bin/bash
# Authors: fcruz, talioto
# Date: 2025-09-08
# Brief description/notes:
#   creates basic directory structure to keep organized the curation folder
#

# Set BASE directory (default to "curation" if no argument provided)
export BASE=${1:-curation}

echo "Creating curation directory structure in: $BASE"

# Define directory structure
directories=(
    "1_pre/in"
    "1_pre/contam_detect"
    "1_pre/wga"
    "2_curated/1_rapid"
    "2_curated/2_decontaminate"
    "2_curated/3_claws"
    "2_curated/4_ear"
    "3_rev/1_rapid"
    "3_rev/2_decontaminate"
    "3_rev/3_claws"
    "3_rev/4_ear"
    "ena/hap1"
    "ena/hap2"
    "ena/primary"
)

# Create directories
for dir in "${directories[@]}"; do
    mkdir -p "$BASE/$dir"
    echo "Created: $BASE/$dir"
done

# Grant Permissions to Team
find "$BASE" -type d -exec chmod 1770 '{}' \;

echo "Directory structure created successfully!"
echo "All directories have been set with permissions 1770 (sticky bit + rwxrwx---)"