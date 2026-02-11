#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <PROJECT_ID> [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --info                Show available data without downloading"
    echo "  --urls                Output download URLs only (properly escaped)"
    echo "  --strategy <list>     Download only specified library strategies (comma-separated)"
    echo ""
    echo "Examples:"
    echo "  $0 PRJEB91173 --info"
    echo "  $0 PRJEB91173 --urls"
    echo "  $0 PRJEB91173 --strategy Hi-C --urls"
    echo "  $0 PRJEB91173 --strategy RNA-Seq"
    echo "  $0 PRJEB91173 --strategy 'WGS,RNA-Seq'"
    echo "  $0 PRJEB91173 --strategy 'Hi-C' --info"
    echo ""
    exit 1
}

# Function to URL encode a string
urlencode() {
    local string="$1"
    local encoded=""
    local pos c o
    
    for (( pos=0 ; pos<${#string} ; pos++ )); do
        c=${string:$pos:1}
        case "$c" in
            [-_.~a-zA-Z0-9/] ) o="$c" ;;
            * ) printf -v o '%%%02X' "'$c"
        esac
        encoded+="$o"
    done
    echo "$encoded"
}

# Check if project ID provided
if [ -z "$1" ]; then
    usage
fi

PROJECT_ID="$1"
shift

# Parse options
INFO_ONLY=false
URLS_ONLY=false
STRATEGY_FILTER=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --info)
            INFO_ONLY=true
            shift
            ;;
        --urls)
            URLS_ONLY=true
            shift
            ;;
        --strategy)
            STRATEGY_FILTER="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT_ID}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_platform,library_strategy,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,submitted_md5,submitted_ftp,submitted_format,sra_ftp,bam_ftp,bam_md5&format=json&download=true&limit=0"

echo "Fetching metadata for project: $PROJECT_ID" >&2

# Download metadata
METADATA_FILE="${PROJECT_ID}_metadata.json"
wget -q -O "$METADATA_FILE" "$API_URL"

if [ ! -s "$METADATA_FILE" ]; then
    echo "Error: Failed to download metadata or project not found" >&2
    exit 1
fi

echo "Metadata downloaded successfully" >&2
echo "" >&2

# Build jq filter for strategy selection
if [ -n "$STRATEGY_FILTER" ]; then
    # Convert comma-separated list to jq array syntax
    IFS=',' read -ra STRATEGIES <<< "$STRATEGY_FILTER"
    JQ_FILTER='['
    for i in "${!STRATEGIES[@]}"; do
        strategy=$(echo "${STRATEGIES[$i]}" | xargs) # trim whitespace
        if [ $i -gt 0 ]; then
            JQ_FILTER+=','
        fi
        JQ_FILTER+="\"$strategy\""
    done
    JQ_FILTER+=']'
    
    JQ_SELECT="select(.library_strategy as \$s | $JQ_FILTER | index(\$s))"
else
    JQ_SELECT="."
fi

# If URLs only, output just the URLs and exit
if [ "$URLS_ONLY" = true ]; then
    jq -r ".[] | $JQ_SELECT | (.fastq_ftp // \"\") + \"|\" + (.submitted_ftp // \"\")" "$METADATA_FILE" | while IFS='|' read -r fastq_ftp submitted_ftp; do
        # Determine which files to use (prefer fastq, fallback to submitted)
        if [ -n "$fastq_ftp" ]; then
            FTP_URLS="$fastq_ftp"
        elif [ -n "$submitted_ftp" ]; then
            FTP_URLS="$submitted_ftp"
        else
            continue
        fi
        
        # Split multiple files (semicolon separated)
        IFS=';' read -ra URLS <<< "$FTP_URLS"
        
        for url in "${URLS[@]}"; do
            # Split path and filename
            dir=$(dirname "$url")
            file=$(basename "$url")
            
            # URL encode the filename
            encoded_file=$(urlencode "$file")
            
            # Output the properly encoded URL
            echo "ftp://${dir}/${encoded_file}"
        done
    done
    exit 0
fi

# If info only, display summary and exit
if [ "$INFO_ONLY" = true ]; then
    echo "=========================================="
    echo "PROJECT SUMMARY: $PROJECT_ID"
    echo "=========================================="
    echo ""
    
    # Count by platform
    echo "INSTRUMENT PLATFORMS:"
    jq -r '.[].instrument_platform' "$METADATA_FILE" | sort | uniq -c | sort -rn
    echo ""
    
    # Count by library strategy
    echo "LIBRARY STRATEGIES:"
    jq -r '.[].library_strategy' "$METADATA_FILE" | sort | uniq -c | sort -rn
    echo ""
    
    # Detailed breakdown
    echo "DETAILED BREAKDOWN:"
    echo "Platform | Strategy |   Runs  | Files Available"
    echo "---------|----------|---------|----------------"
    jq -r 'group_by(.instrument_platform + "|" + .library_strategy) | .[] | 
           .[0].instrument_platform + " | " + 
           .[0].library_strategy + " | " + 
           (length | tostring) + " | " + 
           (if .[0].fastq_ftp != null and .[0].fastq_ftp != "" then "FASTQ" 
            elif .[0].submitted_ftp != null and .[0].submitted_ftp != "" then "Submitted" 
            else "None" end)' "$METADATA_FILE" | column -t -s '|'
    echo ""
    
    # Total files
    TOTAL_RUNS=$(jq 'length' "$METADATA_FILE")
    FASTQ_AVAILABLE=$(jq '[.[] | select(.fastq_ftp != null and .fastq_ftp != "")] | length' "$METADATA_FILE")
    SUBMITTED_AVAILABLE=$(jq '[.[] | select(.submitted_ftp != null and .submitted_ftp != "")] | length' "$METADATA_FILE")
    
    echo "TOTALS:"
    echo "  Total runs: $TOTAL_RUNS"
    echo "  With FASTQ: $FASTQ_AVAILABLE"
    echo "  With submitted files: $SUBMITTED_AVAILABLE"
    echo ""
    echo "=========================================="
    
    if [ -n "$STRATEGY_FILTER" ]; then
        echo ""
        echo "FILTERED BY STRATEGY: $STRATEGY_FILTER"
        echo "=========================================="
        IFS=',' read -ra STRATEGIES <<< "$STRATEGY_FILTER"
        for strategy in "${STRATEGIES[@]}"; do
            strategy=$(echo "$strategy" | xargs) # trim whitespace
            COUNT=$(jq -r --arg strat "$strategy" '[.[] | select(.library_strategy == $strat)] | length' "$METADATA_FILE")
            echo "  $strategy: $COUNT runs"
        done
    fi
    
    exit 0
fi

echo "Filtering for library strategies: $STRATEGY_FILTER" >&2
echo "" >&2

# Parse JSON and download files
jq -r ".[] | $JQ_SELECT | .instrument_platform + \"|\" + .sample_accession + \"|\" + .library_strategy + \"|\" + .run_accession + \"|\" + (.fastq_ftp // \"\") + \"|\" + (.fastq_md5 // \"\") + \"|\" + (.submitted_ftp // \"\") + \"|\" + (.submitted_md5 // \"\")" "$METADATA_FILE" | while IFS='|' read -r platform sample_acc library_strat run_acc fastq_ftp fastq_md5 submitted_ftp submitted_md5; do
    
    # Determine which files to download (prefer fastq, fallback to submitted)
    if [ -n "$fastq_ftp" ] && [ -n "$fastq_md5" ]; then
        FTP_URLS="$fastq_ftp"
        MD5_SUMS="$fastq_md5"
        FILE_TYPE="fastq"
        echo ""
        echo "=========================================="
        echo "Using generated FASTQ files"
    elif [ -n "$submitted_ftp" ] && [ -n "$submitted_md5" ]; then
        FTP_URLS="$submitted_ftp"
        MD5_SUMS="$submitted_md5"
        FILE_TYPE="submitted"
        echo ""
        echo "=========================================="
        echo "⚠ No generated FASTQ available, using submitted files"
    else
        echo ""
        echo "=========================================="
        echo "✗ No files available for run $run_acc"
        echo "=========================================="
        continue
    fi
    
    # Split multiple files (semicolon separated)
    IFS=';' read -ra URLS <<< "$FTP_URLS"
    IFS=';' read -ra MD5S <<< "$MD5_SUMS"
    
    # Download each file
    for i in "${!URLS[@]}"; do
        # Split path and filename for URL encoding
        dir=$(dirname "${URLS[$i]}")
        file=$(basename "${URLS[$i]}")
        
        # URL encode the filename
        encoded_file=$(urlencode "$file")
        
        URL="ftp://${dir}/${encoded_file}"
        EXPECTED_MD5="${MD5S[$i]}"
        ORIGINAL_FILE="$file"  # Use unencoded filename for local file
        
        # Extract read number and extension
        if [[ $ORIGINAL_FILE =~ _([12])\.fastq\.gz$ ]]; then
            READ_NUM="${BASH_REMATCH[1]}"
            EXT=".fastq.gz"
        elif [[ $ORIGINAL_FILE =~ \.([^.]+)$ ]]; then
            # For submitted files, preserve original extension
            READ_NUM=""
            EXT=".${BASH_REMATCH[1]}"
            # Handle compound extensions like .fastq.gz
            if [[ $ORIGINAL_FILE =~ \.([^.]+\.[^.]+)$ ]]; then
                EXT=".${BASH_REMATCH[1]}"
            fi
        else
            READ_NUM=""
            EXT=""
        fi
        
        # Create new filename: platform.library_strategy.sample_accession.run_accession[.READ_NUM].ext
        if [ -n "$READ_NUM" ]; then
            NEW_FILE="${platform}.${library_strat}.${sample_acc}.${run_acc}.${READ_NUM}${EXT}"
        else
            NEW_FILE="${platform}.${library_strat}.${sample_acc}.${run_acc}${EXT}"
        fi
        
        echo "Processing: $NEW_FILE"
        echo "Platform: $platform"
        echo "Sample: $sample_acc"
        echo "Strategy: $library_strat"
        echo "Run: $run_acc"
        echo "Type: $FILE_TYPE"
        
        # Check if renamed file already exists with correct MD5
        if [ -f "$NEW_FILE" ]; then
            EXISTING_MD5=$(md5sum "$NEW_FILE" | awk '{print $1}')
            if [ "$EXISTING_MD5" == "$EXPECTED_MD5" ]; then
                echo "✓ File already exists with correct MD5, skipping"
                echo "=========================================="
                continue
            else
                echo "⚠ File exists but has incorrect MD5, will re-download"
                rm "$NEW_FILE"
            fi
        fi
        
        # Check if original file exists and handle appropriately
        if [ -f "$ORIGINAL_FILE" ]; then
            EXISTING_MD5=$(md5sum "$ORIGINAL_FILE" | awk '{print $1}')
            if [ "$EXISTING_MD5" == "$EXPECTED_MD5" ]; then
                echo "✓ Original file exists with correct MD5, renaming"
                mv "$ORIGINAL_FILE" "$NEW_FILE"
                echo "✓ Renamed to $NEW_FILE"
                echo "=========================================="
                continue
            else
                echo "⚠ Existing file has incorrect MD5, removing and re-downloading"
                rm "$ORIGINAL_FILE"
            fi
        fi
        
        echo "Downloading: $ORIGINAL_FILE"
        echo "URL: $URL"
        echo "Expected MD5: $EXPECTED_MD5"
        echo "=========================================="
        
        # Download file with improved settings
        wget -c --retry-connrefused \
             --waitretry=3 \
             --read-timeout=60 \
             --timeout=30 \
             --tries=0 \
             --passive-ftp \
             "$URL"
        
        # Verify checksum and rename
        if [ -f "$ORIGINAL_FILE" ]; then
            ACTUAL_MD5=$(md5sum "$ORIGINAL_FILE" | awk '{print $1}')
            echo "Actual MD5:   $ACTUAL_MD5"
            
            if [ "$ACTUAL_MD5" == "$EXPECTED_MD5" ]; then
                echo "✓ Checksum verified for $ORIGINAL_FILE"
                mv "$ORIGINAL_FILE" "$NEW_FILE"
                echo "✓ Renamed to $NEW_FILE"
            else
                echo "✗ ERROR: Checksum mismatch for $ORIGINAL_FILE"
                echo "  Renaming to ${ORIGINAL_FILE}.corrupt"
                mv "$ORIGINAL_FILE" "${ORIGINAL_FILE}.corrupt"
            fi
        else
            echo "✗ ERROR: File $ORIGINAL_FILE not found after download"
        fi
    done
done

echo ""
echo "=========================================="
echo "Download complete!"
echo "Metadata saved to: $METADATA_FILE"
echo "=========================================="