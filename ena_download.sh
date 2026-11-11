#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <PROJECT_ID> [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --info                Show available data without downloading"
    echo "  --urls                Output download URLs only (properly escaped)"
    echo "  --strategy <list>     Download only specified library strategies (comma-separated)"
    echo "  -j, --jobs <N>        Download N files in parallel (default: 4)"
    echo ""
    echo "Examples:"
    echo "  $0 PRJEB91173 --info"
    echo "  $0 PRJEB91173 --urls"
    echo "  $0 PRJEB91173 --strategy Hi-C --urls"
    echo "  $0 PRJEB91173 --strategy RNA-Seq"
    echo "  $0 PRJEB91173 --strategy 'WGS,RNA-Seq'"
    echo "  $0 PRJEB91173 --strategy 'Hi-C' --info"
    echo "  $0 PRJEB91173 --jobs 8"
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

# Download a file, trying HTTPS first and falling back to FTP.
# Both protocols serve the same path on ftp.sra.ebi.ac.uk; which one is
# faster/reliable varies by site and over time (FTP data connections stall
# from some networks). -c resumes any partial file, so switching protocol
# or re-running the script never throws away bytes already downloaded.
# Usage: download_file <host/path/encoded_name> <local_file>
DOWNLOAD_ROUNDS=3     # rounds of HTTPS-then-FTP
TRIES_PER_PROTOCOL=5  # wget attempts per protocol per round
download_file() {
    local path="$1" out="$2" round proto
    for (( round=1 ; round<=DOWNLOAD_ROUNDS ; round++ )); do
        for proto in https ftp; do
            echo "Trying ${proto}://${path} (round ${round}/${DOWNLOAD_ROUNDS})"
            if wget -c --retry-connrefused \
                    --waitretry=10 \
                    --read-timeout=60 \
                    --timeout=30 \
                    --tries="$TRIES_PER_PROTOCOL" \
                    --progress="$WGET_PROGRESS" \
                    -O "$out" \
                    "${proto}://${path}"; then
                return 0
            fi
            echo "⚠ ${proto} download failed or stalled, switching protocol"
        done
    done
    return 1
}

# Download one file, verify its MD5 and rename it.
# Usage: process_file <platform> <library_strategy> <sample_acc> <run_acc>
#                     <file_type> <host/path/name> <expected_md5>
# Returns 0 if the file is present and verified at the end, 1 otherwise.
process_file() {
    local platform="$1" library_strat="$2" sample_acc="$3" run_acc="$4"
    local file_type="$5" url="$6" expected_md5="$7"
    local dir file url_path original_file read_num ext new_file actual_md5

    # Split path and filename; URL-encode only the filename
    dir=$(dirname "$url")
    file=$(basename "$url")
    url_path="${dir}/$(urlencode "$file")"  # host/path without protocol
    original_file="$file"  # Use unencoded filename for local file

    # Extract read number and extension
    if [[ $original_file =~ _([12])\.fastq\.gz$ ]]; then
        read_num="${BASH_REMATCH[1]}"
        ext=".fastq.gz"
    elif [[ $original_file =~ \.([^.]+)$ ]]; then
        # For submitted files, preserve original extension
        read_num=""
        ext=".${BASH_REMATCH[1]}"
        # Handle compound extensions like .fastq.gz
        if [[ $original_file =~ \.([^.]+\.[^.]+)$ ]]; then
            ext=".${BASH_REMATCH[1]}"
        fi
    else
        read_num=""
        ext=""
    fi

    # Create new filename: platform.library_strategy.sample_accession.run_accession[.READ_NUM].ext
    if [ -n "$read_num" ]; then
        new_file="${platform}.${library_strat}.${sample_acc}.${run_acc}.${read_num}${ext}"
    else
        new_file="${platform}.${library_strat}.${sample_acc}.${run_acc}${ext}"
    fi

    echo "Processing: $new_file (platform $platform, sample $sample_acc, strategy $library_strat, run $run_acc, type $file_type)"

    # Check if renamed file already exists with correct MD5
    if [ -f "$new_file" ]; then
        if [ "$(md5sum "$new_file" | awk '{print $1}')" == "$expected_md5" ]; then
            echo "✓ File already exists with correct MD5, skipping"
            return 0
        fi
        echo "⚠ File exists but has incorrect MD5, will re-download"
        rm "$new_file"
    fi

    # An existing original file is usually a partial download from an
    # interrupted run: keep it and let wget -c resume it. If it is
    # actually complete, wget finishes immediately; if it is complete
    # but corrupt, the MD5 check below catches it.
    if [ -f "$original_file" ]; then
        echo "Found existing $original_file ($(wc -c < "$original_file" | tr -d ' ') bytes), resuming"
    fi

    echo "Downloading: $original_file (expected MD5 $expected_md5)"

    if ! download_file "$url_path" "$original_file"; then
        echo "✗ ERROR: Download of $original_file failed over both HTTPS and FTP; partial file kept, re-run to resume"
        return 1
    fi

    # Verify checksum and rename
    if [ ! -f "$original_file" ]; then
        echo "✗ ERROR: File $original_file not found after download"
        return 1
    fi
    actual_md5=$(md5sum "$original_file" | awk '{print $1}')
    if [ "$actual_md5" == "$expected_md5" ]; then
        mv "$original_file" "$new_file"
        echo "✓ Checksum verified, renamed to $new_file"
        return 0
    fi
    echo "✗ ERROR: Checksum mismatch for $original_file (got $actual_md5), renaming to ${original_file}.corrupt"
    mv "$original_file" "${original_file}.corrupt"
    return 1
}

# Run process_file with every output line prefixed by the file name, so the
# logs of parallel downloads stay readable when they interleave.
process_file_logged() {
    local name status
    name=$(basename "$6")
    process_file "$@" 2>&1 | awk -v p="[$name] " '{ print p $0; fflush() }'
    status=${PIPESTATUS[0]}
    return "$status"
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
JOBS=4

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
        -j|--jobs)
            JOBS="$2"
            if ! [[ "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
                echo "Error: --jobs needs a positive integer, got '$JOBS'"
                usage
            fi
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

# Settings and functions the parallel workers (bash -c under xargs) need
if [ "$JOBS" -eq 1 ] && [ -t 1 ]; then
    WGET_PROGRESS="bar"
else
    # One line per 32 MB: readable in batch logs and when jobs interleave
    WGET_PROGRESS="dot:giga"
fi
export DOWNLOAD_ROUNDS TRIES_PER_PROTOCOL WGET_PROGRESS
export -f urlencode download_file process_file process_file_logged

echo "Downloading with $JOBS parallel job(s)" >&2

# List every file to download as 7 NUL-separated fields
# (platform, strategy, sample, run, type, url, md5), then hand them to
# xargs, which runs up to $JOBS process_file_logged workers at once.
jq -r ".[] | $JQ_SELECT | .instrument_platform + \"|\" + .sample_accession + \"|\" + .library_strategy + \"|\" + .run_accession + \"|\" + (.fastq_ftp // \"\") + \"|\" + (.fastq_md5 // \"\") + \"|\" + (.submitted_ftp // \"\") + \"|\" + (.submitted_md5 // \"\")" "$METADATA_FILE" | while IFS='|' read -r platform sample_acc library_strat run_acc fastq_ftp fastq_md5 submitted_ftp submitted_md5; do

    # Determine which files to download (prefer fastq, fallback to submitted)
    if [ -n "$fastq_ftp" ] && [ -n "$fastq_md5" ]; then
        FTP_URLS="$fastq_ftp"
        MD5_SUMS="$fastq_md5"
        FILE_TYPE="fastq"
    elif [ -n "$submitted_ftp" ] && [ -n "$submitted_md5" ]; then
        FTP_URLS="$submitted_ftp"
        MD5_SUMS="$submitted_md5"
        FILE_TYPE="submitted"
        echo "⚠ No generated FASTQ available for run $run_acc, using submitted files" >&2
    else
        echo "✗ No files available for run $run_acc" >&2
        continue
    fi

    # Split multiple files (semicolon separated)
    IFS=';' read -ra URLS <<< "$FTP_URLS"
    IFS=';' read -ra MD5S <<< "$MD5_SUMS"

    for i in "${!URLS[@]}"; do
        printf '%s\0' "$platform" "$library_strat" "$sample_acc" "$run_acc" \
            "$FILE_TYPE" "${URLS[$i]}" "${MD5S[$i]}"
    done
done | xargs -0 -n 7 -P "$JOBS" bash -c 'process_file_logged "$@"' _
XARGS_STATUS=$?

echo ""
echo "=========================================="
if [ "$XARGS_STATUS" -eq 0 ]; then
    echo "Download complete!"
else
    echo "✗ Some files failed to download or verify (search the log for ✗); re-run to resume"
fi
echo "Metadata saved to: $METADATA_FILE"
echo "=========================================="
[ "$XARGS_STATUS" -eq 0 ]
