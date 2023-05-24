#!/bin/bash

set -e

# This script implements chunked upload protocol described here
# https://docs.nextcloud.com/server/17/developer_manual/client_apis/WebDAV/chunking.html

# URL finishing in /webdav/ hides the upload interface
##NEXTCLOUD_SERVER=https://b2drop.bsc.es/remote.php/webdav/
NEXTCLOUD_SERVER=https://b2drop.bsc.es/remote.php/dav/

# The destination relative path, ending in a slash
#NEXTCLOUD_PATH=INB_ELIXIR-ES_Hub/BIMCV-PadChest/
NEXTCLOUD_PATH="<yourpath>/"

# You are adviced to create an app password
# There you will get your username
# and the app password to put below
NEXTCLOUD_USER="<youruser>"
NEXTCLOUD_PASS="<yourpassword>"

blockSize=$((1024 * 1024))
#numBlockSize=1
numBlockSize=500
chunkSize=$((blockSize * numBlockSize))

parallelChunks=1
parallelUploads=1
semUploads="$0"

if [ -z "$DEBUG" ] ; then
	PRODFLAGS="-f"
else
	PRODFLAGS="-i"
fi

# This is needed for proper work in parallel downloads
checkChunk=
export checkChunk blockSize numBlockSize chunkSize NEXTCLOUD_SERVER NEXTCLOUD_USER NEXTCLOUD_PASS PRODFLAGS

uploadChunk() {
	local file="$1"
	local skip="$2"
	local tempRemoteFolder="$3"
	local a="$4"
	local b="$5"
	
	local padA="$(printf "%015d" "$a")"
	local padB="$(printf "%015d" "$b")"
		
	#echo curl ${PRODFLAGS} -k -X PUT -H "Content-Length: $((b - a + 1))" -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/"${padA}-${padB}" -T '-'
	if [ -z "$usetemp" ] ; then
		local contentLength=$((b - a + 1))
		dd if="$file" ibs="${blockSize}"c skip="${skip}" count="${numBlockSize}" status=none | curl ${PRODFLAGS} -k -X PUT -H "Content-Length: ${contentLength}" -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/"${padA}-${padB}" -T '-'
		# Check the first transfer worked
		if [ -z "$checkchunk" ] ; then
			checkchunk=1
			local uchunksize="$(curl ${PRODFLAGS} -s -k -X PROPFIND -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/"${padA}-${padB}" | grep -o '<d:getcontentlength>[0-9]*</d:getcontentlength>' | grep -o '[0-9]*')"
			if [ "$uchunksize" -ne "$contentLength" ] ; then
				echo "Avoiding chunked upload due server bug"
				usetemp=1
			fi
		fi
	fi
	
	# Experiment with wget
	#dd if="$file" ibs="${blockSize}"c skip="${skip}" count="${numBlockSize}" of=/tmp/chunktest.bin
	#wget -S --method PUT --body-file=/tmp/chunktest.bin --header="Content-Length: $((b - a + 1))" --user="${NEXTCLOUD_USER}" --password="${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/"${padA}-${padB}"
	
	if [ -n "$usetemp" ] ; then
		local tempchunk="$(mktemp -t NC-XXXXXXXX-UPLOAD.bin)"
		dd if="$file" ibs="${blockSize}"c skip="${skip}" count="${numBlockSize}" of="${tempchunk}" status=none
		curl ${PRODFLAGS} -k -X PUT -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/"${padA}-${padB}" -T "${tempchunk}"
		rm -f "${tempchunk}"
	fi
}
export -f uploadChunk

# When this is called, it is assumed remotePath already exists
uploadFile() {
	local file="$1"
	local remotePath="$2"
	local parallelChunks="$3"
	
	local relFilename="$(basename "$file")"
	local fileSize="$(stat -L -c %s "$file")"
	local fileSize1="$((fileSize - 1))"
	local skip=0
	
	# First, create temp folder for the chunks and assure destination folder exists
	local remoteFolder="${NEXTCLOUD_SERVER}files/${NEXTCLOUD_USER}/${remotePath}"
	local remoteFilename="${remoteFolder}${relFilename}"
	local uploadUUID=$(</proc/sys/kernel/random/uuid)
	local tempRemoteFolder="${NEXTCLOUD_SERVER}uploads/${NEXTCLOUD_USER}/${relFilename}-${uploadUUID}"
	
	echo "Uploading $file into $remoteFilename"
	echo
	echo INFO: If you interrupt the upload, or the upload fails, use next command line to free resources:
	echo curl -X DELETE -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/
	echo
	
	curl ${PRODFLAGS} -k -X MKCOL -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"
	
	local checkchunk=
	# Skip checking uploaded chunk on 0 size files
	if [ "$fileSize" -eq 0 ] ; then
		checkchunk=0
	fi
	local usetemp=
	for offset in $(seq 0 "${chunkSize}" "${fileSize}") ; do
		echo "$offset $fileSize"
		local a="${offset}"
		local b="$((offset + chunkSize -1))"
		if [ "$b" -gt "$fileSize1" ] ; then
			b="$fileSize1"
		fi
		
		if [ "$parallelChunks" = 1 -o "$a" = 0 ] ; then
			# First chunk is uploaded to check bugs
			uploadChunk "$file" "$skip" "$tempRemoteFolder" "$a" "$b"
		else
			sem --jobs "$parallelChunks" --id "$uploadUUID" uploadChunk "$file" "$skip" "$tempRemoteFolder" "$a" "$b"
		fi
		
		skip=$((skip + numBlockSize))
	done
	echo
	echo INFO: If the upload was interrupted, use next command line to free resources:
	echo curl -X DELETE -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${tempRemoteFolder}"/
	echo
	
	if [ "$parallelChunks" -gt 1 ] ; then
		sem --wait --id "$uploadUUID"
	fi
	
	curl ${PRODFLAGS} -k -X MOVE -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" --header "Destination:${remoteFilename}" --header "OC-Total-Length:${fileSize}" "${tempRemoteFolder}"/.file
	echo "* Uploaded $file!"
	echo
}

# This is needed for sem to work
export -f uploadFile

recursiveMkCol() {
	local remotePath="$1"
	
	local -a components
	IFS='/' read -r -a components <<< "${remotePath}"
	
	local remoteFolder="${NEXTCLOUD_SERVER}files/${NEXTCLOUD_USER}/"
	local comp
	for comp in "${components[@]}" ; do
		remoteFolder="${remoteFolder}${comp}/"
		curl -s -S -f -k -X PROPFIND -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${remoteFolder}" >& /dev/null || ( echo "Creating ${remoteFolder}" && curl -s -S ${PRODFLAGS} -k -X MKCOL -u "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" "${remoteFolder}" )
	done
}

uploadEntries() {
	local parallelUploads="$1"
	shift
	local parallelChunks="$1"
	shift
	local remotePath="$1"
	shift
	
	recursiveMkCol "${remotePath}"
	local entry
	for entry in "$@" ; do
		if [ -f "$entry" ] ; then
			if [ "$parallelUploads" = 1 ] ; then
				uploadFile "$entry" "${remotePath}" "$parallelChunks"
				
				echo
				echo "Uploaded $entry"
			else
				sem --jobs "$parallelUploads" --id "${semUploads}" uploadFile "$entry" "${remotePath}" "$parallelChunks"
				
				echo
				echo "Queued $entry"
			fi
		elif [ -d "$entry" ]; then
			echo "Processing directory $entry"
			echo
			
			local newRemotePath="${remotePath}$(basename "$entry")/"
			uploadEntries "${parallelUploads}" "${parallelChunks}" "${newRemotePath}" "${entry}"/*
			
			echo
			echo "Finished processing directory $entry"
		else
			echo "Path $entry does not locally exist"
		fi
	done
}

if [ $# -gt 1 ]; then
	configFile="$1"
	shift
	
	if [ -f "$configFile" ] ; then
		source "$configFile"
	else
		echo "ERROR: Configuration file $configFile not found!" 1>&2
		exit 1
	fi
	
	# Chunk size is computed after the config file has been read
	chunkSize=$((blockSize * numBlockSize))
	
	# Is available parallel uploads feature?
	if type -f sem >& /dev/null ; then
		# Number of parallel upload jobs is controlled
		maxCPUs="$(nproc)"
		if [ "$parallelChunks" -lt 1 ] ; then
			parallelChunks="$maxCPUs"
		elif [ "$parallelChunks" -gt "$maxCPUs" ] ; then
			parallelChunks="$maxCPUs"
		fi
		if [ "$parallelUploads" -lt 1 ] ; then
			parallelUploads="$maxCPUs"
		elif [ "$parallelUploads" -gt "$maxCPUs" ] ; then
			parallelUploads="$maxCPUs"
		fi
		echo "INFO: ${parallelUploads} parallel uploads, ${parallelChunks} parallel chunk uploads will be tried"
		semUploads="$(basename "$0")-$(</proc/sys/kernel/random/uuid)"
	else
		echo "WARNING: Parallel chunk uploads disabled, as `sem` is not available" 1>&2
		# No sem, no parallel uploads
		parallelUploads=1
		parallelChunks=1
	fi
	
	case "${NEXTCLOUD_PATH}" in
		*/)
			# Nothing to be done
			true
			;;
		*)
			# The destination relative path, ending in a slash
			NEXTCLOUD_PATH="${NEXTCLOUD_PATH}/"
	esac

	# The other files are the ones to upload
	declare -a inputs=()
	for input in "$@" ; do
		# Trying to avoid corner cases like '.'
		case "$input" in
			.*)
				input="$(realpath "$input")"
				;;
		esac
		inputs+=("$input")
	done
	uploadEntries "${parallelUploads}" "${parallelChunks}" "${NEXTCLOUD_PATH}" "${inputs[@]}"
	if [ "$parallelUploads" -gt 1 ] ; then
		echo "Waiting for queued and ongoing upload jobs to finish"
		sem --wait --id "${semUploads}"
		echo
		echo "Finished uploads"
	fi
else
	cat <<EOF
Usage: bash $0 {setup_file} {file_or_directory}+
EOF
fi
#curl --progress-bar --verbose --user "${NEXTCLOUD_USER}:${NEXTCLOUD_PASS}" -T "${file}" "${NEXTCLOUD_SERVER}files/${NEXTCLOUD_USER}/${NEXTCLOUD_PATH}"
