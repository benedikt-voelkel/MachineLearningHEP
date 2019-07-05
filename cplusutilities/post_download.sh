#!/bin/bash

function packing()
{
    local output_dir="$1"
    local file_pack="$2"
    echo "===> Write pack to $output_dir"
    mkdir -p $output_dir
    hadd $output_dir/AnalysisResults.root $file_pack
}

function n_job_delay()
{
    local n_max_jobs="$1"
    local sleep_time="$2"
    while true
    do
        n_packing=$(jobs | grep "packing" | grep "Running" | wc -l )
        if (( $n_packing >= $n_max_jobs ))
        then
            sleep $sleep_time
        else
            break
        fi
    done
}

# Some colors
STREAM_START_RED="\033[0;31m"
STREAM_START_GREEN="\033[0;32m"
STREAM_START_YELLOW="\033[0;33m"
STREAM_END_COLOR="\033[0m"

# Top directory where pointing down to the train number
# In that directory the further structure then should be
# .../unmerged/child_N/nnnn/AnalysisResults.root
INPUT_PATH="$1"

# We start from here
CURR_DIR=$(pwd)

# Normal upper bound for merged files
TARGET_PACK_SIZE="1000000"

# Accepted size if one input file is already bigger
MAX_ACCEPTED_INPUT_SIZE="5000000"
# Number of those
MAX_ACCEPTED_BIG_INPUT="25"

# Number of packing jobs
N_PACKING_JOBS=28

MERGED_DIR="merged"

echo "#####"
echo "#####"
echo "MERGING GRID DATA up to target file size of $TARGET_PACK_SIZE kB"
echo "#####"
echo "#####"
echo

[[ ! -d $INPUT_PATH ]] && { echo "$INPUT_PATH is no directory"; return 1; }

# Make it an absolute path...
INPUT_PATH=$(realpath $INPUT_PATH)

# ... and go there
cd $INPUT_PATH

# check for unmerged directory
[[ ! -d "./unmerged" ]] && { echo -e "${STREAM_START_RED}ERROR${STREAM_END_COLOR}: Cannot find unmgered directory"; exit 1; }

# keep the old data savely and produce the merged data as well. This assumes
# there is nothing but the ROOT data from the grid
unmerged_size="$(du -s | awk '{print $1}' )"
free_space="$(df . | grep "/dev" | awk '{print $4}')"

if (( $free_space < $unmerged_size ))
then
    echo -e "${STREAM_START_RED}ERROR${STREAM_END_COLOR}: Not enough disk space left"
fi


# Fail if "merged" directory exists already
[[ -d "$MERGED_DIR" ]] && { echo -e "${STREAM_START_RED}ERROR${STREAM_END_COLOR}: Seems that the merge directory already exists"; exit 1; }


# If we are here, things sould be fine
#mkdir "merged"

# Merge per child so find out childs we have
childs=$(find unmerged -maxdepth 1 -type d -name "child_*" | sort -u)

echo "===> Found childs"
echo "$childs"
echo

# Make the merged dir
mkdir $MERGED_DIR
for c in $childs
do
    c_stripped=${c##unmerged/}
    echo "===> Process $c_stripped"
    # For each child_i there will be a pack_i
    MERGED_CHILD_DIR=$MERGED_DIR/$c_stripped
    mkdir -p $MERGED_CHILD_DIR
    root_files_childs=$(find $c -maxdepth 2 -type f -name "AnalysisResults.root")
    if [[ "$root_files_childs" == "" ]]
    then
        echo -e "${STREAM_START_YELLOW}WARNING${STREAM_END_COLOR}: No ROOT files found in $c"
        continue
    fi
    n_packs="0"
    file_pack=""
    current_size="0"
    n_big_files="0"
    for rfc in $root_files_childs
    do
        next_size=$(du -s $rfc | awk '{print $1}')
        if (( $next_size > $MAX_ACCEPTED_INPUT_SIZE ))
        then
            echo -e "${STREAM_START_RED}ERROR${STREAM_END_COLOR}: File $rfc is bigger than $MAX_ACCEPTED_INPUT_SIZE kB. Not accepted..."
            exit 1
        fi

        if [[ "$file_pack" == "" ]]
        then
            file_pack+="$rfc "
            current_size=$(( $current_size + $next_size ))
            if (( ( $current_size < $MAX_ACCEPTED_INPUT_SIZE ) && ( $current_size > $TARGET_PACK_SIZE ) ))
            then
                n_big_files=$(( $n_big_files + 1 ))
                if (( $n_big_files > $MAX_ACCEPTED_BIG_INPUT ))
                then
                    echo -e "${STREAM_START_RED}ERROR${STREAM_END_COLOR}: More than $MAX_ACCEPTED_BIG_INPUT big files are not accepted"
                    exit 1
                fi
            fi

            continue
        fi
        if (( ($current_size < $TARGET_PACK_SIZE)  && ( $TARGET_PACK_SIZE > ( $next_size + $current_size ) )))
        then
            file_pack+="$rfc "
            current_size=$(( $current_size + $next_size ))
        else
            output_dir="$MERGED_CHILD_DIR/pack_${n_packs}"
            n_job_delay $N_PACKING_JOBS 20
            packing $INPUT_PATH/$output_dir "$file_pack" &
            # Need to add that since it would be skipped otherwise
            file_pack="$rfc "
            current_size="$next_size"
            n_packs="$(( $n_packs + 1 ))"
        fi
    done

    # Handle the last pack
    output_dir="$MERGED_CHILD_DIR/pack_${n_packs}"
    n_job_delay $N_PACKING_JOBS 1
    packing $INPUT_PATH/$output_dir "$file_pack" &

done

n_job_delay 1 20 

echo "DONE"

