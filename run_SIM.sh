#!/bin/bash

# Exit on any error
set -e

# Define variables
EXPNO=$1
INPUT_FILE=$2
TARGET_DIR="Experiments/$EXPNO"
SOURCE_FILE="namelist"   # Change this to your source file name
EXECUTABLE="zoupa"
POST_FILE="./output_postfiles/output_$EXPNO.log"


# Step 1: Check the arguments
if [ -z "$EXPNO" ] || [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 <exp_number> <input_file>" >&2
    exit 1
fi
if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file not found: $INPUT_FILE" >&2
    exit 1
fi

# Step 2: Set the exp number in the input file
# Replaces the first value on the "exp version" line, keeping the rest of the line.
tmp=$(mktemp)
awk -v expno="$EXPNO" '
    /exp version/ { sub(/^[ \t]*[^ \t]+/, expno); found = 1 }
    { print }
    END { if (!found) exit 1 }
' "$INPUT_FILE" > "$tmp" || {
    rm -f "$tmp"
    echo "Error: no 'exp version' line found in $INPUT_FILE" >&2
    exit 1
}
cat "$tmp" > "$INPUT_FILE"   # overwrite contents, keep file permissions
rm -f "$tmp"
echo "Set exp version to $EXPNO in $INPUT_FILE"


# Step 2: Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Step 3: Copy the source file to the directory
cp "$SOURCE_FILE" "$TARGET_DIR/"

#Step 4: logging the experiment 
sh log_experiment.sh $SOURCE_FILE $INPUT_FILE ''

# Step 4: Run the model
nohup ./"$EXECUTABLE" < $INPUT_FILE > "$POST_FILE" &

