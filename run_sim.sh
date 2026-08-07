#!/bin/bash

# Exit on any error
set -e

# Define variables
EXPNO="14"
TARGET_DIR="Experiments/$EXPNO"
SOURCE_FILE="namelistSIM_panarctic"   # Change this to your source file name
EXECUTABLE="zoupa"
POST_FILE="./output_post_files/output_panarctic_$EXPNO"


# Step 2: Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Step 3: Copy the source file to the directory
cp "$SOURCE_FILE" "$TARGET_DIR/"

# Step 4: Run the model
nohup ./"$EXECUTABLE" < input_30day_12hr > "$POST_FILE" &

