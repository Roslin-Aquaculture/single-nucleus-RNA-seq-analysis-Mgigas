#!/bin/bash

# Logfile to store actions
logfile="logfile.txt"

# Clear previous log entries
> $logfile

# Create the final output directory 'star_outputs_collated' in the current working directory
final_dir="./star_outputs_collated"
echo "Creating final output directory: $final_dir"
echo "Creating final output directory: $final_dir" >> $logfile
mkdir -p "$final_dir"

# Check if the final directory was created successfully
if [ $? -eq 0 ]; then
    echo "Directory $final_dir created successfully."
    echo "Directory $final_dir created successfully." >> $logfile
else
    echo "Failed to create directory $final_dir."
    echo "Failed to create directory $final_dir." >> $logfile
    exit 1
fi

# Loop through the directories a1_resultsSolo.out to a8_resultsSolo.out
for dir in a*_resultsSolo.out; do
    # Navigate to the GeneFull/raw subdirectory
    sub_dir="$dir/GeneFull/raw"

    # Extract the prefix (e.g., a1, a2, a3) from the directory name
    prefix=$(echo $dir | cut -d'_' -f1)

    # Create the new subdirectory inside the final output directory with the extracted prefix
    new_dir="$final_dir/${prefix}_subdir"

    echo "Creating directory: $new_dir"
    echo "Creating directory: $new_dir" >> $logfile
    mkdir -p "$new_dir"

    # Check if the directory was created successfully
    if [ $? -eq 0 ]; then
        echo "Directory $new_dir created successfully."
        echo "Directory $new_dir created successfully." >> $logfile
    else
        echo "Failed to create directory $new_dir."
        echo "Failed to create directory $new_dir." >> $logfile
        continue
    fi

    # Files to copy
    files_to_copy=("barcodes.tsv" "features.tsv")

    # Copy each file (barcodes.tsv, features.tsv)
    for file in "${files_to_copy[@]}"; do
        src_file="$sub_dir/$file"
        dest_file="$new_dir/$file"

        echo "Copying $src_file to $dest_file"
        echo "Copying $src_file to $dest_file" >> $logfile
        cp "$src_file" "$dest_file"

        # Check if the file copy was successful
        if [ $? -eq 0 ]; then
            echo "Copied $file successfully."
            echo "Copied $file successfully." >> $logfile
        else
            echo "Failed to copy $file."
            echo "Failed to copy $file." >> $logfile
        fi
    done

    # Copy and rename UniqueAndMult-Rescue.mtx to matrix.mtx
    src_matrix="$sub_dir/UniqueAndMult-Rescue.mtx"
    dest_matrix="$new_dir/matrix.mtx"

    echo "Copying and renaming $src_matrix to $dest_matrix"
    echo "Copying and renaming $src_matrix to $dest_matrix" >> $logfile
    cp "$src_matrix" "$dest_matrix"

    # Check if the matrix file copy was successful
    if [ $? -eq 0 ]; then
        echo "Copied and renamed UniqueAndMult-Rescue.mtx to matrix.mtx successfully."
        echo "Copied and renamed UniqueAndMult-Rescue.mtx to matrix.mtx successfully." >> $logfile
    else
        echo "Failed to copy and rename UniqueAndMult-Rescue.mtx."
        echo "Failed to copy and rename UniqueAndMult-Rescue.mtx." >> $logfile
    fi

    echo "--------------------------------" >> $logfile
done

echo "Script execution complete. Check $logfile for details."
