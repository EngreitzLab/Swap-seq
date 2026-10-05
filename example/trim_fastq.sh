#!/bin/bash

# Create the output directory if it doesn't already exist
mkdir -p trimed_fastq

# Iterate over all .fastq.gz files in the current directory
for file in fastq/*.fastq.gz; do
  # Check to ensure the file does not contain I1, I2, or Undetermined
  if [[ $file != *I1* && $file != *I2* && $file != *Undetermined* ]]; then
    # Define the output file name, keeping the original name
    output_file=$(basename "$file")
    output_path="trimed_fastq/${output_file}"
    
    # Process and trim the file
    zcat "$file" | \
    awk 'NR % 4 == 2 || NR % 4 == 0 { print substr($1, 1, 71) } NR % 4 == 1 || NR % 4 == 3' | \
    gzip > "$output_path"
    
    echo "Processed $file and saved to $output_path"
  else
    echo "Skipping $file (contains I1, I2, or Undetermined)"
  fi
done