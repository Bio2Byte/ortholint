#!/bin/bash
GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium
input_file=$GHOME/clean_achromatium/structures/seqs_to_model.fasta
output_file=$GHOME/clean_achromatium/structures/remaining_seqs_to_model.fasta

# Create a temporary file to store identifiers
temp_file="modelled_structures.txt"

# Extract all identifiers from filenames in the folder
for file in *.pdb; do
    basename "$file" .pdb >> "$temp_file"
done

# Process the input file
awk -v ids_file="$temp_file" '
    BEGIN {
        # Load identifiers into an array
        while ((getline line < ids_file) > 0) {
            ids[line] = 1;
        }
        close(ids_file);
    }
    /^>/ {
        # Check if the current header is in the identifiers
        keep = 1;
        header = substr($0, 2);
        if (header in ids) {
            keep = 0;  # Skip this block
        }
    }
    keep { print }
' "$input_file" > "$output_file"

# Clean up temporary file
rm "$temp_file"

echo "Processing complete. Filtered output written to $output_file."