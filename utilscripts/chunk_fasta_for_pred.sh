#!/bin/bash
GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium
#input_file=$GHOME/clean_achromatium/structures/remaining_seqs_to_model.fasta
input_file=$GHOME/clean_achromatium/structures/seqs_to_model_split_20.fasta
# Define the intervals
intervals=(10 20 30 40 50 60 70 80 90 95 100 105 110 115 120 125 130 135)

# Initialize counters
current_interval_index=0
entry_count=0
#output_file="$GHOME/clean_achromatium/structures/seqs_to_model_split_0_${intervals[$current_interval_index]}.fasta"
output_dir="${GHOME}/clean_achromatium/structures"

awk -v intervals="${intervals[*]}" -v output_dir="$output_dir" '
    BEGIN {
        split(intervals, limits);
        limit = limits[1];  # Initialize first interval limit
        file_idx = 1;       # File counter
        output_file = output_dir "/seqs_to_model_split_20_" file_idx ".fasta";
    }
    /^>/ { 
        # Increment entry count for each new sequence header
        entry_count++;
        if (entry_count > limit) {
            file_idx++;  # Move to next interval
            limit = limits[file_idx];  # Update limit
            output_file = output_dir "/seqs_to_model_split_20_" file_idx ".fasta";  # Update output file
        }
    }
    {
        # Append to the current file
        print > output_file;
    }
' "$input_file"

# Post-processing each split file
for ((i=1; i<=${#intervals[@]}; i++)); do
    output_file="${output_dir}/seqs_to_model_split_20_${i}.fasta"

    # Check if the file exists
    if [[ -f $output_file ]]; then
        echo "Processing file: $output_file"

        # Count headers
        sequence_count=$(grep -c ">" "$output_file")
        echo "Number of sequences: $sequence_count"

        # Calculate prediction time
        prediction_time=$(awk '
            !/^>/ && $0 != "" { 
                len = length($0); 
                length_sq = len ^ 2; 
                total += 0.00038 * length_sq + 0.0551 * len + 4.92
            } 
            END { print total }
        ' "$output_file")
        echo "Predicted time (seconds): $prediction_time"
        echo ""
    fi
done


#n=$1             
# Use awk to skip the first n entries
#awk -v n="$n" '
#    BEGIN { count = 0; skip = 1 }
#    /^>/ { count++ }
#    count > n { skip = 0 }
#    !skip { print }
#' "$input_file" > "$output_file"

#total += 0.00038 * length_sq + 0.0551 * len + 4.92
#total += 1.4856 * exp(0.0079 * length); 