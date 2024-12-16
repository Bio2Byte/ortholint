from Bio import SeqIO
from collections import defaultdict
import glob
import os
import sys

# Input folder containing the nested structure
input_folder = sys.argv[1]#"."
# Output file for the report
report_file = "non_unique_report.txt"

# Use glob to find all .faa files in nested subfolders
fasta_files = glob.glob(f"{input_folder}/**/*.faa", recursive=True)

# Open the report file for writing
with open(report_file, "w") as report:
    # Process each file found by glob
    for fasta_path in fasta_files:
        seq_dict = defaultdict(list)

        # Read the FASTA file
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Extract the first part of the label before the first space
            label = record.id.split(" ")[0]
            seq_dict[label].append(str(record.seq))

        # Count non-unique labels
        non_unique_count = sum(1 for val in seq_dict.values() if len(val) > 1)

        # Write the count to the report
        relative_path = os.path.relpath(fasta_path, input_folder)
        report.write(f"File: {relative_path} - Duplicated labels: {non_unique_count}\n")

print(f"Report saved to {report_file}")
