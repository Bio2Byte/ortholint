import os
import re
import csv
import sys

def extract_data_from_file(file_path):
    """Extract relevant data from an *.iqtree file."""
    data = {
        "file_name": None,
        "sum_len_total_branches": None,
        "sum_len_internal_branches": None,
        "p_compatcness": None,
        "n_seqs": None,
        "n_sites": None,
        "n_constant_sites": None,
        "n_parsimony_informative_sites": None,
        "n_distinct_patterns": None
    }
    
    # Extract the abbreviated file name
    file_name = os.path.basename(file_path).split('_')[-1].replace('.fasta.iqtree', '')
    data["file_name"] = file_name

    # Read the file and parse lines of interest
    with open(file_path, 'r') as f:
        content = f.read()
        # Regex patterns for each value
        patterns = {
            "sum_len_total_branches": r"Total tree length \(sum of branch lengths\): ([\d\.]+)",
            "sum_len_internal_branches": r"Sum of internal branch lengths: ([\d\.]+)",
            "p_compatcness": r"\(([\d\.]+)% of tree length\)",
            "n_seqs": r"Input data: (\d+) sequences",
            "n_sites": r"Input data: \d+ sequences with (\d+) amino-acid sites",
            "n_constant_sites": r"Number of constant sites: (\d+)",
            "n_parsimony_informative_sites": r"Number of parsimony informative sites: (\d+)",
            "n_distinct_patterns": r"Number of distinct site patterns: (\d+)"
        }
        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                data[key] = match.group(1)
    
    return data

def find_iqtree_files(base_path):
    """Find all *.iqtree files recursively starting from base_path."""
    iqtree_files = []
    for root, _, files in os.walk(base_path):
        for file in files:
            if 'work' not in root:
                if file.endswith(".fasta.iqtree"):
                    iqtree_files.append(os.path.join(root, file))
    print(iqtree_files)
    return iqtree_files

def main(base_path, output_csv):
    """Main function to generate the CSV table."""
    iqtree_files = find_iqtree_files(base_path)
    data_rows = []
    
    for file_path in iqtree_files:
        extracted_data = extract_data_from_file(file_path)
        data_rows.append([
            extracted_data["file_name"],
            extracted_data["sum_len_total_branches"],
            extracted_data["sum_len_internal_branches"],
            extracted_data["p_compatcness"],
            extracted_data["n_seqs"],
            extracted_data["n_sites"],
            extracted_data["n_constant_sites"],
            extracted_data["n_parsimony_informative_sites"],
            extracted_data["n_distinct_patterns"]
        ])
    
    # Write to CSV
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([
            "file_name", "sum_len_total_branches", "sum_len_internal_branches",
            "p_compatcness", "n_seqs", "n_sites", "n_constant_sites",
            "n_parsimony_informative_sites", "n_distinct_patterns"
        ])
        writer.writerows(data_rows)
    
    print(f"Data successfully written to {output_csv}")

# Replace 'your_base_path' with the starting directory and 'output.csv' with the desired output file name
if __name__ == "__main__":
    base_path = sys.argv[1]  # Replace with the base directory path
    output_csv = sys.argv[2] +"_tree_stats.csv"    # Replace with the desired CSV file name
    main(base_path, output_csv)
