import os
import csv
from collections import Counter


'''
#move files to appropriate folders
for file in *.pae *.pdb; do
  id=$(echo "$file" | cut -d'_' -f1)
  extension="${file##*.}"
  mkdir -p "$id/$extension"
  mv "$file" "$id/$extension/"
done
'''



# Define the root directory containing the id folders
root_dir = "/scratch/brussel/vo/000/bvo00023/vsc10579/achromatium/clean_achromatium/structures"

complete_ass_l=['bin9','bin51','bin32','Cell7','Cell12','bin0','Cell8','bin2','bin3','bin8','Cell14','bin6','bin33','bin18','Cell26','Cell3']
# Prepare the data structure for 'output'
output_data = []

# Iterate over each "id" folder
for id_folder in os.listdir(root_dir):
    id_path = os.path.join(root_dir, id_folder)
    pdb_path = os.path.join(id_path, "pdb")

    # Skip if it's not a directory or "pdb" folder doesn't exist
    if not os.path.isdir(id_path) or not os.path.isdir(pdb_path):
        continue

    # Get all files in the pdb folder
    files = [f for f in os.listdir(pdb_path) if os.path.isfile(os.path.join(pdb_path, f))]
    
    # Count total occurrences
    occurrence = len(files)

    # Extract second elements from file names
    names = [f.split("_")[1] for f in files if "_" in f]
    
    # Count unique occurrences
    unique_occurrence = len(set(names))
    
    aid = ','.join(sorted(list(set(names))))
    # Append the result for this id
    output_data.append([id_folder, occurrence, unique_occurrence, aid])

# Write results to a TSV file
output_file = root_dir+"/structure_counter.tsv"
with open(output_file, "w", newline="") as tsv_file:
    writer = csv.writer(tsv_file, delimiter="\t")
    writer.writerow(["id", "occurrence", "unique_occurrence", 'assemby_id'])
    writer.writerows(output_data)

print(f"Indexing complete. Output written to {output_file}")
