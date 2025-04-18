import os
import subprocess

# Define the path to the list of subfolder names
list_file = "list.gbk.exclude.plasmid.txt"

# Read the subfolder names from the list file
with open(list_file, 'r') as f:
    subfolders = [line.strip() for line in f.readlines()]

# Loop through each subfolder
for subfolder in subfolders:
    # Construct the path to the candidate folder
    candidate_folder = os.path.join(subfolder.replace(".gbk", ""), "candidate")
    
    # Check if the candidate folder exists
    if os.path.exists(candidate_folder):
        # Find all .faa files in the candidate folder
        for filename in os.listdir(candidate_folder):
            if filename.endswith(".faa"):
                # Construct the full path to the .faa file
                faa_file = os.path.join(candidate_folder, filename)
                
                # Create output folder name by appending "_t4ss" to the .faa filename (without extension)
                output_folder = os.path.join(candidate_folder, filename.replace(".faa", "_t4ss"))
                
                # Construct the macsyfinder command
                macsyfinder_cmd = [
                    "macsyfinder",
                    "--db-type", "ordered_replicon",
                    "--sequence-db", faa_file,
                    "--models", "CONJScan/Chromosome",
                    "all",
                    "-o", output_folder
                ]
                
                # Run the macsyfinder command
                try:
                    print(f"Running MacSyFinder for: {faa_file}")
                    subprocess.run(macsyfinder_cmd, check=True)
                    print(f"Output saved to: {output_folder}")
                except subprocess.CalledProcessError as e:
                    print(f"Error running MacSyFinder for {faa_file}: {e}")
    else:
        print(f"Candidate folder not found: {candidate_folder}")
