import os
import shutil
from math import ceil
import argparse

def main():
    parser = argparse.ArgumentParser(description='folder to n subfolders')
    parser.add_argument('-i', required=True, help='Input folder')
    parser.add_argument('-n', required=True, help='How many subfolders')
    parser.add_argument('-p',  default="subfolder_", help='prefix')
    args = parser.parse_args()
    source_folder = args.i
    n_subfolders = int(args.n)
    output_folder_prefix = args.p
    # Get the list of all files in the source folder
    files = [f for f in os.listdir(source_folder) if os.path.isfile(os.path.join(source_folder, f))]
    
    # Calculate how many files should be in each subfolder
    num_files_per_subfolder = ceil(len(files) / n_subfolders)
    
    # Create the subfolders and distribute the files
    for i in range(n_subfolders):
        subfolder_name = f"{output_folder_prefix}{i+1}"
        subfolder_path = os.path.join(source_folder, subfolder_name)
        os.makedirs(subfolder_path, exist_ok=True)  # Create the subfolder
        
        # Get the files that should go in this subfolder
        start_idx = i * num_files_per_subfolder
        end_idx = start_idx + num_files_per_subfolder
        files_to_move = files[start_idx:end_idx]
        
        # Move the files to the subfolder
        for file_name in files_to_move:
            shutil.move(os.path.join(source_folder, file_name), os.path.join(subfolder_path, file_name))

    print(f"Split folder '{source_folder}' into {n_subfolders} subfolders.")

if __name__ == "__main__":
    main()
    
