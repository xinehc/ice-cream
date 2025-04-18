
import pandas as pd
import argparse
import os

def extract_segment_from_file(segment_info_str):
    segment_start, segment_end = map(int, segment_info_str.split(".."))
    return segment_start, segment_end

def subset_orfs_within_range_icefinder(file_path):
    input_df = pd.read_csv(file_path, sep="\t", header=None, skiprows=7)
    start_values = []
    end_values = []
    for value in input_df[0]:
        parts = value.split('..')
        start_values.append(parts[0].replace("Insert: ", "").strip())
        end_values.append(parts[1].strip() if len(parts) > 1 else None)
    input_df['start'] = start_values
    input_df['end'] = end_values
    input_df['strand'] = input_df[1].map({'+': 1, '-': -1})
    filtered_df = input_df[input_df[5] != "Integrase"]
    filtered_df = filtered_df[['start', 'end', 'strand']]
    filtered_df['long_label'] = "None"
    filtered_df['short_label'] = "None"
    filtered_df['feature_type'] = "other ORFs"
    subset = input_df[input_df[5] == "Integrase"]
    subset = subset[['start', 'end', 'strand']]
    subset['long_label'] = "Integrase"
    subset['short_label'] = "Integrase"
    subset['feature_type'] = "Integrase"
    subset = pd.concat([filtered_df, subset], ignore_index=True)
    return subset

# Function to subset ORFs within range and generate strand and label columns
def subset_familyorfs_within_range(file_path, segment_start, segment_end, start_col, end_col, feature_type,mapping_dict):
    input_df = pd.read_csv(file_path, sep="\t")
    subset = input_df[
        (input_df[start_col] >= segment_start) & (input_df[end_col] <= segment_end)
    ]
    subset['strand'] = subset.apply(lambda x: 1 if x[start_col] < x[end_col] else -1, axis=1)
    subset['long_label'] = subset["model_name"]
    subset["short_label"] = subset["long_label"].map(mapping_dict)
    subset = subset[['start', 'end', 'strand', 'long_label', 'short_label']]
    subset['feature_type'] = feature_type
    return subset

# Function to subset the third input file
def subset_blast_input(file_path, segment_start, segment_end):
    input_df = pd.read_csv(file_path, sep="\t")
    subset = input_df[input_df["query"].str.split("_").str[-2].astype(int) >= segment_start]
    subset = subset[subset["query"].str.split("_").str[-1].astype(int) <= segment_end]
    subset['start'] = subset["query"].str.split("_").str[-2].astype(int)
    subset['end'] = subset["query"].str.split("_").str[-1].astype(int)
    subset['strand'] = subset.apply(lambda x: 1 if x["querystart"] < x["queryend"] else -1, axis=1)
    subset['long_label'] = "ARG:" + subset["New.Subtype"]
    subset['short_label'] = subset['New.New.Subtype']
    subset = subset[['start', 'end', 'strand', 'long_label', 'short_label']]
    subset['feature_type'] = "ARGs"
    return subset

def organize_orfs(range,b,c,d,p):
    db_dir = os.environ.get('DATABASE_FOLDER')
    long_short_label_df = pd.read_csv(os.path.join(db_dir,'long_short_label_plot.txt'), sep="\t")
    mapping_dict = dict(zip(long_short_label_df.long_label, long_short_label_df.short_label))

    segment_start, segment_end = extract_segment_from_file(range)
    
    # Subset ORFs for the second file
    second_subset = subset_familyorfs_within_range(
        b,
        segment_start,
        segment_end,
        "start",
        "end",
        "conjugation",
        mapping_dict
    )
    
    # Subset ORFs for the third file
    third_subset = subset_blast_input(
        c,
        segment_start,
        segment_end
    )
    # for the fourth file
    fourth_subset = subset_orfs_within_range_icefinder(
        d
    )
     
    # Create a feature for DR from the first file
    dr_feature = pd.DataFrame({
        'start': [segment_start, segment_end],
        'end': [segment_start+10, segment_end+10],
        'strand': [1, 1],
        'long_label': ['DR', 'DR'],
        'short_label': ['DR', 'DR'],
        'feature_type': ['DR', 'DR']
    })
    
    # Combine and save ORFs to a single file
    combined_orfs = pd.concat([second_subset, third_subset, dr_feature, fourth_subset], ignore_index=True)
    combined_orfs['start'] = combined_orfs['start'].astype(int)
    combined_orfs = combined_orfs.sort_values(by="start")
    combined_orfs = combined_orfs.drop_duplicates(subset='start', keep='first')
    combined_orfs.to_csv(p+"_combined_orfs.txt", sep="\t", index=False)

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organize ORFs from provided files to generate combined ORFs file.")
    parser.add_argument("-range", help="string for start..end")
    parser.add_argument("-b", help="Path to the icefamily details input file.")
    parser.add_argument("-c", help="Path to the resistance blast input file.")
    parser.add_argument("-d", help="Path to the ice result input file.")
    parser.add_argument("-p", help="prefix for output")
    args = parser.parse_args()
    organize_orfs(args.range,args.b,args.c,args.d,args.p)


