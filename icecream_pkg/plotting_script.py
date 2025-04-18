
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import argparse
import os

def read_features_from_combined_file(filename,long_or_short):
    df = pd.read_csv(filename, sep="\t")
    features = []
    color_mapping = {"conjugation": "#FED800", "ARGs": "#DB4547", "DR": "gray","other ORFs":"#CEFCCC","Integrase":"#BC0FCD"}
    for _, row in df.iterrows():
        start = int(row['start'])
        end = int(row['end'])
        strand = int(row['strand'])
        color = color_mapping.get(row['feature_type'], "gray")
        label = row['long_label'] if long_or_short == 'long' else row['short_label']
        if pd.isna(label):
            label = None
        elif label == "None":
            label = None
        feature = GraphicFeature(
            start=start,
            end=end,
            strand=strand,
            color=color,
            label=label
        )
        features.append(feature)
    return features

def generate_plot(features, segment_start, segment_end, figure_width_scale,inputfile,outputfolder):
    record = GraphicRecord(sequence_length=segment_end-segment_start, first_index=segment_start, features=features)
    ax, _ = record.plot(figure_width=(segment_end-segment_start)/4000*figure_width_scale)
    input_folder, input_filename = os.path.split(inputfile)
    outputname = input_filename.replace("_combined_orfs.txt", "")
    output_filepath = os.path.join(outputfolder, f"{outputname}_graphic_arrangement.jpeg")
    ax.figure.savefig(output_filepath)
    

def plotting_script(combined_orfs_file,outputF,figure_width_scale=1,long_or_short="long"):
    first_input = pd.read_csv(combined_orfs_file, sep="\t")
    first_input['start'] = first_input['start'].astype(int)
    first_input['end'] = first_input['end'].astype(int)
    segment_start = first_input['start'].min()
    segment_end = first_input['end'].max()
    features = read_features_from_combined_file(combined_orfs_file,long_or_short)
    print(features)
    generate_plot(features, segment_start, segment_end, figure_width_scale,combined_orfs_file,outputF)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a plot from combined ORFs file.")
    parser.add_argument('-a',"--combined_orfs_file", help="Path to the combined ORFs file.")
    parser.add_argument('-o',"--outputF", help="Path to the output folder")
    parser.add_argument("--figure_width_scale", type=float, help="Width of the figure to be generated.", default=1)
    parser.add_argument("--long_or_short", choices=["short", "long"], default="long", help="short or long. Some of the genes labels are too long; you can choose to use short name. Default is long.")
    args = parser.parse_args()
    plotting_script(args.combined_orfs_file,args.outputF,args.figure_width_scale,args.long_or_short)
