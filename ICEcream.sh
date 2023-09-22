#!/usr/bin/env bash
set -e

inputF=""
tempF=""
outputF=""
SCRIPT_VERSION="1.10"

remove_trailing_slash() {
    echo "$1" | sed 's:/*$::'
}

display_help() {
    echo "Usage: bash $0 --input <input_folder> --temp <temp_path> --output <output_path>"
    echo ""
    echo "Options:"
    echo "  --input   Path to input folder containing .gbk files"
    echo "  --temp    Path for generating a new TEMPORARY folder"
    echo "  --output  Path for generating a new OUTPUT folder"
	echo "  --version Display script version"
    echo "  --help    Display this help message."
    exit 1
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input)
        inputF=$(remove_trailing_slash "$2")
        shift
        shift
        ;;
        --temp)
        tempF=$(remove_trailing_slash "$2")
        shift
        shift
        ;;
        --output)
        outputF=$(remove_trailing_slash "$2")
        shift
        shift
        ;;
		--help|-h)
        display_help
        ;;
	    --version|-v)
        echo "Script Version: $SCRIPT_VERSION"
        exit 0
        ;;
        *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

if [[ -z $inputF ]] || [[ -z $tempF ]] || [[ -z $outputF ]]; then
    echo "Error: Missing required arguments."
    display_help
fi

folder_not_empty() {
    local folder="$1"
    [[ -d "$folder" && $(ls -A "$folder") ]]
}

if folder_not_empty "${tempF}" || folder_not_empty "${outputF}"; then
    echo 'Error: The temporary folder or the output folder already exist and are not empty, please delete or clear them.'
else
    mkdir -p "${tempF}"
    mkdir -p "${outputF}"
fi


terminal_width=$(tput cols)
welcome_lines=(
    "Welcome to the script of ICEcream"
	"========================"
	"Author: Xiaole (Charlotte) YIN"
	"E-mail: yinlele99@gmail.com"
	"Hong Kong University, Harvard Medical School"
	"Version 1.10, last update in September 2023" 
)
# Calculate padding for the "Welcome" text
printf '=%.0s' $(seq 1 $terminal_width)
echo ""
# Iterate over the lines and print each one centered
for line in "${welcome_lines[@]}"; do
    padding=$(( (terminal_width - ${#line}) / 2 ))
    printf '%*s\n' $((padding + ${#line})) "$line"
done
printf '=%.0s' $(seq 1 $terminal_width)

echo "=== Step 1: locate ICEs"
# run_icefinder
run_icefinder() {
	perl ICEfinder_modified4_yxl.pl -i "$inputF" -t "$tempF" -o "$outputF" > "$tempF"/icefinder.log
    cat "$outputF"/*/*_summary.txt > "$outputF"/icefinder.result.summary.txt
    local ice_count=$(wc -l < "$outputF/icefinder.result.summary.txt")
    echo "=== Step 1.1: icefinder finished, how many ICEs identified: $ice_count"
}
run_icefinder

run_icefinder_amend(){
	grep 'ERROR' "$tempF/icefinder.log" | rev | cut -d '/' -f 1 | rev | sed 's/.ptt file!//g' > "$tempF/icefinder_first_error.txt"
    local potential_ice_count=$(wc -l < "$tempF/icefinder_first_error.txt")
    echo "=== Step 1.1:                     how many ICEs with another chance to be identified after manual annotation: $potential_ice_count"
	if [ "$potential_ice_count" -ne 0 ]; then
		echo "=== Step 1.2: parse, prokka and redo icefinder"
		echo "=== Step 1.2: parse"
		python scripts/GbffParser_YXL_2023.py -i "$inputF" -f .gbk -o "$tempF" > "${tempF}/parser.log"
		mkdir -p "${tempF}/prokka"
		mkdir -p "${tempF}/prokka_out"
		while read -r line; do 
		    cp "${tempF}/${line}/${line}.fna" "${tempF}/prokka"
		done < "${tempF}/icefinder_first_error.txt"
		echo "=== Step 1.2: prokka"
        while read -r line; do 
            prokka --outdir "${tempF}/prokka_out" --force --prefix "$line" --cpus 20 "${tempF}/prokka/${line}.fna" --quiet
        done < "${tempF}/icefinder_first_error.txt" > "${tempF}/prokka.log"
        echo "=== Step 1.2: redo icefinder"
		perl ICEfinder_modified4_yxl.pl -i "${tempF}"/prokka_out -t "${tempF}"/prokka_tmp -o "${tempF}"/prokka_icefinder > "${tempF}/icefinder.log.prokka"
        if [ -f "${tempF}"/prokka_icefinder/*/*_summary.txt ]; then
		    cat "${tempF}"/prokka_icefinder/*/*_summary.txt >> "$outputF/icefinder.result.summary.txt"
			local new_ice_count=$(wc -l < "${tempF}"/prokka_icefinder/*/*_summary.txt)
			echo "=== Step 1.2: finished,how many ICEs identified in this step: $new_ice_count"
			local ice_count=$(wc -l < "$outputF/icefinder.result.summary.txt")
			echo "=== Step 1.2:          how many ICEs identified in total: $ice_count"
		else
		    echo "=== Step 1.2: finished,how many ICEs identified in this step: 0"
			local ice_count=$(wc -l < "$outputF/icefinder.result.summary.txt")
			echo "=== Step 1.2:          how many ICEs identified in total: $ice_count"
		fi
	else
	    echo "=== Step 1.2: skip"
	fi
}
run_icefinder_amend

run_icefamily_identify(){
	echo "=== Step 2: identify ICE families"
    echo "=== Step 2.1: exclude potential plasmid by searching keywords"
	awk '{print $1}' $outputF/icefinder.result.summary.txt > ${tempF}/icefinder.result.summary_ids.txt
	while read -r id; do
	    line="${id}.gbk"
        if ! grep -q "^${line}$" ${tempF}/split.result.summary.noplasmid.txt && grep -q "^${line}$" ${tempF}/split.result.summary.txt; then
            grep "^$id" $outputF/icefinder.result.summary.txt >> $outputF/ICEs_but_truly_plasmids.txt
        else
			echo "=== Step 2.2: studying $line"
			echo "=== Step 2.2: studying $line prodigal"
			local fasta_file="${tempF}/${line}.fa"
			local prodigal_output="${tempF}/${line}.fa_prodigal.faa"
			prodigal -p meta -i "$fasta_file" -a "$prodigal_output" -q >> ${tempF}/icefamily.identification.log
			if [ -f "${fasta_file}a" ]; then
				python ICEfamily_refer/familytools/amendORF.py "${fasta_file}a" "$prodigal_output"
				cat "${fasta_file}a" "${tempF}/${line}.faa_append.faa" > "${tempF}/${line}.fa2.faa"
			else
				python ICEfamily_refer/familytools/amendORF2.py "$fasta_file"
			fi
			echo "=== Step 2.3: studying $line hmmscan"
			hmmscan --tblout "${tempF}/${line}.fa2.faa.icefinder.hmmscan.out" --cpu 20 data/ICE.hmm.db "${tempF}/${line}.fa2.faa" >> ${tempF}/icefamily.identification.log
			grep -v "#" "${tempF}/${line}.fa2.faa.icefinder.hmmscan.out" | tr -s ' ' '\t' > "${tempF}/${line}.fa2.faa.icefinder.hmmscan.modified.out"
			hmmscan --tblout "${tempF}/${line}.fa2.faa.YXL.hmmscan.out" --cpu 20 ICEfamily_refer/conjugation/YXL_26family.hmm "${tempF}/${line}.fa2.faa" >> ${tempF}/icefamily.identification.log
			grep -v "#" "${tempF}/${line}.fa2.faa.YXL.hmmscan.out" | tr -s ' ' '\t' > "${tempF}/${line}.fa2.faa.YXL.hmmscan.modified.out"
			cat "${tempF}/${line}.fa2.faa.icefinder.hmmscan.modified.out" "${tempF}/${line}.fa2.faa.YXL.hmmscan.modified.out" > "${tempF}/${line}.fa2.faa.merge.all.hmmscan.out"
			echo "=== Step 2.4: studying $line align and refer to the ICEcream database"
			mkdir -p "${outputF}/${id}_ICEfamily"
			Rscript ICEfamily_refer/familytools/merge.version7.R "${tempF}/${line}.fa2.faa.merge.all.hmmscan.out" 1e-5 "${outputF}/icefinder.result.summary.txt" "${outputF}" "${tempF}" 
			cat ${outputF}/${id}_ICEfamily/${id}_ICEfamily_result.txt >>${outputF}/icefamily.result.summary.txt
        fi
    done < ${tempF}/icefinder.result.summary_ids.txt 
}
run_icefamily_identify

run_identify_ARG(){
	echo "=== Step 3: annotate ARGs on ICEs"
	while read -r line; do
	    echo "=== Step 3: annotate ARGs on ICEs $line"
	    blastp -db ICEfamily_refer/SARG_20200618_with_multicomponent.fasta -query "${tempF}/${line}.gbk.fa2.faa" -out "${tempF}/${line}.gbk.fa2.faa.tab" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -max_target_seqs 5
        Rscript ICEfamily_refer/familytools/blastx_out_treatment_linux_20230920.R "${tempF}/${line}.gbk.fa2.faa.tab" 70 80 ICEfamily_refer/familytools/structure_20200330_includingmutation.txt 1e-10 "$tempF"
	done < ${tempF}/icefinder.result.summary_ids.txt 
}
run_identify_ARG

draw_in_line(){
	echo "=== Step 4: plot ORFs arrangment of identifed ICEs"
	while read -r line; do
	    echo "=== Step 4: plot ORFs arrangment of identifed ICEs $line"
		read detail file4 start end < <(awk -F'\t' '
		NR==1 {
			for (i=1; i<=NF; i++) {
				col[$i] = i; }
		} 
		NR==2 {
			print $(col["details_file_number"]), $(col["icefinder_folder"]), $(col["ICE_start"]), $(col["ICE_end"]);
		}' $outputF/${line}_ICEfamily/${line}_ICEfamily_result.txt)
		IFS='_' read -ra parts <<< "$detail"
		if [[ ${#parts[@]} -eq 2 ]]; then
		    detail_name="${parts[1]}"
			python ICEfamily_refer/familytools/organize_orfs.py "${start}..${end}" ICEfamily_refer/familytools/long_short_label_plot.txt "${outputF}/${line}_ICEfamily/${line}_classification_summary_details${detail_name}.txt" "${tempF}/extracted_classification_${line}.gbk.fa2.faa.tab.txt" ${outputF}/${line}/${file4} "${outputF}/${line}"
			python ICEfamily_refer/familytools/plotting_script.py "${outputF}/${line}_combined_orfs.txt" "${outputF}"
		elif [[ ${#parts[@]} -gt 2 ]]; then
		    echo "Executing another type of command"
		else
		    echo "Unexpected format for drawing, check the $outputF/${line}_classification_summary_modifi.txt"
		fi		
	done < ${tempF}/icefinder.result.summary_ids.txt 
}
draw_in_line
echo "=== Great! Task finished!"
printf '=%.0s' $(seq 1 $terminal_width)

