#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    echo "Extracts split-aligned fragments from sorted BAM file"
    exit 1
fi

input_dir="$1"
output_dir="$2"
# check input_dir is available
[ -d "$input_dir" ] || { echo "error: input_dir is not available: $input_dir"; exit 1; }
# build the output_dir 
mkdir -p "$output_dir"
mkdir -p "temp_dir"

# process files
for input_bam in "$input_dir"/*.sorted.bam; do
    
    [ -f "$input_bam" ] || continue 

	if [ ! -f "$input_bam" ]; then
   	 echo "Error: Input BAM file '$input_bam' not found"
   	 exit 1
	fi

	if [ ! -f "${input_bam}.bai" ]; then
   	 echo "Creating BAM index..."
    	samtools index "$input_bam"
	fi
	echo "Processing $input_bam..."

	## extract read pairs containing split alignments
	name=$(basename "$input_bam" .sorted.bam)
	output_file1="$output_dir/$name.split_reads.bam"
    
    
    temp_qnames="temp_dir/${name}.qnames.txt"

	echo "Step 1: Extracting qnames of split alignments to file..."

	samtools view "$input_bam" | awk '/SA:Z:/ {print $1}' | sort -u > "$temp_qnames"

    
	if [ -s "$temp_qnames" ]; then
    	echo "Step 2: Extracting read pairs using qname file..."
    	samtools view -h -N "$temp_qnames" "$input_bam" | \
    	samtools view -b - | \
   	    samtools sort -@ 4 -o "$output_file1" -
        samtools index "$output_file1"
    else
        echo "No split alignments found for $name. Creating empty BAM."

        samtools view -H "$input_bam" | samtools view -b - > "$output_file1"
        samtools index "$output_file1"
	fi
    
    echo "finish: $input_bam -> $output_file1"
    
    rm "$temp_qnames"

done

# remove temp dir
rm -rf "temp_dir"

echo "====== All samples finished ======"
