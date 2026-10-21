#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    echo "Extract discordantly aligned fragments from sorted BAM file"
    exit 1
fi

input_dir="$1"
output_dir="$2"

# check input_dir is available
[ -d "$input_dir" ] || { echo "error: input_dir is not available: $input_dir"; exit 1; }

# build the output_dir 
mkdir -p "$output_dir"

# process files
for input_bam in "$input_dir"/*.sorted.bam; do
    # check file is available
    [ -f "$input_bam" ] || continue

    # output file name
    base_name=$(basename "$input_bam" .sorted.bam)
    output_file2="$output_dir/${base_name}.discordant_pairs.bam"
    
    qname_sorted_bam="$input_dir/${base_name}.qname_sorted.bam"
     
    
    sort_order=$(samtools view -H "$input_bam" | grep "^@HD" | grep -o "SO:[a-z]*" | cut -d: -f2 || echo "unknown")
    
    if [ "$sort_order" = "queryname" ]; then
        qname_sorted_bam="$input_bam"
    else
        samtools sort -n -@ 16 "$input_bam" -o "$qname_sorted_bam"
        input_bam="$qname_sorted_bam"
    fi
    
   
    echo "Processing $input_bam..."
    
    ## extract discordant read pairs
    echo "Step 2: Extracting discordant read pairs..."

    # firstly extract properly paired reads, then figure out discordant pairs
    samtools view -h -F 0x2 -F 0x100 -F 0x800 "$input_bam" | \
    awk '
    BEGIN {
        OFS="\t"
        MAX_INSERT_SIZE = 1000
    }
    /^@/ {
        print
        next
    }
    {
        qname = $1
        flag = $2
        rname = $3
        pos = $4
        mapq = $5
        cigar = $6
        rnext = $7
        pnext = $8
        tlen = $9
        
        # if qname is the second read of the pair
        if (qname in data) {
            # process read pair
            split(data[qname], prev, "\t")
            
            chr1 = prev[3]; pos1 = prev[4]; flag1 = prev[2]
            chr2 = rname;   pos2 = pos;     flag2 = flag
            
            is_discordant = 0
            reason = ""
            
            is_r1_reverse = (and(flag1, 16))
            is_r2_reverse = (and(flag2, 16))
            # from different chromosomes
            if(chr1 != "*" && chr2 != "*"){
            	print data[qname] "\tZD:Z:" reason
                print $0 "\tZD:Z:" reason
            }
            
            delete data[qname]
        } else {
            # store the first read
            data[qname] = $0
        }
    }
    ' | samtools view -b - | \
    samtools sort -n -@ 4 -o "$output_file2" -

done

echo "All files processed successfully!"
