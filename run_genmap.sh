
# Generate chromosome sizes file using seqkit
generate_chrom_sizes() {
    chrName=$1
    chrLenght=$2
    chrom_sizes_file=$3	

    echo "Generating chromosome sizes file: $chrom_sizes_file"
    paste $chrName $chrLength > $chrom_sizes_file
}

# Iterate over each fasta file matching *toplevel.fa
for fasta_file in $(find . -type f -name "*toplevel.fa"); do
    # Get the directory of the fasta file
    dir=$(dirname "$fasta_file")
    
    # Get the base name of the fasta file without the extension
    base_name=$(basename "$fasta_file" .fa)
    
    # Create the index file name
    index_file="${dir}/${base_name}_index"
    
    # Generate chromosome sizes file
    chrom_sizes_file="${dir}/${base_name}_chrom_sizes.txt"
    generate_chrom_sizes "${dir}/STAR_index/genomeDir/chrName.txt" "${dir}/STAR_index/genomeDir/chrLenght.txt" "$chrom_sizes_file"

    # if ${dir}/21mers_mappability.wig exists, skip this directory
    if [ -f "${dir}/21mers_mappability.wig" ]; then

    for wig_file in $(find "$dir/" -type f -name "*.wig"); do

	# Get the base name of the WIG file without the extension
        wig_base_name=$(basename "$wig_file" .wig)
        echo $wig_base_name
        
    	chrom_sizes_file="${dir}/${wig_base_name}.chrom.sizes"
	# Create the BigWig file name
        bw_file="${dir}/${wig_base_name}.bw"
        
        # Convert WIG to BigWig
        wigToBigWig "$wig_file" "$chrom_sizes_file" "$bw_file"
        
        echo "Converted $wig_file to $bw_file"
    done
    continue
    fi
    
    # Create the index for genmap
    genmap index --fasta-file "$fasta_file" --index "$index_file"
    
    # Generate mappability tracks (WIG files)
    genmap map -K 21 -E 2 -I "$index_file" -O "${dir}/21mers_mappability" -t -w -bg 
    genmap map -K 28 -E 2 -I "$index_file" -O "${dir}/28mers_mappability" -t -w -bg 
    genmap map -K 35 -E 2 -I "$index_file" -O "${dir}/35mers_mappability" -t -w -bg 
    genmap map -K 50 -E 2 -I "$index_file" -O "${dir}/50mers_mappability" -t -w -bg 
    
    # Convert the resulting WIG files to BigWig files
    for wig_file in $(find "$dir/" -type f -name "*.wig"); do
        # Get the base name of the WIG file without the extension
        wig_base_name=$(basename "$wig_file" .wig)
        
        # Create the BigWig file name
        bw_file="${dir}/${wig_base_name}.bw"
        
        # Convert WIG to BigWig
        wigToBigWig "$wig_file" "$chrom_sizes_file" "$bw_file"
        
        echo "Converted $wig_file to $bw_file"
    done
done
