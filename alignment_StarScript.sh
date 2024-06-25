
    runStar='PATH_TO_STAR_EXECUTABLE_FILE'
    genomeDir='PATH_TO_STAR_INDEX_FILES'
    bam_folder='PATH_TO_OUTPUT_FOLDERS'
 # All fastq files are presuemd to be in a single directory. The format of the paired files should be 
 # Forward read files ending with _fwd.fq.gz and reverese files ending with _rev.fq.gz
 # If the format and arrangmenets of the files do not match please change and rearrange to the specified!
    fastq_folder='PATH_TO_INPUT_FASTQ_FILES'
    work_folder='PATH_TO_WORKING_DIRECTORY'
    N=NUMBER_OF_FILES_TO_PROCESS_SIMULTANEOUSLY
    threadNo=NUMBER_OF_THREADS_ASSIGNED_TO_ANALYZING_EACH_FILE



cd ${work_folder}


r1List=$(ls -1 ${fastq_folder}/*unmapped_fwd.fq.gz | sed s/^.*\\/\//)


ulimit -s unlimited
#### STEP 1: STAR alignment of trimmed fastq files and junctions from the bam file

count=0

for r1 in ${r1List}; do
	name=${r1}
	r2=$(echo ${r1} | sed s/unmapped_fwd/unmapped_rev/)
	f1=$(echo ${r1}|sed s/^.*\\/\//|sed s/\\..*//)
	f2=$(echo ${r2}|sed s/^.*\\/\//|sed s/\\..*//)
	outFolder=$(echo ${f1}|sed s/unmapped_fwd//)
	mkdir $bam_folder/${outFolder}
	echo "STAR alignment for ${outFolder} starts"
        $runStar --genomeDir ${genomeDir}/ \
	        --readFilesIn ${fastq_folder}/${r1} ${fastq_folder}/${r2} \
	        --readFilesCommand zcat \
	        --genomeLoad NoSharedMemory \
	        --outFileNamePrefix ${bam_folder}/${outFolder}/${outFolder} \
        	--outFilterType BySJout \
	        --outFilterMultimapNmax 20 \
	        --alignSJoverhangMin 8 \
	        --alignSJDBoverhangMin 1 \
	        --outFilterMismatchNmax 999 \
       		--outFilterMismatchNoverReadLmax 0.04 \
	        --alignIntronMin 20 \
        	--alignIntronMax 1000000 \
	        --alignMatesGapMax 1000000 \
	        --outSAMunmapped Within \
	        --outSAMattributes All \
	        --outSAMtype BAM SortedByCoordinate \
	        --outSAMheaderHD @HD VN:1.6 SO:coordinate \
	        --sjdbScore 1 --outMultimapperOrder Random \
	        --limitOutSJcollapsed 2000000 \
        	--limitBAMsortRAM 16106127360 \
	        --twopassMode Basic \
        	--runThreadN ${threadNo} &
	 (( ++count % N == 0)) && wait
done


# Indexing bam files
runSamtools=/mnt/data/software/samtools/bin/samtools
N=6

for r1 in ${r1List}; do
	name=${r1}
	r2=$(echo ${r1} | sed s/unmapped_fwd/unmapped_rev/)
	f1=$(echo ${r1}|sed s/^.*\\/\//|sed s/\\..*//)
	f2=$(echo ${r2}|sed s/^.*\\/\//|sed s/\\..*//)
	outFolder=$(echo ${f1}|sed s/unmapped_fwd//)
        $runSamtools index ${bam_folder}/${outFolder}/${outFolder}Aligned.sortedByCoord.out.bam &
	 (( ++count % N == 0)) && wait
done

EF $GTF

#        echo "STAR alignment and junctions for ${f} ends"

