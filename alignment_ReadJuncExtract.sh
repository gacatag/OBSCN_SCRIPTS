    bam_folder='PATH_TO_OUTPUT_FOLDERS'
# All fastq files are presuemd to be in a single directory. The format of the paired files should be 
 # Forward read files ending with _fwd.fq.gz and reverese files ending with _rev.fq.gz
 # If the format and arrangmenets of the files do not match please change and rearrange to the specified!
    fastq_folder='PATH_TO_INPUT_FASTQ_FILES'
    work_folder='PATH_TO_WORKING_DIRECTORY'
    N=NUMBER_OF_FILES_TO_PROCESS_SIMULTANEOUSLY
    runRegtool='PATH_TO_REGTOOL_EXECUTABLE_FILE'

cd ${work_folder}


r1List=$(ls -1 ${fastq_folder}/*unmapped_fwd.fq.gz | sed s/^.*\\/\//)


for r1 in ${r1List}; do

        name=${r1}
    r2=$(echo ${r1} | sed s/unmapped_fwd/unmapped_rev/)
    f1=$(echo ${r1}|sed s/^.*\\/\//|sed s/\\..*//)
    f2=$(echo ${r2}|sed s/^.*\\/\//|sed s/\\..*//)
    outFolder=$(echo ${f1}|sed s/unmapped_fwd//)

    
        ${runRegtool} junctions extract -s RF -o ${bam_folder}/${outFolder}/${outFolder}Aligned_junction.bed ${bam_folder}/${outFolder}/${outFolder}Aligned.sortedByCoord.out.bam &
     (( ++count % N == 0)) && wait
done

