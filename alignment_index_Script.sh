wrkPath='PATH_TO_WORKING_DIRECTORY'
genomeFasta='PATH_TO_GENOME_FASTA_FILE'
annotationGtf='PATH_TO_GENE_ANNOTATION_GTF_FILE'
genomeDir='PATH_TO_INDEX_OUTPUT'
starRun='PATH_TO_STAR_EXECUTABLE_FILE'
threadNo=NUMBER_OF_THREADS

#mkdir ${genomeDir}
cd ${wrkPath}

${starRun} \
--runMode genomeGenerate \
--runThreadN ${threadNo} \
--genomeDir ${genomeDir} \
--genomeFastaFiles ${inFilePath}/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ${inFilePath}/gencode.v41.primary_assembly.annotation.gtf \
--sjdbOverhang 150

echo "Job Finished"
