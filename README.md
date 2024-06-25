These files must be created prior to running any analysis:
F1. bamFiles.txt: Includes list of paths to all the bam files. Each path in a separate line.
F2. batchVec.txt: List of character strings describing the different batches of sequencing.
F3. annoFile.xlsx: It includes one sheet, in which a table with 2 columns exist. The first row must include the columns names 'sample_ID', and 'group'. The first column under 'sample_ID' the name of the samples are listed, which they must be identical to the _EEES.tsv and _IMIS.tsv file names (without the '_EEES.tsv' and '_IMIS.tsv'.

The scripts must be run in the following order:

RNAseq read alignment:
1.ReferenceDownload.sh
2.alignment_index_Script.sh
3.alignment_StarScript.sh

Splicing analysis runs include:
4. buildReference.R
5. interestExonIncAnalysisEEES.R
6. interestIntronRetAnalysisIMIS.R
7. allExonInclusion_ObscnSpecific_DownstreamAnalysis.R
The functions defined in the following files are used in this script:
	7.1. ProcessExons_obscn.R
	7.2. MeasurePSI_obscn.R
	7.3. DifExInc_obscn.R
8. interestBasedOtherAlternativeSplicing_downstream.R


Exon Junction Analysis runs include:
9. alignmnet_ReadJuncExtract.sh
10. allExonInclusion_ObscnSpecific_DownstreamAnalysis_Junctionstudy_globalNorm.R