#!/bin/bash

if [ $# -lt 5 ]; then
	echo "$0 [cultivar] [Lower DP limit for high-confidence variant calling (peak-0.5IQR)] [Upper DP limit for high-confidence variant calling (peak+0.5IQR)] [Upper DP limit for variant calling (peak+1.5IQR)] [Preprocess step (0:start from Trimmomatic 1:use reads after Trimmomatic)]"
	echo "ex. $0 Koshihikari-SBL1 39 59 79 0"
	exit 0
fi

CULTIVAR=$1
DP_LOWER_LIMIT_HC=$2
DP_UPPER_LIMIT_HC=$3
DP_UPPER_LIMIT=$4
PP_STEP=$5

echo "${CULTIVAR} DP_LOWER_LIMIT_HC=[${DP_LOWER_LIMIT_HC}] DP_UPPER_LIMIT_HC=[${DP_UPPER_LIMIT_HC}] DP_UPPER_LIMIT=[${DP_UPPER_LIMIT}] PP_STEP=[${PP_STEP}] start...."
date

# for PP_STEP=0
if [ $PP_STEP -eq 0 ]; then
	READ_DIR=/lustre/home/ykawahar/pft/H27_analysis/KoshihikariRelatedGenome11/$CULTIVAR
# for PP_STEP=1
else
	READ_DIR=/lustre/home/ykawahar/pft/H27_analysis/KoshihikariRelatedGenome11/$CULTIVAR
fi

#===== Parmeters
NCORE=10						# NCORE for alignment
NCORE_VC=10					# NCORE for variant call
MAX_LOOP=20					# max iteration
JAVA_MEM="-Xmx256G -Xms32G"	# Java VM heap 
SAMTOOLS_MEM="24G"				# max RAM for samtools sort


INPUT_FASTQ_PHRED_OFFSET=33	# phred offset value for input FASTQ
ALLOW_VAR_COUNT=0			# max. number of variations (SNP+INDEL) for finishing iteration

#=============================================================================== output file
STAT_OUTFILE=${CULTIVAR}.stat
OUTFILEBASE_1ST=${CULTIVAR}_1st
OUTFILEBASE_FIN=${CULTIVAR}_Final
UNIFIED_ALL_VCF=UnifiedAll.${CULTIVAR}.vcf
UNIFIED_HOM_VCF=UnifiedHom.${CULTIVAR}.vcf

#=============================================================================== DATA
ORG_REF_FASTA=/lustre/home/ykawahar/TropicalJaponicaGenomeAnalysis/000_data/genome/IRGSP-1.0_genome_M_C_unanchored.fa
REF_FASTA_PREFIX=IRGSP-1.0_genome_M_C_unanchored

ORG_REF_GTF=/lustre/home/ykawahar/TropicalJaponicaGenomeAnalysis/000_data/annotation/RAP_MSU_w_cluster_for_Cufflinks.gtf
REF_GTF_PREFIX=RAP_MSU_w_cluster_for_Cufflinks

# path to tools
JRE=/lustre/bio/java/jre1.8.0_05/bin/java
FASTQC=$HOME/tool/NGS_preprocess/FastQC_v0.11.3/fastqc
TRIMMOMATIC=$HOME/tool/NGS_preprocess/Trimmomatic-0.33/trimmomatic-0.33.jar
ADAPTER=/lustre/bio/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa
BWA=$HOME/tool/Alignment/bwa-0.7.12/bwa
SAMTOOLS=$HOME/tool/Others/samtools-1.2/samtools
PICARD=$HOME/tool/Others/picard-tools-1.133/picard.jar
GATK=$HOME/tool/VariantCall/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar
SNPEFF_HOME=$HOME/tool/VariantCall/snpEff-v4.1
SCRIPT_HOME=/lustre/share/pft/script/GenomeAnalysis

#
#=============================================================================== for snpEff
SNPEFF_DB_BASENAME=RAP_MSU_on_IRGSP-1.0
SNPEFF_LOCAL_HOME=`pwd`/snpEff.${CULTIVAR}.$$
#
cp -rf $SNPEFF_HOME $SNPEFF_LOCAL_HOME
rm -rf $SNPEFF_LOCAL_HOME/data

SNPEFF_LOCAL_CONFIG=$SNPEFF_LOCAL_HOME/snpEff.config
config_tmp=config.${CULTIVAR}.$$

cat $SNPEFF_LOCAL_CONFIG | sed -e "s|data\.dir = ./data/|data\.dir = $SNPEFF_LOCAL_HOME/data/|" > $config_tmp
echo "${SNPEFF_DB_BASENAME}.genome : Oryza sativa Nipponbare IRGSP-1.0" >> $config_tmp
mv -f $config_tmp $SNPEFF_LOCAL_CONFIG
#
#=============================================================================== FastQC and Trimmomatic

if [ $PP_STEP -eq 0 ]; then
	mkdir -p $CULTIVAR
	cd $CULTIVAR
	WD=`pwd`;

	READ01=$READ_DIR/${CULTIVAR}_r1.fastq.bz2
	READ02=$READ_DIR/${CULTIVAR}_r2.fastq.bz2
	if [ ! -f $READ01 ]; then
		echo "[ERROR] $READ01 was not found."
		exit -1
	fi

	echo "Original reads:";
	echo "$READ01";
	echo "$READ02";

	READ1=${CULTIVAR}_r1.trimmed_pe.fastq.gz
	READ2=${CULTIVAR}_r2.trimmed_pe.fastq.gz

	mkdir FastQC_before_Trimmomatic
	echo "$FASTQC --nogroup --threads $NCORE --outdir FastQC_before_Trimmomatic --format fastq $READ01 $READ02"
	$FASTQC --nogroup --threads $NCORE --outdir FastQC_before_Trimmomatic --format fastq $READ01 $READ02

	echo "$JRE $JAVA_MEM -jar $TRIMMOMATIC PE -threads $NCORE -phred${INPUT_FASTQ_PHRED_OFFSET} $READ01 $READ02 $READ1 ${CULTIVAR}_r1.trimmed_se.fastq.gz $READ2 ${CULTIVAR}_r2.trimmed_se.fastq.gz ILLUMINACLIP:$ADAPTER:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 TOPHRED33"
	$JRE $JAVA_MEM -jar $TRIMMOMATIC PE -threads $NCORE -phred${INPUT_FASTQ_PHRED_OFFSET} $READ01 $READ02 $READ1 ${CULTIVAR}_r1.trimmed_se.fastq.gz $READ2 ${CULTIVAR}_r2.trimmed_se.fastq.gz ILLUMINACLIP:$ADAPTER:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 TOPHRED33

	mkdir FastQC_after_Trimmomatic
	echo "$FASTQC --nogroup --threads $NCORE --outdir FastQC_after_Trimmomatic --format fastq $READ1 $READ2"
	$FASTQC --nogroup --threads $NCORE --outdir FastQC_after_Trimmomatic --format fastq $READ1 $READ2

else	# use calculated data
	cd $CULTIVAR
	WD=`pwd`;

	READ1=${CULTIVAR}_r1.trimmed_pe.fastq.gz
	READ2=${CULTIVAR}_r2.trimmed_pe.fastq.gz
	if [ ! -f $READ01 ]; then
		echo "[ERROR] $READ01 was not found."
		exit -1
	fi
fi
#
echo "Preprocessed reads:";
echo "$READ1";
echo "$READ2";


#=============================================================================== prepare
if [ -f $STAT_OUTFILE ]; then
	rm -f $STAT_OUTFILE
fi
#
ln -s $ORG_REF_FASTA ${REF_FASTA_PREFIX}.${CULTIVAR}.00.fa
#
#=============================================================================== Loop start
NEXT_I=0
final_round=0
while :
do
	ROUND=$(printf "%02d" $NEXT_I)
	CORRECTED_REF_FASTA=${REF_FASTA_PREFIX}.${CULTIVAR}.${ROUND}.fa

	echo "Round $ROUND analysis based on $CORRECTED_REF_FASTA"
	date
#
	#=============================================================================== bwa, samtools and picard (1)
	# Add some indexing to reference files
	COM_BWA_INDEX="$BWA index $CORRECTED_REF_FASTA"
	echo $COM_BWA_INDEX; $COM_BWA_INDEX
	COM_SAMTOOLS_FAIDX="$SAMTOOLS faidx $CORRECTED_REF_FASTA"
	echo $COM_SAMTOOLS_FAIDX; $COM_SAMTOOLS_FAIDX
	COM_PICARD_CREATE_DICT="$JRE $JAVA_MEM -jar $PICARD CreateSequenceDictionary MAX_RECORDS_IN_RAM=10000000 REFERENCE=$CORRECTED_REF_FASTA OUTPUT=${REF_FASTA_PREFIX}.${CULTIVAR}.${ROUND}.dict"
	echo $COM_PICARD_CREATE_DICT; $COM_PICARD_CREATE_DICT

	# Alignment by BWA MEM
	RAW_SAM=raw_alignment.sam
	echo "$BWA mem -t $NCORE -M -T 30 -R \"@RG\tID:${CULTIVAR}\tSM:${CULTIVAR}\tPL:illumina\tLB:TruSeqV3\" $CORRECTED_REF_FASTA $READ1 $READ2 > $RAW_SAM"
	$BWA mem -t $NCORE -M -T 30 -R "@RG\tID:${CULTIVAR}\tSM:${CULTIVAR}\tPL:illumina\tLB:TruSeqV3" $CORRECTED_REF_FASTA $READ1 $READ2 > $RAW_SAM

	# convert SAM2BAM
	RAW_BAM=raw_alignment.${CULTIVAR}.${ROUND}.bam
	COM_CONVERT_SAM2BAM="$SAMTOOLS view -Sb -@ $NCORE -o $RAW_BAM $RAW_SAM"
	echo $COM_CONVERT_SAM2BAM; $COM_CONVERT_SAM2BAM

	TMP_SORT_BAM=raw_alignment.${CULTIVAR}.${ROUND}.sorted.tmp.bam
	COM_SORT_RAW_BAM="$SAMTOOLS sort -T ${CULTIVAR} -m $SAMTOOLS_MEM -@ $NCORE -o $TMP_SORT_BAM $RAW_BAM"
	echo $COM_SORT_RAW_BAM; $COM_SORT_RAW_BAM
	mv $TMP_SORT_BAM $RAW_BAM

	COM_INDEX_RAW_BAM="$SAMTOOLS index $RAW_BAM"
	echo $COM_INDEX_RAW_BAM; $COM_INDEX_RAW_BAM
	rm -f $RAW_SAM

	# extract high mapQ alignment
	HQ_BAM=hq_alignment.bam
	COM_EXTRACT_HQ_BAM="$SAMTOOLS view -b -@ $NCORE -q 20 -F 0x100 -o $HQ_BAM $RAW_BAM"
	echo $COM_EXTRACT_HQ_BAM; $COM_EXTRACT_HQ_BAM

	# nsort hq alignment
	NSORT_HQ_BAM=nsorted_hq_alignment.bam
	COM_NSORT_HQ_BAM="$SAMTOOLS sort -T ${CULTIVAR} -m $SAMTOOLS_MEM -@ $NCORE -n -o $NSORT_HQ_BAM $HQ_BAM"
	echo $COM_NSORT_HQ_BAM; $COM_NSORT_HQ_BAM
	rm -f $HQ_BAM

	# fixmate nsort hq alignment
	FIXMATE_NSORT_HQ_BAM=fixmate_nsorted_hq_alignment.bam
	COM_FIXMATE_NSORT_HQ_BAM="$SAMTOOLS fixmate -r -O bam $NSORT_HQ_BAM $FIXMATE_NSORT_HQ_BAM"
	echo $COM_FIXMATE_NSORT_HQ_BAM; $COM_FIXMATE_NSORT_HQ_BAM
	rm -f $NSORT_HQ_BAM

	# psort fixmate hq alignment
	FIXMATE_HQ_BAM=fixmate_hq_alignment.bam
	COM_PSORT_FIXMATE_HQ_BAM="$SAMTOOLS sort -T ${CULTIVAR} -m $SAMTOOLS_MEM -@ $NCORE -o $FIXMATE_HQ_BAM $FIXMATE_NSORT_HQ_BAM"
	echo $COM_PSORT_FIXMATE_HQ_BAM; $COM_PSORT_FIXMATE_HQ_BAM
	rm -f $FIXMATE_NSORT_HQ_BAM

	# extract properly mapped pe reads
	HQ_PE_BAM=hq_pe_alignment.bam
	COM_EXTRACT_HQ_PE_BAM="$SAMTOOLS view -b -@ $NCORE -F 0x10C -f 0x2 -o $HQ_PE_BAM $FIXMATE_HQ_BAM"
	echo $COM_EXTRACT_HQ_PE_BAM; $COM_EXTRACT_HQ_PE_BAM
	rm -f $FIXMATE_HQ_BAM

	# Mark duplicate reads (optical duplicates could bias variant detection by adding excessive coverage depth at a variant locus
	RMDUP_MATRIX=rmdup.matrix
	RMDUP_BAM=rmdup_alignment.bam
	COM_RMDUP="$JRE $JAVA_MEM -jar $PICARD MarkDuplicates MAX_RECORDS_IN_RAM=10000000 INPUT=$HQ_PE_BAM OUTPUT=$RMDUP_BAM METRICS_FILE=$RMDUP_MATRIX REMOVE_DUPLICATES=true"
	echo $COM_RMDUP; $COM_RMDUP
	COM_INDEX_RMDUP_BAM="$SAMTOOLS index $RMDUP_BAM"
	echo $COM_INDEX_RMDUP_BAM; $COM_INDEX_RMDUP_BAM
	rm -f $HQ_PE_BAM $RMDUP_MATRIX

	#=========== gatk
	##Improve mapping by local realignment (indels) 
	# Identify regions for indel local realignment of the selected chromosome
	TARGET_INTERVALS_LIST=target_intervals.list
	COM_GATK_TARGETCREATOR="$JRE $JAVA_MEM -jar $GATK -T RealignerTargetCreator --num_threads $NCORE_VC -I $RMDUP_BAM -R $CORRECTED_REF_FASTA -o $TARGET_INTERVALS_LIST"
	echo $COM_GATK_TARGETCREATOR; $COM_GATK_TARGETCREATOR
	
	# Perform indel local realignment of the selected chromosome
	REALN_BAM_PREFIX=realigned_alignment
	REALN_BAM=${REALN_BAM_PREFIX}.bam
	COM_GATK_REALIGN="$JRE $JAVA_MEM -jar $GATK -T IndelRealigner -I $RMDUP_BAM -R $CORRECTED_REF_FASTA -targetIntervals $TARGET_INTERVALS_LIST --consensusDeterminationModel USE_READS -o $REALN_BAM"
	echo $COM_GATK_REALIGN; $COM_GATK_REALIGN
	rm -f $RMDUP_BAM ${RMDUP_BAM}.bai $TARGET_INTERVALS_LIST 

	#=========== bwa, samtools and picard (2)
	# nsort realigned BAM
	NSORT_REALN_BAM=nsorted_realigned_alignment.bam
	COM_NSORT_REALN_BAM="$SAMTOOLS sort -T ${CULTIVAR} -m $SAMTOOLS_MEM -@ $NCORE -n -o $NSORT_REALN_BAM $REALN_BAM"
	echo $COM_NSORT_REALN_BAM; $COM_NSORT_REALN_BAM
	rm -f $REALN_BAM ${REALN_BAM_PREFIX}.bai 

	# fixmate nsort realigned BAM
	FIXMATE_NSORT_REALN_BAM=fixmate_nsorted_realigned_alignment.bam
	COM_FIXMATE_NSORT_REALN_BAM="$SAMTOOLS fixmate -r -O bam $NSORT_REALN_BAM $FIXMATE_NSORT_REALN_BAM"
	echo $COM_FIXMATE_NSORT_REALN_BAM; $COM_FIXMATE_NSORT_REALN_BAM
	rm -f $NSORT_REALN_BAM

	# psort fixmate realigned BAM
	FIXMATE_REALN_BAM=fixmate_realigned_alignment.bam
	COM_PSORT_FIXMATE_REALN_BAM="$SAMTOOLS sort -T ${CULTIVAR}.${ROUND} -m $SAMTOOLS_MEM -@ $NCORE -o $FIXMATE_REALN_BAM $FIXMATE_NSORT_REALN_BAM"
	echo $COM_PSORT_FIXMATE_REALN_BAM; $COM_PSORT_FIXMATE_REALN_BAM
	COM_INDEX_FIXMATE_REALN_BAM="$SAMTOOLS index $FIXMATE_REALN_BAM"
	echo $COM_INDEX_FIXMATE_REALN_BAM; $COM_INDEX_FIXMATE_REALN_BAM
	rm -f $FIXMATE_NSORT_REALN_BAM


	#========== HaplotypeCaller 1st
	HC_1st_VCF=variant_info_by_HC_1st.vcf
	FILTERED_HC_1st_VCF=variant_info_by_HC_1st.filtered.vcf

	# Call variants in your sequence data (diploid genome => HaplotypeCaller)
	COM_HC_1st="$JRE $JAVA_MEM -jar $GATK -T HaplotypeCaller --num_cpu_threads_per_data_thread $NCORE_VC --min_base_quality_score 20 -I $FIXMATE_REALN_BAM -R $CORRECTED_REF_FASTA -o $HC_1st_VCF"
	echo $COM_HC_1st; $COM_HC_1st

	# COM_SELECT_VAR_HC_1st="$JRE $JAVA_MEM -jar $GATK -T SelectVariants --num_threads $NCORE_VC --reference_sequence $CORRECTED_REF_FASTA --variant $HC_1st_VCF --out $FILTERED_HC_1st_VCF --selectexpressions \"QUAL >= 2000 && DP >= $DP_LOWER_LIMIT_HC && DP <= $DP_UPPER_LIMIT_HC && QD >= 30.0 && FS < 10.0 && MQ >= 40.0\""
	echo "$JRE $JAVA_MEM -jar $GATK -T SelectVariants --num_threads $NCORE_VC --reference_sequence $CORRECTED_REF_FASTA --variant $HC_1st_VCF --out $FILTERED_HC_1st_VCF --selectexpressions \"QUAL >= 2000 && DP >= $DP_LOWER_LIMIT_HC && DP <= $DP_UPPER_LIMIT_HC && QD >= 30.0 && FS < 10.0 && MQ >= 40.0\""
	$JRE $JAVA_MEM -jar $GATK -T SelectVariants --num_threads $NCORE_VC --reference_sequence $CORRECTED_REF_FASTA --variant $HC_1st_VCF --out $FILTERED_HC_1st_VCF --selectexpressions "QUAL >= 2000 && DP >= $DP_LOWER_LIMIT_HC && DP <= $DP_UPPER_LIMIT_HC && QD >= 30.0 && FS < 10.0 && MQ >= 40.0"
	# echo $COM_SELECT_VAR_HC_1st; $COM_SELECT_VAR_HC_1st
	rm -f $HC_1st_VCF ${HC_1st_VCF}.idx

	# Analyze patterns of covariation in the sequence dataset 
	BQSR_RECAL_TABLE=BQSR_recal_data.table
	COM_BQSR="$JRE $JAVA_MEM -jar $GATK -T BaseRecalibrator --num_cpu_threads_per_data_thread $NCORE_VC --reference_sequence $CORRECTED_REF_FASTA -I $FIXMATE_REALN_BAM -knownSites $FILTERED_HC_1st_VCF --out $BQSR_RECAL_TABLE"
	echo $COM_BQSR; $COM_BQSR
	rm -f $FILTERED_HC_1st_VCF ${FILTERED_HC_1st_VCF}.idx

	# Apply the recalibration to your sequence data 
	BQSR_BAM=final_alignment.${CULTIVAR}.${ROUND}.bam
	COM_PRINT_READS="$JRE $JAVA_MEM -jar $GATK -T PrintReads --num_cpu_threads_per_data_thread $NCORE_VC --reference_sequence $CORRECTED_REF_FASTA -I $FIXMATE_REALN_BAM -BQSR $BQSR_RECAL_TABLE --out $BQSR_BAM"
	echo $COM_PRINT_READS; $COM_PRINT_READS
	# COM_INDEX_BQSR_BAM="$SAMTOOLS index $BQSR_BAM"
	# echo $COM_INDEX_BQSR_BAM; $COM_INDEX_BQSR_BAM
	rm -f $FIXMATE_REALN_BAM ${FIXMATE_REALN_BAM}.bai $BQSR_RECAL_TABLE
	
	#========== HaplotypeCaller 2nd
	HC_2nd_VCF=variant_info_by_HC_2nd.vcf
	FILTERED_HC_2nd_VCF=variant_info_by_HC_2nd.filtered.vcf
	COM_HC_2nd="$JRE $JAVA_MEM -jar $GATK -T HaplotypeCaller --num_cpu_threads_per_data_thread $NCORE_VC --min_base_quality_score 20 -I $BQSR_BAM -R $CORRECTED_REF_FASTA -o $HC_2nd_VCF"
	echo $COM_HC_2nd; $COM_HC_2nd

	echo "$JRE $JAVA_MEM -jar $GATK -T VariantFiltration --reference_sequence $CORRECTED_REF_FASTA --variant $HC_2nd_VCF --out $FILTERED_HC_2nd_VCF --clusterSize 3 --clusterWindowSize 10 --filterExpression \"QUAL < 200 || DP < 3.0 || DP > $DP_UPPER_LIMIT || QD < 20.0 || FS >= 40.0 || MQ < 30.0\" --filterName \"FILTER\" --genotypeFilterExpression \"GQ < 20.0\" --genotypeFilterName \"FILTER\""
	$JRE $JAVA_MEM -jar $GATK -T VariantFiltration --reference_sequence $CORRECTED_REF_FASTA --variant $HC_2nd_VCF --out $FILTERED_HC_2nd_VCF --clusterSize 3 --clusterWindowSize 10 --filterExpression "QUAL < 200 || DP < 3.0 || DP > $DP_UPPER_LIMIT || QD < 20.0 || FS >= 40.0 || MQ < 30.0" --filterName "FILTER" --genotypeFilterExpression "GQ < 20.0" --genotypeFilterName "FILTER"
	# COM_VAR_FILTRATION="$JRE $JAVA_MEM -jar $GATK -T VariantFiltration --reference_sequence $CORRECTED_REF_FASTA --variant $HC_2nd_VCF --out $FILTERED_HC_2nd_VCF --clusterSize 3 --clusterWindowSize 10 --filterExpression \"QUAL < 200 || DP < 3.0 || DP > $DP_UPPER_LIMIT || QD < 20.0 || FS >= 40.0 || MQ < 30.0\" --filterName \"FILTER\" --genotypeFilterExpression \"GQ < 20.0\" --genotypeFilterName \"FILTER\""
	# echo $COM_VAR_FILTRATION; $COM_VAR_FILTRATION
	rm -f $HC_2nd_VCF ${HC_2nd_VCF}.idx

	#=============== Statistics and snpEff for 1st round
	if [ $NEXT_I -eq 0 ]; then
		COM_CALC_COV_DP_ROUND0="$SCRIPT_HOME/calc_coverage_depth_from_BAM.pl --ref_file $CORRECTED_REF_FASTA --bam_file $BQSR_BAM --outfile_base $OUTFILEBASE_1ST"
		echo $COM_CALC_COV_DP_ROUND0; $COM_CALC_COV_DP_ROUND0

		# COM_FLAGSTAT_ROUND0="$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_1ST}_flagstat.txt"
		# echo $COM_FLAGSTAT_ROUND0; $COM_FLAGSTAT_ROUND0
		echo "$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_1ST}_flagstat.txt"
		$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_1ST}_flagstat.txt

		## snpEff
		if [ -d $SNPEFF_LOCAL_HOME/data ]; then
			rm -rf $SNPEFF_LOCAL_HOME/data
		fi

		ORG_PROTEIN_FASTA=protein.${CULTIVAR}.${ROUND}.fa
		# COM_MAKE_AA_SEQ_ROUND0="$SCRIPT_HOME/residueFromGtf.pl $ORG_REF_FASTA $ORG_REF_GTF > $ORG_PROTEIN_FASTA"
		# echo $COM_MAKE_AA_SEQ_ROUND0; $COM_MAKE_AA_SEQ_ROUND0
		echo "$SCRIPT_HOME/residueFromGtf.pl $ORG_REF_FASTA $ORG_REF_GTF > $ORG_PROTEIN_FASTA"
		$SCRIPT_HOME/residueFromGtf.pl $ORG_REF_FASTA $ORG_REF_GTF > $ORG_PROTEIN_FASTA


		mkdir -p $SNPEFF_LOCAL_HOME/data/genomes
		mkdir $SNPEFF_LOCAL_HOME/data/$SNPEFF_DB_BASENAME
		ln -s $ORG_REF_FASTA $SNPEFF_LOCAL_HOME/data/genomes/${SNPEFF_DB_BASENAME}.fa
		ln -s $ORG_REF_GTF $SNPEFF_LOCAL_HOME/data/${SNPEFF_DB_BASENAME}/genes.gtf
		ln -s $WD/$ORG_PROTEIN_FASTA $SNPEFF_LOCAL_HOME/data/${SNPEFF_DB_BASENAME}/protein.fa

		COM_BUILD_SNPEFF_DB_ROUND0="$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar build -gtf22 -v $SNPEFF_DB_BASENAME"
		echo $COM_BUILD_SNPEFF_DB_ROUND0; $COM_BUILD_SNPEFF_DB_ROUND0

		FILTERED_HC_2nd_SNPEFF_VCF=variant_info_by_HC_2nd.filtered.snpEff.vcf
		# COM_SNPEFF_ROUND0="$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME $FILTERED_HC_2nd_VCF > $FILTERED_HC_2nd_SNPEFF_VCF"
		# echo $COM_SNPEFF_ROUND0; $COM_SNPEFF_ROUND0
		echo "$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME $FILTERED_HC_2nd_VCF > $FILTERED_HC_2nd_SNPEFF_VCF"
		$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME $FILTERED_HC_2nd_VCF > $FILTERED_HC_2nd_SNPEFF_VCF

		mv -f snpEff_genes.txt snpEff_genes_1st.${CULTIVAR}.txt
		mv -f snpEff_summary.html snpEff_summary_1st.${CULTIVAR}.html

		ALL_VCF_WO_SNPEFF=variant_info_All.${CULTIVAR}.${ROUND}.vcf
		# COM_EXCLUDE_FILTERED_RECORD_FROM_FILTERED_HC_2nd_VCF="$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_VCF > $ALL_VCF_WO_SNPEFF"
		# echo $COM_EXCLUDE_FILTERED_RECORD_FROM_FILTERED_HC_2nd_VCF; $COM_EXCLUDE_FILTERED_RECORD_FROM_FILTERED_HC_2nd_VCF
		echo "$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_VCF > $ALL_VCF_WO_SNPEFF"
		$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_VCF > $ALL_VCF_WO_SNPEFF

		cp -f $ALL_VCF_WO_SNPEFF $UNIFIED_ALL_VCF

		ALL_VCF=variant_info_All.snpEff.${CULTIVAR}.${ROUND}.vcf
		# COM_EXCLUDE_FILTERED_RECORD="$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_SNPEFF_VCF > $ALL_VCF"
		# echo $COM_EXCLUDE_FILTERED_RECORD; $COM_EXCLUDE_FILTERED_RECORD
		echo "$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_SNPEFF_VCF > $ALL_VCF"
		$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_SNPEFF_VCF > $ALL_VCF
		rm -f $FILTERED_HC_2nd_SNPEFF_VCF ${FILTERED_HC_2nd_SNPEFF_VCF}.idx

		#===== make RefConf VCF
		# REFCONF_VCF=variant_info_by_HC_2nd.RefConf.vcf
		# FILTERED_REFCONF_VCF=variant_info_by_HC_2nd.filtered.RefConf.vcf
		# COM_HC_REFCONF="$JRE $JAVA_MEM -jar $GATK -T HaplotypeCaller --num_cpu_threads_per_data_thread $NCORE_VC --emitRefConfidence BP_RESOLUTION --min_base_quality_score 20 -I $BQSR_BAM -R $CORRECTED_REF_FASTA -o $REFCONF_VCF"
		# echo $COM_HC_REFCONF; $COM_HC_REFCONF

		# echo "$JRE $JAVA_MEM -jar $GATK -T VariantFiltration --reference_sequence $CORRECTED_REF_FASTA --variant $REFCONF_VCF --out $FILTERED_REFCONF_VCF --genotypeFilterExpression \"DP < 3.0 || DP > $DP_UPPER_LIMIT || GQ < 20.0\" --genotypeFilterName \"FILTER\""
		# $JRE $JAVA_MEM -jar $GATK -T VariantFiltration --reference_sequence $CORRECTED_REF_FASTA --variant $REFCONF_VCF --out $FILTERED_REFCONF_VCF --genotypeFilterExpression "DP < 3.0 || DP > $DP_UPPER_LIMIT || GQ < 20.0" --genotypeFilterName "FILTER"
		### COM_VAR_FILTRATION_REFCONF="$JRE $JAVA_MEM -jar $GATK -T VariantFiltration --reference_sequence $CORRECTED_REF_FASTA --variant $REFCONF_VCF --out $FILTERED_REFCONF_VCF --genotypeFilterExpression \"DP < 3.0 || DP > $DP_UPPER_LIMIT || GQ < 20.0\" --genotypeFilterName \"FILTER\""
		### echo $COM_VAR_FILTRATION_REFCONF; $COM_VAR_FILTRATION_REFCONF
		# rm -f $REFCONF_VCF ${REFCONF_VCF}.idx
		#=====

		HOM_VCF_WO_SNPEFF=variant_info_Hom.${CULTIVAR}.${ROUND}.vcf
		echo "$JRE $JAVA_MEM -jar $GATK -T SelectVariants --reference_sequence $CORRECTED_REF_FASTA --variant $ALL_VCF_WO_SNPEFF --out $HOM_VCF_WO_SNPEFF --excludeFiltered -select \"vc.getGenotype(\'$CULTIVAR\').isHomVar()\""
		$JRE $JAVA_MEM -jar $GATK -T SelectVariants --reference_sequence $CORRECTED_REF_FASTA --variant $ALL_VCF_WO_SNPEFF --out $HOM_VCF_WO_SNPEFF --excludeFiltered -select "vc.getGenotype('$CULTIVAR').isHomVar()"
		# COM_SELECT_HOM_VAR_FROM_VCF_WO_SNPEFF="$JRE $JAVA_MEM -jar $GATK -T SelectVariants --reference_sequence $CORRECTED_REF_FASTA --variant $ALL_VCF_WO_SNPEFF --out $HOM_VCF_WO_SNPEFF --excludeFiltered -select \"vc.getGenotype(\'$CULTIVAR\').isHomVar()\""
		# echo $COM_SELECT_HOM_VAR_FROM_VCF_WO_SNPEFF; $COM_SELECT_HOM_VAR_FROM_VCF_WO_SNPEFF

		HOM_VCF=variant_info_Hom.snpEff.${CULTIVAR}.${ROUND}.vcf

	else
		ALL_VCF=variant_info_All.${CULTIVAR}.${ROUND}.vcf
		# COM_EXCLUDE_FILTERED_RECORD="$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_VCF > $ALL_VCF"
		# echo $COM_EXCLUDE_FILTERED_RECORD; $COM_EXCLUDE_FILTERED_RECORD
		echo "$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_VCF > $ALL_VCF"
		$SCRIPT_HOME/exclude_FILTER_record_from_VCF.pl $FILTERED_HC_2nd_VCF > $ALL_VCF
		rm -f $FILTERED_HC_2nd_VCF ${FILTERED_HC_2nd_VCF}.idx

		HOM_VCF=variant_info_Hom.${CULTIVAR}.${ROUND}.vcf
	fi

	rm -f $FILTERED_HC_2nd_VCF ${FILTERED_HC_2nd_VCF}.idx

	echo "$JRE $JAVA_MEM -jar $GATK -T SelectVariants --reference_sequence $CORRECTED_REF_FASTA --variant $ALL_VCF --out $HOM_VCF --excludeFiltered -select \"vc.getGenotype(\'$CULTIVAR\').isHomVar()\""
	$JRE $JAVA_MEM -jar $GATK -T SelectVariants --reference_sequence $CORRECTED_REF_FASTA --variant $ALL_VCF --out $HOM_VCF --excludeFiltered -select "vc.getGenotype('$CULTIVAR').isHomVar()"
	# COM_SELECT_HOM_VAR="$JRE $JAVA_MEM -jar $GATK -T SelectVariants --reference_sequence $CORRECTED_REF_FASTA --variant $ALL_VCF --out $HOM_VCF --excludeFiltered -select \"vc.getGenotype(\'$CULTIVAR\').isHomVar()\""
	# echo $COM_SELECT_HOM_VAR; $COM_SELECT_HOM_VAR

	if [ $NEXT_I -eq 0 ]; then
		cp -f $HOM_VCF_WO_SNPEFF $UNIFIED_HOM_VCF
	fi

#=================== union VCF file for after Round 1
	if [ $NEXT_I -ne 0 ]; then
		mv -f $UNIFIED_ALL_VCF ${UNIFIED_ALL_VCF}.tmp
		mv -f $UNIFIED_HOM_VCF ${UNIFIED_HOM_VCF}.tmp

		# COM_UNION_ALL_VCF="$SCRIPT_HOME/unionVcf.pl ${UNIFIED_ALL_VCF}.tmp $ALL_VCF > $UNIFIED_ALL_VCF"
		# echo $COM_UNION_ALL_VCF; $COM_UNION_ALL_VCF
		echo "$SCRIPT_HOME/unionVcf.pl ${UNIFIED_ALL_VCF}.tmp $ALL_VCF > $UNIFIED_ALL_VCF"
		$SCRIPT_HOME/unionVcf.pl ${UNIFIED_ALL_VCF}.tmp $ALL_VCF > $UNIFIED_ALL_VCF
		# COM_UNION_HOM_VCF="$SCRIPT_HOME/unionVcf.pl ${UNIFIED_HOM_VCF}.tmp $HOM_VCF > $UNIFIED_HOM_VCF"
		# echo $COM_UNION_HOM_VCF; $COM_UNION_HOM_VCF
		echo "$SCRIPT_HOME/unionVcf.pl ${UNIFIED_HOM_VCF}.tmp $HOM_VCF > $UNIFIED_HOM_VCF"
		$SCRIPT_HOME/unionVcf.pl ${UNIFIED_HOM_VCF}.tmp $HOM_VCF > $UNIFIED_HOM_VCF

		rm -f ${UNIFIED_ALL_VCF}.tmp
		rm -f ${UNIFIED_HOM_VCF}.tmp

		if [ $NEXT_I -ne 1 ]; then
			# rm -f $BQSR_BAM ${BQSR_BAM}.bai
			# rm -f $RAW_BAM ${RAW_BAM}.bai
			rm -f ${BQSR_BAM}.bai
			rm -f ${RAW_BAM}.bai
		fi
	fi

#=============================================================================== Statistics
	HOM_SNP_COUNT=`grep -v ^# $HOM_VCF | awk '{print $4$5}' | grep '^.\{2\}$' | wc -l | awk '{print $1}'`
	HOM_INDEL_COUNT=`grep -v ^# $HOM_VCF | awk '{print $4$5}' | grep '.\{3\}' | wc -l | awk '{print $1}'`
	echo -e "${CULTIVAR}\t${ROUND}\t${HOM_VCF}\t${HOM_SNP_COUNT}\t${HOM_INDEL_COUNT}" >> $STAT_OUTFILE
	final_round=$NEXT_I

#=============================================================================== snpEff (last)
	HOM_VAR_COUNT=`grep -v ^# $HOM_VCF | wc -l | awk '{print $1}'`
	# Meet the finishing criteria
	TMP_MAX_LOOP=`expr $MAX_LOOP - 1`
	if [ $NEXT_I -ge $TMP_MAX_LOOP ]; then
		cp $UNIFIED_ALL_VCF ${UNIFIED_ALL_VCF}.wo_snpEff
		cp $UNIFIED_HOM_VCF ${UNIFIED_HOM_VCF}.wo_snpEff
		
		mv -f $UNIFIED_ALL_VCF ${UNIFIED_ALL_VCF}.tmp
		mv -f $UNIFIED_HOM_VCF ${UNIFIED_HOM_VCF}.tmp

		# COM_SNPEFF_ALL_ROUND_FIN="$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_ALL_VCF}.tmp > $UNIFIED_ALL_VCF"
		# echo $COM_SNPEFF_ALL_ROUND_FIN; $COM_SNPEFF_ALL_ROUND_FIN
		echo "$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_ALL_VCF}.tmp > $UNIFIED_ALL_VCF"
		$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_ALL_VCF}.tmp > $UNIFIED_ALL_VCF
		mv -f snpEff_genes.txt snpEff_genes_UnifiedAll.${CULTIVAR}.txt
		mv -f snpEff_summary.html snpEff_summary_UnifiedAll.${CULTIVAR}.html

		# COM_SNPEFF_HOM_ROUND_FIN="$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_HOM_VCF}.tmp > $UNIFIED_HOM_VCF"
		# echo $COM_SNPEFF_HOM_ROUND_FIN; $COM_SNPEFF_HOM_ROUND_FIN
		echo "$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_HOM_VCF}.tmp > $UNIFIED_HOM_VCF"
		$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_HOM_VCF}.tmp > $UNIFIED_HOM_VCF
		mv -f snpEff_genes.txt snpEff_genes_UnifiedHom.${CULTIVAR}.txt
		mv -f snpEff_summary.html snpEff_summary_UnifiedHom.${CULTIVAR}.html

		rm -f ${UNIFIED_ALL_VCF}.tmp
		rm -f ${UNIFIED_HOM_VCF}.tmp

		NEXT_I=`expr $NEXT_I + 1`
		NEXT_ROUND=$(printf "%02d" $NEXT_I)

		NEXT_CORRECTED_GTF=${REF_GTF_PREFIX}.${CULTIVAR}.${NEXT_ROUND}.gtf
		# COM_MAKE_FINAL_GTF="$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF"
		# echo $COM_MAKE_FINAL_GTF; $COM_MAKE_FINAL_GTF
		echo "$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF"
		$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF

		NEXT_CORRECTED_REF_FASTA=${REF_FASTA_PREFIX}.${CULTIVAR}.${NEXT_ROUND}.fa
		#COM_MAKE_FINAL_FASTA="$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA"
		# echo $COM_MAKE_FINAL_FASTA; $COM_MAKE_FINAL_FASTA
		echo "$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA"
		$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA

		NEXT_CORRECTED_PROTEIN_FASTA=protein.${CULTIVAR}.${NEXT_ROUND}.fa 
		# COM_MAKE_FINAL_AA_SEQ="$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA"
		# echo $COM_MAKE_FINAL_AA_SEQ; $COM_MAKE_FINAL_AA_SEQ
		echo "$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA"
		$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA

		COM_CALC_COV_DP_ROUND_FIN="$SCRIPT_HOME/calc_coverage_depth_from_BAM.pl --ref_file $CORRECTED_REF_FASTA --bam_file $BQSR_BAM --outfile_base $OUTFILEBASE_FIN"
		echo $COM_CALC_COV_DP_ROUND_FIN; $COM_CALC_COV_DP_ROUND_FIN

		# COM_FLAGSTAT_ROUND_FIN="$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_FIN}_flagstat.txt"
		# echo $COM_FLAGSTAT_ROUND_FIN; $COM_FLAGSTAT_ROUND_FIN
		echo "$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_FIN}_flagstat.txt"
		$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_FIN}_flagstat.txt

		rm -rf $SNPEFF_LOCAL_HOME
	break
	fi

	# Meet the finishing criteria
	if [ $HOM_VAR_COUNT -le $ALLOW_VAR_COUNT ]; then
		cp $UNIFIED_ALL_VCF ${UNIFIED_ALL_VCF}.wo_snpEff
		cp $UNIFIED_HOM_VCF ${UNIFIED_HOM_VCF}.wo_snpEff

		mv -f $UNIFIED_ALL_VCF ${UNIFIED_ALL_VCF}.tmp
		mv -f $UNIFIED_HOM_VCF ${UNIFIED_HOM_VCF}.tmp

		# COM_SNPEFF_ALL_ROUND_FIN="$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_ALL_VCF}.tmp > $UNIFIED_ALL_VCF"
		# echo $COM_SNPEFF_ALL_ROUND_FIN; $COM_SNPEFF_ALL_ROUND_FIN
		echo "$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_ALL_VCF}.tmp > $UNIFIED_ALL_VCF"
		$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_ALL_VCF}.tmp > $UNIFIED_ALL_VCF
		mv -f snpEff_genes.txt snpEff_genes_UnifiedAll.${CULTIVAR}.txt
		mv -f snpEff_summary.html snpEff_summary_UnifiedAll.${CULTIVAR}.html

		# COM_SNPEFF_HOM_ROUND_FIN="$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_HOM_VCF}.tmp > $UNIFIED_HOM_VCF"
		# echo $COM_SNPEFF_HOM_ROUND_FIN; $COM_SNPEFF_HOM_ROUND_FIN
		echo "$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_HOM_VCF}.tmp > $UNIFIED_HOM_VCF"
		$JRE $JAVA_MEM -jar $SNPEFF_LOCAL_HOME/snpEff.jar -v $SNPEFF_DB_BASENAME ${UNIFIED_HOM_VCF}.tmp > $UNIFIED_HOM_VCF
		mv -f snpEff_genes.txt snpEff_genes_UnifiedHom.${CULTIVAR}.txt
		mv -f snpEff_summary.html snpEff_summary_UnifiedHom.${CULTIVAR}.html

		rm -f ${UNIFIED_ALL_VCF}.tmp
		rm -f ${UNIFIED_HOM_VCF}.tmp

		NEXT_I=`expr $NEXT_I + 1`
		NEXT_ROUND=$(printf "%02d" $NEXT_I)

		NEXT_CORRECTED_GTF=${REF_GTF_PREFIX}.${CULTIVAR}.${NEXT_ROUND}.gtf
		# COM_MAKE_FINAL_GTF="$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF"
		# echo $COM_MAKE_FINAL_GTF; $COM_MAKE_FINAL_GTF
		echo "$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF"
		$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF

		NEXT_CORRECTED_REF_FASTA=${REF_FASTA_PREFIX}.${CULTIVAR}.${NEXT_ROUND}.fa
		# COM_MAKE_FINAL_FASTA="$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA"
		# echo $COM_MAKE_FINAL_FASTA; $COM_MAKE_FINAL_FASTA
		echo "$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA"
		$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA

		NEXT_CORRECTED_PROTEIN_FASTA=protein.${CULTIVAR}.${NEXT_ROUND}.fa 
		# COM_MAKE_FINAL_AA_SEQ="$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA"
		# echo $COM_MAKE_FINAL_AA_SEQ; $COM_MAKE_FINAL_AA_SEQ
		echo "$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA"
		$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA

		COM_CALC_COV_DP_ROUND_FIN="$SCRIPT_HOME/calc_coverage_depth_from_BAM.pl --ref_file $CORRECTED_REF_FASTA --bam_file $BQSR_BAM --outfile_base $OUTFILEBASE_FIN"
		echo $COM_CALC_COV_DP_ROUND_FIN; $COM_CALC_COV_DP_ROUND_FIN

		# COM_FLAGSTAT_ROUND_FIN="$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_FIN}_flagstat.txt"
		# echo $COM_FLAGSTAT_ROUND_FIN; $COM_FLAGSTAT_ROUND_FIN
		echo "$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_FIN}_flagstat.txt"
		$SAMTOOLS flagstat $BQSR_BAM > ${OUTFILEBASE_FIN}_flagstat.txt

		rm -rf $SNPEFF_LOCAL_HOME
	break
	fi


#=================== update data
	NEXT_I=`expr $NEXT_I + 1`
	NEXT_ROUND=$(printf "%02d" $NEXT_I)

	NEXT_CORRECTED_REF_FASTA=${REF_FASTA_PREFIX}.${CULTIVAR}.${NEXT_ROUND}.fa
	# COM_MAKE_NEXT_REF_FASTA="$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA"
	# echo $COM_MAKE_NEXT_REF_FASTA; $COM_MAKE_NEXT_REF_FASTA
	echo "$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA"
	$SCRIPT_HOME/updateFasta.pl $CORRECTED_REF_FASTA $HOM_VCF > $NEXT_CORRECTED_REF_FASTA

	if [ $NEXT_I -ne 2 ]; then
		rm -f ${CORRECTED_REF_FASTA}* ${REF_FASTA_PREFIX}.${CULTIVAR}.${ROUND}.dict
	fi

	if [ $NEXT_I -eq 1 ]; then
		NEXT_CORRECTED_GTF=${REF_GTF_PREFIX}.${CULTIVAR}.${NEXT_ROUND}.gtf
		# COM_MAKE_FINAL_GTF="$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF"
		# echo $COM_MAKE_FINAL_GTF; $COM_MAKE_FINAL_GTF
		echo "$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF"
		$SCRIPT_HOME/updateGtf.pl $ORG_REF_GTF $UNIFIED_HOM_VCF > $NEXT_CORRECTED_GTF

		NEXT_CORRECTED_PROTEIN_FASTA=protein.${CULTIVAR}.${NEXT_ROUND}.fa 
		# COM_MAKE_FINAL_AA_SEQ="$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA"
		# echo $COM_MAKE_FINAL_AA_SEQ; $COM_MAKE_FINAL_AA_SEQ
		echo "$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA"
		$SCRIPT_HOME/residueFromGtf.pl $NEXT_CORRECTED_REF_FASTA $NEXT_CORRECTED_GTF > $NEXT_CORRECTED_PROTEIN_FASTA
	fi
done

#=============================================================================== Statistics final
REPLACED_COUNT=`expr $final_round + 1`
echo "${REPLACED_COUNT} rounds of variant replacement were performed." >> $STAT_OUTFILE

#
#=============================================================================== file delete
##### rm -f *.bai
##### rm -f *.idx
##### rm -f *.dict

##### rm -f *.amb
##### rm -f *.ann
##### rm -f *.bwt
##### rm -f *.fai
##### rm -f *.pac
##### rm -f *.sa

echo "Completed!!"
date
#
