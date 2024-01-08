#!/bin/bash

## Setting input parameters from cenotetaker3.py

CENOTE_SCRIPTS=$1
original_contigs=$2
run_title=$3
PROPHAGE=$4
CPU=$5
VERSION=$6
ANNOTATION_MODE=$7
TEMPLATE_FILE=$8
READS=$9
circ_length_cutoff=${10}
linear_length_cutoff=${11}
CIRC_MINIMUM_DOMAINS=${12}
LIN_MINIMUM_DOMAINS=${13}
HALL_TYPE=${14}
C_DBS=${15}
HMM_DBS=${16}
WRAP=${17}
CALLER=${18}
HHSUITE_TOOL=${19} ##
ISO_SOURCE=${20}
COLLECT_DATE=${21}
META_TYPE=${22}
SRR=${23}
SRX=${24}
BIOSAMP=${25}
PRJ=${26}
ASSEMBLER=${27}
MOL_TYPE=${28}
DATA_SOURCE=${29}


### check HHsuite DBs
if [ -s ${C_DBS}/hhsearch_DBs/NCBI_CD/NCBI_CD_a3m.ffdata ] ; then
	CD_HHSUITE="${C_DBS}/hhsearch_DBs/NCBI_CD/NCBI_CD"
else
	CD_HHSUITE=""
fi
if [ -s ${C_DBS}/hhsearch_DBs/pfam_32_db/pfam_a3m.ffdata ] ; then
	PFAM_HHSUITE="${C_DBS}/hhsearch_DBs/pfam_32_db/pfam"
else
	PFAM_HHSUITE=""
fi
if [ -s ${C_DBS}/hhsearch_DBs/pdb70/pdb70_a3m.ffdata ] ; then
	PDB_HHSUITE="${C_DBS}/hhsearch_DBs/pdb70/pdb70"
else
	PDB_HHSUITE=""
fi

HHSUITE_DB_STR=""
if [ -n "$CD_HHSUITE" ] ; then
	HHSUITE_DB_STR="${HHSUITE_DB_STR}-d ${CD_HHSUITE} "
fi
if [ -n "$PFAM_HHSUITE" ] ; then
	HHSUITE_DB_STR="${HHSUITE_DB_STR}-d ${PFAM_HHSUITE} "
fi
if [ -n "$PDB_HHSUITE" ] ; then
	HHSUITE_DB_STR="${HHSUITE_DB_STR}-d ${PDB_HHSUITE}"
fi
if [ ! -n "$HHSUITE_DB_STR" ] ; then
	echo "HHsuite databases not found at ${C_DBS}/hhsearch_DBs"
	echo "$HHSUITE_TOOL will not be run"
	HHSUITE_TOOL="none"
fi


### echo colors
Color_Off='\033[0m'       # Text Reset
# Bold
BBlack='\033[1;30m'       # Black
BRed='\033[1;31m'         # Red
BGreen='\033[1;32m'       # Green
BYellow='\033[1;33m'      # Yellow
BBlue='\033[1;34m'        # Blue
BPurple='\033[1;35m'      # Purple
BCyan='\033[1;36m'        # Cyan
BWhite='\033[1;37m'       # White


### if read provided, check that each file exists

if [ "${READS}" != "none" ] ; then
	for READ_FILE in $READS ; do
		if [ -s $READ_FILE ] ; then
			echo $READ_FILE
		else
			echo "$READ_FILE not found."
			echo "exiting"
			exit
		fi
	done
fi

MDYT=$( date +"%m-%d-%y---%T" )
echo -e "${BBlack}time update: configuring run directory  ${MDYT}${Color_Off}"


if [ ! -d "${run_title}/ct_processing" ]; then
	mkdir ${run_title}/ct_processing
fi

TEMP_DIR="${run_title}/ct_processing"

### setting universal minimum length
if [ $circ_length_cutoff -gt $linear_length_cutoff ] ; then
	LENGTH_MINIMUM=$linear_length_cutoff
else
	LENGTH_MINIMUM=$circ_length_cutoff

fi


### save all the parameters in the run_arguments.txt file then print to terminal
echo "@@@@@@@@@@@@@@@@@@@@@@@@@" >> ${run_title}/run_arguments.txt
echo "Your specified arguments:" >> ${run_title}/run_arguments.txt
echo "Cenote-Taker version:              $VERSION" >> ${run_title}/run_arguments.txt
echo "original contigs:                  $original_contigs" >> ${run_title}/run_arguments.txt
echo "title of this run:                 $run_title" >> ${run_title}/run_arguments.txt
echo "Prune prophages?                   $PROPHAGE" >> ${run_title}/run_arguments.txt
echo "CPUs used for run:                 $CPU" >> ${run_title}/run_arguments.txt
echo "Annotation only?                   $ANNOTATION_MODE" >> ${run_title}/run_arguments.txt
echo "minimum circular contig length:    $circ_length_cutoff" >> ${run_title}/run_arguments.txt
echo "minimum linear contig length:      $linear_length_cutoff" >> ${run_title}/run_arguments.txt
echo "virus hallmark type(s) to count:   $HALL_TYPE" >> ${run_title}/run_arguments.txt
echo "min. viral hallmarks for linear:   $LIN_MINIMUM_DOMAINS" >> ${run_title}/run_arguments.txt
echo "min. viral hallmarks for circular: $CIRC_MINIMUM_DOMAINS" >> ${run_title}/run_arguments.txt
echo "Wrap contigs?                      $WRAP" >> ${run_title}/run_arguments.txt
echo "HMM db version                     $HMM_DBS" >> ${run_title}/run_arguments.txt
echo "ORF Caller:                        $CALLER" >> ${run_title}/run_arguments.txt
echo "Cenote DBs directory:              $C_DBS" >> ${run_title}/run_arguments.txt
echo "Cenote scripts directory:          $CENOTE_SCRIPTS" >> ${run_title}/run_arguments.txt
echo "Template file:                     $TEMPLATE_FILE" >> ${run_title}/run_arguments.txt
echo "read file(s):                      $READS" >> ${run_title}/run_arguments.txt
echo "HHsuite tool:                      $HHSUITE_TOOL" >> ${run_title}/run_arguments.txt

cat ${run_title}/run_arguments.txt

echo " "


### filtering input contigs by minimum length and renaming for cenote-taker
#- input:
#-- ${original_contigs}
#- output:
#-- ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta

if [ -s ${original_contigs} ] ; then
	
	seqkit seq --quiet -g -m $LENGTH_MINIMUM $original_contigs |\
	  seqkit replace --quiet -p '^' -r ${run_title}_{nr}@#@# |\
	  sed 's/@#@#/ /g' > ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta

else
	echo "${original_contigs} not found"
	exit
fi


### split contigs into equal parts for prodigal ORF calling, call ORFs
#- input:
#-- ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta
#- output:
#-- ${TEMP_DIR}/contig_name_map.tsv
#-- ${TEMP_DIR}/ORF_orig_contigs/split/*.faa (1 or more)

if [ -s ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta ] ; then
	if [ ! -d "${TEMP_DIR}/ORF_orig_contigs" ]; then
		mkdir ${TEMP_DIR}/ORF_orig_contigs
	fi

	# table with ct name and input name in separate columns
	TABQ=$'\t'
	grep -F ">" ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta |\
	  sed "s/ /\t/" | sed 's/>//g' > ${TEMP_DIR}/contig_name_map.tsv

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BRed}time update: running pyrodigal on all contigs  ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/pyrodigal_gv_runner.py ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta\
	  ${TEMP_DIR}/ORF_orig_contigs $CPU $CALLER

	seqkit split --quiet -j $CPU -p $CPU -O ${TEMP_DIR}/ORF_orig_contigs/split ${TEMP_DIR}/ORF_orig_contigs/pyrodigal_gv_AAs.prod.faa

else
	echo "couldn't find ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta"
	echo "exiting"
	exit
fi

### set minimum hallmark genes
if [ $CIRC_MINIMUM_DOMAINS -gt $LIN_MINIMUM_DOMAINS ] ; then
	HALLMARK_MINIMUM=$LIN_MINIMUM_DOMAINS
else
	HALLMARK_MINIMUM=$CIRC_MINIMUM_DOMAINS
fi

### run pyhmmer on prodigal ORF files, virion DB and rep DB
#- input: -#
#-- ${TEMP_DIR}/ORF_orig_contigs/split/*.faa (1 or more)
#- output: -#
#-- ${TEMP_DIR}/orig_pyhmmer_virion/contig_hit_count.tsv
#---- 	fields
#---- 	(contig	count)
#-- ${TEMP_DIR}/orig_pyhmmer_virion/pyhmmer_report_AAs.tsv
#--- 	fields
#--- 	(ORFquery	contig	target	evalue	pvalue)
#-- ${TEMP_DIR}/orig_pyhmmer_rep/contig_hit_count.tsv
#---- 	fields
#---- 	(contig	count)
#-- ${TEMP_DIR}/orig_pyhmmer_rep/pyhmmer_report_AAs.tsv
#----	fields
#----	(ORFquery	contig	target	evalue	pvalue)
#-- ${TEMP_DIR}/contigs_to_keep.txt
#-- ${TEMP_DIR}/hallmarks_per_orig_contigs.tsv
#-- ${TEMP_DIR}/hallmarks_for_keepcontigs1.txt

SPLIT_ORIG_AAs=$( find ${TEMP_DIR}/ORF_orig_contigs/split -type f -name "*.faa" )

if [ -n "$SPLIT_ORIG_AAs" ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BRed}time update: running pyhmmer on all ORFs  ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/ORF_orig_contigs/split ${TEMP_DIR}/orig_pyhmmer_virion\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/Virion_HMMs.h3m $CPU 1e-7 0.1

	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/ORF_orig_contigs/split ${TEMP_DIR}/orig_pyhmmer_rep\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/DNA_rep_HMMs.h3m $CPU 1e-7 0.1

	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/ORF_orig_contigs/split ${TEMP_DIR}/orig_pyhmmer_rdrp\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/RDRP_HMMs.h3m $CPU 1e-7 0.8

	python ${CENOTE_SCRIPTS}/python_modules/combine_hallmark_counts.py ${TEMP_DIR}/orig_pyhmmer_virion\
	  ${TEMP_DIR}/orig_pyhmmer_rep ${TEMP_DIR}/orig_pyhmmer_rdrp ${HALLMARK_MINIMUM} "${HALL_TYPE}" ${TEMP_DIR} ${TEMP_DIR}/contig_name_map.tsv

else
	echo "couldn't find prodigal AA seqs in ${TEMP_DIR}/ORF_orig_contigs/split"
fi


### grabbing contigs with minimum marker gene number
#- input: -#
#-- ${TEMP_DIR}/contigs_to_keep.txt
#-- ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta
#- output: -#
#-- ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta

if [ -s ${TEMP_DIR}/contigs_to_keep.txt ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BYellow}grabbing contigs with minimum marker gene number ${MDYT}${Color_Off}"

	KEEP_COUNT=$( cat ${TEMP_DIR}/contigs_to_keep.txt | wc -l )

	if [ $KEEP_COUNT -ge 1 ] ; then

		seqkit grep --quiet -f ${TEMP_DIR}/contigs_to_keep.txt\
		  ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta > ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta

	else
		echo "no contigs with ${KEEP_COUNT} hallmark gene(s) were detected. Exiting."
		exit
	fi
else
	echo "couldn't find ${TEMP_DIR}/contigs_to_keep.txt. Exiting."
	exit
fi

### detecting DTRs and ITRs. Trimming DTRs
#- input: -#
#-- ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta
#-- ${TEMP_DIR}/hallmarks_per_orig_contigs.tsv
#-- variables: $circ_length_cutoff $linear_length_cutoff $CIRC_MINIMUM_DOMAINS $WRAP (boolean)
#- output: -#
#-- ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta
#-- ${TEMP_DIR}/hallmark_contigs_terminal_repeat_summary.tsv
#-- ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv
#----	fields
#----	(contig	in_length_contig	out_length_contig	dtr_seq	itr_seq)
#-- ${TEMP_DIR}/contigs_over_threshold.txt
#-- ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta

if [ -s ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BYellow}detecting DTRs and ITRs ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/terminal_repeats.py ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta\
	${TEMP_DIR}/hallmarks_per_orig_contigs.tsv $circ_length_cutoff $linear_length_cutoff $CIRC_MINIMUM_DOMAINS\
	$LIN_MINIMUM_DOMAINS ${TEMP_DIR} $WRAP

	if [ -s ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta ] && [ -s ${TEMP_DIR}/contigs_over_threshold.txt ] ; then
		
		seqkit grep --quiet -f ${TEMP_DIR}/contigs_over_threshold.txt\
		  ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta > ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta

	else
		echo "couldn't find contigs_over_threshold.txt. No virus contigs over set thresholds."
	fi

else
	echo "couldn't find hallmark contigs at ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta"
	exit
fi

### Rotating DTR contigs
#- input: -#
#-- ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta
#-- ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv
#- output: -#
#-- ${TEMP_DIR}/oriented_hallmark_contigs.fasta

if [ -s ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta ] && [ -s ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv ] ; then
	if [ "$WRAP" == "True" ] ; then

		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BYellow}rotating DTR contigs ${MDYT}${Color_Off}"

		if [ ! -d "${TEMP_DIR}/rotation" ]; then
			mkdir ${TEMP_DIR}/rotation
		fi

		awk '{OFS=FS="\t"}{ if (NR != 1 && $4 != "NA") {print $1, $3} }' ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv |\
		  while read DTR_CONTIG TRIM_LENGTH ; do
			seqkit grep --quiet -p "${DTR_CONTIG}" ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta |\
			  prodigal -a ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa -i /dev/stdin -c -p meta -q >/dev/null 2>&1
			FWD_GENES=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa | sed 's/ # /	/g' |\
				awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $4}}' | wc -l )
			REV_GENES=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa | sed 's/ # /	/g' |\
				awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == -1) {print $4}}' | wc -l )
			if [ $FWD_GENES -ge $REV_GENES ] && [ $FWD_GENES -ge 1 ]; then
				START_BASE=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa | sed 's/ # /	/g' |\
					awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $2, ($3-$2)}}' |\
					sort -rg -k2,2 | head -n1 | cut -f1 )
				seqkit grep --quiet -p "${DTR_CONTIG}" ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta |\
				  seqkit restart --quiet -i ${START_BASE} > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
			elif [ $REV_GENES -ge 1 ]; then
				seqkit grep --quiet -p "${DTR_CONTIG}" ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta |\
				  seqkit seq --quiet -t DNA -r -p > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.fna
				prodigal -a ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.faa -i ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.fna -p meta -q >/dev/null 2>&1
				RC_FWD_GENES=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.faa | sed 's/ # /	/g' |\
					awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $4}}' | wc -l )
				if [ $RC_FWD_GENES -ge 1 ] ; then 
					START_BASE=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.faa | sed 's/ # /	/g' |\
						awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $2, ($3-$2)}}' |\
						sort -rg -k2,2 | head -n1 | cut -f1 )
					cat ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.fna |\
					  seqkit restart --quiet -i ${START_BASE} > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
				else
					echo "Can't find suitable ORF to set rotation of ${DTR_CONTIG} and will remain unrotated"
					seqkit grep --quiet -p "${DTR_CONTIG}" ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
				fi
			else
				echo "Can't find suitable ORF to set rotation of $nucl_fa and will remain unrotated"
				seqkit grep --quiet -p "${DTR_CONTIG}" ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
			fi
		done

		awk '{OFS=FS="\t"}{ if (NR != 1 && $4 == "NA") {print $1} }' ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv > ${TEMP_DIR}/hallmark_contigs_no_DTRS.txt

		if [ -s ${TEMP_DIR}/hallmark_contigs_no_DTRS.txt ] ; then
			seqkit grep --quiet -f ${TEMP_DIR}/hallmark_contigs_no_DTRS.txt\
			  ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta > ${TEMP_DIR}/hallmark_contigs_no_DTRS.fasta
		fi

		ALL_DTRS=$( find ${TEMP_DIR}/rotation -type f -name "*.rotate.fasta" )

		if [ -n "$ALL_DTRS" ] ; then
			echo "$ALL_DTRS" | while read SEQ ; do
				cat $SEQ
			done > ${TEMP_DIR}/hallmark_contigs_with_DTRS.rotated.fasta
		fi

		if [ -s ${TEMP_DIR}/hallmark_contigs_with_DTRS.rotated.fasta ] ; then
			cat ${TEMP_DIR}/hallmark_contigs_with_DTRS.rotated.fasta >> ${TEMP_DIR}/oriented_hallmark_contigs.fasta
		fi

		if [ -s ${TEMP_DIR}/hallmark_contigs_no_DTRS.fasta ] ; then
			cat ${TEMP_DIR}/hallmark_contigs_no_DTRS.fasta >> ${TEMP_DIR}/oriented_hallmark_contigs.fasta
		fi


	else
		echo "not wrapping"
		
		cp ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta ${TEMP_DIR}/oriented_hallmark_contigs.fasta
	fi
else
	echo "couldn't find threshold contigs with processed terminal repeats at ${TEMP_DIR}/threshold_trimmed_TRs_contigs.fasta"
fi


if [ ! -d ${TEMP_DIR}/reORF ]; then
	mkdir ${TEMP_DIR}/reORF
fi

if [ ! -d ${TEMP_DIR}/hallmark_tax ]; then
	mkdir ${TEMP_DIR}/hallmark_tax
fi

### evaluate ORF caller argument

if [ $CALLER == "prodigal" ] || [ $CALLER == "prodigal-gv" ] || [ $CALLER == "phanotate" ] ; then
	#- input: -#
	#-- ${TEMP_DIR}/contigs_to_keep.txt
	#- output: -#
	#-- ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt
	#-- ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt
	echo "forcing final ORF calls to be $CALLER"

	if [ $CALLER == "prodigal" ] || [ $CALLER == "prodigal-gv" ] ; then

		cp ${TEMP_DIR}/contigs_to_keep.txt ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt

	elif [ $CALLER == "phanotate" ] ; then

		cp ${TEMP_DIR}/contigs_to_keep.txt ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt

	fi

else

	### blastp-style mmseqs hallmark genes for taxonomy
	#- input: -#
	#-- ${TEMP_DIR}/hallmarks_for_keepcontigs1.txt
	#-- ${TEMP_DIR}/ORF_orig_contigs/split/*.faa (1 or more)
	#- output: -#
	#-- ${TEMP_DIR}/hallmark_tax/orig_hallmark_genes.faa
	#-- ${TEMP_DIR}/hallmark_tax/orig_hallmarks_align.tsv


	if [ -s ${TEMP_DIR}/hallmarks_for_keepcontigs1.txt ] && [ -n "$SPLIT_ORIG_AAs" ] ; then

		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BGreen}mmseqs of original hallmark genes for taxonomy calls ${MDYT}${Color_Off}"



		seqkit grep --quiet -f ${TEMP_DIR}/hallmarks_for_keepcontigs1.txt\
		  $SPLIT_ORIG_AAs > ${TEMP_DIR}/hallmark_tax/orig_hallmark_genes.faa

		if [ -s ${TEMP_DIR}/hallmark_tax/orig_hallmark_genes.faa ] ; then
			mmseqs createdb ${TEMP_DIR}/hallmark_tax/orig_hallmark_genes.faa ${TEMP_DIR}/hallmark_tax/orig_hallmark_genesDB -v 1

			mmseqs search ${TEMP_DIR}/hallmark_tax/orig_hallmark_genesDB\
			  ${C_DBS}/mmseqs_DBs/refseq_virus_prot_taxDB\
			  ${TEMP_DIR}/hallmark_tax/orig_hallmarks_resDB ${TEMP_DIR}/hallmark_tax/tmp -v 1 --start-sens 1 --sens-steps 3 -s 7

			mmseqs convertalis ${TEMP_DIR}/hallmark_tax/orig_hallmark_genesDB\
			  ${C_DBS}/mmseqs_DBs/refseq_virus_prot_taxDB\
			  ${TEMP_DIR}/hallmark_tax/orig_hallmarks_resDB ${TEMP_DIR}/hallmark_tax/orig_hallmarks_align.tsv\
			  --format-output query,target,pident,alnlen,evalue,theader,taxlineage -v 1

		else
			echo "couldn't find ${TEMP_DIR}/hallmark_tax/orig_hallmark_genes.faa"
			exit
		fi

	else
		echo "couldn't find ${TEMP_DIR}/contigs_to_keep.txt or $SPLIT_ORIG_AAs for mmseqs hallmarks"
		exit
	fi

	### parse taxonomy on hallmark gene mmseqs2 search and decide final ORF caller
	#- input: -#
	#-- ${TEMP_DIR}/hallmark_tax/orig_hallmarks_align.tsv
	#-- ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv
	#- output: -#
	#-- ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt
	#-- ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt
	#-- ${TEMP_DIR}/hallmark_tax/orf_caller_each_seq.tsv
	#----	fields
	#----	(contig	out_length_contig	query	target	pident	alnlen	evalue	theader	taxlineage	ORFcaller	pos	Note)

	if [ -e ${TEMP_DIR}/hallmark_tax/orig_hallmarks_align.tsv ] && [ -s ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BGreen}choosing ORF caller for each sequence ${MDYT}${Color_Off}"

		python ${CENOTE_SCRIPTS}/python_modules/orfcaller_decision1.py ${TEMP_DIR}/hallmark_tax/orig_hallmarks_align.tsv\
		  ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv ${TEMP_DIR}/hallmark_tax

	else
		echo "couldn't find ${TEMP_DIR}/hallmark_tax/orig_hallmarks_align.tsv or ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv"
		exit
	fi
fi

## redo ORF calls for everything. Some need phanotate, some were rotated
if [ -s ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt ] || [ -s ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BBlue}redoing ORF calls for each sequence ${MDYT}${Color_Off}"


	## adding contigs that had no hits in mmseqs search to list of contigs that need prodigal ORF calling
	if [ -s ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt ] ; then
		cat ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt >> ${TEMP_DIR}/hallmark_tax/taxed_seqs1.txt
	fi

	if [ -s ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt ] ; then
		cat ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt >> ${TEMP_DIR}/hallmark_tax/taxed_seqs1.txt
	fi
	if [ -s ${TEMP_DIR}/hallmark_tax/taxed_seqs1.txt ] ; then
		grep -v -f ${TEMP_DIR}/hallmark_tax/taxed_seqs1.txt\
			${TEMP_DIR}/contigs_to_keep.txt >> ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt
	fi


else
	echo "couldn't find prodigal_seqs1.txt or phanotate_seqs1.txt"
	exit
fi


### Prodigal for hallmark, rotated seqs
#- input: -#
#-- ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt
#-- ${TEMP_DIR}/oriented_hallmark_contigs.fasta
#- output: -#
#-- ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.faa
#-- ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.gff
#-- ${TEMP_DIR}/reORF/reORFcalled_all.faa (append)
#-- ${TEMP_DIR}/reORF/contig_gcodes1.txt

if [ -s ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt ] && [ -s ${TEMP_DIR}/oriented_hallmark_contigs.fasta ]; then


	seqkit grep --quiet -f ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt\
	  ${TEMP_DIR}/oriented_hallmark_contigs.fasta > ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.fna

	python ${CENOTE_SCRIPTS}/python_modules/pyrodigal_gv_runner.py ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.fna\
	  ${TEMP_DIR}/reORF $CPU $CALLER

	if [ -s ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.faa ] ; then
		cat ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.faa >> ${TEMP_DIR}/reORF/reORFcalled_all.faa

		## extract genetic code from prodigal files.
		cat ${TEMP_DIR}/hallmark_tax/prodigal_seqs1.txt | while read SEQ ; do 
			GCODE=$( grep -A1 "\"${SEQ}\"" ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.gff | tail -n1 |\
			  sed 's/.*transl_table=\([0-9]\{1,2\}\).*/\1/' )
			echo -e "${SEQ}\t${GCODE}"
		done > ${TEMP_DIR}/reORF/contig_gcodes1.txt

	else
		echo "can't find ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.faa"
	fi


else
	echo "no prodigal list. OK."
fi


### phanotate for hallmark, rotated seqs
#- input: -#
#-- ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt
#-- ${TEMP_DIR}/oriented_hallmark_contigs.fasta
#- output: -#
#-- ${TEMP_DIR}/reORF/phan_split/*bed (1 or more)
#-- ${TEMP_DIR}/reORF/phan_split/*.faa (1 or more)
#-- ${TEMP_DIR}/reORF/reORFcalled_all.faa (append)

if [ -s ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt ] ; then

	if [ ! -d ${TEMP_DIR}/reORF/phan_split ]; then
		mkdir ${TEMP_DIR}/reORF/phan_split
	fi

	seqkit grep --quiet -f ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt ${TEMP_DIR}/oriented_hallmark_contigs.fasta |\
	  seqkit split --quiet -j $CPU -p $CPU -O ${TEMP_DIR}/reORF/phan_split

	SPLIT_PHAN_CONTIGS=$( find ${TEMP_DIR}/reORF/phan_split -type f -name "*.fasta" )

	if [ -n "$SPLIT_PHAN_CONTIGS" ] ; then
		## this part actually just calls the coordinates

		PHAN_SEQS_L=$( cat ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt | wc -l )
		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BBlue}time update: running phanotate for re-ORF call on ${PHAN_SEQS_L} seqs ${MDYT}${Color_Off}"

		echo "$SPLIT_PHAN_CONTIGS" | sed 's/.fasta//g' |\
		  xargs -n 1 -I {} -P $CPU phanotate.py -f tabular {}.fasta -o {}.phan_genes.bad_fmt.tsv\
		    >${TEMP_DIR}/reORF/phan_split/phan_log.txt 2>&1


		SPLIT_PHAN_TABS=$( find ${TEMP_DIR}/reORF/phan_split -type f -name "*.phan_genes.bad_fmt.tsv" )

		for PHAN_TSV in $SPLIT_PHAN_TABS; do
			#echo $PHAN_TSV
			awk '{OFS=FS="\t"}{ if ($1 !~ /^#/) { if ($2>$1) {print $4, ($1-1), $2, $4"_"NR, $5, $3} else {print $4, ($2-1), $1, $4"_"NR, $5, $3}}}' ${PHAN_TSV} > ${PHAN_TSV%.phan_genes.bad_fmt.tsv}.phan_genes.bed
			rm ${PHAN_TSV}
		done


	else
		echo "can't find split phanotate contigs"

	fi

	PHAN_TABS=$( find ${TEMP_DIR}/reORF/phan_split -type f -name "*.phan_genes.bed" )
	if [ -n "$PHAN_TABS" ] ; then

		## this part extracts the gene seqs
		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BBlue}time update: running bedtools to extract phanotate ORF calls ${MDYT}${Color_Off}"


		echo "$PHAN_TABS" | sed 's/.phan_genes.bed//g' |\
		  xargs -n 1 -I {} -P $CPU bedtools getfasta -fi {}.fasta -bed {}.phan_genes.bed -fo {}.phan_genes.fasta -s -nameOnly

		## this part translates the gene seqs
		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BBlue}time update: running seqkit translate to on phanotate ORF calls ${MDYT}${Color_Off}"

		echo "$PHAN_TABS" | sed 's/.phan_genes.bed//g' |\
		  xargs -n 1 -I {} -P $CPU seqkit translate --quiet -x -T 11 {}.phan_genes.fasta -o {}.faa >/dev/null 2>&1 


		echo "$PHAN_TABS" | sed 's/.phan_genes.bed//g' | while read PHAN ; do
			sed 's/(.*//g' ${PHAN}.faa
		done >> ${TEMP_DIR}/reORF/reORFcalled_all.faa

	fi

else
	echo "no phanotate list. OK."
fi


### split ORFs of hallmark contigs
#- input: -#
#-- ${TEMP_DIR}/reORF/reORFcalled_all.faa
#- output: -#
#-- ${TEMP_DIR}/reORF_pyhmmer1_split/*.faa (1 or more)

if [ -s ${TEMP_DIR}/reORF/reORFcalled_all.faa ] ; then

	if [ ! -d ${TEMP_DIR}/reORF_pyhmmer1_split ]; then
		mkdir ${TEMP_DIR}/reORF_pyhmmer1_split
	fi

	if [ ! -d ${TEMP_DIR}/reORF_pyhmmer2_split ]; then
		mkdir ${TEMP_DIR}/reORF_pyhmmer2_split
	fi
	
	if [ ! -d ${TEMP_DIR}/reORF_mmseqs_combined ]; then
		mkdir ${TEMP_DIR}/reORF_mmseqs_combined
	fi
	

	seqkit split --quiet -j $CPU -p $CPU -O ${TEMP_DIR}/reORF_pyhmmer1_split ${TEMP_DIR}/reORF/reORFcalled_all.faa

else
	echo "couldn't find  ${TEMP_DIR}/reORF/reORFcalled_all.faa for splitting"
	exit
fi


### pyhmmer1 (virion and rep hallmarks)
#- input: -#
#-- ${TEMP_DIR}/reORF_pyhmmer1_split/*.faa (1 or more)
#- output: -#
#-- ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt
#-- ${TEMP_DIR}/reORF_pyhmmer2_split/*.no1.faa (1 or more)

SPLIT_REORF_AAs=$( find ${TEMP_DIR}/reORF_pyhmmer1_split -type f -name "*.faa" )

if [ -n "$SPLIT_REORF_AAs" ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BCyan}time update: running pyhmmer hallmark dbs on reORFs ${MDYT}${Color_Off}"

	if [ ! -d ${TEMP_DIR}/virion_reORF_pyhmmer ]; then
		mkdir ${TEMP_DIR}/virion_reORF_pyhmmer
	fi
	
	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/reORF_pyhmmer1_split ${TEMP_DIR}/virion_reORF_pyhmmer\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/Virion_HMMs.h3m $CPU 1e-5 0.1

	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/reORF_pyhmmer1_split ${TEMP_DIR}/rep_reORF_pyhmmer\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/DNA_rep_HMMs.h3m $CPU 1e-5 0.1

	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/reORF_pyhmmer1_split ${TEMP_DIR}/rdrp_reORF_pyhmmer\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/RDRP_HMMs.h3m $CPU 1e-5 0.8

	if [ -s ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv ]\
	  && [ -s ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv ]\
	  && [ -s ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv ]; then
		awk FNR!=1 ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
		  ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
		  ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv |\
		  cut -f1 > ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt

	elif [ -s ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv ]\
	  && [ -s ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
		awk FNR!=1 ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
		  ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv |\
		  cut -f1 > ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt

	elif [ -s ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv ]\
	  && [ -s ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
		awk FNR!=1 ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
		  ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv |\
		  cut -f1 > ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt

	elif [ -s ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
		awk FNR!=1 ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv |\
		  cut -f1 > ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt

	elif [ -s ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
		awk FNR!=1 ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv |\
		  cut -f1 > ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt

	elif [ -s ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
		awk FNR!=1 ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv |\
		  cut -f1 > ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt
	fi

	if [ -s ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt ] ; then
		echo "$SPLIT_REORF_AAs" | while read AA ; do
			BASE_AA=$( basename $AA )
			seqkit grep --quiet -j $CPU -v -f ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt\
			  $AA > ${TEMP_DIR}/reORF_pyhmmer2_split/${BASE_AA%.faa}.no1.faa
		done

	else
		echo "$SPLIT_REORF_AAs" | while read AA ; do
			BASE_AA=$( basename $AA )
			cp $AA ${TEMP_DIR}/reORF_pyhmmer2_split/${BASE_AA%.faa}.no1.faa
		done
	fi


else
	echo "couldn't find prodigal AA seqs in ${TEMP_DIR}/reORF_pyhmmer1_split"
fi

## pyhmmer2 (other virus HMMs)
#- input: -#
#-- ${TEMP_DIR}/reORF_pyhmmer2_split/*.no1.faa (1 or more)
#- output: -#
#-- ${TEMP_DIR}/comm_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/comm_reORF_pyhmmer/hit_this_round1.txt
#-- ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa 

SECOND_REORF_AAs=$( find ${TEMP_DIR}/reORF_pyhmmer2_split -type f ! -size 0 -name "*.no1.faa" )

if [ -n "$SECOND_REORF_AAs" ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BCyan}time update: running pyhmmer additional annotation HMMs on reORFs ${MDYT}${Color_Off}"


	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/reORF_pyhmmer2_split ${TEMP_DIR}/comm_reORF_pyhmmer\
	  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/Useful_Annotation_HMMs.h3m $CPU 1e-6 0.1


	if [ -s ${TEMP_DIR}/comm_reORF_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
		tail -n+2 ${TEMP_DIR}/comm_reORF_pyhmmer/pyhmmer_report_AAs.tsv | cut -f1 > ${TEMP_DIR}/comm_reORF_pyhmmer/hit_this_round1.txt

		echo "$SECOND_REORF_AAs" | while read AA ; do
			seqkit grep --quiet -j $CPU -v -f ${TEMP_DIR}/comm_reORF_pyhmmer/hit_this_round1.txt $AA >> ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa
		done

	else
		echo "$SECOND_REORF_AAs" | while read AA ; do
			cat $AA
		done > ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa
	fi


else
	echo "couldn't find prodigal AA seqs in ${TEMP_DIR}/reORF_pyhmmer2_split"
fi




### mmseqs2 with CDD profiles
#- input: -#
#-- ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa
#- output: -#
#-- ${TEMP_DIR}/reORF_mmseqs_combined/no2_seqs_CDD.tsv
#-- ${TEMP_DIR}/reORF_mmseqs_combined/summary_no2_AAs_vs_CDD.besthit.tsv
#----	fields
#----	(query	target	sequence_identity	align_length	evalue	bitscore	cdd_num	cdd_accession	shortname	description	other_num)

if [ -s ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BCyan}time update: running mmseqs2 against CDD with ORFs not annotated with HMMs ${MDYT}${Color_Off}"

	mmseqs createdb ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2DB -v 1

	mmseqs search ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2DB\
	  ${C_DBS}/mmseqs_DBs/CDD ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2_vs_CDD_resDB\
	  ${TEMP_DIR}/reORF_mmseqs_combined/tmp -s 4 -v 1

	mmseqs convertalis ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2DB\
	  ${C_DBS}/mmseqs_DBs/CDD\
	  ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2_vs_CDD_resDB ${TEMP_DIR}/reORF_mmseqs_combined/no2_seqs_CDD.tsv\
	  --format-output query,target,pident,alnlen,evalue,bits -v 1

	python ${CENOTE_SCRIPTS}/python_modules/parse_mmseqs_cdd_results1.py ${TEMP_DIR}/reORF_mmseqs_combined/no2_seqs_CDD.tsv\
	  ${C_DBS}/mmseqs_DBs/cddid_all.tbl ${TEMP_DIR}/reORF_mmseqs_combined

else
	echo "couldn't find ${TEMP_DIR}/reORF_mmseqs_combined/all_AA_seqs.no2.faa for mmseqs CDD"
fi


### Assess all the gene annotations to make annotation table and determine virus chunks
#- input: -#
#-- ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv
#-- ${TEMP_DIR}/reORF/phan_split/*bed (1 or more)
#-- ${TEMP_DIR}/reORF/pyrodigal_gv_AAs.prod.gff
#-- ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/comm_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv
#-- ${TEMP_DIR}/reORF_mmseqs_combined/summary_no2_AAs_vs_CDD.besthit.tsv
#-- ${HALL_TYPE} (argument)
#- output: -#
#-- ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.tsv
#----	fields
#----	(contig	gene_name	gene_orient	gene_start	gene_stop	contig_length	dtr_seq	evidence_acession	evidence_description	Evidence_source	vscore_category)
#-- ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.hallmarks.bed
#-- ${TEMP_DIR}/assess_prune/prune_figures/*chunks.tsv (0 or more)
#----	fields
#----	(contig	left_cutoff	right_cutoff	chunk_number)
#-- ${TEMP_DIR}/assess_prune/prune_figures/*figures.pdf (0 or more)

if [ -s ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv ] && [ -s ${C_DBS}/viral_cdds_and_pfams_191028.txt ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BCyan}time update: assessing each gene on all contigs and scoring contigs for virusness  ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/assess_virus_genes1.py ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv\
	  ${TEMP_DIR}/reORF/phan_split ${TEMP_DIR}/reORF ${TEMP_DIR}/virion_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
	  ${TEMP_DIR}/comm_reORF_pyhmmer/pyhmmer_report_AAs.tsv ${TEMP_DIR}/rep_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
	  ${TEMP_DIR}/rdrp_reORF_pyhmmer/pyhmmer_report_AAs.tsv\
	  ${TEMP_DIR}/reORF_mmseqs_combined/summary_no2_AAs_vs_CDD.besthit.tsv ${C_DBS}/viral_cdds_and_pfams_191028.txt \
	  ${TEMP_DIR}/assess_prune "${HALL_TYPE}" ${PROPHAGE}

else
	echo "couldn't start assess and prune script"
fi

### find coordinates of viruses within input contigs
#- input: -#
#--  ${TEMP_DIR}/assess_prune/prune_figures/*chunks.tsv
#--  ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.hallmarks.bed
#- output: -#
#--  ${TEMP_DIR}/assess_prune/indiv_seqs/*.bed
#--  ${TEMP_DIR}/assess_prune/indiv_seqs/*viruses.tsv

CHUNK_FILES=$( find ${TEMP_DIR}/assess_prune/prune_figures -type f -name "*chunks.tsv" )

if [ -s ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.hallmarks.bed ] && [ -n "$CHUNK_FILES" ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BPurple}time update: pulling out virus parts of contigs >= 10 kb ${MDYT}${Color_Off}"

	if [ ! -d ${TEMP_DIR}/assess_prune/indiv_seqs ]; then
		mkdir ${TEMP_DIR}/assess_prune/indiv_seqs
	fi

	echo ""
	echo "$CHUNK_FILES" | while read SEQ_CHUNK ; do

		B_SEQ=$( basename $SEQ_CHUNK )

		tail -n+2 $SEQ_CHUNK > ${TEMP_DIR}/assess_prune/indiv_seqs/${B_SEQ%.tsv}.bed

		bedtools intersect -c -a ${TEMP_DIR}/assess_prune/indiv_seqs/${B_SEQ%.tsv}.bed\
		  -b ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.hallmarks.bed |\
		  awk -v minh="$LIN_MINIMUM_DOMAINS" '{OFS=FS="\t"}{if ($5 >= minh) {print}}' > ${TEMP_DIR}/assess_prune/indiv_seqs/${B_SEQ%.chunks.tsv}.viruses.tsv

	done
else
	echo "couldn't find chunk files"
fi


### adjust viruses based on pruning assessment
#- input: -#
#--  ${TEMP_DIR}/assess_prune/indiv_seqs/*viruses.tsv
#--  ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.tsv
#--  ${TEMP_DIR}/oriented_hallmark_contigs.fasta
#- output: -#
#--  ${TEMP_DIR}/contig_gene_annotation_summary.pruned.tsv
#--  ${TEMP_DIR}/prune_coords.bed
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta
#--  ${TEMP_DIR}/hypothetical_proteins.after_chunk.txt

if [ -s ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.tsv ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BPurple}time update: reconfiguring gene/contig coordinates after prune ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/adjust_viruses1.py ${TEMP_DIR}/assess_prune/indiv_seqs\
	  ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.tsv ${TEMP_DIR} ${PROPHAGE}


	seqkit subseq --quiet -j $CPU --bed ${TEMP_DIR}/prune_coords.bed ${TEMP_DIR}/oriented_hallmark_contigs.fasta |\
	  sed 's/>.* />/g' > ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta


else
	echo "couldn't find ${TEMP_DIR}/assess_prune/contig_gene_annotation_summary.tsv for adjust viruses"

fi

### pyhmmer3 (selected PHROGS HMMS)
#- input: -#
#--  ${TEMP_DIR}/hypothetical_proteins.after_chunk.txt
#--  ${TEMP_DIR}/reORF/reORFcalled_all.faa
#- output: -#
#--  ${TEMP_DIR}/phrogs_pyhmmer/pyhmmer_report_AAs.tsv
#--  ${TEMP_DIR}/phrogs_pyhmmer/hit_this_round1.txt

if [ -s ${TEMP_DIR}/hypothetical_proteins.after_chunk.txt ]; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BRed}time update: running pyhmmer on PHROGs HMMs on reORFs ${MDYT}${Color_Off}"

	if [ ! -d ${TEMP_DIR}/reORF_phrogs_split ] ; then
		mkdir ${TEMP_DIR}/reORF_phrogs_split
	fi


	seqkit grep --quiet -f ${TEMP_DIR}/hypothetical_proteins.after_chunk.txt ${TEMP_DIR}/reORF/reORFcalled_all.faa |\
	  seqkit split --quiet -j $CPU -p $CPU -O ${TEMP_DIR}/reORF_phrogs_split

	PHROGS_FASTAs=$( find ${TEMP_DIR}/reORF_phrogs_split -type f ! -size 0 -name "*.fasta" )

	if [ -n "$PHROGS_FASTAs" ] ; then
		for PFAST in $PHROGS_FASTAs ; do
			mv $PFAST ${PFAST%.fasta}.faa
		done
	fi

	PHROGS_AAs=$( find ${TEMP_DIR}/reORF_phrogs_split -type f ! -size 0 -name "*.faa" )

	if [ -n "$PHROGS_AAs" ] ; then


		python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/reORF_phrogs_split ${TEMP_DIR}/phrogs_pyhmmer\
		  ${C_DBS}/hmmscan_DBs/${HMM_DBS}/phrogs_for_ct.h3m $CPU 1e-4 0.1

		if [ -s ${TEMP_DIR}/phrogs_pyhmmer/pyhmmer_report_AAs.tsv ] ; then
			tail -n+2 ${TEMP_DIR}/phrogs_pyhmmer/pyhmmer_report_AAs.tsv | cut -f1 > ${TEMP_DIR}/phrogs_pyhmmer/hit_this_round1.txt

			echo "$PHROGS_AAs" | while read AA ; do
				seqkit grep --quiet -j $CPU -v -f ${TEMP_DIR}/phrogs_pyhmmer/hit_this_round1.txt $AA >> ${TEMP_DIR}/phrogs_pyhmmer/all_AA_seqs.no_phrogs.faa
			done

		else
			if [ ! -d ${TEMP_DIR}/phrogs_pyhmmer ] ; then
				mkdir ${TEMP_DIR}/phrogs_pyhmmer
			fi

			echo "$PHROGS_AAs" | while read AA ; do
				cat $AA
			done > ${TEMP_DIR}/phrogs_pyhmmer/all_AA_seqs.no_phrogs.faa

		fi

	else
		echo "couldn't find seqs for phrogs HMMs"
	fi
else
	echo "couldn't find hypothetical proteins for phrogs HMMscan"
fi


### hhblits/hhsearch (installed DBs)
#- input: -#
#--  ${TEMP_DIR}/phrogs_pyhmmer/all_AA_seqs.no_phrogs.faa
#- output: -#
#--  ${TEMP_DIR}/hhpred/hhpred_report_AAs.tsv

if  [[ $HHSUITE_TOOL = "hhsearch" ]] || [[ $HHSUITE_TOOL = "hhblits" ]] ; then
	if [ -s ${TEMP_DIR}/phrogs_pyhmmer/all_AA_seqs.no_phrogs.faa ] ; then

		MDYT=$( date +"%m-%d-%y---%T" )
		echo -e "${BRed}time update: hhsuite search of hypothetical proteins ${MDYT}${Color_Off}"

		if [ ! -d ${TEMP_DIR}/hhpred ]; then
			mkdir ${TEMP_DIR}/hhpred
		fi

		if [ ! -d ${TEMP_DIR}/hhpred/AA_files ]; then
			mkdir ${TEMP_DIR}/hhpred/AA_files
		fi

		seqkit split --quiet -j $CPU -s 1 -O ${TEMP_DIR}/hhpred/AA_files ${TEMP_DIR}/phrogs_pyhmmer/all_AA_seqs.no_phrogs.faa

		HH_AAs=$( find ${TEMP_DIR}/hhpred/AA_files -type f -name "*.faa" )

		if [ -n "$HH_AAs" ] && [[ $HHSUITE_TOOL = "hhblits" ]] ; then
			echo "$HH_AAs" | sed 's/.faa//g' |\
			  xargs -n 1 -I {} -P $CPU hhblits -i {}.faa ${HHSUITE_DB_STR} -o {}.out.hhr\
			  -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >${TEMP_DIR}/hhpred/hhblits_logs.txt 2>&1

		elif [ -n "$HH_AAs" ] && [[ $HHSUITE_TOOL = "hhsearch" ]] ; then
			echo "$HH_AAs" | sed 's/.faa//g' |\
			  xargs -n 1 -I {} -P $CPU hhsearch -i {}.faa ${HHSUITE_DB_STR} -o {}.out.hhr\
			  -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >${TEMP_DIR}/hhpred/hhblits_logs.txt 2>&1

		else
			echo "AA seqs for hhpred not found"
		fi

		###parse annoying hhpred output files
		python ${CENOTE_SCRIPTS}/python_modules/hhpred_to_table.py ${TEMP_DIR}/hhpred/AA_files ${TEMP_DIR}/hhpred


	else
		echo "no list of proteins for hhpred"
	fi
else
	echo "not running hhsearch/hhblits"
fi

### tRNAscan-SE
#- input: -#
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta
#- output: -#
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.tRNAscan.tsv
#----	fields
#----	(seq_name	trna_number	begin	end	type	anticodon	intron_begin	intron_end	score	note)

if [ -s ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BRed}time update: tRNAscan-SE on virus contigs ${MDYT}${Color_Off}"

	tRNAscan-SE -Q -G -o ${TEMP_DIR}/oriented_hallmark_contigs.pruned.tRNAscan.tsv --brief\
	  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta >${TEMP_DIR}/trnascan_se_logs.txt 2>&1

	## python script to remove overlapping genes and replace them with tRNAs
else
	echo "couldn't find oriented pruned contigs for tRNAscan-SE"
fi


### Mapping reads to virus contigs
#- input: -#
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta
#--  Fastq files in $READS variable
#- output: -#
#--  ${TEMP_DIR}/mapping_reads/oriented_hallmark_contigs.pruned.coverage.tsv

if [ -s ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta ] && [ "${READS}" != "none" ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BYellow}time update: mapping reads to virus contigs ${MDYT}${Color_Off}"

	if [ ! -d ${TEMP_DIR}/mapping_reads ]; then
		mkdir ${TEMP_DIR}/mapping_reads
	fi

	minimap2 -t $CPU -ax sr ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta ${READS} |\
	  samtools sort - | samtools coverage -o ${TEMP_DIR}/mapping_reads/oriented_hallmark_contigs.pruned.coverage.tsv -

else
	echo "not mapping reads"
fi


### redo hallmark taxonomy on reORF viruses/chunks
#- input: -#
#--  ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt
#--  ${TEMP_DIR}/reORF/reORFcalled_all.faa
#--  ${TEMP_DIR}/contig_gene_annotation_summary.pruned.tsv
#- output: -#
#--  ${TEMP_DIR}/final_taxonomy/hallmark_proteins.faa
#--  ${TEMP_DIR}/final_taxonomy/hallmark_proteins_align.tsv
#--  ${TEMP_DIR}/final_taxonomy/virus_taxonomy_summary.tsv
#----	fields
#----	(contig	chunk_name	taxon	taxonomy_hierarchy	taxon_level	avg_hallmark_AAI_to_ref)

if [ ! -d ${TEMP_DIR}/final_taxonomy ]; then
	mkdir ${TEMP_DIR}/final_taxonomy
fi

if [ -s ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BYellow}time update: reassessing taxonomy on final virus seqs with mmseqs2 ${MDYT}${Color_Off}"


	seqkit grep --quiet -f ${TEMP_DIR}/virion_reORF_pyhmmer/hit_this_round1.txt\
	  ${TEMP_DIR}/reORF/reORFcalled_all.faa > ${TEMP_DIR}/final_taxonomy/hallmark_proteins.faa

	if [ -s ${TEMP_DIR}/final_taxonomy/hallmark_proteins.faa ] ; then
		mmseqs createdb ${TEMP_DIR}/final_taxonomy/hallmark_proteins.faa ${TEMP_DIR}/final_taxonomy/hallmark_proteinsDB -v 1

		mmseqs search ${TEMP_DIR}/final_taxonomy/hallmark_proteinsDB\
		  ${C_DBS}/mmseqs_DBs/refseq_virus_prot_taxDB\
		  ${TEMP_DIR}/final_taxonomy/hallmark_proteins_resDB ${TEMP_DIR}/final_taxonomy/tmp -v 1 --start-sens 1 --sens-steps 3 -s 7

		mmseqs convertalis ${TEMP_DIR}/final_taxonomy/hallmark_proteinsDB\
		  ${C_DBS}/mmseqs_DBs/refseq_virus_prot_taxDB\
		  ${TEMP_DIR}/final_taxonomy/hallmark_proteins_resDB ${TEMP_DIR}/final_taxonomy/hallmark_proteins_align.tsv\
		  --format-output query,target,pident,alnlen,evalue,theader,taxlineage -v 1

	else
		echo "couldn't find ${TEMP_DIR}/final_taxonomy/hallmark_proteins.faa"
	fi
else
	touch ${TEMP_DIR}/final_taxonomy/hallmark_proteins_align.tsv
fi

if [ -e ${TEMP_DIR}/final_taxonomy/hallmark_proteins_align.tsv ] ; then

	##python script to merge these results with gene annotation table, find best hit, decide taxon
	python ${CENOTE_SCRIPTS}/python_modules/vote_taxonomy.py ${TEMP_DIR}/final_taxonomy/hallmark_proteins_align.tsv\
	  ${TEMP_DIR}/contig_gene_annotation_summary.pruned.tsv ${TEMP_DIR}/final_taxonomy

else
	echo "couldn't find mmseqs hallmark tax table for final taxonomy assessment"
fi


#-#- blastn-style mmseqs2 taxonomy for species level

### Make sequin-formatted fsa file
#- input: -#
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta
#--  ${TEMP_DIR}/final_taxonomy/virus_taxonomy_summary.tsv
#--  ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv
#--  ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt
#--  ${TEMP_DIR}/reORF/contig_gcodes1.txt
#- output: -#
#--  ${run_title}/sequin_and_genome_maps/*.fsa (1 or more)

if [ -s ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta ] &&\
   [ -s ${TEMP_DIR}/final_taxonomy/virus_taxonomy_summary.tsv ] &&\
   [ -s ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BYellow}time update: Making genome map and sequin files ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/make_sequin_fsas.py ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta\
	  ${TEMP_DIR}/final_taxonomy/virus_taxonomy_summary.tsv\
	  ${TEMP_DIR}/threshold_contigs_terminal_repeat_summary.tsv ${TEMP_DIR} ${run_title}/sequin_and_genome_maps\
	  ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt ${TEMP_DIR}/reORF/contig_gcodes1.txt $ISO_SOURCE\
	  $COLLECT_DATE $META_TYPE $SRR $SRX $BIOSAMP $PRJ $MOL_TYPE $DATA_SOURCE

else

	echo "couldn't find files to make fsa's"
fi

### Make sequin-formatted tbl annotation file
#- input: -#
#--  ${TEMP_DIR}/contig_gene_annotation_summary.pruned.tsv
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.tRNAscan.tsv
#--  ${TEMP_DIR}/phrogs_pyhmmer/pyhmmer_report_AAs.tsv
#--  ${TEMP_DIR}/hhpred/hhpred_report_AAs.tsv
#- output: -#
#--  ${run_title}/sequin_and_genome_maps/*.tbl
#--  ${run_title}/final_ORF_list.txt
#--  ${run_title}/final_genes_to_contigs_annotation_summary.tsv


if [ -s ${TEMP_DIR}/contig_gene_annotation_summary.pruned.tsv ] ; then
	#-#-  remove overlapping tRNAs/genes and replace them with tRNAs

	## sequin tbl
	python ${CENOTE_SCRIPTS}/python_modules/make_sequin_tbls.py ${TEMP_DIR}/contig_gene_annotation_summary.pruned.tsv\
	  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.tRNAscan.tsv ${TEMP_DIR}/phrogs_pyhmmer/pyhmmer_report_AAs.tsv\
	  ${TEMP_DIR}/hhpred/hhpred_report_AAs.tsv ${run_title}/sequin_and_genome_maps

else
	echo "couldn't find annotation file for tbl generation"

fi

### Make sequin-formatted comment (cmt) file
#- input: -#
#--  ${run_title}/sequin_and_genome_maps/*.fsa
#--  ${TEMP_DIR}/mapping_reads/oriented_hallmark_contigs.pruned.coverage.tsv
#- output: -#
#--  ${run_title}/sequin_and_genome_maps/*.cmt

FSA_FILES=$( find ${run_title}/sequin_and_genome_maps -type f -name "*fsa" )

if [ -n "$FSA_FILES" ] ; then
	for REC in $FSA_FILES ; do
		if [ -s ${TEMP_DIR}/mapping_reads/oriented_hallmark_contigs.pruned.coverage.tsv ] ; then
			COVERAGE=$( awk -v SEQNAME="${REC%.fsa}" '{OFS=FS="\t"}{ if ($1 == SEQNAME) {print $7}}' \
				${TEMP_DIR}/mapping_reads/oriented_hallmark_contigs.pruned.coverage.tsv | head -n1 )

		else
			COVERAGE=1
		fi

		echo "StructuredCommentPrefix	##Genome-Assembly-Data-START##" > ${REC%.fsa}.cmt
		echo "Assembly Method	${ASSEMBLER}" >> ${REC%.fsa}.cmt
		echo "Genome Coverage	"$COVERAGE"x" >> ${REC%.fsa}.cmt
		echo "Sequencing Technology	Illumina" >> ${REC%.fsa}.cmt
		echo "Annotation Pipeline	Cenote-Taker 3" >> ${REC%.fsa}.cmt
		echo "URL	https://github.com/mtisza1/Cenote-Taker3" >> ${REC%.fsa}.cmt	
	done
fi

### Merge final virus seqs
#- input: -#
#--  ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta
#- output: -#
#--  ${run_title}/${run_title}_virus_sequences.fna

if [ -s ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta ] ; then
	seqkit replace --quiet -p "\s.+" ${TEMP_DIR}/oriented_hallmark_contigs.pruned.fasta > ${run_title}/${run_title}_virus_sequences.fna

else
	echo "couldn't find file for final virus seqs"
fi

### Merge final translated ORFs
#- input: -#
#--  ${run_title}/final_ORF_list.txt
#--  ${TEMP_DIR}/reORF/reORFcalled_all.faa
#- output: -#
#--  ${run_title}/${run_title}_virus_AA.faa

if [ -s ${run_title}/final_ORF_list.txt ] ; then
	seqkit grep --quiet -f ${run_title}/final_ORF_list.txt\
	  ${TEMP_DIR}/reORF/reORFcalled_all.faa | seqkit replace --quiet -p "\s.+" > ${run_title}/${run_title}_virus_AA.faa

	rm ${run_title}/final_ORF_list.txt
else
	echo "couldn't find list of final ORFs"
fi

### Run tbl2asn
#- input: -#
#--  ${run_title}/sequin_and_genome_maps/*.fsa
#--  ${run_title}/sequin_and_genome_maps/*.tbl
#--  ${run_title}/sequin_and_genome_maps/*.cmt
#- output: -#
#--  ${run_title}/sequin_and_genome_maps/*.gbf
#--  ${run_title}/sequin_and_genome_maps/*.sqn
#--  ${run_title}/sequin_and_genome_maps/*.val

if [ -s ${TEMPLATE_FILE} ] ; then
	tbl2asn -V vb -t ${TEMPLATE_FILE} -X C -p ${run_title}/sequin_and_genome_maps >\
	  ${run_title}/sequin_and_genome_maps/tbl2asn.log 2>&1
else
	echo "could not find template file for tbl2asn"
fi

### Make run summary files
#- input: -#
#--  ${run_title}/final_genes_to_contigs_annotation_summary.tsv
#--  ${TEMP_DIR}/contig_name_map.tsv
#--  ${TEMP_DIR}/final_taxonomy/virus_taxonomy_summary.tsv
#--  ${run_title}/sequin_and_genome_maps/*.fsa
#--  ${TEMP_DIR}/reORF/contig_gcodes1.txt
#--  ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt
#- output: -#
#--  ${run_title}/${run_title}_virus_summary.tsv
#----	fields
#----	(contig	input_name	organism	virus_seq_length	end_feature	gene_count	virion_hallmark_count
#----	 rep_hallmark_count	virion_hallmark_genes	rep_hallmark_genes	taxonomy_hierarchy	ORF_caller)
#--  ${run_title}/${run_title}_prune_summary.tsv
#----	fields
#----	(contig	contig_length	chunk_length	chunk_name	chunk_start	chunk_stop)

if [ -s ${run_title}/final_genes_to_contigs_annotation_summary.tsv ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo -e "${BYellow}time update: Making virus summary table ${MDYT}${Color_Off}"

	python ${CENOTE_SCRIPTS}/python_modules/virus_summary.py ${TEMP_DIR}/contig_name_map.tsv \
	  ${run_title}/final_genes_to_contigs_annotation_summary.tsv ${TEMP_DIR}/final_taxonomy/virus_taxonomy_summary.tsv \
	  ${run_title}/sequin_and_genome_maps ${run_title} ${TEMP_DIR}/reORF/contig_gcodes1.txt\
	  ${TEMP_DIR}/hallmark_tax/phanotate_seqs1.txt $CALLER

	python ${CENOTE_SCRIPTS}/python_modules/summary_statement.py\
	  ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta $LENGTH_MINIMUM\
	  ${run_title}/${run_title}_virus_summary.tsv ${run_title}/${run_title}_prune_summary.tsv\
	  ${run_title}/final_genes_to_contigs_annotation_summary.tsv $ANNOTATION_MODE

else
	echo "couldn't find files to make run summary"

fi


# gtf/gff

MDYT=$( date +"%m-%d-%y---%T" )
echo -e "${BYellow}Cenote-Taker finishing now ${MDYT}${Color_Off}"

echo -e "${BPurple}output: ${run_title}${Color_Off}"





