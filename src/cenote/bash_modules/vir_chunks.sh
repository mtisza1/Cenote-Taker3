#!/bin/bash

HM_BED=$1
PRUNE_DIR=$2
LIN_HM=$3
INDIV_DIR=$4

CHUNK_FILES=$( find ${PRUNE_DIR}/ -type f -name "*chunks.tsv" )

if [ -s ${HM_BED} ] && [ -n "$CHUNK_FILES" ] ; then

	#MDYT=$( date +"%m-%d-%y---%T" )
	#echo -e "${BPurple}time update: pulling out virus parts of contigs >= 10 kb ${MDYT}${Color_Off}"
	echo ""
	if [ ! -d ${INDIV_DIR} ]; then
		mkdir ${INDIV_DIR}
	fi

	echo "$CHUNK_FILES" | while read SEQ_CHUNK ; do
		B_SEQ=$( basename $SEQ_CHUNK )

		tail -n+2 $SEQ_CHUNK > ${INDIV_DIR}/${B_SEQ%.tsv}.bed

		bedtools intersect -c -a ${INDIV_DIR}/${B_SEQ%.tsv}.bed\
		  -b ${HM_BED} |\
		  awk -v minh="$LIN_HM"\
		    '{OFS=FS="\t"}{if ($5 >= minh) {print}}' > ${INDIV_DIR}/${B_SEQ%.chunks.tsv}.viruses.tsv

	done
else
	echo "couldn't find chunk files"
fi

touch ${INDIV_DIR}/.finished_vir_chunks