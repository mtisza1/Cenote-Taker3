#!/bin/bash

SEQ_LIST=$1
PROD_GFFS=$2
BIG_GFF=$3
GCODES=$4

cat $PROD_GFFS > $BIG_GFF

if [ -s $SEQ_LIST ] ; then
	cat $SEQ_LIST | while read SEQ ; do 
		GCODE=$( grep -A1 "\"${SEQ}\"" $BIG_GFF | tail -n1 |\
		  sed 's/.*transl_table=\([0-9]\{1,2\}\).*/\1/' )
		echo -e "${SEQ}\t${GCODE}"
	done > $GCODES
else
	touch $GCODES
fi