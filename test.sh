#!/usr/bin/env bash
# Author: Dominik Harmim <harmim6@gmail.com>

OUT='vid'
SOURCE="$OUT.cpp"
BUILD_OPTIONS='-std=c++14'
#BUILD_OPTIONS='-std=c++14 -D DEBUG'
RUN_OPTIONS='-q'

if [[ "$#" -ne 1 ]]; then
	echo "Expecting one argument: $0 altitudes" >&2
	exit 1
fi

if ! [[ "$1" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
	echo 'Invalid format of altitudes. Expecting "^[0-9]+(,[0-9]+)*$".' >&2
	exit 1
fi

ALTS_COUNT=`echo "$1" | grep -o ',' | wc -l | xargs`
if [[ "$ALTS_COUNT" -le 1 ]]; then NP=1
else
	ALTS_COUNT_ROUNDED=`echo "l($ALTS_COUNT) / l(2)" | bc -l`
	ALTS_COUNT_ROUNDED=`python -c "from math import ceil; \
		print int(ceil($ALTS_COUNT_ROUNDED))"`
	NP=`echo "scale=0; 2^$ALTS_COUNT_ROUNDED / 2" | bc -l`
fi

if [[ `uname -s` = "Darwin" ]]; then
	export OMPI_MCA_btl='self,tcp'
	export PMIX_MCA_gds='^ds12'
	RUN_OPTIONS="$RUN_OPTIONS --oversubscribe"
else
	BUILD_OPTIONS="$BUILD_OPTIONS --prefix /usr/local/share/OpenMPI"
	RUN_OPTIONS="$RUN_OPTIONS --prefix /usr/local/share/OpenMPI"
fi

mpic++ ${BUILD_OPTIONS} -o "$OUT" "$SOURCE"
CODE="$?"
if [[ "$CODE" -ne 0 ]]; then exit "$CODE"; fi

mpirun ${RUN_OPTIONS} -np "$NP" "$OUT" "$1"
CODE="$?"
rm -f "$OUT"
if [[ "$CODE" -ne 0 ]]; then exit "$CODE"; fi
