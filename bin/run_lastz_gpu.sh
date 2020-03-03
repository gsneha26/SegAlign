#!/bin/bash

if [ $# -eq 1 ] || [ $# -eq 0 ] 
then
	wga $1
else
	ref=$(basename -s .fa $1)
	query=$(basename -s .fa $2)
	
	FOLDER1=$PWD/output_$RANDOM/
	FOLDER=$FOLDER1/data_${ref}_${query}_$RANDOM/
	mkdir $FOLDER1
	mkdir $FOLDER
	
	cd $FOLDER
	mkdir "ref"
	faSplit byname $1 "ref"/
	cd "ref"
	for i in *;
	do
	    echo "faToTwoBit $i $(basename -s .fa $i).2bit" >> cmd.sh
	done
	refChr=$(wc -l cmd.sh | cut -d " " -f 1)
	parallel --jobs=$refChr < cmd.sh
	rm cmd.sh *.fa
	
	cd $FOLDER
	mkdir "query"
	faSplit byname $2 "query"/
	cd "query"
	for i in *;
	do
	    echo "faToTwoBit $i $(basename -s .fa $i).2bit" >> cmd.sh
	done
	queryChr=$(wc -l cmd.sh | cut -d " " -f 1)
	parallel --jobs=$queryChr < cmd.sh
	rm cmd.sh *.fa
	
  cd $FOLDER1
	wga $1 $2 $FOLDER 
  rm -rf $FOLDER
fi
