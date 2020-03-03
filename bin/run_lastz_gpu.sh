#!/bin/bash

if [ $# -eq 1 ] || [ $# -eq 0 ] 
then
	wga $1
else
	refPath=$(readlink -f $1)
	queryPath=$(readlink -f $2)
	ref=$(basename -s .fa $1)
	query=$(basename -s .fa $2)
	
	FOLDER1=$PWD/output_$RANDOM/
	FOLDER=$FOLDER1/data_${ref}_${query}_$RANDOM/
	mkdir $FOLDER1
	mkdir $FOLDER
	
	cd $FOLDER
	mkdir "ref"
	faSplit byname $refPath "ref"/
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
	faSplit byname $queryPath "query"/
	cd "query"
	for i in *;
	do
	    echo "faToTwoBit $i $(basename -s .fa $i).2bit" >> cmd.sh
	done
	queryChr=$(wc -l cmd.sh | cut -d " " -f 1)
	parallel --jobs=$queryChr < cmd.sh
	rm cmd.sh *.fa
	
  cd $FOLDER1
  if [ $# -eq 3 ]
  then
      wga $refPath $queryPath $FOLDER $3 
  else
      wga $refPath $queryPath $FOLDER  
  fi

  rm *.segments
  rm -rf $FOLDER
fi
