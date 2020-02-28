#!/bin/bash

# Matthew Lee
# Winter 2020
# JagResearch
# MultiMutant

###### ARGUMENTS ######

# $1 = location of pipeline
# $2 = number of instances to create
# $3 = {pdbID}.{chainID}

###### NOTES ######

# this script is intended to be used inside the output
# folder of a successful invocation of multiMutant

###### SCRIPT ######

#list of subdirectories containing output of multiMutant
directories=$(ls -d */)
#total number of pipelines created
progress = 0
#number of pipelines to be running at once
target = $2

while ${#directories[@]} > progress;
do
	for progress < target;
	do
		#copy pipeline into subdirectory
		cp $1 directories[progress]
		
		#cd into subdirectory
		cd directories[progress]

		#copy input files into pipeline input folder
		cp *_em.pdb WT_C/data/$3.cur.in
		#start pipeline in separate process
		./invokePipeline_WT.sh #todo

		#then in that process copy output to top level directory and delete pipeline

		#return to parent dir
		cd ..
	done

	target *= 2

	#wait 1 minute
done


# have cleanup function that checks after x minutes to see if output is in a folder where
# Kinari pipeline was run, if so, then delete the pipeline from that folder

