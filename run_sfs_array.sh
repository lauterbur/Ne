#!/bin/bash

Ns="10 20 50 100 200 500 1000"
#Ns=10
#ns="2 5 10 20 50 100"
ns="200 500 1000"
#ns="10"
#Ts="1 10 100 1000"

#for i in $Ns
#do
#        for j in $ns
#        do
#		for k in $Ts
#		do
i=$1
j=$2
#k=$3
	                scale=8 
			THETA=$( bc -l <<<"$i*4*0.000001" )
#               	 echo ${THETA}
               		RHO=$( bc -l <<<"$i*4*0.000001/1000" )
			echo ${i}-${j}
#			TIME=$( bc -l <<<"$k/2000" )
#			# TIME is relative to the ancestral N (*2 for diploid), number of generations after event (bottleneck) to sample
			RELSIZE=$( bc -l <<<"$i/1000" )
			# RELSIZE is size change between ancestral and new N
#			./SFS_code/bin/sfs_code 1 10001 -A -I -a N -n ${j} -N ${i} -t ${THETA} -r ${RHO} -L 1 5000 > ${i}-${j}/${i}-${j}-sfs.out
			./SFS_code/bin/convertSFS_CODE ${i}-${j}/${i}-${j}-sfs.out --ms T 0 F ${i}-${j}/${i}-${j}-sfs-msformat.txt
	                cat ${i}-${j}/${i}-${j}-sfs-msformat.txt | ./msdir/sample_stats > ${i}-${j}/${i}-${j}-samplestats-sfs.txt
#                done
#        done
#done

