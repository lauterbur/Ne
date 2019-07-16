#!/bin/bash

Ns="10 20 50 100 200 500 1000"
#Ns=10
ns="200 500 1000"

i=$1
j=$2
                        scale=8 
                        THETA=$( bc -l <<<"$i*4*0.000001" )
#                        echo ${THETA}
                        RHO=$( bc -l <<<"$i*4*0.000001/1000" )
                        echo ${i}-${j}-${k}

	                hapn=$((j * 2))
	                printf "./msdir/ms %s 10000 -t %s -r %s 5000" ${hapn} ${THETA} ${RHO} > ${i}-${j}/${i}-${j}-slim.out
        	        printf "\nSLiM %s 10000\n" "${hapn}" >> ${i}-${j}/${i}-${j}-slim.out # tells sample_stats how many iterations
                	printf "SLiM %s 10000\n" "${hapn}"
	                printf "000 001 002\n\n" >> ${i}-${j}/${i}-${j}-slim.out
			for l in {1..10001}
                	do
                        	echo $l
                       		./SLiM/bin/slim ${i}-${j}/${i}-${j}-run_slim.txt
	                        printf "./msdir/ms %s %s -t %s -r %s 5000\n" ${hapn} ${ITERS} ${THETA} ${RHO}
        	        done

	                cat ${i}-${j}/${i}-${j}-slim.out | ./msdir/sample_stats > ${i}-${j}/${i}-${j}-samplestats-slim.txt



