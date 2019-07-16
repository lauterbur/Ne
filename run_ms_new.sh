#!/bin/bash

Ns="10 20 50 100 200 500 1000"
#ns="2 5 10 20 50 100"
ns="200 500 1000"
#Ns="10"
#ns="2"

for i in $Ns
do
        for j in $ns
        do
		THETA=$( bc -l <<<"$i*4*0.000001*5000" )
#		echo ${THETA}
		RHO=$( bc -l <<<"$i*4*0.000001*5000/1000" )
#		echo ${RHO}
		SAMPLES=$( bc -l <<<"$j*2" )
		./msdir/ms ${SAMPLES} 10001 -t ${THETA} -r ${RHO} 5000 > ${i}-${j}/${i}-${j}-ms-temp.out
		printf './msdir/ms %s 10001 -t %s -r %s 5000' ${SAMPLES} ${THETA} ${RHO}
#                ITERS=0
#                while [ ${ITERS} -lt 10001 ]
#                do
#                        ITERS="$(grep 'positions' ${i}-${j}/${i}-${j}-ms-temp.out | wc -l)"
#                        echo ${ITERS}
#                        printf "\n" >> ${i}-${j}/${i}-${j}-ms-temp.out
#			./msdir/ms ${SAMPLES} 1000 -t ${THETA} -r ${RHO} 5000 >> ${i}-${j}/${i}-${j}-ms-temp.out
#			wc -l ${i}-${j}/${i}-${j}-samplestats-ms.txt
		LINES="$(grep '//' ${i}-${j}/${i}-${j}-ms-temp.out | wc -l)"
		printf './msdir/ms %s %s -t %s -r %s 5000\n' ${SAMPLES} ${LINES} ${THETA} ${RHO} > ${i}-${j}/${i}-${j}-ms.out
		tail -n +2 ${i}-${j}/${i}-${j}-ms-temp.out >> ${i}-${j}/${i}-${j}-ms.out
		cat ${i}-${j}/${i}-${j}-ms.out | ./msdir/sample_stats > ${i}-${j}/${i}-${j}-samplestats-ms.txt
                done
	done
done
