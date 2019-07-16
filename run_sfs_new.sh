#!/bin/bash

Ns="10 20 50 100 200 500 1000"
#Ns=2000 5000 10000"
#Ns=10
#ns="2 5 10 20 50 100"
ns="200 500 1000"
#ns="10"

for i in $Ns
do
        for j in $ns
        do
                scale=8 
		THETA=$( bc -l <<<"$i*4*0.000001" )
#                echo ${THETA}
                RHO=$( bc -l <<<"$i*4*0.000001/1000" )
                SAMPLES=$( bc -l <<<"$ns*2" )
		echo ${i}-${j}
                ./SFS_Code/bin/sfs_code 1 10001 -A -I -a N -n ${j} -N ${i} -t ${THETA} -r ${RHO} -L 1 5000 > ${i}-${j}/${i}-${j}-sfs.out
#                ITERS=0
#                while [ ${ITERS} -lt 10001 ]
#                do
#                        printf "\n" >> ${i}-${j}/${i}-${j}-sfs.out
#	                ./SFS_Code/bin/convertSFS_CODE ${i}-${j}/${i}-${j}-sfs.out --ms T 0 F ${i}-${j}/${i}-${j}-sfs-msformat.txt
#	                cat ${i}-${j}/${i}-${j}-sfs-msformat.txt | ./msdir/sample_stats > ${i}-${j}/${i}-${j}-samplestats-sfs.txt
#                        ITERS="$(grep 'positions' ${i}-${j}/${i}-${j}-sfs-msformat.txt | wc -l)"
#                        echo ${ITERS}
#			echo ${i}-${j}
 #                       ./SFS_Code/bin/sfs_code 1 10001 -A -a N -n ${j} -N ${i} -t ${THETA} -r ${RHO} -L 1 5000 > ${i}-${j}/${i}-${j}-sfs-temp.out
#			tail -n +3 ${i}-${j}/${i}-${j}-sfs-temp.out >> ${i}-${j}/${i}-${j}-sfs.out
#                done
        done
done
