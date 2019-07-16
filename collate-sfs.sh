#!/bin/bash

Ns="10 20 50 100 200 500 1000"
# 2000 5000 10000"
#Ns="10"
#ns="2 5 10 20 50 100"
ns="200 500 1000"
#ns="10"

printf "Ne,n,ss,pi\n" > collate-all-samplestats-sfs.csv

for i in $Ns
do
        hapN=$((i * 2))
        for j in $ns
        do
                OUT="${i}-${j}"
#               echo ${OUT}
                lines=`wc -l ${OUT}/${OUT}-samplestats-sfs.txt | awk '{print $1; exit}'`
                segsites=`cat ${OUT}/${OUT}-samplestats-sfs.txt | cut -f 4`
                segsites=( $segsites )
                pisites=`cat ${OUT}/${OUT}-samplestats-sfs.txt | cut -f 2`
                pisites=( $pisites )
                for l in $(seq 1 $lines)
                do
#                       ss=`echo ${segsites} | awk -v var=$l '{print $var}'`
                        ss=`echo ${segsites[$l]}`
#                       echo ${segsites}
#                       echo ${ss}
#                       pi=`echo ${pisites} | awk -v var=$l '{print $var}'`
                        pi=`echo ${pisites[$l]}`
#                       echo ${pisites}
#                       echo ${pi}
                        printf "${i},${j},${ss},${pi}\n" >> collate-all-samplestats-sfs.csv
#                       echo "${i},${j},${ss},${pi}\n"
                done
        done
done

