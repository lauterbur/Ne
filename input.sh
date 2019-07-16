#!/bin/bash

Ns="10 20 50 100 200 500 1000"
# 2000 5000 10000"
#Ns="10"
#ns="2 5 10 20 50 100"
ns="200 500 1000"
#ns="10"
#Ts="1 10 100 1000"

#printf "Ne,n,ss,pi\n" > collate-all-samplestats-ms.csv

for i in $Ns
do
        for j in $ns
        do
#                for k in $Ts
#		do
		printf "%d %d\n" $i $j >> input.txt
#                done
        done
done

