#!/bin/bash

Ns="10 20 50 100 200 500 1000"
#Ns="10"
#ns="2 5 10 20 50 100"
ns="200 500 1000"
#ns="2"
for i in $Ns
do
	BURN=$( bc -l <<<"$i"*10 )
        for j in $ns
        do
                mkdir ${i}-${j}
		echo "initialize() {
        initializeMutationRate(1e-6); //mutation rate (vary this?)
        initializeMutationType('m1', 0.5, 'f', 0.0); // single mutation, co-dominant, no selection
        initializeGenomicElementType('g1', m1, 1.0); // single type of chromosomal region (in this case, non-coding because m1 not selected on and hopefully that's where SNPs are)
        initializeGenomicElement(g1, 0, 4999); //5000 bp regions
        initializeRecombinationRate(1e-9);
}
1 {
        sim.addSubpop('p1', ${i}); // single population of Ns diploid individuals
}

$BURN late() { // after 5*2N burn-in generations, sample
//        cat('generation: ' + sim.generation + '\n');
//        cat(sim.mutations.size() + ' mutations in population\n');
        allIndividuals = sim.subpopulations.individuals;
	if (${j} > ${i}) // if more samples than individuals
	       sampledIndividuals = sample(allIndividuals, ${j}, T); // ns individuals, sampled with replacement
	else
		sampledIndividuals = sample(allIndividuals, ${j}, F); // ns individuals, sampled without replacement
        g = sampledIndividuals.genomes;
//        m = sortBy(unique(g.mutations), 'position');
//        cat(size(m) + ' segsites\n');
//        if (size(m) >= 1)
//          cat('bam');
        sampledIndividuals.genomes.outputMS(filePath='${i}-${j}/${i}-${j}-slim.out', append = T); // MS-style output
        sim.simulationFinished();
}" > "${i}-${j}/${i}-${j}-run_slim.txt"
	done
done
