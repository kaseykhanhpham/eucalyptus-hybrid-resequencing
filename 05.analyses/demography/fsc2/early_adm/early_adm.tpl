//Parameters for the coalescence simulation program : fastsimcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NGLOB$
NCORD$
//Haploid samples sizes 
40
20
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0	0	  	
0  	0
//Migration matrix 1
0       MIG1$
MIG2$    0
//Migration matrix 2
0	0   
0       0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
TADM$ 0 0 0 1 0 1
6000 0 1 1 NANC 0 2 absoluteResize
//Number of independent loci [chromosome] 
11 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 4.80e-7 OUTEXP
