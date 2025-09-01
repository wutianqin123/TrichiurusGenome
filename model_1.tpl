//Number of population samples (demes)
4 samples to simulate
//Population effective sizes (number of genes) 
NPOP0$
NPOP1$
NPOP2$
NPOP3$
//Sample sizes
20
9
8
11
//Growth rates: negative growth implies populations expansion
0
0
0
0
//Number of migration matrices: 0 implies no migration between demes
2
//migration matrix 0
0.000	M01$	M02$	0.000
M10$	0.000	0.000	0.000
M20$	0.000	0.000	0.000
0.000	0.000	0.000	0.000
//migration matrix 1
0	0	0	0
0	0	0	0
0	0	0	0
0	0	0	0
//historical event: time , source, sink, migrants, new size, new growth rate, migr. matrix
3 historical event
TSP$	2	1	1	1	0	1
TSP1$	1	0	1	1	0	1
TSP2$	0	3	1	1	0	1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut. rate + optional parameters
FREQ    1       0       2.5e-8
