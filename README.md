# Stable_Nanopores_Generation-
All possible stable nanopores are created using Redelmeir's algorithm and based on their geometrical properties searching is also enabled

stable_nanopores_parallel.m generates a .mat file that contains data of all possible stable nanopores of a particular size, n. 
'n' is the number of triangles in the polyiamond or number of atoms that are missing in the nanopore.
It generates spreadsheets that sort and record geometrical properties (major axis, minor axis, shape factor) of these stable nanopores.
It also creates an xyz file that saves the coordinates of the nanopores.
Apart from these the validity of the algorithm can be checked using the code itself by generating fixed and free polyiamonds and matching their number with the OEIS data.
All the other .m files given in the main folder are functions that are used by stable_nanopores.m
The word parallel denotes to the use of parallel computing. The program requires Mapping toolbox and Parallel Computing toolbox from MATLAB. 

The Supplementary_Codes folder contains program files that can be used to extract and visualize further details of the data generated by stable_nanopores_parallel.m
figs_min_max.m creates figures of nanopore rims with maximum and minimum properties along with the properties.
min_max_xyz creates xyz files of nanopores with maximum and minimum properties. 
The functions min_max_plotrimedge, xyz_A.m and xyz_num.m are improvised versions of functions from the main folder.
Polyiamond_polyhex.m plots polyiamonds and polyhexes in figures and helps in their visualization.
Before using the files in Supplementary_Codes folder, make sure that the folder contains the data generated by stable_nanopores_parallel.m