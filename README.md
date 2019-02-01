# EdMot_pro: Edge enhancement approach for motif-aware community detection 
This is the Matlab (protected) program of EdMot.
You are suggested to use a 64 bit Matlab with version higher than R2014b. Thank you.

Usage:
result = EdMot(uG,str,Nrun,ty)

% Input: uG, node adjacency matrix, in sparse format
%        str, name of the graph partitioning ('Louvain', 'SC','AP',NMF);
%        Nrun, number of runs
%        ty, true cluster labels
% Output：result,consisting of NMI, F-score and modularity as well as the corresponding standard deviations.
