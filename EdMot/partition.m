function labels = partition(A,ks,str,~)
%partitioning using traditional mrthods
%Input: A, the adjacency matrix;
%       ks, number of clusters;
%       str, partitioning method, which can be 'Louvain','NMF', etc.
if strcmp(str,'Louvain')
    [COM,~] = cluster_jl(A);
    labels = (COM.COM{1,1})';    
end

if strcmp(str,'AP')
    [i,j,s] = find(A);
    S = [i,j,s];
    p = median(s); 
    [labels,~,~,~]=apcluster(S,p);
end

if strcmp(str,'Ncut')
    [NcutDiscrete,~,~] = ncutW(A,ks);
    [~,labels] = find(NcutDiscrete);
end


 if strcmp(str,'SC')
     [labels, ~] = KmeansCluster(A, ks);
 end
 
if strcmp(str,'NMF')
     labels = NMF_getLabel(A,ks);
end

end

