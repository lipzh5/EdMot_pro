%EdMot: An Edge Enhancement Approach for Motif-aware Community Detection
%Pre-partitioning, can be replaced with other methods
load('data\pol-book.mat');  str = 'SC';
% str, partitioning method, which can be 'Louvain','SC','AP','NMF', etc.
[LCC0, lcc_inds0, ci0, sizes0] = LargestConnectedComponent(sparse(uG)); 
%[LCC0, lcc_inds0, ci0, sizes0] = LargestConnectedComponent(G);
ty_lcc = ty(lcc_inds0);
ks = numel(unique(ty_lcc)); % number of clusters
A = sparse(LCC0);
W = MotifAdjacency(A,'m4');
[LCC, lcc_inds, ci, sizes] = LargestConnectedComponent(W); 
%Pre-partitioning, can be replaced with other methods
%****************
label = partition(LCC,ks,str);
Ulabel = unique(label); %label may be not sequential numbers 
C = numel(unique(label)); %number of communities in the LCC
C_nodes = cell(C,1);
rid = [];
cid = [];
for i=1:C
    idx = find(label==(Ulabel(i)));
    C_nodes{i} = idx;
    if numel(idx)>=2
        perm = nchoosek(idx,2);
        rid = [rid;perm(:,1);perm(:,2)];
        cid = [cid;perm(:,2);perm(:,1)]; 
    end
end
num = numel(rid);
s = ones(num,1);
n = numel(label);
atomA = sparse(rid,cid,s,n,n);
A(lcc_inds,lcc_inds) = A(lcc_inds,lcc_inds)+atomA; %将LCC部分赋值，使其结构原子化，一定会保持原始结构

%Final partitioning, can be repalced 
com = partition(A,ks,str);
if find(com==0)
   com = com+1; 
end
%*********************************
F1 = F1Over(ty_lcc, com); 
JC = JCOver(ty_lcc, com);
nmi_score = NMI(ty_lcc, com);
pty = Purity(ty_lcc,com);
Q = modul(LCC0,com);
%mod=COM.MOD;
result = [nmi_score,pty,F1,Q];