function [B,ind_sel] = topK(A,K)
%Choose top-K components to pre-partition
% Input: A, the matrix; K, the parameter
% Output: B, the selected components from A
%         ind_sel, selected inds from A
[ci, sizes] = components(A);
[size_sorted, I] = sort(sizes,'descend');
II = I(1:K); %select top K 
idx = [];
for i=1:K
    tmp = find(ci==II(i));
    idx = [idx;tmp];
end

B = A(idx,idx);
ind_sel = idx;
end

