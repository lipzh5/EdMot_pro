function [label_result,ac,pty,f1,Q,link_number, process_record] = active(selected_type, label_type, maxIter, batch, display_result, n, M, nClass, labels_real)
if nargin == 0
    selected_type = 'entropy-diff'; %random, entropy entropy-diff
    label_type = 'fully-real'; %query, hubs, fully-real fully
end
if nargin < 3
    maxIter = 60;
    batch = 1;
end
if nargin < 5
    display_result = 1;
end


M_source = M;
%show_matrix(M_source,'source');
% options = [];
% options.maxIter = 1000;
% options.alpha = 1;
%[U,V, nmf_iter] = GNMF(M, nClass, zeros(n), options);
[V, nmf_iter] = NMF_ALSE(M, nClass);
% for large network
% option.nRepeat = 20;
% option.maxIter = 1000;
% option.alpha = 0;
% [U_final, V_final, nIter_final, objhistory_final] = GNMF_S(M, nClass, [], option);
edge_visited = zeros(n);
node_visited = zeros(1,n);
link_selected_record = zeros(maxIter,2);
trend = zeros(maxIter,1);

S=convert3(V);
en = entropy_vector(S);
[~,label_result]=max(S,[],2);
label_original = label_result;
%ac(1)=compute_NMI(labels_real,label_result);
ac(1) = NMI(labels_real,label_result);
pty(1) = Purity(labels_real,label_result);
f1(1) = F1Over(labels_real,label_result);
Q(1) = modul(M,label_result);
%disp([datestr(now) '  iter:0; ac:' num2str(ac(1))]);
global link_accu;
link_accu = 0;
link_number(1) = link_accu;



hubs = select_hub(label_result, en, nClass);

if  display_result
        show_matrix(M,num2str(0));
    end

iter = 1;
while iter <= maxIter
    if iter == 15
%        disp(['hehe']);
    end
    %disp([datestr(now) '  iter:' num2str(iter)]);
    if strcmp(selected_type, 'baseline')
        select_process = tic;
        [p1, p2, edge_visited] = select_link_baseline(M, edge_visited, batch);
        time_select = toc(select_process);
    elseif ~strcmp(selected_type, 'random')
        select_process = tic;
        [p1, p2, edge_visited, node_visited] = select_link(selected_type, M, label_result, S, en, edge_visited, node_visited, batch);
        time_select = toc(select_process);
    elseif strcmp(selected_type, 'random')
        [p1, p2, edge_visited] = select_link_random(M, edge_visited, batch);
    end
    if length(p1)==0 || length(p2) == 0
        break;
    end
    %[p1, p2, edge_visited] = select_link_random(M, edge_visited);
    %link_selected_record(iter,:) = [p1, p2];
    M_old = M;
    label_process = tic;
    [M, hub1, hub2, edge_visited] = label_link(label_type, M, p1, p2, hubs, labels_real, label_result, edge_visited);
    time_label = toc(label_process);
    %[U,V, nmf_iter err] = GNMF(M, nClass, zeros(n), options, U, V);
    [V, nmf_iter err] = NMF_ALSE(M, nClass, V);
    
%     option.nRepeat = 1;
%     option.maxIter = 1000;
%     option.alpha = 0;
%     nmf_process = tic;
%     [U_final, V_final, nIter_final, objhistory_final] = GNMF_S(M, nClass, [], option, 1, V_final);
%     time_nmf = toc(nmf_process);
    
    [index_x, index_y] = find((M - M')~=0);
    
    %nmf_iter
    %err
    
    S=convert3(V);
    en = entropy_vector(S);
    [~,label_result]=max(S,[],2);
    %ac(iter+1)=compute_NMI(labels_real,label_result);
    ac(iter+1)=NMI(labels_real,label_result);

    link_number(iter+1)=link_accu;
    trend_temp = ac(iter+1) - ac(iter);
    if trend_temp < 0
        trend_temp = 0;
    end
    if trend_temp < 0
        %pause;
    end
    trend(iter) = trend_temp;
    
    %show_matrix(M,['iter:' num2str(iter) '  trend:' num2str(trend_temp)]);
    %pause
    process_record = [];
%     process_record(iter).label_original = label_original;
%     process_record(iter).M_old = M_old;
%     process_record(iter).M_new = M;
%     process_record(iter).M_change = M - M_old;
%     process_record(iter).V = V_final;
%     process_record(iter).label_result = label_result;
%     process_record(iter).en = en;
%     process_record(iter).ac = ac(iter);
%     process_record(iter).trend = trend_temp;
%     process_record(iter).link_selected = [p1, p2];
%     process_record(iter).link_accu = link_accu;

    if mod(iter, 20) == 0 && display_result
        show_matrix(M,num2str(iter));
    end
    
    
    hubs = select_hub(label_result, en, nClass);
    
    
    if iter > 4
         accu_trend = sum(trend(iter-4:iter));
    else
        accu_trend = sum(trend(1:iter));
    end
    
    %disp([datestr(now) '  iter:' num2str(iter) ' ; ac:' num2str(ac(iter+1)) ' edges:' num2str(sum(sum(M))/2) 'time_select' num2str(time_select) 'time_label' num2str(time_label) 'time_nmf' num2str(time_nmf) 'trend:' num2str(accu_trend)]);

    
%     if accu_trend < 0.001 && iter >5
%         break;
%     end
    
    iter = iter +1;
    
end

if display_result
        show_matrix(M,num2str(iter-1));
    end

M_change = M - M_source;
%show_matrix(M,'end');
%show_matrix(M_change,'change');
edge_add = sum(sum(M_change > 0))/2;
edge_remove = sum(sum(M_change < 0))/2;

%save result.mat
if display_result
    show_result(M_source, M, ac);
end

%figure;plot(ac);
%pause
end

function [point1_set, point2_set, edge_visited] = select_link_baseline(adj, edge_visited, batch)
[x,y,~]=find(edge_visited==0);
edge_number = length(x);
if (edge_number > 0)
    selected_index = randi(edge_number, batch, 1);
    point1_set = x(selected_index);
    point2_set = y(selected_index);
    
    for i = 1:length(point1_set)
        point1 = point1_set(i);
        point2 = point2_set(i);
        edge_visited(point1, point2) = 1;
        edge_visited(point2, point1) = 1;
    end
else
    point1_set = [];
    point2_set = [];
end
end

function [point1_set, point2_set, edge_visited] = select_link_random(adj, edge_visited, batch)
[x,y,~]=find(adj-edge_visited>0);
edge_number = length(x);
if (edge_number > 0)
    selected_index = randi(edge_number, batch, 1);
    point1_set = x(selected_index);
    point2_set = y(selected_index);
    
    for i = 1:length(point1_set)
        point1 = point1_set(i);
        point2 = point2_set(i);
        edge_visited(point1, point2) = 1;
        edge_visited(point2, point1) = 1;
    end
else
    point1_set = [];
    point2_set = [];
end
end

function [point1_set, point2_set, edge_visited, node_visited] = select_link(selected_type,adj, label_result, S, en, edge_visited, node_visited, batch)
if nargin < 8
    batch = 1;
end
node_number = length(en);
edge_number = sum(sum(adj>0));
link_result = zeros(edge_number, 5);
index = 1;
for k=1:node_number
    for j=k+1:node_number
        if edge_visited(k,j)==0 && adj(k,j)~=0
            
            if strcmp(selected_type, 'entropy-diff')
                if label_result(k) ~= label_result(j)
                    link_en = en(k) + en(j);
                    link_diff = S(k) .* S(j)';
                    link_result(index,:) = [k, j, en(k), en(j), link_en];
                    index = index + 1;
                end
            elseif strcmp(selected_type, 'entropy')
                link_en = en(k) + en(j);
                link_diff = S(k) .* S(j)';
                link_result(index,:) = [k, j, en(k), en(j), link_en];
                index = index + 1;
            end
        end
    end
end
link_result = link_result(1:index-1,:);

if size(link_result,1) == 0
    point1_set = [];
    point2_set = [];
else
    link_en_vector = link_result(:,5);
    %[~, idx] = max(link_en_vector);
    [maxvalue, idx] = sort(link_en_vector, 'descend');
    idx = idx(1:min(batch, length(idx)));
    
    p = link_result(idx,1:2);
    point1_set = p(:,1);
    point2_set = p(:,2);
    
    for t = 1:length(point1_set)
        point1 = point1_set(t);
        point2 = point2_set(t);
        node_visited(point1) = node_visited(point1)+1;
        node_visited(point2) = node_visited(point2)+1;
        
        edge_visited(point1, point2) = 1;
        edge_visited(point2, point1) = 1;
    end
end
end

function en = link_entropy(en1 , en2)
en = en1 + en2;
end

function M = label_selected_link(M, p1_set, p2_set, label_real)
global link_accu;
for t = 1: length(p1_set)
end
end

function [M, hub1, hub2, edge_visited] = label_link(label_type, M, p1_set, p2_set, hubs, label_real, label_result, edge_visited,weight, type)
if nargin < 9
    weight = 0;
end
if nargin < 10
    type = 'result';
end
if strcmp(label_type, 'fully-real')
    type = 'real';
end

global link_accu;

guess_correct = 0;
guess_wrong = 0;
hub1 = [];
hub2 = [];

for t = 1: length(p1_set)
    p1 = p1_set(t);
    p2 = p2_set(t);
    %decide whether p1 and p2 in the same community
    if label_real(p1) ~= label_real(p2)
        M(p1, p2) = 0;
        M(p2, p1) = 0;
    else
        M(p1, p2) = 1;
        M(p2, p1) = 1;
    end
    link_accu = link_accu + 1;
    
    if ~strcmp(label_type, 'query')
        hub_selected = 0;
        %connect p1 to a hub
        for hub = hubs'
            if label_real(hub) == label_real(p1)
                M(p1, hub) = 1;
                M(hub, p1) = 1;
                edge_visited(p1, hub) = 1;
                edge_visited(hub, p1) = 1;
                hub_selected = hub;
                link_accu = link_accu + 1;
                break;
            end
        end
        hub1(t) = hub_selected;
        
        %disconect p1 with other neighbour which have different communties with
        %p1's hubs
        if hub_selected  ~= 0 && (strcmp(label_type, 'fully') || strcmp(label_type, 'fully-real'))
            neighbour_p1 = setdiff(find(M(p1,:)>0), p2);
            for q=neighbour_p1
                cmp_real = (label_real(hub) ~= label_real(q));
                cmp_result = (label_result(hub) ~= label_result(q));
                if strcmp(type, 'real')
                    cmp = cmp_real;
                    link_accu = link_accu + 1;
                elseif strcmp(type,'result')
                    cmp = cmp_result;
                    %disp(['disconnect ' num2str(p1) ' and ' num2str(q) ':' num2str(cmp_real == cmp_result)]);
                    guess_correct = guess_correct + (cmp_real == cmp_result);
                    guess_wrong = guess_wrong + ( 1 - (cmp_real == cmp_result));
                end
                
                if cmp
                    M(q, p1)=M(q, p1)*weight;
                    M(p1,q)=M(p1,q)*weight;
                    
                    
                end
            end
        end
    end
    
    if ~strcmp(label_type, 'query')
        %connect p2 to a hub
        hub_selected = 0;
        for hub = hubs'
            if label_real(hub) == label_real(p2)
                M(p2, hub) = 1;
                M(hub, p2) = 1;
                edge_visited(p2, hub) = 1;
                edge_visited(hub, p2) = 1;
                hub_selected = hub;
                link_accu = link_accu + 1;
                break;
            end
        end
        hub2(t) = hub_selected;
        
        %disconect p2 with other neighbour which have different communties with
        %p2's hubs
        if hub_selected  ~= 0 && (strcmp(label_type, 'fully') || strcmp(label_type, 'fully-real'))
            neighbour_p2 = setdiff(find(M(p2,:)>0), p1);
            for q=neighbour_p2
                cmp_real = (label_real(hub) ~= label_real(q));
                cmp_result = (label_result(hub) ~= label_result(q));
                if strcmp(type, 'real')
                    cmp = cmp_real;
                    link_accu = link_accu + 1;
                elseif strcmp(type,'result')
                    cmp = cmp_result;
                    %disp(['disconnect ' num2str(p2) ' and ' num2str(q) ':' num2str(cmp_real == cmp_result)]);
                    guess_correct = guess_correct + (cmp_real == cmp_result);
                    guess_wrong = guess_wrong + ( 1 - (cmp_real == cmp_result));
                end
                
                if cmp
                    M(q, p2)=M(q, p2)*weight;
                    M(p2, q)=M(p2, q)*weight;
                end
            end
        end
    end
end
guess_total = guess_correct + guess_wrong;
if strcmp(type,'result')
    % disp(['correct guess:' num2str(guess_correct) ' ; wrong guess:' num2str(guess_wrong)]);
end
end

function hubs = select_hub( label_result, en, nClass)
hubs_rate = 3;
hubs = [];
for idx = 1:nClass;
    member = (label_result==idx);
    member_en = member .* en;
    v_max = ones(nClass,1)./nClass;
    member_en(member == 0) =  -10*sum(v_max.*log2(v_max));
    [~, I] = sort(member_en);
    hubs = [hubs; I(1:hubs_rate)];
end
end

function S=convert3(X)
X = max(full(X), 1e-10);
b=sum(X,2);
c=size(X,2);
D=repmat(b, [1,c]);
S=X./D;
end

function show_result(source_adj, result_adj,  ac)
figure;
subplot(2,2,1,'Color',[1 1 1]);imagesc(source_adj);
subplot(2,2,2,'Color',[1 1 1]);imagesc(result_adj);
subplot(2,2,3,'Color',[1 1 1]);imagesc(result_adj - source_adj);
subplot(2,2,4);plot(ac);
figure
show_matrix(result_adj - source_adj,'diff')
end



function en = entropy_vector(labelmatrix)
n = size(labelmatrix,1);
en = zeros(n,1);
for i=1:n
    v = labelmatrix(i,:);
    v(v==0) = 1e-10;
    en(i,1) =  -sum(v.*log2(v));
end
end