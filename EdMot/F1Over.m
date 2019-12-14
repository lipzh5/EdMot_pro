function [ ff, F1rec, rec, pre ] = F1Over( y, cy )
% This is a function to calculate F1-score

ff=0;
if(size(y, 2)>1 && size(y, 1)>1)
    if(size(cy, 1)==1)
        cy= cy';
    end
    if(size(cy, 2)==1)
        n = size(y,1);
        m = max(cy);
        F1rec= zeros(n,m);
        pre= zeros(n,m);
        rec= zeros(n,m);
        for i=1:n
            for j=1:m
                [F1rec(i,j), rec(i,j), pre(i,j)]= F1(y(i,:),cy==j);
            end
        end
        re= max(F1rec);
        si= zeros(1, m);
        for i=1:m
            si(i)= sum(cy(i,:));
        end
        ff= sum(re.*si)/sum(si);
    end
    if(size(cy, 2)~=1)
        n = size(y,1);
        m = size(cy,1);
        F1rec= zeros(n,m);
        pre= zeros(n,m);
        rec= zeros(n,m);
        for i=1:n
            for j=1:m
                [F1rec(i,j), rec(i,j), pre(i,j)]= F1(y(i,:),cy(j,:));
            end
        end
        re= max(F1rec);
        si= zeros(1, m);
        for i=1:m
            si(i)= sum(cy(i,:));
        end
        ff= sum(re.*si)/sum(si);
    end
end
if(size(y, 1)==1)
    y= y';
end
if(size(cy, 1)==1)
    cy= cy';
end
if(size(y, 2)==1)
        if(size(cy, 2)==1)
            n = max(y);
            m = max(cy);
            F1rec= zeros(n,m);
            pre= zeros(n,m);
            rec= zeros(n,m);
            for i=1:n
                for j=1:m
                    [F1rec(i,j), rec(i,j), pre(i,j)]= F1(y==i,cy==j);
                end
            end
            re= max(F1rec);
            si= zeros(1, m);
            for i=1:m
                si(i)= sum(cy==i);
            end
            ff= sum(re.*si)/sum(si);
        end
        if(size(cy, 2)~=1)
            n = max(y);
            m = size(cy,1);
            F1rec= zeros(n,m);
            pre= zeros(n,m);
            rec= zeros(n,m);
            for i=1:n
                for j=1:m
                    [F1rec(i,j), rec(i,j), pre(i,j)]= F1(y==i,cy(j,:));
                end
            end
            re= max(F1rec);
            si= zeros(1, m);
            for i=1:m
                si(i)= sum(cy(i,:));
            end
            ff= sum(re.*si)/sum(si);
        end
end

end

