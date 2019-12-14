function [ ff, F1rec, rec, pre ] = JCOver( y, cy )
% This is a function to calculate JC-score

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
                JCrec(i,j)= JC(y(i,:),cy==j);
            end
        end
        re1= max(JCrec);
        re2= max(JCrec,[],2);
        ff= (sum(re1,2)/size(re1,2)+sum(re2,1)/size(re2,1))/2;
    end
    if(size(cy, 2)~=1)
        n = size(y,1);
        m = size(cy,1);
        F1rec= zeros(n,m);
        pre= zeros(n,m);
        rec= zeros(n,m);
        for i=1:n
            for j=1:m
                JCrec(i,j)= JC(y(i,:),cy(j,:));
            end
        end
        re1= max(JCrec);
        re2= max(JCrec,[],2);
        ff= (sum(re1,2)/size(re1,2)+sum(re2,1)/size(re2,1))/2;
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
                    JCrec(i,j)= JC(y==i,cy==j);
                end
            end
            re1= max(JCrec);
            re2= max(JCrec,[],2);
            ff= (sum(re1,2)/size(re1,2)+sum(re2,1)/size(re2,1))/2;
        end
        if(size(cy, 2)~=1)
            n = max(y);
            m = size(cy,1);
            F1rec= zeros(n,m);
            pre= zeros(n,m);
            rec= zeros(n,m);
            for i=1:n
                for j=1:m
                    JCrec(i,j)= JC(y==i,cy(j,:));
                end
            end
            re1= max(JCrec);
            re2= max(JCrec,[],2);
            ff= (sum(re1,2)/size(re1,2)+sum(re2,1)/size(re2,1))/2;
        end
end

end