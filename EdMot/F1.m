function [ F1score ,rec, pre] = F1( y, cy )
% This is a function to calculate F1-score

if(size(y, 2)==1)
    y=y';
end
if(size(cy, 2)==1)
    cy= cy';
end
F1score= 0;
pre= sum((y+cy)==2);
rec= pre/sum(cy);
pre= pre/sum(y);
if(pre+rec>0)
          F1score= 2*pre*rec/(pre+rec);
end
end

