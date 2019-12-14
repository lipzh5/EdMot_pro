function [ JCscore] = JC( y, cy )
% This is a function to calculate JC-score


if(size(y, 2)==1)
    y=y';
end
if(size(cy, 2)==1)
    cy= cy';
end
JCscore= 0;
pre= sum((y+cy)==2);
summ= sum((y+cy)>0);
if(summ>0)
          JCscore= pre/summ;
end
end

