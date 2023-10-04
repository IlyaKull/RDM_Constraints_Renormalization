function [p] = pTr_Inds(x,inds)

p=zeros(inds(end)/size(inds,2));

for j=1:size(inds,2)
    p = p + x(inds(:,j),inds(:,j));
end
end