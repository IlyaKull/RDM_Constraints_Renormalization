function [styles]= makeLineStyles(colors,markers,lines)

if nargin==0
    colors='kbrgmc';
    markers='ox+d^p';
    lines={'-','--','-.'};
end

styles={};
for i=1:numel(markers)
    for j=1:numel(colors)
        for k=1:numel(lines)
            styles=[styles,[markers(i),colors(j),lines{k}]];
        end
    end 
end

styles = styles(randperm(numel(styles)));