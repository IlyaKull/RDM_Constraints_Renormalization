function [S,GFX,tableOut] = plotFromDotOutFiles(modelStr,GFX,DlinSpec,Sin,tableIn)
if nargin < 4 
    [S] = retrieveDataFromOUTfile(modelStr);
else
    S=Sin;
end

if nargin <5
    tableOut=[];
else
    tableOut = [tableIn,nan(size(tableIn,1),1)];
end

switch modelStr
    case 'heisenberg3'
        d=3;
    otherwise
        d=2;
end
DsInLegend=[];
for j=1:size(S,2)
%      mrkSize = 5+ (S(j).maxIter - 2000)/1000;
    mrkSize=10;
   sE= plot(GFX.axes.energy,S(j).n+2,S(j).deltaE,[DlinSpec{S(j).D},'*-.'],'MarkerSize',mrkSize,'MarkerFaceColor',DlinSpec{S(j).D});
     tableOut=[tableOut; real([S(j).n+2,S(j).D,S(j).deltaE]),5];
   
   sV= plot(GFX.axes.nVars, numVars(d,S(j).n,S(j).D),S(j).deltaE ,[DlinSpec{S(j).D},'*-.'],'MarkerSize',mrkSize,'MarkerFaceColor',DlinSpec{S(j).D});
   
   if ~ismember( S(j).D,DsInLegend)
    GFX.legends.energy = [GFX.legends.energy, sprintf('SCSMF D=%d',S(j).D)]; %,S(j).maxIter,S(j).resid(1))];
    DsInLegend=[DsInLegend,S(j).D];
    GFX.subsets.energy=[GFX.subsets.energy,sE];
   end
       
end
