function [Py,avAsymm] = projectToPSDcone(y,PD)
 
[yCell,numVarsY] = readVecToCell(y,PD.yinds);
Py=y;

asymm=zeros(numVarsY,1);
switch PD.parallelize
    case {'process','thread'}
        parfor i=1:numVarsY
            [yCell{i} ,asymm(i)] = projectOneMatrix(yCell{i}); % overwrite yCell
        end
    case 'single'
        for i=1:numVarsY
            [yCell{i} ,asymm(i)] = projectOneMatrix(yCell{i}); % overwrite yCell
        end
end
for  i=1:numVarsY 
    Py(PD.yinds(i,2):PD.yinds(i,3))=reshape(yCell{i},[],1);
end
avAsymm = mean(asymm);