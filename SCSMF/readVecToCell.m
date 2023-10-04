function [outCell,numVars] = readVecToCell(varVec,varInds)
numVars=varInds(end,1);
outCell=cell(numVars,1);
for i=1:numVars 
    readInds=varInds(i,2):varInds(i,3);
    outCell{i}=reshape(varVec(readInds), varInds(i,4),varInds(i,4));
end