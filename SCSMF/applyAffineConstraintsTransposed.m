function [x]=applyAffineConstraintsTransposed(y,PD,PAR,LxI,IxR,inds_L,inds_R)
 
[yCell] = readVecToCell(y,PD.yinds);

%% compute NEGATIVE OF output entries
numVarsX=PD.xinds(end,1);
k0m2=PD.k0-2; %index shift constant for use in parfor
% initialize xCell (final output)
xCell=cell(numVarsX,1);
for i=1:numVarsX
   PD.xCell{i}=nan(PD.xinds(i,4));
end

xCell{1}=  trace(yCell{1});
xCell{2}= pTr_Inds(yCell{1},PD.inds_L0)...
         -pTr_Inds(yCell{1},PD.inds_R0);
     
xCell{3}= PD.V0xI*yCell{1}*PD.V0xI' ...
        - pTr_Inds(yCell{2},inds_L);
    
xCell{4}= PD.IxV0*yCell{1}* PD.IxV0' ...
         -pTr_Inds(yCell{2},inds_R);

if PAR %parfor version 

    xCellEven=cell(PD.n-1,1);
    xCellOdd=cell(PD.n-1,1);
    for i=5:numVarsX
        xCellEven{i}=nan(PD.xinds(i,4));
        xCellOdd{i}=nan(PD.xinds(i,4));
    end

    % we need two copies of y, one shifted by one for the parallelization
    yCell_kPlusOne=cell(PD.n-1,1);
    for k=PD.k0:PD.n-1
        yCell_kPlusOne{k-k0m2}=yCell{k-k0m2+1};
    end
    switch PD.parallelize
        case 'process'
            parfor k=PD.k0:PD.n-1 
                    xCellEven{k}= LxI.Value{k} * yCell{k-k0m2} * LxI.Value{k}' ...
                                  -pTr_Inds(yCell_kPlusOne{k-k0m2},inds_L);

                    xCellOdd{k}= IxR.Value{k} * yCell{k-k0m2} * IxR.Value{k}'  ...
                                    -pTr_Inds(yCell_kPlusOne{k-k0m2},inds_R);
            end
        case 'thread'
            parfor k=PD.k0:PD.n-1 
                    xCellEven{k}= LxI{k} * yCell{k-k0m2} * LxI{k}' ...
                                  -pTr_Inds(yCell_kPlusOne{k-k0m2},inds_L);

                    xCellOdd{k}= IxR{k} * yCell{k-k0m2}* IxR{k}'  ...
                                    -pTr_Inds(yCell_kPlusOne{k-k0m2},inds_R);
            end
    end
    
    for k=PD.k0:PD.n-1
        xCell([5+2*(k-PD.k0), (5+2*(k-PD.k0)+1)]) = {xCellEven{k};xCellOdd{k}};
    end
    clear xCellEven xCellOdd

else %single thread version
    for k=PD.k0:PD.n-1
        xCell{5+2*(k-PD.k0)} = LxI{k} * yCell{k-k0m2} * LxI{k}' ...
                          -pTr_Inds(yCell{k-k0m2+1},inds_L);

        xCell{5+2*(k-PD.k0)+1} = IxR{k} * yCell{k-k0m2}* IxR{k}'  ...
                            -pTr_Inds(yCell{k-k0m2+1},inds_R);
    end
end
    
%% reshape x into a vector
x=nan(PD.xinds(end,3),1); 
for i=1:numVarsX 
    x(PD.xinds(i,2):PD.xinds(i,3))=(-1).* reshape(xCell{i},[],1); %note the -1! 
end