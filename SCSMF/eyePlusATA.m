function [x]=eyePlusATA(x,PD,PAR,LxI,IxR,inds_L,inds_R)
[xCell,numVarsX] = readVecToCell(x,PD.xinds);
numVarsY=PD.yinds(end,1);

yCell=cell(numVarsY,1);
k0m2=PD.k0-2;
d=PD.d;
%% apply A

yCell{1}=  speye(PD.yinds(1,4)) * xCell{1} + ...
                    + kron(speye(PD.d),xCell{2}) - kron(xCell{2},speye(PD.d)) + ...
                     PD.V0xI'*xCell{3}*PD.V0xI + PD.IxV0'*xCell{4}*PD.IxV0 ;
if PAR
    % prepare copies of gamma vars
    k=PD.k0:PD.n-1;
    gammaLk(k)=xCell(2*(k-PD.k0)+5);
    gammaRk(k)=xCell(2*(k-PD.k0)+1+5);
    gammaLkMinusOne(k)=xCell(2*(k-PD.k0)+5-2);
    gammaRkMinusOne(k)=xCell(2*(k-PD.k0)+1+5-2);
    
    switch PD.parallelize
        case 'process'
            parfor k=PD.k0:PD.n-1
                yCell{k-k0m2} =  LxI.Value{k}' * gammaLk{k} * LxI.Value{k} + ...
                            IxR.Value{k}' * gammaRk{k} * IxR.Value{k}...
                            - kron(speye(d),gammaLkMinusOne{k}) ...
                            - kron(gammaRkMinusOne{k},speye(d));
            end
        case 'thread'
            parfor k=PD.k0:PD.n-1
                yCell{k-k0m2} =  LxI{k}' * gammaLk{k} * LxI{k} + ...
                            IxR{k}' * gammaRk{k} * IxR{k}...
                            - kron(speye(d),gammaLkMinusOne{k}) ...
                            - kron(gammaRkMinusOne{k},speye(d));
            end
     end
    
    yCell{numVarsY} =  - kron(speye(PD.d),gammaLk{PD.n-1}) ...
                       - kron(gammaRk{PD.n-1},speye(PD.d));
else %single process
    for k=PD.k0:PD.n-1
        yCell{k-PD.k0+2} =  LxI{k}' * xCell{2*(k-PD.k0)+5} * LxI{k} + ...         %gammaL{k0} is stored at x{2*(k-k0)+5}
                    IxR{k}' * xCell{2*(k-PD.k0)+1+5} * IxR{k}...               %gammaR{k0} is stored at x{2*(k-k0)+1+5}
                    - kron(speye(PD.d),xCell{2*(k-PD.k0)+5-2}) ...                
                    - kron(xCell{2*(k-PD.k0)+1+5-2} ,speye(PD.d));                
    end
    
    yCell{numVarsY} =  - kron(speye(PD.d), xCell{2*(PD.n-1-PD.k0)+5}) ...
                       - kron(xCell{2*(PD.n-1-PD.k0)+1+5},speye(PD.d));
end
%% apply A^transpose


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

                    xCellOdd{k}= IxR.Value{k} * yCell{k-k0m2}* IxR.Value{k}'  ...
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
for i=1:numVarsX 
    x(PD.xinds(i,2):PD.xinds(i,3))= x(PD.xinds(i,2):PD.xinds(i,3)) + reshape(xCell{i},[],1);  
end

