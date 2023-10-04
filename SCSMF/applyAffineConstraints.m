function [y]=applyAffineConstraints(x,PD,PAR,LxI,IxR)
[xCell] = readVecToCell(x,PD.xinds);
numVarsY=PD.yinds(end,1);
yCell=cell(numVarsY,1);

 
% switch PD.parallelize
%     case 'process'
%         LxI=PD.LxIconst.Value;
%         IxR=PD.IxRconst.Value;
%     case {'thread','single'}
%         LxI=PD.LxI;
%         IxR=PD.IxR;
%     otherwise
%         error('PD.parallelize parameter should be eiter process,thread or single')
% end
%  

%% compute output entries
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
    clear xCell

    % because of stupid parfor rules:
        k0m2=PD.k0-2;
        d=PD.d;
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
%% reshape y into a vector
y=nan(PD.yinds(end,3),1); 
for i=1:numVarsY 
    y(PD.yinds(i,2):PD.yinds(i,3))= (-1).*reshape(yCell{i},[],1); %note the -1
end