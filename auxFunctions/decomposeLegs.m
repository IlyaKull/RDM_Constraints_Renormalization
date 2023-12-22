function [AmultiDim] = decomposeLegs(A2Dim, dimsOut, dimsIn) 
% decompose a many-body operator into a multidimensional array 
% out legs: d1 d2 d3 d4 (dimsOut)
%           |  |  |  |
%           AAAAAAAAAA
%            |  |  |
% in legs:   D1 D2 D3   (dimsIn)
%
sizeA = size(A2Dim);
assert(length(sizeA)==2,'Input must be 2-dim')
assert(prod(dimsIn)==sizeA(2),'dimsIn does not match size(A)(2)')
assert(prod(dimsOut)==sizeA(1),'dimsOut does not match size(A)(1)')

AmultiDim = reshape(A2Dim,[dimsOut(end:-1:1),dimsIn(end:-1:1)]);
AmultiDim = permute(AmultiDim,[length(dimsOut):-1:1, length(dimsOut) + (length(dimsIn):-1:1)]);


% to do the reverse operation: from multilidim array to matrix
% reverse the order of the 'Out' and 'In' indices:
% A2Dim = reshape(permute(AmultiDim,[ndOut:-1:1, ndOut +(ndIn:-1:1)]),prod(dimsOut),prod(dimsIn))

% here is a test
% Lo = rand([2 1]);
% Oo = rand([3 1]);
% Ro = rand([4 1]);
% dimOut=2*3*4;
% Li = rand([5 1]);
% Oi = rand([6 1]);
% Ri = rand([7 1]);
% dimIn=5*6*7;
% 
% 
% Ti=tensor(Li,Oi,Ri);
% To=tensor(Lo,Oo,Ro);
% T=To*Ti.';
% Td=decomposeLegs(T,[2 3 4],[5 6 7]);
% 
% Tdr=reshape(permute(Td,[3 2 1, 6 5 4]),dimOut,dimIn);
% diffff=sum(abs(Tdr-T),'all')
% 
% Ti2=tensor(Li,Ri)*(Oi'*Oi);
% 
% T2=To*Ti2.';
% TxOi = ncon({Td,Oi},{[-1 -2 -3 -4 1 -5],[1]});
% TxOi2D=reshape(permute(TxOi,[3 2 1 5 4]),dimOut,dimIn/6);
% difff= norm( TxOi2D - T2)
% 
% Ti3=tensor(Li,Oi)*(Ri'*Ri);
% 
% T3=To*Ti3.';
% TxRi = ncon({Td,Ri},{[-1 -2 -3 -4 -5 1],[1]});
% TxRi2D=reshape(permute(TxRi,[3 2 1 5 4]),dimOut,dimIn/7);
% difff= norm( TxRi2D - T3)
