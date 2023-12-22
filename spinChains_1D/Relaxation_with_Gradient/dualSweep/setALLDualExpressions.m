function [expr] = setALLDualExpressions(d,n,k0,H,CGmaps,vars) 
 
if isstruct(H)
    suppH=H.suppH; % the number of sites H acts on 
    d = H.d; % the dimention of each site
    H=H.Hamiltonian; % the Hamiltonian itself
else % H is a d^2 x d^2 matrix
    suppH=2;
    d=sqrt(size(H,1));
end
assert(size(H,1) == d^suppH,'dimentinos of H inconsisten with support H and d')
assert(k0+1 >= suppH,'k0+1 must be >= support H')


expr=cell(n,1);
 
V0xI=tensor(CGmaps.V0,eye(d)); 
IxV0=tensor(eye(d),CGmaps.V0); 

% expr{k0-1} vars: gamma{k0-1},alpha
expr{k0-1} = tensor(H,eye(d^(k0-(suppH-1)))) ...
        + V0xI' * vars.gammaL{k0-1} * V0xI ...
        + IxV0' * vars.gammaR{k0-1} * IxV0 ...
        + tensor(eye(d),vars.alpha) ...
        - tensor(vars.alpha,eye(d)) ;
        

LxI=tensor(CGmaps.L{k0},eye(d));
IxR=tensor(eye(d),CGmaps.R{k0});

% expr{k0} vars: gamma{k0}, gamma{k0-1}
expr{k0} = LxI'*vars.gammaL{k0}*LxI + IxR'*vars.gammaR{k0}*IxR ...
         - tensor(eye(d),vars.gammaL{k0-1}) - tensor(vars.gammaR{k0-1},eye(d)) ;

% compatibility constraints k>k0
for k=k0+1:n-1
    LxI=tensor(CGmaps.L{k},eye(d));
    IxR=tensor(eye(d),CGmaps.R{k});

    % expr{k} vars: gamma{k}, gamma{k-1}
    expr{k} = LxI'*vars.gammaL{k}*LxI + IxR'*vars.gammaR{k}*IxR ...
            - tensor(eye(d),vars.gammaL{k-1}) - tensor(vars.gammaR{k-1},eye(d));
end

% expr{n} vars: gamma{n-1}
expr{n} = -tensor(eye(d),vars.gammaL{n-1}) - tensor(vars.gammaR{n-1},eye(d)) ;


 