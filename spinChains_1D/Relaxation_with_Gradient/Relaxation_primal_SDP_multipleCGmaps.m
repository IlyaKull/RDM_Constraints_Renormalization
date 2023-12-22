function [ Elower ]=Relaxation_primal_SDP_multipleCGmaps(n,D,k0,H,CGmaps,sdpOps)

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

if k0>n 
   Elower=nan;
   warning('%d= k0 > n = %d, no coarse-graining possible for D=%d and n=%d. aborting!!!',k0,n,D,n);
   return
end


% when multiple CGmaps (sets of) are provided, rho is constrained by
% multiple chains of omegas
OmegaS = cell(numel(CGmaps),1);
constraints=[];

for j=1:numel(CGmaps)
    
    if j==1
        [rho,OmegaS{j}]=declarePrimalVars(n,d,D,k0);
        constraints = [constraints, rho  >= 0];
        constraints= [constraints,...
                 pTr(rho,1,[d,d^(k0-1),d]) == pTr(rho,3,[d,d^(k0-1),d]) ];
        constraints= [constraints, trace(rho) == 1];
    else
        [~,OmegaS{j}]=declarePrimalVars(n,d,D,k0);
    end
    
    for l=k0:n
        constraints = [constraints, OmegaS{j}{l} >= 0];
    end


    % compatibility constraints k=k0
    V0xI=tensor(CGmaps(j).V0,eye(d)); 
    IxV0=tensor(eye(d),CGmaps(j).V0); 
    constraints = [constraints,...
                    V0xI * rho * V0xI' == pTr(OmegaS{j}{k0},1,[d,D^2,d]),...
                    IxV0 * rho * IxV0' == pTr(OmegaS{j}{k0},3,[d,D^2,d])...
                   ];

    % compatibility constraints k>k0
    for k=k0:n-1
        LxI=tensor(CGmaps(j).L{k},eye(d));
        IxR=tensor(eye(d),CGmaps(j).R{k});
       constraints = [constraints,...
                     LxI * OmegaS{j}{k} * LxI' == pTr(OmegaS{j}{k+1},1,[d,D^2,d]),...
                     IxR * OmegaS{j}{k} * IxR' == pTr(OmegaS{j}{k+1},3,[d,D^2,d])...
                     ];
    end
end
%% objective function (evaluate 2-body hamiltonian on rho{k0-1}
objectiveFunc=trace( rho * tensor(H,eye(d^(k0-1))) );


%% optimize
optimize(constraints,objectiveFunc,sdpOps);

Elower=value(objectiveFunc);

 
