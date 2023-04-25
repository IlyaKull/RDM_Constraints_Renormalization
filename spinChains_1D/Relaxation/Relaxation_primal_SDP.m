function [ Elower ]=Relaxation_primal_SDP(n,D,k0,H,CGmaps,sdpOps)

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
 
%% sdp variables
[rho,omega]=declarePrimalVars(n,d,D,k0);
%% constraints
constraints=[];
% positivity
constraints = [constraints, rho  >= 0];
for l=k0:n
    constraints = [constraints, omega{l} >= 0];
end

% local compatability constraint for rho
constraints= [constraints,...
             pTr(rho,1,[d,d^(k0-1),d]) == pTr(rho,3,[d,d^(k0-1),d]) ];
% normalization 
constraints= [constraints, trace(rho) == 1];

% compatibility constraints k=k0
V0xI=tensor(CGmaps.V0,eye(d)); 
IxV0=tensor(eye(d),CGmaps.V0); 
constraints = [constraints,...
                V0xI * rho * V0xI' == pTr(omega{k0},1,[d,D^2,d]),...
                IxV0 * rho * IxV0' == pTr(omega{k0},3,[d,D^2,d])...
               ];
 
% compatibility constraints k>k0
for k=k0:n-1
    LxI=tensor(CGmaps.L{k},eye(d));
    IxR=tensor(eye(d),CGmaps.R{k});
   constraints = [constraints,...
                 LxI * omega{k} * LxI' == pTr(omega{k+1},1,[d,D^2,d]),...
                 IxR * omega{k} * IxR' == pTr(omega{k+1},3,[d,D^2,d])...
                 ];
end

%% objective function (evaluate 2-body hamiltonian on rho{k0-1}
objectiveFunc=trace( rho * tensor(H,eye(d^(k0-1))) );


%% optimize
optimize(constraints,objectiveFunc,sdpOps);

Elower=value(objectiveFunc);

 
