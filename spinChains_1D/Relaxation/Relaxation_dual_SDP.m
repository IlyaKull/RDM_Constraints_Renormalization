function [ Elower,sdpVars,solInfo] = Relaxation_dual_SDP(n,D,k0,H,CGmaps, sdpOps)

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

%% 
sdpVars = declareDualVars(n,d,D,k0);
 
[constraints] = setDualConstraints(d,n,k0,H,CGmaps,sdpVars,suppH);

objectiveFunc = - sdpVars.epsil; %minus sign for maximization

solInfo = optimize(constraints,objectiveFunc,sdpOps);
Elower = -value(objectiveFunc);
  
end

 
