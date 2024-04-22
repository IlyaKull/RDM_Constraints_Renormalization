function [ ElocTIRig,ElocTI] = LTI_SDP_DUAL(N,H,sdpOps)
%solve the locally translation invariant problem for a state of N spins

% support different input formats for Hamiltonian term H
if isstruct(H)
    suppH=H.suppH; % the number of sites H acts on 
    d = H.d; % the dimention of each site
    H=H.Hamiltonian; % the Hamiltonian itself
else % H is a d^2 x d^2 matrix
    suppH=2;
    d=sqrt(size(H,1));
end
 
 
alpha=sdpvar(d^(N-1)); %first state is on k0+1 sites. Alpha is the local transl. inv. dual variable

sdpvar epsil % dual variable of the tr(rho)==1 constraint

expr = kron(H,speye(d^(N-suppH))) ...
    + dontUseKron(alpha,d^(N-1),d,'IX') ...
    - dontUseKron(alpha,d^(N-1),d,'XI') ...
    - epsil*speye(d^(N)); 

constraints = [ expr >= 0 ];
 
objectiveFunc = - epsil; %minus sign for maximization

 optimize(constraints,objectiveFunc,sdpOps);
% %%
ElocTI=-value(objectiveFunc);
 
expr = value(expr);

% we can obtain a certified result by making sure that the dual variables obtained from the solver are feasible (to machine precision)
minEig = min(real(eig( expr )));
% the corrected energy value is then:
ElocTIRig=value(epsil) + minEig;
 
