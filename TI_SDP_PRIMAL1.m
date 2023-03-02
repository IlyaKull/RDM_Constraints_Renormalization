function [ Elower,sol,omega,k0] = TI_SDP_PRIMAL1(n,D,H,vuMPS,methodOps,yalmipOps)

d=sqrt(size(H,1));
 
%% sdp variables
[rho,omega,dims,k0]=declarePrimalVars(n,d,D);
if k0>=n 
   Elower=nan;
   sol.solvertime=nan;
   sol.problem='k0>=n, skipped this n';
   sol.info=[];
   warning('k0>=n, skipped this n');
   omega=[];
   return
end

% keep only last rho:
for k=1:k0-2
    rho{k}=[];
end
%% isometries from mps
 [V0,L,R] = CoarseGrainingMaps(n,k0,vuMPS,methodOps);
%% positivity constraints
constraints=[];

constraints = [constraints,...
                   rho{k0-1} >= 0];


    for l=k0:n
        constraints = [constraints,...
                       omega{l} >= 0];
    end

%% local compatability constraint for rho{k0-1}
constraints= [constraints,...
             pTr(rho{k0-1},1,dims{k0-1}) == pTr(rho{k0-1},3,dims{k0-1}) ];
%% normalization 
constraints= [constraints,...
              trace(rho{k0-1}) == 1];
%% objective function (evaluate 2-body hamiltonian on rho{k0-1}
objectiveFunc=trace(...
                   rho{k0-1} * tensor(H,eye(prod(dims{k0-1})/d^2)) ...
                   );

%% compatibility constraints k=k0
V0xI=tensor(V0,eye(d)); 
IxV0=tensor(eye(d),V0); 
constraints = [constraints,...
                V0xI * rho{k0-1} * V0xI' == pTr(omega{k0},1,dims{k0}),...
                IxV0 * rho{k0-1} * IxV0' == pTr(omega{k0},3,dims{k0})...
               ];
 
%% compatibility constraints k>k0
for k=k0:n-1
    LxI=tensor(L{k},eye(d));
    IxR=tensor(eye(d),R{k});
   constraints = [constraints,...
                 LxI * omega{k} * LxI' == pTr(omega{k+1},1,dims{k+1}),...
                 IxR * omega{k} * IxR' == pTr(omega{k+1},3,dims{k+1})];
end



%% optimize
%   
sol = optimize(constraints,objectiveFunc,yalmipOps);
%%
Elower=value(objectiveFunc);

 
