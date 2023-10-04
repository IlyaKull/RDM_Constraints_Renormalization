function [energyRig,violation]=FeasibleDualPtSCSmod(u,PD)

xindsInU=Uinds(PD);
[xCell,numVarsX] = readVecToCell(u(xindsInU)/u(end)/PD.scale_sigma ,PD.xinds); % x=u(xindsInU)/u(end)/PD.scale_sigma;
epsil= -1*xCell{1}; %%%% !!! note -1
alpha= xCell{2};
betaL= xCell{3};
betaR= xCell{4};
gammaL=cell(PD.n-1,1);
gammaR=cell(PD.n-1,1);
j=5;k=PD.k0;
while j < numVarsX
    gammaL{k}=xCell{j};
    gammaR{k}=xCell{j+1};
    j=j+2;k=k+1;
end
[energyRig,violation]=FeasiblePtFromDualVar(alpha,betaL,betaR,gammaL,gammaR,epsil,...
                  PD.H,PD.V0xI,PD.IxV0,PD.LxI,PD.IxR,PD.k0,PD.n,PD.d,PD.D,PD.methodOps.whichMaps);


end

function [E,violation]=FeasiblePtFromDualVar(alpha,betaL,betaR,gammaL,gammaR,epsil,H,V0xI,IxV0,LxI,IxR,k0,n,d,D,isoOrMPS)
 
%%
MinEig=nan(n-1,1);
correctedMinEig=nan(n-1,1);
dim=D^2*d;

switch isoOrMPS
    case {'mpsLeft','mix'}
        Wl=0; % all the weight goes to the right because the R CG map is an isometry
    otherwise
        Wl=0.5;
end

%%
% the last constraint is:
%    -tensor(eye(d),gammaL{n-1}) - tensor(gammaR{n-1},eye(d)) >= 0

MinEig(n-1)=min(real(eig(...
                          -tensor(eye(d),gammaL{n-1}) - tensor(gammaR{n-1},eye(d))...
                          )));
% to satisfy the constraint we remove the minimum eigenvalue from the LHS
% by updating gammaL\R{n-1}
    gammaL{n-1}= gammaL{n-1} + Wl*eye(dim)*MinEig(n-1);
    gammaR{n-1}= gammaR{n-1} + (1-Wl)*eye(dim)*MinEig(n-1);
correctedMinEig(n-1)=min(real(eig(...
                          -tensor(eye(d),gammaL{n-1}) - tensor(gammaR{n-1},eye(d))...
                          )));
            

% second to last constraint and onwards
for k=n-1:-1:k0+1
     
    % the constraint is:
    % LxI'*gammaL{k}*LxI + IxR'*gammaR{k}*IxR - tensor(eye(d),gammaL{k-1}) - tensor(gammaR{k-1},eye(d)) >= 0      

    MinEig(k-1)=min(real(eig(...
      LxI{k}'*gammaL{k}*LxI{k} + IxR{k}'*gammaR{k}*IxR{k} - tensor(eye(d),gammaL{k-1}) - tensor(gammaR{k-1},eye(d))...
                      )));

        gammaL{k-1}= gammaL{k-1} + Wl*eye(dim)*MinEig(k-1);
        gammaR{k-1}= gammaR{k-1} + (1-Wl)*eye(dim)*MinEig(k-1);
    correctedMinEig(k-1)=min(real(eig(...
      LxI{k}'*gammaL{k}*LxI{k} + IxR{k}'*gammaR{k}*IxR{k} - tensor(eye(d),gammaL{k-1}) - tensor(gammaR{k-1},eye(d))...
                      )));


 end

% second constraint involves betaL\R
% LxI'*gammaL{k0}*LxI + IxR'*gammaR{k0}*IxR - tensor(eye(d),betaL) - tensor(betaR,eye(d))>= 0

MinEig(k0-1)=min(real(eig(...
    LxI{k0}'*gammaL{k0}*LxI{k0} + IxR{k0}'*gammaR{k0}*IxR{k0} - tensor(eye(d),betaL) - tensor(betaR,eye(d))...
                    )));

    betaL= betaL + Wl*eye(dim)*MinEig(k0-1);
    betaR= betaR + (1-Wl)*eye(dim)*MinEig(k0-1);
correctedMinEig(k0-1)=min(real(eig(...
    LxI{k0}'*gammaL{k0}*LxI{k0} + IxR{k0}'*gammaR{k0}*IxR{k0} - tensor(eye(d),betaL) - tensor(betaR,eye(d))...
                    )));

% the first constraint
% tensor(H,eye(d^(k0-1))) + V0xI'*betaL*V0xI + IxV0'*betaR*IxV0 + tensor(eye(d),alpha) - tensor(alpha,eye(d)) - epsil*eye(d^(k0+1)) >= 0 

 
 
Eshift=min(real(eig(...
     tensor(H,eye(d^(k0-1))) + V0xI'*betaL*V0xI + IxV0'*betaR*IxV0 + tensor(eye(d),alpha) - tensor(alpha,eye(d)) - epsil*eye(d^(k0+1))...
                    )));


E= epsil+Eshift;  

FirstConstrCorrectedMinEig=min(real(eig(...
    tensor(H,eye(d^(k0-1))) + V0xI'*betaL*V0xI + IxV0'*betaR*IxV0 + tensor(eye(d),alpha) - tensor(alpha,eye(d)) - E*eye(d^(k0+1))...
                    )));
                
% only the negative corrected eigs violate the PSD conditions
violation= sum(correctedMinEig(correctedMinEig<0)) + (FirstConstrCorrectedMinEig<0)*FirstConstrCorrectedMinEig;
 % adding E to epsil will make the first inequality hold which makes the
 % updated dual point (alpha,betaL,betaR,gammaL,gammaR,E) feasible
 
end
