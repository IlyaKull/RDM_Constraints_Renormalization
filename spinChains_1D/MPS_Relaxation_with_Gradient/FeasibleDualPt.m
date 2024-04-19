function [E]=FeasibleDualPt(n,D,k0,H,sdpVars,CGmaps, CGops)
                       
% takes the dual variables returned by the solver (approximately feasible)
% and corrects them to make the point feasible. The energy returned is thus
% a rigorous lower bound.

% the correction is done by starting at the last constraint and adding an
% identity term to gammaL/R{n-1} to make the inequality hold. then
% correcting the constraint before by adding an identity to gammaL/R{n-2}
% and so on. The correction can be done with different weight (Wl) between gammaL/R

suppH=H.suppH;
d=H.d;
H=H.Hamiltonian;

%% obtain the values from the desicion variables
alpha=value(sdpVars.alpha);
betaL=value(sdpVars.betaL);
betaR=value(sdpVars.betaR);
gammaL=cell(n-1,1);
gammaR=cell(n-1,1);
for j=k0:n-1
    gammaL{j}=value(sdpVars.gammaL{j});
    gammaR{j}=value(sdpVars.gammaR{j});
end
epsil=value(sdpVars.epsil);
%%
MinEig=nan(n-1,1);
correctedMinEig=nan(n-1,1);
dim=D^2*d;

switch CGops.whichMaps
    case {'mpsLeft','mixLeft'}
        Wl=0; % all the weight goes to the right because the R CG map is an isometry
    otherwise
        Wl=0.5;
end

%%
if n>k0
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
        LxI=tensor(CGmaps.L{k},eye(d));
        IxR=tensor(eye(d),CGmaps.R{k});
        % the constraint is:
        % LxI'*gammaL{k}*LxI + IxR'*gammaR{k}*IxR - tensor(eye(d),gammaL{k-1}) - tensor(gammaR{k-1},eye(d)) >= 0      

        MinEig(k-1)=min(real(eig(...
          LxI'*gammaL{k}*LxI + IxR'*gammaR{k}*IxR - tensor(eye(d),gammaL{k-1}) - tensor(gammaR{k-1},eye(d))...
                          )));

            gammaL{k-1}= gammaL{k-1} + Wl*eye(dim)*MinEig(k-1);
            gammaR{k-1}= gammaR{k-1} + (1-Wl)*eye(dim)*MinEig(k-1);
        correctedMinEig(k-1)=min(real(eig(...
          LxI'*gammaL{k}*LxI + IxR'*gammaR{k}*IxR - tensor(eye(d),gammaL{k-1}) - tensor(gammaR{k-1},eye(d))...
                          )));


     end

    % second constraint involves betaL\R
    LxI=tensor(CGmaps.L{k0},eye(d));
    IxR=tensor(eye(d),CGmaps.R{k0});

    % LxI'*gammaL{k0}*LxI + IxR'*gammaR{k0}*IxR - tensor(eye(d),betaL) - tensor(betaR,eye(d))>= 0


    MinEig(k0-1)=min(real(eig(...
        LxI'*gammaL{k0}*LxI + IxR'*gammaR{k0}*IxR - tensor(eye(d),betaL) - tensor(betaR,eye(d))...
                        )));


        betaL= betaL + Wl*eye(dim)*MinEig(k0-1);
        betaR= betaR + (1-Wl)*eye(dim)*MinEig(k0-1);
    correctedMinEig(k0-1)=min(real(eig(...
        LxI'*gammaL{k0}*LxI + IxR'*gammaR{k0}*IxR - tensor(eye(d),betaL) - tensor(betaR,eye(d))...
                        )));

else %k0==n
    % second constraint involves betaL\R
 MinEig(k0-1)=min(real(eig(...
          - tensor(eye(d),betaL) - tensor(betaR,eye(d))...
                        )));


        betaL= betaL + Wl*eye(dim)*MinEig(k0-1);
        betaR= betaR + (1-Wl)*eye(dim)*MinEig(k0-1);
    correctedMinEig(k0-1)=min(real(eig(...
          - tensor(eye(d),betaL) - tensor(betaR,eye(d))...
                        )));
% the first constraint
% tensor(H,eye(d^(k0-1))) + V0xI'*betaL*V0xI + IxV0'*betaR*IxV0 + tensor(eye(d),alpha) - tensor(alpha,eye(d)) - epsil*eye(d^(k0+1)) >= 0 
end 

 V0xI=kron(CGmaps.V0,speye(d));
 IxV0=kron(speye(d),CGmaps.V0);
 
Eshift=min(real(eig(...
     tensor(H,eye(d^(k0-(suppH-1)))) + V0xI'*betaL*V0xI + IxV0'*betaR*IxV0 + tensor(eye(d),alpha) - tensor(alpha,eye(d)) - epsil*eye(d^(k0+1))...
                   )));

E= epsil+Eshift;  
                



